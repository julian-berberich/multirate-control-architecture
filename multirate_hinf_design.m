%% MULTIRATE_HINF_SINGLE
%
% Single-point robust L2-gain study of a multirate sampled-data loop with two
% independent, asynchronous zero-order holds and a polytopic uncertainty in
% the slow actuator. Three quantities are computed at the same dwell-time
% bounds, with the same relaxation, so the numbers are comparable:
%
%   (1) gamma for the GIVEN gains K1 = [-4 0], K2 = [0 -2], via the primal
%       (analysis) formulation: X(tau1,tau2) parameterized as a spline.
%   (2) gamma for the same given gains via the DUAL formulation, obtained by
%       substituting L_Ji := Ki * Y, which stays LINEAR in Y because Ki is a
%       constant. Y = X^{-1} is parameterized as a spline.
%   (3) gamma for DESIGNED clock-scheduled gains: same dual LMIs, L_Ji free.
%
% TWO BUILT-IN CONSISTENCY CHECKS
%   * (3) <= (2) is RIGOROUS. Same parameterization, and the fixed-gain
%     feasible set is literally a subset of the design one (L tied to Y vs.
%     L free). A violation is a bug, not conservatism.
%   * (1) vs (2) is a CROSS-FORMULATION check on the same quantity through
%     independent code paths. They are not identical at a fixed mesh: a
%     spline's inverse is not a spline, so the two parameterize different
%     function classes. They should agree closely and converge together
%     under knot refinement -- see the refinement check at the end.
%
% MODEL
%   Plant    xdot = Ap x + B1 u1 + B2(beta) u2 + Ep d,   x in R^2
%   Loop 1 (fast) actuates x1, loop 2 (slow) actuates x2, the disturbance
%   enters the slowly actuated channel, and beta in [beta_min, beta_max] is
%   the uncertain slow-actuator effectiveness.
%
%   Augmented state z = [x; u1_hold; u2_hold]:
%       flow    zdot = A(beta) z + Be d,   e = C z + De d
%       jump i  z^+  = A_Ji z + B_Ji * (K_Ji z)
%
%   There are NO jump performance channels: x^+ = x^-, so a process
%   disturbance has nothing to enter through at a reset and no signal is
%   emitted there. Consequently the jump LMI coincides with the pure
%   stability jump condition -- gamma enters only the flow LMI, while the
%   jump conditions are what tie gamma to the sampling rates.
%
% NOTATION
%   beta = uncertain slow-actuator effectiveness (plant parameter).
%   gamma = certified L2 gain. (The symbol nu, used in the sweep script for
%   the sampling budget, does not appear here.)
%
% NUMERICAL CONVENTIONS (carried over from the sweep script)
%   * cfg.tol is the LMI strictness margin and cfg.resTol the residual
%     acceptance threshold, with |resTol| < tol. An ACCEPTED solution then
%     satisfies M <= -(tol-|resTol|)*I, i.e. genuine strict feasibility.
%   * NO regularization is placed on L when gamma is the reported quantity.
%     A penalty rho*sum(L.^2) with rho = 1e-4 contributes O(1)-O(100) to the
%     objective for these control-point counts and would bias gamma upward.
%     L is bounded by a box instead. A gain penalty IS used in the optional
%     fixed-gamma design, where gamma is prescribed rather than reported.
%   * The bounds on Y (resp. X) are loose numerical guards only. A tight
%     lower bound such as Y >= 1e-2*I caps X = Y^{-1} <= 100*I, which is a
%     real restriction on the certificate and can produce false infeasibility.
%     The eigenvalue range actually attained is reported.
%   * Knots place an EXACT breakpoint at T_min on each axis. The convex-hull
%     argument certifies an LMI on a patch only if it is imposed at every
%     control point of that patch, so selecting control points by a test such
%     as "grid value >= T_min" is unsound unless T_min is a knot.
%   * Segment counts are chosen adaptively so every segment satisfies
%     ||A||*h <= alpha: required resolution scales with ||A||*T, not with T.
%
% Requires: YALMIP and an SDP solver (MOSEK by default).

clear; close all; clc;
rng(1);

%% ------------------------------------------------------------------------
%  Options
%  ------------------------------------------------------------------------

DO_SPEC_DESIGN = false;    % after minimizing, re-solve at a prescribed gamma
gamma_spec     = NaN;      % target for that re-solve (set a number to use)

DO_VERIFY      = true;     % dense-grid re-check of every certificate
DO_SIMULATE    = true;     % closed-loop simulation + empirical gain
DO_REFINEMENT  = true;     % re-solve at halved alpha and report the change

%% ------------------------------------------------------------------------
%  Plant, uncertainty, given gains
%  ------------------------------------------------------------------------

cfg.Ap    = [ 2.0  0.50;
             -0.5  0.25];          % open loop unstable
cfg.B1    = [1; 0];                % fast loop actuates x1
cfg.B2fun = @(b) [0; b];           % slow loop, uncertain effectiveness
cfg.Ep    = [0; 1];                % disturbance enters x2

cfg.beta_min = 0.6;
cfg.beta_max = 1.0;

% Gains from the original stability example, designed WITHOUT any robustness
% or performance requirement. Rows act on x; padded with zeros for the two
% hold states of the augmented model.
K1 = [-4  0];
K2 = [ 0 -2];

%% ------------------------------------------------------------------------
%  Dwell-time bounds and performance weights
%  ------------------------------------------------------------------------

cfg.T1min = 0.03;   cfg.T1max = 0.05;      % fast loop
cfg.T2min = 0.30;   cfg.T2max = 0.50;      % slow loop

cfg.w1  = 1.0;      % weight on x1
cfg.w2  = 1.0;      % weight on x2
cfg.rho = 0.1;      % weight on control effort

%% ------------------------------------------------------------------------
%  Relaxation and numerics
%  ------------------------------------------------------------------------

cfg.d1 = 2;   cfg.d2 = 2;          % Bezier degree per axis
cfg.alpha = 0.40;                  % require ||A|| * h <= alpha
cfg.s_max = 8;                     % cap on segments per region

cfg.tol    = 1e-5;                 % LMI strictness margin
cfg.resTol = -1e-6;                % residual acceptance (|resTol| < tol)
cfg.Ymin   = 1e-6;   cfg.Ymax = 1e6;   % loose guards on Y (dual)
cfg.Xmin   = 1e-6;   cfg.Xmax = 1e6;   % loose guards on X (primal)
cfg.Lmax   = 5e2;                  % box on L (no penalty when gamma is free)
cfg.regL   = 1e-3;                 % gain penalty, ONLY in fixed-gamma design

cfg.options = sdpsettings('solver','mosek','verbose',0);
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-10;
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-10;

%% ------------------------------------------------------------------------
%  Augmented model
%  ------------------------------------------------------------------------

n      = size(cfg.Ap,1) + 2;
Averts = {augA(cfg, cfg.beta_min), augA(cfg, cfg.beta_max)};

Be = [cfg.Ep; 0; 0];
C  = diag([cfg.w1 cfg.w2 cfg.rho cfg.rho]);
De = zeros(size(C,1), size(Be,2));
nd = size(Be,2);
ne = size(C,1);

B_J1 = [0;0;1;0];   A_J1 = diag([1 1 0 1]);
B_J2 = [0;0;0;1];   A_J2 = diag([1 1 1 0]);

Kf1 = [K1, 0, 0];                  % gains on the augmented state
Kf2 = [K2, 0, 0];
AK1 = A_J1 + B_J1*Kf1;             % fixed-gain closed-loop jump maps
AK2 = A_J2 + B_J2*Kf2;

mdl = struct('n',n,'Averts',{Averts},'Be',Be,'C',C,'De',De,'nd',nd,'ne',ne, ...
             'A_J1',A_J1,'A_J2',A_J2,'B_J1',B_J1,'B_J2',B_J2, ...
             'Kf1',Kf1,'Kf2',Kf2,'AK1',AK1,'AK2',AK2);

Anorm    = max(norm(Averts{1}), norm(Averts{2}));
cfg.hmax = cfg.alpha / Anorm;

%% ------------------------------------------------------------------------
%  Knots
%  ------------------------------------------------------------------------

[kn1, idxJ1, N1, s1] = make_axis(cfg.T1min, cfg.T1max, cfg.d1, cfg);
[kn2, idxJ2, N2, s2] = make_axis(cfg.T2min, cfg.T2max, cfg.d2, cfg);

fprintf('=====================================================================\n');
fprintf(' MULTIRATE ROBUST L2-GAIN, SINGLE OPERATING POINT\n');
fprintf('=====================================================================\n');
fprintf('  dwell times    T1 in [%.3f, %.3f] s,  T2 in [%.3f, %.3f] s\n', ...
        cfg.T1min, cfg.T1max, cfg.T2min, cfg.T2max);
fprintf('  uncertainty    beta in [%.2f, %.2f], %d vertices\n', ...
        cfg.beta_min, cfg.beta_max, numel(Averts));
fprintf('  ||A|| = %.4f  ->  max segment width h = %.5f s\n', Anorm, cfg.hmax);
fprintf('  axis 1: %d+%d segments, degree %d, %d control points\n', s1(1), s1(2), cfg.d1, N1);
fprintf('  axis 2: %d+%d segments, degree %d, %d control points\n', s2(1), s2(2), cfg.d2, N2);
fprintf('  tensor grid: %d x %d = %d matrix blocks of size %d\n\n', N1, N2, N1*N2, n);

kn = struct('k1',kn1,'k2',kn2,'i1',idxJ1,'i2',idxJ2,'N1',N1,'N2',N2);

%% ------------------------------------------------------------------------
%  (1) Given gains, PRIMAL formulation (Theorem dissipativity)
%  ------------------------------------------------------------------------

fprintf('--- (1) given gains, primal formulation ------------------------------\n');
[gFixP, Xv, dFixP] = solve_primal(cfg, mdl, kn, NaN);
report_solve(gFixP, dFixP, cfg);

%% ------------------------------------------------------------------------
%  (2) Given gains, DUAL formulation (L_Ji := Ki * Y, linear in Y)
%  ------------------------------------------------------------------------

fprintf('--- (2) given gains, dual formulation --------------------------------\n');
[gFixD, YvF, L1vF, L2vF, dFixD] = solve_dual(cfg, mdl, kn, true, NaN);
report_solve(gFixD, dFixD, cfg);

%% ------------------------------------------------------------------------
%  (3) Designed clock-scheduled gains (dual, L free)
%  ------------------------------------------------------------------------

fprintf('--- (3) designed gains, dual formulation -----------------------------\n');
[gDes, Yv, L1v, L2v, dDes] = solve_dual(cfg, mdl, kn, false, NaN);
report_solve(gDes, dDes, cfg);

%% ------------------------------------------------------------------------
%  Optional: re-solve the design at a prescribed gamma, minimizing gains
%  ------------------------------------------------------------------------

if DO_SPEC_DESIGN && isfinite(gamma_spec)
    if ~isnan(gDes) && gamma_spec < gDes
        fprintf(['--- spec design SKIPPED: gamma_spec = %.4f is below the\n' ...
                 '    achievable minimum %.4f.\n\n'], gamma_spec, gDes);
    else
        fprintf('--- spec design at gamma = %.4f -------------------------------\n', gamma_spec);
        [gSpec, Yv2, L1v2, L2v2, dSpec] = solve_dual(cfg, mdl, kn, false, gamma_spec^2);
        report_solve(gSpec, dSpec, cfg);
        if ~isnan(gSpec)
            Yv = Yv2;  L1v = L1v2;  L2v = L2v2;   % verify/simulate this one
        end
    end
end

%% ------------------------------------------------------------------------
%  Consistency checks
%  ------------------------------------------------------------------------

fprintf('--- consistency checks -----------------------------------------------\n');

if ~isnan(gDes) && ~isnan(gFixD)
    okMono = gDes <= gFixD*(1 + 1e-6);
    fprintf('  design <= fixed (same dual parameterization) : %s   %.4f vs %.4f\n', ...
            tf2str(okMono), gDes, gFixD);
    if ~okMono
        warning(['Design gamma exceeds fixed-gain gamma under the SAME ' ...
                 'parameterization. The fixed-gain feasible set is a subset ' ...
                 'of the design one, so this indicates a bug, not conservatism.']);
    end
end

if ~isnan(gFixP) && ~isnan(gFixD)
    relgap = abs(gFixP - gFixD) / max(gFixP, gFixD);
    fprintf('  primal vs dual, same fixed gains             : %.2f%% relative gap\n', 100*relgap);
    fprintf('    (not expected to be zero: X-spline and Y-spline are different\n');
    fprintf('     function classes; they converge under refinement)\n');
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  Recovered clock-scheduled gains
%  ------------------------------------------------------------------------

Yat  = @(t1,t2) eval_bezier2D(Yv,  t1, t2, kn1, kn2, cfg.d1, cfg.d2);
L1at = @(t1,t2) eval_bezier2D(L1v, t1, t2, kn1, kn2, cfg.d1, cfg.d2);
L2at = @(t1,t2) eval_bezier2D(L2v, t1, t2, kn1, kn2, cfg.d1, cfg.d2);
K1at = @(t1,t2) L1at(t1,t2) / Yat(t1,t2);
K2at = @(t1,t2) L2at(t1,t2) / Yat(t1,t2);

if ~isnan(gDes)
    fprintf('--- recovered gains (design) -----------------------------------------\n');
    fprintf('  fixed  K_J1 = [%7.3f %7.3f %7.3f %7.3f]  (given, constant)\n', Kf1);
    fprintf('  design K_J1(T1min, 0     ) = [%7.3f %7.3f %7.3f %7.3f]\n', K1at(cfg.T1min,0));
    fprintf('  design K_J1(T1max, T2max ) = [%7.3f %7.3f %7.3f %7.3f]\n', K1at(cfg.T1max,cfg.T2max));
    fprintf('  fixed  K_J2 = [%7.3f %7.3f %7.3f %7.3f]  (given, constant)\n', Kf2);
    fprintf('  design K_J2(0,      T2min) = [%7.3f %7.3f %7.3f %7.3f]\n', K2at(0,cfg.T2min));
    fprintf('  design K_J2(T1max,  T2max) = [%7.3f %7.3f %7.3f %7.3f]\n\n', K2at(cfg.T1max,cfg.T2max));
end

%% ------------------------------------------------------------------------
%  Dense-grid verification
%
%  Both certificates are re-checked in the ORIGINAL coordinates X, using the
%  recovered gains. The convex-hull argument already certifies them exactly;
%  this is a safeguard against construction errors (misaligned knots, wrong
%  index sets) and is cheap enough to always run.
%  ------------------------------------------------------------------------

if DO_VERIFY
    if ~isnan(gFixP)
        XatP = @(t1,t2) eval_bezier2D(Xv, t1, t2, kn1, kn2, cfg.d1, cfg.d2);
        verify_certificate('(1) given gains, primal', XatP, ...
            @(t1,t2) AK1, @(t1,t2) AK2, cfg, mdl, kn, gFixP^2);
    end
    if ~isnan(gFixD)
        YatF  = @(t1,t2) eval_bezier2D(YvF, t1, t2, kn1, kn2, cfg.d1, cfg.d2);
        verify_certificate('(2) given gains, dual', @(t1,t2) inv(YatF(t1,t2)), ...
            @(t1,t2) AK1, @(t1,t2) AK2, cfg, mdl, kn, gFixD^2);
    end
    if ~isnan(gDes)
        verify_certificate('(3) designed gains, dual', @(t1,t2) inv(Yat(t1,t2)), ...
            @(t1,t2) A_J1 + B_J1*K1at(t1,t2), ...
            @(t1,t2) A_J2 + B_J2*K2at(t1,t2), cfg, mdl, kn, gDes^2);
    end
end

%% ------------------------------------------------------------------------
%  Closed-loop simulation and empirical gain
%
%  Randomized sampling intervals, windowed sinusoidal disturbances over a
%  frequency grid, at both uncertainty vertices. The empirical ratio is a
%  LOWER bound on the true induced gain and must sit below gamma.
%  ------------------------------------------------------------------------

if DO_SIMULATE
    Tfinal = 8;  freqs = [0.2 0.5 1 2 5];
    fprintf('--- empirical gain estimates (lower bounds) --------------------------\n');
    if ~isnan(gFixP) || ~isnan(gFixD)
        rF = empirical_gain(cfg, mdl, @(t1,t2) Kf1, @(t1,t2) Kf2, Tfinal, freqs);
        fprintf('  given gains    : %.4f   (certified %.4f / %.4f)\n', ...
                rF, gFixP, gFixD);
    end
    if ~isnan(gDes)
        rD = empirical_gain(cfg, mdl, K1at, K2at, Tfinal, freqs);
        fprintf('  designed gains : %.4f   (certified %.4f)\n', rD, gDes);
    end
    fprintf('\n');

    dfun = @(t) (t <= 2) .* sin(1.0*t);
    [tF, ZF] = simulate_cl(Averts{1}, Be, A_J1, A_J2, B_J1, B_J2, ...
                           @(t1,t2) Kf1, @(t1,t2) Kf2, cfg, Tfinal, dfun);
    if ~isnan(gDes)
        [tD, ZD] = simulate_cl(Averts{1}, Be, A_J1, A_J2, B_J1, B_J2, ...
                               K1at, K2at, cfg, Tfinal, dfun);
    end

    figure('Name','States','Color','w');
    plot(tF, ZF(:,1), 'k--', 'LineWidth',2); hold on
    plot(tF, ZF(:,2), 'k-',  'LineWidth',2);
    if ~isnan(gDes)
        plot(tD, ZD(:,1), 'b--', 'LineWidth',2);
        plot(tD, ZD(:,2), 'b-',  'LineWidth',2);
        legend({'$x_1$ given','$x_2$ given','$x_1$ design','$x_2$ design'}, ...
               'Interpreter','latex','FontSize',12);
    else
        legend({'$x_1$ given','$x_2$ given'},'Interpreter','latex','FontSize',12);
    end
    grid on
    xlabel('time $t$','Interpreter','latex','FontSize',18);
    ylabel('states','Interpreter','latex','FontSize',18);
    set(gca,'FontSize',13);

    figure('Name','Inputs','Color','w');
    plot(tF, ZF(:,3), 'k--', 'LineWidth',2); hold on
    plot(tF, ZF(:,4), 'k-',  'LineWidth',2);
    if ~isnan(gDes)
        plot(tD, ZD(:,3), 'b--', 'LineWidth',2);
        plot(tD, ZD(:,4), 'b-',  'LineWidth',2);
    end
    grid on
    xlabel('time $t$','Interpreter','latex','FontSize',18);
    ylabel('inputs','Interpreter','latex','FontSize',18);
    set(gca,'FontSize',13);
end

%% ------------------------------------------------------------------------
%  Refinement check
%
%  cfg.alpha is a heuristic, not a convergence guarantee. Halving it should
%  move gamma only marginally; the primal and dual fixed-gain values should
%  also move CLOSER to each other.
%  ------------------------------------------------------------------------

if DO_REFINEMENT
    fprintf('--- refinement check (alpha halved) ----------------------------------\n');
    cfg2 = cfg;  cfg2.alpha = cfg.alpha/2;  cfg2.hmax = cfg2.alpha/Anorm;
    [k1b, j1b, N1b, s1b] = make_axis(cfg.T1min, cfg.T1max, cfg.d1, cfg2);
    [k2b, j2b, N2b, s2b] = make_axis(cfg.T2min, cfg.T2max, cfg.d2, cfg2);
    knb = struct('k1',k1b,'k2',k2b,'i1',j1b,'i2',j2b,'N1',N1b,'N2',N2b);
    fprintf('  axis 1: %d+%d segments (%d pts), axis 2: %d+%d segments (%d pts)\n', ...
            s1b(1), s1b(2), N1b, s2b(1), s2b(2), N2b);

    gFixP2 = solve_primal(cfg2, mdl, knb, NaN);
    gFixD2 = solve_dual(cfg2, mdl, knb, true,  NaN);
    gDes2  = solve_dual(cfg2, mdl, knb, false, NaN);

    prt = @(nm,a,b) fprintf('  %-28s %8.4f -> %8.4f   (%+.2f%%)\n', ...
                            nm, a, b, 100*(b-a)/max(a,eps));
    prt('given gains, primal', gFixP, gFixP2);
    prt('given gains, dual',   gFixD, gFixD2);
    prt('designed gains',      gDes,  gDes2);
    if ~isnan(gFixP2) && ~isnan(gFixD2)
        fprintf('  primal-vs-dual gap: %.2f%% -> %.2f%%\n', ...
                100*abs(gFixP-gFixD)/max(gFixP,gFixD), ...
                100*abs(gFixP2-gFixD2)/max(gFixP2,gFixD2));
    end
    fprintf('\n');
end

%% ------------------------------------------------------------------------
%  Summary
%  ------------------------------------------------------------------------

fprintf('=====================================================================\n');
fprintf(' SUMMARY   (T1 in [%.3f,%.3f], T2 in [%.3f,%.3f], beta in [%.2f,%.2f])\n', ...
        cfg.T1min, cfg.T1max, cfg.T2min, cfg.T2max, cfg.beta_min, cfg.beta_max);
fprintf('---------------------------------------------------------------------\n');
fprintf('  given gains K1=[%g %g], K2=[%g %g]\n', K1, K2);
fprintf('    primal formulation                  gamma = %s\n', n2s(gFixP));
fprintf('    dual formulation                    gamma = %s\n', n2s(gFixD));
fprintf('  designed clock-scheduled gains        gamma = %s\n', n2s(gDes));
if ~isnan(gDes) && ~isnan(gFixD)
    fprintf('  improvement from design               %.2fx\n', gFixD/gDes);
end
fprintf('=====================================================================\n');

save('hinf_single_results.mat','cfg','mdl','kn','gFixP','gFixD','gDes', ...
     'Xv','Yv','L1v','L2v','YvF','K1','K2');


%% ========================================================================
%  Local functions -- model and knots
%  ========================================================================

function A = augA(cfg, beta)
    np = size(cfg.Ap,1);
    A  = [cfg.Ap, cfg.B1, cfg.B2fun(beta);
          zeros(2, np+2)];
end

function [knots, idxJump, N, s] = make_axis(Tmin, Tmax, d, cfg)
% Adaptive segment counts plus a knot vector with an EXACT breakpoint at Tmin.
% Resolution scales with ||A||*T, so each region is refined until every
% segment satisfies ||A||*h <= alpha.
    need_pre  = Tmin          / cfg.hmax;
    need_post = (Tmax - Tmin) / cfg.hmax;
    s_pre  = min(cfg.s_max, max(1, ceil(need_pre)));
    s_post = min(cfg.s_max, max(1, ceil(need_post)));
    if (need_pre > cfg.s_max) || (need_post > cfg.s_max)
        warning('make_axis:capped', ...
            ['Segment cap s_max = %d bound on [%.4f, %.4f]; an infeasible ' ...
             'result here is INCONCLUSIVE.'], cfg.s_max, Tmin, Tmax);
    end
    [knots, idxJump] = build_knots(Tmin, Tmax, s_pre, s_post, d);
    N = (s_pre + s_post)*d + 1;
    s = [s_pre s_post];
end

function [knots, idxJump] = build_knots(Tmin, Tmax, s_pre, s_post, d)
% Knot vector on [0,Tmax] with an exact breakpoint at Tmin, so the
% jump-constrained control points form a union of whole Bezier patches.
    assert(Tmin > 0 && Tmax > Tmin, 'Require 0 < Tmin < Tmax.');
    pre   = linspace(0, Tmin, s_pre + 1);
    post  = linspace(Tmin, Tmax, s_post + 1);
    knots = [pre, post(2:end)];
    idxJump = (s_pre*d + 1) : ((s_pre + s_post)*d + 1);
end

%% ========================================================================
%  Local functions -- solvers
%  ========================================================================

function [gamma, Xval, dg] = solve_primal(cfg, mdl, kn, gam2fix)
% Given gains, primal formulation (Theorem dissipativity):
%   [ Xdot + A'X + XA + C'C   X Be + C'De ]
%   [        *               -g^2 I + De'De ] < 0        per uncertainty vertex
%   A_Ji^K' X^{0,i} A_Ji^K - X < 0                       on T_i
% Linear in (X, gamma^2).

    yalmip('clear');
    gamma = NaN;  Xval = {};  dg = nan(1,4);

    n = mdl.n;  nd = mdl.nd;
    d1 = cfg.d1;  d2 = cfg.d2;
    m1 = numel(kn.k1) - 1;   m2 = numel(kn.k2) - 1;

    Xc = cell(kn.N1, kn.N2);
    for p = 1:kn.N1
        for q = 1:kn.N2
            Xc{p,q} = sdpvar(n, n, 'symmetric');
        end
    end

    if isnan(gam2fix)
        gam2 = sdpvar(1,1);  minimizing = true;
    else
        gam2 = gam2fix;      minimizing = false;
    end

    Fc = {};
    if minimizing, Fc{end+1} = (gam2 >= 0); end

    for p = 1:kn.N1
        for q = 1:kn.N2
            Fc{end+1} = (Xc{p,q} >= cfg.Xmin*eye(n));
            Fc{end+1} = (Xc{p,q} <= cfg.Xmax*eye(n));
        end
    end

    Ifl = eye(n + nd);
    for v = 1:numel(mdl.Averts)
        A = mdl.Averts{v};
        for seg1 = 1:m1
            h1 = kn.k1(seg1+1) - kn.k1(seg1);
            for seg2 = 1:m2
                h2 = kn.k2(seg2+1) - kn.k2(seg2);
                [Cloc, dC] = bezier_patch_derivative(Xc, seg1, seg2, d1, d2, h1, h2);
                for i = 0:d1
                    for j = 0:d2
                        Xt  = Cloc{i+1,j+1};
                        dXt = dC{i+1,j+1};
                        M = [ dXt + A'*Xt + Xt*A + mdl.C'*mdl.C, ...
                                    Xt*mdl.Be + mdl.C'*mdl.De ;
                              mdl.Be'*Xt + mdl.De'*mdl.C, ...
                                   -gam2*eye(nd) + mdl.De'*mdl.De ];
                        Fc{end+1} = (M <= -cfg.tol*Ifl);
                    end
                end
            end
        end
    end

    for p = kn.i1
        for q = 1:kn.N2
            Fc{end+1} = (mdl.AK1'*Xc{1,q}*mdl.AK1 - Xc{p,q} <= -cfg.tol*eye(n));
        end
    end
    for q = kn.i2
        for p = 1:kn.N1
            Fc{end+1} = (mdl.AK2'*Xc{p,1}*mdl.AK2 - Xc{p,q} <= -cfg.tol*eye(n));
        end
    end

    F = [Fc{:}];
    if minimizing, obj = gam2; else, obj = 0; end

    try
        sol = optimize(F, obj, cfg.options);
    catch ME
        warning('solve_primal:solverError','%s', ME.message);
        dg(3) = -99;  return
    end

    dg(3) = sol.problem;
    if sol.problem ~= 0, return; end
    dg(1) = min(check(F));
    if dg(1) <= cfg.resTol, return; end

    if minimizing, gamma = sqrt(max(value(gam2),0)); else, gamma = sqrt(gam2fix); end
    Xval = cellfun(@value, Xc, 'UniformOutput', false);
    [dg(2), dg(4)] = eig_range(Xval);
end

function [gamma, Yval, L1val, L2val, dg] = solve_dual(cfg, mdl, kn, fixedGains, gam2fix)
% Dual (inverse-Lyapunov) formulation, Proposition 2:
%   [ -Ydot + A Y + Y A'   Be      Y C' ]
%   [        *            -g^2 I   De'  ] < 0             per uncertainty vertex
%   [        *              *      -I   ]
%   [ -Y                (A_Ji Y + B_Ji L_Ji)' ]
%   [  A_Ji Y + B_Ji L_Ji   -Y^{0,i}          ] < 0       on T_i
%
% fixedGains = true  : substitute L_Ji := Ki * Y. Since Ki is a constant this
%                      is LINEAR in Y, and the jump block becomes AKi*Y, whose
%                      Schur complement is exactly the primal jump condition
%                      with X = Y^{-1}. The feasible set is a SUBSET of the
%                      design one, so the resulting gamma is an upper bound
%                      for the designed gamma.
% fixedGains = false : L_Ji free (design).

    yalmip('clear');
    gamma = NaN;  Yval = {};  L1val = {};  L2val = {};  dg = nan(1,4);

    n = mdl.n;  nd = mdl.nd;  ne = mdl.ne;
    d1 = cfg.d1;  d2 = cfg.d2;
    m1 = numel(kn.k1) - 1;   m2 = numel(kn.k2) - 1;

    Yc = cell(kn.N1, kn.N2);
    for p = 1:kn.N1
        for q = 1:kn.N2
            Yc{p,q} = sdpvar(n, n, 'symmetric');
        end
    end

    L1c = cell(kn.N1, kn.N2);
    L2c = cell(kn.N1, kn.N2);
    for p = 1:kn.N1
        for q = 1:kn.N2
            L1c{p,q} = zeros(1,n);
            L2c{p,q} = zeros(1,n);
        end
    end
    if fixedGains
        for p = kn.i1
            for q = 1:kn.N2, L1c{p,q} = mdl.Kf1 * Yc{p,q}; end
        end
        for q = kn.i2
            for p = 1:kn.N1, L2c{p,q} = mdl.Kf2 * Yc{p,q}; end
        end
    else
        for p = kn.i1
            for q = 1:kn.N2, L1c{p,q} = sdpvar(1, n, 'full'); end
        end
        for q = kn.i2
            for p = 1:kn.N1, L2c{p,q} = sdpvar(1, n, 'full'); end
        end
    end

    if isnan(gam2fix)
        gam2 = sdpvar(1,1);  minimizing = true;
    else
        gam2 = gam2fix;      minimizing = false;
    end

    Fc = {};
    if minimizing, Fc{end+1} = (gam2 >= 0); end

    for p = 1:kn.N1
        for q = 1:kn.N2
            Fc{end+1} = (Yc{p,q} >= cfg.Ymin*eye(n));
            Fc{end+1} = (Yc{p,q} <= cfg.Ymax*eye(n));
        end
    end

    Ifl = eye(n + nd + ne);
    for v = 1:numel(mdl.Averts)
        A = mdl.Averts{v};
        for seg1 = 1:m1
            h1 = kn.k1(seg1+1) - kn.k1(seg1);
            for seg2 = 1:m2
                h2 = kn.k2(seg2+1) - kn.k2(seg2);
                [Cloc, dC] = bezier_patch_derivative(Yc, seg1, seg2, d1, d2, h1, h2);
                for i = 0:d1
                    for j = 0:d2
                        Yt  = Cloc{i+1,j+1};
                        dYt = dC{i+1,j+1};
                        M = [ -dYt + A*Yt + Yt*A',  mdl.Be,          Yt*mdl.C' ;
                               mdl.Be',            -gam2*eye(nd),    mdl.De'   ;
                               mdl.C*Yt,            mdl.De,         -eye(ne)  ];
                        Fc{end+1} = (M <= -cfg.tol*Ifl);
                    end
                end
            end
        end
    end

    Ijp = eye(2*n);
    for p = kn.i1
        for q = 1:kn.N2
            M = mdl.A_J1*Yc{p,q} + mdl.B_J1*L1c{p,q};
            Fc{end+1} = ([-Yc{p,q}, M'; M, -Yc{1,q}] <= -cfg.tol*Ijp);
        end
    end
    for q = kn.i2
        for p = 1:kn.N1
            M = mdl.A_J2*Yc{p,q} + mdl.B_J2*L2c{p,q};
            Fc{end+1} = ([-Yc{p,q}, M'; M, -Yc{p,1}] <= -cfg.tol*Ijp);
        end
    end

    gainpen = 0;
    if ~fixedGains
        for p = kn.i1
            for q = 1:kn.N2
                Fc{end+1} = (L1c{p,q} >= -cfg.Lmax);
                Fc{end+1} = (L1c{p,q} <=  cfg.Lmax);
                gainpen = gainpen + sum(L1c{p,q}.^2);
            end
        end
        for q = kn.i2
            for p = 1:kn.N1
                Fc{end+1} = (L2c{p,q} >= -cfg.Lmax);
                Fc{end+1} = (L2c{p,q} <=  cfg.Lmax);
                gainpen = gainpen + sum(L2c{p,q}.^2);
            end
        end
    end

    F = [Fc{:}];
    % A gain penalty is used ONLY when gamma is prescribed. When gamma is the
    % reported quantity the penalty would bias it upward (it reaches O(1)-O(100)
    % for these control-point counts), so the objective is gamma^2 alone.
    if minimizing
        obj = gam2;
    else
        obj = cfg.regL * gainpen;
    end

    try
        sol = optimize(F, obj, cfg.options);
    catch ME
        warning('solve_dual:solverError','%s', ME.message);
        dg(3) = -99;  return
    end

    dg(3) = sol.problem;
    if sol.problem ~= 0, return; end
    dg(1) = min(check(F));
    if dg(1) <= cfg.resTol, return; end

    if minimizing, gamma = sqrt(max(value(gam2),0)); else, gamma = sqrt(gam2fix); end
    Yval  = cellfun(@value, Yc, 'UniformOutput', false);
    L1val = cellfun(@(M) val_or_zero(M,n), L1c, 'UniformOutput', false);
    L2val = cellfun(@(M) val_or_zero(M,n), L2c, 'UniformOutput', false);
    [dg(2), dg(4)] = eig_range(Yval);
end

function V = val_or_zero(M, n)
    if isnumeric(M), V = M; else, V = value(M); end
    if isempty(V), V = zeros(1,n); end
end

function [lo, hi] = eig_range(Cval)
    lo = inf; hi = -inf;
    for p = 1:size(Cval,1)
        for q = 1:size(Cval,2)
            e = eig(Cval{p,q});
            lo = min(lo, min(e));  hi = max(hi, max(e));
        end
    end
end

%% ========================================================================
%  Local functions -- Bezier machinery
%  ========================================================================

function Ce = bezier_elevate1D(C, n)
% Exact Bezier degree elevation, degree n -> degree n+1.
    Ce = cell(1, n+2);
    for k = 0:n+1
        term = 0;
        if k >= 1, term = term + (k/(n+1))     * C{k};   end
        if k <= n, term = term + (1 - k/(n+1)) * C{k+1}; end
        Ce{k+1} = term;
    end
end

function [Cloc, dC] = bezier_patch_derivative(Cctrl, seg1, seg2, d1, d2, h1, h2)
% Control points of patch (seg1,seg2) and the diagonal directional derivative
% dC/dtau1 + dC/dtau2, exactly degree-elevated back to (d1,d2) so it lives in
% the same Bernstein basis as the patch. h1,h2 are this patch's own widths.
    base1 = (seg1-1)*d1;   base2 = (seg2-1)*d2;

    Cloc = cell(d1+1, d2+1);
    for i = 0:d1
        for j = 0:d2
            Cloc{i+1,j+1} = Cctrl{base1+i+1, base2+j+1};
        end
    end

    Q1 = cell(d1, d2+1);
    for i = 0:d1-1
        for j = 0:d2
            Q1{i+1,j+1} = d1*(Cloc{i+2,j+1} - Cloc{i+1,j+1});
        end
    end
    Q1e = cell(d1+1, d2+1);
    for j = 0:d2
        col = cell(1,d1);
        for i = 0:d1-1, col{i+1} = Q1{i+1,j+1}; end
        colE = bezier_elevate1D(col, d1-1);
        for i = 0:d1, Q1e{i+1,j+1} = colE{i+1}; end
    end

    Q2 = cell(d1+1, d2);
    for i = 0:d1
        for j = 0:d2-1
            Q2{i+1,j+1} = d2*(Cloc{i+1,j+2} - Cloc{i+1,j+1});
        end
    end
    Q2e = cell(d1+1, d2+1);
    for i = 0:d1
        row = cell(1,d2);
        for j = 0:d2-1, row{j+1} = Q2{i+1,j+1}; end
        rowE = bezier_elevate1D(row, d2-1);
        for j = 0:d2, Q2e{i+1,j+1} = rowE{j+1}; end
    end

    dC = cell(d1+1, d2+1);
    for i = 0:d1
        for j = 0:d2
            dC{i+1,j+1} = Q1e{i+1,j+1}/h1 + Q2e{i+1,j+1}/h2;
        end
    end
end

function [seg, s] = locate_segment(t, knots)
    m = numel(knots) - 1;
    t = min(max(t, knots(1)), knots(end));
    seg = m;
    for k = 1:m
        if t <= knots(k+1) + 1e-12, seg = k; break; end
    end
    s = min(max((t - knots(seg))/(knots(seg+1) - knots(seg)), 0), 1);
end

function val = eval_bezier2D(Cval, t1, t2, knots1, knots2, d1, d2)
    [seg1, s1] = locate_segment(t1, knots1);
    [seg2, s2] = locate_segment(t2, knots2);
    base1 = (seg1-1)*d1;   base2 = (seg2-1)*d2;
    val = 0;
    for i = 0:d1
        bi = nchoosek(d1,i) * s1^i * (1-s1)^(d1-i);
        for j = 0:d2
            bj = nchoosek(d2,j) * s2^j * (1-s2)^(d2-j);
            val = val + bi*bj*Cval{base1+i+1, base2+j+1};
        end
    end
end

function pts = segment_interior_grid(knots, nPerSeg)
% Points strictly inside each segment, so finite differences never straddle a
% knot (where the composite surface is only C^0).
    pts = [];
    for k = 1:numel(knots)-1
        e = linspace(knots(k), knots(k+1), nPerSeg + 2);
        pts = [pts, e(2:end-1)]; %#ok<AGROW>
    end
end

%% ========================================================================
%  Local functions -- verification, simulation, reporting
%  ========================================================================

function verify_certificate(name, Xat, AK1at, AK2at, cfg, mdl, kn, gam2)
% Dense-grid re-check in the ORIGINAL coordinates X. Xat, AK1at, AK2at are
% function handles, so the same routine covers the primal certificate, the
% dual fixed-gain certificate (X = Y^{-1}) and the designed clock-scheduled
% controller (AKi depends on the clocks).

    nPerSeg = 6;
    g1 = segment_interior_grid(kn.k1, nPerSeg);
    g2 = segment_interior_grid(kn.k2, nPerSeg);
    e  = 1e-6 * min([diff(kn.k1), diff(kn.k2)]);

    n = mdl.n;  nd = mdl.nd;
    wPos = inf; wFlow = -inf; wJump = -inf;

    for t1 = g1
        for t2 = g2
            X  = Xat(t1,t2);
            dX = (Xat(t1+e,t2+e) - Xat(t1-e,t2-e)) / (2*e);
            wPos = min(wPos, min(eig((X+X')/2)));

            for v = 1:numel(mdl.Averts)
                A = mdl.Averts{v};
                M = [ dX + A'*X + X*A + mdl.C'*mdl.C, ...
                            X*mdl.Be + mdl.C'*mdl.De ;
                      mdl.Be'*X + mdl.De'*mdl.C, ...
                           -gam2*eye(nd) + mdl.De'*mdl.De ];
                wFlow = max(wFlow, max(eig((M+M')/2)));
            end

            if t1 >= cfg.T1min
                AKq = AK1at(t1,t2);
                M = AKq'*Xat(0,t2)*AKq - X;
                wJump = max(wJump, max(eig((M+M')/2)));
            end
            if t2 >= cfg.T2min
                AKq = AK2at(t1,t2);
                M = AKq'*Xat(t1,0)*AKq - X;
                wJump = max(wJump, max(eig((M+M')/2)));
            end
        end
    end

    fprintf('--- verification: %s\n', name);
    fprintf('    grid %d x %d, %d vertices\n', numel(g1), numel(g2), numel(mdl.Averts));
    fprintf('    min eig X                 : %+.3e  (want > 0)\n', wPos);
    fprintf('    max eig flow BRL          : %+.3e  (want < 0)\n', wFlow);
    fprintf('    max eig closed-loop jump  : %+.3e  (want < 0)\n', wJump);
    if wPos > 0 && wFlow < 0 && wJump < 0
        fprintf('    VERIFIED on the grid.\n\n');
    else
        warning('verify_certificate:failed','%s failed grid verification.', name);
        fprintf('\n');
    end
end

function t = draw_sampling_times(Tmin, Tmax, Tfinal)
% Inter-sample intervals drawn i.i.d. uniformly from [Tmin,Tmax]. Clocks start
% at zero and the first instant is one full interval later, so every reset
% occurs at a certified clock value.
    t = [];  tk = 0;
    while true
        tk = tk + Tmin + (Tmax - Tmin)*rand();
        if tk > Tfinal, break; end
        t(end+1) = tk; %#ok<AGROW>
    end
end

function [tvec, Z] = simulate_cl(A, Be, A_J1, A_J2, B_J1, B_J2, ...
                                 K1at, K2at, cfg, Tfinal, dfun)
    t_f = draw_sampling_times(cfg.T1min, cfg.T1max, Tfinal);
    t_s = draw_sampling_times(cfg.T2min, cfg.T2max, Tfinal);
    events = unique([t_f(:); t_s(:)]);
    events = events(events <= Tfinal);
    isFast = ismember(events, t_f(:));
    isSlow = ismember(events, t_s(:));

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
    n = size(A,1);
    z = zeros(n,1);  tvec = 0;  Z = z.';  tcur = 0;
    lastFast = 0;  lastSlow = 0;

    for k = 1:numel(events)+1
        if k <= numel(events), te = events(k); else, te = Tfinal; end
        if te > tcur
            [tt, zz] = ode45(@(t,zz) A*zz + Be*dfun(t), [tcur te], z, opts);
            tvec = [tvec; tt]; Z = [Z; zz]; %#ok<AGROW>
            z = zz(end,:).';  tcur = te;
        end
        if k <= numel(events)
            tau1 = te - lastFast;  tau2 = te - lastSlow;
            if isFast(k)
                z = A_J1*z + B_J1*(K1at(tau1,tau2)*z);  lastFast = te;
            end
            if isSlow(k)
                z = A_J2*z + B_J2*(K2at(tau1,tau2)*z);  lastSlow = te;
            end
            tvec = [tvec; te]; Z = [Z; z.']; %#ok<AGROW>
        end
    end
end

function worst = empirical_gain(cfg, mdl, K1at, K2at, Tfinal, freqs)
% ||e||_L2 / ||d||_L2 over a frequency grid and both uncertainty vertices.
% Energy integrals are carried as extra ODE states. LOWER bound on the gain.
    worst = 0;
    n = mdl.n;
    for v = 1:numel(mdl.Averts)
        A = mdl.Averts{v};
        for om = freqs
            dfun = @(t) (t <= 2) .* sin(om*t);
            t_f = draw_sampling_times(cfg.T1min, cfg.T1max, Tfinal);
            t_s = draw_sampling_times(cfg.T2min, cfg.T2max, Tfinal);
            events = unique([t_f(:); t_s(:)]);
            events = events(events <= Tfinal);
            isFast = ismember(events, t_f(:));
            isSlow = ismember(events, t_s(:));

            opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
            odef = @(t,s) [ A*s(1:n) + mdl.Be*dfun(t);
                            norm(mdl.C*s(1:n) + mdl.De*dfun(t))^2;
                            norm(dfun(t))^2 ];
            xi = [zeros(n,1); 0; 0];  tcur = 0;
            lastFast = 0;  lastSlow = 0;

            for k = 1:numel(events)+1
                if k <= numel(events), te = events(k); else, te = Tfinal; end
                if te > tcur
                    [~, ss] = ode45(odef, [tcur te], xi, opts);
                    xi = ss(end,:).';  tcur = te;
                end
                if k <= numel(events)
                    tau1 = te - lastFast;  tau2 = te - lastSlow;
                    z = xi(1:n);
                    if isFast(k)
                        z = mdl.A_J1*z + mdl.B_J1*(K1at(tau1,tau2)*z);  lastFast = te;
                    end
                    if isSlow(k)
                        z = mdl.A_J2*z + mdl.B_J2*(K2at(tau1,tau2)*z);  lastSlow = te;
                    end
                    xi(1:n) = z;
                end
            end
            worst = max(worst, sqrt(xi(n+1)/max(xi(n+2), eps)));
        end
    end
end

function report_solve(gamma, dg, cfg)
    if isnan(gamma)
        if dg(3) == -99
            fprintf('    solver exception.\n\n');
        elseif isnan(dg(1))
            fprintf('    NOT certified (solver problem code %g).\n\n', dg(3));
        else
            fprintf(['    NOT certified: solved but residual %.3e <= resTol %.1e.\n' ...
                     '    Loosen cfg.resTol before reading this as infeasible.\n\n'], ...
                    dg(1), cfg.resTol);
        end
        return
    end
    fprintf('    gamma = %.4f    worst residual %.3e\n', gamma, dg(1));
    fprintf('    eig range of the certificate: [%.3e, %.3e]\n', dg(2), dg(4));
    if dg(2) < 10*min(cfg.Ymin, cfg.Xmin)
        fprintf('    WARNING: lower guard nearly active; loosen it.\n');
    end
    if dg(4) > 0.1*max(cfg.Ymax, cfg.Xmax)
        fprintf('    WARNING: upper guard nearly active; loosen it.\n');
    end
    fprintf('\n');
end

function s = n2s(g)
    if isnan(g), s = '   not certified'; else, s = sprintf('%15.4f', g); end
end

function s = tf2str(b)
    if b, s = 'OK  '; else, s = 'FAIL'; end
end

%% MULTIRATE_HINF_RATE_TRADEOFF
%
% Sweeps the robust multirate L2-gain design over the sampling-rate plane
% and produces the rate-allocation figures:
%
%   Fig 1 - iso-gamma contours of the certified L2 gain over (T1_max,T2_max),
%           the uncertified region, a family of sampling-budget curves
%           1/T1_max + 1/T2_max = nu, the gamma-minimizing point on each,
%           and the locus of those optima (allocation frontier).
%   Fig 2 - gamma* along each budget curve.
%
% The design solved at each point is the one in multirate_hinf_design.m:
% robust state-feedback synthesis with clock-scheduled reset gains,
% flow-only performance channels, and a polytopic uncertainty beta in the
% slow actuator.
%
% NOTATION
%   beta = uncertain slow-actuator effectiveness (the plant parameter)
%   nu   = sampling budget 1/T1_max + 1/T2_max  [Hz]
%   These are deliberately different symbols; do not conflate them.
%
% -------------------------------------------------------------------------
% RANGES
%   Both axes run up to 1 s, so the box includes the region where the
%   nominally "fast" loop is slower than the nominally "slow" one. The
%   T1 = T2 diagonal is drawn for reference. Grids are LOGARITHMIC: the
%   interesting behaviour is concentrated at small T.
%
% APERIODICITY CONVENTION
%   T_min_i = kappa * T_max_i with kappa fixed, so relative sampling jitter
%   is constant across the sweep and changes in gamma are attributable to
%   the sampling RATE, not to a changing degree of aperiodicity.
%   kappa = 0.6 reproduces the baseline [0.03,0.05], [0.30,0.50].
%
% ADAPTIVE KNOT REFINEMENT
%   Spline resolution needed on axis i scales with ||A||*T_i, not with T_i
%   alone. Segment counts are chosen per axis and per grid point so every
%   segment satisfies ||A||*h <= alpha. Without this, points would be
%   reported "not certified" because the RELAXATION ran out of resolution
%   rather than because no controller exists.
%
% BUDGET CURVES
%   nu = 1/T1_max + 1/T2_max is the total guaranteed sampling rate (loop i
%   is sampled at least once every T_max_i). LARGE nu = generous budget.
%   The baseline (0.05, 0.50) sits at nu = 22 Hz.
%
%   Restricting to the equality is without loss of generality provided
%   gamma* decreases monotonically toward small T, in which case the optimum
%   over 1/T1 + 1/T2 <= nu lies on the boundary. Because T_min_i =
%   kappa*T_max_i shifts rather than nests the admissible sampling sets,
%   this monotonicity is not automatic -- read it off the contour map rather
%   than assuming it.
%
% NUMERICAL CONVENTIONS
%   The strictness margin cfg.tol and the residual acceptance threshold
%   cfg.resTol satisfy |resTol| < tol, so an ACCEPTED solution satisfies
%   M <= -(tol-|resTol|)*I, i.e. genuine strict feasibility. The count of
%   points rejected by the residual test alone is reported: if it is
%   non-zero, loosen cfg.resTol before believing a hole in the contour map.
%
%   No regularization is placed on L. The objective is gamma^2 alone. A
%   penalty rho*sum(L.^2) with rho = 1e-4 contributes O(1)-O(100) to the
%   objective for the control-point counts used here and would bias the
%   reported gamma upward; L is bounded by its box instead.
%
% Requires: YALMIP and an SDP solver (MOSEK by default).

clear; close all; clc;

%% ------------------------------------------------------------------------
%  Sweep configuration
%  ------------------------------------------------------------------------

cfg.kappa = 0.6;                          % T_min_i = kappa * T_max_i

n1 = 15;  n2 = 15;                        % grid resolution
T1_grid = logspace(log10(0.030), log10(1.0), n1);
T2_grid = logspace(log10(0.030), log10(1.0), n2);

nu_list = [22 14 9 6 4];                  % sampling budgets [Hz], generous -> tight
n_curve = 15;                              % solves per budget curve

T1_base = 0.05;   T2_base = 0.50;         % baseline design point (nu = 22 Hz)
gamma_cap = 5e2;                          % larger values -> "not certified"

USE_PARFOR = false;
SAVE_FILE  = 'rate_tradeoff_results.mat';

%% ------------------------------------------------------------------------
%  Plant, uncertainty, performance weights
%  ------------------------------------------------------------------------

cfg.Ap    = [ 2.0  0.50;
             -0.5  0.25];
cfg.B1    = [1; 0];                       % fast loop actuates x1
cfg.B2fun = @(b) [0; b];                  % slow loop, uncertain effectiveness
cfg.Ep    = [0; 1];                       % disturbance enters x2

cfg.beta_min = 0.6;                       % uncertain actuator effectiveness
cfg.beta_max = 1.0;

cfg.w1  = 1.0;
cfg.w2  = 1.0;
cfg.rho = 0.1;

%% ------------------------------------------------------------------------
%  Relaxation, numerics, solver
%  ------------------------------------------------------------------------

cfg.d1 = 2;   cfg.d2 = 2;                 % Bezier degree per axis

cfg.alpha = 0.40;                         % require ||A|| * h <= alpha
cfg.s_max = 5;                            % cap on segments per region

cfg.tol    = 1e-5;                        % LMI strictness margin
cfg.resTol = -1e-6;                       % residual acceptance (|resTol| < tol)
cfg.Ymin   = 1e-6;                        % Y >= Ymin*I   (positivity only)
cfg.Ymax   = 1e6;                         % Y <= Ymax*I   (numerical guard)
cfg.Lmax   = 5e2;                         % box on L      (no penalty on L)

Anorm = max(norm(augA(cfg, cfg.beta_min)), norm(augA(cfg, cfg.beta_max)));
cfg.hmax = cfg.alpha / Anorm;
fprintf('||A|| = %.4f  ->  max segment width h = %.5f s\n', Anorm, cfg.hmax);

cfg.options = sdpsettings('solver','mosek','verbose',0);
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-9;
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-9;
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-9;

% diagnostics column layout returned by solve_gamma
%   1 s1tot | 2 s2tot | 3 capped | 4 resid | 5 solver problem | 6 min eig Y
%   7 max eig Y
DCOL = struct('s1',1,'s2',2,'capped',3,'resid',4,'prob',5,'ymin',6,'ymax',7);

%% ------------------------------------------------------------------------
%  Grid sweep
%  ------------------------------------------------------------------------

nTot = n1*n2;
fprintf('Grid sweep: %d x %d = %d designs\n', n1, n2, nTot);

Gvec = nan(nTot,1);
Dvec = nan(nTot,7);
tSweep = tic;

if USE_PARFOR
    parfor idx = 1:nTot
        [i1, i2] = ind2sub([n1 n2], idx);
        [g, dg] = solve_gamma(T1_grid(i1), T2_grid(i2), cfg);
        Gvec(idx) = g;  Dvec(idx,:) = dg;
    end
else
    for idx = 1:nTot
        [i1, i2] = ind2sub([n1 n2], idx);
        [g, dg] = solve_gamma(T1_grid(i1), T2_grid(i2), cfg);
        Gvec(idx) = g;  Dvec(idx,:) = dg;
        if mod(idx, max(1,round(nTot/20))) == 0 || idx == nTot
            el = toc(tSweep);
            fprintf('  %3d/%3d   elapsed %6.1fs   eta %6.1fs\n', ...
                    idx, nTot, el, el*(nTot-idx)/idx);
        end
    end
end

% G indexed (T2, T1) so contour(T1_grid, T2_grid, G) works directly.
G       = reshape(Gvec, [n1 n2]).';
Capped  = reshape(Dvec(:,DCOL.capped), [n1 n2]).' > 0;
Resid   = reshape(Dvec(:,DCOL.resid ), [n1 n2]).';
Prob    = reshape(Dvec(:,DCOL.prob  ), [n1 n2]).';
G(G > gamma_cap) = NaN;

report_diagnostics(G, Capped, Resid, Prob, Dvec, DCOL, cfg, toc(tSweep));

%% ------------------------------------------------------------------------
%  Budget curves
%  ------------------------------------------------------------------------

curves = struct('nu',{},'T1',{},'T2',{},'g',{}, ...
                'T1opt',{},'T2opt',{},'gopt',{},'atEdge',{}, ...
                'T1full',{},'T2full',{});

for b = 1:numel(nu_list)
    nu = nu_list(b);

    lo = max(min(T1_grid), 1.001/nu);
    hi = max(T1_grid);
    if lo >= hi
        fprintf('Budget nu = %g Hz: asymptote outside the box, skipped.\n', nu);
        continue
    end

    T1f = logspace(log10(lo), log10(hi), 400);
    T2f = 1 ./ (nu - 1./T1f);
    ok  = isfinite(T2f) & T2f > 0 & ...
          T2f >= min(T2_grid) & T2f <= max(T2_grid) & ...
          T1f >= min(T1_grid) & T1f <= max(T1_grid);
    T1f = T1f(ok);  T2f = T2f(ok);

    if numel(T1f) < 2
        fprintf('Budget nu = %g Hz: curve outside the box, skipped.\n', nu);
        continue
    end

    sel = unique(round(linspace(1, numel(T1f), n_curve)));
    T1c = T1f(sel);  T2c = T2f(sel);

    fprintf('Budget nu = %g Hz: %d designs\n', nu, numel(T1c));
    gc = nan(size(T1c));
    for k = 1:numel(T1c)
        gc(k) = solve_gamma(T1c(k), T2c(k), cfg);
    end
    gc(gc > gamma_cap) = NaN;

    [gopt, kopt] = min(gc);
    atEdge = false;
    if isempty(gopt) || isnan(gopt)
        T1o = NaN; T2o = NaN; gopt = NaN;
        fprintf('  no certified point on this curve\n');
    else
        T1o = T1c(kopt);  T2o = T2c(kopt);
        atEdge = (kopt == 1) || (kopt == numel(T1c));
        fprintf('  optimum: T1 = %.4f (%.2f Hz), T2 = %.4f (%.2f Hz), gamma = %.4f\n', ...
                T1o, 1/T1o, T2o, 1/T2o, gopt);
        if atEdge
            fprintf(['  WARNING: optimum sits at a curve ENDPOINT -- the curve is\n' ...
                     '           clipped by the plotted box, so this is not a\n' ...
                     '           genuine interior optimum. Widen the grid.\n']);
        end
    end

    curves(end+1) = struct('nu',nu,'T1',T1c,'T2',T2c,'g',gc, ...
                           'T1opt',T1o,'T2opt',T2o,'gopt',gopt, ...
                           'atEdge',atEdge,'T1full',T1f,'T2full',T2f); %#ok<SAGROW>
end

g_base = solve_gamma(T1_base, T2_base, cfg);
fprintf('\nBaseline (%.3f, %.3f): gamma = %s   [nu = %.2f Hz]\n\n', ...
        T1_base, T2_base, num2str(g_base,'%.4f'), 1/T1_base + 1/T2_base);

save(SAVE_FILE, 'T1_grid','T2_grid','G','Capped','Resid','Prob','Dvec', ...
     'curves','T1_base','T2_base','g_base','nu_list','cfg');

%% ------------------------------------------------------------------------
%  Figure 1: contour map with the budget-curve family
%  ------------------------------------------------------------------------

figure('Name','Rate allocation','Color','w','Position',[80 80 780 600]);
ax = axes; hold(ax,'on');
set(ax,'Color',[0.86 0.86 0.86]);          % shows through where G is NaN

gFin = G(~isnan(G));
if isempty(gFin)
    error('Nothing certified anywhere; widen the grid or relax the setup.');
end
if numel(unique(gFin)) < 2
    lev = unique(gFin) * [0.9 1.1];        % degenerate case guard
else
    lev = logspace(log10(min(gFin)), log10(max(gFin)), 14);
end

contourf(T1_grid, T2_grid, G, lev, 'LineColor','none');
set(ax,'ColorScale','log','XScale','log','YScale','log');
colormap(ax, parula);
cb = colorbar; title(cb,'$\gamma$','Interpreter','latex','FontSize',14);

[cc, hc] = contour(T1_grid, T2_grid, G, lev, 'k-', 'LineWidth', 0.7);
clabel(cc, hc, 'FontSize', 9, 'LabelSpacing', 300, 'Color','k');

% mark uncertified points whose segment cap bound (inconclusive, not a boundary)
[TT1, TT2] = meshgrid(T1_grid, T2_grid);
susp = isnan(G) & Capped;
if any(susp(:))
    plot(TT1(susp), TT2(susp), 'x', 'Color',[0.8 0 0], 'MarkerSize',8,'LineWidth',1.5);
end

% T1 = T2 diagonal: below it the "fast" loop is actually the slower one
%dg = [max(min(T1_grid),min(T2_grid)), min(max(T1_grid),max(T2_grid))];
%plot(dg, dg, ':', 'Color',[0.25 0.25 0.25], 'LineWidth', 1.6);

cmapB = lines(max(numel(curves),1));
hLeg = gobjects(0); sLeg = {};
for b = 1:numel(curves)
    cv = curves(b);
    plot(cv.T1full, cv.T2full, 'w-', 'LineWidth', 3.5);
    h = plot(cv.T1full, cv.T2full, '--', 'Color', cmapB(b,:), 'LineWidth', 2);
    hLeg(end+1) = h;                                          %#ok<SAGROW>
    sLeg{end+1} = sprintf('$\\nu = %g$', cv.nu);           %#ok<SAGROW>
    if ~isnan(cv.T1opt)
        plot(cv.T1opt, cv.T2opt, 'p', 'MarkerSize', 17, ...
             'MarkerFaceColor', cmapB(b,:), 'MarkerEdgeColor','k','LineWidth',1);
    end
end

% allocation frontier: locus of optima as the budget tightens
if ~isempty(curves)
    fT1 = [curves.T1opt];  fT2 = [curves.T2opt];
    ok  = ~isnan(fT1) & ~[curves.atEdge];
    if nnz(ok) >= 2
        %plot(fT1(ok), fT2(ok), 'k-', 'LineWidth', 2.4);
    end
end

plot(T1_base, T2_base, 'o', 'MarkerSize', 10, 'MarkerFaceColor','w', ...
     'MarkerEdgeColor','k','LineWidth',1.5);
text(T1_base+0.004, T2_base+0.04, '  baseline', 'Interpreter','latex', ...
     'FontSize',12,'VerticalAlignment','top');

xlabel('sampling time bound $\bar{T}_1$','Interpreter','latex','FontSize',18);
ylabel('sampling time bound $\bar{T}_2$','Interpreter','latex','FontSize',18);
if ~isempty(hLeg)
    legend(hLeg, sLeg, 'Interpreter','latex','Location','southwest','FontSize',11);
end
set(ax,'FontSize',13,'Layer','top','Box','on');
axis([min(T1_grid) max(T1_grid) min(T2_grid) max(T2_grid)]);
hold(ax,'off');

%% ------------------------------------------------------------------------
%  Figure 2: gamma along each budget curve
%  ------------------------------------------------------------------------

figure('Name','Gain along budget curves','Color','w','Position',[900 80 640 460]);
hold on
hL = gobjects(0); sL = {};
for b = 1:numel(curves)
    cv = curves(b);
    h = plot(1./cv.T1, cv.g, '-o', 'Color', cmapB(b,:), 'LineWidth', 2, ...
             'MarkerFaceColor','w','MarkerSize',5);
    hL(end+1) = h;                                            %#ok<SAGROW>
    sL{end+1} = sprintf('$\\nu = %g$ Hz', cv.nu);             %#ok<SAGROW>
    if ~isnan(cv.gopt)
        plot(1/cv.T1opt, cv.gopt, 'p', 'MarkerSize', 16, ...
             'MarkerFaceColor', cmapB(b,:), 'MarkerEdgeColor','k','LineWidth',1);
    end
end
grid on; set(gca,'XScale','log','YScale','log','FontSize',13);
xlabel('sampling rate bound $1/\bar{T}_1$ [Hz]','Interpreter','latex','FontSize',18);
ylabel('$\gamma^\star$','Interpreter','latex','FontSize',18);
if ~isempty(hL)
    legend(hL, sL, 'Interpreter','latex','Location','best','FontSize',11);
end
hold off

save('rate_tradeoff_results.mat')


%% ========================================================================
%  Local functions
%  ========================================================================

function A = augA(cfg, beta)
% Augmented flow matrix z = [x; u1_hold; u2_hold] for a given effectiveness.
    np = size(cfg.Ap,1);
    A  = [cfg.Ap, cfg.B1, cfg.B2fun(beta);
          zeros(2, np+2)];
end

function [s_pre, s_post, capped] = choose_segments(Tmax, cfg)
% Adaptive segment counts. The pre-region [0, kappa*Tmax] and post-region
% [kappa*Tmax, Tmax] are refined independently so every segment width
% satisfies ||A||*h <= alpha. Counts are capped at cfg.s_max; "capped" flags
% that the requested resolution was not delivered, in which case an
% uncertified result is INCONCLUSIVE rather than a feasibility boundary.

    need_pre  = cfg.kappa     * Tmax / cfg.hmax;
    need_post = (1-cfg.kappa) * Tmax / cfg.hmax;

    s_pre  = min(cfg.s_max, max(1, ceil(need_pre)));
    s_post = min(cfg.s_max, max(1, ceil(need_post)));

    capped = (need_pre > cfg.s_max) || (need_post > cfg.s_max);
end

function [gamma, dg] = solve_gamma(T1max, T2max, cfg)
% Robust multirate L2-gain design for one pair of upper dwell-time bounds.
%   gamma : minimized certified gain, NaN if nothing is certified
%   dg    : [s1tot s2tot capped resid problem minEigY maxEigY]

    yalmip('clear');                       % keep the model table from growing

    gamma = NaN;
    dg    = nan(1,7);

    T1min = cfg.kappa * T1max;
    T2min = cfg.kappa * T2max;

    [s1_pre, s1_post, cap1] = choose_segments(T1max, cfg);
    [s2_pre, s2_post, cap2] = choose_segments(T2max, cfg);
    dg(1) = s1_pre + s1_post;
    dg(2) = s2_pre + s2_post;
    dg(3) = double(cap1 || cap2);

    % ---- augmented model ------------------------------------------------
    n      = size(cfg.Ap,1) + 2;
    Averts = {augA(cfg, cfg.beta_min), augA(cfg, cfg.beta_max)};

    Be = [cfg.Ep; 0; 0];
    C  = diag([cfg.w1 cfg.w2 cfg.rho cfg.rho]);
    De = zeros(size(C,1), size(Be,2));
    nd = size(Be,2);
    ne = size(C,1);

    B_J1 = [0;0;1;0];   A_J1 = diag([1 1 0 1]);
    B_J2 = [0;0;0;1];   A_J2 = diag([1 1 1 0]);

    % ---- knots (exact breakpoint at T_min on each axis) -----------------
    d1 = cfg.d1;  d2 = cfg.d2;
    [knots1, idxJump1] = build_knots(T1min, T1max, s1_pre, s1_post, d1);
    [knots2, idxJump2] = build_knots(T2min, T2max, s2_pre, s2_post, d2);

    m1 = numel(knots1) - 1;   N1 = m1*d1 + 1;
    m2 = numel(knots2) - 1;   N2 = m2*d2 + 1;

    % ---- variables ------------------------------------------------------
    Yctrl  = cell(N1,N2);
    L1ctrl = cell(N1,N2);
    L2ctrl = cell(N1,N2);
    for p = 1:N1
        for q = 1:N2
            Yctrl{p,q}  = sdpvar(n, n, 'symmetric');
            L1ctrl{p,q} = zeros(1,n);
            L2ctrl{p,q} = zeros(1,n);
        end
    end
    % L_Ji is only evaluated for tau_i >= T_min_i, whose patches use exactly
    % the control points in idxJump_i; the rest stay at zero.
    for p = idxJump1
        for q = 1:N2, L1ctrl{p,q} = sdpvar(1, n, 'full'); end
    end
    for q = idxJump2
        for p = 1:N1, L2ctrl{p,q} = sdpvar(1, n, 'full'); end
    end

    gam2 = sdpvar(1,1);

    % Constraints accumulate in a cell and are concatenated once: repeated
    % F = [F, c] is quadratic in the number of constraints.
    Fc = {};
    Fc{end+1} = (gam2 >= 0);

    % ---- positivity / numerical guards on Y -----------------------------
    for p = 1:N1
        for q = 1:N2
            Fc{end+1} = (Yctrl{p,q} >= cfg.Ymin*eye(n));
            Fc{end+1} = (Yctrl{p,q} <= cfg.Ymax*eye(n));
        end
    end

    % ---- flow LMIs, one family per uncertainty vertex --------------------
    %   [ -Ydot + A Y + Y A'   Be        Y C' ]
    %   [        *            -g^2 I     De'  ]  < 0
    %   [        *              *        -I   ]
    Ifl = eye(n + nd + ne);
    for v = 1:numel(Averts)
        A = Averts{v};
        for seg1 = 1:m1
            h1 = knots1(seg1+1) - knots1(seg1);
            for seg2 = 1:m2
                h2 = knots2(seg2+1) - knots2(seg2);
                [Cloc, dC] = bezier_patch_derivative(Yctrl, seg1, seg2, ...
                                                     d1, d2, h1, h2);
                for i = 0:d1
                    for j = 0:d2
                        Yt  = Cloc{i+1,j+1};
                        dYt = dC{i+1,j+1};
                        M = [ -dYt + A*Yt + Yt*A',  Be,             Yt*C'    ;
                               Be',                -gam2*eye(nd),   De'      ;
                               C*Yt,                De,            -eye(ne) ];
                        Fc{end+1} = (M <= -cfg.tol*Ifl);
                    end
                end
            end
        end
    end

    % ---- jump LMIs -------------------------------------------------------
    %   [ -Y                    (A_Ji Y + B_Ji L_Ji)' ]
    %   [  A_Ji Y + B_Ji L_Ji   -Y^{0,i}              ]  < 0
    Ijp = eye(2*n);
    for p = idxJump1
        for q = 1:N2
            M = A_J1*Yctrl{p,q} + B_J1*L1ctrl{p,q};
            Fc{end+1} = ([-Yctrl{p,q}, M'; M, -Yctrl{1,q}] <= -cfg.tol*Ijp);
        end
    end
    for q = idxJump2
        for p = 1:N1
            M = A_J2*Yctrl{p,q} + B_J2*L2ctrl{p,q};
            Fc{end+1} = ([-Yctrl{p,q}, M'; M, -Yctrl{p,1}] <= -cfg.tol*Ijp);
        end
    end

    % ---- box on L (no penalty: a penalty would bias gamma) ---------------
    for p = idxJump1
        for q = 1:N2
            Fc{end+1} = (L1ctrl{p,q} >= -cfg.Lmax);
            Fc{end+1} = (L1ctrl{p,q} <=  cfg.Lmax);
        end
    end
    for q = idxJump2
        for p = 1:N1
            Fc{end+1} = (L2ctrl{p,q} >= -cfg.Lmax);
            Fc{end+1} = (L2ctrl{p,q} <=  cfg.Lmax);
        end
    end

    F = [Fc{:}];

    % ---- solve -----------------------------------------------------------
    try
        sol = optimize(F, gam2, cfg.options);
    catch ME
        warning('solve_gamma:solverError', ...
                'Solver failed at (%.4f, %.4f): %s', T1max, T2max, ME.message);
        dg(5) = -99;
        return
    end

    dg(5) = sol.problem;
    if sol.problem ~= 0
        return
    end

    resid = min(check(F));
    dg(4) = resid;
    if resid <= cfg.resTol
        return                              % solved, but residuals not trusted
    end

    gamma = sqrt(max(value(gam2), 0));

    % eigenvalue range of Y over the control points, to expose whether the
    % Ymin/Ymax guards are anywhere near active
    emin = inf; emax = -inf;
    for p = 1:N1
        for q = 1:N2
            e = eig(value(Yctrl{p,q}));
            emin = min(emin, min(e));  emax = max(emax, max(e));
        end
    end
    dg(6) = emin;  dg(7) = emax;
end

function [knots, idxJump] = build_knots(Tmin, Tmax, s_pre, s_post, d)
% Knot vector on [0,Tmax] with an EXACT breakpoint at Tmin, so the
% jump-constrained control points form a union of whole Bezier patches.
% This is mandatory: selecting control points by a numerical test such as
% "grid value >= Tmin" leaves part of the admissible clock range governed by
% unconstrained control points, producing certificates that are not valid.
    assert(Tmin > 0 && Tmax > Tmin, 'Require 0 < Tmin < Tmax.');
    pre   = linspace(0, Tmin, s_pre + 1);
    post  = linspace(Tmin, Tmax, s_post + 1);
    knots = [pre, post(2:end)];
    idxJump = (s_pre*d + 1) : ((s_pre + s_post)*d + 1);
end

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
% Control points of patch (seg1,seg2) and the diagonal directional
% derivative dC/dtau1 + dC/dtau2, exactly degree-elevated back to (d1,d2) so
% it lives in the same Bernstein basis as the patch itself. h1,h2 are this
% patch's own segment widths (the knot vectors are non-uniform).
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

function report_diagnostics(G, Capped, Resid, Prob, Dvec, DCOL, cfg, elapsed)
% Print everything needed to tell a genuine feasibility boundary from a
% numerical artefact.
    nTot  = numel(G);
    nCert = sum(~isnan(G(:)));
    fprintf('\nGrid done in %.1f s. Certified at %d of %d points.\n', ...
            elapsed, nCert, nTot);

    nSolvedButRejected = sum(Prob(:) == 0 & isnan(G(:)) & ~isnan(Resid(:)));
    if nSolvedButRejected > 0
        fprintf(['NOTE: %d point(s) SOLVED but were rejected by the residual\n' ...
                 '      test (resTol = %.1e). Worst rejected residual = %.3e.\n' ...
                 '      Loosen cfg.resTol before reading these as infeasible.\n'], ...
                nSolvedButRejected, cfg.resTol, ...
                min(Resid(Prob == 0 & isnan(G))));
    end

    nSusp = sum(isnan(G(:)) & Capped(:));
    if nSusp > 0
        fprintf(['NOTE: %d uncertified point(s) hit the segment cap s_max = %d\n' ...
                 '      (marked with red crosses). These are INCONCLUSIVE.\n'], ...
                nSusp, cfg.s_max);
    end

    nErr = sum(Dvec(:,DCOL.prob) == -99);
    if nErr > 0
        fprintf('NOTE: %d point(s) raised a solver exception.\n', nErr);
    end

    fprintf('Segments used: axis 1 in [%d,%d], axis 2 in [%d,%d] (cap %d)\n', ...
            min(Dvec(:,DCOL.s1)), max(Dvec(:,DCOL.s1)), ...
            min(Dvec(:,DCOL.s2)), max(Dvec(:,DCOL.s2)), cfg.s_max);

    ey = Dvec(~isnan(Dvec(:,DCOL.ymin)), [DCOL.ymin DCOL.ymax]);
    if ~isempty(ey)
        fprintf('eig(Y) over certified points: min %.3e, max %.3e\n', ...
                min(ey(:,1)), max(ey(:,2)));
        if min(ey(:,1)) < 10*cfg.Ymin
            fprintf('      WARNING: Ymin guard (%.1e) is nearly active.\n', cfg.Ymin);
        end
        if max(ey(:,2)) > 0.1*cfg.Ymax
            fprintf('      WARNING: Ymax guard (%.1e) is nearly active.\n', cfg.Ymax);
        end
    end
    fprintf('\n');
end

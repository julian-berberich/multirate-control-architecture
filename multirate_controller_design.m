%% MULTIRATE_CONTROLLER_DESIGN
%
% Synthesis of clock-scheduled reset gains for a multirate sampled-data
% control loop, via the dual (inverse-Lyapunov) LMI formulation.
%
% SETUP
%   Continuous-time plant       xdot = A x + B1 u1 + B2 u2
%   with two independent, asynchronous zero-order-hold loops:
%     - "fast" loop, inter-sample times in [Tf_min, Tf_max]
%     - "slow" loop, inter-sample times in [Ts_min, Ts_max]
%   Unlike multirate_stability_analysis.m, the gains are DECISION VARIABLES.
%
% FORMULATION
%   To keep the problem linear, the gains are not parameterized directly.
%   Instead we solve for
%       Y(tau1,tau2) = X(tau1,tau2)^{-1}          (inverse Lyapunov matrix)
%       L_i(tau1,tau2), i = 1,2                   (auxiliary gain factors)
%   as tensor-product composite Bezier matrix functions, and recover
%       K_i(tau1,tau2) = L_i(tau1,tau2) Y(tau1,tau2)^{-1}.
%   The jump conditions enter as Schur complements, jointly affine in
%   (Y, L_i), so the convex-hull property of the Bernstein basis still
%   certifies them exactly on each patch.
%
% FLOW CONDITION
%   In this zero-order-hold representation the augmented flow zdot = A_full z
%   is autonomous: the inputs enter the state only through resets, so there
%   is no continuous-time control channel. The flow condition of the
%   underlying proposition is therefore instantiated with B = 0, L = 0 and
%   reduces to a Lyapunov differential inequality on Y. All design freedom
%   resides in the two jump conditions.
%
% CLOCK DEPENDENCE OF THE GAINS
%   As written, K_i is allowed to depend on the full clock vector
%   (tau1, tau2). To restrict K_i to its own clock only, set
%   RESTRICT_OWN_CLOCK = true below; the recovered gains then depend on
%   tau_i alone, which is the more implementable structure and lets the two
%   variants be compared numerically.
%
% Requires: YALMIP and an SDP solver (MOSEK by default).
%
% -------------------------------------------------------------------------
% KNOT PLACEMENT (important)
%
%   The jump LMIs must hold only for tau_i in [T_min_i, T_max_i]. The
%   convex-hull argument certifies an LMI on a patch only if the LMI is
%   imposed at EVERY control point of that patch. Selecting control points
%   by a numerical test such as "grid value >= T_min" therefore yields an
%   exact restriction to [T_min, T_max] only when T_min coincides with a
%   segment boundary; otherwise part of the admissible clock range is
%   governed by control points that carry no constraint, and the resulting
%   design is not certified (it can report feasible for a closed loop that
%   diverges in simulation).
%
%   The knot vectors below therefore place an EXACT breakpoint at T_min on
%   each axis. The jump-constrained index sets are then read off directly
%   from the segment structure and are unions of whole patches by
%   construction, for any ratio T_max / T_min.
% -------------------------------------------------------------------------

clear; close all; clc;

rng(1);   % reproducible sampling sequences in the simulation

RESTRICT_OWN_CLOCK = false;   % true: K_i depends on tau_i only

%% ------------------------------------------------------------------------
%  Plant
%  ------------------------------------------------------------------------

A  = [ 2.0  0.50;
      -0.5  0.25];        % open loop unstable

B1 = [1; 0];              % actuated by the fast loop
B2 = [0; 1];              % actuated by the slow loop

n = size(A,1);

%% ------------------------------------------------------------------------
%  Dwell-time bounds
%  ------------------------------------------------------------------------

h = 0;

Tf_min = 0.03;   Tf_max = 0.05+h;     % fast loop inter-sample range
Ts_min = 0.30;   Ts_max = 0.50+h;     % slow loop inter-sample range

%% Simulation time 
Tfinal = 20;

%% ------------------------------------------------------------------------
%  Augmented impulsive representation (open loop reset maps)
%  ------------------------------------------------------------------------

A_full = [A       B1  B2;
          zeros(2,4)     ];
n_full = size(A_full,1);

B_J1 = [0;0;1;0];        % injects the refreshed fast hold
B_J2 = [0;0;0;1];        % injects the refreshed slow hold

A_J1_ol = diag([1 1 0 1]);   % fast reset, structural part
A_J2_ol = diag([1 1 1 0]);   % slow reset, structural part

%% ------------------------------------------------------------------------
%  Composite Bezier parameterization of Y, L1, L2
%  ------------------------------------------------------------------------

tol = 1e-4;      % strict-inequality margin used in all LMIs
rho = 1e-3;      % weight on gain magnitude in the objective

d1 = 2;  m1_pre = 1;  m1_post = 2;   % fast axis: degree, segments before/after Tf_min
d2 = 2;  m2_pre = 1;  m2_post = 2;   % slow axis: degree, segments before/after Ts_min

[knots1, idxJump1] = build_knots(Tf_min, Tf_max, m1_pre, m1_post, d1);
[knots2, idxJump2] = build_knots(Ts_min, Ts_max, m2_pre, m2_post, d2);

m1 = numel(knots1) - 1;   nctrl1 = m1*d1 + 1;
m2 = numel(knots2) - 1;   nctrl2 = m2*d2 + 1;

Yctrl  = cell(nctrl1, nctrl2);
L1ctrl = cell(nctrl1, nctrl2);
L2ctrl = cell(nctrl1, nctrl2);

for p = 1:nctrl1
    for q = 1:nctrl2
        Yctrl{p,q} = sdpvar(n_full, n_full, 'symmetric');
    end
end

% L1 is only ever evaluated for tau1 in [Tf_min,Tf_max], whose patches use
% exclusively the control points in idxJump1; the remaining rows are
% therefore irrelevant and are held at zero rather than parameterized.
% The same applies to L2 along axis 2.
for p = 1:nctrl1
    for q = 1:nctrl2
        L1ctrl{p,q} = zeros(1, n_full);
        L2ctrl{p,q} = zeros(1, n_full);
    end
end

for p = idxJump1
    for q = 1:nctrl2
        if RESTRICT_OWN_CLOCK && q > 1
            L1ctrl{p,q} = L1ctrl{p,1};          % no dependence on tau2
        else
            L1ctrl{p,q} = sdpvar(1, n_full, 'full');
        end
    end
end

for q = idxJump2
    for p = 1:nctrl1
        if RESTRICT_OWN_CLOCK && p > 1
            L2ctrl{p,q} = L2ctrl{1,q};          % no dependence on tau1
        else
            L2ctrl{p,q} = sdpvar(1, n_full, 'full');
        end
    end
end

I = eye(n_full);
F = [];

%% ------------------------------------------------------------------------
%  Positive definiteness, reduced to the faces of the clock box
%
%  Instead of imposing Y > 0 on the whole box, it suffices to impose it on
%  the faces {tau_i = 0}. Note the argument differs from the primal one: the
%  dual-to-primal congruence needs Y nonsingular, which is exactly what is to
%  be shown, so the primal proof cannot simply be transcribed. Working
%  directly in the dual coordinates,
%
%     d/ds [ exp(-A s) Y(tau + s*1) exp(-A' s) ]
%         = exp(-A s) [ dY/ds - A Y - Y A' ] exp(-A' s)  >  0
%
%  by the flow LMI, so that quantity is INCREASING along the diagonal.
%  Propagating BACKWARD by rho = min_i tau_i lands on a face and gives
%     Y(tau) >= exp(A rho) Y(tau - rho*1) exp(A' rho) > 0.
%  The jump LMIs are not needed for this direction.
%
%  In the tensor basis the face {tau_1 = 0} is exactly the composite Bezier
%  curve carried by the FIRST ROW of control points, and {tau_2 = 0} by the
%  FIRST COLUMN, so the convex-hull property makes this restriction exact.
%  The constraint count drops from prod_i s_i*(d_i+1) to nctrl1+nctrl2-1.
%  ------------------------------------------------------------------------

for q = 1:nctrl2
    F = [F, Yctrl{1,q} >= tol*I];          % face tau_1 = 0
end
for p = 2:nctrl1                            % face tau_2 = 0 (corner already done)
    F = [F, Yctrl{p,1} >= tol*I];
end

%% ------------------------------------------------------------------------
%  Flow LMIs
%     A_full Y + Y A_full' - dY/dtau1 - dY/dtau2 < 0
%  on [0,Tf_max] x [0,Ts_max]
%  ------------------------------------------------------------------------

for seg1 = 1:m1
    h1 = knots1(seg1+1) - knots1(seg1);
    for seg2 = 1:m2
        h2 = knots2(seg2+1) - knots2(seg2);

        [Cloc, dC] = bezier_patch_derivative(Yctrl, seg1, seg2, d1, d2, h1, h2);

        for i = 0:d1
            for j = 0:d2
                Yt  = Cloc{i+1,j+1};
                dYt = dC{i+1,j+1};

                F = [F, A_full*Yt + Yt*A_full' - dYt <= -tol*I];
            end
        end
    end
end

%% ------------------------------------------------------------------------
%  Jump LMIs (design), Schur-complement form
%
%     [ -Y                      (A_J_i Y + B_J_i L_i)' ]
%     [  A_J_i Y + B_J_i L_i    -Y^{0,i}               ] < 0
%
%  Y^{0,i} denotes Y evaluated at tau_i = 0, exactly the boundary curve
%  carried by the first row (i = 1) or first column (i = 2) of control
%  points. Since the Bernstein weights of the omitted axis sum to one,
%  imposing the inequality at every control-point pair of a patch certifies
%  it on the whole patch.
%  ------------------------------------------------------------------------

I2 = eye(2*n_full);

for p = idxJump1
    for q = 1:nctrl2
        Yt = Yctrl{p,q};
        Y0 = Yctrl{1,q};                          % Y(0, tau2)
        M  = A_J1_ol*Yt + B_J1*L1ctrl{p,q};

        F = [F, [-Yt, M'; M, -Y0] <= -tol*I2];
    end
end

for p = 1:nctrl1
    for q = idxJump2
        Yt = Yctrl{p,q};
        Y0 = Yctrl{p,1};                          % Y(tau1, 0)
        M  = A_J2_ol*Yt + B_J2*L2ctrl{p,q};

        F = [F, [-Yt, M'; M, -Y0] <= -tol*I2];
    end
end

%% ------------------------------------------------------------------------
%  Normalization, conditioning, objective
%  ------------------------------------------------------------------------

F = [F, trace(Yctrl{1,1}) == n_full];     % removes the scale invariance

obj = 0;
for p = 1:nctrl1
    for q = 1:nctrl2
        F   = [F, -100*I <= Yctrl{p,q} <= 100*I];
        obj = obj + sum(sum(Yctrl{p,q}.^2));
    end
end

for p = idxJump1
    for q = 1:nctrl2
        if isa(L1ctrl{p,q}, 'sdpvar')
            F   = [F, -500 <= L1ctrl{p,q} <= 500];
            obj = obj + rho*sum(L1ctrl{p,q}.^2);
        end
    end
end

for q = idxJump2
    for p = 1:nctrl1
        if isa(L2ctrl{p,q}, 'sdpvar')
            F   = [F, -500 <= L2ctrl{p,q} <= 500];
            obj = obj + rho*sum(L2ctrl{p,q}.^2);
        end
    end
end

%% ------------------------------------------------------------------------
%  Solve
%  ------------------------------------------------------------------------

options = sdpsettings('solver','mosek','verbose',0);
options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-10;
options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-10;

sol = optimize(F, obj, options);

fprintf('Solver status: %s\n', sol.info);

if sol.problem ~= 0
    disp('Infeasible or solver failure -- no controller obtained.');
    return
end

res = min(check(F));
fprintf('Worst LMI residual at the solution: %.3e\n', res);
if res < -1e-5
    warning('Reported feasible but residuals are violated; solve is inconclusive.');
end

Yctrl_val  = cellfun(@value, Yctrl, 'UniformOutput', false);
L1ctrl_val = cellfun(@(M) value_or_zero(M, n_full), L1ctrl, 'UniformOutput', false);
L2ctrl_val = cellfun(@(M) value_or_zero(M, n_full), L2ctrl, 'UniformOutput', false);

%% ------------------------------------------------------------------------
%  Recovered gains at the corners of their admissible clock domains
%  ------------------------------------------------------------------------

K1at = @(t1,t2) eval_bezier2D(L1ctrl_val, t1, t2, knots1, knots2, d1, d2) / ...
                eval_bezier2D(Yctrl_val,  t1, t2, knots1, knots2, d1, d2);
K2at = @(t1,t2) eval_bezier2D(L2ctrl_val, t1, t2, knots1, knots2, d1, d2) / ...
                eval_bezier2D(Yctrl_val,  t1, t2, knots1, knots2, d1, d2);

fprintf('\nFast gain K1(tau1,tau2), tau1 in [Tf_min,Tf_max]:\n');
fprintf('  at (Tf_min, 0)      : '); fprintf('%9.4f', K1at(Tf_min, 0));      fprintf('\n');
fprintf('  at (Tf_max, Ts_max) : '); fprintf('%9.4f', K1at(Tf_max, Ts_max)); fprintf('\n');

fprintf('\nSlow gain K2(tau1,tau2), tau2 in [Ts_min,Ts_max]:\n');
fprintf('  at (0, Ts_min)      : '); fprintf('%9.4f', K2at(0, Ts_min));      fprintf('\n');
fprintf('  at (Tf_max, Ts_max) : '); fprintf('%9.4f', K2at(Tf_max, Ts_max)); fprintf('\n');

%% ------------------------------------------------------------------------
%  Independent verification on a dense grid
%
%  The convex-hull argument already certifies the design exactly. This grid
%  check is a safeguard against construction errors (e.g. misaligned knots,
%  wrong index sets) and is cheap enough to always run. Sample points are
%  placed strictly inside each segment so the finite-difference stencil
%  never straddles a knot, where Y is only C^0.
%  ------------------------------------------------------------------------

verify_design(Yctrl_val, L1ctrl_val, L2ctrl_val, knots1, knots2, d1, d2, ...
              A_full, A_J1_ol, A_J2_ol, B_J1, B_J2, Tf_min, Ts_min);

%% ------------------------------------------------------------------------
%  Closed-loop simulation with randomized sampling
%
%  Each inter-sample interval is drawn uniformly from its admissible range,
%  so the simulation exercises the aperiodic dwell-time set the LMIs cover.
%  Clocks start at zero and the first reset of each loop occurs after one
%  full drawn interval, so every gain evaluation uses a certified clock
%  value. The scheduled gains are evaluated at the actual clock values
%  immediately before each reset.
%  ------------------------------------------------------------------------

x0     = [2; -1];
z0     = [x0; 0; 0];

t_f = draw_sampling_times(Tf_min, Tf_max, Tfinal);
t_s = draw_sampling_times(Ts_min, Ts_max, Tfinal);

[tsim, Zsim] = simulate_scheduled(A_full, A_J1_ol, A_J2_ol, B_J1, B_J2, ...
                                  K1at, K2at, t_f, t_s, Tfinal, z0);

figure('Name','States (designed controller)');
plot(tsim, Zsim(:,1), 'k--', 'LineWidth', 2); hold on
plot(tsim, Zsim(:,2), 'b',   'LineWidth', 2); grid on
xlabel('time $t$','Interpreter','latex','FontSize',20);
ylabel('states','Interpreter','latex','FontSize',20);
legend('State $x_1$','State $x_2$','Interpreter','latex','FontSize',12);
set(gca,'FontSize',12);

figure('Name','Inputs (designed controller)');
plot(tsim, Zsim(:,3), 'k--', 'LineWidth', 2); hold on
plot(tsim, Zsim(:,4), 'b',   'LineWidth', 2); grid on
xlabel('time $t$','Interpreter','latex','FontSize',20);
ylabel('inputs','Interpreter','latex','FontSize',20);
legend('Input $\bar{u}_1$','Input $\bar{u}_2$','Interpreter','latex', ...
       'FontSize',12,'Location','SouthEast');
set(gca,'FontSize',12);

fprintf('Final state norm: %.3e\n', norm(Zsim(end,1:n)));


%% ========================================================================
%  Local functions
%  ========================================================================

function [knots, idxJump] = build_knots(Tmin, Tmax, m_pre, m_post, d)
% Composite Bezier knot vector on [0,Tmax] with an exact breakpoint at Tmin.
%
%   m_pre    segments cover [0, Tmin]     (no jump constraint)
%   m_post   segments cover [Tmin, Tmax]  (jump constrained)
%
% idxJump lists the global control-point indices belonging to the segments
% in [Tmin, Tmax]. Because segments share boundary control points, this is
% exactly (m_pre*d + 1) : (m_pre + m_post)*d + 1, a union of whole patches.

    assert(Tmin > 0 && Tmax > Tmin, 'Require 0 < Tmin < Tmax.');

    pre   = linspace(0, Tmin, m_pre + 1);
    post  = linspace(Tmin, Tmax, m_post + 1);
    knots = [pre, post(2:end)];

    nctrl   = (m_pre + m_post)*d + 1;
    idxJump = (m_pre*d + 1) : nctrl;
end

function V = value_or_zero(M, n_full)
% value() for decision variables, passthrough for the fixed zero blocks.
    if isa(M, 'sdpvar')
        V = value(M);
    else
        V = zeros(1, n_full);
    end
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
% Control points of patch (seg1,seg2) together with the directional
% derivative dC/dtau1 + dC/dtau2, exactly degree-elevated back to (d1,d2)
% so that it lives in the same Bernstein basis as the patch itself.
% h1, h2 are this patch's own segment widths (knots may be non-uniform).

    base1 = (seg1-1)*d1;
    base2 = (seg2-1)*d2;

    Cloc = cell(d1+1, d2+1);
    for i = 0:d1
        for j = 0:d2
            Cloc{i+1,j+1} = Cctrl{base1+i+1, base2+j+1};
        end
    end

    % derivative along axis 1: degree (d1-1, d2), elevated back to d1
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

    % derivative along axis 2: degree (d1, d2-1), elevated back to d2
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
% Segment of a possibly non-uniform knot vector containing t, and the local
% Bernstein parameter s in [0,1] within that segment.
    m = numel(knots) - 1;
    t = min(max(t, knots(1)), knots(end));
    seg = m;
    for k = 1:m
        if t <= knots(k+1) + 1e-12
            seg = k;
            break
        end
    end
    s = (knots(seg+1) - knots(seg));
    s = min(max((t - knots(seg))/s, 0), 1);
end

function val = eval_bezier2D(Cctrl_val, t1, t2, knots1, knots2, d1, d2)
% Evaluate a numeric tensor-product composite Bezier surface at (t1,t2).
    [seg1, s1] = locate_segment(t1, knots1);
    [seg2, s2] = locate_segment(t2, knots2);

    base1 = (seg1-1)*d1;
    base2 = (seg2-1)*d2;

    val = 0;
    for i = 0:d1
        bi = nchoosek(d1,i) * s1^i * (1-s1)^(d1-i);
        for j = 0:d2
            bj = nchoosek(d2,j) * s2^j * (1-s2)^(d2-j);
            val = val + bi*bj*Cctrl_val{base1+i+1, base2+j+1};
        end
    end
end

function pts = segment_interior_grid(knots, nPerSeg)
% Sample points strictly inside each segment, so that finite differences
% never straddle a knot (where the composite surface is only C^0).
    pts = [];
    for k = 1:numel(knots)-1
        e = linspace(knots(k), knots(k+1), nPerSeg + 2);
        pts = [pts, e(2:end-1)]; %#ok<AGROW>
    end
end

function verify_design(Yctrl_val, L1ctrl_val, L2ctrl_val, knots1, knots2, ...
                       d1, d2, A_full, A_J1_ol, A_J2_ol, B_J1, B_J2, ...
                       Tf_min, Ts_min)
% Dense-grid re-check of the closed-loop conditions in the ORIGINAL (primal)
% coordinates X = Y^{-1}, using the recovered gains. This tests the
% synthesis result rather than merely re-evaluating the LMIs that were
% handed to the solver.

    nPerSeg = 6;
    g1 = segment_interior_grid(knots1, nPerSeg);
    g2 = segment_interior_grid(knots2, nPerSeg);

    hmin = min([diff(knots1), diff(knots2)]);
    e    = 1e-6 * hmin;

    Yat  = @(t1,t2) eval_bezier2D(Yctrl_val,  t1, t2, knots1, knots2, d1, d2);
    L1at = @(t1,t2) eval_bezier2D(L1ctrl_val, t1, t2, knots1, knots2, d1, d2);
    L2at = @(t1,t2) eval_bezier2D(L2ctrl_val, t1, t2, knots1, knots2, d1, d2);

    worstPos = inf; worstFlow = -inf; worstJump = -inf;

    for t1 = g1
        for t2 = g2
            Y  = Yat(t1,t2);
            X  = eye(size(Y)) / Y;
            dY = (Yat(t1+e,t2+e) - Yat(t1-e,t2-e)) / (2*e);
            dX = -X*dY*X;

            worstPos  = min(worstPos,  min(eig(X)));
            worstFlow = max(worstFlow, max(eig(dX + A_full'*X + X*A_full)));

            if t1 >= Tf_min
                K1  = L1at(t1,t2) / Y;
                AK1 = A_J1_ol + B_J1*K1;
                X0  = eye(size(Y)) / Yat(0,t2);
                worstJump = max(worstJump, max(eig(AK1'*X0*AK1 - X)));
            end
            if t2 >= Ts_min
                K2  = L2at(t1,t2) / Y;
                AK2 = A_J2_ol + B_J2*K2;
                X0  = eye(size(Y)) / Yat(t1,0);
                worstJump = max(worstJump, max(eig(AK2'*X0*AK2 - X)));
            end
        end
    end

    fprintf('\nGrid verification in primal coordinates (%d x %d points)\n', ...
            numel(g1), numel(g2));
    fprintf('  min eig X = Y^{-1}            : %+.3e  (want > 0)\n', worstPos);
    fprintf('  max eig flow LMI              : %+.3e  (want < 0)\n', worstFlow);
    fprintf('  max eig closed-loop jump LMIs : %+.3e  (want < 0)\n', worstJump);

    if worstPos > 0 && worstFlow < 0 && worstJump < 0
        fprintf('  design verified on the grid.\n\n');
    else
        warning('Grid verification failed -- design is not certified.');
    end
end

function t = draw_sampling_times(Tmin, Tmax, Tfinal)
% Sampling instants with inter-sample intervals drawn i.i.d. uniformly from
% [Tmin, Tmax]. The clock starts at zero; the first instant is one full
% interval later, so every reset occurs at an admissible clock value.
    t  = [];
    tk = 0;
    while true
        tk = tk + Tmin + (Tmax - Tmin)*rand();
        if tk > Tfinal, break; end
        t(end+1) = tk; %#ok<AGROW>
    end
end

function [tvec, Z] = simulate_scheduled(A_full, A_J1_ol, A_J2_ol, B_J1, B_J2, ...
                                        K1at, K2at, t_f, t_s, Tfinal, z0)
% Simulate the impulsive closed loop with clock-scheduled gains, evaluated
% at the actual clock values immediately before each reset. Integration runs
% in the augmented coordinates z = [x; u1_hold; u2_hold], so the simulated
% model is exactly the model the LMIs certify and the held inputs are read
% directly off the state. Duplicated time stamps at reset instants render
% the inputs as proper staircases.
%
% Exactly simultaneous fast and slow resets have probability zero under
% randomized sampling; should they occur they are applied fast-then-slow.

    events = unique([t_f(:); t_s(:)]);
    events = events(events <= Tfinal);
    isFast = ismember(events, t_f(:));
    isSlow = ismember(events, t_s(:));

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

    z = z0(:);
    tvec = 0; Z = z.';
    tcur = 0;
    lastFast = 0; lastSlow = 0;

    for k = 1:numel(events)+1
        if k <= numel(events), t_event = events(k); else, t_event = Tfinal; end

        if t_event > tcur
            [tt, zz] = ode45(@(t,zz) A_full*zz, [tcur t_event], z, opts);
            tvec = [tvec; tt]; Z = [Z; zz]; %#ok<AGROW>
            z = zz(end,:).';
            tcur = t_event;
        end

        if k <= numel(events)
            tau1 = t_event - lastFast;
            tau2 = t_event - lastSlow;

            if isFast(k)
                z = A_J1_ol*z + B_J1*(K1at(tau1,tau2)*z);
                lastFast = t_event;
            end
            if isSlow(k)
                z = A_J2_ol*z + B_J2*(K2at(tau1,tau2)*z);
                lastSlow = t_event;
            end

            tvec = [tvec; t_event]; Z = [Z; z.']; %#ok<AGROW>
        end
    end
end
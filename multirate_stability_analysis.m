%% MULTIRATE_STABILITY_ANALYSIS
%
% Stability analysis of a multirate sampled-data control loop via a
% clock-dependent Lyapunov certificate parameterized as a tensor-product
% composite Bezier (B-spline) matrix function.
%
% SETUP
%   Continuous-time plant       xdot = A x + B1 u1 + B2 u2
%   with two independent, asynchronous zero-order-hold loops:
%     - "fast" loop: u1 <- K1 x, inter-sample times in [Tf_min, Tf_max]
%     - "slow" loop: u2 <- K2 x, inter-sample times in [Ts_min, Ts_max]
%   The gains K1, K2 are FIXED here; see multirate_controller_design.m
%   for the synthesis counterpart.
%
% IMPULSIVE REFORMULATION
%   On the augmented state z = [x; u1_hold; u2_hold] the closed loop is
%     flow:       zdot = A_full z          (plant evolves, holds frozen)
%     fast reset: z^+  = A_J1 z            (u1_hold <- K1 x)
%     slow reset: z^+  = A_J2 z            (u2_hold <- K2 x)
%
% CERTIFICATE
%   P(tau1,tau2) depends on the two clocks tau_i = time since the last
%   reset of loop i, and is written in a tensor-product Bernstein basis.
%   By the convex-hull property of that basis, imposing a linear matrix
%   inequality at every control point of a patch certifies it on the
%   entire patch -- exactly, with no gridding.
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
%   certificate is not valid (it can report feasible for a closed loop that
%   diverges in simulation).
%
%   The knot vectors below therefore place an EXACT breakpoint at T_min on
%   each axis. The jump-constrained index sets are then read off directly
%   from the segment structure and are unions of whole patches by
%   construction, for any ratio T_max / T_min.
% -------------------------------------------------------------------------

clear; close all; clc;

rng(1);   % reproducible sampling sequences in the simulation

%% ------------------------------------------------------------------------
%  Plant
%  ------------------------------------------------------------------------

A  = [ 2.0  0.50;
      -0.5  0.25];        % open loop unstable

B1 = [1; 0];              % actuated by the fast loop
B2 = [0; 1];              % actuated by the slow loop

n = size(A,1);

%% ------------------------------------------------------------------------
%  Dwell-time bounds and fixed controller gains
%  ------------------------------------------------------------------------

h = 0;

Tf_min = 0.03;   Tf_max = 0.05+h;     % fast loop inter-sample range
Ts_min = 0.30;   Ts_max = 0.50+h;     % slow loop inter-sample range

K1 = -[4 0];     % fast loop gain (uses x1)
K2 = -[0 2];     % slow loop gain (uses x2)

%% ------------------------------------------------------------------------
%  Augmented impulsive representation of the closed loop
%  ------------------------------------------------------------------------

A_full = [A       B1  B2;
          zeros(2,4)     ];
n_full = size(A_full,1);

B_J1 = [0;0;1;0];        % injects the refreshed fast hold
B_J2 = [0;0;0;1];        % injects the refreshed slow hold

A_J1_ol = diag([1 1 0 1]);   % fast reset, structural part
A_J2_ol = diag([1 1 1 0]);   % slow reset, structural part

A_J1 = A_J1_ol + B_J1*[K1 0 0];   % fast reset, closed loop
A_J2 = A_J2_ol + B_J2*[K2 0 0];   % slow reset, closed loop

%% ------------------------------------------------------------------------
%  Composite Bezier parameterization of P(tau1,tau2)
%  ------------------------------------------------------------------------

tol = 1e-4;      % strict-inequality margin used in all LMIs

d1 = 2;  m1_pre = 1;  m1_post = 2;   % fast axis: degree, segments before/after Tf_min
d2 = 2;  m2_pre = 1;  m2_post = 2;   % slow axis: degree, segments before/after Ts_min

[knots1, idxJump1] = build_knots(Tf_min, Tf_max, m1_pre, m1_post, d1);
[knots2, idxJump2] = build_knots(Ts_min, Ts_max, m2_pre, m2_post, d2);

m1 = numel(knots1) - 1;   nctrl1 = m1*d1 + 1;
m2 = numel(knots2) - 1;   nctrl2 = m2*d2 + 1;

Pctrl = cell(nctrl1, nctrl2);
for p = 1:nctrl1
    for q = 1:nctrl2
        Pctrl{p,q} = sdpvar(n_full, n_full, 'symmetric');
    end
end

I = eye(n_full);
F = [];

%% ------------------------------------------------------------------------
%  Positive definiteness, reduced to the faces of the clock box
%
%  Instead of imposing P > 0 on the whole box, it suffices to impose it on
%  the faces {tau_i = 0}. Given the flow and jump inequalities below,
%  positivity on the faces propagates to the whole box: from any tau, flow
%  forward along the diagonal by delta = min_i(Tmax_i - tau_i); since
%  d/ds[exp(A's) P exp(As)] < 0 along the flow, P(tau) dominates the value at
%  the endpoint, where some clock has reached its upper bound and the jump
%  inequality applies, landing on a face.
%
%  In the tensor basis the face {tau_1 = 0} is exactly the composite Bezier
%  curve carried by the FIRST ROW of control points, and {tau_2 = 0} by the
%  FIRST COLUMN, so the convex-hull property makes this restriction exact.
%  The constraint count drops from prod_i s_i*(d_i+1) to nctrl1+nctrl2-1
%  (81 -> 13 with the default settings).
%  ------------------------------------------------------------------------

for q = 1:nctrl2
    F = [F, Pctrl{1,q} >= tol*I];          % face tau_1 = 0
end
for p = 2:nctrl1                            % face tau_2 = 0 (corner already done)
    F = [F, Pctrl{p,1} >= tol*I];
end

%% ------------------------------------------------------------------------
%  Flow LMIs
%     dP/dtau1 + dP/dtau2 + A_full' P + P A_full < 0
%  on [0,Tf_max] x [0,Ts_max]
%  ------------------------------------------------------------------------

for seg1 = 1:m1
    h1 = knots1(seg1+1) - knots1(seg1);
    for seg2 = 1:m2
        h2 = knots2(seg2+1) - knots2(seg2);

        [Cloc, dC] = bezier_patch_derivative(Pctrl, seg1, seg2, d1, d2, h1, h2);

        for i = 0:d1
            for j = 0:d2
                Pt  = Cloc{i+1,j+1};
                dPt = dC{i+1,j+1};

                F = [F, dPt + A_full'*Pt + Pt*A_full <= -tol*I];
            end
        end
    end
end

%% ------------------------------------------------------------------------
%  Jump LMIs
%     fast:  A_J1' P(0,tau2) A_J1 - P(tau1,tau2) < 0,  tau1 in [Tf_min,Tf_max]
%     slow:  A_J2' P(tau1,0) A_J2 - P(tau1,tau2) < 0,  tau2 in [Ts_min,Ts_max]
%
%  P(0,tau2) is exactly the boundary curve carried by the first row of
%  control points, P(tau1,0) by the first column. Since the Bernstein
%  weights of the omitted axis sum to one, imposing the inequality at every
%  control-point pair of a patch certifies it on the whole patch.
%  ------------------------------------------------------------------------

for p = idxJump1
    for q = 1:nctrl2
        F = [F, A_J1'*Pctrl{1,q}*A_J1 - Pctrl{p,q} <= -tol*I];
    end
end

for p = 1:nctrl1
    for q = idxJump2
        F = [F, A_J2'*Pctrl{p,1}*A_J2 - Pctrl{p,q} <= -tol*I];
    end
end

%% ------------------------------------------------------------------------
%  Normalization, conditioning, objective
%  ------------------------------------------------------------------------

F = [F, trace(Pctrl{1,1}) == n_full];     % removes the scale invariance

for p = 1:nctrl1
    for q = 1:nctrl2
        F = [F, -100*I <= Pctrl{p,q} <= 100*I];
    end
end

obj = 0;
for p = 1:nctrl1
    for q = 1:nctrl2
        obj = obj + sum(sum(Pctrl{p,q}.^2));
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
    disp('Infeasible or solver failure -- no certificate obtained.');
    return
end

res = min(check(F));
fprintf('Worst LMI residual at the solution: %.3e\n', res);
if res < -1e-5
    warning('Reported feasible but residuals are violated; solve is inconclusive.');
end

Pctrl_val = cellfun(@value, Pctrl, 'UniformOutput', false);

%% ------------------------------------------------------------------------
%  Independent verification on a dense grid
%
%  The convex-hull argument already certifies the LMIs exactly. This grid
%  check is a safeguard against construction errors (e.g. misaligned knots,
%  wrong index sets) and is cheap enough to always run. Sample points are
%  placed strictly inside each segment so the finite-difference stencil
%  never straddles a knot, where P is only C^0.
%
%  It also confirms the positivity reduction directly: P > 0 is imposed only
%  on the two faces, but is checked here over the whole box.
%  ------------------------------------------------------------------------

verify_analysis(Pctrl_val, knots1, knots2, d1, d2, A_full, A_J1, A_J2, ...
                Tf_min, Ts_min);

%% ------------------------------------------------------------------------
%  Closed-loop simulation with randomized sampling
%
%  Each inter-sample interval is drawn uniformly from its admissible range,
%  so the simulation exercises the aperiodic dwell-time set the LMIs cover.
%  Clocks start at zero and the first reset of each loop occurs after one
%  full drawn interval, consistent with the range dwell-time condition.
%  ------------------------------------------------------------------------

Tfinal = 5;
x0     = [2; -1];
z0     = [x0; 0; 0];

t_f = draw_sampling_times(Tf_min, Tf_max, Tfinal);
t_s = draw_sampling_times(Ts_min, Ts_max, Tfinal);

[tsim, Zsim] = simulate_fixed_gain(A_full, A_J1, A_J2, t_f, t_s, Tfinal, z0);

figure('Name','States');
plot(tsim, Zsim(:,1), 'k--', 'LineWidth', 2); hold on
plot(tsim, Zsim(:,2), 'b',   'LineWidth', 2); grid on
xlabel('time $t$','Interpreter','latex','FontSize',20);
ylabel('states','Interpreter','latex','FontSize',20);
legend('State $x_1$','State $x_2$','Interpreter','latex','FontSize',12);
set(gca,'FontSize',12);

figure('Name','Inputs');
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

function verify_analysis(Pctrl_val, knots1, knots2, d1, d2, A_full, ...
                         A_J1, A_J2, Tf_min, Ts_min)
% Dense-grid re-check of every LMI, independent of the solver's report.

    nPerSeg = 6;
    g1 = segment_interior_grid(knots1, nPerSeg);
    g2 = segment_interior_grid(knots2, nPerSeg);

    hmin = min([diff(knots1), diff(knots2)]);
    e    = 1e-6 * hmin;                 % stays inside the current segment

    Pat = @(t1,t2) eval_bezier2D(Pctrl_val, t1, t2, knots1, knots2, d1, d2);

    worstPos = inf; worstFlow = -inf; worstJump = -inf;

    for t1 = g1
        for t2 = g2
            P  = Pat(t1,t2);
            dP = (Pat(t1+e,t2+e) - Pat(t1-e,t2-e)) / (2*e);

            worstPos  = min(worstPos,  min(eig(P)));
            worstFlow = max(worstFlow, max(eig(dP + A_full'*P + P*A_full)));

            if t1 >= Tf_min
                M = A_J1'*Pat(0,t2)*A_J1 - P;
                worstJump = max(worstJump, max(eig(M)));
            end
            if t2 >= Ts_min
                M = A_J2'*Pat(t1,0)*A_J2 - P;
                worstJump = max(worstJump, max(eig(M)));
            end
        end
    end

    fprintf('\nGrid verification (%d x %d points)\n', numel(g1), numel(g2));
    fprintf('  min eig P                     : %+.3e  (want > 0)\n', worstPos);
    fprintf('  max eig flow LMI              : %+.3e  (want < 0)\n', worstFlow);
    fprintf('  max eig jump LMIs             : %+.3e  (want < 0)\n', worstJump);

    if worstPos > 0 && worstFlow < 0 && worstJump < 0
        fprintf('  certificate verified on the grid.\n\n');
    else
        warning('Grid verification failed -- certificate is not valid.');
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

function [tvec, Z] = simulate_fixed_gain(A_full, A_J1, A_J2, t_f, t_s, Tfinal, z0)
% Simulate the impulsive closed loop with constant gains. Integration runs
% in the augmented coordinates z = [x; u1_hold; u2_hold], so the simulated
% model is exactly the model the LMIs certify and the held inputs are read
% directly off the state. Duplicated time stamps at reset instants render
% the inputs as proper staircases.

    events = unique([t_f(:); t_s(:)]);
    events = events(events <= Tfinal);
    isFast = ismember(events, t_f(:));
    isSlow = ismember(events, t_s(:));

    opts = odeset('RelTol',1e-8,'AbsTol',1e-10);

    z = z0(:);
    tvec = 0; Z = z.';
    tcur = 0;

    for k = 1:numel(events)+1
        if k <= numel(events), t_event = events(k); else, t_event = Tfinal; end

        if t_event > tcur
            [tt, zz] = ode45(@(t,zz) A_full*zz, [tcur t_event], z, opts);
            tvec = [tvec; tt]; Z = [Z; zz]; %#ok<AGROW>
            z = zz(end,:).';
            tcur = t_event;
        end

        if k <= numel(events)
            if isFast(k), z = A_J1*z; end
            if isSlow(k), z = A_J2*z; end
            tvec = [tvec; t_event]; Z = [Z; z.']; %#ok<AGROW>
        end
    end
end
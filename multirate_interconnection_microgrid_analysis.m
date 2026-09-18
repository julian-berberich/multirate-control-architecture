%% MULTIRATE_INTERCONNECTION_EXAMPLE
%
% Seven impulsive subsystems in a star, each sampled on its own clock, two of
% them sharing a clock and additionally exchanging information at those
% instants through a hardwired bus.
%
% -------------------------------------------------------------------------
% SUBSYSTEM MODEL  (notation of the paper)
%
%   Physical state  xi_j = (V_j, I_j, v_j): PCC voltage, filter current,
%   integrator state. Sampled-data implementation adds the zero-order-hold
%   register uhat_j, giving the impulsive state
%
%       x_j = (V_j, I_j, v_j, uhat_j) in R^4 .
%
%   FLOW, t not in T:      xdot_j = A_j x_j + B_e,j d_j
%                          e_j    = C_j x_j + D_e,j d_j
%
%       A_j = [   0      1/C_j    0     0   ]      B_e,j = [ 1/C_j ]
%             [ -1/L_j  -R_j/L_j  0   1/L_j ]              [   0   ]
%             [  -1        0      0     0   ]              [   0   ]
%             [   0        0      0     0   ]              [   0   ]
%
%       C_j = [1 0 0 0],   D_e,j = Rs,   B_j = 0,  D_j = 0.
%
%   The control does not act during flow: it enters only through the hold
%   state uhat_j, which is part of x_j. This matches the standing assumption
%   B_j = 0, D_j = 0 of the interconnection theorem.
%
%   JUMP, t in T_i with I_j = {i(j)} a single clock per subsystem:
%
%       x_j^+     = A_J,j x_j^- + B_J,j u_J,j + B_J,e,j d_J,j
%       e_J,j     = C_J,j x_j^-
%
%       A_J,j = diag(1,1,1,0),   B_J,j = e_4,   D_J,j = 0,  D_J,e,j = 0,
%       B_J,e,j = kappa*e_4  and  C_J,j = [0 1 0 0]   for j in {2,5},
%       B_J,e,j and C_J,j EMPTY                        otherwise.
%
%   Componentwise the reset reads
%
%       V_j^+ = V_j^-,  I_j^+ = I_j^-,  v_j^+ = v_j^-,
%       uhat_j^+ = K_J,j(theta^-) x_j^-  +  kappa * d_J,j ,
%                  \___ local feedback __/   \_ shared bus _/
%
%   so only the hold register changes; the physical state is untouched.
%
% -------------------------------------------------------------------------
% WHAT kappa IS
%
%   kappa is the single nonzero entry of B_J,e,j. It is SYSTEM DATA, like the
%   line resistances, NOT a controller: subsystems 2 and 5 are wired to a
%   common analog current-sharing bus whose gain is fixed in hardware. The
%   controller u_J,j = K_J,j(theta^-) x_j^- is purely local, exactly as in the
%   paper's controller parameterization. The neighbour's current reaches the
%   hold register only through the fixed path B_J,e,j, and the interconnection
%   theorem accounts for it through the jump supply rate P_J,j. That is what
%   the jump port is for; kappa cannot be folded into K_J,j because d_J,j is
%   not a function of x_j. Setting kappa = 0 makes B_J,e,j = 0 and removes the
%   jump coupling entirely.
%
% -------------------------------------------------------------------------
% WHY THE FEEDTHROUGH Rs IS THERE
%
%   Passivity on the flow port with NO feedthrough makes the (2,2) block of
%   the flow LMI identically zero, which forces the exact pinning
%       X B_e = sig*C'   i.e.  X(:,1) = sig*C_j*e_1 ,
%   and that makes strict dissipativity IMPOSSIBLE for a sampled-data
%   subsystem:
%     (i)  jumps leave the physical state unchanged, so A_J^K e_1 = e_1 + B_J k_1
%          and the (1,1) entry of the jump condition is k_1^2 (X^0)_44 >= 0.
%          Hence k_1 = 0: the reset cannot use the passivity output at all.
%     (ii) the pinning also kills every X_1k, so on the (I,v,uhat) block the
%          flow condition reduces to dXh + Ah'Xh + Xh Ah + eta Xh <= 0.
%     (iii) dim ker(Ah) = 2 and dim ker(Gh - I) = 2 with Gh the reset on that
%          block. Two 2-dimensional subspaces of R^3 always intersect, and for
%          w in the intersection the flow forces
%              w' Xh(theta) w <= exp(-eta*theta) w' Xh(0) w
%          while the jump forces w' Xh(theta) w >= w' Xh(0) w. Contradiction
%          for every eta > 0.
%
%   A series resistance Rs > 0 at the PCC, e_j = V_j + Rs*d_j, gives
%   D_e,j = Rs and turns the (2,2) block into -2*sig*Rs < 0. The pinning
%   becomes a penalised mismatch instead of an equality and the chain never
%   starts. The interconnection condition is unaffected, since it depends only
%   on M and the supply rates, and well-posedness det(I + Rs*L) ~= 0 is
%   automatic because L >= 0.
%
%   With Rs > 0 no load conductance is needed: the model below is the
%   reference's verbatim, with the loads left in the disturbance. Rs is the
%   ONLY addition to the plant. All LMIs are imposed STRICTLY, so the
%   dissipativity theorem applies exactly as stated.
%
% -------------------------------------------------------------------------
% SUPPLY RATES
%   Flow: P_j = [0 sig; sig 0], sig = 10, all j. With lambda_j = 1,
%         [M;I]'Z(lambda)[M;I] = -2*sig*L <= 0, structural.
%   Jump: P_J,j = [nu sigJ; sigJ 0] for j in {2,5}, empty otherwise.
%         The (2,2) block of the jump condition is kappa^2 (X^0)_44 - nu, so
%         nu MUST be positive; the interconnection needs nu <= sigJ. Taking
%         nu = sigJ makes [M_J;I]'Z_J[M_J;I] = 2(nu-sigJ)*L2 = 0 and leaves
%         the largest local margin. sigJ enters the OFF-diagonal independently
%         of kappa and must be small: sigJ = 5 is infeasible, 0.1 is not.
%
%   Note the sign: d = M e with M = -L, since d_j = sum_k (e_k - e_j)/R_jk.
%
% MODEL SOURCE
%   M. Tucci, S. Riverso, G. Ferrari-Trecate, "Line-independent plug-and-play
%   controllers for voltage stabilization in DC microgrids", IEEE TCST 26(3),
%   1115-1123, 2018; converter and line parameters from Table 2 of the
%   technical report arXiv:1609.02456. Added here: the sampling ranges (the
%   reference designs in continuous time), the PCC resistance Rs, and the
%   replication of three units to reach seven.
%
% Requires: YALMIP and an SDP solver (MOSEK by default).

clear; close all; clc;
rng(1);

%% ------------------------------------------------------------------------
%  Configuration
%  ------------------------------------------------------------------------

cfg.sig   = 10;         % flow-port passivity scaling (common to all units)
cfg.Rs    = 0.05;       % PCC series resistance -> feedthrough D_e,j = Rs
cfg.sigJ  = 0.1;        % jump-port scaling
cfg.nu    = 0.1;        % = sigJ, see header
cfg.kappa = 0.2;        % gain of the shared current-sharing bus (system data)
cfg.eta   = 0.02;       % guaranteed decay rate of the storage

cfg.d     = 2;          % Bezier degree
cfg.alpha = 0.40;       % adaptive rule: ||A||*h <= alpha
cfg.s_max = 8;          % cap on segments per region

cfg.tol   = 1e-6;       % strictness margin, flow LMI
cfg.tolJ  = 1e-6;       % strictness margin, jump LMI
cfg.resTol= -1e-7;      % residual acceptance
cfg.Ymin  = 1e-8;       % loose positivity guard  (inactive at the solution)
cfg.Ymax  = 1e4;        % loose boundedness guard (inactive at the solution)

cfg.options = sdpsettings('solver','mosek','verbose',0);
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-10;
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS   = 1e-10;
cfg.options.mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS   = 1e-10;

%% ------------------------------------------------------------------------
%  Subsystems, clocks, topology
%  ------------------------------------------------------------------------

N = 7;  n = 4;
Ccap = 2.20;                               % capacitance, common to all units

%          j =    1     2     3     4     5     6     7
R_all   = [ 0.2   0.3   0.1   0.5   0.3   0.1   0.5 ];
L_all   = [ 1.8   2.0   2.2   3.0   2.0   2.2   3.0 ];
clkOf   = [   1     2     3     4     2     5     6 ];   % G_2, G_5 share clock 2
Tmin_c  = [ 0.4   0.5   0.8   1.5   1.0   1.8 ];
Tmax_c  = [ 0.6   0.7   1.2   2.0   1.4   2.4 ];
nClk    = numel(Tmin_c);
pair    = [2 5];
assert(clkOf(pair(1)) == clkOf(pair(2)), 'Jump-coupled units must share a clock.');

R1j = [ 0.05  0.07  0.03  0.05  0.07  0.03 ];        % hub-to-leaf lines
Lap = zeros(N);
for k = 1:N-1
    wgt = 1/R1j(k);
    Lap(1,1)     = Lap(1,1) + wgt;
    Lap(k+1,k+1) = wgt;
    Lap(1,k+1)   = -wgt;
    Lap(k+1,1)   = -wgt;
end

Cy   = [1 0 0 0];                 % flow-port output   e_j = V_j + Rs*d_j
CJ   = [0 1 0 0];                 % jump-port output   e_J,j = I_j
A_J  = diag([1 1 1 0]);
B_J  = [0;0;0;1];
BJe  = cfg.kappa*[0;0;0;1];       % B_J,e,j : the shared bus, system data

A = cell(1,N);  Be = cell(1,N);
for j = 1:N
    Rj = R_all(j);  Lj = L_all(j);
    Ahat = [  0      1/Ccap   0 ;
             -1/Lj  -Rj/Lj    0 ;
             -1      0        0 ];
    Bhat = [0; 1/Lj; 0];
    A{j}  = [ Ahat, Bhat ; zeros(1,4) ];
    Be{j} = [1/Ccap; 0; 0; 0];
end

%% ------------------------------------------------------------------------
%  Report
%  ------------------------------------------------------------------------

fprintf('=====================================================================\n');
fprintf(' MULTIRATE INTERCONNECTION: %d subsystems, %d clocks\n', N, nClk);
fprintf('=====================================================================\n');
fprintf('  sig = %g, Rs = %g, sigJ = nu = %g, kappa = %g, eta = %g\n', ...
        cfg.sig, cfg.Rs, cfg.sigJ, cfg.kappa, cfg.eta);
fprintf('  loads left in the disturbance (no conductance): reference model verbatim\n');
fprintf('  open-loop eigenvalues of the four distinct types:\n');
for j = 1:4
    ev = eig(A{j}(1:3,1:3));
    fprintf('    G_%d : ', j); fprintf('%+.4f%+.4fi  ', [real(ev).'; imag(ev).']);
    fprintf('\n');
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  Local design: one SDP per (dynamics, clock) pair
%  ------------------------------------------------------------------------

Kfun = cell(1,N);  ok = false(1,N);
fprintf('--- local design -----------------------------------------------------\n');
for j = 1:N
    if j == pair(2) && ok(pair(1))
        Kfun{j} = Kfun{pair(1)};  ok(j) = true;
        fprintf('  G_%d : same problem as G_%d\n', j, pair(1));
        continue
    end
    ci = clkOf(j);  hasJ = ismember(j, pair);
    [Kfun{j}, info] = design_local(A{j}, Be{j}, Cy, CJ, A_J, B_J, BJe, ...
                                   Tmin_c(ci), Tmax_c(ci), cfg, hasJ);
    ok(j) = ~isempty(Kfun{j});
    fprintf('  G_%d : clock %d, T in [%.1f,%.1f], seg %d+%d, %d ctrl pts', ...
            j, ci, Tmin_c(ci), Tmax_c(ci), info.spre, info.spost, info.N);
    if ok(j)
        fprintf('  ->  feasible (resid %.1e)\n', info.resid);
    else
        fprintf('  ->  INFEASIBLE (problem %g)\n', info.problem);
    end
end
if ~all(ok)
    error(['Local design failed. After a parameter change the first knobs are ' ...
           'cfg.Rs (raise) and cfg.sigJ (lower).']);
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  Interconnection conditions
%  ------------------------------------------------------------------------

fprintf('--- interconnection conditions ---------------------------------------\n');
eL  = eig((Lap+Lap.')/2);
eL2 = eig([1 -1; -1 1]);
eWP = eig(eye(N) + cfg.Rs*Lap);
fprintf('  flow  -2*sig*L <= 0        : min eig(L)  = %+.3e   (need >= 0)\n', min(eL));
fprintf('  jump  2(nu-sigJ)*L2 <= 0   : nu - sigJ   = %+.3e   (need <= 0)\n', cfg.nu-cfg.sigJ);
fprintf('                               min eig(L2) = %+.3e   (need >= 0)\n', min(eL2));
fprintf('  well-posed det(I+Rs*L)~=0  : min eig     = %+.4f    (need > 0)\n', min(eWP));
if ~(min(eL) >= -1e-10 && (cfg.nu-cfg.sigJ) <= 1e-12 && min(eL2) >= -1e-10 && min(eWP) > 0)
    error('Interconnection condition violated.');
end
fprintf('  all hold  ->  the interconnection is exponentially stable.\n\n');

%% ------------------------------------------------------------------------
%  Dense-grid verification of each local certificate (primal coordinates)
%  ------------------------------------------------------------------------

fprintf('--- grid verification (X = Y^{-1}, recovered gains) ------------------\n');
for j = 1:N
    if j == pair(2), continue; end
    ci = clkOf(j);
    v = verify_local(Kfun{j}, A{j}, Be{j}, Cy, CJ, A_J, B_J, BJe, ...
                     Tmin_c(ci), cfg, ismember(j,pair));
    fprintf('  G_%d : min eig X = %+.3e (>0)   flow = %+.3e (<0)   jump = %+.3e (<0)\n', ...
            j, v.pos, v.flow, v.jump);
    if ~(v.pos > 0 && v.flow < 0 && v.jump < 0)
        warning('verify:failed','G_%d failed grid verification.', j);
    end
end
fprintf('\n');

%% ------------------------------------------------------------------------
%  Closed-loop and open-loop simulation
%
%  The feedthrough makes the coupling implicit: e = V + Rs*d and d = -L e
%  give d = -Ltil*V with Ltil = L*(I + Rs*L)^{-1}, which is again PSD.
%  ------------------------------------------------------------------------

Tfinal = 200;                                   % ms
Ltil   = Lap/(eye(N) + cfg.Rs*Lap);

Aglob = zeros(N*n);
for j = 1:N
    rows = (j-1)*n + (1:n);
    Aglob(rows,rows) = A{j};
    for k = 1:N
        cols = (k-1)*n + (1:n);
        Aglob(rows,cols) = Aglob(rows,cols) - Ltil(j,k)*Be{j}*Cy;
    end
end

Z0 = zeros(N*n,1);
V0 = [ 1.0  -0.6   0.8  -0.4   0.5  -0.9   0.7 ];
I0 = [-0.5   0.3  -0.2   0.6  -0.4   0.2  -0.3 ];
for j = 1:N
    Z0((j-1)*n + (1:3)) = [V0(j); I0(j); 0];
end

[tsim, Zsim] = simulate_network(Aglob, A_J, B_J, CJ, cfg.kappa, Kfun, ...
                                clkOf, pair, Tmin_c, Tmax_c, Tfinal, Z0, n, N);
xnorm = zeros(numel(tsim), N);  uhold = zeros(numel(tsim), N);
for j = 1:N
    idx = (j-1)*n + (1:n);
    xnorm(:,j) = vecnorm(Zsim(:,idx(1:3)), 2, 2);
    uhold(:,j) = Zsim(:,idx(4));
end

% open loop: u = 0, so the holds stay at zero and there are no resets
[tol_, Zol] = ode45(@(t,z) Aglob*z, [0 Tfinal], Z0, ...
                    odeset('RelTol',1e-9,'AbsTol',1e-12));
xnorm_ol = zeros(numel(tol_), N);
for j = 1:N
    xnorm_ol(:,j) = vecnorm(Zol(:,(j-1)*n + (1:3)), 2, 2);
end

sel  = tsim > 5;
pfit = polyfit(tsim(sel), log(max(max(xnorm(sel,:),[],2),1e-300)), 1);
evOL = eig(Aglob);

fprintf('--- simulation -------------------------------------------------------\n');
fprintf('  OPEN LOOP (u = 0)\n');
fprintf('    max Re eig = %+.2e, %d eigenvalues at the origin, geometric\n', ...
        max(real(evOL)), sum(abs(evOL) < 1e-9));
fprintf('    multiplicity %d: the zero eigenvalue is DEFECTIVE, so generic\n', ...
        N*n - rank(Aglob,1e-9));
fprintf('    initial conditions grow linearly. The holds start at zero here,\n');
fprintf('    which does not excite those modes, so the response stays bounded --\n');
fprintf('    V and I decay while the integrator states settle at nonzero\n');
fprintf('    constants. The origin is NOT asymptotically stable.\n');
fprintf('    max_j ||x_j||:  %.3e at t=0   ->   %.3e at t=%g\n', ...
        max(xnorm_ol(1,:)), max(xnorm_ol(end,:)), Tfinal);
fprintf('  CLOSED LOOP\n');
fprintf('    max_j ||x_j||:  %.3e at t=0   ->   %.3e at t=%g\n', ...
        max(xnorm(1,:)), max(xnorm(end,:)), Tfinal);
fprintf('    empirical decay rate %.4f /ms   (certificate uses eta = %g)\n', ...
        -pfit(1), cfg.eta);
fprintf('  improvement at t=%g: factor %.0f\n\n', ...
        Tfinal, max(xnorm_ol(end,:))/max(xnorm(end,:)));

%% ------------------------------------------------------------------------
%  Plots
%  ------------------------------------------------------------------------

col = lines(N);
lgd = arrayfun(@(j) sprintf('$G_%d$',j), 1:N, 'UniformOutput', false);

ylo   = max(min([xnorm(:); xnorm_ol(:)]), 1e-14);
yhi   = max([xnorm(:); xnorm_ol(:)]);
ylim_ = [10^(floor(log10(ylo))), 10^(ceil(log10(yhi)))];

figure('Name','Open-loop states','Color','w','Position',[80 500 620 420]);
hold on
for j = 1:N, plot(tol_, max(xnorm_ol(:,j),1e-14), 'LineWidth', 1.4, 'Color', col(j,:)); end
set(gca,'YScale','log','FontSize',12); grid on; ylim(ylim_); xlim([0 Tfinal]);
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$\|x_j\|$','Interpreter','latex','FontSize',16);
legend(lgd,'Interpreter','latex','FontSize',11,'Location','southwest');
hold off

figure('Name','Closed-loop states','Color','w','Position',[720 500 620 420]);
hold on
for j = 1:N, plot(tsim, max(xnorm(:,j),1e-14), 'LineWidth', 1.4, 'Color', col(j,:)); end
set(gca,'YScale','log','FontSize',12); grid on; ylim(ylim_); xlim([0 Tfinal]);
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$\|x_j\|$','Interpreter','latex','FontSize',16);
legend(lgd,'Interpreter','latex','FontSize',11,'Location','southwest');
hold off

figure('Name','Closed-loop inputs','Color','w','Position',[400 20 620 420]);
hold on
for j = 1:N, plot(tsim, uhold(:,j), 'LineWidth', 1.4, 'Color', col(j,:)); end
grid on; set(gca,'FontSize',12); xlim([0 Tfinal]);
xlabel('$t$','Interpreter','latex','FontSize',16);
ylabel('$\hat{u}_j$','Interpreter','latex','FontSize',16);
legend(lgd,'Interpreter','latex','FontSize',11,'Location','northeast');
hold off


%% ========================================================================
%  Local functions
%  ========================================================================

function [Kfun, info] = design_local(A, Be, Cy, CJ, A_J, B_J, BJe, ...
                                     Tmin, Tmax, cfg, hasJump)
% Dual-variable design of a clock-dependent reset gain rendering one
% subsystem strictly dissipative w.r.t. its flow port and, if hasJump, its
% jump port. One clock, so the certificate is a univariate spline.
%
%   flow :  [ -Ydot + A Y + Y A' + eta Y     Be - sig*Y*Cy' ]  <  0
%           [            *                     -2*sig*Rs    ]
%   jump :  [ -Y   -sigJ*Y*CJ'   (A_J Y + B_J L)' ]
%           [  *      -nu             BJe'        ]  <  0
%           [  *       *              -Y0         ]
% (without a jump port the middle row and column are absent). Both are
% imposed STRICTLY, so the dissipativity theorem applies as stated.

    yalmip('clear');
    Kfun = [];
    info = struct('spre',NaN,'spost',NaN,'N',NaN,'resid',NaN,'problem',NaN);

    n = size(A,1);  d = cfg.d;

    hmax  = cfg.alpha / norm(A);
    spre  = min(cfg.s_max, max(1, ceil(Tmin/hmax)));
    spost = min(cfg.s_max, max(1, ceil((Tmax-Tmin)/hmax)));
    knots = [linspace(0,Tmin,spre+1), linspace(Tmin,Tmax,spost+1)];
    knots(spre+2) = [];                       % drop the duplicated Tmin
    m  = spre + spost;
    Np = m*d + 1;
    idxJ = (spre*d + 1) : Np;                 % control points with theta >= Tmin
    info.spre = spre;  info.spost = spost;  info.N = Np;

    Yc = cell(1,Np);  Lc = cell(1,Np);
    for p = 1:Np
        Yc{p} = sdpvar(n,n,'symmetric');
        Lc{p} = sdpvar(1,n,'full');
    end

    F = {};
    for p = 1:Np
        F{end+1} = (Yc{p} >= cfg.Ymin*eye(n));                        %#ok<*AGROW>
        F{end+1} = (Yc{p} <= cfg.Ymax*eye(n));
    end

    for seg = 1:m
        h = knots(seg+1) - knots(seg);
        [Yl, dY] = bezier_seg_derivative(Yc, seg, d, h);
        for i = 0:d
            Yt  = Yl{i+1};
            off = Be - cfg.sig*Yt*Cy.';
            M = [ -dY{i+1} + A*Yt + Yt*A.' + cfg.eta*Yt,  off ;
                   off.',                                -2*cfg.sig*cfg.Rs ];
            F{end+1} = (M <= -cfg.tol*eye(n+1));
        end
    end

    Y0 = Yc{1};
    for p = idxJ
        Mj = A_J*Yc{p} + B_J*Lc{p};
        if hasJump
            F{end+1} = ([ -Yc{p},                 -cfg.sigJ*Yc{p}*CJ.',  Mj.'  ;
                          -cfg.sigJ*CJ*Yc{p},     -cfg.nu,               BJe.' ;
                           Mj,                     BJe,                 -Y0   ] ...
                        <= -cfg.tolJ*eye(2*n+1));
        else
            F{end+1} = ([ -Yc{p},  Mj.' ;
                           Mj,    -Y0  ] <= -cfg.tolJ*eye(2*n));
        end
    end

    obj = 0;
    for p = idxJ, obj = obj + sum(Lc{p}.^2); end

    try
        sol = optimize([F{:}], obj, cfg.options);
    catch ME
        warning('design_local:solver','%s', ME.message);
        info.problem = -99;  return
    end

    info.problem = sol.problem;
    if sol.problem ~= 0, return; end
    info.resid = min(check([F{:}]));
    if info.resid <= cfg.resTol, return; end

    Yv = cellfun(@value, Yc, 'UniformOutput', false);
    Lv = cellfun(@value, Lc, 'UniformOutput', false);
    Kfun = struct('K',     @(th) eval_bezier1D(Lv, th, knots, d) ...
                              / eval_bezier1D(Yv, th, knots, d), ...
                  'Y',     @(th) eval_bezier1D(Yv, th, knots, d), ...
                  'knots', knots, 'd', d);
end

function v = verify_local(S, A, Be, Cy, CJ, A_J, B_J, BJe, Tmin, cfg, hasJump)
% Dense-grid re-check of the local certificate in the ORIGINAL coordinates.
% Sample points sit strictly inside each segment, so the finite-difference
% stencil never straddles a knot (the spline is only C^0 there).
    knots = S.knots;  d = S.d;  n = size(A,1);
    th = [];
    for k = 1:numel(knots)-1
        e = linspace(knots(k), knots(k+1), 9);
        th = [th, e(2:end-1)];                                        %#ok<AGROW>
    end
    e = 1e-7;
    v = struct('pos',inf,'flow',-inf,'jump',-inf);
    for t = th
        Y  = S.Y(t);   X = eye(n)/Y;
        dY = (S.Y(min(t+e,knots(end))) - S.Y(max(t-e,0)))/(2*e);
        dX = -X*dY*X;
        v.pos = min(v.pos, min(eig((X+X.')/2)));

        off = X*Be - cfg.sig*Cy.';
        M = [ dX + A.'*X + X*A + cfg.eta*X,  off ;
              off.',                        -2*cfg.sig*cfg.Rs ];
        v.flow = max(v.flow, max(eig((M+M.')/2)));

        if t >= Tmin - 1e-12
            K  = S.K(t);  AK = A_J + B_J*K;  X0 = eye(n)/S.Y(0);
            if hasJump
                Mj = [ AK.'*X0*AK - X,                    AK.'*X0*BJe - cfg.sigJ*CJ.' ;
                      (AK.'*X0*BJe - cfg.sigJ*CJ.').',    BJe.'*X0*BJe - cfg.nu       ];
            else
                Mj = AK.'*X0*AK - X;
            end
            v.jump = max(v.jump, max(eig((Mj+Mj.')/2)));
        end
    end
end

function Ce = bezier_elevate1D(C, n)
% Exact Bezier degree elevation, degree n -> degree n+1.
    Ce = cell(1,n+2);
    for k = 0:n+1
        term = 0;
        if k >= 1, term = term + (k/(n+1))     * C{k};   end
        if k <= n, term = term + (1 - k/(n+1)) * C{k+1}; end
        Ce{k+1} = term;
    end
end

function [Cloc, dC] = bezier_seg_derivative(Cctrl, seg, d, h)
% Control points of segment seg and the derivative d/dtheta, degree-elevated
% back to d so that it lives in the same Bernstein basis.
    base = (seg-1)*d;
    Cloc = cell(1,d+1);
    for i = 0:d, Cloc{i+1} = Cctrl{base+i+1}; end
    Q = cell(1,d);
    for i = 0:d-1, Q{i+1} = d*(Cloc{i+2} - Cloc{i+1}); end
    Qe = bezier_elevate1D(Q, d-1);
    dC = cell(1,d+1);
    for i = 0:d, dC{i+1} = Qe{i+1}/h; end
end

function val = eval_bezier1D(Cval, t, knots, d)
    m = numel(knots) - 1;
    t = min(max(t, knots(1)), knots(end));
    seg = m;
    for k = 1:m
        if t <= knots(k+1) + 1e-12, seg = k; break; end
    end
    s = min(max((t - knots(seg))/(knots(seg+1) - knots(seg)), 0), 1);
    base = (seg-1)*d;
    val = 0;
    for i = 0:d
        val = val + nchoosek(d,i)*s^i*(1-s)^(d-i)*Cval{base+i+1};
    end
end

function t = draw_ticks(Tmin, Tmax, Tfinal)
% Inter-sample intervals drawn uniformly from [Tmin,Tmax]. Clocks start at
% zero, the first tick one full interval later, so every reset happens at a
% certified clock value.
    t = [];  tk = 0;
    while true
        tk = tk + Tmin + (Tmax - Tmin)*rand();
        if tk > Tfinal, break; end
        t(end+1) = tk;                                                %#ok<AGROW>
    end
end

function [tvec, Z] = simulate_network(Aglob, A_J, B_J, CJ, kappa, Kfun, ...
                                      clkOf, pair, Tmin_c, Tmax_c, Tfinal, Z0, n, N)
% Impulsive closed loop with independent, aperiodic clocks. At a tick of
% clock i every subsystem on that clock resets; the jump-coupled pair also
% exchanges currents at that instant through the shared bus.

    nClk  = numel(Tmin_c);
    ticks = cell(1,nClk);
    for i = 1:nClk, ticks{i} = draw_ticks(Tmin_c(i), Tmax_c(i), Tfinal); end

    events = unique([ticks{:}]);
    events = events(events <= Tfinal);
    opts   = odeset('RelTol',1e-9,'AbsTol',1e-12);

    Z = Z0(:).';  tvec = 0;  tcur = 0;  last = zeros(1,nClk);

    for k = 1:numel(events)+1
        if k <= numel(events), te = events(k); else, te = Tfinal; end

        if te > tcur
            [tt, zz] = ode45(@(t,z) Aglob*z, [tcur te], Z(end,:).', opts);
            tvec = [tvec; tt(2:end)];  Z = [Z; zz(2:end,:)];          %#ok<AGROW>
            tcur = te;
        end
        if k > numel(events), break; end

        z = Z(end,:).';
        for i = 1:nClk
            if ~any(abs(ticks{i} - te) < 1e-12), continue; end
            th  = te - last(i);
            mem = find(clkOf == i);
            if isequal(sort(mem(:).'), sort(pair(:).'))
                a = pair(1);  b = pair(2);
                ia = (a-1)*n + (1:n);  ib = (b-1)*n + (1:n);
                Ka = Kfun{a}.K(th);
                za = z(ia);  zb = z(ib);
                z(ia) = A_J*za + B_J*(Ka*za) + kappa*B_J*(CJ*zb - CJ*za);
                z(ib) = A_J*zb + B_J*(Ka*zb) + kappa*B_J*(CJ*za - CJ*zb);
            else
                for j = mem
                    ij = (j-1)*n + (1:n);
                    z(ij) = A_J*z(ij) + B_J*(Kfun{j}.K(th)*z(ij));
                end
            end
            last(i) = te;
        end
        tvec = [tvec; te];  Z = [Z; z.'];                             %#ok<AGROW>
    end
end
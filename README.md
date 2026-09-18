# Multirate Control Architecture

MATLAB code accompanying

> J. Berberich, *A Robust Control Architecture for Multirate Systems and their Interconnections*.

The scripts reproduce the three numerical examples of the paper: stabilizing
multirate control, robust multirate $\mathcal{H}_\infty$-control, and stability
analysis of an interconnection of impulsive systems.

---

## What the code does

Multirate systems are modeled as impulsive systems with several independent
jump-time sequences, one clock $\theta_i$ per sequence. Stability and
dissipativity are certified by a Lyapunov-like function

$$V(t) = x(t)^\top \mathcal{X}(\theta(t))\, x(t),$$

whose matrix depends on all clocks. The resulting conditions are differential
LMIs on the clock box $[0,\bar T_1]\times\cdots\times[0,\bar T_q]$.

The differential LMIs are relaxed into a finite-dimensional SDP by writing
$\mathcal{X}$ (or its inverse $\mathcal{Y}$, for synthesis) in a
**tensor-product composite Bézier basis**. By the convex-hull property of the
Bernstein basis, imposing an LMI at every control point of a patch certifies it
on the whole patch — exactly, with no gridding.

Two implementation points are load-bearing and are worth knowing before
modifying the code:

**Knot placement.** The jump LMIs must hold only for
$\tau_i \in [\underline T_i, \bar T_i]$. Since the convex-hull argument
certifies an LMI on a patch only if it is imposed at *every* control point of
that patch, the restriction is exact only when $\underline T_i$ coincides with a
segment boundary. All scripts therefore place a knot at $\underline T_i$ on each
axis. Selecting control points instead by a numerical test such as
`grid value >= T_min` yields certificates that are not valid.

**Adaptive refinement.** The spline resolution needed on axis $i$ scales with
$\|A\|\,\bar T_i$, not with $\bar T_i$ alone. The scripts that sweep over
sampling rates choose segment counts per axis so that every segment satisfies
$\|A\|\,h \le \alpha$. Without this, points get reported as "not certified"
because the relaxation ran out of resolution rather than because no controller
exists.

---

## Requirements

- MATLAB (tested with R2023b or later)
- [YALMIP](https://yalmip.github.io/)
- An SDP solver. All scripts default to **MOSEK**; to use a different solver,
  edit the `sdpsettings` call near the top of each file. SeDuMi and SDPT3 work
  but are slower and may need looser tolerances.

No toolboxes beyond base MATLAB are required. Each script is self-contained:
all helper routines (Bézier evaluation, degree elevation, knot construction,
verification, simulation) are local functions at the end of the file.

---

## Files

| File | Paper section | Runtime |
|---|---|---|
| `multirate_stability_analysis.m` | 5.A | seconds |
| `multirate_controller_design.m` | 5.A | seconds |
| `multirate_hinf_design_once.m` | 5.B | ~1 min |
| `multirate_hinf_design_tradeoff.m` | 5.B, Fig. 5 | ~1–2 h |
| `multirate_interconnection_microgrid_analysis.m` | 5.C | ~1 min |

Runtimes are indicative, measured with MOSEK on a desktop machine.

### `multirate_stability_analysis.m`

Stability analysis (Theorem 1) for a two-rate sampled-data loop with **fixed**
gains $K_1 = -[4\ \ 0]$, $K_2 = -[0\ \ 2]$. Solves the differential LMIs for a
clock-dependent Lyapunov matrix $\mathcal{X}(\tau_1,\tau_2)$, re-checks the
certificate on a dense grid in the original coordinates, and simulates the
closed loop with inter-sample times drawn uniformly at random from the
range dwell-time intervals.

Set `h > 0` to enlarge both $\bar T_i$; the LMIs become infeasible for the fixed
gains, which motivates the synthesis script.

### `multirate_controller_design.m`

Synthesis of clock-scheduled reset gains (Proposition 2) via the dual
formulation: solve for $\mathcal{Y} = \mathcal{X}^{-1}$ and auxiliary factors
$L_i$, then recover $K_i = L_i\mathcal{Y}^{-1}$. The jump conditions enter as
Schur complements and stay jointly affine in $(\mathcal{Y}, L_i)$.

In the zero-order-hold representation the augmented flow is autonomous, so the
flow condition is instantiated with $B = 0$, $L = 0$ and reduces to a Lyapunov
differential inequality on $\mathcal{Y}$. **All design freedom sits in the two
jump conditions.**

Setting `RESTRICT_OWN_CLOCK = true` constrains each $K_i$ to depend on its own
clock $\tau_i$ only — the more implementable structure — which lets the two
variants be compared numerically.

### `multirate_hinf_design_once.m`

Robust $\mathcal{L}_2$-gain study at one pair of dwell-time bounds, with a
polytopic uncertainty $\beta \in [0.6, 1]$ in the slow actuator. Computes three
quantities under the same relaxation:

1. $\gamma$ for the given gains, **primal** formulation ($\mathcal{X}$ splined);
2. $\gamma$ for the same gains, **dual** formulation, substituting
   $L_{\mathrm{J},i} := K_i\mathcal{Y}$ — linear in $\mathcal{Y}$ because $K_i$
   is constant;
3. $\gamma$ for **designed** clock-scheduled gains, same dual LMIs, $L$ free.

Two consistency checks are built in. (3) $\le$ (2) is **rigorous**: same
parameterization, and the fixed-gain feasible set is a subset of the design one,
so a violation is a bug, not conservatism. (1) versus (2) is a
cross-formulation check on the same quantity through independent code paths;
the two need not agree at a fixed mesh, since a spline's inverse is not a
spline, but they should converge under refinement — which the built-in
refinement check (`DO_REFINEMENT`) exercises by halving `cfg.alpha`.

Toggles: `DO_VERIFY`, `DO_SIMULATE`, `DO_REFINEMENT`, and `DO_SPEC_DESIGN`
(re-solve at a prescribed $\gamma$, minimizing gain magnitude).

### `multirate_hinf_design_tradeoff.m`

Sweeps the robust design over the sampling-rate plane and produces the
rate-allocation figure. Both axes run over $[0.03, 1]$ on a logarithmic grid,
so the box includes the region where the nominally "fast" loop is slower than
the nominally "slow" one; the $\bar T_1 = \bar T_2$ diagonal is drawn for
reference.

Relative jitter is held fixed, $\underline T_i = \kappa\,\bar T_i$ with
$\kappa = 0.6$, so changes in $\gamma^\star$ are attributable to the sampling
**rate** rather than to a changing degree of aperiodicity.

Outputs iso-$\gamma$ contours, a family of sampling-budget curves
$1/\bar T_1 + 1/\bar T_2 = \nu$ for $\nu \in \{4,6,9,14,22\}$ Hz, the
$\gamma^\star$-minimizing point on each, and the locus of those optima. Results
are saved to `rate_tradeoff_results.mat`.

This is the long-running script. Set `USE_PARFOR = true` if the Parallel
Computing Toolbox is available, or reduce `n1`, `n2` and `n_curve` for a
coarser but much faster sweep.

> **Note on symbols.** In this script `beta` is the uncertain plant parameter
> and `nu` the sampling budget. They are deliberately distinct.

### `multirate_interconnection_microgrid_analysis.m`

Compositional stability analysis (Theorem 5) of a DC microgrid with $N = 7$
distributed generation units in a star topology. Each unit is a single-rate
impulsive system with its own clock, except $G_2$ and $G_5$, which share a
clock and additionally exchange information at those instants through a
hardwired bus — so $q = 6$ clocks in total.

The script designs local controllers rendering each unit strictly passive
(Proposition 3), checks the two coupling conditions of Theorem 5, re-verifies
each local certificate on a dense grid, and simulates the network with all six
clocks independent and aperiodic.

Both coupling conditions hold **structurally**, for any connected topology and
any line resistances, with $\lambda_j = 1$ and no search required:

- flow: $-2\bar\sigma\mathcal{L} \preceq 0$, since $\mathcal{L}\succeq 0$;
- jump: $2(\nu - \bar\sigma_{\mathrm{J}})\mathcal{L}_2 \preceq 0$, satisfied with
  equality for $\nu = \bar\sigma_{\mathrm{J}}$.

Two modelling points deserve attention if you change parameters.

**The series resistance `cfg.Rs` is not cosmetic.** With no feedthrough, the
$(2,2)$ block of the passivity LMI vanishes, which forces the exact pinning
$\mathcal{X}B_{\mathrm d} = \bar\sigma C^\top$. Since jumps leave the plant state
unchanged, the $(1,1)$ entry of the jump condition then equals
$k_1^2(\mathcal{X}^{0})_{44} \ge 0$, so the reset cannot depend on the passivity
output at all, and strict dissipativity becomes structurally impossible. Setting
$R_{\mathrm s} > 0$ makes the pinning a penalized mismatch instead of an equality
and the chain never starts.

**The jump-port scaling `cfg.sigJ` must be small.** It enters the off-diagonal of
the jump LMI independently of `cfg.kappa`, so it costs margin whatever $\kappa$
is: $\bar\sigma_{\mathrm J} = 5$ is infeasible for every $\kappa$ tried, while
$0.1$ leaves comfortable margin. And $\nu > 0$ is *forced* — the $(2,2)$ block of
the jump condition is $\kappa^2(\mathcal{X}^{0})_{44} - \nu$, so a jump port
always requires an input-strict supply rate.

If a local design turns out infeasible after a parameter change, the first two
knobs are `cfg.Rs` (raise) and `cfg.sigJ` (lower).

---

## Reproducing the figures

| Paper figure | Script | Generated figures |
|---|---|---|
| Fig. 3 (states, inputs) | `multirate_stability_analysis.m` | `States`, `Inputs` |
| Fig. 4 (rate allocation) | `multirate_hinf_design_tradeoff.m` | `Rate allocation`, `Gain along budget curves` |
| Fig. 6 (open/closed loop) | `multirate_interconnection_microgrid_analysis.m` | `Open-loop states`, `Closed-loop states` |

Simulations draw inter-sample times at random. The scripts that simulate call
`rng` with a fixed seed where reproducibility matters; qualitative behaviour is
insensitive to the seed, exact trajectories are not.

---

## Reading the output

Each script prints a diagnostic block before its figures. Worth reading rather
than skipping, since it distinguishes failure modes that look identical on a
plot:

- **infeasible** — the LMIs have no solution at this relaxation;
- **solved but residual-rejected** — the solver converged but
  `min(check(F))` fell below `cfg.resTol`; the worst residual is printed, and
  the threshold should be loosened before reading this as infeasible;
- **segment-cap-limited** — the adaptive refinement hit `cfg.s_max`, so an
  infeasible result is *inconclusive* rather than a feasibility boundary.

The strictness margin `cfg.tol` and the residual threshold `cfg.resTol` satisfy
$|\texttt{resTol}| < \texttt{tol}$, so an accepted solution satisfies
$M \preceq -(\texttt{tol} - |\texttt{resTol}|)I$, i.e. genuine strict
feasibility.

Where $\gamma$ is the reported quantity, **no regularization is placed on $L$**.
A penalty $\rho\sum\|L\|^2$ with $\rho = 10^{-4}$ contributes $O(1)$–$O(10^2)$
to the objective for these control-point counts and would bias $\gamma$ upward;
$L$ is bounded by a box instead. A gain penalty *is* used in the fixed-$\gamma$
design, where $\gamma$ is prescribed rather than measured.

---

## Citing

```bibtex
@article{berberich2026multirate,
  author  = {Berberich, Julian},
  title   = {A Robust Control Architecture for Multirate Systems
             and their Interconnections},
  year    = {2026}
}
```

## License

MIT — see `LICENSE`.

## Acknowledgment

The code as well as this Readme file were generated with the assistance of Claude (Opus 5) and subsequently
reviewed and tested by the author, who maintains full responsibility for its
content.

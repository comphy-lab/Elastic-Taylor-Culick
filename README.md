# Elastic Taylor--Culick

Taylor--Culick retraction cases in Basilisk C, organised using the CoMPhy
project layout. There are two geometries -- axisymmetric hole opening and a
planar semi-infinite sheet -- sharing one constitutive path: the scalar
2D/axisymmetric log-conformation solver from
[MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow).

## Repository layout

```
├── simulationCases/ - Active Basilisk entry points and generated case folders
│   ├── TaylorCulick.c - Axisymmetric elastic and viscoelastic hole-opening case
│   └── TaylorCulickPlanar.c - Planar semi-infinite sheet, same constitutive path
├── src-local/ - Stable rheology headers and typed parameter parser
│   ├── log-conform-viscoelastic-scalar-2D.h - Upstream scalar log-conformation solver
│   ├── two-phaseVE.h - Per-phase VE material-property coupling
│   ├── parse_params.h - key=value file loader
│   └── params.h - typed runtime accessors
├── postProcess/tip_to_csv.py - Tip log to t,tstar,x_tip,v_tip,v_over_VTC CSV
├── scripts/params.sh - Shared shell parameter helpers
├── default.params - Purely elastic axisymmetric default (lambda1 = 1e30)
├── default-viscoelastic.params - Finite-relaxation axisymmetric default
├── default-planar.params - Planar Newtonian default (Oh_SB = 0.1)
├── sweep.params - Two-combination sweep definition
├── runSimulation.sh - Single-case compile/run driver
├── runParameterSweep.sh - Deterministic Cartesian sweep driver
├── AGENTS.md - Repository guidance
└── README.md - Project documentation
```

Generated case directories such as simulationCases/c1000/ are ignored and
contain the copied source, case.params, executable, logs, dumps, and snapshots
for that run.

## Requirements

- Basilisk C with qcc available in PATH
- bash, awk, and standard POSIX utilities

The local source snapshot was copied from MultiRheoFlow commit
4695a434f0750c8476c9e094bb8561384092299a on 2026-08-06. The scalar solver
matches upstream blob cea628698212c4172780b35ec8949fbfb1bc4570 exactly.
The companion two-phaseVE.h was whitespace-normalized from upstream blob
537123b1390f798ad7262321f7505b8c83efa1db.

## Single-case runs

Compile and run the elastic default:

```bash
bash runSimulation.sh --input default.params
```

Run the finite-relaxation case:

```bash
bash runSimulation.sh --input default-viscoelastic.params
```

Compile without running:

```bash
bash runSimulation.sh --input default.params --compile-only
```

Run the planar case:

```bash
bash runSimulation.sh --case simulationCases/TaylorCulickPlanar.c \
  --input default-planar.params
```

The executable receives case.params as argv[1]; this keeps each output
directory self-describing and restartable.

`--outdir DIR` puts the run directory somewhere other than
`simulationCases/c<CaseNo>/`, which is what you want when run data must stay
out of the checkout. `--openmp` adds `-fopenmp`, so the run honours
`OMP_NUM_THREADS`.

## Planar case

`simulationCases/TaylorCulickPlanar.c` retracts a semi-infinite sheet of
**full** thickness `h0` that is symmetric about the midplane `y = 0`. The free
edge is closed by a semicircular cap of radius `h0/2` centred on the midplane
at `x = xtip0 + h0/2`, so the interface meets the midplane at `x = xtip0` and
retracts towards `+x`. The bottom boundary keeps Basilisk's default symmetry
condition, which is exactly the midplane mirror condition; left, right and top
are open.

It is the same code as the axisymmetric case apart from the geometry: no
`axi.h`, a planar initial condition and boundary conditions, no `AThTh`
conformation component (the solver already guards it with `#if AXI`), and no
`2*pi*y` weight in the kinetic-energy integral. The elastic and viscoelastic
capability is unchanged -- `G1`, `lambda1`, `G2`, `lambda2` behave exactly as
in the axisymmetric case, and `G1 = 0` gives the Newtonian limit by parameter
rather than by a separate stripped source file.

`Ldomain` must be large enough that the sheet stays effectively semi-infinite:
the edge travels at most `sqrt(2)*t`, so keep `Ldomain` well above
`xtip0 + sqrt(2)*tmax`.

### Non-dimensionalisation

Lengths are scaled with the **full** thickness `h0`, densities with the liquid
density, and stresses with `sigma/h0`. With `rho1 = sigma = h0 = 1`:

- `mu1` is the Ohnesorge number `Oh = mu/sqrt(rho*sigma*h0)`;
- the Taylor--Culick speed is `V_TC = sqrt(2*sigma/(rho*h0)) = sqrt(2)`;
- the capillary time is `sqrt(rho*h0^3/sigma) = 1`.

Savva & Bush (*JFM* **626**, 2009) put the **half**-thickness in their
Ohnesorge number,

```
Oh_SB = mu/sqrt(2*h0*rho*sigma) = Oh/sqrt(2),
```

so a case quoted at `Oh_SB` is run here with `mu1 = sqrt(2)*Oh_SB` -- the
factor of `sqrt(2)` is easy to lose. Their viscous time is
`tau_vis = mu*h0/(2*sigma) = mu1/2` and their reduced time is
`t* = t/tau_vis`. Both clocks are written to the tip file.

### Tip diagnostic

The planar case measures the retraction in-code and writes `c<CaseNo>-tip.dat`
with columns `t`, `tstar`, `x_tip`, `x_tip_vof`, `x_tip_global`:

- `x_tip` is the smallest `x` over the reconstructed VOF facets inside the
  midplane band `y < h0/10`, i.e. the interface position on the midplane;
- `x_tip_vof` is the independent estimate `integral (1-f) dx` along the bottom
  row of cells;
- `x_tip_global` is the same minimum taken over the whole interface.

All three are sub-cell accurate. `x_tip` and `x_tip_vof` agree while the
midplane is crossed exactly once, so a growing gap between them flags rim
pinch-off or an entrained bubble rather than a genuine tip motion.

Differentiate offline:

```bash
python3 postProcess/tip_to_csv.py simulationCases/c1000
```

which writes `tip_velocity.csv` with `t,tstar,x_tip,v_tip,v_over_VTC`. The
speed is a local first-order least-squares slope over `--window` samples,
which is far less noisy than a two-point difference at the sub-cell scale of
the VOF reconstruction.

## Phase rheology

two-phaseVE.h maps the parameters below to the scalar solver:

- f = 1: rho1, mu1, G1, lambda1
- f = 0: rho2, mu2, G2, lambda2

Use a finite positive lambda for a viscoelastic liquid. Use lambda = 1e30
for the upstream purely elastic limit. Both phases can therefore be assigned
independently, including an elastic solid and a finite-relaxation liquid.
- `TOLelastic` is the phase-fraction cutoff used at mixed cells; it defaults to
  `1e-2` and should remain small compared with one. The runtime range is
  `0 <= TOLelastic < 0.5`, which keeps at least one phase active in every
  mixed cell.

## Parameter sweep

Preview the checked-in two-case sweep:

```bash
bash runParameterSweep.sh --dry-run
```

CASE_START..CASE_END must equal the Cartesian product generated by all
SWEEP_* entries. Each generated case is executed through runSimulation.sh,
so the single-case contract remains the only compilation path.

## Validation boundary

The local checks establish parser, runner, and compilation contracts. They do
not by themselves verify Taylor--Culick convergence or validate the physical
model against independent data. Those claims require separate refinement and
external-comparison evidence.

/**
# Elastic Taylor--Culick retraction

Axisymmetric two-phase Taylor--Culick retraction using the CoMPhy scalar
log-conformation solver.  The same case supports a finite-relaxation
viscoelastic phase and an elastic phase by assigning independent `G1`,
`lambda1`, `G2`, and `lambda2` values in the parameter file.

The constitutive implementation is the stable scalar 2D/axisymmetric solver
from [MultiRheoFlow](https://github.com/comphy-lab/MultiRheoFlow), while
`two-phaseVE.h` maps the phase properties onto its `Gp` and `lambda` fields.

## Phase convention

- `f = 1`: dense retracting phase, with `rho1`, `mu1`, `G1`, `lambda1`.
- `f = 0`: surrounding phase, with `rho2`, `mu2`, `G2`, `lambda2`.
- `lambda = 1e30`: upstream solver convention for the purely elastic limit.

## Runtime parameters

Parameters are loaded from a `key=value` file through `src-local/params.h`.
The root runner passes the copied `case.params` file to the executable.
*/

#include "axi.h"
#include "navier-stokes/centered.h"

#define FILTERED
#include "log-conform-viscoelastic-scalar-2D.h"
#include "two-phaseVE.h"

#include "navier-stokes/conserving.h"
#include "tension.h"
#define PARSE_PARAMS_IMPLEMENTATION
#include "params.h"
#undef PARSE_PARAMS_IMPLEMENTATION

/**
## Numerical controls

The tolerances are deliberately explicit so that a refinement study can
change them in one place without changing the constitutive model.
*/
#define FERR 1e-3
#define VELERR 1e-6
#define AERR 1e-4
#define KERR 1e-6

/**
## Outer boundary conditions

The top and right boundaries are open outflow boundaries.  In `axi.h`, the
bottom boundary is the symmetry axis; the local scalar solver applies the
axis condition to the off-diagonal conformation component.
*/
u.n[top] = neumann(0.);
p[top] = dirichlet(0.);
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

int CaseNo, MAXlevel, MINlevel;
double tmax, Ldomain;
char dumpFile[128], logFile[128], snapshotFile[160];

/**
## Periodic-event increments need a non-zero static initialiser

`qcc` registers events inside `_init_solver()`, which runs before any
statement of `main()`.  `init_event()` classifies each event expression
there by calling it twice and watching whether `i` or `t` change; an
expression that leaves both unchanged is taken to be a *condition* rather
than an *increment*, and that classification is never redone.

So `event writing_files (t = 0.; t += tsnap)` with `tsnap` assigned only
inside `main()` is registered while `tsnap == 0.`, is misclassified as a
condition, and fires exactly once, at `t = 0`.  Nothing warns: the run
completes normally, writes its final `dump`, and leaves
`intermediate/` holding a single snapshot.

A non-zero file-scope initialiser makes the classification correct.  The
runtime value is still what gets used, because `init_event()` re-runs at
`iter == 0` -- after `main()` -- so `tsnap` stays settable from
`case.params`.
*/
double tsnap = 0.25;

/**
### main()

Loads runtime parameters, assigns both phase rheologies, and starts the
Basilisk event loop.
*/
int main (int argc, char const * argv[])
{
  params_init_from_argv(argc, argv);

  CaseNo = param_int("CaseNo", 1000);
  MAXlevel = param_int("MAXlevel", 10);
  MINlevel = param_int("MINlevel", max(6, MAXlevel - 4));
  Ldomain = param_double("Ldomain", 100.);
  tmax = param_double("tmax", 25.);
  tsnap = param_double("tsnap", 0.25);
  dtmax = param_double("dtmax", 1e-5);

  rho1 = param_double("rho1", 1.);
  mu1 = param_double("mu1", 5e-2);
  rho2 = param_double("rho2", 1e-3);
  mu2 = param_double("mu2", 1e-5);

  G1 = param_double("G1", param_double("Ec", 1.));
  lambda1 = param_double("lambda1", param_double("De", 1e30));
  G2 = param_double("G2", 0.);
  lambda2 = param_double("lambda2", 0.);
  TOLelastic = param_double("TOLelastic", 1e-2);

  if (CaseNo < 1000 || MAXlevel < 1 || MAXlevel > 20 ||
      MINlevel < 1 || MINlevel > MAXlevel || Ldomain <= 0. ||
      tmax <= 0. || tsnap <= 0. || dtmax <= 0. || dtmax > tmax ||
      rho1 <= 0. || rho2 <= 0. || mu1 < 0. || mu2 < 0. ||
      G1 < 0. || G2 < 0. || lambda1 < 0. || lambda2 < 0. ||
      TOLelastic < 0. || TOLelastic >= 0.5) {
    fprintf(ferr, "ERROR: invalid runtime parameters.\n");
    return 1;
  }

  L0 = Ldomain;
  X0 = 0.;
  Y0 = 0.;
  init_grid(1 << MINlevel);

  if (system("mkdir -p intermediate") != 0) {
    fprintf(ferr, "ERROR: unable to create intermediate output directory.\n");
    return 1;
  }

  sprintf(dumpFile, "dump");
  sprintf(logFile, "c%d-log", CaseNo);

  f.sigma = 1.;
  TOLERANCE = 1e-4;
  CFL = 0.5;

  /**
  ## The timestep cap must be written to `DT`, not to `dtmax`

  `navier-stokes/centered.h` contains

  ~~~literatec
  event set_dtmax (i++,last) dtmax = DT;
  ~~~

  so `dtmax` is re-read from the global `DT` at the top of *every*
  timestep.  Assigning `dtmax` here would bind the first step only; from
  `i = 1` onwards the cap silently reverts to `DT`.  The failure is
  resolution-dependent -- a coarse run survives and looks converged while
  the same script blows up at higher `MAXlevel` -- which is exactly
  backwards from how a refinement study should behave.  Setting `DT` makes
  the cap persistent, which is what `dtmax=` in the parameter file is
  meant to express.

  ## The capillary constraint is `tension.h`'s job, and it does it

  `tension.h` contributes its own `stability` event that narrows `dtmax`
  to the capillary limit

  $$\Delta t_\sigma = \sqrt{\rho_m\Delta_{min}^3/(\pi\sigma)}$$

  using the *actual* finest cell carrying an interface, not a worst-case
  guess.  That narrowing was checked against this solver stack with an
  instrumented `tension.h`: at `i = 0` it sees finite `dmin`, `rho_m` and
  `sigma` and does apply the cap.  So the constraint is reported below
  rather than re-imposed here -- clamping `DT` to the `MAXlevel` value
  would pin every run to the finest cell the adaptation *could* produce
  even while the interface still sits on coarser cells.

  `dt_sigma` is printed at `MAXlevel` as the worst case the run can reach,
  which is what to compare an achieved `dt` against when a high-resolution
  run misbehaves.
  */
  const double rho_m = (rho1 + rho2)/2.;
  const double Delta_min = Ldomain/(1 << MAXlevel);
  const double dt_sigma = sqrt(rho_m*cube(Delta_min)/(pi*f.sigma));
  DT = dtmax;

  if (pid() == 0) {
    fprintf(ferr,
            "CaseNo=%d MAXlevel=%d MINlevel=%d Ldomain=%g "
            "tmax=%g dtmax_requested=%g dt_sigma=%g DT=%g\n",
            CaseNo, MAXlevel, MINlevel, Ldomain, tmax, dtmax,
            dt_sigma, DT);
    fprintf(ferr,
            "phase1: rho=%g mu=%g G=%g lambda=%g; "
            "phase2: rho=%g mu=%g G=%g lambda=%g\n",
            rho1, mu1, G1, lambda1, rho2, mu2, G2, lambda2);
    fprintf(ferr, "TOLelastic=%g\n", TOLelastic);
  }

  run();
}

/**
## Initial interface

The initial condition is the circular cap and retracting sheet used by the
original Taylor--Culick cases.
*/
event init (t = 0)
{
  const double hole0 = 1.;
  const double h0 = 1.;

  if (!restore(file = dumpFile)) {
    /**
    The `- Delta` terms are load-bearing.  `refine()` evaluates its
    condition at cell *centres*, so a plain `x < h0/2 + pad` test never
    fires on the initial grid whenever `Ldomain/2^MINlevel` exceeds about
    `h0`: the first column of cells has its centre outside the band it is
    meant to resolve, nothing is refined, and `fraction()` then writes one
    smeared value across the whole sheet.  The run still starts and still
    produces a plausible-looking retraction, but the initial interface is
    wrong -- and nothing in the output says so.  Testing the cell's inner
    edge instead of its centre makes the pre-refinement independent of
    `MINlevel`.

    The band is resolved along its whole radial extent first, because the
    flat part of the sheet is just as under-resolved as the rim is; the
    second pass then takes the rim neighbourhood to `MAXlevel`.  The
    padding is a fraction of `h0` rather than a fixed `0.025` so that it
    stays proportionate to the sheet thickness.
    */
    const double pad = 0.1*h0;
    refine(x - Delta < h0/2. + pad && level < MAXlevel - 3);
    refine(x - Delta < h0/2. + pad &&
           y - Delta < hole0 + 2.*h0 && level < MAXlevel);
    fraction(f, y < hole0 + h0/2.
             ? sq(h0/2.) - (sq(x) + sq(y - h0/2. - hole0))
             : h0/2. - x);

    /**
    ## Fail loudly on a wrong initial condition

    The smeared-initial-condition failure above costs one wasted campaign
    precisely because it is silent, so the initial state is checked against
    its analytic value before the first timestep rather than trusted.

    In `axi.h` the metric is `cm = y`, so `sum f dv()` is the liquid volume
    divided by `2*pi`.  For the flat sheet `0 < x < h0/2`,
    `hole0 + h0/2 < y < L0`, that is
    `(h0/2)(L0^2 - (hole0 + h0/2)^2)/2`; the semicircular rim of radius
    `h0/2` centred on `(0, hole0 + h0/2)` contributes its area times its
    centroid radius, `(pi h0^2/8)(hole0 + h0/2)`.

    The tolerance is loose enough that ordinary VOF discretisation error
    never trips it and tight enough that the smearing failure -- which was
    wrong by a factor of ninety -- cannot get past.
    */
    const double yrim = hole0 + h0/2.;
    const double expected = (h0/2.)*(sq(L0) - sq(yrim))/2.
      + (pi*sq(h0)/8.)*yrim;
    double vol = 0.;
    foreach (reduction(+:vol))
      vol += f[]*dv();

    if (pid() == 0)
      fprintf(ferr, "initial liquid volume/(2 pi) = %g (expected %g, "
              "relative error %.3g)\n",
              vol, expected, fabs(vol - expected)/expected);
    if (fabs(vol - expected) > 0.05*expected) {
      fprintf(ferr, "ERROR: initial volume fraction is wrong; the sheet is "
              "probably unresolved on the initial grid.\n");
      return 1;
    }
  }
}

/**
## Adaptive mesh refinement

The interface, velocity, curvature, and all axisymmetric conformation
components participate in the wavelet error estimate.
*/
scalar KAPPA[];

event adapt_mesh (i++)
{
  curvature(f, KAPPA);
  adapt_wavelet((scalar *) {f, u.x, u.y, KAPPA,
                            A11, A12, A22, AThTh},
                (double[]) {FERR, VELERR, VELERR, KERR,
                             AERR, AERR, AERR, AERR},
                MAXlevel);
}

/**
## Restart and snapshots
*/
event writing_files (t = 0.; t += tsnap)
{
  dump(file = dumpFile);
  sprintf(snapshotFile, "intermediate/snapshot-%5.4f", t);
  dump(file = snapshotFile);
}

/**
## Progress log

The logged kinetic energy is the axisymmetric integral.  A run is stopped
early only for the same two deterministic failure conditions used by the
legacy case: energy blow-up or complete decay after the initial transient.
*/
event log_writing (i++)
{
  double ke = 0.;
  int stop = 0;
  int dump_state = 0;
  foreach (reduction(+:ke))
    ke += (2.*pi*y)*(0.5*rho(f[])*
                     (sq(u.x[]) + sq(u.y[])))*sq(Delta);

  if (pid() == 0) {
    FILE * fp = fopen(logFile, i == 0 ? "w" : "a");
    if (!fp) {
      fprintf(ferr, "ERROR: cannot open %s\n", logFile);
      stop = 1;
    }
    else {
      if (i == 0) {
        fprintf(fp, "CaseNo %d, MAXlevel %d, G1 %g, lambda1 %g\n",
                CaseNo, MAXlevel, G1, lambda1);
        fprintf(fp, "i dt t ke\n");
      }
      fprintf(fp, "%d %.8e %.8e %.8e\n", i, dt, t, ke);
      fclose(fp);
      fprintf(ferr, "%d %.8e %.8e %.8e\n", i, dt, t, ke);

      if (ke > 1e2 && i > 10) {
        fprintf(ferr, "ERROR: kinetic energy blew up.\n");
        stop = 1;
        dump_state = 1;
      }
      if (ke < 1e-8 && i > 10) {
        fprintf(ferr, "Kinetic energy decayed below the stopping threshold.\n");
        stop = 1;
        dump_state = 1;
      }
    }
  }

  mpi_all_reduce(stop, MPI_INT, MPI_MAX);
  mpi_all_reduce(dump_state, MPI_INT, MPI_MAX);
  if (dump_state)
    dump(file = dumpFile);
  assert(ke > -1e-10);
  return stop;
}

/**
## Completion
*/
event stop_simulation (t = tmax)
{
  if (pid() == 0)
    fprintf(ferr, "Case %d complete at t=%g.\n", CaseNo, t);
  return 1;
}

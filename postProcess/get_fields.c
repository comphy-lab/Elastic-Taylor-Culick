/**
# Field snapshot on a regular grid

Restores a single snapshot and interpolates `f`, `u.x`, `u.y`, `KAPPA`, and
`PHI` (the viscous dissipation rate per unit volume, raw -- not log-scaled;
the plotting layer decides the display transform) (plus, when compiled with
`-DVISCOELASTIC=1`, the scalar log-conformation and stress fields
`A11 A12 A22 T11 T12 T22` from
`src-local/log-conform-viscoelastic-scalar-2D.h`) onto a regular
`nx` x `ny` grid over `[xmin,xmax] x [ymin,ymax]`, and writes one row per
grid point to stdout.

`PHI = 2 mu E:E` with `E = sym(grad u)` the strain-rate tensor, i.e.
`2 mu (E_xx^2 + E_yy^2 + 2 E_xy^2)` for planar 2D -- the planar analogue of
the axisymmetric `D11^2+D22^2+D33^2+2*D13^2` invariant in
`comphy-lab/DropsAtLubis`'s `getData.c`, and the same reconstruction
`Taylor-Culick-FEM/postProcess/animate_planar.py` uses for the co-moving FEM
cross-check -- plotted there with `hot_r` and `LogNorm` on the raw value,
which is the convention this column is meant to feed. `mu(f) =
clamp(f,0,1)*(mu1-mu2) + mu2` matches `two-phase-generic.h`'s own macro;
`mu1`/`mu2` are runtime case parameters, not part of the dump, so they are
required arguments here.

A field name absent from the restored dump (e.g. the viscoelastic fields
against a Newtonian-case snapshot, or vice versa) is read back as exactly
zero by Basilisk's `restore()` rather than erroring -- compile the matching
variant for the case a snapshot came from, or the corresponding output
columns are silently all-zero.

Usage: `get_fields snapshot-file xmin xmax ymin ymax ny mu1 mu2`
Columns: `x y f ux uy kappa phi` (+ `A11 A12 A22 T11 T12 T22` if
VISCOELASTIC).
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];
scalar KAPPA[];
scalar PHI[];
#if VISCOELASTIC
scalar A11[], A12[], A22[], T11[], T12[], T22[];
#endif

int main(int argc, char const * argv[])
{
  if (argc < 9) {
    fprintf(stderr,
            "usage: %s snapshot-file xmin xmax ymin ymax ny mu1 mu2\n",
            argv[0]);
    return 1;
  }
  double xmin = atof(argv[2]), xmax = atof(argv[3]);
  double ymin = atof(argv[4]), ymax = atof(argv[5]);
  int ny = atoi(argv[6]);
  double mu1 = atof(argv[7]), mu2 = atof(argv[8]);
  if (!isfinite(xmin) || !isfinite(xmax) || xmax <= xmin ||
      !isfinite(ymin) || !isfinite(ymax) || ymax <= ymin || ny <= 0 ||
      !isfinite(mu1) || !isfinite(mu2)) {
    fprintf(stderr,
            "%s: require finite bounds with xmax > xmin, ymax > ymin, "
            "ny > 0, and finite viscosities\n", argv[0]);
    return 1;
  }
  if (!restore(file = argv[1])) {
    fprintf(stderr, "%s: cannot restore '%s'\n", argv[0], argv[1]);
    return 1;
  }
  f.prolongation = fraction_refine;

  foreach() {
    double exx = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    double eyy = (u.y[0,1] - u.y[0,-1])/(2.*Delta);
    double exy = 0.5*((u.x[0,1] - u.x[0,-1]) + (u.y[1,0] - u.y[-1,0]))/
                 (2.*Delta);
    double mu = clamp(f[], 0., 1.)*(mu1 - mu2) + mu2;
    PHI[] = 2.*mu*(sq(exx) + sq(eyy) + 2.*sq(exy));
  }

  double dy = (ymax - ymin)/ny;
  int nx = (int)((xmax - xmin)/dy);
  if (nx < 1)
    nx = 1;
  double dx = (xmax - xmin)/nx;

  fprintf(stdout, "# nx %d ny %d\n", nx, ny);

  for (int i = 0; i < nx; i++) {
    double x = xmin + dx*(i + 0.5);
    for (int j = 0; j < ny; j++) {
      double y = ymin + dy*(j + 0.5);
      fprintf(stdout, "%g %g %g %g %g %g %g",
              x, y, interpolate(f, x, y),
              interpolate(u.x, x, y), interpolate(u.y, x, y),
              interpolate(KAPPA, x, y), interpolate(PHI, x, y));
#if VISCOELASTIC
      fprintf(stdout, " %g %g %g %g %g %g",
              interpolate(A11, x, y), interpolate(A12, x, y),
              interpolate(A22, x, y), interpolate(T11, x, y),
              interpolate(T12, x, y), interpolate(T22, x, y));
#endif
      fputc('\n', stdout);
    }
  }
  return 0;
}

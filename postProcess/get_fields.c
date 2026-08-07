/**
# Field snapshot on a regular grid

Restores a single snapshot and interpolates `f`, `u.x`, `u.y`, `KAPPA` (plus,
when compiled with `-DVISCOELASTIC=1`, the scalar log-conformation and stress
fields `A11 A12 A22 T11 T12 T22` from
`src-local/log-conform-viscoelastic-scalar-2D.h`) onto a regular
`nx` x `ny` grid over `[xmin,xmax] x [ymin,ymax]`, and writes one row per
grid point to stdout.

A field name absent from the restored dump (e.g. the viscoelastic fields
against a Newtonian-case snapshot, or vice versa) is read back as exactly
zero by Basilisk's `restore()` rather than erroring -- compile the matching
variant for the case a snapshot came from, or the corresponding output
columns are silently all-zero.

Usage: `get_fields snapshot-file xmin xmax ymin ymax ny`
Columns: `x y f ux uy kappa` (+ `A11 A12 A22 T11 T12 T22` if VISCOELASTIC).
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];
scalar KAPPA[];
#if VISCOELASTIC
scalar A11[], A12[], A22[], T11[], T12[], T22[];
#endif

int main(int argc, char const * argv[])
{
  if (argc < 7) {
    fprintf(stderr, "usage: %s snapshot-file xmin xmax ymin ymax ny\n",
            argv[0]);
    return 1;
  }
  double xmin = atof(argv[2]), xmax = atof(argv[3]);
  double ymin = atof(argv[4]), ymax = atof(argv[5]);
  int ny = atoi(argv[6]);
  if (!restore(file = argv[1])) {
    fprintf(stderr, "%s: cannot restore '%s'\n", argv[0], argv[1]);
    return 1;
  }
  f.prolongation = fraction_refine;

  double dy = (ymax - ymin)/ny;
  int nx = (int)((xmax - xmin)/dy);
  double dx = (xmax - xmin)/nx;

  fprintf(stdout, "# nx %d ny %d\n", nx, ny);

  for (int i = 0; i < nx; i++) {
    double x = xmin + dx*(i + 0.5);
    for (int j = 0; j < ny; j++) {
      double y = ymin + dy*(j + 0.5);
      fprintf(stdout, "%g %g %g %g %g %g",
              x, y, interpolate(f, x, y),
              interpolate(u.x, x, y), interpolate(u.y, x, y),
              interpolate(KAPPA, x, y));
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

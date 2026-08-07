/**
# Snapshot-based tip position, mirrored from the in-code diagnostic

Restores a single snapshot and recomputes the same three tip-position
estimates as the `tip_output` event in `TaylorCulickPlanar.c` /
`TaylorCulickPlanarNewtonian.c` (`x_tip`, `x_tip_vof`, `x_tip_global`), plus
the local velocity `u.x` at the `x_tip` cell.  This is deliberately a
*post-hoc, independent* recomputation from the dump alone -- it does not read
the in-code `c<N>-tip.dat` log -- so a mismatch between this and
`tip_to_csv.py`'s output at the same time flags a bug in one of the two
paths rather than being assumed correct by construction.

Usage: `get_tip snapshot-file [h0]` (h0 defaults to 1, matching the case
default; pass the case's actual `h0` for a non-default sheet thickness).
Prints one line: `t x_tip x_tip_vof x_tip_global u_tip_x` to stdout.
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];

int main(int argc, char const * argv[])
{
  if (argc < 2) {
    fprintf(stderr, "usage: %s snapshot-file [h0]\n", argv[0]);
    return 1;
  }
  double h0 = argc > 2 ? atof(argv[2]) : 1.;
  if (!restore(file = argv[1])) {
    fprintf(stderr, "%s: cannot restore '%s'\n", argv[0], argv[1]);
    return 1;
  }
  f.prolongation = fraction_refine;

  double xtip = HUGE, xtipglobal = HUGE, xvof = 0., utip = 0.;
  const double band = h0/10.;

  foreach() {
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);
      coord segment[2];
      if (facets(n, alpha, segment) == 2)
        for (int k = 0; k < 2; k++) {
          double xf = x + segment[k].x*Delta;
          double yf = y + segment[k].y*Delta;
          if (xf < xtipglobal)
            xtipglobal = xf;
          if (yf < band && xf < xtip) {
            xtip = xf;
            utip = u.x[];
          }
        }
    }
    if (y < 0.75*Delta)
      xvof += (1. - clamp(f[], 0., 1.))*Delta;
  }

  fprintf(stdout, "%.8e %.8e %.8e %.8e %.8e\n",
          t, xtip, xvof, xtipglobal, utip);
  return 0;
}

/**
# Interface facets from a snapshot

Restores a single `dump()`-written snapshot and writes the VOF interface of
`f` as line segments (`output_facets()`, `fractions.h`) to stdout: pairs of
`x y` per segment endpoint, blank-line separated, ready for `plot`.

Usage: `get_facets snapshot-file`
*/
#include "navier-stokes/centered.h"
#include "fractions.h"

scalar f[];

int main(int argc, char const * argv[])
{
  if (argc < 2) {
    fprintf(stderr, "usage: %s snapshot-file\n", argv[0]);
    return 1;
  }
  if (!restore(file = argv[1])) {
    fprintf(stderr, "%s: cannot restore '%s'\n", argv[0], argv[1]);
    return 1;
  }
  f.prolongation = fraction_refine;
  output_facets(f, stdout);
  return 0;
}

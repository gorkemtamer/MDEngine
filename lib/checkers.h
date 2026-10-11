#ifndef CHECKERS_H
#define CHECKERS_H

#include "xyzio.h"

// ---------------------------------------------------------------------
// Checks for unphysical geometries. Each check takes a frame only and
// returns 0 if the frame passes, nonzero otherwise. They print nothing
// on a normal failure, so they can be called in loops; the caller
// decides what to report.
// ---------------------------------------------------------------------

// Checks every pair of particles against the minimum distance
// MIN_DIST_FACTOR * (element_radius(i) + element_radius(j)), in reduced
// units (see elements.h / elements.c).
// First makes sure the boundary conditions are defined (f->periodic):
//   1  periodic: distances use the minimum-image convention (nearest copy
//      of j in the periodic images of the box)
//   0  open: plain distances
// Returns  0  all pairs are far enough apart,
//         >0  number of pairs that are too close,
//         -1  invalid input, printed to stderr: f->periodic not 0/1,
//             non-positive box lengths with PBC, or a label without a
//             radius in elements.c.
int CheckMinDistance(const XYZFrame *f);

#endif

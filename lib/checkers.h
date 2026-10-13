#ifndef CHECKERS_H
#define CHECKERS_H

#include "xyzio.h"

// ---------------------------------------------------------------------
// Checks for unphysical geometries. Each check takes a frame only and
// returns 0 if the frame passes, nonzero otherwise. They print nothing
// on a normal failure, so they can be called in loops; the caller
// decides what to report.
// ---------------------------------------------------------------------

// 1 if particle i is closer than the minimum distance
// to any particle j < i, 0 otherwise. Used when placing particles one by one.
int TooCloseToEarlier(const XYZFrame& f, int i);

// Checks every pair of particles against the minimum distance
// MIN_DIST_FACTOR * (sigma_i/2 + sigma_j/2), in reduced
// units (see elements.h / elements.cpp).
// Distances depend on f.boundary():
//   Boundary::Periodic  minimum-image convention (nearest copy of j in the
//                       periodic images of the box)
//   Boundary::Open      plain distances
// The XYZFrame constructor guarantees a valid boundary and positive box.
// Returns  0  all pairs are far enough apart,
//         >0  number of pairs that are too close,

int CheckMinDistance(const XYZFrame& f);

#endif

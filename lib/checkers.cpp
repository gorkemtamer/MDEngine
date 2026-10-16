// checkers.c
//
// Assert-style checks for unphysical properties of a geometry.
// Meant to be run on input/output and, later, at the end of each iteration.

#include "checkers.h"
#include "elements.h"
#include "vec3.h"

//Checks during geometry generation
int TooCloseToEarlier(const XYZFrame& f, int i){
    const double ri = 0.5 * f.element(i).sigma;
    for (int j = 0; j < i; ++j){
        double rj = 0.5 * f.element(j).sigma;
        Vec3 d = f.MinImage(f.pos(j) - f.pos(i));
        double dmin = PLACE_DIST_FACTOR * (ri + rj);
        if (norm2(d) < dmin * dmin) return 1;
    }
    return 0;
}

// Plain O(n^2) loop over all pairs.
// Checks during simulation
// TODO: cell lists if n gets large. TRY TO AVOID THIS
int CheckMinDistance(const XYZFrame& f) {    
    int too_close = 0;
    for (int i = 0; i < f.n(); ++i) {
        double ri = 0.5 * f.element(i).sigma;
        for (int j = i + 1; j < f.n(); ++j) {
            double rj = 0.5 * f.element(j).sigma;
            Vec3 d = f.MinImage(f.pos(i)-f.pos(j));
            double dmin = MIN_DIST_FACTOR * (ri + rj);    // fraction of sigma_ij
            if (norm2(d) < dmin * dmin) ++too_close;   // compare squares, no sqrt
        }
    }
    return too_close;
}

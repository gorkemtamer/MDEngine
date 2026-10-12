// checkers.c
//
// Assert-style checks for unphysical properties of a geometry.
// Meant to be run on input/output and, later, at the end of each iteration.

#include <stdio.h>
#include <math.h>
#include "checkers.h"
#include "elements.h"
#include "vec3.h"

// Plain O(n^2) loop over all pairs.
// TODO: cell lists if n gets large.
int CheckMinDistance(const XYZFrame *f) {
    const Vec3 L = f->box;
    
    int too_close = 0;

    for (int i = 0; i < f->n; ++i) {
        double ri = element_radius(f->label(i));
        if (ri < 0) {
            fprintf(stderr, "CheckMinDistance: no radius for label '%s' (particle %d), add it to elements.cpp\n",
                    f->label(i), i);
            return -1;
        }
        for (int j = i + 1; j < f->n; ++j) {
            double rj = element_radius(f->label(j));
            if (rj < 0) continue;      // reported when the outer loop reaches j
            Vec3 d = Vec3{ f->x[i], f->y[i], f->z[i]} - Vec3{ f->x[j], f->y[j], f->z[j]};
            if (f->boundary == Boundary::Periodic) d = min_image(d, L);
            double dmin = MIN_DIST_FACTOR * (ri + rj);    // fraction of sigma_ij
            if (norm2(d) < dmin * dmin) ++too_close;   // compare squares, no sqrt
        }
    }
    return too_close;
}

// checkers.c
//
// Assert-style checks for unphysical properties of a geometry.
// Meant to be run on input/output and, later, at the end of each iteration.

#include <stdio.h>
#include <math.h>
#include "checkers.h"
#include "elements.h"

// Minimum-image convention: shift d by whole box lengths so that it lies
// in [-L/2, L/2], i.e. the distance to the nearest periodic copy.
static double min_image(double d, double L) {
    return d - L * round(d / L);
}

// Plain O(n^2) loop over all pairs.
// TODO: cell lists if n gets large.
int CheckMinDistance(const XYZFrame *f) {
    // boundary conditions must be defined before any distance makes sense
    if (f->periodic != 0 && f->periodic != 1) {
        fprintf(stderr, "CheckMinDistance: boundary conditions not set (periodic = %d, "
                        "must be 1 for PBC or 0 for OBC)\n", f->periodic);
        return -1;
    }
    const double *L = f->lattice;
    if (f->periodic && (L[0] <= 0 || L[1] <= 0 || L[2] <= 0)) {
        fprintf(stderr, "CheckMinDistance: periodic boundaries need positive box lengths "
                        "(got %g %g %g)\n", L[0], L[1], L[2]);
        return -1;
    }

    int too_close = 0;

    for (int i = 0; i < f->n; ++i) {
        double ri = element_radius(f->label[i]);
        if (ri < 0) {
            fprintf(stderr, "CheckMinDistance: no radius for label '%s' (particle %d), add it to elements.c\n",
                    f->label[i], i);
            return -1;
        }
        for (int j = i + 1; j < f->n; ++j) {
            double rj = element_radius(f->label[j]);
            if (rj < 0) continue;      // reported when the outer loop reaches j
            double dx = f->x[i] - f->x[j];
            double dy = f->y[i] - f->y[j];
            double dz = f->z[i] - f->z[j];
            if (f->periodic) {
                dx = min_image(dx, L[0]);
                dy = min_image(dy, L[1]);
                dz = min_image(dz, L[2]);
            }
            double dmin = MIN_DIST_FACTOR * (ri + rj);    // fraction of sigma_ij
            if (dx*dx + dy*dy + dz*dz < dmin * dmin) ++too_close;   // compare squares, no sqrt
        }
    }
    return too_close;
}

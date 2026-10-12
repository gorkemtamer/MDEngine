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
    // boundary conditions must be defined before any distance makes sense
    if (f->periodic != 0 && f->periodic != 1) {
        fprintf(stderr, "CheckMinDistance: boundary conditions not set (periodic = %d, "
                        "must be 1 for PBC or 0 for OBC)\n", f->periodic);
        return -1;
    }
    const Vec3 L{f->lattice[0], f->lattice[1], f->lattice[2]};
    if (f->periodic && (L.x <= 0 || L.y <= 0 || L.z <= 0)) {
        fprintf(stderr, "CheckMinDistance: periodic boundaries need positive box lengths "
                        "(got %g %g %g)\n", L.x, L.y, L.z);
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
            Vec3 d = Vec3{ f->x[i], f->y[i], f->z[i]} - Vec3{ f->x[j], f->y[j], f->z[j]};
            if (f->periodic) d = min_image(d, L);
            double dmin = MIN_DIST_FACTOR * (ri + rj);    // fraction of sigma_ij
            if (norm2(d) < dmin * dmin) ++too_close;   // compare squares, no sqrt
        }
    }
    return too_close;
}

// xyz_io.c
//
// Definition of XYZ Frame and helper functions to read and write XYZ Frame 
// to a straight CSV file
// XYZ Frame structure also straightforward:
// each column gives the coordinates and velocities of the particles in the 
// system
//
// TODO: later make this a library
//
// Test:   make xyzio_test && bin/xyzio_test   (writes into output/)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "xyzio.h"
#include "checkers.h"
#include "paths.h"
#include "elements.h"
#include "params.h"
#include "vec3.h"

// Allocate a frame with n rows. All columns are allocated; caller fills them.
XYZFrame* xyz_alloc(int n) {
    XYZFrame *f = static_cast<XYZFrame*>(malloc(sizeof(XYZFrame)));
    if (!f) return NULL;

    f->n = n;
    f->comment[0] = '\0';
    f->lattice[0] = f->lattice[1] = f->lattice[2] = 0.0;
    f->periodic = -1;      // caller must choose PBC (1) or OBC (0)
    f->label = static_cast<char (*)[LABEL_LEN]>(malloc(n * sizeof(*f->label)));
    f->x  = static_cast<double*>(malloc(n * sizeof(double)));
    f->y  = static_cast<double*>(malloc(n * sizeof(double)));
    f->z  = static_cast<double*>(malloc(n * sizeof(double)));
    f->vx = static_cast<double*>(malloc(n * sizeof(double)));
    f->vy = static_cast<double*>(malloc(n * sizeof(double)));
    f->vz = static_cast<double*>(malloc(n * sizeof(double)));
    f->fx = static_cast<double*>(calloc(n, sizeof(double)));
    f->fy = static_cast<double*>(calloc(n, sizeof(double)));
    f->fz = static_cast<double*>(calloc(n, sizeof(double)));
    f->epot = 0.0;
    f->forces_valid = 0;

    if (!f->label || !f->x || !f->y || !f->z || !f->vx || !f->vy || !f->vz ||
        !f->fx || !f->fy || !f->fz) {
        fprintf(stderr, "xyz_alloc: allocation failed for n=%d\n", n);
        xyz_free(f); return NULL;
    }
    return f;
}

void xyz_free(XYZFrame *f) {
    if (!f) return;
    free(f->label);
    free(f->x); free(f->y); free(f->z);
    free(f->vx); free(f->vy); free(f->vz);
    free(f->fx); free(f->fy); free(f->fz);
    free(f);
}

// ---------------------------------------------------------------------
// Write an XYZ file:
//   line 1: n
//   line 2: comment
//   line 3: Lattice Lx Ly Lz periodic P     (P = 1: periodic, 0: open)
//   then label x y z vx vy vz, one row per particle.
// ---------------------------------------------------------------------
int xyz_write(const char *path, const XYZFrame *f) {
    int bad = CheckMinDistance(f);
    if (bad != 0) {
        if (bad > 0) fprintf(stderr, "xyz_write: %d particle pair(s) too close, not writing %s\n", bad, path);
        return -1;
    }

    FILE *fp = fopen(path, "w");
    if (!fp) { perror("xyz_write: fopen"); return -1; }

    fprintf(fp, "%d\n", f->n);
    fprintf(fp, "%s\n", f->comment);
    fprintf(fp, "Lattice %.6f %.6f %.6f periodic %d\n",
            f->lattice[0], f->lattice[1], f->lattice[2], f->periodic);

    for (int i = 0; i < f->n; ++i) {
        fprintf(fp, "%s %.6f %.6f %.6f %.6f %.6f %.6f\n",
                f->label[i], f->x[i], f->y[i], f->z[i],
                f->vx[i], f->vy[i], f->vz[i]);
    }

    fclose(fp);
    return 0;
}

// ---------------------------------------------------------------------
// Place n particles uniformly at random in the unit box [0,1)^3, scale the
// coordinates and the box by `scale`, and write the result to `path`.
// Labels are "Ar". Velocities are drawn from the Maxwell-Boltzmann
// distribution at DEFAULT_T (params.h), then the center-of-mass velocity
// is subtracted so the system does not drift as a whole.
// CheckMinDistance is also implemented to avoid overlaps.
// Returns the frame (caller frees it with xyz_free), or NULL on failure.
//
// ---------------------------------------------------------------------
XYZFrame* xyz_random(const char *path, int n, double scale, Rng *rng) {
    if (n <= 0 || scale <= 0) {
        fprintf(stderr, "xyz_random: need n > 0 and scale > 0 (got n=%d, scale=%g)\n", n, scale);
        return NULL;
    }

    XYZFrame *f = xyz_alloc(n);
    if (!f) return NULL;

    snprintf(f->comment, sizeof(f->comment), "random geometry, n=%d, scale=%g", n, scale);
    f->lattice[0] = f->lattice[1] = f->lattice[2] = scale;   // 1 1 1 box, scaled
    f->periodic = 1;      //(now PBC revert back to 0) open boundaries: the box is only the region particles are placed in

    // Rejection sampling: place particle i, check the first i+1 particles
    // (by temporarily setting f->n), redraw i if it is too close to any
    // earlier one. Gives up if the box is too dense.
    const int max_tries = 10000;
    for (int i = 0; i < n; ++i) {
        snprintf(f->label[i], LABEL_LEN, "Ar");
        f->n = i + 1;
        int tries = 0;
        do {
            if (++tries > max_tries) {
                fprintf(stderr, "xyz_random: could not place particle %d after %d tries, "
                                "box too small for n=%d (increase scale)\n", i, max_tries, n);
                f->n = n; xyz_free(f);
                return NULL;
            }
            f->x[i] = scale * rng_uniform(rng);
            f->y[i] = scale * rng_uniform(rng);
            f->z[i] = scale * rng_uniform(rng);
        } while (CheckMinDistance(f) != 0);
    }

    f->n = n;

    // Maxwell-Boltzmann velocities at DEFAULT_T, using each particle's mass
    Vec3 p{}; // total momentum
    double mtot = 0; // total mass
    for (int i = 0; i < n; ++i) {
        double m = element_mass(f->label[i]);
        f->vx[i] = rng_MaxwellBoltzmann(rng, DEFAULT_T, m);
        f->vy[i] = rng_MaxwellBoltzmann(rng, DEFAULT_T, m);
        f->vz[i] = rng_MaxwellBoltzmann(rng, DEFAULT_T, m);
        p += m * Vec3{ f->vx[i], f->vy[i], f->vz[i] };
        mtot += m;
    }

    // remove the center-of-mass velocity: total momentum becomes zero
    const Vec3 vcm = p / mtot;
    for (int i = 0; i < n; ++i) {
        f->vx[i] -= vcm.x;
        f->vy[i] -= vcm.y;
        f->vz[i] -= vcm.z;
    }

    if (xyz_write(path, f) != 0) { xyz_free(f); return NULL; }
    return f;
}

// ---------------------------------------------------------------------
// Read an extended XYZ file into a freshly allocated XYZFrame.
// Returns NULL on failure (bad file, malformed header, etc).
// ---------------------------------------------------------------------
XYZFrame* xyz_read(const char *path) {
    FILE *fp = fopen(path, "r");
    if (!fp) { perror("xyz_read: fopen"); return NULL; }

    char line[512];
    int n;

    // line 1: particle count
    if (!fgets(line, sizeof(line), fp) || sscanf(line, "%d", &n) != 1) {
        fprintf(stderr, "xyz_read: failed to parse particle count\n");
        fclose(fp);
        return NULL;
    }

    XYZFrame *f = xyz_alloc(n);
    if (!f) { fclose(fp); return NULL; }

    // line 2: comment (keep the whole line, minus trailing newline)
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "xyz_read: failed to read comment line\n");
        xyz_free(f); fclose(fp);
        return NULL;
    }
    line[strcspn(line, "\n")] = '\0';
    size_t len = strlen(line);
    if (len >= sizeof(f->comment)) len = sizeof(f->comment) - 1;  // truncate long comments
    memcpy(f->comment, line, len);
    f->comment[len] = '\0';

    // line 3: box lengths and boundary conditions, required
    if (!fgets(line, sizeof(line), fp) ||
        sscanf(line, "Lattice %lf %lf %lf periodic %d",
               &f->lattice[0], &f->lattice[1], &f->lattice[2], &f->periodic) != 4 ||
        f->lattice[0] <= 0 || f->lattice[1] <= 0 || f->lattice[2] <= 0 ||
        (f->periodic != 0 && f->periodic != 1)) {
        fprintf(stderr, "xyz_read: line 3 must be 'Lattice Lx Ly Lz periodic P' "
                        "with positive lengths and P = 0 (open) or 1 (periodic)\n");
        xyz_free(f); fclose(fp);
        return NULL;
    }

    // remaining lines: label x y z vx vy vz
    for (int i = 0; i < n; ++i) {
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "xyz_read: unexpected EOF at row %d\n", i);
            xyz_free(f); fclose(fp);
            return NULL;
        }
        int ok = sscanf(line, "%7s %lf %lf %lf %lf %lf %lf",
                         f->label[i], &f->x[i], &f->y[i], &f->z[i],
                         &f->vx[i], &f->vy[i], &f->vz[i]);
        if (ok != 7) {
            fprintf(stderr, "xyz_read: malformed row %d (parsed %d fields)\n", i, ok);
            xyz_free(f); fclose(fp);
            return NULL;
        }
    }

    fclose(fp);

    int bad = CheckMinDistance(f);
    if (bad != 0) {
        if (bad > 0) fprintf(stderr, "xyz_read: %d particle pair(s) too close in %s\n", bad, path);
        xyz_free(f);
        return NULL;
    }
    return f;
}

#ifdef XYZIO_TEST
// print header, lattice and all rows of a frame
static void print_frame(const XYZFrame *f) {
    printf("Read back %d rows, comment: \"%s\"\n", f->n, f->comment);
    printf("Lattice: %g %g %g, %s boundaries\n", f->lattice[0], f->lattice[1], f->lattice[2],
           f->periodic ? "periodic" : "open");
    printf("%-4s %10s %10s %10s %10s %10s %10s\n",
           "lbl", "x", "y", "z", "vx", "vy", "vz");
    for (int i = 0; i < f->n; ++i) {
        printf("%-4s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               f->label[i], f->x[i], f->y[i], f->z[i],
               f->vx[i], f->vy[i], f->vz[i]);
    }
}

// ---------------------------------------------------------------------
// Test 1: build a small frame in memory, write it, read it back, compare.
// Test 2: generate a random frame, read it back, print it.
// ---------------------------------------------------------------------
int main(void) {
    printf("== Test 1: write / read round-trip ==\n");
    const char *path = OUTPUT_DIR "test.xyz";
    int n = 4;

    XYZFrame *out = xyz_alloc(n);
    snprintf(out->comment, sizeof(out->comment), "test frame, n=%d", n);
    out->lattice[0] = out->lattice[1] = out->lattice[2] = 10.0;
    out->periodic = 1;

    for (int i = 0; i < n; ++i) {
        snprintf(out->label[i], LABEL_LEN, "Ar");
        out->x[i]  = i * 1.1;
        out->y[i]  = i * 2.2;
        out->z[i]  = i * 3.3;
        out->vx[i] = 0.1 * i;
        out->vy[i] = -0.1 * i;
        out->vz[i] = 0.05 * i;
    }

    if (xyz_write(path, out) != 0) { xyz_free(out); return 1; }
    printf("Wrote %s\n", path);

    XYZFrame *in = xyz_read(path);
    if (!in) { xyz_free(out); return 1; }

    print_frame(in);

    // Sanity check: compare with a small tolerance, not exact equality.
    // %.6f only keeps 6 digits after the decimal point, which is not always
    // enough to recover the exact same double bit pattern you started with
    // (e.g. 3*1.1 can come back a few ULPs off after a text round-trip).
    // That's expected and fine for this file format; use a tolerance here,
    // and never compare floating-point values with == in general.
    const double tol = 1e-6;
    int mismatch = 0;
    for (int i = 0; i < n; ++i) {
        if (fabs(out->x[i] - in->x[i]) > tol || fabs(out->vz[i] - in->vz[i]) > tol) mismatch = 1;
    }
    printf(mismatch ? "MISMATCH after round-trip!\n" : "Round-trip OK (within 1e-6 tolerance).\n");

    xyz_free(out);
    xyz_free(in);

    printf("\n== Test 2: random geometry ==\n");
    const char *rpath = OUTPUT_DIR "random_test.xyz";
    Rng rng;
    rng_seed(&rng, 42);      // fixed seed: same geometry every run

    XYZFrame *gen = xyz_random(rpath, 10, 5.0, &rng);
    if (!gen) return 1;
    xyz_free(gen);              // only the file is needed from here on
    printf("Wrote %s\n", rpath);

    XYZFrame *rin = xyz_read(rpath);
    if (!rin) return 1;
    print_frame(rin);
    xyz_free(rin);

    return 0;
}

#endif
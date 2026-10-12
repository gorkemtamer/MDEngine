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

void xyz_free(XYZFrame *f) {
    delete f; //free vectors!!!    
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
    fprintf(fp, "%s\n", f->comment.c_str());
    fprintf(fp, "Lattice %.6f %.6f %.6f periodic %d\n",
            f->box.x, f->box.y, f->box.z, boundary_flag(f->boundary));

    for (int i = 0; i < f->n; ++i) {
        fprintf(fp, "%s %.6f %.6f %.6f %.6f %.6f %.6f\n",
                f->label(i), f->x[i], f->y[i], f->z[i],
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

    // 1 1 1 box, scaled, PBC
    XYZFrame *f = new XYZFrame(n, Vec3{scale, scale, scale},Boundary::Periodic);

    char buf[256];
    snprintf(buf, sizeof(buf),"random geometry, n=%d, scale=%g", n, scale);
    f->comment = buf;
    const int argon = element_index("Ar");

    // Rejection sampling: place particle i, check the first i+1 particles
    // (by temporarily setting f->n), redraw i if it is too close to any
    // earlier one. Gives up if the box is too dense.
    const int max_tries = 10000;
    for (int i = 0; i < n; ++i) {
        f->type[i] = argon;
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
        double m = element_mass(f->label(i));
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

    // line 2: comment (keep the whole line, minus trailing newline)
    if (!fgets(line, sizeof(line), fp)) {
        fprintf(stderr, "xyz_read: failed to read comment line\n");
        fclose(fp);
        return NULL;
    }
    line[strcspn(line, "\n")] = '\0';
    std::string comment = line;

    // line 3: box lengths and boundary conditions, required
    Vec3 box;
    int flag;
    Boundary bc;
    if (!fgets(line, sizeof(line), fp) ||
        sscanf(line, "Lattice %lf %lf %lf periodic %d",
            &box.x, &box.y, &box.z, &flag) != 4 ||
            box.x <= 0 || box.y <= 0 || box.z <= 0 ||
            !boundary_from_flag(flag, bc)) {
        fprintf(stderr, "xyz_read: line 3 must be 'Lattice Lx Ly Lz periodic P' "
                        "with positive lengths and P = 0 (open) or 1 (periodic)\n");
        fclose(fp);
        return NULL;
    }

    XYZFrame *f = new XYZFrame(n, box, bc); //box and flag already checked: cannot throw
    f-> comment = comment;

    // remaining lines: label x y z vx vy vz
    for (int i = 0; i < n; ++i) {
        if (!fgets(line, sizeof(line), fp)) {
            fprintf(stderr, "xyz_read: unexpected EOF at row %d\n", i);
            xyz_free(f); fclose(fp);
            return NULL;
        }
        char lbl[LABEL_LEN];
        int ok = sscanf(line, "%7s %lf %lf %lf %lf %lf %lf",
                         lbl, &f->x[i], &f->y[i], &f->z[i],
                         &f->vx[i], &f->vy[i], &f->vz[i]);
        if (ok != 7) {
            fprintf(stderr, "xyz_read: malformed row %d (parsed %d fields)\n", i, ok);
            xyz_free(f); fclose(fp);
            return NULL;
        }
        f->type[i] = element_index(lbl);
        if (f->type[i]<0)
        {
            fprintf(stderr, "xyz_read: unknown element '%s' in row %d, add it to elements.cpp\n",lbl,i);
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
    printf("Read back %d rows, comment: \"%s\"\n", f->n, f->comment.c_str());
    printf("Lattice: %g %g %g, %s boundaries\n", f->box.x, f->box.y, f->box.z,
           boundary_name(f->boundary));
    printf("%-4s %10s %10s %10s %10s %10s %10s\n",
           "lbl", "x", "y", "z", "vx", "vy", "vz");
    for (int i = 0; i < f->n; ++i) {
        printf("%-4s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               f->label(i), f->x[i], f->y[i], f->z[i],
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

    XYZFrame *out = new XYZFrame(n, Vec3{10.,10.,10.},Boundary::Periodic);
    char buf[64];
    snprintf(buf, sizeof(buf), "test frame, n=%d", n);
    out->comment = buf;

    for (int i = 0; i < n; ++i) {
        out->type[i] = element_index("Ar");
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
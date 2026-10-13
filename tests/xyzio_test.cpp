#include <cmath>
#include <cstdio>
#include <exception>
#include <format>
#include "xyzio.h"
#include "elements.h"
#include "params.h"
#include "paths.h"

// print header, lattice and all rows of a frame
static void print_frame(const XYZFrame& f) {
    printf("Read back %d rows, comment: \"%s\"\n", f.n(), f.comment.c_str());
    printf("Lattice: %g %g %g, %s boundaries\n", f.box().x, f.box().y, f.box().z,
           boundary_name(f.boundary()));
    printf("%-4s %10s %10s %10s %10s %10s %10s\n",
           "lbl", "x", "y", "z", "vx", "vy", "vz");
    for (int i = 0; i < f.n(); ++i) {
        Vec3 r = f.pos(i), v = f.vel(i);
        printf("%-4s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               f.label(i), r.x, r.y, r.z, v.x, v.y, v.z);
    }
}

// ---------------------------------------------------------------------
// Test 1: build a small frame in memory, write it, read it back, compare.
// Test 2: generate a random frame, write it, read it back, print it.
// ---------------------------------------------------------------------
int main(void) {
    try {
        printf("== Test 1: write / read round-trip ==\n");
        const char *path = OUTPUT_DIR "test.xyz";
        int n = 4;

        XYZFrame out(n, Vec3{ 10., 10., 10. }, Boundary::Periodic);
        out.comment = std::format("test frame, n={}", n);
        for (int i = 0; i < n; ++i) {
            out.set_element(i, element_index("Ar"));
            out.set_pos(i, Vec3{ i * 1.1, i * 2.2, i * 3.3 });
            out.set_vel(i, Vec3{ 0.1 * i, -0.1 * i, 0.05 * i });
        }

        xyz_write(path, out);
        printf("Wrote test.xyz\n");

        XYZFrame in = xyz_read(path);
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
            if (std::fabs(out.pos(i).x - in.pos(i).x) > tol ||
                std::fabs(out.vel(i).z - in.vel(i).z) > tol) mismatch = 1;
        }
        printf(mismatch ? "MISMATCH after round-trip!\n" : "Round-trip OK (within 1e-6 tolerance).\n");

        printf("\n== Test 2: random geometry ==\n");
        const char *rpath = OUTPUT_DIR "random_test.xyz";
        Rng rng;
        rng_seed(&rng, 42);      // fixed seed: same geometry every run

        XYZFrame gen = RandomFrame(10, 5.0, Boundary::Periodic, DEFAULT_T, rng);
        xyz_write(rpath, gen);
        printf("Wrote random_test.xyz\n" );

        XYZFrame rin = xyz_read(rpath);
        print_frame(rin);
    } catch (const std::exception& e) {
        fprintf(stderr, "xyzio_test: %s\n", e.what());
        return 1;
    }
    return 0;
}
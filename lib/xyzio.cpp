// Test:   make xyzio_test && bin/xyzio_test   (writes into output/)

#include <cstdio>
#include <cmath>
#include <format>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "xyzio.h"
#include "checkers.h"
#include "paths.h"
#include "elements.h"
#include "params.h"
#include "vec3.h"


namespace {
[[noreturn]] void read_error(const std::string& path, const std::string& what)
{
    throw std::runtime_error("xyz_read: " + path + ": " + what);
}
}

// ---------------------------------------------------------------------
// Read an extended XYZ file into a freshly allocated XYZFrame.
// Throws runtime_error on failure (bad file, malformed header, etc).
// ---------------------------------------------------------------------
XYZFrame xyz_read(const std::string& path){
    std::ifstream in(path);
    if (!in) read_error(path, "cannot open file");

    std::string line;
    //line 1: particle count
    int n;
    if (!std::getline(in, line) || !(std::istringstream(line) >> n) || n < 0)
        read_error(path, "line 1 must be the particle count");

    //line 2: comment, kept as it is
    std::string comment;
    if (!std::getline(in, comment)) read_error(path, "missing comment line");

    //line 3: Lattice Lx Ly Lz periodic P
    Vec3 box;
    int flag;
    Boundary bc;
    std::string kw_lattice, kw_periodic;
    if (!std::getline(in, line) || !(std::istringstream(line) >> kw_lattice >> box.x >> box.y >> box.z >> kw_periodic >> flag) ||
        kw_lattice != "Lattice" || kw_periodic != "periodic" ||
        box.x <= 0 || box.y <= 0 || box.z <= 0 || !boundary_from_flag(flag, bc))
        read_error(path, "line 3 must be 'Lattice Lx Ly Lz periodic P' "
                         "with positive lengths and P = 0 (open) or 1 (periodic)");
    XYZFrame f(n, box, bc);
    f.comment = comment;

    //one row per particle: label x y z vx vy vz
    for (int i = 0; i < n; ++i) {
        if (!std::getline(in, line))
            read_error(path, "file ends after " + std::to_string(i) + " of " + std::to_string(n) + " particles");
        std::string label;
        Vec3 r, v;
        if (!(std::istringstream(line) >> label >> r.x >> r.y >> r.z >> v.x >> v.y >> v.z))
            read_error(path, "malformed row for particle " + std::to_string(i));
        int k = element_index(label.c_str());
        if (k < 0)
            read_error(path, "unknown element '" + label + "' (particle " + std::to_string(i) +
                             "), add it to elements.cpp");
        f.set_element(i, k);
        f.set_pos(i, r);
        f.set_vel(i, v);
    }

    if (int bad = CheckMinDistance(f); bad > 0)
        read_error(path, std::to_string(bad) + " particle pair(s) too close");
    return f;
}

// ---------------------------------------------------------------------
// Write an XYZ file:
//   line 1: n
//   line 2: comment
//   line 3: Lattice Lx Ly Lz periodic P     (P = 1: periodic, 0: open)
//   then label x y z vx vy vz, one row per particle.
// ---------------------------------------------------------------------
void xyz_write(const std::string& path, const XYZFrame& f){
    if (int bad = CheckMinDistance(f); bad > 0)
        throw std::runtime_error("xyz_write: " + path + ": " + std::to_string(bad) +
                                 " particle pair(s) too close, not written");

    std::ofstream out(path);
    if (!out) throw std::runtime_error("xyz_write: cannot open " + path + " for writing");

    const Vec3 L = f.box();
    out << f.n() << '\n'
        << f.comment << '\n'
        << std::format("Lattice {:.6f} {:.6f} {:.6f} periodic {}\n",
                       L.x, L.y, L.z, boundary_flag(f.boundary()));
    for (int i = 0; i < f.n(); ++i) {
        Vec3 r = f.pos(i), v = f.vel(i);
        out << std::format("{} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}\n",
                           f.label(i), r.x, r.y, r.z, v.x, v.y, v.z);
    }
    if (!out) throw std::runtime_error("xyz_write: error while writing " + path);
}

// ---------------------------------------------------------------------
// Place n particles uniformly at random in [0, scale)^3, at least the minimum
// distance apart; Maxwell-Boltzmann velocities at temperature T with the
// center-of-mass velocity removed. Labels are "Ar". Nothing is written to disk.
// Throws std::invalid_argument for n <= 0 or scale <= 0, std::runtime_error
// if the box is too small for n particles.
// ---------------------------------------------------------------------
XYZFrame RandomFrame(int n, double scale, Boundary bc, double T, Rng& rng) {
    if (n <= 0 || scale <= 0)
        throw std::invalid_argument(std::format("RandomFrame: need n > 0 and scale > 0 (got n={}, scale={:g})", n, scale));
    
    // 1 1 1 box, scaled 
    XYZFrame f(n, Vec3{scale,scale,scale},bc);

    f.comment = std::format("random geometry, n={}, scale={:g}", n, scale);
    
    const int argon = element_index("Ar");
    for (int i = 0; i < n; ++i) f.set_element(i,argon);

    // Rejection sampling: place particle i, redraw it while it is too close
    // to any earlier one (TooCloseToEarlier). Gives up if the box is too dense.
    const int max_tries = 10000;
    for (int i = 0; i < n; ++i) {
        int tries = 0;
        do {
            if(++tries > max_tries) 
                throw std::runtime_error(std::format("RandomFrame: could not place particle {} after {} tries, box too small for n={} (increase scale)", i, max_tries, n));
            Vec3 r{ scale*rng_uniform(&rng), scale*rng_uniform(&rng), scale*rng_uniform(&rng)};
            f.set_pos(i, r);
        } while (TooCloseToEarlier(f,i));
    }

    // Maxwell-Boltzmann velocities at temperature T, using each particle's mass
    for (int i = 0; i < n; ++i) {
        double m = f.mass(i);
        Vec3 v{
            rng_MaxwellBoltzmann(&rng, T, m),
            rng_MaxwellBoltzmann(&rng, T, m),
            rng_MaxwellBoltzmann(&rng, T, m)};
        f.set_vel(i,v);
    }

    // remove the center-of-mass velocity: total momentum becomes zero
    const Vec3 vcm = f.TotalMomentum() / f.TotalMass();
    for (int i = 0; i < n; ++i) f.set_vel(i,f.vel(i) - vcm);

    return f;
}

#ifdef XYZIO_TEST
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
        printf("Wrote %s\n", path);

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
        printf("Wrote %s\n", rpath);

        XYZFrame rin = xyz_read(rpath);
        print_frame(rin);
    } catch (const std::exception& e) {
        fprintf(stderr, "xyzio_test: %s\n", e.what());
        return 1;
    }
    return 0;
}

#endif

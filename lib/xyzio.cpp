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
            Vec3 r{ scale*rng.uniform(), scale*rng.uniform(), scale*rng.uniform()};
            f.set_pos(i, r);
        } while (TooCloseToEarlier(f,i));
    }

    // Maxwell-Boltzmann velocities at temperature T, using each particle's mass
    for (int i = 0; i < n; ++i) {
        double m = f.mass(i);
        Vec3 v{
            rng.MaxwellBoltzmann(T, m),
            rng.MaxwellBoltzmann(T, m),
            rng.MaxwellBoltzmann(T, m)};
        f.set_vel(i,v);
    }

    // remove the center-of-mass velocity: total momentum becomes zero
    const Vec3 vcm = f.TotalMomentum() / f.TotalMass();
    for (int i = 0; i < n; ++i) f.set_vel(i,f.vel(i) - vcm);

    return f;
}
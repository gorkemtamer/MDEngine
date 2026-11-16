// wca_wall.cpp
//
// WCA wall forces and energy

#include <format>
#include <stdexcept>
#include "wca_wall.h"

namespace {
constexpr double WCA_CUTOFF = 1.122462048309373;   // 2^(1/6): minimum of LJ, in units of s
}

WCAWall::WCAWall(double eps) : eps_(eps) {
    if (eps <= 0)
        throw std::invalid_argument(std::format("WCAWall: eps must be positive (got {:g})", eps));
}

double WCAWall::AddForces(XYZFrame& f) const {
    if (f.boundary() == Boundary::Periodic)
        throw std::runtime_error("WCAWall: only for open boundaries (the walls are the faces of the box)");

    const Vec3 L = f.box();
    double U = 0.0;
    for (int i = 0; i < f.n(); ++i) {
        const double s = 0.5 * f.element(i).sigma;
        const double rc = WCA_CUTOFF * s;
        const Vec3 r = f.pos(i);

        // one wall at distance z; returns the force along the wall normal
        // (pointing into the box) and adds the energy to U
        auto wall = [&](double z) {
            if (z <= 0)
                throw std::runtime_error(std::format("WCAWall: particle {} is outside the box", i));
            if (z >= rc) return 0.0;
            const double q = s * s / (z * z);
            const double sr6 = q * q * q;
            U += 4 * eps_ * (sr6 * sr6 - sr6) + eps_;
            return 24 * eps_ * (2 * sr6 * sr6 - sr6) / z;
        };
        // low wall pushes towards +, high wall towards -
        const Vec3 F{ wall(r.x) - wall(L.x - r.x),
                      wall(r.y) - wall(L.y - r.y),
                      wall(r.z) - wall(L.z - r.z) };
        f.add_force(i, F);
    }
    return U;
}
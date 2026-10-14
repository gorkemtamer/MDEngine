// harmonic_trap.cpp
//
// Harmonic trap forces and energy

#include <format>
#include <stdexcept>
#include "harmonic_trap.h"

HarmonicTrap::HarmonicTrap(double k, Vec3 center) : k_(k), center_(center) {
    if (k <= 0)
        throw std::invalid_argument(std::format("HarmonicTrap: k must be positive (got {:g})", k));
}

// same k for all particles -maybe create something more general later-
// one potential center: expectation is cluster formation
double HarmonicTrap::AddForces(XYZFrame& f) const {
    if (f.boundary() == Boundary::Periodic)
        throw std::runtime_error("HarmonicTrap: only for open boundaries");

    double U = 0.0;
    for (int i = 0; i < f.n(); ++i) {
        Vec3 d = f.pos(i) - center_;
        U += 0.5 * k_ * norm2(d);
        f.add_force(i, -k_ * d);
    }
    return U;
}
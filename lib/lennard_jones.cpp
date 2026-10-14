// lennard_jones.cpp
//
// Lennard-Jones forces and energy (naive all-pairs)

#include <cmath>
#include <format>
#include <stdexcept>
#include "lennard_jones.h"
#include "elements.h"

LennardJones::LennardJones(double cutoff)
    : cutoff_(cutoff)
    , ntypes_(element_count())
    , table_(element_count() * element_count())
{
    if (cutoff <= 0)
        throw std::invalid_argument(std::format("LennardJones: cutoff must be positive (got {:g})", cutoff));
    
    // energy of the plain LJ potential at the cutoff, per unit 4eps_ij
    const double sc6 = std::pow(1./cutoff_,6);    
    shift_ = sc6 * sc6 - sc6;

    for (int a = 0; a < ntypes_; ++a)
        for (int b = 0; b < ntypes_; ++b){
            Pair& p = table_[a * ntypes_ + b];
            p.sigma = 0.5 * (element(a).sigma + element(b).sigma);   // Lorentz
            p.eps   = std::sqrt(element(a).eps * element(b).eps);    // Berthelot
            double rc = cutoff_ * p.sigma;
            p.rc2 = rc * rc;    
        }
}

double LennardJones::AddForces(XYZFrame& f) const {
    const int n = f.n();

    // minimum image only sees the nearest copy, so the cutoff sphere must fit in half the box
    if (f.boundary() == Boundary::Periodic) {
        double sig_max = 0;
        for (int i = 0; i < n; ++i) sig_max = std::fmax(sig_max, f.element(i).sigma);
        const double rc_max = cutoff_ * sig_max;
        const double L[3] = { f.box().x, f.box().y, f.box().z };
        for (int k = 0; k < 3; ++k)
            if (L[k] < 2.0 * rc_max)
                throw std::runtime_error(std::format(
                    "LennardJones: box length {:g} is smaller than 2 * cutoff = {:g}", L[k], 2.0 * rc_max));
    }
    double U = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Vec3 d = f.MinImage(f.pos(i) - f.pos(j));
            double r2 = norm2(d);

            const Pair& p = pair(f.type(i), f.type(j));
            if (r2 >= p.rc2) continue;
            if (r2 == 0.0)
                throw std::runtime_error(std::format(
                    "LennardJones: particles {} and {} are at the same position", i, j));

            double s2  = p.sigma * p.sigma / r2;
            double s6  = s2 * s2 * s2;
            double s12 = s6 * s6;

            U += 4.0 * p.eps * (s12 - s6 - shift_);

            // F_i = -dV/dr * (r_i - r_j)/r = 24 eps (2 s12 - s6) / r^2 * (r_i - r_j)
            double fr = 24.0 * p.eps * (2.0 * s12 - s6) / r2;
            Vec3 F = fr * d;        // force on i from j
            f.add_force(i, F);
            f.add_force(j, -F);     // Newton's third law
        }
    }
    return U;
}
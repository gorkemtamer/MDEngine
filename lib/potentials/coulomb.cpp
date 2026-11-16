// coulomb.cpp
//
// Coulomb forces and energy (naive all-pairs)

#include <cmath>
#include <format>
#include <stdexcept>
#include <vector>
#include "coulomb.h"
#include "units.h"

Coulomb::Coulomb(double relative_permittivity)
    : eps_r_(relative_permittivity),
      prefactor_(COULOMB_K / relative_permittivity)
{
    if (relative_permittivity <= 0.) throw std::invalid_argument(std::format(
        "Coulomb: relative permittivity must be positive (got {:g})", relative_permittivity));
}

double Coulomb::AddForces(XYZFrame& f) const {
    if (f.boundary() == Boundary::Periodic)
        throw std::runtime_error("Coulomb: only open boundaries, for periodic not implemented yet.");

    //check charges
    const int n = f.n();
    std::vector<double> q(n);
    for (int i = 0; i < n; ++i) q[i] = f.element(i).charge;

    // now compute
    double U = 0.;
    for (int i = 0; i < n; ++i) {
        if (q[i] == 0) continue;    //neutral: no Coulomb term
        for (int j = i + 1; j < n; ++j){
            if (q[j]==0.) continue;
            Vec3 d = f.pos(i) - f.pos(j);
            double r2 = norm2(d);
            if (r2 == 0.0) throw std::runtime_error(std::format(
                    "Coulomb: particles {} and {} are at the same position", i, j));
            double r = std::sqrt(r2);
            double kqq = prefactor_ * q[i] * q[j];

            U += kqq / r;

            //  F_i = -dV/dr * (r_i - r_j)/r = K q_i q_j / r^3 * (r_i - r_j)
            Vec3 F = (kqq / (r2 * r)) * d; //force on i from j
            f.add_force(i, F);
            f.add_force(j, -F); //Newton's third law
        }
    }
    return U;    
}
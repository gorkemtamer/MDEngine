// integrator.c
//
// Lennard-Jones forces (naive all-pairs) and the velocity-Verlet integrator.

#include <cstdio>
#include <cmath>
#include <format>
#include <vector>
#include <stdexcept>
#include "integrator.h"
#include "elements.h"
#include "params.h"
#include "vec3.h"

void ComputeForces(XYZFrame& f) {
    const int n = f.n();
    const double L[3] = {f.box().x, f.box().y, f.box().z};

    // per-particle sigma and epsilon, looked up once instead of once per pair
    std::vector<double> sig(n), eps(n);
    double sig_max = 0;
    for (int i = 0; i < n; ++i) {
        sig[i] = f.element(i).sigma;
        eps[i] = f.element(i).eps;
        if (sig[i] > sig_max) sig_max = sig[i];
    }

    // minimum image only sees the nearest copy, so the cutoff sphere must fit in half the box
    if (f.boundary() == Boundary::Periodic) {
        double rc_max = LJ_CUTOFF * sig_max;
        for (int k = 0; k < 3; ++k)
            if (L[k] < 2.0 * rc_max) throw std::runtime_error(std::format(
                    "ComputeForces: box length {:g} is smaller than 2 * cutoff = {:g}", 
                    L[k], 2.0 * rc_max));            
    }

    f.clear_forces();
    // energy of the plain LJ potential at the cutoff, per unit 4*eps_ij
    const double sc6 = std::pow(1.0 / LJ_CUTOFF, 6);
    const double shift = sc6 * sc6 - sc6;

    double U = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Vec3 d = f.MinImage(f.pos(i) - f.pos(j));
            double r2 = norm2(d);

            double s  = 0.5 * (sig[i] + sig[j]);      // Lorentz
            double e  = std::sqrt(eps[i] * eps[j]);        // Berthelot
            double rc = LJ_CUTOFF * s;
            if (r2 >= rc * rc) continue;
            if (r2 == 0.0) throw std::runtime_error(std::format(
                    "ComputeForces: particles {} and {} are at the same position", i, j));

            double s2  = s * s / r2;
            double s6  = s2 * s2 * s2;
            double s12 = s6 * s6;

            U += 4.0 * e * (s12 - s6 - shift);

            // F_i = -dV/dr * (r_i - r_j)/r = 24 eps (2 s12 - s6) / r^2 * (r_i - r_j)
            double fr = 24.0 * e * (2.0 * s12 - s6) / r2;
            Vec3 F = fr * d; //force on i from j
            f.add_force(i,F);
            f.add_force(j,-F); //Newton's third law
        }
    }
    f.set_potential_energy(U);
}

void VelocityVerletStep(XYZFrame& f, double dt) {
    if (dt <= 0) throw std::invalid_argument(std::format("VelocityVerletStep: time step must be positive (got {:g})", dt));

    // F(t): only needed on the first step, afterwards it is kept from the previous step
    if (!f.forces_valid()) ComputeForces(f);

    // half kick + drift
    for (int i = 0; i < f.n(); ++i) {
        double h = 0.5 * dt / f.mass(i);
        Vec3 v = f.vel(i) + h * f.force(i);
        f.set_vel(i,v);
        f.set_pos(i,f.Wrap(f.pos(i)+dt*v));
    }

    // F(t + dt)
    ComputeForces(f);

    // second half kick
    for (int i = 0; i < f.n(); ++i) {
        double h = 0.5 * dt / f.mass(i);
        f.set_vel(i,f.vel(i)+h*f.force(i));
    }
}

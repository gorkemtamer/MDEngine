// integrator.c
//
// Lennard-Jones forces (naive all-pairs) and the velocity-Verlet integrator.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integrator.h"
#include "elements.h"
#include "params.h"
#include "vec3.h"

int ComputeForces(XYZFrame *f) {
    const int n = f->n;
    const double L[3] = {f -> box.x, f-> box.y, f-> box.z};
    const Vec3 box{ L[0], L[1], L[2]};

    // per-particle sigma and epsilon, looked up once instead of once per pair
    double *sig = static_cast<double*>(malloc(n * sizeof(double)));
    double *eps = static_cast<double*>(malloc(n * sizeof(double)));
    if (!sig || !eps) {
        fprintf(stderr, "ComputeForces: allocation failed for n=%d\n", n);
        free(sig); free(eps);
        return -1;
    }
    double sig_max = 0;
    for (int i = 0; i < n; ++i) {
        double r = element_radius(f->label(i));
        double e = element_epsilon(f->label(i));
        if (r < 0 || e < 0) {
            fprintf(stderr, "ComputeForces: no LJ parameters for label '%s' (particle %d)\n",
                    f->label(i), i);
            free(sig); free(eps);
            return -1;
        }
        sig[i] = 2.0 * r;
        eps[i] = e;
        if (sig[i] > sig_max) sig_max = sig[i];
    }

    // minimum image only sees the nearest copy, so the cutoff sphere must fit in half the box
    if (f-> boundary == Boundary::Periodic) {
        double rc_max = LJ_CUTOFF * sig_max;
        for (int k = 0; k < 3; ++k)
            if (L[k] < 2.0 * rc_max) {
                fprintf(stderr, "ComputeForces: box length %g is smaller than 2 * cutoff = %g\n",
                        L[k], 2.0 * rc_max);
                free(sig); free(eps);
                return -1;
            }
    }

    for (int i = 0; i < n; ++i) f->fx[i] = f->fy[i] = f->fz[i] = 0.0;

    // energy of the plain LJ potential at the cutoff, per unit 4*eps_ij
    const double sc6 = pow(1.0 / LJ_CUTOFF, 6);
    const double shift = sc6 * sc6 - sc6;

    double U = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            Vec3 d = Vec3{ f->x[i], f->y[i], f->z[i] } - Vec3{ f->x[j], f->y[j], f->z[j] };
            if (f-> boundary == Boundary::Periodic) d = min_image(d,box);
            double r2 = norm2(d);

            double s  = 0.5 * (sig[i] + sig[j]);      // Lorentz
            double e  = sqrt(eps[i] * eps[j]);        // Berthelot
            double rc = LJ_CUTOFF * s;
            if (r2 >= rc * rc) continue;
            if (r2 == 0.0) {
                fprintf(stderr, "ComputeForces: particles %d and %d are at the same position\n", i, j);
                free(sig); free(eps);
                return -1;
            }

            double s2  = s * s / r2;
            double s6  = s2 * s2 * s2;
            double s12 = s6 * s6;

            U += 4.0 * e * (s12 - s6 - shift);

            // F_i = -dV/dr * (r_i - r_j)/r = 24 eps (2 s12 - s6) / r^2 * (r_i - r_j)
            double fr = 24.0 * e * (2.0 * s12 - s6) / r2;
            Vec3 F = fr * d; //force on i from j
            f->fx[i] += F.x;  f->fx[j] -= F.x; // Newton's third law:
            f->fy[i] += F.y;  f->fy[j] -= F.y; // force on j is minus
            f->fz[i] += F.z;  f->fz[j] -= F.z; // the force on i
        }
    }

    free(sig); free(eps);
    f->epot = U;
    f->forces_valid = 1;
    return 0;
}

// wrap x into [0, L)
static double wrap(double x, double L) {
    return x - L * floor(x / L);
}

XYZFrame* VelocityVerletStep(XYZFrame *f, double dt) {
    if (dt <= 0) {
        fprintf(stderr, "VelocityVerletStep: time step must be positive (got %g)\n", dt);
        return NULL;
    }

    // F(t): only needed on the first step, afterwards it is kept from the previous step
    if (!f->forces_valid && ComputeForces(f) != 0) return NULL;

    // half kick + drift
    for (int i = 0; i < f->n; ++i) {
        double h = 0.5 * dt / element_mass(f->label(i));
        f->vx[i] += h * f->fx[i];
        f->vy[i] += h * f->fy[i];
        f->vz[i] += h * f->fz[i];
        f->x[i] += dt * f->vx[i];
        f->y[i] += dt * f->vy[i];
        f->z[i] += dt * f->vz[i];
        if (f-> boundary == Boundary::Periodic) {
            f->x[i] = wrap(f->x[i], f-> box.x);
            f->y[i] = wrap(f->y[i], f-> box.y);
            f->z[i] = wrap(f->z[i], f-> box.z);
        }
    }
    f->forces_valid = 0;

    // F(t + dt)
    if (ComputeForces(f) != 0) return NULL;

    // second half kick
    for (int i = 0; i < f->n; ++i) {
        double h = 0.5 * dt / element_mass(f->label(i));
        f->vx[i] += h * f->fx[i];
        f->vy[i] += h * f->fy[i];
        f->vz[i] += h * f->fz[i];
    }
    return f;
}

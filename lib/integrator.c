// integrator.c
//
// Lennard-Jones forces (naive all-pairs) and the velocity-Verlet integrator.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "integrator.h"
#include "elements.h"
#include "params.h"

// same minimum-image convention as in checkers.c
static double min_image(double d, double L) {
    return d - L * round(d / L);
}

int ComputeForces(XYZFrame *f) {
    const int n = f->n;
    const double *L = f->lattice;

    if (f->periodic != 0 && f->periodic != 1) {
        fprintf(stderr, "ComputeForces: boundary conditions not set (periodic = %d)\n", f->periodic);
        return -1;
    }

    // per-particle sigma and epsilon, looked up once instead of once per pair
    double *sig = malloc(n * sizeof(double));
    double *eps = malloc(n * sizeof(double));
    if (!sig || !eps) {
        fprintf(stderr, "ComputeForces: allocation failed for n=%d\n", n);
        free(sig); free(eps);
        return -1;
    }
    double sig_max = 0;
    for (int i = 0; i < n; ++i) {
        double r = element_radius(f->label[i]);
        double e = element_epsilon(f->label[i]);
        if (r < 0 || e < 0) {
            fprintf(stderr, "ComputeForces: no LJ parameters for label '%s' (particle %d)\n",
                    f->label[i], i);
            free(sig); free(eps);
            return -1;
        }
        sig[i] = 2.0 * r;
        eps[i] = e;
        if (sig[i] > sig_max) sig_max = sig[i];
    }

    // minimum image only sees the nearest copy, so the cutoff sphere must fit in half the box
    if (f->periodic) {
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
            double dx = f->x[i] - f->x[j];
            double dy = f->y[i] - f->y[j];
            double dz = f->z[i] - f->z[j];
            if (f->periodic) {
                dx = min_image(dx, L[0]);
                dy = min_image(dy, L[1]);
                dz = min_image(dz, L[2]);
            }
            double r2 = dx*dx + dy*dy + dz*dz;

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
            f->fx[i] += fr * dx;  f->fx[j] -= fr * dx;   // Newton's third law:
            f->fy[i] += fr * dy;  f->fy[j] -= fr * dy;   // force on j is minus
            f->fz[i] += fr * dz;  f->fz[j] -= fr * dz;   // the force on i
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
        double h = 0.5 * dt / element_mass(f->label[i]);
        f->vx[i] += h * f->fx[i];
        f->vy[i] += h * f->fy[i];
        f->vz[i] += h * f->fz[i];
        f->x[i] += dt * f->vx[i];
        f->y[i] += dt * f->vy[i];
        f->z[i] += dt * f->vz[i];
        if (f->periodic) {
            f->x[i] = wrap(f->x[i], f->lattice[0]);
            f->y[i] = wrap(f->y[i], f->lattice[1]);
            f->z[i] = wrap(f->z[i], f->lattice[2]);
        }
    }
    f->forces_valid = 0;

    // F(t + dt)
    if (ComputeForces(f) != 0) return NULL;

    // second half kick
    for (int i = 0; i < f->n; ++i) {
        double h = 0.5 * dt / element_mass(f->label[i]);
        f->vx[i] += h * f->fx[i];
        f->vy[i] += h * f->fy[i];
        f->vz[i] += h * f->fz[i];
    }
    return f;
}

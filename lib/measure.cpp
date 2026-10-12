// measure.c
//
// Observables of a single frame.

#include <math.h>
#include "measure.h"
#include "integrator.h"
#include "elements.h"

double KineticEnergy(const XYZFrame *f) {
    double ke = 0.0;
    for (int i = 0; i < f->n; ++i)
        ke += 0.5 * element_mass(f->label[i]) *
              (f->vx[i]*f->vx[i] + f->vy[i]*f->vy[i] + f->vz[i]*f->vz[i]);
    return ke;
}

double PotentialEnergy(XYZFrame *f) {
    // one LJ implementation for forces and energy: reuse ComputeForces
    if (!f->forces_valid && ComputeForces(f) != 0) return NAN;
    return f->epot;
}

// measure.c
//
// Observables of a single frame.

#include <math.h>
#include "measure.h"
#include "integrator.h"
#include "elements.h"

double PotentialEnergy(XYZFrame *f) {
    // one LJ implementation for forces and energy: reuse ComputeForces
    if (!f->forces_valid() && ComputeForces(f) != 0) return NAN;
    return f->PotentialEnergy();
}

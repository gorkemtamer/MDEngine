// measure.c
//
// Observables of a single frame.

#include <cmath>
#include "measure.h"
#include "integrator.h"
#include "elements.h"

double PotentialEnergy(XYZFrame& f) {
    // one LJ implementation for forces and energy: reuse ComputeForces
    if (!f.forces_valid()) ComputeForces(f) ;
    return f.PotentialEnergy();
}

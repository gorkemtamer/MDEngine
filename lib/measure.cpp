// measure.cpp
//
// Observables of a single frame.

#include "measure.h"

double PotentialEnergy(XYZFrame& f, const ForceField& ff) {
    // forces and energy come from the same computation: reuse it
    if (!f.forces_valid()) ff.Compute(f);
    return f.PotentialEnergy();
}

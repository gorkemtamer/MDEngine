#ifndef MEASURE_H
#define MEASURE_H

#include "xyzio.h"
#include "forcefield.h"

// ---------------------------------------------------------------------
// Observables measured on a single frame, in reduced units
// (see reduced_units.md).
// ---------------------------------------------------------------------


// Total potential energy U of the frame from all potentials of the force
// field ff (the same ones that give the forces).
// Uses f.PotentialEnergy() if the forces are up to date (f.forces_valid()),
// which is always the case right after a VelocityVerletStep; otherwise it calls
// ff.Compute first, so f's forces are filled as a side effect.
// Errors of the force field propagate.
double PotentialEnergy(XYZFrame& f, const ForceField& ff);

#endif

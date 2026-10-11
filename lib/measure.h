#ifndef MEASURE_H
#define MEASURE_H

#include "xyzio.h"

// ---------------------------------------------------------------------
// Observables measured on a single frame, in reduced units
// (see reduced_units.md).
// ---------------------------------------------------------------------

// Total kinetic energy  K = sum_i 1/2 m_i v_i^2.
double KineticEnergy(const XYZFrame *f);

// Total potential energy U of the frame (Lennard-Jones, same potential as
// the forces: see integrator.h).
// Uses f->epot if the forces are up to date (f->forces_valid), which is
// always the case right after a VelocityVerletStep; otherwise it calls
// ComputeForces first, so f's forces are filled as a side effect.
// Returns NAN on error (message printed by ComputeForces).
double PotentialEnergy(XYZFrame *f);

#endif

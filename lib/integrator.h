#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "xyzio.h"
#include "forcefield.h"

// ---------------------------------------------------------------------
// Time integration, in reduced units (see reduced_units.md).
// The forces come from a ForceField (forcefield.h).
// ---------------------------------------------------------------------

// Advances the frame by one velocity-Verlet step of length dt:
//     v(t + dt/2) = v(t)        + dt/2 * F(t) / m
//     r(t + dt)   = r(t)        + dt   * v(t + dt/2)
//     F(t + dt)   = forces at the new positions
//     v(t + dt)   = v(t + dt/2) + dt/2 * F(t + dt) / m
// F(t) is reused from the previous step when f.forces_valid() is true,
// so each step costs one force evaluation. With periodic boundaries the
// positions are wrapped back into the box.
// Updates f in place.
// Throws std::invalid_argument for dt<=0; errors from the force field propagate.
void VelocityVerletStep(XYZFrame& f, const ForceField& ff, double dt);

#endif

#ifdef INTEGRATOR_TEST



#endif
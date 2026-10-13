#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "xyzio.h"

// ---------------------------------------------------------------------
// Forces and time integration, in reduced units (see reduced_units.md).
//
// Potential: Lennard-Jones between every pair, truncated at
// LJ_CUTOFF * sigma_ij (params.h) and shifted so the energy is zero at
// the cutoff. Unlike pairs use the Lorentz-Berthelot rules
//     sigma_ij = (sigma_i + sigma_j) / 2,   eps_ij = sqrt(eps_i * eps_j).
// Periodic frames use the minimum-image convention.
// ---------------------------------------------------------------------

// Computes all forces and the total potential energy
// stores them in the frame (add_force, set_potential_energy)
// Naive O(n^2) loop over pairs j > i; each pair force is added to i and
// subtracted from j (Newton's third law), so every pair is computed once.
// throws std::runtime_error when the box is narrower than
// 2xcutoff or two particles coincide
void ComputeForces(XYZFrame& f); 

// Advances the frame by one velocity-Verlet step of length dt:
//     v(t + dt/2) = v(t)        + dt/2 * F(t) / m
//     r(t + dt)   = r(t)        + dt   * v(t + dt/2)
//     F(t + dt)   = forces at the new positions
//     v(t + dt)   = v(t + dt/2) + dt/2 * F(t + dt) / m
// F(t) is reused from the previous step when f.forces_valid() is true,
// so each step costs one force evaluation. With periodic boundaries the
// positions are wrapped back into the box.
// Updates f in place.
// Throws std::invalid_argument for dt<=0;errors from ComputeForces propogate.
void VelocityVerletStep(XYZFrame& f, double dt);

#endif

#ifdef INTEGRATOR_TEST



#endif
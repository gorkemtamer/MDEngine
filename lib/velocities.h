#pragma once
#include "frame.h"
#include "rng.h"

// ---------------------------------------------------------------------
// Initial velocities, in reduced units (see reduced_units.md).
// Used by the frame builders (RandomFrame, LatticeFrame) and after an
// energy minimization, which leaves the velocities untouched.
// ---------------------------------------------------------------------

// Maxwell-Boltzmann velocities at temperature T: each component of v_i is
// Gaussian with variance T / m_i (m_i from the particle's element), drawn
// in the order x, y, z of particle 0, 1, ... from rng. The center-of-mass
// velocity is removed afterwards, so <K> = (3N - 3) T / 2.
// The temperature is not rescaled to exactly T: K fluctuates as in the
// canonical ensemble (a thermostat brings it to T).
// Positions and forces are not touched. T = 0 gives zero velocities.
// Throws std::invalid_argument for T < 0 or a particle without an element.
void AssignVelocities(XYZFrame& f, double T, Rng& rng);

// Subtracts the center-of-mass velocity from every particle: the total
// momentum becomes zero. Nothing happens for an empty frame.
void RemoveCenterOfMassVelocity(XYZFrame& f);


// Removes the rigid rotation about the center of mass R: with
// L = sum m_i (r_i - R) x v_i and the inertia tensor I about R,
//     omega = I^-1 L,     v_i -= omega x (r_i - R)
// afterwards L = 0. The total momentum is not changed (call
// RemoveCenterOfMassVelocity as well). Both together leave N_f = 3N - 6
// degrees of freedom, <K> = (3N - 6) T / 2 after AssignVelocities.
// Linear configurations (all particles on one line, e.g. two particles)
// can only rotate about the perpendicular axes: omega = L / sum m_i |r_i - R|^2.
// Only meaningful when L is conserved: open boundaries without walls
// (a HarmonicTrap conserves L only about its center).
// Nothing happens for fewer than two particles.
// Throws std::runtime_error for periodic frames (no angular momentum).
void RemoveAngularMomentum(XYZFrame& f);
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
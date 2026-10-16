#pragma once
#include "potential.h"

// ---------------------------------------------------------------------
// WCA walls: the six faces of the box (x = 0, x = Lx, y = 0, ...) repel
// the particles with the purely repulsive part of Lennard-Jones. For a
// particle at distance z from a wall:
//     U(z) = 4 eps [(s/z)^12 - (s/z)^6] + eps   for z < 2^(1/6) s,  else 0
// U and F go smoothly to zero at the cutoff (WCA: shifted at the LJ minimum).
// s = sigma_i / 2: Lorentz rule with a point-like wall (sigma_wall = 0), so
// the particle's surface touches the wall at z ~ sigma_i / 2.
// A one-body potential like HarmonicTrap: it closes an open system (the box
// becomes a container) and breaks momentum conservation on purpose.
// Only for open boundaries: in a periodic box the faces are not walls.
// ---------------------------------------------------------------------

class WCAWall : public Potential {
public:
    // eps: wall strength (reduced units: eps_Ar); throws std::invalid_argument for eps <= 0
    explicit WCAWall(double eps = 1.0);

    // throws std::runtime_error for periodic frames and for particles
    // outside the box (on or beyond a wall: time step too large?)
    double AddForces(XYZFrame& f) const override;

    double eps() const { return eps_; }

private:
    double eps_;
};
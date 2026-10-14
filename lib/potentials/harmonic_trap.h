#pragma once
#include "potential.h"

// ---------------------------------------------------------------------
// Harmonic trap: every particle is pulled towards a fixed center c,
//     V = sum_i 1/2 k |r_i - c|^2,      F_i = -k (r_i - c)
// A one-body potential (no pairs): it keeps an open cluster from
// evaporating. It breaks momentum conservation on purpose (the trap is
// an external force). Each particle oscillates with omega = sqrt(k / m_i).
// Only for open boundaries: in a periodic box a particle leaving one side
// would jump to the other, and with it its trap energy.
// ---------------------------------------------------------------------

class HarmonicTrap : public Potential {
public:
    // k: sprin constant (reduced units: eps / sigma^2),
    // center: trap center
    // throws std::invalid_argument for k<=0
    HarmonicTrap(double k, Vec3 center);

    // throws std::runtime_error for periodic frames
    double AddForces(XYZFrame& f) const override;

    double k() const {return k_;}
    Vec3 center() const {return center_;}

private:
    double k_;
    Vec3 center_;
};
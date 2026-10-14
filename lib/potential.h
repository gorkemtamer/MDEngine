#pragma once
#include "frame.h"

// ---------------------------------------------------------------------
// One contribution to the interactions of a system. 
// Several potentials can act on the same frame;
// a ForceField holds them and adds up their contributions.
// ---------------------------------------------------------------------

class Potential {
public:
    virtual ~Potential() = default;
    // Adds this potential's forces to f (f.add_force) and returns its
    // potential energy. Must not clear the forces already in f: other
    // potentials contribute to the same frame (ForceField clears them).
    virtual double AddForces(XYZFrame& f) const = 0;
};

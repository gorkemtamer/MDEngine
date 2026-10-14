#pragma once
#include <vector>
#include "potential.h"
#include "params.h"

// ---------------------------------------------------------------------
// Lennard-Jones between every pair of particles:
//     V(r) = 4 eps_ij [ (sigma_ij/r)^12 - (sigma_ij/r)^6 ] - V_shift
// truncated at cutoff * sigma_ij and shifted so the energy is zero at the
// cutoff (the force still jumps there). Unlike pairs use the
// Lorentz-Berthelot rules
//     sigma_ij = (sigma_i + sigma_j) / 2,   eps_ij = sqrt(eps_i * eps_j),
// with sigma_i, eps_i from the element table (elements.h).
// Periodic frames use the minimum-image convention.
// Naive O(n^2) loop over pairs j > i; each pair force is added to i and
// subtracted from j (Newton's third law).
// ---------------------------------------------------------------------
class LennardJones : public Potential {
public:
    // cutoff in units of sigma_ij.
    // Throws std::invalid_argument for cutoff <= 0.
    explicit LennardJones(double cutoff = LJ_CUTOFF);

    // Throws std::runtime_error if a periodic box is narrower
    // than 2 * cutoff, or two particles are at the same position.
    double AddForces(XYZFrame& f) const override;

    double cutoff() const {return cutoff_;}

private:
    // parameters of one element pair, computed once in the constructor
    struct Pair {
        double sigma; // sigma_ij
        double eps; // eps_ij
        double rc2; // (cutoff * sigma_ij)^2
    };

    const Pair& pair(int a, int b) const {return table_[a * ntypes_ + b];}

    double cutoff_; 
    double shift_;  // LJ energy at the cutoff per 4eps_ij
    int ntypes_;  //#elements in element table  
    std::vector<Pair> table_; // ntypes_ x ntypes_, indexed by element index
};
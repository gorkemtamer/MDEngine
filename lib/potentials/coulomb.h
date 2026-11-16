#pragma once
#include "potential.h"

// ---------------------------------------------------------------------
// Coulomb interaction between charged particles:
//     V = sum_{i<j} K q_i q_j / (eps_r r_ij)
// with charges q_i from the element table (in e), r in sigma_Ar and
// K = COULOMB_K (units.h, ~394). eps_r is the relative permittivity of a
// surrounding medium (1 = vacuum; e.g. ~78 for water screens by 1/78).
// Naive: every pair, no cutoff (1/r decays too slowly to cut off), O(n^2).
// Pairs with a neutral particle are skipped, so neutral atoms cost nothing.
// Only for open boundaries: with periodic copies the sum over all images
// needs Ewald summation (not implemented).
// ---------------------------------------------------------------------

class Coulomb : public Potential {
public:
    //Throws std::invalid_argument for eps_r<=0
    explicit Coulomb(double relative_premittivity = 1.);

    //Throw std::runtime_error for periodic frames, 
    //or if two particles coincide.
    double AddForces(XYZFrame& f) const override;

    double relative_permittivity() const { return eps_r_;}

private:
    double eps_r_;
    double prefactor_; //COULOMB_K / eps_r
};
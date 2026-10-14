#pragma once
#include <numbers>

// ---------------------------------------------------------------------
// The reduced (Lennard-Jones) unit system, with argon as the reference
// particle: sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1, k_B = 1.
// Charges are in units of the elementary charge e.
// See reduced_units.md for all derived units and conversions.
// ---------------------------------------------------------------------

// reference values = 1 reduced unit
inline constexpr double SIGMA_REF = 3.418;     // sigma of Ar in Angstrom  -> length unit
inline constexpr double EPS_REF   = 124.0;     // epsilon/k_B of Ar in K   -> energy unit
inline constexpr double MASS_REF  = 39.948;    // mass of Ar in u          -> mass unit

// SI constants (exact in the 2019 SI, except epsilon_0: CODATA 2018)
inline constexpr double ELEMENTARY_CHARGE   = 1.602176634e-19;   // C
inline constexpr double VACUUM_PERMITTIVITY = 8.8541878128e-12;  // F/m
inline constexpr double BOLTZMANN           = 1.380649e-23;      // J/K
inline constexpr double ANGSTROM            = 1e-10;             // m

// Coulomb constant in reduced units: two charges q_i, q_j (in e) at
// distance r (in sigma_Ar) have the energy  U = COULOMB_K * q_i q_j / r
// (in epsilon_Ar).  COULOMB_K = e^2 / (4 pi eps_0 sigma_Ar eps_Ar) ~ 394.3,
// i.e. two unit charges at 1 sigma interact ~400 times stronger than two
// argon atoms at their LJ minimum.
inline constexpr double COULOMB_K =
    ELEMENTARY_CHARGE * ELEMENTARY_CHARGE /
    (4.0 * std::numbers::pi * VACUUM_PERMITTIVITY * (SIGMA_REF * ANGSTROM) * (EPS_REF * BOLTZMANN));

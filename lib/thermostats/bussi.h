#pragma once
#include "thermostat.h"
#include "rng.h"

// ---------------------------------------------------------------------
// Stochastic velocity rescaling (Bussi, Donadio, Parrinello,
// J. Chem. Phys. 126, 014101 (2007)). All velocities are scaled by one
// factor alpha, chosen so that the kinetic energy K relaxes to its target
// N_f T / 2 with time constant tau and fluctuates as in the canonical
// ensemble. With c = exp(-dt/tau) and N_f Gaussian numbers R_1..R_Nf:
//     K_new = K + (1-c) (K_target sum R_i^2 / N_f - K)
//               + 2 R_1 sqrt(c (1-c) K K_target / N_f)
//     alpha = sqrt(K_new / K)
// tau -> infinity: plain NVE.  tau = 0: new canonical K every step.
// Typical: tau ~ 0.1 (some 20-100 steps).
//
// Degrees of freedom N_f = 3N - removed_dof:
//   removed_dof = 3 (default): total momentum conserved and zero (free or
//                 interacting particles without external forces)
//   removed_dof = 0: external forces (WCAWall, HarmonicTrap) change the
//                 total momentum, so all 3N velocity components are free
// rng must live at least as long as the thermostat.
// ---------------------------------------------------------------------
class Bussi : public Thermostat {
public:
    // Throws std::invalid_argument for T <= 0 or tau < 0.
    Bussi(double T, double tau, Rng& rng, int removed_dof = 3);
    // Forbids a temporary generator, e.g. Bussi(T, tau, Rng(1)).
    Bussi(double, double, Rng&&, int = 3) = delete;

    void Apply(XYZFrame& f, double dt) override;

    double tau() const { return tau_; }
    int removed_dof() const {return removed_dof_;}

private:
    double tau_;
    Rng& rng_;
    int removed_dof_;
};
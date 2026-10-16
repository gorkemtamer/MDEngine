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
// N_f = 3N - 3: the frames start with zero total momentum, and a common
// scaling factor keeps it zero.
// The thermostat does not own the random generator, it only refers to it:
// rng must live at least as long as the thermostat.
// ---------------------------------------------------------------------
class Bussi : public Thermostat {
public:
    // Throws std::invalid_argument for T <= 0 or tau < 0.
    Bussi(double T, double tau, Rng& rng);
    // Forbids a temporary generator, e.g. Bussi(T, tau, Rng(1)).
    Bussi(double, double, Rng&&) = delete;

    void Apply(XYZFrame& f, double dt) override;

    double tau() const { return tau_; }

private:
    double tau_;
    Rng& rng_;
};
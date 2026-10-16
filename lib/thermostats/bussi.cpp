// bussi.cpp
#include <cmath>
#include <format>
#include <stdexcept>
#include "thermostats/bussi.h"

Bussi::Bussi(double T, double tau, Rng& rng) : Thermostat(T), tau_(tau), rng_(rng) {
    if (tau < 0) throw std::invalid_argument(std::format("Bussi: tau must be >= 0 (got {:g})", tau));
}

void Bussi::Apply(XYZFrame& f, double dt) {
    const int nf = 3 * f.n() - 3;          // degrees of freedom, momentum removed
    if (nf <= 0) throw std::invalid_argument("Bussi: needs at least 2 particles");
    const double K = f.KineticEnergy();
    if (K <= 0) throw std::runtime_error("Bussi: kinetic energy is zero, nothing to rescale");

    const double K_target = 0.5 * nf * T_;
    const double c = (tau_ > 0) ? std::exp(-dt / tau_) : 0.0;

    // N_f Gaussian numbers: R_1 separately, plus the sum of all squares
    const double r1 = rng_.normal();
    double sum = r1 * r1;
    for (int i = 1; i < nf; ++i) {
        double r = rng_.normal();
        sum += r * r;
    }

    // K_new >= 0 always (it is a sum of squares, see the paper's appendix)
    const double K_new = K + (1 - c) * (K_target * sum / nf - K)
                       + 2 * r1 * std::sqrt(c * (1 - c) * K * K_target / nf);
    const double alpha = std::sqrt(K_new / K);
    for (int i = 0; i < f.n(); ++i) f.set_vel(i, alpha * f.vel(i));

    energy_added_ += K_new - K;
}
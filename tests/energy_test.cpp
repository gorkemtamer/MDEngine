// energy_test.cpp
//
// Energy conservation of the velocity-Verlet integrator.
// The same start (lattice, T = 1) is run for the same short time twice,
// with time steps dt and dt/2:
//   1. the energy stays close to its start value (largest deviation small)
//   2. halving dt makes the RMS deviation about 4 times smaller: velocity
//      Verlet is a second-order method, its energy error scales like dt^2.
//      A first-order (wrong) integrator gives only about 2.
// Why a short run and the RMS: the potential is shifted in energy but the
// force jumps at the cutoff. Every pair crossing the cutoff adds a small
// error that does not scale like dt^2. Over long runs these add up and make
// the ratio noisy; over 0.5 tau the ratio stays within 3.4 - 5.3 for 50
// different starts (checked when writing this test).

#include <cmath>
#include <exception>
#include "integrator.h"
#include "measure.h"
#include "xyzio.h"
#include "test_helpers.h"

struct EnergyError {
    double max;   // largest |E(t) - E(0)| per particle
    double rms;   // root mean square of E(t) - E(0) per particle
};

// Runs the lattice start for time t_end with step dt.
static EnergyError energy_error(double dt, double t_end) {
    Rng rng(5);                                   // same start for every dt
    XYZFrame f = LatticeFrame(5, 1.2, 1.0, rng);
    const double E0 = f.KineticEnergy() + PotentialEnergy(f);
    const int steps = static_cast<int>(std::lround(t_end / dt));

    double max = 0, sum2 = 0;
    for (int s = 0; s < steps; ++s) {
        VelocityVerletStep(f, dt);
        double d = f.KineticEnergy() + PotentialEnergy(f) - E0;
        max = std::fmax(max, std::fabs(d));
        sum2 += d * d;
    }
    return { max / f.n(), std::sqrt(sum2 / steps) / f.n() };
}

int main(){
    try {
        const double t_end = 0.5;
        const EnergyError e1 = energy_error(0.005, t_end);
        const EnergyError e2 = energy_error(0.0025, t_end);
        const double ratio = e1.rms / e2.rms;

        check(e1.max < 1e-3, "dt = 0.005:  max |E - E0| per particle", e1.max);
        check(e2.max < 1e-3, "dt = 0.0025: max |E - E0| per particle", e2.max);
        check(ratio > 3.0 && ratio < 6.0, "RMS error ratio dt / (dt/2), expected ~4", ratio);
    } catch(const std::exception& e) {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
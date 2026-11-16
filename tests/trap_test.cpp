// trap_test.cpp
//
// One particle in a harmonic trap is a harmonic oscillator with a known
// exact solution. Released at rest at distance A from the center:
//     x(t) = A cos(omega t),   omega = sqrt(k / m)
// This tests the trap and the integrator together, against an exact answer:
//   1. after one full period the particle is back at x = A
//   2. the position error over the run shrinks ~4x when dt is halved
//      (second-order integrator)
//   3. a periodic frame is rejected (the trap is only for open boundaries)

#include <cmath>
#include <exception>
#include <memory>
#include <numbers>
#include "integrators/velocity_verlet.h"
#include "xyzio.h"
#include "potentials/harmonic_trap.h"
#include "test_helpers.h"

// Largest |x(t) - A cos(omega t)| over one period, for time step dt.
static double max_position_error(double dt, double k, double A) {
    ForceField ff;
    ff.Add(std::make_unique<HarmonicTrap>(k, Vec3{ 5, 5, 5 }));
    XYZFrame f(1, Vec3{ 10, 10, 10 }, Boundary::Open);
    f.set_element(0, element_index("Kr"));          // m != 1, so the mass matters
    f.set_pos(0, Vec3{ 5 + A, 5, 5 });
    f.set_vel(0, Vec3{ 0, 0, 0 });

    const double omega = std::sqrt(k / f.mass(0));
    const double period = 2 * std::numbers::pi / omega;
    const int steps = static_cast<int>(std::lround(period / dt));

    VelocityVerlet vv(ff, dt);
    double err = 0;
    for (int s = 1; s <= steps; ++s) {
        vv.Step(f);
        double exact = A * std::cos(omega * s * dt);
        err = max_nan(err, std::fabs((f.pos(0).x - 5) - exact));
    }
    return err;
}

int main() {
    try {
        const double k = 3.0, A = 1.0;
        const double e1 = max_position_error(0.01, k, A);
        const double e2 = max_position_error(0.005, k, A);

        check(e1 < 1e-3, "dt = 0.01:  max |x - A cos(wt)| over one period", e1);
        check(e2 < 1e-3, "dt = 0.005: max |x - A cos(wt)| over one period", e2);
        check(e1 / e2 > 3.5 && e1 / e2 < 4.5, "error ratio dt / (dt/2), expected 4", e1 / e2);

        // periodic frames are rejected
        ForceField ff;
        ff.Add(std::make_unique<HarmonicTrap>(k, Vec3{ 5, 5, 5 }));
        XYZFrame p(1, Vec3{ 10, 10, 10 }, Boundary::Periodic);
        p.set_element(0, element_index("Ar"));
        bool thrown = false;
        try { ff.Compute(p); } catch (const std::runtime_error&) { thrown = true; }
        check(thrown, "periodic frame is rejected");
    } catch (const std::exception& e) {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
// reversibility_test.cpp
//
// Velocity Verlet is time-reversible: run N steps, flip every velocity,
// run N steps again, and the system retraces its path back to the start.
// At the end every particle must be at its start position, with its start
// velocity reversed. Checked for a periodic and an open frame.
//
// In exact arithmetic the return is perfect. In floating point, the
// rounding errors (~1e-16 per step) are amplified by the chaotic dynamics:
// for the dense periodic lattice roughly x100 every 250 steps. After 500
// steps they are ~1e-11, after 1000 ~1e-7, and after 2000 the reversal
// fails completely (checked when writing this test). So the run is kept
// at 500 steps with a tolerance of 1e-8. A non-reversible integrator
// misses the start by many orders of magnitude more.

#include <cmath>
#include <exception>
#include <vector>
#include "integrators/velocity_verlet.h"
#include "xyzio.h"
#include "test_helpers.h"

struct ReturnError {
    double pos;   // largest distance from the start position
    double vel;   // largest |v_end + v_start|
};

// N steps forward, flip velocities, N steps forward.
static ReturnError run_and_return(XYZFrame& f, const ForceField& ff, double dt, int steps) {
    VelocityVerlet vv(ff, dt);
    const int n = f.n();
    std::vector<Vec3> r0(n), v0(n);
    for (int i = 0; i < n; ++i) {
        r0[i] = f.pos(i);
        v0[i] = f.vel(i);
    }

    for (int s = 0; s < steps; ++s) vv.Step(f);
    for (int i = 0; i < n; ++i) f.set_vel(i, -f.vel(i));
    for (int s = 0; s < steps; ++s) vv.Step(f);

    ReturnError e{ 0, 0 };
    for (int i = 0; i < n; ++i) {
        // MinImage: a particle wrapped across a periodic box is still at its start
        e.pos = std::fmax(e.pos, norm(f.MinImage(f.pos(i) - r0[i])));
        e.vel = std::fmax(e.vel, norm(f.vel(i) + v0[i]));
    }
    return e;
}

int main() {
    try {
        const double dt = 0.005;
        const int steps = 500;
        const double tol = 1e-8;
        const ForceField ff = LJForceField();
        Rng rng(7);

        XYZFrame p = LatticeFrame(5, 1.2, 1.0, rng);
        ReturnError e = run_and_return(p, ff, dt, steps);
        check(e.pos < tol, "periodic: max distance from start position", e.pos);
        check(e.vel < tol, "periodic: max |v_end + v_start|", e.vel);

        XYZFrame o = RandomFrame(30, 6.0, Boundary::Open, 1.0, rng);
        e = run_and_return(o, ff, dt, steps);
        check(e.pos < tol, "open: max distance from start position", e.pos);
        check(e.vel < tol, "open: max |v_end + v_start|", e.vel);
    } catch (const std::exception& e) {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
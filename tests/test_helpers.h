#pragma once

/*
    Shared helpers for physics tests.
    Each test program calls check() for every condition
    it tests and ends with return test_result().
    It is either 0 meaning all checks passed, or 1 otherwise.
    Ctest only looks at the exit code
*/

#include <cstdio>
#include <stdexcept>
#include "frame.h"
#include "rng.h"

// #failed checks in this test program
inline int g_failures = 0;

// Print one line per check: PASS/FAIL, what was checked
// and the measured value for information
inline void check(bool ok, const char *what, double value){
    std::printf("%s %-50s %.6e\n", ok ? "PASS" : "FAIL", what, value);
    if (!ok) ++g_failures;
}

// Same, for checks without a number (e.g. "sequence matches")
inline void check(bool ok, const char *what){
    std::printf("%s %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++g_failures;
}

// Summary line and exit code for main()
inline int test_result(){
    if (g_failures == 0) {
        std::printf("all checks passed\n");
        return 0;
    }
    std::printf("%d check(s) FAILED\n", g_failures);
    return 1;
}

// A calm starting geometry: k*k*k particles of element `label` on a simple
// cubic lattice with spacing a (in sigma), in a periodic box of length k*a.
// Maxwell-Boltzmann velocities at temperature T, center-of-mass velocity
// removed. Much gentler than RandomFrame, so energy tests can be strict.
inline XYZFrame LatticeFrame(int k, double a, double T, Rng& rng, const char *label = "Ar") {
    const int el = element_index(label);
    if (el < 0) throw std::invalid_argument("LatticeFrame: unknown element");

    XYZFrame f(k * k * k, Vec3{ k * a, k * a, k * a }, Boundary::Periodic);
    int i = 0;
    for (int ix = 0; ix < k; ++ix)
        for (int iy = 0; iy < k; ++iy)
            for (int iz = 0; iz < k; ++iz, ++i) {
                f.set_element(i, el);
                f.set_pos(i, Vec3{ (ix + 0.5) * a, (iy + 0.5) * a, (iz + 0.5) * a });
            }

    for (i = 0; i < f.n(); ++i) {
        double m = f.mass(i);
        Vec3 v{ rng.MaxwellBoltzmann(T, m), rng.MaxwellBoltzmann(T, m), rng.MaxwellBoltzmann(T, m) };
        f.set_vel(i, v);
    }
    const Vec3 vcm = f.TotalMomentum() / f.TotalMass();
    for (i = 0; i < f.n(); ++i) f.set_vel(i, f.vel(i) - vcm);
    return f;
}




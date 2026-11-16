#pragma once

/*
    Shared helpers for physics tests.
    Each test program calls check() for every condition
    it tests and ends with return test_result().
    It is either 0 meaning all checks passed, or 1 otherwise.
    Ctest only looks at the exit code
*/

#include <cstdio>
#include <cmath>
#include <memory>
#include "xyzio.h"   // LatticeFrame
#include "forcefield.h"
#include "potentials/lennard_jones.h"

// #failed checks in this test program
inline int g_failures = 0;

// Print one line per check: PASS/FAIL, what was checked
// and the measured value for information.
// A NaN value always fails, whatever ok says.
inline void check(bool ok, const char *what, double value){
    if (std::isnan(value)) ok = false;
    std::printf("%s %-50s %.6e\n", ok ? "PASS" : "FAIL", what, value);
    if (!ok) ++g_failures;
}

// Same, for checks without a number (e.g. "sequence matches")
inline void check(bool ok, const char *what){
    std::printf("%s %s\n", ok ? "PASS" : "FAIL", what);
    if (!ok) ++g_failures;
}

// Larger of a and b for running maxima of errors, like std::fmax, but NaN
// if either is NaN. std::fmax(err, NaN) returns err: a NaN error (e.g. from
// 0/0) would vanish from the maximum and "err < tol" would pass.
// Use this instead of std::fmax in the tests.
inline double max_nan(double a, double b){
    return (std::isnan(a) || std::isnan(b)) ? NAN : std::fmax(a, b);
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

// The interactions used by the physics tests: Lennard-Jones only, default cutoff.
inline ForceField LJForceField() {
    ForceField ff;
    ff.Add(std::make_unique<LennardJones>());
    return ff;
}

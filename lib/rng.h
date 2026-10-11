#ifndef RNG_H
#define RNG_H

#include <stdint.h>

// ---------------------------------------------------------------------
// xoshiro256** pseudo-random number generator (Blackman & Vigna, 2018).
// Fast, good statistical quality, period 2^256 - 1.
// Not cryptographically secure (not needed for simulations).
//
// The state lives in a struct, so you can have independent streams
// (e.g. one per thread) and the same seed always gives the same sequence.
// ---------------------------------------------------------------------
typedef struct {
    uint64_t s[4];
    int    has_spare;   // rng_normal makes normals in pairs; the second one
    double spare;       // is kept here for the next call
} Rng;

// Initialize the state from a single 64-bit seed (any value, 0 included).
void rng_seed(Rng *r, uint64_t seed);

// Next raw 64-bit random number.
uint64_t rng_u64(Rng *r);

// Uniform double in [0, 1): 0 possible, 1 never.
double rng_uniform(Rng *r);

// Standard normal (Gaussian) number: mean 0, variance 1.
double rng_normal(Rng *r);

// One Cartesian velocity component from the Maxwell-Boltzmann distribution:
// Gaussian with mean 0 and variance kT / mass (reduced units, k_B = 1).
// Call it three times for vx, vy, vz.
double rng_MaxwellBoltzmann(Rng *r, double kT, double mass);

#endif

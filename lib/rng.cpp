// rng.c
//
// xoshiro256** generator, adapted from the public-domain reference
// implementation: https://prng.di.unimi.it/xoshiro256starstar.c
// Seeding uses splitmix64, as recommended by the authors, so that even
// simple seeds (0, 1, 2, ...) give well-mixed, nonzero starting states.

#include <cmath>
#include "rng.h"

static inline uint64_t rotl(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

// splitmix64: only used to expand one 64-bit seed into the 256-bit state
static uint64_t splitmix64(uint64_t *x) {
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void rng_seed(Rng *r, uint64_t seed) {
    for (int i = 0; i < 4; ++i) r->s[i] = splitmix64(&seed);
    r->has_spare = 0;   // a reseeded stream must not reuse an old normal
}

uint64_t rng_u64(Rng *r) {
    uint64_t *s = r->s;
    const uint64_t result = rotl(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl(s[3], 45);

    return result;
}

// Top 53 bits -> double in [0, 1). 53 bits = full precision of a double,
// so every value is k / 2^53 with k = 0 .. 2^53 - 1.
double rng_uniform(Rng *r) {
    return (rng_u64(r) >> 11) * 0x1.0p-53;
}

// Box-Muller transform: two uniforms -> two independent standard normals,
//   R = sqrt(-2 ln u1),  z0 = R cos(2 pi u2),  z1 = R sin(2 pi u2).
// z1 is stored and returned by the next call.
double rng_normal(Rng *r) {
    if (r->has_spare) {
        r->has_spare = 0;
        return r->spare;
    }
    const double two_pi = 6.283185307179586;
    double u1 = 1.0 - rng_uniform(r);     // (0, 1]: avoids log(0)
    double u2 = rng_uniform(r);
    double R = std::sqrt(-2.0 * std::log(u1));

    r->spare = R * std::sin(two_pi * u2);
    r->has_spare = 1;
    return R * std::cos(two_pi * u2);
}

double rng_MaxwellBoltzmann(Rng *r, double kT, double mass) {
    return std::sqrt(kT / mass) * rng_normal(r);
}

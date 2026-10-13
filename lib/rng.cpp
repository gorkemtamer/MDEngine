// rng.cpp
//
// xoshiro256** generator, adapted from the public-domain reference
// implementation: https://prng.di.unimi.it/xoshiro256starstar.c
// Seeding uses splitmix64, as recommended by the authors, so that even
// simple seeds (0, 1, 2, ...) give well-mixed, nonzero starting states.

#include <cmath>
#include "rng.h"

namespace {

inline std::uint64_t rot1(std::uint64_t x, int k){
    return (x << k) | (x >> (64 - k));    
}

//splitmix64: only used to expand one 64-bit seed into the
//256-bit state
std::uint64_t splitmix64(std::uint64_t& x){
    std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

} //namespace

void Rng::seed(std::uint64_t s) {
    for (int i = 0; i < 4; ++i) s_[i] = splitmix64(s);
    has_spare_ = false;   // a reseeded stream must not reuse an old normal
}

std::uint64_t Rng::u64() {
    const std::uint64_t result = rot1(s_[1] * 5, 7) * 9;
    const std::uint64_t t = s_[1] << 17;

    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rot1(s_[3], 45);

    return result;
}

// Top 53 bits -> double in [0, 1). 53 bits = full precision of a double,
// so every value is k / 2^53 with k = 0 .. 2^53 - 1.
double Rng::uniform() {
    return (u64() >> 11) * 0x1.0p-53;
}

// Box-Muller transform: two uniforms -> two independent standard normals,
//   R = sqrt(-2 ln u1),  z0 = R cos(2 pi u2),  z1 = R sin(2 pi u2).
// z1 is stored and returned by the next call.
double Rng::normal() {
    if (has_spare_) {
        has_spare_ = false;
        return spare_;
    }
    const double two_pi = 6.283185307179586;
    double u1 = 1.0 - uniform();     // (0, 1]: avoids log(0)
    double u2 = uniform();
    double R = std::sqrt(-2.0 * std::log(u1));

    spare_ = R * std::sin(two_pi * u2);
    has_spare_ = true;
    return R * std::cos(two_pi * u2);
}

double Rng::MaxwellBoltzmann(double kT, double mass) {
    return std::sqrt(kT / mass) * normal();
}
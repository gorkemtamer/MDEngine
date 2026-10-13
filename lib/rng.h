#ifndef RNG_H
#define RNG_H

#include <cstdint>

// ---------------------------------------------------------------------
// xoshiro256** pseudo-random number generator (Blackman & Vigna, 2018).
// Fast, good statistical quality, period 2^256 - 1.
// Not cryptographically secure (not needed for simulations).
//
// Each Rng object is an independent stream (e.g. one thread);
// the same seed always gives the same sequence
// Rng rng(42);
// double u = rng.uniform();
// ---------------------------------------------------------------------
class Rng{
public:
    //Seeded from a single 64-bit value (any value, 0 included)
    explicit Rng(std::uint64_t seed) {this->seed(seed);}

    //A copy would continue the same sequence: 
    //Two "independent" streams that are in fact identical.
    //Pass generators by reference instead.
    Rng(const Rng&) = delete;
    Rng& operator=(const Rng&) = delete;

    //Restart the sequence from a new seed
    void seed(std::uint64_t s);

    //Next raw 64-bit random number
    std::uint64_t u64();

    //Uniform double in [0,1): 0 possible, 1 never.
    double uniform();
    //Standard normal (Gaussian) number: mean 0, variance 1.
    double normal();
    //One cartesian velocity component from the Maxwell-Boltzmann distribution:
    //Gaussian with mean 0 and variance kT/mass (reduced units, kB=1)
    //Call it three times for vx,vy,vz
    double MaxwellBoltzmann(double kT, double mass);
private:
    std::uint64_t s_[4];
    bool has_spare_ = false; //normal() makes normals in pairs; the second one
    double spare_ = 0.0; //is kept here for the next call
};

#endif

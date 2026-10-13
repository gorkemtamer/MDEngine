/*
    rng_test.cpp
    Tests for the random number generator
    1. xoshiro256** reference sequence, state {1,2,3,4}
    as produced by the published reference implementation 
    (prng.di.unimi.it)
    2. uniform(): every value in [0,1), mean 1/2
    3. normal(): mean 0, variance 1
    4. MaxwellBoltzmann(kT,m): variance kT/m
    5. same seed -> same sequence;
    seed() restarts it, including the stored second Box-Muller value
    Seeds are fixed to ensure reproducibility
*/

#include <array>
#include <cmath>
#include <cstdint>
#include <exception>
#include "rng.h"
#include "test_helpers.h"

int main(){
    try{
        //1. reference sequence
        Rng ref(std::array<std::uint64_t, 4>{ 1, 2, 3, 4 });
        const std::uint64_t expected[] = { 11520ULL, 0ULL, 1509978240ULL, 1215971899390074240ULL };
        bool same = true;
        for (std::uint64_t e : expected) same = same && (ref.u64() == e);
        check(same, "xoshiro256** reference sequence");

        const int N = 1000000;
        Rng rng(2026);
        
        //2.uniform:range and mean (std error: sqrt(1/12 / N))
        // remember monte carlo sampling
        double sum = 0, lo=1, hi=0;
        for (int i = 0; i < N; i++){
            double u =rng.uniform();
            sum += u;
            lo = std::fmin(lo,u);
            hi = std::fmax(hi,u);
        }
        check(lo>=0.0 && hi < 1., "uniform(): all values in [0,1)");
        check(std::fabs(sum / N - 0.5) < 5 * std::sqrt(1./12 / N),
        "uniform(): mean - 1/2", sum / N - 0.5);

        //3.normal: mean (std error 1/sqrt(N)) and variance (std error sqrt(2/N))
        double s1 = 0, s2 = 0;
        for (int i = 0; i < N; ++i){
            double z = rng.normal();
            s1 += z;
            s2 += z * z;
        }
        double mean = s1/N , var = s2/N - mean * mean;
        check(std::fabs(mean) < 5 / std::sqrt(N), "normal(): mean", mean);
        check(std::fabs(var - 1.0) < 5 * std::sqrt(2.0 / N), "normal(): variance - 1", var - 1.0);

        //4.Maxwell-Boltzmann component: variance kT/m
        const double kT = 1.5, m = 2.0;
        s2 = 0;
        for (int i = 0; i < N; ++i) {
            double v = rng.MaxwellBoltzmann(kT, m);
            s2 += v * v;
        }
        check(std::fabs(s2 / N / (kT / m) - 1.0) < 5 * std::sqrt(2.0 / N),
              "MaxwellBoltzmann(): variance / (kT/m) - 1", s2 / N / (kT / m) - 1.0);

        // 5. reproducibility
        Rng a(7), b(7);
        same = true;
        for (int i = 0; i < 1000; ++i) same = same && (a.uniform() == b.uniform());
        check(same, "same seed gives the same sequence");

        Rng c(11), fresh(11);
        c.normal();                 // leaves a stored second value behind
        c.seed(11);                 // must discard it
        check(c.normal() == fresh.normal(), "seed() restarts the sequence (incl. Box-Muller spare)");
    } catch (const std::exception& e){
        std::printf("FAIL exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
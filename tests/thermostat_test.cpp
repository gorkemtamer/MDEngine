// thermostat_test.cpp
//
// The Bussi thermostat must sample the canonical ensemble. There the
// kinetic energy K has a known distribution, whatever the potential:
// mean N_f T / 2 and variance N_f T^2 / 2  (N_f = 3N - 3).
//   1. free particles (empty force field): mean and variance of K.
//      Only 8 particles, so N_f = 21 differs clearly from 3N = 24.
//   2. Lennard-Jones lattice: mean of K, and H = K + U - energy_added()
//      is conserved (the thermostat's bookkeeping is right)
//   3. free particles in a walled box: N_f = 3N (removed_dof = 0), since
//      the walls change the total momentum
//   4. bad parameters are rejected

#include <cmath>
#include <exception>
#include <stdexcept>
#include "xyzio.h"
#include "rng.h"
#include "integrators/velocity_verlet.h"
#include "thermostats/bussi.h"
#include "potentials/wca_wall.h"
#include "test_helpers.h"

struct KStats { double mean, var, drift; };

// Equilibrates for n_eq steps, then measures K over n_steps steps, and the
// largest deviation of H = K + U - energy_added() per particle.
static KStats sample(XYZFrame& f, const ForceField& ff, Bussi& th, double dt, int n_eq, int n_steps) {
    VelocityVerlet vv(ff, dt);
    for (int s = 0; s < n_eq; ++s) { vv.Step(f); th.Apply(f, dt); }

    const double H0 = f.KineticEnergy() + f.PotentialEnergy(ff) - th.energy_added();
    double sum = 0, sum2 = 0, drift = 0;
    for (int s = 0; s < n_steps; ++s) {
        vv.Step(f);
        th.Apply(f, dt);
        const double K = f.KineticEnergy();
        sum += K;
        sum2 += K * K;
        drift = std::fmax(drift, std::fabs(K + f.PotentialEnergy(ff) - th.energy_added() - H0));
    }
    const double mean = sum / n_steps;
    return { mean, sum2 / n_steps - mean * mean, drift / f.n() };
}

int main() {
    try {
        const double T = 1.5, tau = 0.1, dt = 0.005;
        Rng rng(7);

        // 1. free particles: only the thermostat changes K, many cheap steps
        XYZFrame f = LatticeFrame(2, 1.2, 0.5, rng);   // starts too cold (T = 0.5)
        double nf = 3.0 * f.n() - 3.0;
        ForceField none;
        Bussi th(T, tau, rng);
        KStats s = sample(f, none, th, dt, 1000, 200000);
        check(std::fabs(s.mean / (0.5 * nf * T) - 1) < 0.02, "free: <K> / (N_f T/2), expected 1", s.mean / (0.5 * nf * T));
        check(std::fabs(s.var / (0.5 * nf * T * T) - 1) < 0.1, "free: var(K) / (N_f T^2/2), expected 1", s.var / (0.5 * nf * T * T));

        // 2. Lennard-Jones lattice
        XYZFrame g = LatticeFrame(5, 1.2, 0.5, rng);
        nf = 3.0 * g.n() - 3.0;
        const ForceField lj = LJForceField();
        Bussi th_lj(T, tau, rng);
        s = sample(g, lj, th_lj, dt, 1000, 4000);
        check(std::fabs(s.mean / (0.5 * nf * T) - 1) < 0.03, "LJ: <K> / (N_f T/2), expected 1", s.mean / (0.5 * nf * T));
        check(s.drift < 1e-2, "LJ: max |H - H0| per particle, H = K + U - added", s.drift);

        // 3. free particles in a walled box: the walls change the total
        //    momentum, so N_f = 3N (removed_dof = 0)
        XYZFrame b = RandomFrame(8, 3.0, Boundary::Open, 0.5, rng, 0.6);
        nf = 3.0 * b.n();
        ForceField walls;
        walls.Add(std::make_unique<WCAWall>());
        Bussi th_box(T, tau, rng, 0);
        s = sample(b, walls, th_box, dt, 1000, 200000);
        check(std::fabs(s.mean / (0.5 * nf * T) - 1) < 0.02, "walls: <K> / (3N T/2), expected 1", s.mean / (0.5 * nf * T));
        check(std::fabs(s.var / (0.5 * nf * T * T) - 1) < 0.1, "walls: var(K) / (3N T^2/2), expected 1", s.var / (0.5 * nf * T * T));

        // 4. bad parameters
        bool thrown = false;
        try { Bussi bad(0.0, tau, rng); } catch (const std::invalid_argument&) { thrown = true; }
        check(thrown, "T = 0 is rejected");
        thrown = false;
        try { Bussi bad(T, -1.0, rng); } catch (const std::invalid_argument&) { thrown = true; }
        check(thrown, "tau < 0 is rejected");
        thrown = false;
        try { Bussi bad(T, tau, rng, -1); } catch (const std::invalid_argument&) { thrown = true; }
        check(thrown, "removed_dof < 0 is rejected");
    } catch (const std::exception& e) {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
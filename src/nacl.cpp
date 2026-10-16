#include <cstdio>
#include <filesystem>
#include <format>
#include <stdexcept>
#include "xyzio.h"
#include "rng.h"
#include "paths.h"
#include "params.h"
#include "integrator.h"
#include "checkers.h"
#include "measure.h"
#include "forcefield.h"
#include "elements.h"
#include "potentials/lennard_jones.h"
#include "potentials/coulomb.h"

/* Same as RunMD in main.cpp, but writes and prints only every `every`-th
   step: Na-Cl needs a small time step, so many steps.
   TODO: merge with main.cpp's RunMD (Phase 6, VelocityVerlet class). */
static void RunMD(XYZFrame& geom, const ForceField& ff, const std::string& dir,
                  double dt, int n_steps, int every) {
    std::filesystem::create_directories(dir);
    xyz_write(dir + "initialconfig.xyz", geom);

    double ekin = geom.KineticEnergy();
    double epot = PotentialEnergy(geom, ff);
    printf("%6s %8s %12s %12s %12s\n", "step", "t", "Ekin", "Epot", "Etot");
    printf("%6d %8.3f %12.4f %12.4f %12.4f\n", 0, 0.0, ekin, epot, ekin + epot);

    for (int step = 1; step <= n_steps; ++step) {
        VelocityVerletStep(geom, ff, dt);
        if (CheckMinDistance(geom) != 0) throw std::runtime_error("minDistance failed at step " + std::to_string(step));
        if (step % every != 0) continue;

        geom.comment = std::format("MD step {}, t = {:g}", step, step * dt);
        xyz_write(dir + std::format("step{:05d}.xyz", step), geom);
        ekin = geom.KineticEnergy();
        epot = geom.PotentialEnergy();
        printf("%6d %8.3f %12.4f %12.4f %12.4f\n", step, step * dt, ekin, epot, ekin + epot);
    }
}

/* Rock-salt (NaCl) crystal of nx*ny*nz ions with nearest-neighbour distance a,
   open boundaries: a simple cubic lattice whose sites alternate between Na
   and Cl ((ix+iy+iz) even: Na). Maxwell-Boltzmann velocities at T with the
   center-of-mass velocity removed, as in LatticeFrame. */
static XYZFrame RockSaltFrame(int nx, int ny, int nz, double a, double T, Rng& rng) {
    const int na = element_index("Na"), cl = element_index("Cl");
    XYZFrame f(nx * ny * nz, Vec3{ nx * a, ny * a, nz * a }, Boundary::Open);
    f.comment = std::format("rock salt {}x{}x{}, a={:g}", nx, ny, nz, a);
    int i = 0;
    for (int ix = 0; ix < nx; ++ix)
        for (int iy = 0; iy < ny; ++iy)
            for (int iz = 0; iz < nz; ++iz, ++i) {
                f.set_element(i, (ix + iy + iz) % 2 == 0 ? na : cl);
                f.set_pos(i, Vec3{ (ix + 0.5) * a, (iy + 0.5) * a, (iz + 0.5) * a });
            }

    for (i = 0; i < f.n(); ++i) {
        double m = f.mass(i);
        f.set_vel(i, Vec3{ rng.MaxwellBoltzmann(T, m), rng.MaxwellBoltzmann(T, m), rng.MaxwellBoltzmann(T, m) });
    }
    const Vec3 vcm = f.TotalMomentum() / f.TotalMass();
    for (i = 0; i < f.n(); ++i) f.set_vel(i, f.vel(i) - vcm);
    return f;
}

int main() {
    try {
    /* interactions: LJ (short-range repulsion) + Coulomb, open boundaries only */
    ForceField ff;
    ff.Add(std::make_unique<LennardJones>());
    ff.Add(std::make_unique<Coulomb>());

    const double dt      = 0.001;  /* Na-Cl vibration period ~0.036 tau, see reduced_units.md */
    const int    n_steps = 5000;   /* 5 tau */
    const int    every   = 10;     /* write/print every 10th step: 500 frames */

    /* 1. random: 100 ions in a 10x10x10 box, Na and Cl in turn (neutral),
          frames in output/NaClRandom/ */
    Rng random_rng(42);  /* fixed seed */
    XYZFrame geom = RandomFrame(100, 10.0, Boundary::Open, DEFAULT_T, random_rng);
    const int na = element_index("Na"), cl = element_index("Cl");
    for (int i = 0; i < geom.n(); ++i) geom.set_element(i, i % 2 == 0 ? na : cl);
    printf("NaClRandom: %d ions, frames in output/NaClRandom/\n", geom.n());
    RunMD(geom, ff, std::string(OUTPUT_DIR) + "NaClRandom/", dt, n_steps, every);
    printf("--------------------------------\n");

    /* 2. rock-salt crystal 5x5x4 = 100 ions (50 Na, 50 Cl), a = 0.8 sigma
          (~2.7 A, real NaCl: 2.82 A), frames in output/NaClLattice/ */
    Rng lattice_rng(1);  /* fixed seed */
    XYZFrame crystal = RockSaltFrame(5, 5, 4, 0.8, DEFAULT_T, lattice_rng);
    printf("NaClLattice: %d ions, frames in output/NaClLattice/\n", crystal.n());
    RunMD(crystal, ff, std::string(OUTPUT_DIR) + "NaClLattice/", dt, n_steps, every);
    printf("--------------------------------\n");
    return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
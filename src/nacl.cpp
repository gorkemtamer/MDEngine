#include <cstdio>
#include <format>
#include <stdexcept>
#include <string>
#include "xyzio.h"
#include "rng.h"
#include "paths.h"
#include "params.h"
#include "elements.h"
#include "velocities.h"
#include "potentials/lennard_jones.h"
#include "potentials/coulomb.h"
#include "potentials/wca_wall.h"
#include "simulation.h"

/* LJ (short-range repulsion) + Coulomb, open boundaries only */
static ForceField ff_cluster() {
    ForceField ff;
    ff.Add(std::make_unique<LennardJones>());
    ff.Add(std::make_unique<Coulomb>());
    return ff;
}

/* Full pipeline for one start geometry, frames and energies in
   output/<name>/equilibration/ and output/<name>/production/:
   minimize, Maxwell-Boltzmann velocities at T = 1, NVT (Bussi, tau = 0.1),
   then NVE. Na-Cl needs a small time step (vibration period ~0.036 tau),
   so frames are written every 10th step and energies printed every 500th.
   removed_dof: 6 for a free cluster, 0 inside walls (see simulation.h). */
static void Run(const std::string& name, XYZFrame start, ForceField ff, int removed_dof) {
    Simulation sim(std::move(start), std::move(ff), 11, removed_dof);

    const MinimizeResult m = sim.Minimize({ .force_tol = 1e-2, .max_steps = 5000 });
    printf("%s: %d ions, minimized in %d steps: U = %.4f, max |F| = %.2e%s\n", name.c_str(), sim.frame().n(),
           m.steps, m.energy, m.max_force, m.converged ? "" : " (not converged)");

    sim.AssignVelocities(DEFAULT_T);
    const std::string dir = std::string(OUTPUT_DIR) + name + "/";
    printf("equilibration (NVT), frames in %sequilibration/\n", dir.c_str());
    sim.Equilibrate({ .dt = 0.001, .steps = 10000, .ensemble = Ensemble::NVT, .T = DEFAULT_T, .tau = 0.1,
                      .write_every = 10, .dir = dir + "equilibration/", .print_every = 500 });
    printf("production (NVE), frames in %sproduction/\n", dir.c_str());
    sim.Produce({ .dt = 0.001, .steps = 3000, .write_every = 10, .dir = dir + "production/", .print_every = 500 });
    printf("--------------------------------\n");
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

    AssignVelocities(f, T, rng);
    return f;
}

int main() {
    try {

    /* 1. random: 100 ions in a 10x10x10 box, Na and Cl in turn (neutral) */
    Rng random_rng(42);  /* fixed seed */
    XYZFrame geom = RandomFrame(100, 10.0, Boundary::Open, DEFAULT_T, random_rng);
    const int na = element_index("Na"), cl = element_index("Cl");
    for (int i = 0; i < geom.n(); ++i) geom.set_element(i, i % 2 == 0 ? na : cl);
    Run("NaClRandom", std::move(geom), ff_cluster(), 6);

    /* 2. rock-salt crystal 5x5x4 = 100 ions (50 Na, 50 Cl), a = 0.8 sigma
          (~2.7 A, real NaCl: 2.82 A) */
    Rng lattice_rng(1);  /* fixed seed */
    Run("NaClLattice", RockSaltFrame(5, 5, 4, 0.8, DEFAULT_T, lattice_rng), ff_cluster(), 6);

    /* 3. closed box: random start as in 1., inside WCA walls (the faces of
          the 10x10x10 box), placed at least 0.75 from the walls (beyond the
          wall cutoff 0.72 of Cl) */
    ForceField ff_box = ff_cluster();
    ff_box.Add(std::make_unique<WCAWall>());
    Rng box_rng(42);  /* fixed seed */
    XYZFrame boxed = RandomFrame(100, 10.0, Boundary::Open, DEFAULT_T, box_rng, 0.75);
    for (int i = 0; i < boxed.n(); ++i) boxed.set_element(i, i % 2 == 0 ? na : cl);
    Run("NaClBox", std::move(boxed), std::move(ff_box), 0);   /* walls: N_f = 3N */

    return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}
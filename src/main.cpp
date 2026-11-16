#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <filesystem>
#include <format>
#include <stdexcept>
#include "xyzio.h"
#include "rng.h"
#include "paths.h"
#include "params.h"
#include "simulation.h"
#include "potentials/lennard_jones.h"

/* Lennard-Jones only, default cutoff (one per simulation: a Simulation owns its force field) */
static ForceField LJ() {
    ForceField ff;
    ff.Add(std::make_unique<LennardJones>());
    return ff;
}

int main() {
    try {    
    /* checking things for now*/
    printf("--------------------------------\n");
    printf("Read the test input now:\n");
    XYZFrame input = xyz_read(INPUT_DIR "test.xyz");
    printf("Number of particles: %d\n", input.n());
    printf("Read the comment: %s\n",input.comment.c_str());
    printf("Box (Lattice): %f %f %f, %s boundaries\n", input.box().x, input.box().y, input.box().z,
         boundary_name(input.boundary()));
    printf("Now the coord and velocities:\n");
    printf("%-4s %10s %10s %10s %10s %10s %10s\n",
           "lbl", "x", "y", "z", "vx", "vy", "vz");
    for (int i = 0; i < input.n(); ++i) {
        Vec3 r = input.pos(i), v = input.vel(i);
        printf("%-4s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
               input.label(i), r.x, r.y, r.z, v.x, v.y, v.z);
    }
    printf("--------------------------------\n");
    printf("Testing Velocity-Verlet scheme for 1000 timesteps\n");

    /* NVE, every step written to the stage directory and printed */
    MDParams nve{ .dt = 0.005, .steps = 1000, .write_every = 1, .print_every = 1 };

    /* 1. random initial geometry, frames in output/RandomFrame/ */
    const int    n_particles = 20;
    const double box         = 6.0;    /* box edge length in sigma */
    Rng md_rng(42);  /* fixed seed: same trajectory every run */

    XYZFrame geom = RandomFrame(n_particles, box, Boundary::Periodic, DEFAULT_T, md_rng);
    xyz_write(INPUT_DIR "initialconfig.xyz", geom);
    printf("RandomFrame: %d particles, frames in output/RandomFrame/, "
           "initial geometry also in input/initialconfig.xyz\n", geom.n());
    Simulation random(std::move(geom), LJ(), 42, 3);   /* periodic: N_f = 3N - 3 */
    nve.dir = OUTPUT_DIR "RandomFrame/";
    random.Produce(nve);
    printf("--------------------------------\n");

    /* 2. lattice initial geometry (same as in tests/force_test.cpp),
          frames in output/LatticeFrame/ */
    Rng lattice_rng(1);  /* fixed seed */

    XYZFrame lattice = LatticeFrame(5, 1.2, DEFAULT_T, lattice_rng);
    printf("LatticeFrame: %d particles on a 5x5x5 simple cubic lattice (spacing 1.2), "
           "frames in output/LatticeFrame/\n", lattice.n());
    Simulation lat(std::move(lattice), LJ(), 1, 3);
    nve.dir = OUTPUT_DIR "LatticeFrame/";
    lat.Produce(nve);
    printf("--------------------------------\n");
    return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}


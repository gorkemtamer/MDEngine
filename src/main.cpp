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
#include "integrator.h"
#include "checkers.h"
#include "measure.h"
#include "forcefield.h"
#include "potentials/lennard_jones.h"

/* Runs n_steps velocity-Verlet steps of length dt on geom, checking the
   minimum distance after every step. Writes the start geometry
   (initialconfig.xyz) and every step (stepNNNN.xyz) into dir, which is
   created if needed, and prints the energies. 
   Later to be implemented as the main scheme.
   */
static void RunMD(XYZFrame& geom, const ForceField& ff, const std::string& dir, double dt, int n_steps) {
    std::filesystem::create_directories(dir);
    xyz_write(dir + "initialconfig.xyz", geom);

    double ekin = geom.KineticEnergy();
    double epot = PotentialEnergy(geom, ff);
    printf("%5s %8s %12s %12s %12s\n", "step", "t", "Ekin", "Epot", "Etot");
    printf("%5d %8.3f %12.6f %12.6f %12.6f\n", 0, 0.0, ekin, epot, ekin + epot);

    for (int step = 1; step <= n_steps; ++step) {
        /* updates geom in place, so it already is the new geometry */
        VelocityVerletStep(geom, ff, dt);
        if (CheckMinDistance(geom) != 0) throw std::runtime_error("minDistance failed at step " + std::to_string(step));
        geom.comment = std::format("MD step {}, t = {:g}", step, step * dt);
        xyz_write(dir + std::format("step{:04d}.xyz", step), geom);

        ekin = geom.KineticEnergy();
        epot = geom.PotentialEnergy();      /* forces are fresh after the step: no recomputation */
        printf("%5d %8.3f %12.6f %12.6f %12.6f\n", step, step * dt, ekin, epot, ekin + epot);
    }
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
    printf("Testing Velocity-Verlet scheme for 100 timesteps\n");

    /* interactions: Lennard-Jones only, default cutoff */
    ForceField ff;
    ff.Add(std::make_unique<LennardJones>());

    const double dt      = 0.005;  /* time step in tau */
    const int    n_steps = 1000;

    /* 1. random initial geometry, frames in output/RandomFrame/ */
    const int    n_particles = 20;
    const double box         = 6.0;    /* box edge length in sigma */
    Rng md_rng(42);  /* fixed seed: same trajectory every run */

    XYZFrame geom = RandomFrame(n_particles, box, Boundary::Periodic, DEFAULT_T, md_rng);
    xyz_write(INPUT_DIR "initialconfig.xyz", geom);
    printf("RandomFrame: %d particles, frames in output/RandomFrame/, "
           "initial geometry also in input/initialconfig.xyz\n", geom.n());
    RunMD(geom, ff, std::string(OUTPUT_DIR) + "RandomFrame/", dt, n_steps);
    printf("--------------------------------\n");

    /* 2. lattice initial geometry (same as in tests/force_test.cpp),
          frames in output/LatticeFrame/ */
    Rng lattice_rng(1);  /* fixed seed */

    XYZFrame lattice = LatticeFrame(5, 1.2, DEFAULT_T, lattice_rng);
    printf("LatticeFrame: %d particles on a 5x5x5 simple cubic lattice (spacing 1.2), "
           "frames in output/LatticeFrame/\n", lattice.n());
    RunMD(lattice, ff, std::string(OUTPUT_DIR) + "LatticeFrame/", dt, n_steps);
    printf("--------------------------------\n");
    return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <format>
#include <stdexcept>
#include "xyzio.h"
#include "rng.h"
#include "paths.h"
#include "params.h"
#include "integrator.h"
#include "checkers.h"
#include "measure.h"

/* this part will be moved later inside*/
double EPS = 1.;
double SIGMA = 1.;

double LD_pot(double r, double e, double s){
    return 4*e * (std::pow(s/r , 12) - std::pow(s/r , 6));
}

int main() {
    try {
        /* random number generator: fix the seed (e.g. 12345) for reproducible runs,
       time(nullptr) gives a different sequence every second */
    Rng rng;
    rng_seed(&rng, (uint64_t)time(nullptr));
    
    /* checking things for now*/
    printf("--------------------------------\n");
    printf("Checking everything:\n");
    printf("PI is %f, meaning math.h is called\n", M_PI);
    printf("LD parameters are %f and %f, meaning globals are working\n", EPS, SIGMA);
    printf("Check LD potential it should be 0: %f \n", LD_pot(1,EPS,SIGMA));
    printf("Some random number between 0 and 1: %f \n", rng_uniform(&rng));
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

    /* simple MD loop: random initial geometry (open boundaries), then
       velocity-Verlet steps, checking the geometry and saving every frame */
    const int    n_particles = 20;
    const double box         = 6.0;    /* box edge length in sigma */
    const double dt          = 0.005;  /* time step in tau */
    const int    n_steps     = 1000;

    Rng md_rng;
    rng_seed(&md_rng, 42);             /* fixed seed: same trajectory every run */

    XYZFrame geom = RandomFrame(n_particles, box, Boundary::Periodic, DEFAULT_T, md_rng);
    xyz_write(OUTPUT_DIR "initialconfig.xyz", geom);
    xyz_write(INPUT_DIR "initialconfig.xyz", geom);    
    printf("Initial geometry written to output/ and input/initialconfig.xyz\n");

    double ekin = geom.KineticEnergy();
    double epot = PotentialEnergy(geom);
    printf("%5s %8s %12s %12s %12s\n", "step", "t", "Ekin", "Epot", "Etot");
    printf("%5d %8.3f %12.6f %12.6f %12.6f\n", 0, 0.0, ekin, epot, ekin + epot);

    for (int step = 1; step <= n_steps; ++step) {
        /* updates geom in place, so it already is the new geometry */
        VelocityVerletStep(geom, dt);
        if (CheckMinDistance(geom) != 0) throw std::runtime_error("minDistance failed at step " + std::to_string(step));
        geom.comment = std::format("MD step {}, t = {:g}", step, step * dt);
        xyz_write(std::string(OUTPUT_DIR) + std::format("step{:04d}.xyz",step),geom);

        ekin = geom.KineticEnergy();
        epot = geom.PotentialEnergy();      /* forces are fresh after the step: no recomputation */
        printf("%5d %8.3f %12.6f %12.6f %12.6f\n", step, step * dt, ekin, epot, ekin + epot);
    }

    printf("--------------------------------\n");
    return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}


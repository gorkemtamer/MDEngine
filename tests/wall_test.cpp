// wall_test.cpp
//
// WCA walls (the faces of the box) keep the particles inside a closed box.
//   1. one Cl against the formula: at z = s = sigma_i/2 from a wall
//      U = eps and F = 24 eps / s, pointing into the box; beyond the cutoff
//      2^(1/6) s nothing
//   2. Lennard-Jones gas in a walled box (NVE): energy conserved while the
//      particles bounce off the walls, error ~4x smaller for dt/2
//   3. rejected: periodic frames, particles outside the box, eps <= 0

#include <cmath>
#include <exception>
#include <memory>
#include <stdexcept>
#include "xyzio.h"
#include "rng.h"
#include "measure.h"
#include "integrators/velocity_verlet.h"
#include "potentials/wca_wall.h"
#include "test_helpers.h"

// Largest |E - E0| per particle over a run of length t_end with step dt.
static double max_energy_error(double dt, double t_end) {
    Rng rng(4);                                   // same start for every dt
    XYZFrame f = RandomFrame(30, 5.0, Boundary::Open, 2.0, rng, 0.6);
    ForceField ff = LJForceField();
    ff.Add(std::make_unique<WCAWall>());
    VelocityVerlet vv(ff, dt);
    const double E0 = f.KineticEnergy() + PotentialEnergy(f, ff);
    const int steps = static_cast<int>(std::lround(t_end / dt));
    double err = 0;
    for (int s = 0; s < steps; ++s) {
        vv.Step(f);
        err = std::fmax(err, std::fabs(f.KineticEnergy() + f.PotentialEnergy() - E0));
    }
    return err / f.n();
}

int main() {
    try {
        ForceField walls;
        walls.Add(std::make_unique<WCAWall>(2.0));     // eps = 2

        // 1. one Cl (not Ar: sigma_i must enter) at z = s from the wall x = 0,
        //    then from the wall y = Ly, then just beyond the cutoff
        XYZFrame one(1, Vec3{ 6, 6, 6 }, Boundary::Open);
        one.set_element(0, element_index("Cl"));
        const double s = 0.5 * one.element(0).sigma;  // 0.644
        const double F0 = 24 * 2.0 / s;                // force at z = s, eps = 2
        one.set_pos(0, Vec3{ s, 3, 3 });
        double U = PotentialEnergy(one, walls);        // expected eps = 2
        Vec3 F = one.force(0);                         // expected F0 in +x
        check(std::fabs(U - 2) < 1e-12 && norm(F - Vec3{ F0, 0, 0 }) < 1e-10,
              "low wall: U and F vs. formula", std::fabs(U - 2) + norm(F - Vec3{ F0, 0, 0 }));
        one.set_pos(0, Vec3{ 3, 6 - s, 3 });
        U = PotentialEnergy(one, walls);
        F = one.force(0);
        check(std::fabs(U - 2) < 1e-12 && norm(F - Vec3{ 0, -F0, 0 }) < 1e-10,
              "high wall: U and F vs. formula", std::fabs(U - 2) + norm(F - Vec3{ 0, -F0, 0 }));
        one.set_pos(0, Vec3{ 3, 3, 1.123 * s });      // beyond 2^(1/6) s = 1.1225 s
        U = PotentialEnergy(one, walls);
        check(U == 0 && norm(one.force(0)) == 0, "beyond the cutoff: no energy, no force");

        // 2. LJ gas bouncing off the walls (T = 2: many wall hits in 5 tau).
        //    The error comes mostly from close LJ pairs of the random start
        //    (same size without walls), hence the loose bound; the ratio is
        //    the real check.
        const double e1 = max_energy_error(0.005, 5.0);
        const double e2 = max_energy_error(0.0025, 5.0);
        check(e1 < 5e-2, "LJ + walls, dt = 0.005:  max |E - E0| per particle", e1);
        check(e2 < 5e-2, "LJ + walls, dt = 0.0025: max |E - E0| per particle", e2);
        check(e1 / e2 > 3.0 && e1 / e2 < 6.0, "error ratio dt / (dt/2), expected ~4", e1 / e2);

        // 3. rejected
        bool thrown = false;
        XYZFrame p(1, Vec3{ 6, 6, 6 }, Boundary::Periodic);
        p.set_element(0, element_index("Ar"));
        p.set_pos(0, Vec3{ 3, 3, 3 });
        try { walls.Compute(p); } catch (const std::runtime_error&) { thrown = true; }
        check(thrown, "periodic frame is rejected");

        thrown = false;
        one.set_pos(0, Vec3{ 3, 3, 6.1 });
        try { walls.Compute(one); } catch (const std::runtime_error&) { thrown = true; }
        check(thrown, "particle outside the box is rejected");

        thrown = false;
        try { WCAWall bad(0.0); } catch (const std::invalid_argument&) { thrown = true; }
        check(thrown, "eps = 0 is rejected");
    } catch (const std::exception& e) {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
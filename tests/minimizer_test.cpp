// minimizer_test.cpp
//
// Steepest descent against systems with a known minimum:
//   1. particles in a harmonic trap end at the center (U = 0)
//   2. an Ar-Kr pair ends at r = 2^(1/6) sigma_ij with U = -eps_ij (+ shift)
//   3. four Ar atoms end as a regular tetrahedron: 6 pairs at the LJ minimum
//   4. periodic LJ lattice: U goes down, positions stay in the box
//   5. velocities are not touched, forces and energy belong to the result
//   6. step limit reached: converged = false; step length grows by 1.2 up to max_step
//   7. bad parameters and Coulomb alone are rejected

#include <cmath>
#include <exception>
#include <memory>
#include <stdexcept>
#include "minimizers/steepest_descent.h"
#include "potentials/harmonic_trap.h"
#include "potentials/coulomb.h"
#include "xyzio.h"
#include "rng.h"
#include "test_helpers.h"

// LJ pair energy of elements a, b at distance r, truncated and shifted as in LennardJones
static double lj_pair(const char* a, const char* b, double r) {
    const Element& ea = element(element_index(a));
    const Element& eb = element(element_index(b));
    const double s = 0.5 * (ea.sigma + eb.sigma), e = std::sqrt(ea.eps * eb.eps);
    auto lj = [&](double x) { double q = std::pow(s / x, 6); return 4 * e * (q * q - q); };
    return lj(r) - lj(LJ_CUTOFF * s);
}

int main() {
    try {
        const double tol = 1e-6;
        const double rmin_factor = std::pow(2.0, 1.0 / 6.0);

        // 1. three particles in a trap (no interaction between them)
        {
            ForceField ff;
            ff.Add(std::make_unique<HarmonicTrap>(2.0, Vec3{ 5, 5, 5 }));
            XYZFrame f(3, Vec3{ 10, 10, 10 }, Boundary::Open);
            const Vec3 start[3] = { { 6, 5.5, 4 }, { 3, 5, 5 }, { 5, 7.5, 5.2 } };
            for (int i = 0; i < 3; ++i) { f.set_element(i, element_index("Ar")); f.set_pos(i, start[i]); }
            SteepestDescent sd(ff, tol, 10000);
            MinimizeResult r = sd.Minimize(f);
            double dist = 0;
            for (int i = 0; i < 3; ++i) dist = max_nan(dist, norm(f.pos(i) - Vec3{ 5, 5, 5 }));
            check(r.converged, "trap: converged (max |F| < tol)", r.max_force);
            check(dist < tol / 2.0, "trap: max distance from the center", dist);   // |F| = k |r - c|
            check(std::fabs(r.energy) < 1e-12, "trap: U at the end, expected 0", r.energy);
        }

        // 2. Ar-Kr pair
        {
            const ForceField ff = LJForceField();
            XYZFrame f(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            f.set_element(0, element_index("Ar"));
            f.set_element(1, element_index("Kr"));
            f.set_pos(0, Vec3{ 4, 5, 5 });
            f.set_pos(1, Vec3{ 5.6, 5.3, 4.9 });                  // r = 1.63, beyond the minimum
            SteepestDescent sd(ff, tol, 10000);
            MinimizeResult r = sd.Minimize(f);
            const double s_ij = 0.5 * (element(element_index("Ar")).sigma + element(element_index("Kr")).sigma);
            const double r_min = rmin_factor * s_ij;
            const double d = norm(f.pos(1) - f.pos(0));
            check(r.converged, "Ar-Kr: converged", r.max_force);
            check(std::fabs(d - r_min) < 1e-6, "Ar-Kr: |r - 2^(1/6) sigma_ij|", std::fabs(d - r_min));
            check(std::fabs(r.energy - lj_pair("Ar", "Kr", r_min)) < 1e-12, "Ar-Kr: U vs. -eps_ij (shifted)",
                  std::fabs(r.energy - lj_pair("Ar", "Kr", r_min)));
        }

        // 3. Ar4 from a distorted tetrahedron
        {
            const ForceField ff = LJForceField();
            XYZFrame f(4, Vec3{ 10, 10, 10 }, Boundary::Open);
            const Vec3 start[4] = { { 5, 5, 5 }, { 6.3, 5.1, 4.9 }, { 5.5, 6.0, 5.2 }, { 5.6, 5.4, 6.1 } };
            for (int i = 0; i < 4; ++i) { f.set_element(i, element_index("Ar")); f.set_pos(i, start[i]); }
            SteepestDescent sd(ff, tol, 10000);
            MinimizeResult r = sd.Minimize(f);
            double derr = 0;
            for (int i = 0; i < 4; ++i)
                for (int j = i + 1; j < 4; ++j) derr = max_nan(derr, std::fabs(norm(f.pos(j) - f.pos(i)) - rmin_factor));
            const double U_exp = 6 * lj_pair("Ar", "Ar", rmin_factor);
            check(r.converged, "Ar4: converged", r.max_force);
            check(derr < 1e-6, "Ar4: all 6 distances at 2^(1/6), max error", derr);
            check(std::fabs(r.energy - U_exp) < 1e-12, "Ar4: U vs. 6 * pair minimum", std::fabs(r.energy - U_exp));
        }

        // 4. + 5. periodic LJ lattice with thermal velocities, randomly displaced
        {
            const ForceField ff = LJForceField();
            Rng rng(3);
            XYZFrame f = LatticeFrame(5, 1.2, 1.0, rng);
            // shifted by half a spacing: the lattice sites lie on the box faces,
            // so particles cross them during the minimization
            const Vec3 shift{ 0.6, 0.6, 0.6 };
            for (int i = 0; i < f.n(); ++i)
                f.set_pos(i, f.Wrap(f.pos(i) - shift + 0.1 * Vec3{ rng.uniform() - 0.5, rng.uniform() - 0.5, rng.uniform() - 0.5 }));
            const XYZFrame before = f;
            const double U0 = f.PotentialEnergy(ff);
            SteepestDescent sd(ff, 1e-3, 20000);
            MinimizeResult r = sd.Minimize(f);
            bool inside = true, same_v = true;
            for (int i = 0; i < f.n(); ++i) {
                Vec3 p = f.pos(i);
                inside = inside && p.x >= 0 && p.x < f.box().x && p.y >= 0 && p.y < f.box().y && p.z >= 0 && p.z < f.box().z;
                same_v = same_v && norm(f.vel(i) - before.vel(i)) == 0.0;
            }
            check(r.converged, "periodic lattice: converged (tol 1e-3)", r.max_force);
            check(r.energy < U0, "periodic lattice: U lowered", r.energy - U0);
            check(inside, "periodic lattice: positions stay in the box");
            check(same_v, "velocities not touched");
            const double U_stored = r.energy;
            ff.Compute(f);
            check(U_stored == f.PotentialEnergy(ff) && r.max_force == MaxForce(f),
                  "result energy and forces belong to the final positions");
        }

        // 6. step limit
        {
            const ForceField ff = LJForceField();
            XYZFrame f(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            f.set_element(0, element_index("Ar"));
            f.set_element(1, element_index("Ar"));
            f.set_pos(0, Vec3{ 4, 5, 5 });
            f.set_pos(1, Vec3{ 5.5, 5, 5 });
            SteepestDescent sd(ff, tol, 3);
            MinimizeResult r = sd.Minimize(f);
            check(!r.converged && r.steps == 3, "step limit: stops after 3 steps, not converged");
        }
        {
            // one particle far from a weak trap center: every step is downhill,
            // so it moves by h = 0.01, 0.012, ... up to 0.1 per step
            ForceField ff;
            ff.Add(std::make_unique<HarmonicTrap>(0.1, Vec3{ 5, 5, 5 }));
            XYZFrame f(1, Vec3{ 10, 10, 10 }, Boundary::Open);
            f.set_element(0, element_index("Ar"));
            f.set_pos(0, Vec3{ 9, 5, 5 });
            SteepestDescent sd(ff, tol, 30, 0.01, 0.1);
            sd.Minimize(f);
            double x = 9, h = 0.01;
            for (int s = 0; s < 30; ++s) { x -= h; h = std::fmin(1.2 * h, 0.1); }
            check(std::fabs(f.pos(0).x - x) < 1e-12, "step length: 30 accepted steps, growth 1.2, cap 0.1",
                  std::fabs(f.pos(0).x - x));
        }

        // 7. rejections
        {
            const ForceField ff = LJForceField();
            auto throws = [](auto make) {
                try { make(); } catch (const std::invalid_argument&) { return true; }
                return false;
            };
            check(throws([&] { SteepestDescent(ff, 0.0, 10); }), "force_tol = 0 rejected");
            check(throws([&] { SteepestDescent(ff, 1e-3, 0); }), "max_steps = 0 rejected");
            check(throws([&] { SteepestDescent(ff, 1e-3, 10, 0.0); }), "initial step = 0 rejected");
            check(throws([&] { SteepestDescent(ff, 1e-3, 10, 0.1, 0.05); }), "max step < initial step rejected");

            ForceField coulomb_only;
            coulomb_only.Add(std::make_unique<Coulomb>());
            XYZFrame f(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            f.set_element(0, element_index("Na"));
            f.set_element(1, element_index("Cl"));
            f.set_pos(0, Vec3{ 4, 5, 5 });
            f.set_pos(1, Vec3{ 6, 5, 5 });
            bool thrown = false;
            SteepestDescent sd(coulomb_only, 1e-3, 100);
            try { sd.Minimize(f); } catch (const std::runtime_error&) { thrown = true; }
            check(thrown, "Coulomb alone rejected (no minimum)");
        }
    } catch (const std::exception& e) {
        std::printf("FAIL unexpected exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
// velocity_test.cpp
//
// AssignVelocities against the Maxwell-Boltzmann distribution, for a
// mixture of light and heavy particles:
//   1. total momentum zero after the center-of-mass velocity is removed
//   2. equipartition per element: <m v^2 / 2> = 3/2 T for Ar and for Kr
//      (fails if the mass is ignored)
//   3. variance of one component per element: <v_x^2> = T / m
//   4. positions untouched, forces stay valid
//   5. T = 0 gives zero velocities; T < 0 and particles without an
//      element are rejected; an empty frame is fine

#include <cmath>
#include <exception>
#include <stdexcept>
#include "velocities.h"
#include "test_helpers.h"

int main() {
    try {
        const double T = 1.7;
        const int n = 20000;                            // half Ar (m = 1), half Kr (m = 2.1)
        const int ar = element_index("Ar"), kr = element_index("Kr");
        XYZFrame f(n, Vec3{ 10, 10, 10 }, Boundary::Open);
        for (int i = 0; i < n; ++i) {
            f.set_element(i, i % 2 == 0 ? ar : kr);
            f.set_pos(i, Vec3{ 0.1 * (i % 100), 0.1 * (i / 100 % 100), 5 });
        }
        const XYZFrame before = f;

        Rng rng(11);
        AssignVelocities(f, T, rng);

        // 1. momentum, relative to the sum of all |p_i| (rounding only)
        double p_abs = 0;
        for (int i = 0; i < n; ++i) p_abs += f.mass(i) * norm(f.vel(i));
        const double p_rel = norm(f.TotalMomentum()) / p_abs;
        check(p_rel < 1e-12, "|total momentum| / sum |p_i|", p_rel);

        // 2. + 3. per element: kinetic energy per particle and <v_x^2>
        double k[2] = { 0, 0 }, vx2[2] = { 0, 0 };
        for (int i = 0; i < n; ++i) {
            const int e = (f.type(i) == ar) ? 0 : 1;
            k[e] += 0.5 * f.mass(i) * norm2(f.vel(i));
            vx2[e] += f.vel(i).x * f.vel(i).x;
        }
        const int half = n / 2;
        // relative standard deviation of a mean of 3*half (resp. half) squared Gaussians: sqrt(2/(3*half)), sqrt(2/half)
        const double kAr = k[0] / half / (1.5 * T), kKr = k[1] / half / (1.5 * T);
        const double vAr = vx2[0] / half / (T / element(ar).mass), vKr = vx2[1] / half / (T / element(kr).mass);
        check(std::fabs(kAr - 1) < 5 * std::sqrt(2.0 / (3 * half)), "Ar: <m v^2/2> / (3/2 T) - 1", kAr - 1);
        check(std::fabs(kKr - 1) < 5 * std::sqrt(2.0 / (3 * half)), "Kr: <m v^2/2> / (3/2 T) - 1", kKr - 1);
        check(std::fabs(vAr - 1) < 5 * std::sqrt(2.0 / half), "Ar: <v_x^2> / (T/m) - 1", vAr - 1);
        check(std::fabs(vKr - 1) < 5 * std::sqrt(2.0 / half), "Kr: <v_x^2> / (T/m) - 1", vKr - 1);

        // 4. positions untouched, forces still valid
        {
            XYZFrame g(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            g.set_element(0, ar); g.set_element(1, kr);
            g.set_pos(0, Vec3{ 4, 5, 5 }); g.set_pos(1, Vec3{ 5.2, 5, 5 });
            const ForceField ff = LJForceField();
            ff.Compute(g);
            AssignVelocities(g, T, rng);
            bool same_pos = true;
            for (int i = 0; i < f.n(); ++i) same_pos = same_pos && norm(f.pos(i) - before.pos(i)) == 0.0;
            check(same_pos && g.forces_valid(), "positions untouched, forces still valid");
        }

        // 5. T = 0, rejections, empty frame
        {
            AssignVelocities(f, 0.0, rng);
            double vmax = 0;
            for (int i = 0; i < n; ++i) vmax = std::fmax(vmax, norm(f.vel(i)));
            check(vmax == 0.0, "T = 0: all velocities zero", vmax);

            bool thrown = false;
            try { AssignVelocities(f, -1.0, rng); } catch (const std::invalid_argument&) { thrown = true; }
            check(thrown, "T < 0 rejected");

            XYZFrame u(3, Vec3{ 10, 10, 10 }, Boundary::Open);   // elements not set
            u.set_element(0, ar);
            thrown = false;
            try { AssignVelocities(u, T, rng); } catch (const std::invalid_argument&) { thrown = true; }
            check(thrown, "particle without an element rejected");

            XYZFrame empty(0, Vec3{ 10, 10, 10 }, Boundary::Open);
            AssignVelocities(empty, T, rng);
            check(true, "empty frame: no error");
        }
    } catch (const std::exception& e) {
        std::printf("FAIL unexpected exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
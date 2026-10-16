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
// RemoveAngularMomentum:
//   6. a rigid rotation + translation keeps only the translation
//      (3D cluster and a linear chain, where the inertia tensor is singular)
//   7. random cluster: L = 0 afterwards, momentum and positions unchanged
//   8. <K> = (3N - 6) T / 2 after both removals (N_f of an open cluster)
//   9. periodic frames rejected; one particle unchanged


#include <cmath>
#include <exception>
#include <stdexcept>
#include "velocities.h"
#include "test_helpers.h"
#include "xyzio.h"

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
            for (int i = 0; i < n; ++i) vmax = max_nan(vmax, norm(f.vel(i)));
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

        // 6. rigid rotation omega0 about the center of mass plus a translation u
        {
            const Vec3 omega0{ 0.3, -1.1, 0.7 }, u{ 0.2, 0.1, -0.4 };
            auto rigid = [&](XYZFrame& g) {
                const Vec3 R = g.CenterOfMass();
                for (int i = 0; i < g.n(); ++i) g.set_vel(i, u + cross(omega0, g.pos(i) - R));
            };
            // a sum, not std::fmax: fmax ignores NaN, a sum keeps it (and NaN < tol is false)
            auto sum_dev = [&](const XYZFrame& g) {
                double d = 0;
                for (int i = 0; i < g.n(); ++i) d += norm(g.vel(i) - u);
                return d;
            };
            XYZFrame c(5, Vec3{ 10, 10, 10 }, Boundary::Open);          // 3D cluster, mixed masses
            const Vec3 pc[5] = { { 4, 5, 5 }, { 5.2, 5.1, 4.8 }, { 4.6, 6.0, 5.3 }, { 4.9, 5.4, 6.2 }, { 5.8, 4.3, 5.6 } };
            for (int i = 0; i < 5; ++i) { c.set_element(i, i % 2 == 0 ? ar : kr); c.set_pos(i, pc[i]); }
            rigid(c);
            RemoveAngularMomentum(c);
            check(sum_dev(c) < 1e-12, "rigid rotation + translation: sum |v_i - u|", sum_dev(c));

            XYZFrame l(3, Vec3{ 10, 10, 10 }, Boundary::Open);          // linear chain along (1, 2, -1)
            for (int i = 0; i < 3; ++i) {
                l.set_element(i, i == 1 ? kr : ar);
                l.set_pos(i, Vec3{ 3, 3, 6 } + (0.9 * i * i + 1.1 * i) * Vec3{ 1, 2, -1 } / std::sqrt(6.0));
            }
            rigid(l);       // the part of omega0 along the chain moves nothing
            RemoveAngularMomentum(l);
            check(sum_dev(l) < 1e-12, "linear chain: sum |v_i - u|", sum_dev(l));
        }

        // 7. random open cluster with mixed masses
        {
            Rng r7(5);
            XYZFrame g = RandomFrame(40, 6.0, Boundary::Open, 1.0, r7);
            for (int i = 0; i < g.n(); i += 3) g.set_element(i, kr);
            AssignVelocities(g, 1.0, r7);
            const XYZFrame g0 = g;
            const Vec3 R = g.CenterOfMass();
            double l_abs = 0;
            for (int i = 0; i < g.n(); ++i) l_abs += g.mass(i) * norm(cross(g.pos(i) - R, g.vel(i)));
            RemoveAngularMomentum(g);
            const double l_rel = norm(g.AngularMomentum()) / l_abs;
            check(l_rel < 1e-12, "random cluster: |L| / sum |l_i| after removal", l_rel);
            const double dp = norm(g.TotalMomentum() - g0.TotalMomentum());
            bool same_pos = true;
            for (int i = 0; i < g.n(); ++i) same_pos = same_pos && norm(g.pos(i) - g0.pos(i)) == 0.0;
            check(dp < 1e-13 && same_pos, "random cluster: momentum and positions unchanged", dp);
        }

        // 8. N_f = 3N - 6: average K over many small clusters
        {
            Rng r8(8);
            const int frames = 20000, n5 = 5;
            double ksum = 0;
            for (int s = 0; s < frames; ++s) {
                XYZFrame g = RandomFrame(n5, 3.0, Boundary::Open, T, r8);
                g.set_element(0, kr);
                AssignVelocities(g, T, r8);
                RemoveAngularMomentum(g);
                ksum += g.KineticEnergy();
            }
            const double ratio = ksum / frames / ((3 * n5 - 6) * T / 2);
            // per frame K has relative standard deviation sqrt(2 / N_f); 5 sigma of the mean
            check(std::fabs(ratio - 1) < 5 * std::sqrt(2.0 / (3 * n5 - 6) / frames), "<K> / ((3N - 6) T / 2) - 1", ratio - 1);
        }

        // 9. rejections, one particle
        {
            XYZFrame per(2, Vec3{ 10, 10, 10 }, Boundary::Periodic);
            per.set_element(0, ar); per.set_element(1, ar);
            bool thrown = false;
            try { RemoveAngularMomentum(per); } catch (const std::runtime_error&) { thrown = true; }
            check(thrown, "RemoveAngularMomentum: periodic frame rejected");

            XYZFrame one(1, Vec3{ 10, 10, 10 }, Boundary::Open);
            one.set_element(0, ar);
            one.set_pos(0, Vec3{ 2, 3, 4 });
            one.set_vel(0, Vec3{ 1, 2, 3 });
            RemoveAngularMomentum(one);
            check(norm(one.vel(0) - Vec3{ 1, 2, 3 }) == 0.0, "RemoveAngularMomentum: one particle unchanged");
        }
    } catch (const std::exception& e) {
        std::printf("FAIL unexpected exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
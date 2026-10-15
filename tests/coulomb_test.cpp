// coulomb_test.cpp
//
// Coulomb potential (lib/potentials/coulomb.h) against formulas:
//   1. one Na+ Cl- pair: U = -K / r, force attractive with |F| = K / r^2
//   2. Na+, Ar, Cl- in a row: the neutral Ar between them changes nothing
//      (pairs with a neutral particle are skipped, the loop must go on)
//   3. relative permittivity eps_r: energy and force scale by 1/eps_r
//   4. a periodic frame is rejected (would need Ewald summation)
// Forces vs. -dU/dr and Newton's third law with charges are checked in
// force_test and newton_test.

#include <cmath>
#include <exception>
#include <memory>
#include "measure.h"
#include "units.h"
#include "potentials/coulomb.h"
#include "test_helpers.h"

// open frame with the given elements, all on the x axis
static XYZFrame line_frame(std::initializer_list<const char*> labels, std::initializer_list<double> xs) {
    XYZFrame f(static_cast<int>(labels.size()), Vec3{ 20, 20, 20 }, Boundary::Open);
    int i = 0;
    for (const char *l : labels) f.set_element(i++, element_index(l));
    i = 0;
    for (double x : xs) f.set_pos(i++, Vec3{ x, 10, 10 });
    return f;
}

int main() {
    try {
        ForceField ff;
        ff.Add(std::make_unique<Coulomb>());
        const double r = 0.9;

        // 1. one ion pair
        XYZFrame pair = line_frame({ "Na", "Cl" }, { 10, 10 + r });
        double U = PotentialEnergy(pair, ff);
        double expected = -COULOMB_K / r;
        check(std::fabs(U - expected) < 1e-12 * std::fabs(expected),
              "Na-Cl pair energy vs. -K/r (relative error)", std::fabs(U - expected) / std::fabs(expected));
        // Na at the smaller x is pulled towards Cl: force in +x direction
        double Fx = pair.force(0).x, Fexp = COULOMB_K / (r * r);
        check(std::fabs(Fx - Fexp) < 1e-12 * Fexp && norm(pair.force(0) + pair.force(1)) < 1e-12 * Fexp,
              "Na-Cl force: attractive, |F| = K/r^2 (relative error)", std::fabs(Fx - Fexp) / Fexp);

        // 2. a neutral particle in between changes nothing
        XYZFrame three = line_frame({ "Na", "Ar", "Cl" }, { 10, 10 + r / 2, 10 + r });
        U = PotentialEnergy(three, ff);
        check(std::fabs(U - expected) < 1e-12 * std::fabs(expected) && norm(three.force(1)) == 0.0,
              "Na, Ar, Cl: Ar neither changes the energy nor feels a force", std::fabs(U - expected) / std::fabs(expected));

        // 3. relative permittivity: everything scaled by 1/eps_r
        const double eps_r = 78.4;
        ForceField water;
        water.Add(std::make_unique<Coulomb>(eps_r));
        XYZFrame pw = line_frame({ "Na", "Cl" }, { 10, 10 + r });
        double Uw = PotentialEnergy(pw, water);
        double errU = std::fabs(Uw * eps_r - expected) / std::fabs(expected);
        double errF = std::fabs(pw.force(0).x * eps_r - Fexp) / Fexp;
        check(std::fmax(errU, errF) < 1e-12, "eps_r = 78.4: energy and force scaled by 1/eps_r", std::fmax(errU, errF));

        // 4. periodic frames are rejected
        XYZFrame p(2, Vec3{ 20, 20, 20 }, Boundary::Periodic);
        p.set_element(0, element_index("Na"));
        p.set_element(1, element_index("Cl"));
        p.set_pos(1, Vec3{ 1, 0, 0 });
        bool thrown = false;
        try { ff.Compute(p); } catch (const std::runtime_error&) { thrown = true; }
        check(thrown, "periodic frame is rejected");
    } catch (const std::exception& e) {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
// newton_test.cpp
//
// Newton's third law (force on j from i = minus force on i from j) and
// its consequences:
//   1. the forces of every frame add up to zero
//   2. the total momentum stays constant during a run
//   3. with open boundaries, the total angular momentum also stays constant.
//      This additionally needs every pair force to point along the line
//      between the two particles (a central force). Periodic boundaries
//      break rotational symmetry, so there it is not conserved.
// Velocity Verlet keeps both momenta constant up to rounding, so the
// tolerances only allow for rounding errors.
// So idea is simple, check the momentum conservation

#include <cmath>
#include <exception>
#include "integrator.h"
#include "xyzio.h"
#include "test_helpers.h"

// |sum_i F_i| / sum_i |F_i|: zero if all forces cancel
static double force_sum_error(XYZFrame& f) {
    ComputeForces(f);
    Vec3 sum{};
    double scale = 0;
    for (int i = 0; i < f.n(); ++i) {
        sum += f.force(i);
        scale += norm(f.force(i));
    }
    return norm(sum) / scale;
}

// total angular momentum about the origin, sum_i m_i r_i x v_i
static Vec3 angular_momentum(const XYZFrame& f){
    Vec3 L{};
    for (int i = 0; i < f.n(); ++i) L += f.mass(i) * cross(f.pos(i), f.vel(i));
    return L;
}

// sum_i m_i |v_i|: the size against which a change in momentum
// is measured
static double momentum_scale(const XYZFrame& f){
    double s = 0;
    for (int i = 0; i < f.n(); ++i) s += f.mass(i) * norm(f.vel(i));
    return s;
}

// sum_i m_i |r_i x v_i|: same for the angular momentum
static double angular_scale(const XYZFrame& f){
    double s = 0;
    for (int i = 0; i < f.n(); ++i) s += f.mass(i) * norm(cross(f.pos(i),f.vel(i)));
    return s;
}

int main() {
    try {
        const double dt = 0.005;
        const int steps = 500;
        Rng rng(3);
        // perturbed lattice with Ar, Kr, Ne (as in force_test)
        XYZFrame p = LatticeFrame(5, 1.2, 1., rng);
        const int el[3] = {element_index("Ar"), element_index("Kr"), element_index("Ne")};
        for (int i=0; i < p.n(); ++i){
            p.set_element(i, el[i%3]);
            Vec3 d{0.2 * rng.uniform() - 0.1, 0.2 * rng.uniform() - 0.1, 0.2 * rng.uniform() - 0.1};
            p.set_pos(i, p.pos(i) + d);
        } 
        double err = force_sum_error(p);
        check( err < 1e-12, "periodic: |sum F| / |F|", err);
                Vec3 P0 = p.TotalMomentum();
        double dP = 0, scale = momentum_scale(p);
        for (int s = 0; s < steps; ++s) {
            VelocityVerletStep(p, dt);
            dP = std::fmax(dP, norm(p.TotalMomentum() - P0));
        }
        check(dP / scale < 1e-12, "periodic: max |P(t) - P(0)| / sum m|v| (500 steps)", dP / scale);

        // open: random geometry
        XYZFrame o = RandomFrame(30, 6.0, Boundary::Open, 1.0, rng);
        err = force_sum_error(o);
        check(err < 1e-12, "open: |sum F| / sum |F|", err);

        P0 = o.TotalMomentum();
        Vec3 L0 = angular_momentum(o);
        double dL = 0;
        dP = 0;
        scale = momentum_scale(o);
        double lscale = angular_scale(o);
        for (int s = 0; s < steps; ++s) {
            VelocityVerletStep(o, dt);
            dP = std::fmax(dP, norm(o.TotalMomentum() - P0));
            dL = std::fmax(dL, norm(angular_momentum(o) - L0));
        }
        check(dP / scale < 1e-12, "open: max |P(t) - P(0)| / sum m|v| (500 steps)", dP / scale);
        check(dL / lscale < 1e-12, "open: max |L(t) - L(0)| / sum m|r x v| (500 steps)", dL / lscale);
    }
    catch(const std::exception& e)
    {
        std::printf("FAIL  exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
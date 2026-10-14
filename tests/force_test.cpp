/*
force_test.cpp

Checks that the computed forces are the negative
gradient of the potential energy, F_i = -dU/dr_i
by comparing every force component with a central
finite difference of the energy
    F_num = -( U(r_i + h e_k) - U(r_i - h e_k))/(2h)
Run for two frames:
 - periodic: perturbed lattice with a mix of Ar, Kr
 and Ne (also tests the mixing rules)
 - open: a random geometry from RandomFrame
The FD check only tests that forces and energy agree with
each other; a wrong parameter would be wrong in both
and still pass. So one Ar-Kr pair energy is also compared 
with LJ formula, evaluated independently here.
*/

#include <cmath>
#include <exception>
#include <vector>
#include "integrator.h"
#include "params.h"
#include "measure.h"
#include "xyzio.h"
#include "test_helpers.h"

// component k (0, 1, 2) of a vector
static double& comp(Vec3& v, int k)       { return k == 0 ? v.x : (k == 1 ? v.y : v.z); }
static double  comp(const Vec3& v, int k) { return k == 0 ? v.x : (k == 1 ? v.y : v.z); }

// Largest |F_num - F| over all particles and components, divided by the
// largest |F| (so the result is a relative error).
static double max_force_error(XYZFrame& f, const ForceField& ff, double h) {
    ff.Compute(f);
    std::vector<Vec3> F(f.n());
    double fmax = 0;
    for (int i = 0; i < f.n(); ++i) {
        F[i] = f.force(i);
        fmax = std::fmax(fmax, std::fmax(std::fabs(F[i].x), std::fmax(std::fabs(F[i].y), std::fabs(F[i].z))));
    }
    double err = 0;
    for (int i = 0; i < f.n(); ++i) {
        const Vec3 r0 = f.pos(i);
        for (int k = 0; k < 3; ++k) {
            Vec3 r = r0;
            comp(r, k) = comp(r, k) + h;  f.set_pos(i, r);  double up = PotentialEnergy(f, ff);
            comp(r, k) = comp(r0, k) - h; f.set_pos(i, r);  double um = PotentialEnergy(f, ff);
            f.set_pos(i, r0);
            double f_num = -(up - um) / (2 * h);
            err = std::fmax(err, std::fabs(f_num - comp(F[i], k)));
        }
    }
    return err / fmax;
}

int main(){
    try{
        const double h = 1e-6; //FD step (sigma)
        const double tol = 1e-6; // tol by rel error
        const ForceField ff = LJForceField();

        //periodic: 5x5x5 lattice (spacing 1.2), every particle shifted randomly
        //by up to 0.1 sigma per direction (a perfect lattice has zero forces by
        //symmetry), elements Ar,Kr,Ne in turn
        Rng rng(1);
        XYZFrame p = LatticeFrame(5,1.2,1.,rng);
        const int el[3] = { element_index("Ar"), element_index("Kr"), element_index("Ne") };
        for (int i = 0; i < p.n(); ++i) {
            p.set_element(i, el[i % 3]);
            Vec3 d{ 0.2 * rng.uniform() - 0.1, 0.2 * rng.uniform() - 0.1, 0.2 * rng.uniform() - 0.1 };
            p.set_pos(i, p.pos(i) + d);
        }
        double err = max_force_error(p, ff, h);
        check(err < tol, "periodic, Ar/Kr/Ne: max |F_num - F| / max |F|", err);

        // open: random geometry
        XYZFrame o = RandomFrame(30, 6.0, Boundary::Open, 1.0, rng);
        err = max_force_error(o, ff, h);
        check(err < tol, "open, Ar: max |F_num - F| / max |F|", err);

        // one Ar-Kr pair at distance r: U = 4 eps_ij [(s/r)^12 - (s/r)^6] - shift,
        // sigma_ij = (sigma_Ar + sigma_Kr)/2, eps_ij = sqrt(eps_Ar eps_Kr), shift = value at the cutoff
        const Element& ar = element(element_index("Ar"));
        const Element& kr = element(element_index("Kr"));
        const double s_ij = 0.5 * (ar.sigma + kr.sigma);
        const double e_ij = std::sqrt(ar.eps * kr.eps);
        const double r = 1.1 * s_ij;
        auto lj = [&](double x) { double q = std::pow(s_ij / x, 6); return 4 * e_ij * (q * q - q); };
        const double expected = lj(r) - lj(LJ_CUTOFF * s_ij);

        XYZFrame pair(2, Vec3{ 10, 10, 10 }, Boundary::Open);
        pair.set_element(0, element_index("Ar"));
        pair.set_element(1, element_index("Kr"));
        pair.set_pos(0, Vec3{ 4, 5, 5 });
        pair.set_pos(1, Vec3{ 4 + r, 5, 5 });
        double U = PotentialEnergy(pair, ff);
        check(std::fabs(U - expected) < 1e-12 * std::fabs(expected),
              "Ar-Kr pair energy vs. LJ formula (relative error)", std::fabs(U - expected) / std::fabs(expected));
    } catch(const std::exception& e){
        std::printf("FAIL exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}
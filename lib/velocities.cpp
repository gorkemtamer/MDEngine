//velocities.cpp
#include <format>
#include <stdexcept>
#include "velocities.h"

void AssignVelocities(XYZFrame& f, double T, Rng& rng) {
    if (T < 0) throw std::invalid_argument(std::format("AssignVelocities: temperature must not be negative (got {:g})", T));
    for (int i = 0; i < f.n(); ++i)
        if (f.type(i) < 0) throw std::invalid_argument(std::format("AssignVelocities: particle {} has no element (mass unknown)", i));

    for (int i = 0; i < f.n(); ++i){
        double m = f.mass(i);
        Vec3 v{rng.MaxwellBoltzmann(T, m),
               rng.MaxwellBoltzmann(T, m),
               rng.MaxwellBoltzmann(T, m)};
        f.set_vel(i,v);
    }
    RemoveCenterOfMassVelocity(f);
}

void RemoveCenterOfMassVelocity(XYZFrame& f){
    if (f.n() == 0) return;     // no mass: the center-of-mass velocity is undefined
    const Vec3 vcm = f.TotalMomentum() / f.TotalMass();
    for (int i = 0; i < f.n(); ++i) f.set_vel(i, f.vel(i) - vcm);
}


void RemoveAngularMomentum(XYZFrame& f) {
    if (f.boundary() == Boundary::Periodic)
        throw std::runtime_error("RemoveAngularMomentum: only for open boundaries (no angular momentum in a periodic box)");
    if (f.n() < 2) return;

    const Vec3 R = f.CenterOfMass();
    const Vec3 L = f.AngularMomentum();

    // inertia tensor about R: I = sum m (|r|^2 1 - r r^T), symmetric
    double xx = 0, yy = 0, zz = 0, xy = 0, xz = 0, yz = 0;
    for (int i = 0; i < f.n(); ++i) {
        const Vec3 r = f.pos(i) - R;
        const double m = f.mass(i);
        xx += m * r.x * r.x;  yy += m * r.y * r.y;  zz += m * r.z * r.z;
        xy += m * r.x * r.y;  xz += m * r.x * r.z;  yz += m * r.y * r.z;
    }
    const double mr2 = xx + yy + zz;                     // sum m |r|^2
    if (mr2 == 0) return;                                // all particles at R: no rotation
    const double Ixx = mr2 - xx, Iyy = mr2 - yy, Izz = mr2 - zz;
    const double Ixy = -xy, Ixz = -xz, Iyz = -yz;

    // omega = I^-1 L via the adjugate (cofactors of a symmetric matrix)
    const double cxx = Iyy * Izz - Iyz * Iyz, cyy = Ixx * Izz - Ixz * Ixz, czz = Ixx * Iyy - Ixy * Ixy;
    const double cxy = Ixz * Iyz - Ixy * Izz, cxz = Ixy * Iyz - Ixz * Iyy, cyz = Ixy * Ixz - Ixx * Iyz;
    const double det = Ixx * cxx + Ixy * cxy + Ixz * cxz;

    Vec3 omega;
    if (det > 1e-12 * mr2 * mr2 * mr2) {
        omega = Vec3{ cxx * L.x + cxy * L.y + cxz * L.z,
                      cxy * L.x + cyy * L.y + cyz * L.z,
                      cxz * L.x + cyz * L.y + czz * L.z } / det;
    } else {
        // linear: I = mr2 (1 - u u^T) for the axis u, and L is perpendicular to u
        omega = L / mr2;
    }
    for (int i = 0; i < f.n(); ++i) f.set_vel(i, f.vel(i) - cross(omega, f.pos(i) - R));
}
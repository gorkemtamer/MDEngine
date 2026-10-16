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
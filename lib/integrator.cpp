// integrator.cpp
//
// The velocity-Verlet integrator.

#include <format>
#include <stdexcept>
#include "integrator.h"
#include "vec3.h"

void VelocityVerletStep(XYZFrame& f, const ForceField& ff, double dt) {
    if (dt <= 0) throw std::invalid_argument(std::format("VelocityVerletStep: time step must be positive (got {:g})", dt));
    ff.CheckNotCollapsing();   // TEMPORARY, see forcefield.h

    // F(t): only needed on the first step, afterwards it is kept from the previous step
    if (!f.forces_valid()) ff.Compute(f);

    // half kick + drift
    for (int i = 0; i < f.n(); ++i) {
        double h = 0.5 * dt / f.mass(i);
        Vec3 v = f.vel(i) + h * f.force(i);
        f.set_vel(i,v);
        f.set_pos(i,f.Wrap(f.pos(i)+dt*v));
    }

    // F(t + dt)
    ff.Compute(f);

    // second half kick
    for (int i = 0; i < f.n(); ++i) {
        double h = 0.5 * dt / f.mass(i);
        f.set_vel(i,f.vel(i)+h*f.force(i));
    }
}

// velocity_verlet.cpp
#include "integrators/velocity_verlet.h"
#include "vec3.h"

void VelocityVerlet::Step(XYZFrame& f) {
    ff_.CheckNotCollapsing();   // TEMPORARY, see forcefield.h

    // F(t): only needed on the first step, afterwards it is kept from the previous step
    if (!f.forces_valid()) ff_.Compute(f);

    // half kick + drift
    for (int i = 0; i < f.n(); ++i) {
        double h = 0.5 * dt_ / f.mass(i);
        Vec3 v = f.vel(i) + h * f.force(i);
        f.set_vel(i,v);
        f.set_pos(i,f.Wrap(f.pos(i)+dt_*v));
    }

    // F(t + dt)
    ff_.Compute(f);

    // second half kick
    for (int i = 0; i < f.n(); ++i) {
        double h = 0.5 * dt_ / f.mass(i);
        f.set_vel(i,f.vel(i)+h*f.force(i));
    }
}
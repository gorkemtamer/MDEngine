#pragma once
#include "integrator.h"

// ---------------------------------------------------------------------
// Velocity-Verlet: second order, symplectic, time-reversible.
// One step of length dt:
//     v(t + dt/2) = v(t)        + dt/2 * F(t) / m
//     r(t + dt)   = r(t)        + dt   * v(t + dt/2)
//     F(t + dt)   = forces at the new positions
//     v(t + dt)   = v(t + dt/2) + dt/2 * F(t + dt) / m
// F(t) is reused from the previous step when f.forces_valid() is true,
// so each step costs one force evaluation. With periodic boundaries the
// positions are wrapped back into the box.
// ---------------------------------------------------------------------
class VelocityVerlet : public Integrator {
public:
    using Integrator::Integrator;     // same constructors as Integrator

    void Step(XYZFrame& f) override;
    int order() const override { return 2; }
};
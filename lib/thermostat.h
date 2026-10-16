#pragma once
#include "frame.h"

// ---------------------------------------------------------------------
// Interface for thermostats, in reduced units (see reduced_units.md).
// A thermostat couples the system to a heat bath at temperature T. It is
// applied after every integrator step and changes only the velocities, so
// the forces of the frame stay valid:
//     vv.Step(f);
//     th.Apply(f, vv.dt());
// Implementations live in lib/thermostats/.
// ---------------------------------------------------------------------
class Thermostat {
public:
    // Throws std::invalid_argument for T <= 0.
    explicit Thermostat(double T);
    virtual ~Thermostat() = default;

    // Couples f to the heat bath for one time step dt (the integrator's).
    virtual void Apply(XYZFrame& f, double dt) = 0;

    double temperature() const { return T_; }

    // Kinetic energy the thermostat has added so far (negative if it removed
    // energy). H = K + U - energy_added() is conserved by the exact dynamics:
    // its drift measures the integration error, as the total energy does
    // without a thermostat.
    double energy_added() const { return energy_added_; }

protected:
    double T_;
    double energy_added_ = 0.0;
};
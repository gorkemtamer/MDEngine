#pragma once
#include "frame.h"
#include "forcefield.h"

// ---------------------------------------------------------------------
// Interface for time integrators, in reduced units (see reduced_units.md).
// An integrator advances a frame by one time step dt, with the forces
// from a ForceField. Implementations live in lib/integrators/, e.g.
//     VelocityVerlet vv(ff, 0.005);
//     for (...) vv.Step(frame);
//
// The integrator does not own the force field, it only refers to it:
// ff must live at least as long as the integrator.
// Same design as the Potential class
// ---------------------------------------------------------------------

class Integrator {
public:
    //Throw std::invalid_argument for dt <= 0.
    Integrator(const ForceField& ff, double dt);
    
    // Forbids a temporary force field, e.g. VelocityVerlet(LJForceField(), dt):
    // it would be destroyed at the end of the line and leave ff_ dangling.   
    Integrator(ForceField&&, double) = delete;
    virtual ~Integrator() = default;

    //Advances f by one time step, in place. Not const:
    //integrators may keep a state between steps (e.g. earlier accelerations,
    //a random generator).
    //Errors from the force field propogate.
    virtual void Step(XYZFrame& f) = 0;

    // Order p of the method: the global error scales like dt^p, so halving
    // dt makes it 2^p times smaller (4 for velocity-Verlet).
    virtual int order() const = 0;

    double dt() const {return dt_;}
    const ForceField& force_field() const { return ff_;}

protected:
    const ForceField& ff_;
    double dt_;
};


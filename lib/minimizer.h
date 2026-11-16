#pragma once
#include "frame.h"
#include "forcefield.h"

// ---------------------------------------------
/*
Interface for energy minimizers, in reduced units (see reduced_units.md).
A minimizer moves the particles of a frame downhill in the potential
energy U of a ForceField until the largest force is below force_tol
(a local minimum, or the step limit). Velocities are not touched.
Implementations live in lib/minimizers/
    SteepestDescent sd(ff, 1e-4, 10000);
    MinimizeResult r = sd.Minimize(frame);
    if (!r.converged) ...
The minimizer does not own the force field, it only refers to it:
ff must live at least as long as the minimizer.
Same design as the Integrator class.
*/
// ---------------------------------------------

//What a minimization reached.
struct MinimizeResult {
    bool converged; //true if max_force < force_tol
    int steps; // #accepted+rejected steps
    double energy; //U at the end
    double max_force; //largest |F_i| at the end
};

class Minimizer {
public: 
    //force_tol: convergence criterion, largest |F_i| (eps/sigma)
    //max_steps: upper limit of steps
    //Throws std::invalid_argument for force_tol <=0 or max_steps <=0
    Minimizer(const ForceField& ff, double force_tol, int max_steps);

    //Forbids a temporary force field (would leave ff_ dangling).
    Minimizer(ForceField&&, double, int) = delete;
    virtual ~Minimizer() = default;

    //Minimizes U of f in place. On return the forces and the energy of f
    //belong to its positions. Not converging is not an error:
    //check the result. Errors of the force field propagate.
    virtual MinimizeResult Minimize(XYZFrame& f) = 0;

    double force_tol() const { return force_tol_; }
    int max_steps() const { return max_steps_; }
    const ForceField& force_field() const { return ff_; }
protected:
    const ForceField& ff_;
    double force_tol_;
    int max_steps_;
};

// Largest |F_i| of f; the forces must be up to date.
double MaxForce(const XYZFrame& f);


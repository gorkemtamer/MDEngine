#pragma once
#include "minimizer.h"

// ---------------------------------------------------------------------
// Steepest descent with an adaptive step (as in GROMACS). Each step moves
// every particle along its force, the largest one by the step length h:
//     r_i' = r_i + h * F_i / max_j |F_j|
// If U(r') < U(r) the step is accepted and h grows (h *= 1.2, at most
// max_step); otherwise it is rejected and h shrinks (h *= 0.2).
// Robust and simple, but slow near the minimum (zig-zag in narrow valleys).
// Stops when max |F_i| < force_tol (converged), after max_steps steps, or
// when h has shrunk below 1e-12 sigma (U cannot be lowered any further at
// double precision; not converged).
// max_step caps the displacement per step, which keeps particles from
// jumping through a neighbour or out of a walled box (WCAWall throws).
// With periodic boundaries the positions are wrapped back into the box.
// ---------------------------------------------------------------------

class SteepestDescent : public Minimizer {
public: 
    // initial_step, max_step: displacement of the particle with the
    // largest force, in sigma.
    // Throws std::invalid_argument for initial_step <= 0 or
    // max_step < initial_step, and for bad force_tol or max_steps.
    SteepestDescent(const ForceField& ff, double force_tol, int max_steps,
                    double initial_step = 0.01, double max_step = 0.1);
    SteepestDescent(ForceField&&, double, int, double = 0.01, double = 0.1) = delete;

    // Throws std::runtime_error for a Coulomb-only force field (see
    // ForceField::CheckNotCollapsing).
    MinimizeResult Minimize(XYZFrame& f) override;

    double initial_step() const { return initial_step_; }
    double max_step() const { return max_step_; }

private:
    double initial_step_;
    double max_step_;
};

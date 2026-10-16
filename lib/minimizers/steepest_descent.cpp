// steepest_descent.cpp
#include <algorithm>
#include <format>
#include <stdexcept>
#include "minimizers/steepest_descent.h"

SteepestDescent::SteepestDescent(const ForceField& ff, double force_tol, int max_steps,
                                 double initial_step, double max_step)
    : Minimizer(ff, force_tol, max_steps),
      initial_step_(initial_step),
      max_step_(max_step)
{
    if (initial_step <= 0)
        throw std::invalid_argument(std::format("SteepestDescent: initial step must be positive (got {:g})", initial_step));
    if (max_step < initial_step)
        throw std::invalid_argument(std::format("SteepestDescent: max step {:g} smaller than initial step {:g}", max_step, initial_step));
}

MinimizeResult SteepestDescent::Minimize(XYZFrame& f) {
    ff_.CheckNotCollapsing();   // TEMPORARY, see forcefield.h: Coulomb alone has no minimum

    double U = f.PotentialEnergy(ff_);
    double Fmax = MaxForce(f);
    double h = initial_step_;
    int step = 0;

    while (Fmax >= force_tol_ && step < max_steps_ && h > 1e-12) {
        ++step;
        // trial positions on a copy, so a rejected step leaves f unchanged
        XYZFrame trial = f;
        for (int i = 0; i < f.n(); ++i)
            trial.set_pos(i, trial.Wrap(f.pos(i) + (h / Fmax) * f.force(i)));

        const double U_trial = trial.PotentialEnergy(ff_);
        if (U_trial < U) {                      // downhill: accept, try a longer step
            f = std::move(trial);
            U = U_trial;
            Fmax = MaxForce(f);
            h = std::min(1.2 * h, max_step_);
        } else {                                // uphill: reject, shorter step
            h *= 0.2;
        }
    }
    return { Fmax < force_tol_, step, U, Fmax };
}
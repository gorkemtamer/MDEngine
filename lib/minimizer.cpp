// minimizer.cpp
//
// The part common to all minimizers.

#include <cmath>
#include <format>
#include <stdexcept>
#include "minimizer.h"

Minimizer::Minimizer(const ForceField& ff, double force_tol, int max_steps)
    : ff_(ff),
      force_tol_(force_tol),
      max_steps_(max_steps)
{
    if (force_tol <= 0) throw std::invalid_argument(std::format("Minimizer: force tolerance must be positive (got {:g})", force_tol));
    if (max_steps <= 0) throw std::invalid_argument(std::format("Minimizer: max_steps must be positive (got {})", max_steps));
}

double MaxForce(const XYZFrame& f) {
    double m2 = 0.0;
    for (int i = 0; i < f.n(); ++i) m2 = std::fmax(m2, norm2(f.force(i)));
    return std::sqrt(m2);
}
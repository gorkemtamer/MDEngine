// integrator.cpp
//
// The part common to all integrators.

#include <format>
#include <stdexcept>
#include "integrator.h"

Integrator::Integrator(const ForceField& ff, double dt) 
    : ff_(ff), 
      dt_(dt) 
{
    if (dt <= 0) throw std::invalid_argument(std::format("Integrator: time step must be positive (got {:g})", dt));
}

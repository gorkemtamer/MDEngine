// forcefield.cpp
#include <stdexcept>
#include "forcefield.h"

void ForceField::Add(std::unique_ptr<Potential> p) {
    if (!p) throw std::invalid_argument("ForceField::Add: empty potential");
    potentials_.push_back(std::move(p));
}

void ForceField::Compute(XYZFrame& f) const {
    f.clear_forces();
    double U = 0.0;
    for (const auto& p : potentials_) U += p->AddForces(f);
    f.set_potential_energy(U);
}
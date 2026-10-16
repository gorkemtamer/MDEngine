// forcefield.cpp
#include <stdexcept>
#include "forcefield.h"
#include "potentials/coulomb.h"         // only for the temporary CheckNotCollapsing
#include "potentials/lennard_jones.h"

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

// TEMPORARY: to be replaced by the planned distance / step-size checks
// (see ROADMAP.md). dynamic_cast asks "is this potential a Coulomb / a
// LennardJones?"; it returns nullptr if not.
void ForceField::CheckNotCollapsing() const {
    bool coulomb = false, repulsion = false;
    for (const auto& p : potentials_) {
        if (dynamic_cast<const Coulomb*>(p.get()))      coulomb = true;
        if (dynamic_cast<const LennardJones*>(p.get())) repulsion = true;
    }
    if (coulomb && !repulsion)
        throw std::runtime_error("ForceField: Coulomb only system collapses "
                                 "(Coulomb needs a short-range repulsion, e.g. LennardJones)");
}

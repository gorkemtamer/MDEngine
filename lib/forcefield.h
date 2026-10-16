#pragma once
#include <memory>
#include <vector>
#include "potential.h"

// ---------------------------------------------------------------------
// All interactions of a simulation (FF): a list of potentials whose
// forces and energies add up, e.g. Lennard-Jones + a trap.
//     ForceField ff;
//     ff.Add(std::make_unique<LennardJones>());
//     ff.Compute(frame);
// An empty force field gives zero forces (free particles).
// ---------------------------------------------------------------------
class ForceField {
public:
    // Adds a potential; the force field owns it from then on.
    void Add(std::unique_ptr<Potential> p);

    // Clears the forces of f, lets every potential add its forces and
    // energy, and stores the total energy in f (which marks the forces as
    // valid). Errors of the potentials propagate.
    void Compute(XYZFrame& f) const;

    // TEMPORARY safety check (see ROADMAP.md, Loose ends): throws
    // std::runtime_error if the force field contains Coulomb but no
    // short-range repulsion (LennardJones). Such a system collapses:
    // opposite charges fall into each other and the integrator ejects them.
    // Called by VelocityVerletStep: evaluating energies and forces of a
    // fixed configuration with Coulomb alone is still allowed.
    void CheckNotCollapsing() const;

private:
    std::vector<std::unique_ptr<Potential>> potentials_;
};
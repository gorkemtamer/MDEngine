#pragma once
#include <cstdint>
#include <string>
#include "frame.h"
#include "forcefield.h"
#include "rng.h"
#include "minimizer.h"

// ---------------------------------------------------------------------
// One complete MD run as a sequence of stages, in reduced units (see
// reduced_units.md):
//     Simulation sim(std::move(frame), std::move(ff), 42, 3);
//     sim.Minimize({ .force_tol = 1e-3 });
//     sim.AssignVelocities(1.0);
//     sim.Equilibrate({ .dt = 0.005, .steps = 2000, .ensemble = Ensemble::NVT, .T = 1.0 });
//     sim.Produce({ .dt = 0.005, .steps = 10000, .write_every = 10, .dir = OUTPUT_DIR "run/" });
// The stages are called in order by the program; each one continues from
// the frame the previous one left. A pre-production run is a short Produce
// without output.
//
// The simulation owns the frame, the force field and one random generator
// (velocities and thermostats). It builds the integrator (velocity-Verlet),
// the thermostat (Bussi) and the minimizer (steepest descent) from the
// parameters of each stage. Not copyable or movable: the integrator and
// thermostat refer to its members.
// ---------------------------------------------------------------------

enum class Ensemble {NVE, NVT, NPT};

//Energy minimization (SteepestDescent, see minimizers/steepest_descent.h)
struct MinimizeParams {
    double force_tol = 1e-3;     // converged when max |F_i| < force_tol (eps / sigma)
    int max_steps = 10000;
    double initial_step = 0.01;  // sigma
    double max_step = 0.1;       // sigma
};

// One MD stage (equilibration or production).
struct MDParams {
    double dt = 0.005;                  // time step (tau)
    int steps = 1000;                   // number of steps
    Ensemble ensemble = Ensemble::NVE;  // NVE: velocity-Verlet only; NVT: + Bussi; NPT: not yet (barostat, Phase 8)
    double T = 1.0;                     // NVT: target temperature
    double tau = 0.1;                   // NVT: coupling time of the thermostat
    int write_every = 0;                // frames and energies to dir every n-th step (0: nothing written)
    std::string dir = "";               // output directory, created if needed; required if write_every > 0
    int print_every = 0;                // energies to stdout every n-th step (0: nothing printed)
    double max_move = 0.1;              // largest allowed displacement of a particle per step (sigma)
};

class Simulation {
public:
    // seed: of the random generator for velocities and thermostats.
    // removed_dof: degrees of freedom removed from the 3N velocities, so
    // N_f = 3N - removed_dof (temperature, thermostat):
    //   3: total momentum conserved (periodic, or open without external forces)
    //   6: open without walls or trap: AssignVelocities also removes the
    //      angular momentum
    //   0: walls (WCAWall) or a trap (HarmonicTrap)
    // Throws std::invalid_argument for removed_dof not 0, 3 or 6, for 6 with
    // periodic boundaries, or for an empty frame.
    Simulation(XYZFrame frame, ForceField ff, std::uint64_t seed, int removed_dof);
    Simulation(const Simulation&) = delete;
    Simulation& operator=(const Simulation&) = delete;

    // Minimizes the potential energy (positions only). Not converging is
    // not an error: check the result. Throws std::runtime_error if the
    // minimization collapsed, i.e. particles ended closer than the minimum
    // distance (CheckMinDistance): a force field without a minimum, e.g.
    // Coulomb without a short-range repulsion.
    MinimizeResult Minimize(const MinimizeParams& p);

    // Maxwell-Boltzmann velocities at T, momentum removed (and angular
    // momentum for removed_dof = 6). Throws std::invalid_argument for T < 0.
    void AssignVelocities(double T);

    // MD stages; the same loop, named for the reader (and for later
    // differences, e.g. measurements during production).
    // Each step: integrator, then thermostat (NVT); energies and frames
    // written / printed every n-th step, including step 0 (the start).
    // Throws std::invalid_argument for bad parameters (NPT not implemented
    // yet), std::runtime_error if a particle moves more than max_move in one
    // step ("time step too large") or a written frame has too close particles.
    void Equilibrate(const MDParams& p) { RunMD(p, "equilibration"); }
    void Produce(const MDParams& p) { RunMD(p, "production"); }

    const XYZFrame& frame() const { return frame_; }
    const ForceField& force_field() const { return ff_; }
    int dof() const { return 3 * frame_.n() - removed_dof_; }   // N_f
    double Temperature() const { return 2.0 * frame_.KineticEnergy() / dof(); }

private:
    void RunMD(const MDParams& p, const std::string& stage);

    XYZFrame frame_;
    ForceField ff_;
    Rng rng_;
    int removed_dof_;
};
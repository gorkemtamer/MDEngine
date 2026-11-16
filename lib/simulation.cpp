// simulation.cpp
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <vector>
#include "simulation.h"
#include "velocities.h"
#include "xyzio.h"
#include "checkers.h"
#include "integrators/velocity_verlet.h"
#include "thermostats/bussi.h"
#include "minimizers/steepest_descent.h"

Simulation::Simulation(XYZFrame frame, ForceField ff, std::uint64_t seed, int removed_dof)
    : frame_(std::move(frame)),
      ff_(std::move(ff)),
      rng_(seed),
      removed_dof_(removed_dof)
{
    if (removed_dof != 0 && removed_dof != 3 && removed_dof != 6)
        throw std::invalid_argument(std::format("Simulation: removed_dof must be 0, 3 or 6 (got {})", removed_dof));
    if (removed_dof == 6 && frame_.boundary() == Boundary::Periodic)
        throw std::invalid_argument("Simulation: removed_dof = 6 (angular momentum) needs open boundaries");
    if (dof() <= 0)
        throw std::invalid_argument(std::format("Simulation: no degrees of freedom left ({} particles, removed_dof = {})",
                                                frame_.n(), removed_dof));
}

MinimizeResult Simulation::Minimize(const MinimizeParams& p) {
    SteepestDescent sd(ff_, p.force_tol, p.max_steps, p.initial_step, p.max_step);
    const MinimizeResult r = sd.Minimize(frame_);
    if (int bad = CheckMinDistance(frame_); bad > 0)
        throw std::runtime_error(std::format("Simulation: minimization collapsed: {} particle pair(s) too close "
                                             "(no minimum? e.g. Coulomb without a short-range repulsion)", bad));
    return r;
}

void Simulation::AssignVelocities(double T) {
    ::AssignVelocities(frame_, T, rng_);        // :: the free function, not this method
    if (removed_dof_ == 6) RemoveAngularMomentum(frame_);
}

void Simulation::RunMD(const MDParams& p, const std::string& stage) {
    if (p.steps < 0) throw std::invalid_argument(std::format("Simulation: {}: steps must not be negative", stage));
    if (p.write_every < 0 || p.print_every < 0)
        throw std::invalid_argument(std::format("Simulation: {}: write_every and print_every must not be negative", stage));
    if (p.write_every > 0 && p.dir.empty())
        throw std::invalid_argument(std::format("Simulation: {}: write_every > 0 needs an output directory", stage));
    if (p.max_move <= 0) throw std::invalid_argument(std::format("Simulation: {}: max_move must be positive", stage));
    if (p.ensemble == Ensemble::NPT)
        throw std::invalid_argument(std::format("Simulation: {}: NPT is not implemented yet (needs a barostat)", stage));

    VelocityVerlet vv(ff_, p.dt);
    std::unique_ptr<Bussi> th;                  // only for NVT
    if (p.ensemble == Ensemble::NVT) th = std::make_unique<Bussi>(p.T, p.tau, rng_, removed_dof_);

    std::ofstream log;
    if (p.write_every > 0) {
        std::filesystem::create_directories(p.dir);
        log.open(p.dir + "energies.dat");
        if (!log) throw std::runtime_error("Simulation: cannot open " + p.dir + "energies.dat");
        log << std::format("# {}: dt = {:g}, N_f = {}\n", stage, p.dt, dof());
    }
    const std::string header = std::format("{:>6} {:>8} {:>10} {:>14} {:>14} {:>14}", "step", "t", "T", "Ekin", "Epot", "H");
    if (log.is_open()) log << '#' << header.substr(1) << '\n';
    if (p.print_every > 0) std::printf("%s\n", header.c_str());

    // energies of the current frame: one row of the log / the printout
    auto report = [&](int step) {
        const bool write = p.write_every > 0 && step % p.write_every == 0;
        const bool print = p.print_every > 0 && step % p.print_every == 0;
        if (!write && !print) return;
        const double K = frame_.KineticEnergy();
        const double U = frame_.PotentialEnergy(ff_);
        const double H = K + U - (th ? th->energy_added() : 0.0);   // conserved: drift = integration error
        const std::string row = std::format("{:6d} {:8.3f} {:10.6f} {:14.6f} {:14.6f} {:14.6f}",
                                            step, step * p.dt, 2.0 * K / dof(), K, U, H);
        if (print) std::printf("%s\n", row.c_str());
        if (write) {
            log << row << '\n';
            frame_.comment = std::format("{} step {}, t = {:g}", stage, step, step * p.dt);
            xyz_write(p.dir + std::format("step{:06d}.xyz", step), frame_);
        }
    };

    report(0);
    std::vector<Vec3> old(frame_.n());
    for (int step = 1; step <= p.steps; ++step) {
        for (int i = 0; i < frame_.n(); ++i) old[i] = frame_.pos(i);
        vv.Step(frame_);
        if (th) th->Apply(frame_, p.dt);

        // a particle jumping far in one step means the time step is too large
        // (e.g. ejected from a collision); MinImage undoes the periodic wrap
        for (int i = 0; i < frame_.n(); ++i)
            if (norm2(frame_.MinImage(frame_.pos(i) - old[i])) > p.max_move * p.max_move)
                throw std::runtime_error(std::format("Simulation: {} step {}: particle {} moved more than {:g} sigma "
                                                     "(time step too large?)", stage, step, i, p.max_move));
        report(step);
    }
}
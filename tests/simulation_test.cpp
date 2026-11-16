// simulation_test.cpp
//
// The Simulation class only wires existing parts together, so it is tested
// against the same parts called by hand:
//   1. NVE stage: frame bit-identical to a plain velocity-Verlet loop
//   2. AssignVelocities + NVT stage: bit-identical to the same calls by hand
//      with the same seed (one generator, velocities first, then the thermostat)
//   3. Minimize: same result as SteepestDescent called directly
//   4. AssignVelocities: P = 0; L = 0 for removed_dof = 6; N_f = 3N - removed_dof
//   5. output: frames and energies every write_every steps, step 0 included
//   6. a too large time step is stopped by the step-size check
//   6b. collapse without a minimum (Na+ and Cl- with Coulomb alone): stopped
//       by the step-size check in MD, by the distance check after Minimize
//   7. bad parameters rejected (removed_dof, NPT, write_every without dir, ...)

#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include "simulation.h"
#include "velocities.h"
#include "paths.h"
#include "integrators/velocity_verlet.h"
#include "thermostats/bussi.h"
#include "minimizers/steepest_descent.h"
#include "potentials/coulomb.h"
#include "test_helpers.h"

// true if both frames have bit-identical positions and velocities
static bool same_state(const XYZFrame& a, const XYZFrame& b) {
    for (int i = 0; i < a.n(); ++i)
        if (norm(a.pos(i) - b.pos(i)) != 0.0 || norm(a.vel(i) - b.vel(i)) != 0.0) return false;
    return a.n() == b.n();
}

// open Ar cluster: 3x3x3 lattice with spacing 1.12, perturbed, in a 20^3 box
static XYZFrame cluster(Rng& rng) {
    XYZFrame f(27, Vec3{ 20, 20, 20 }, Boundary::Open);
    int i = 0;
    for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
            for (int z = 0; z < 3; ++z, ++i) {
                f.set_element(i, element_index(i % 4 == 0 ? "Kr" : "Ar"));
                f.set_pos(i, Vec3{ 9.0 + 1.12 * x, 9.0 + 1.12 * y, 9.0 + 1.12 * z }
                                 + 0.05 * Vec3{ rng.uniform() - 0.5, rng.uniform() - 0.5, rng.uniform() - 0.5 });
            }
    return f;
}

int main() {
    try {
        // 1. NVE
        {
            Rng r(5);
            XYZFrame start = LatticeFrame(5, 1.2, 1.0, r);
            XYZFrame by_hand = start;
            Simulation sim(start, LJForceField(), 1, 3);
            sim.Produce({ .dt = 0.004, .steps = 200 });

            const ForceField ff = LJForceField();
            VelocityVerlet vv(ff, 0.004);
            for (int s = 0; s < 200; ++s) vv.Step(by_hand);
            check(same_state(sim.frame(), by_hand), "NVE stage = 200 velocity-Verlet steps (bit-identical)");
        }

        // 2. velocities + NVT on a free cluster (removed_dof = 6): the simulation's
        //    generator is seeded with 7 and used first for the velocities, then by
        //    the thermostat; the hand-made run does the same with its own Rng(7)
        {
            Rng r(5);
            XYZFrame start = cluster(r);
            XYZFrame by_hand = start;
            Simulation sim(start, LJForceField(), 7, 6);
            sim.AssignVelocities(0.5);
            sim.Equilibrate({ .dt = 0.004, .steps = 200, .ensemble = Ensemble::NVT, .T = 0.3, .tau = 0.05 });

            const ForceField ff = LJForceField();
            Rng rng(7);
            AssignVelocities(by_hand, 0.5, rng);
            RemoveAngularMomentum(by_hand);
            VelocityVerlet vv(ff, 0.004);
            Bussi th(0.3, 0.05, rng, 6);
            for (int s = 0; s < 200; ++s) { vv.Step(by_hand); th.Apply(by_hand, 0.004); }
            check(same_state(sim.frame(), by_hand), "velocities + NVT stage = by hand (bit-identical)");
        }

        // 3. Minimize
        {
            Rng r(3);
            XYZFrame start = cluster(r);
            XYZFrame by_hand = start;
            Simulation sim(start, LJForceField(), 1, 6);
            const MinimizeResult a = sim.Minimize({ .force_tol = 1e-4, .max_steps = 5000, .initial_step = 0.02, .max_step = 0.05 });

            const ForceField ff = LJForceField();
            SteepestDescent sd(ff, 1e-4, 5000, 0.02, 0.05);
            const MinimizeResult b = sd.Minimize(by_hand);
            check(a.converged && a.steps == b.steps && a.energy == b.energy && same_state(sim.frame(), by_hand),
                  "Minimize = SteepestDescent (converged, same steps, energy, positions)", a.max_force);
        }

        // 4. velocities and degrees of freedom
        {
            Rng r(3);
            Simulation free_cluster(cluster(r), LJForceField(), 2, 6);
            free_cluster.AssignVelocities(1.0);
            const double P = norm(free_cluster.frame().TotalMomentum());
            const double L = norm(free_cluster.frame().AngularMomentum());
            check(P < 1e-12 && L < 1e-12, "removed_dof = 6: P = 0 and L = 0 after AssignVelocities", P + L);

            Simulation periodic(LatticeFrame(5, 1.2, 1.0, r), LJForceField(), 2, 3);
            periodic.AssignVelocities(1.0);
            check(norm(periodic.frame().TotalMomentum()) < 1e-12 && periodic.dof() == 3 * 125 - 3,
                  "removed_dof = 3: P = 0, N_f = 3N - 3");
            const double T = periodic.Temperature();
            check(std::fabs(T - 2.0 * periodic.frame().KineticEnergy() / 372) < 1e-15, "Temperature() = 2K / N_f", T);
        }

        // 5. output files: 20 steps, written every 5th: steps 0, 5, 10, 15, 20
        {
            const std::string dir = std::string(OUTPUT_DIR) + "simulation_test/";
            std::filesystem::remove_all(dir);
            Rng r(5);
            Simulation sim(LatticeFrame(5, 1.2, 1.0, r), LJForceField(), 1, 3);
            sim.Produce({ .dt = 0.005, .steps = 20, .write_every = 5, .dir = dir });

            int frames = 0;
            for (const auto& e : std::filesystem::directory_iterator(dir))
                if (e.path().extension() == ".xyz") ++frames;
            const bool all = std::filesystem::exists(dir + "step000000.xyz") && std::filesystem::exists(dir + "step000020.xyz");
            check(frames == 5 && all, "5 frames written (steps 0 to 20, every 5th)", frames);

            std::ifstream log(dir + "energies.dat");
            std::string line;
            int comments = 0, rows = 0;
            while (std::getline(log, line)) (line[0] == '#' ? comments : rows)++;
            check(comments == 2 && rows == 5, "energies.dat: 2 header lines, 5 rows", rows);

            XYZFrame last = xyz_read(dir + "step000020.xyz");
            check(norm(last.pos(7) - sim.frame().pos(7)) < 1e-6, "last frame on disk = final state (6 digits)",
                  norm(last.pos(7) - sim.frame().pos(7)));
        }

        // 6. step-size check: two Ar atoms far inside the repulsive core, large dt
        {
            XYZFrame pair(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            pair.set_element(0, element_index("Ar"));
            pair.set_element(1, element_index("Ar"));
            pair.set_pos(0, Vec3{ 5, 5, 5 });
            pair.set_pos(1, Vec3{ 5.8, 5, 5 });               // F = 760: first step moves ~0.15 sigma
            Simulation sim(pair, LJForceField(), 1, 3);
            bool thrown = false;
            try { sim.Produce({ .dt = 0.02, .steps = 10 }); }
            catch (const std::runtime_error& e) { thrown = std::string(e.what()).find("time step too large") != std::string::npos; }
            check(thrown, "time step too large: particle moves > 0.1 sigma, stopped");

            Simulation calm(pair, LJForceField(), 1, 3);
            thrown = false;
            try { calm.Produce({ .dt = 0.001, .steps = 10 }); } catch (const std::runtime_error&) { thrown = true; }
            check(!thrown, "same pair with dt = 0.001: no error");
        }

        // 6b. Coulomb alone: opposite charges fall into each other
        {
            XYZFrame ions(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            ions.set_element(0, element_index("Na"));
            ions.set_element(1, element_index("Cl"));
            ions.set_pos(0, Vec3{ 4, 5, 5 });
            ions.set_pos(1, Vec3{ 6, 5, 5 });
            auto coulomb_only = [] { ForceField ff; ff.Add(std::make_unique<Coulomb>()); return ff; };
            auto message_has = [](auto run, const char* text) {
                try { run(); } catch (const std::runtime_error& e) { return std::string(e.what()).find(text) != std::string::npos; }
                return false;
            };
            Simulation md(ions, coulomb_only(), 1, 3);
            check(message_has([&] { md.Produce({ .dt = 0.001, .steps = 100000 }); }, "time step too large"),
                  "Coulomb alone, MD: stopped by the step-size check");
            Simulation mini(ions, coulomb_only(), 1, 3);
            check(message_has([&] { mini.Minimize({ .force_tol = 1e-3, .max_steps = 100000 }); }, "collapsed"),
                  "Coulomb alone, Minimize: collapse detected");
        }

        // 7. rejections
        {
            auto throws = [](auto run) {
                try { run(); } catch (const std::invalid_argument&) { return true; }
                return false;
            };
            Rng r(5);
            const XYZFrame lat = LatticeFrame(5, 1.2, 1.0, r);
            check(throws([&] { Simulation(lat, LJForceField(), 1, 2); }), "removed_dof = 2 rejected");
            check(throws([&] { Simulation(lat, LJForceField(), 1, 6); }), "removed_dof = 6 with periodic boundaries rejected");
            XYZFrame two(2, Vec3{ 10, 10, 10 }, Boundary::Open);
            check(throws([&] { Simulation(two, LJForceField(), 1, 6); }), "no degrees of freedom left (N = 2, removed 6) rejected");

            Simulation sim(lat, LJForceField(), 1, 3);
            check(throws([&] { sim.Equilibrate({ .ensemble = Ensemble::NPT }); }), "NPT rejected (not implemented)");
            check(throws([&] { sim.Produce({ .steps = -1 }); }), "negative steps rejected");
            check(throws([&] { sim.Produce({ .write_every = 5 }); }), "write_every without a directory rejected");
            check(throws([&] { sim.Produce({ .dt = 0.0 }); }), "dt = 0 rejected");
            check(throws([&] { sim.Produce({ .max_move = 0.0 }); }), "max_move = 0 rejected");
            check(throws([&] { sim.AssignVelocities(-1.0); }), "T < 0 rejected");
        }
    } catch (const std::exception& e) {
        std::printf("FAIL unexpected exception: %s\n", e.what());
        return 1;
    }
    return test_result();
}

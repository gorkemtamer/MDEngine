Personal project for a simple MDEngine.

--will be expanded--

build instructions: (check makefile)
 - make            (builds everything into bin/)
 - bin/main
 - bin/viewer [file.xyz]
 - bin/xyzio_test

layout:
 - src/     programs (main.c, viewer.c)
 - lib/     shared code (xyzio, rng, checkers, elements, paths)
 - bin/     executables, created by make
 - input/   input .xyz files
 - output/  generated files

units: everything is in reduced (Lennard-Jones) units, see reduced_units.md


Starting point: Lennard-Jones potential and hard spheres.

Time integrator: Velocity-verlet

-----------------------------------------------------
7.10.2026:

- First version of the viewer for a single snapshot (mostly done with the help of claude)
- Random number generator for uniform (xoshiro256*), normal (Box-Muller method) and Maxwell-Boltzmann distirbutions for a randomized initial configuration
- Simple implementation of atomic radii to avoid overlapping hard spheres by CheckMinDistance and a 
reduced atomic radii list (normalized to Ar radius)
- xyzio readers and writers (test for read, write and random also implemented)
- simplified build through makefile
- due to Maxwell-Boltzmann distribution mass and temperature are also introduced prematurely
- naive implemented Velocity-Verlet integrator
hardcoded to a LJ potential (epot is computed directly from dU/dr analytically)
-----------------------------------------------------
12.10.2026:

- Switch to cpp in order to avoid overhead in `struct`.
- Vec3 helper functions to collect operations to one place
- for PBC: nearest image distance
- current implementation mostly implemented assuming PBC, OBC will come later
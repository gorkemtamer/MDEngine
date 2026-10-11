#ifndef PARAMS_H
#define PARAMS_H

// ---------------------------------------------------------------------
// Default simulation parameters, in reduced units (see reduced_units.md).
// Placeholders until they can be set from an input file.
// ---------------------------------------------------------------------

// Temperature T* = k_B T / epsilon_Ar.  1.0 corresponds to 124 K for argon.
#define DEFAULT_T 1.0

// Lennard-Jones cutoff in units of sigma_ij: pairs farther apart than
// LJ_CUTOFF * sigma_ij do not interact. 2.5 is the usual choice.
// With periodic boundaries the box must be at least 2 * cutoff wide.
#define LJ_CUTOFF 2.5

#endif

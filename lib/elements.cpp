// elements.c
//
// Per-element data: Lennard-Jones sigma (minimum-distance check, viewer
// sphere size, forces), Lennard-Jones epsilon (forces) and atomic mass. The table keeps literature units
// (Angstrom, Kelvin, atomic mass units u) so values are easy to compare,
// and the functions convert to reduced units with argon as the reference:
// sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1. See reduced_units.md for the full bookkeeping.
//
// Values are approximate:
//   sigma, epsilon, noble gases - Hirschfelder, Curtiss & Bird (viscosity fits)
//   sigma, epsilon, H, C, N, O  - typical OPLS-AA values (depend on chemistry!)
//   mass                - IUPAC standard atomic weights
// Same elements as the color table in viewer.c; add new ones here.

#include <string.h>
#include "elements.h"

#define SIGMA_REF 3.418     // sigma of Ar in Angstrom = 1 reduced length unit
#define EPS_REF   124.0     // epsilon/k_B of Ar in K  = 1 reduced energy unit
#define MASS_REF  39.948    // mass of Ar in u         = 1 reduced mass unit

typedef struct {
    const char *label;
    double sigma;   // Angstrom
    double eps;     // epsilon / k_B in K
    double mass;    // u
} Element;

static const Element ELEMENTS[] = {
    //  label  sigma   eps      mass
    { "H",  2.500,  15.1,    1.008  },
    { "C",  3.500,  33.2,   12.011  },
    { "N",  3.250,  85.5,   14.007  },
    { "O",  2.960, 105.7,   15.999  },
    { "He", 2.576,  10.22,   4.0026 },
    { "Ne", 2.789,  35.7,   20.180  },
    { "Ar", 3.418, 124.0,   39.948  },
    { "Kr", 3.610, 190.0,   83.798  },
    { "Xe", 4.055, 229.0,  131.293  },
};

static const Element *find(const char *label) {
    for (size_t k = 0; k < sizeof(ELEMENTS) / sizeof(ELEMENTS[0]); ++k)
        if (strcmp(ELEMENTS[k].label, label) == 0) return &ELEMENTS[k];
    return NULL;
}

double element_radius(const char *label) {
    const Element *e = find(label);
    return e ? 0.5 * e->sigma / SIGMA_REF : -1.0;
}

double element_mass(const char *label) {
    const Element *e = find(label);
    return e ? e->mass / MASS_REF : -1.0;
}

double element_epsilon(const char *label) {
    const Element *e = find(label);
    return e ? e->eps / EPS_REF : -1.0;
}

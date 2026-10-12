// elements.cpp
//
// Per-element data: Lennard-Jones sigma (minimum-distance check, viewer
// sphere size, forces), Lennard-Jones epsilon (forces) and atomic mass.
// The source table keeps literature units (Angstrom, Kelvin, atomic mass
// units u) so values are easy to compare; it is converted once, at compile
// time, to reduced units with argon as the reference:
// sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1. See reduced_units.md for the full bookkeeping.
//
// Values are approximate:
//   sigma, epsilon, noble gases - Hirschfelder, Curtiss & Bird (viscosity fits)
//   sigma, epsilon, H, C, N, O  - typical OPLS-AA values (depend on chemistry!)
//   mass                - IUPAC standard atomic weights
// Same elements as the color table in viewer.cpp; add new ones here.

#include <array>
#include <cassert>
#include <cstring>
#include "elements.h"

constexpr double SIGMA_REF = 3.418;     // sigma of Ar in Angstrom = 1 reduced length unit
constexpr double EPS_REF   = 124.0;     // epsilon/k_B of Ar in K  = 1 reduced energy unit
constexpr double MASS_REF  = 39.948;    // mass of Ar in u         = 1 reduced mass unit

namespace {

// literature units, as in the papers
struct ElementSource {
    const char *label;
    double sigma;   // Angstrom
    double eps;     // epsilon / k_B in K
    double mass;    // u
};

constexpr ElementSource SOURCE[] = {
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

constexpr int N_ELEMENTS = sizeof(SOURCE) / sizeof(SOURCE[0]);

// the same table in reduced units, computed by the compiler
constexpr std::array<Element, N_ELEMENTS> make_table() {
    std::array<Element, N_ELEMENTS> t{};
    for (int k = 0; k < N_ELEMENTS; ++k)
        t[k] = { SOURCE[k].label,
                 SOURCE[k].sigma / SIGMA_REF,
                 SOURCE[k].eps   / EPS_REF,
                 SOURCE[k].mass  / MASS_REF };
    return t;
}

constexpr std::array<Element, N_ELEMENTS> TABLE = make_table();

} // namespace

int element_index(const char *label) {
    for (int k = 0; k < N_ELEMENTS; ++k)
        if (std::strcmp(TABLE[k].label, label) == 0) return k;
    return -1;
}

const Element& element(int index) {
    assert(index >= 0 && index < N_ELEMENTS);
    return TABLE[index];
}

int element_count() {
    return N_ELEMENTS;
}

double element_radius(const char *label) {
    int k = element_index(label);
    return k >= 0 ? 0.5 * TABLE[k].sigma : -1.0;
}

double element_mass(const char *label) {
    int k = element_index(label);
    return k >= 0 ? TABLE[k].mass : -1.0;
}

double element_epsilon(const char *label) {
    int k = element_index(label);
    return k >= 0 ? TABLE[k].eps : -1.0;
}

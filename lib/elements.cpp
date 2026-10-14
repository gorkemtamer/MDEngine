// elements.cpp
//
// Per-element data: Lennard-Jones sigma (minimum-distance check, viewer
// sphere size, forces), Lennard-Jones epsilon (forces), atomic mass and
// charge (Coulomb).
// The source table keeps literature units (Angstrom, Kelvin, atomic mass
// units u) so values are easy to compare; it is converted once, at compile
// time, to reduced units with argon as the reference (units.h):
// sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1. See reduced_units.md for the full bookkeeping.
//
// Values are approximate:
//   sigma, epsilon, noble gases - Hirschfelder, Curtiss & Bird (viscosity fits)
//   sigma, epsilon, H, C, N, O  - typical OPLS-AA values (depend on chemistry!)
//   sigma, epsilon, Na+, Cl-    - Smith & Dang (1994), with
//                                 1 kcal/mol = 503.22 K
//   mass                - IUPAC standard atomic weights
//   charge              - in units of e; 0 for neutral atoms
// New elements go at the end of SOURCE, so existing indices do not change.
// Same elements as the color table in viewer.cpp; add new ones here.

#include <array>
#include <cassert>
#include <cstring>
#include "elements.h"
#include "units.h"

namespace {

// literature units, as in the papers
struct ElementSource {
    const char *label;
    double sigma;   // Angstrom
    double eps;     // epsilon / k_B in K
    double mass;    // u
    double charge;  // e
};

constexpr ElementSource SOURCE[] = {
    //  label  sigma   eps      mass     charge
    { "H",  2.500,  15.1,    1.008,    0.0 },
    { "C",  3.500,  33.2,   12.011,    0.0 },
    { "N",  3.250,  85.5,   14.007,    0.0 },
    { "O",  2.960, 105.7,   15.999,    0.0 },
    { "He", 2.576,  10.22,   4.0026,   0.0 },
    { "Ne", 2.789,  35.7,   20.180,    0.0 },
    { "Ar", 3.418, 124.0,   39.948,    0.0 },
    { "Kr", 3.610, 190.0,   83.798,    0.0 },
    { "Xe", 4.055, 229.0,  131.293,    0.0 },
    { "Na", 2.350,  65.42,  22.990,   +1.0 },   // Na+ ion
    { "Cl", 4.400,  50.32,  35.45,    -1.0 },   // Cl- ion
};

constexpr int N_ELEMENTS = sizeof(SOURCE) / sizeof(SOURCE[0]);

// the same table in reduced units, computed by the compiler
constexpr std::array<Element, N_ELEMENTS> make_table() {
    std::array<Element, N_ELEMENTS> t{};
    for (int k = 0; k < N_ELEMENTS; ++k)
        t[k] = { SOURCE[k].label,
                 SOURCE[k].sigma / SIGMA_REF,
                 SOURCE[k].eps   / EPS_REF,
                 SOURCE[k].mass  / MASS_REF,
                 SOURCE[k].charge };
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

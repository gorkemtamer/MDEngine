#ifndef ELEMENTS_H
#define ELEMENTS_H

// ---------------------------------------------------------------------
// Per-element data in reduced (Lennard-Jones) units with argon as the
// reference: sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1; charges in units of e
// (see units.h, elements.cpp and reduced_units.md).
//
// Look a label up once with element_index(), then use element(index):
//     int k = element_index("Ar");          // -1 if unknown
//     double m = element(k).mass;
//
// radius_i = sigma_i / 2, so radius_i + radius_j = sigma_ij (Lorentz rule).
// Two particles i, j are "too close" when
//     d_ij < MIN_DIST_FACTOR * (radius_i + radius_j)
// ---------------------------------------------------------------------

// Fraction of sigma_ij below which a pair counts as unphysical.
// At r = 0.8 sigma the LJ energy is already ~ +43 epsilon.
inline constexpr double MIN_DIST_FACTOR = 0.8;

// One element, already converted to reduced units.
struct Element {
    const char *label;
    double sigma;   // LJ sigma_i / sigma_Ar
    double eps;     // LJ epsilon_i / epsilon_Ar
    double mass;    // m_i / m_Ar
    double charge;  // q_i / e (0 for neutral atoms)
};

// Index of `label` in the element table, or -1 if the label is unknown.
int element_index(const char *label);

// Element number `index` (0 <= index < element_count()).
const Element& element(int index);

// Number of elements in the table.
int element_count();

#endif

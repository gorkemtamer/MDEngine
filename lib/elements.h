#ifndef ELEMENTS_H
#define ELEMENTS_H

// ---------------------------------------------------------------------
// Per-element data in reduced (Lennard-Jones) units with argon as the
// reference: sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1
// (see elements.cpp and reduced_units.md).
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
#define MIN_DIST_FACTOR 0.8

// One element, already converted to reduced units.
struct Element {
    const char *label;
    double sigma;   // LJ sigma_i / sigma_Ar
    double eps;     // LJ epsilon_i / epsilon_Ar
    double mass;    // m_i / m_Ar
};

// Index of `label` in the element table, or -1 if the label is unknown.
int element_index(const char *label);

// Element number `index` (0 <= index < element_count()).
const Element& element(int index);

// Number of elements in the table.
int element_count();

// ---------------------------------------------------------------------
// Label-based lookups, kept so existing code compiles unchanged.
// Each call searches the table by name; prefer element_index() + element()
// in loops. To be removed once frames store element indices (Phase 3).
// All return -1.0 if the label is unknown.
// ---------------------------------------------------------------------
double element_radius(const char *label);    // sigma_i / 2
double element_mass(const char *label);      // m_i
double element_epsilon(const char *label);   // epsilon_i

#endif

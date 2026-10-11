#ifndef ELEMENTS_H
#define ELEMENTS_H

// ---------------------------------------------------------------------
// Per-element data in reduced (Lennard-Jones) units with argon as the
// reference: sigma_Ar = 1, epsilon_Ar = 1, m_Ar = 1
// (see elements.c and reduced_units.md).
//
// radius_i = sigma_i / 2, so radius_i + radius_j = sigma_ij (Lorentz rule).
// Two particles i, j are "too close" when
//     d_ij < MIN_DIST_FACTOR * (element_radius(label_i) + element_radius(label_j))
// ---------------------------------------------------------------------

// Fraction of sigma_ij below which a pair counts as unphysical.
// At r = 0.8 sigma the LJ energy is already ~ +43 epsilon.
#define MIN_DIST_FACTOR 0.8

// Radius sigma_i / 2 of the element `label` in reduced units,
// or -1.0 if the label is unknown.
double element_radius(const char *label);

// Mass m_i / m_Ar of the element `label` in reduced units,
// or -1.0 if the label is unknown.
double element_mass(const char *label);

// LJ well depth epsilon_i / epsilon_Ar of the element `label` in reduced
// units, or -1.0 if the label is unknown.
double element_epsilon(const char *label);

#endif

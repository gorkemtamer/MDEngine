# Reduced units

MDEngine works in **reduced (Lennard-Jones) units**. Every quantity is stored
as a plain number measured in units built from three reference constants:
a length σ, an energy ε and a mass m. This file records where those
constants come from and lists the values used in the code, for bookkeeping.

## Why reduced units

For the Lennard-Jones potential

    V(r) = 4 ε [ (σ/r)^12 − (σ/r)^6 ]

choosing σ, ε and m as the units makes the equations of motion free of
physical constants: σ = ε = m = 1 and $k_B$ = 1. This has three advantages:

- the numbers in the simulation are of order 1, which avoids round-off
  problems with values like 10⁻²¹ J or 10⁻¹⁰ m;
- one simulation describes every substance with the same potential shape.
  Results convert to argon, krypton, ... afterwards by multiplying with
  the right reference values;
- results can be compared directly with the large body of LJ literature,
  which is written in these units (T\*, ρ\*, P\*, …).

A star (T\*, ρ\*, …) marks a reduced quantity in this file. The code
drops the star: every number in the code and in `.xyz` files is reduced.

## Reference particle: argon

The reference constants are those of argon. These are the classic
Hirschfelder, Curtiss & Bird viscosity-fit values, chosen so that σ and ε
belong to the same fit:

| Base unit | Symbol | Value | Defined in code |
|---|---|---|---|
| length | σ_Ar | 3.418 Å | `SIGMA_REF` in `lib/units.h` |
| energy | ε_Ar | ε/k_B = 124.0 K | `EPS_REF` in `lib/units.h` |
| mass | m_Ar | 39.948 u | `MASS_REF` in `lib/units.h` |
| Boltzmann constant | k_B | 1 | implicit |
| charge | e | 1.602 × 10⁻¹⁹ C | charges are given in units of e (not reduced) |

If you change the reference (e.g. to σ = 3.405 Å, ε/k_B = 119.8 K, another
common Ar set), update the constants in `lib/units.h` and every table below.
The Coulomb constant below is computed from them automatically.

## Derived units

All other units follow from σ, ε and m. SI values are for the argon
reference above.

| Quantity | Reduced unit | SI value | Convert reduced → physical |
|---|---|---|---|
| length | σ | 3.418 × 10⁻¹⁰ m = 3.418 Å | x = x\* · σ |
| energy | ε | 1.712 × 10⁻²¹ J = 0.01069 eV = 1.031 kJ/mol = 0.2464 kcal/mol | E = E\* · ε |
| mass | m | 6.634 × 10⁻²⁶ kg | M = M\* · m |
| temperature | ε / k_B | 124.0 K | T = T\* · 124.0 K |
| time | τ = σ √(m/ε) | 2.128 × 10⁻¹² s = 2.128 ps | t = t\* · τ |
| velocity | σ / τ | 160.6 m/s | v = v\* · σ/τ |
| force | ε / σ | 5.009 × 10⁻¹² N | F = F\* · ε/σ |
| pressure | ε / σ³ | 4.287 × 10⁷ Pa = 42.87 MPa | P = P\* · ε/σ³ |
| number density | 1 / σ³ | 2.504 × 10²⁸ m⁻³ | n = ρ\* / σ³ |
| mass density | m / σ³ | 1661 kg/m³ | ρ = ρ\* · m/σ³ |

Charges and the Coulomb constant: charges stay in units of e, so the
Coulomb energy of two charges q_i, q_j (in e) at distance r\* (in σ) is

    U* = K* · q_i q_j / r*,    K* = e² / (4π ε₀ σ_Ar ε_Ar) = 394.26

(`COULOMB_K` in `lib/units.h`, computed from the SI constants). Two unit
charges at 1 σ therefore interact about 400 times more strongly than two
argon atoms at their LJ minimum. Ionic systems move correspondingly faster:
a Na⁺–Cl⁻ pair vibrates with a period of about 0.036 τ (argon: ~1 τ), so
they need a time step of about 0.001 τ or smaller.

Handy numbers:

| | Reduced | Physical (Ar) |
|---|---|---|
| typical MD time step | Δt\* = 0.005 | 10.6 fs |
| 1 fs | 4.70 × 10⁻⁴ τ | |
| LJ triple point | T\* ≈ 0.69, liquid ρ\* ≈ 0.85 | ≈ 86 K |
| LJ critical point | T\* ≈ 1.32, ρ\* ≈ 0.32 | ≈ 164 K |
| energy at r = 0.8 σ | V ≈ +43 ε | used for `MIN_DIST_FACTOR` |

The triple and critical points are for the full (untruncated) LJ fluid;
a cut-off potential shifts them slightly.

## Elements

σ, ε and mass are written in `lib/elements.cpp` in Å, K and u, and
converted once, at compile time, to reduced units. `element(element_index(label))`
returns an `Element` with `.sigma` = σ_i / σ_Ar, `.eps` = ε_i / ε_Ar,
`.mass` = m_i / m_Ar and `.charge` = q_i / e. For particle i of a frame, `f->element(i)` gives the
same `Element` and `f->mass(i)` its mass. The particle radius used by the
minimum-distance check and the viewer is σ_i / 2.

| Element | σ (Å) | σ\* | radius\* = σ\*/2 | mass (u) | m\* | ε/k_B (K) | ε\* | q (e) |
|---|---|---|---|---|---|---|---|---|
| H  | 2.500 | 0.7314 | 0.3657 | 1.008   | 0.0252 | 15.1  | 0.122 | 0 |
| C  | 3.500 | 1.0240 | 0.5120 | 12.011  | 0.3007 | 33.2  | 0.268 | 0 |
| N  | 3.250 | 0.9508 | 0.4754 | 14.007  | 0.3506 | 85.5  | 0.690 | 0 |
| O  | 2.960 | 0.8660 | 0.4330 | 15.999  | 0.4005 | 105.7 | 0.852 | 0 |
| He | 2.576 | 0.7537 | 0.3768 | 4.0026  | 0.1002 | 10.22 | 0.0824 | 0 |
| Ne | 2.789 | 0.8160 | 0.4080 | 20.180  | 0.5052 | 35.7  | 0.288 | 0 |
| **Ar** | **3.418** | **1** | **0.5** | **39.948** | **1** | **124.0** | **1** | **0** |
| Kr | 3.610 | 1.0562 | 0.5281 | 83.798  | 2.0977 | 190   | 1.532 | 0 |
| Xe | 4.055 | 1.1864 | 0.5932 | 131.293 | 3.2866 | 229   | 1.847 | 0 |
| Na | 2.350 | 0.6875 | 0.3438 | 22.990  | 0.5755 | 65.42 | 0.528 | +1 |
| Cl | 4.400 | 1.2873 | 0.6437 | 35.45   | 0.8874 | 50.32 | 0.406 | −1 |

Sources and caveats:

- **Noble gases:** σ, ε from Hirschfelder, Curtiss & Bird (viscosity fits).
- **H, C, N, O:** OPLS-AA values for aliphatic H, sp³ C, amide N and carbonyl
  O. Converted with 1 kcal/mol = 503.2 K. These atoms' parameters
  depend strongly on their chemical environment, so treat them as rough
  placeholders.
- **Na⁺, Cl⁻:** σ, ε from Smith & Dang (1994), converted with
  1 kcal/mol = 503.22 K. In `.xyz` files they are labelled `Na` and `Cl`.
- **Masses:** IUPAC standard atomic weights.
- **Charges:** per element; 0 for all neutral atoms.

## Mixing rules for unlike pairs

For a pair i–j of different elements, the Lorentz–Berthelot rules apply:

    σ_ij = (σ_i + σ_j) / 2        (Lorentz)
    ε_ij = √(ε_i ε_j)             (Berthelot)

Because a particle's radius is σ_i / 2, the sum of two radii is σ_ij.

## How quantities are stored in the code

| Where | Quantity | Unit |
|---|---|---|
| `XYZFrame::pos(i)` and `.xyz` columns 2–4 | positions | σ |
| `XYZFrame::vel(i)` and `.xyz` columns 5–7 | velocities | σ/τ |
| `XYZFrame::box()` and `Lattice` line | box lengths | σ |
| `xyz_random(..., scale, ...)` | box edge length | σ |
| particle radius `0.5 * f->element(i).sigma` | σ_i / 2 | σ |
| `XYZFrame::mass(i)` | m_i | m_Ar |
| `rng_MaxwellBoltzmann(r, kT, mass)` | kT, mass → velocity component | ε, m_Ar → σ/τ |
| `MIN_DIST_FACTOR` (= 0.8) | minimum pair distance / σ_ij | dimensionless |
| viewer sphere radius | σ_i / 2 | σ |
| `EPS`, `SIGMA` in `src/main.cpp` | LJ parameters | ε, σ (both 1) |
| `DEFAULT_T` in `lib/params.h` (= 1.0) | temperature used by `xyz_random` | ε/k_B (1.0 = 124 K) |
| `LJ_CUTOFF` in `lib/params.h` (= 2.5) | LJ cutoff radius / σ_ij | dimensionless |
| `XYZFrame::force(i)` | forces | ε/σ |
| `XYZFrame::KineticEnergy()`, `XYZFrame::PotentialEnergy()`, `PotentialEnergy(f)` (`lib/measure.h`) | energies | ε |
| `XYZFrame::TotalMomentum()` | total momentum | m_Ar σ/τ |
| `VelocityVerletStep(f, dt)` | time step | τ (typical 0.005 ≈ 10.6 fs) |
| `Element.sigma`, `.eps`, `.mass`, `.charge` (`element(k)`) | σ_i, ε_i, m_i, q_i | σ_Ar, ε_Ar, m_Ar, e |
| `COULOMB_K` in `lib/units.h` (= 394.26) | Coulomb constant e²/(4πε₀σε) | ε σ / e² |

Example: a Maxwell–Boltzmann velocity for argon at 100 K is
`rng_MaxwellBoltzmann(&rng, 100.0 / 124.0, element(element_index("Ar")).mass)`,
i.e. kT\* = 0.806 and m\* = 1.

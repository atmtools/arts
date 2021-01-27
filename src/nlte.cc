/* Copyright (C) 2018
   Richard Larsson
                            
   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */
   
/**
 * @file nlte.h
 * @author Richard Larsson
 * @date 2018-03-07
 * 
 * @brief Deep calculations for NLTE
 */

#include "nlte.h"
#include "interpolation_lagrange.h"

void statistical_equilibrium_equation(MatrixView A,
                                      ConstVectorView Aij,
                                      ConstVectorView Bij,
                                      ConstVectorView Bji,
                                      ConstVectorView Cij,
                                      ConstVectorView Cji,
                                      ConstVectorView Jij,
                                      const ArrayOfIndex& upper,
                                      const ArrayOfIndex& lower) {
  const Index nlines = Aij.nelem();

  A = 0.0;
  for (Index iline = 0; iline < nlines; iline++) {
    const Index i = upper[iline];
    const Index j = lower[iline];

    A(j, j) -= Bji[iline] * Jij[iline] + Cji[iline];
    A(i, i) -= Aij[iline] + Bij[iline] * Jij[iline] + Cij[iline];

    A(j, i) += Aij[iline] + Bij[iline] * Jij[iline] + Cij[iline];
    A(i, j) += Bji[iline] * Jij[iline] + Cji[iline];
  }
}

void dampened_statistical_equilibrium_equation(
    MatrixView A,
    ConstVectorView x,
    ConstVectorView Aij,
    ConstVectorView Bij,
    ConstVectorView Bji,
    ConstVectorView Cij,
    ConstVectorView Cji,
    ConstVectorView Jij,
    ConstVectorView Lambda,
    const ArrayOfIndex& upper,
    const ArrayOfIndex& lower,
    const Numeric& total_number_count) {
  const Index nlines = Aij.nelem();

  A = 0.0;
  for (Index iline = 0; iline < nlines; iline++) {
    const Index i = upper[iline];
    const Index j = lower[iline];

    const Numeric Source =
        total_number_count *
        (x[i] * Aij[iline] / (x[j] * Bji[iline] - x[i] * Bij[iline]));

    A(j, j) -= Bji[iline] * (Jij[iline] - Lambda[iline] * Source) + Cji[iline];
    A(i, i) -= Aij[iline] * (1.0 - Lambda[iline]) +
               Bij[iline] * (Jij[iline] - Lambda[iline] * Source) + Cij[iline];

    A(j, i) += Aij[iline] * (1.0 - Lambda[iline]) +
               Bij[iline] * (Jij[iline] - Lambda[iline] * Source) + Cij[iline];
    A(i, j) += Bji[iline] * (Jij[iline] - Lambda[iline] * Source) + Cji[iline];
  }
}

void set_constant_statistical_equilibrium_matrix(MatrixView A,
                                                 VectorView x,
                                                 const Numeric& sem_ratio,
                                                 const Index row) {
  A(row, joker) = 1.0;
  x[row] = sem_ratio;
}

Vector createAij(const ArrayOfArrayOfAbsorptionLines& abs_lines) {
  // Size of problem
  const Index n = nelem(abs_lines);
  Vector Aij(n);
  
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        Aij[i] = band.A(k); 
        i++;
      }
    }
  }
  return Aij;
}

Vector createBij(const ArrayOfArrayOfAbsorptionLines& abs_lines) {
  extern const Numeric PLANCK_CONST, SPEED_OF_LIGHT;
  const static Numeric c0 =
  2.0 * PLANCK_CONST / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
  
  // Size of problem
  const Index n = nelem(abs_lines);
  Vector Bij(n);
  
  // Base equation for single state:  B21 = A21 c^2 / 2 h f^3  (nb. SI, don't use this without checking your need)
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        Bij[i] = band.A(k) / (c0 * band.F0(k) * band.F0(k) * band.F0(k));
        i++;
      }
    }
  }
  return Bij;
}

Vector createBji(const Vector& Bij, const ArrayOfArrayOfAbsorptionLines& abs_lines) {
  // Size of problem
  const Index n = Bij.nelem();
  Vector Bji(n);

  // Base equation for single state:  B12 = B21 g2 / g1
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines();k++) {
        Bji[i] = Bij[i] * band.g_upp(k) / band.g_low(k);
        i++;
      }
    }
  }
  return Bji;
}

Vector createCji(const Vector& Cij,
                 const ArrayOfArrayOfAbsorptionLines& abs_lines,
                 const Numeric& T) {
  const Index n = nelem(abs_lines);
  Vector Cji(n);
  setCji(Cji, Cij, abs_lines, T);
  return Cji;
}

void setCji(Vector& Cji,
            const Vector& Cij,
            const ArrayOfArrayOfAbsorptionLines& abs_lines,
            const Numeric& T) {
  extern const Numeric PLANCK_CONST, BOLTZMAN_CONST;
  const static Numeric c0 = -PLANCK_CONST / BOLTZMAN_CONST;
  const Numeric constant = c0 / T;

  // Base equation for single state:  C12 = C21 exp(-hf / kT) g2 / g1
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        Cji[i] = Cij[i] * exp(constant * band.F0(k)) * band.g_upp(k) / band.g_low(k);
        i++;
      }
    }
  }
}

void nlte_collision_factorsCalcFromCoeffs(
    Vector& Cij,
    Vector& Cji,
    const ArrayOfArrayOfAbsorptionLines& abs_lines,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfGriddedField1& collision_coefficients,
    const ArrayOfQuantumIdentifier& collision_line_identifiers,
    const SpeciesAuxData& isotopologue_ratios,
    const ConstVectorView vmr,
    const Numeric& T,
    const Numeric& P) {
  extern const Numeric BOLTZMAN_CONST;

  // size of problem
  const Index nspec = abs_species.nelem();
  const Index ntrans = collision_line_identifiers.nelem();

  // reset Cij for summing later
  Cij = 0;

  // For all species
  for (Index i = 0; i < nspec; i++) {
    // Compute the number density noting that free_electrons will behave differently
    const Numeric numden =
        vmr[i] * (abs_species[i][0].SpeciesNameMain() == "free_electrons"
                      ? 1.0
                      : P / (BOLTZMAN_CONST * T));
    
    for (Index j = 0; j < ntrans; j++) {
      Index iline=0;
      for (auto& lines: abs_lines) {
        for (auto& band: lines) {
          const Numeric isot_ratio =
          isotopologue_ratios.getParam(band.Species(), band.Isotopologue())[0]
          .data[0];
          for (Index k=0; k<band.NumLines(); k++) {
            
            const auto& transition = collision_line_identifiers[j];
            const auto& gf1 = collision_coefficients[i][j];
            
            if (Absorption::id_in_line(band, transition, k)) {
              // Standard linear ARTS interpolation
              const FixedLagrangeInterpolation<1> lag(0, T, gf1.get_numeric_grid(0), false);
              const auto itw = interpweights(lag);
              
              Cij[iline] += interp(gf1.data, itw, lag) * numden * isot_ratio;
              iline++;
              break;
            }
            iline++;
          }
        }
      }
    }
  }

  // Compute the reverse
  setCji(Cji, Cij, abs_lines, T);
}

void nlte_positions_in_statistical_equilibrium_matrix(
    ArrayOfIndex& upper,
    ArrayOfIndex& lower,
    const ArrayOfArrayOfAbsorptionLines& abs_lines,
    const EnergyLevelMap& nlte_field) {
  const Index nl = nelem(abs_lines), nq = nlte_field.Levels().nelem();

  upper = ArrayOfIndex(nl, -1);
  lower = ArrayOfIndex(nl, -1);

  Index i=0;
  for (auto& lines: abs_lines) {
    for (const AbsorptionLines& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        for (Index iq = 0; iq < nq; iq++) {
          if (Absorption::id_in_line_lower(band, nlte_field.Levels()[iq], i))
            lower[i] = iq;
          if (Absorption::id_in_line_upper(band, nlte_field.Levels()[iq], i))
            upper[i] = iq;
        }
        i++;
      }
    }
  }
  
  i = 0;
  for (Index il = 0; il < nl; il++)
    if (upper[il] < 0 or lower[il] < 0) i++;
  ARTS_USER_ERROR_IF (i > 1,
        "Must set upper and lower levels completely for all but one level");
}

Index find_first_unique_in_lower(const ArrayOfIndex& upper,
                                 const ArrayOfIndex& lower) noexcept {
  for (const Index& l : lower) {
    if (std::find(upper.cbegin(), upper.cend(), l) == upper.cend())
      return l;
  }
  return upper.nelem() - 1;
}

void check_collision_line_identifiers(const ArrayOfQuantumIdentifier& collision_line_identifiers) {
  auto p = std::find_if(collision_line_identifiers.cbegin(), collision_line_identifiers.cend(), 
                        [spec=collision_line_identifiers.front().Species(), isot=collision_line_identifiers.front().Isotopologue()]
                        (auto& x) {
                          return
                          spec not_eq x.Species() or 
                          isot not_eq x.Isotopologue() or 
                          x.Type() not_eq QuantumIdentifier::TRANSITION;});
  ARTS_USER_ERROR_IF (p not_eq collision_line_identifiers.cend(),
    *p, "\n"
    "does not match the requirements for a line identifier\n"
    "Your list of species is:\n",
    collision_line_identifiers, "\n"
    "This contains more than one isotopologue or it contains some non-transition type identifiers.\n"
    "It will therefore fail in current code.  You can only input transitions, and a single isotopologue.\n")
}

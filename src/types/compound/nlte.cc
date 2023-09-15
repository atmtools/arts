/**
 * @file nlte.h
 * @author Richard Larsson
 * @date 2018-03-07
 * 
 * @brief Deep calculations for NLTE
 */

#include "nlte.h"
#include "arts_constants.h"
#include "interp.h"
#include "quantum_numbers.h"

std::ostream& operator<<(std::ostream& os, const VibrationalEnergyLevels& vib) {
  bool any = false;
  os << '{';
  for (auto& a: vib) {
    if (any) os << ',' << ' ';
    any = true;
    os << std::quoted(var_string(a.first)) << ": " << a.second;
  }
  return os << '}';
}

void statistical_equilibrium_equation(MatrixView A,
                                      const ConstVectorView& Aij,
                                      const ConstVectorView& Bij,
                                      const ConstVectorView& Bji,
                                      const ConstVectorView& Cij,
                                      const ConstVectorView& Cji,
                                      const ConstVectorView& Jij,
                                      const ArrayOfIndex& upper,
                                      const ArrayOfIndex& lower) {
  const Index nlines = Aij.size();

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
    const ConstVectorView& x,
    const ConstVectorView& Aij,
    const ConstVectorView& Bij,
    const ConstVectorView& Bji,
    const ConstVectorView& Cij,
    const ConstVectorView& Cji,
    const ConstVectorView& Jij,
    const ConstVectorView& Lambda,
    const ArrayOfIndex& upper,
    const ArrayOfIndex& lower,
    const Numeric& total_number_count) {
  const Index nlines = Aij.size();

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
  const Index n = size(abs_lines);
  Vector Aij(n);
  
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        Aij[i] = band.lines[k].A; 
        i++;
      }
    }
  }
  return Aij;
}

Vector createBij(const ArrayOfArrayOfAbsorptionLines& abs_lines) {
  constexpr Numeric c0 = 2.0 * Constant::h / Math::pow2(Constant::c);
  
  // Size of problem
  const Index n = size(abs_lines);
  Vector Bij(n);
  
  // Base equation for single state:  B21 = A21 c^2 / 2 h f^3  (nb. SI, don't use this without checking your need)
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        Bij[i] = band.lines[k].A / (c0 * Math::pow3(band.lines[k].F0));
        i++;
      }
    }
  }
  return Bij;
}

Vector createBji(const Vector& Bij, const ArrayOfArrayOfAbsorptionLines& abs_lines) {
  // Size of problem
  const Index n = Bij.size();
  Vector Bji(n);

  // Base equation for single state:  B12 = B21 g2 / g1
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines();k++) {
        Bji[i] = Bij[i] * band.lines[k].gupp / band.lines[k].glow;
        i++;
      }
    }
  }
  return Bji;
}

Vector createCji(const Vector& Cij,
                 const ArrayOfArrayOfAbsorptionLines& abs_lines,
                 const Numeric& T) {
  const Index n = size(abs_lines);
  Vector Cji(n);
  setCji(Cji, Cij, abs_lines, T);
  return Cji;
}

void setCji(Vector& Cji,
            const Vector& Cij,
            const ArrayOfArrayOfAbsorptionLines& abs_lines,
            const Numeric& T) {
  static constexpr Numeric c0 = -Constant::planck_constant / Constant::boltzmann_constant;
  const Numeric constant = c0 / T;

  // Base equation for single state:  C12 = C21 exp(-hf / kT) g2 / g1
  Index i=0;
  for (auto& lines: abs_lines) {
    for (auto& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        Cji[i] = Cij[i] * exp(constant * band.lines[k].F0) * band.lines[k].gupp / band.lines[k].glow;
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
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const ConstVectorView& vmr,
    const Numeric& T,
    const Numeric& P) {
  // size of problem
  const Index nspec = abs_species.size();
  const Index ntrans = collision_line_identifiers.size();

  // reset Cij for summing later
  Cij = 0;

  // For all species
  for (Index i = 0; i < nspec; i++) {
    // Compute the number density noting that free_electrons will behave differently
    const Numeric numden =
        vmr[i] * (abs_species[i].FreeElectrons() ? 1.0 : P / (Constant::k * T));
    
    for (Index j = 0; j < ntrans; j++) {
      Index iline=0;
      for (auto& lines: abs_lines) {
        for (auto& band: lines) {
          const Numeric isot_ratio = isotopologue_ratios[band.Isotopologue()];
          for (Index k=0; k<band.NumLines(); k++) {
            
            const auto& transition = collision_line_identifiers[j];
            const auto& gf1 = collision_coefficients[i][j];
            const Quantum::Number::StateMatch lt(transition, band.lines[k].localquanta, band.quantumidentity);
            
            if (lt == Quantum::Number::StateMatchType::Full) {
              // Standard linear ARTS interpolation
              const FixedLagrangeInterpolation<1> lag(0, T, gf1.get_numeric_grid(0));
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
    const ArrayOfQuantumIdentifier& nlte_qid) {
  const Index nl = size(abs_lines), nq = nlte_qid.size();

  upper = ArrayOfIndex(nl, -1);
  lower = ArrayOfIndex(nl, -1);
  
  Index i=0;
  for (auto& lines: abs_lines) {
    for (const AbsorptionLines& band: lines) {
      for (Index k=0; k<band.NumLines(); k++) {
        for (Index iq = 0; iq < nq; iq++) {
          const Quantum::Number::StateMatch lt(nlte_qid[iq], band.lines[k].localquanta, band.quantumidentity);
          if (lt == Quantum::Number::StateMatchType::Level and lt.low)
            lower[i] = iq;
          if (lt == Quantum::Number::StateMatchType::Level and lt.upp)
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
                                 const ArrayOfIndex& lower) ARTS_NOEXCEPT {
  for (const Index& l : lower) {
    if (std::find(upper.cbegin(), upper.cend(), l) == upper.cend())
      return l;
  }
  return upper.size() - 1;
}

void check_collision_line_identifiers(const ArrayOfQuantumIdentifier& collision_line_identifiers) {
  auto p =
      std::find_if(collision_line_identifiers.cbegin(),
                   collision_line_identifiers.cend(),
                   [isot = collision_line_identifiers.front().Isotopologue()](
                       auto& x) { return isot not_eq x.Isotopologue(); });
  ARTS_USER_ERROR_IF (p not_eq collision_line_identifiers.cend(),
    *p, "\n"
    "does not match the requirements for a line identifier\n"
    "Your list of species is:\n",
    collision_line_identifiers, "\n"
    "This contains more than one isotopologue or it contains some non-transition type identifiers.\n"
    "It will therefore fail in current code.  You can only input transitions, and a single isotopologue.\n")
}

std::pair<Numeric, Numeric>
VibrationalEnergyLevels::lower_upper(const Key &key) const {
  return {at(key.LowerLevel()), at(key.UpperLevel())};
}

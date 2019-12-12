/* Copyright 2018, Richard Larsson
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/**
 * @file linemixing.cc 
 * @author Richard Larsson
 * @date 2018-08-06
 * 
 * @brief Line mixing calculation implementation
 * 
 * This file contains only experimental code to test the F90-routines.
 * It might evolve to replace it at some point but that is beyond present ability
 */

#include "linemixing.h"
#include "abs_species_tags.h"
#include "complex.h"
#include "lin_alg.h"
#include "linefunctions.h"
#include "linescaling.h"
#include "sorting.h"
#include "species_info.h"
#include "wigner_functions.h"

inline Numeric getB0(const SpeciesTag& main) {
  if (main.IsSpecies("CO2"))
    return Conversion::kaycm2freq(0.39021);  // Herzberg 1966
  else if (main.IsSpecies("O2")) {
    if (main.IsIsotopologue("66"))
      return 43100.44276e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
    else if (main.IsIsotopologue("67"))
      return Conversion::kaycm2freq(1.35);
    else if (main.IsIsotopologue("68"))
      return 40707.38657e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
    else if (main.IsIsotopologue("88"))
      return 38313.72938e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
  } else if (main.IsSpecies("O2")) {
    return -1e99;
  } else if (main.IsSpecies("CH4"))
    return Conversion::kaycm2freq(5.2);

  throw std::runtime_error(
      "Error in getB0: Unsupported main species/isotopologue.  Check your broadening data");
}

inline AdiabaticFactor adiabatic_factor(const SpeciesTag& main,
                                        const SpeciesTag& collider) {
  if (main.IsSpecies("CO2") and collider.IsSpecies("N2"))
    return AdiabaticFactor({2.2e-10}, AdiabaticFactor::Type::Hartmann);
  else if (main.IsSpecies("CO2") and collider.IsSpecies("O2"))
    return AdiabaticFactor({2.4e-10}, AdiabaticFactor::Type::Hartmann);
  else if (main.IsSpecies("O2")) {
    return AdiabaticFactor({0.545e-10}, AdiabaticFactor::Type::Hartmann);
  }

  throw std::runtime_error(
      "Error in adiabatic_factor: Unsupported species pairs, main-collider.  Check your broadening data");
}

inline BasisRate basis_rate(const SpeciesTag& main,
                            const SpeciesTag& collider,
                            const Numeric& T,
                            const Numeric& T0) {
  if (main.IsSpecies("CO2") and collider.IsSpecies("N2"))
    return BasisRate({Conversion::hitran2arts_broadening(0.0180) *
                          std::pow(T0 / T, 0.85),       // Hz/Pa
                      0.81 * std::pow(T0 / T, 0.0152),  // unitless
                      0.008},                           // unitless
                     BasisRate::Type::Hartmann);
  else if (main.IsSpecies("CO2") and collider.IsSpecies("O2"))
    return BasisRate({Conversion::hitran2arts_broadening(0.0168) *
                          std::pow(T0 / T, 0.50),       // Hz/Pa
                      0.82 * std::pow(T0 / T, -0.091),  // unitless
                      0.007},                           // unitless
                     BasisRate::Type::Hartmann);
  else if (main.IsSpecies("O2"))
    return BasisRate({-1e99, -1e99, -1e99}, BasisRate::Type::Hartmann);

  throw std::runtime_error(
      "Error in basis_rate: Unsupported species pairs, main-collider.  Check your broadening data");
}

enum class Species { CO2, O2_66 };

Matrix relaxation_matrix_calculations(const AbsorptionLines& band,
                                      const Vector& population,
                                      const SpeciesTag& collider,
                                      const Numeric& collider_vmr,
                                      const Numeric& T,
                                      const Index& size) try {
  const Index n = band.NumLines();
  Matrix W(n, n);

  auto main = SpeciesTag(band.SpeciesName());
  
  const BasisRate br = basis_rate(main, collider, T, band.T0());
  const AdiabaticFactor af = adiabatic_factor(main, collider);
  const Numeric B0 = getB0(main);
  const ArrayOfArrayOfSpeciesTag pseudo_species(
      {ArrayOfSpeciesTag(1, collider), ArrayOfSpeciesTag(1, main)});
  const Vector pseudo_vmrs({1, 0});
  
  Vector vmrs = LineShape::vmrs(pseudo_vmrs, pseudo_species, band.QuantumIdentity(), band.BroadeningSpecies(), band.Self(), band.Bath(), band.LineShapeType());

  Species spec;
  if (main.IsSpecies("CO2"))
    spec = Species::CO2;
  else if (main.IsSpecies("O2") and main.IsIsotopologue("66"))
    spec = Species::O2_66;
  else
    throw "Unsupported species";

#pragma omp parallel for schedule( \
    guided, 1) if (DO_FAST_WIGNER && !arts_omp_in_parallel())
  for (Index i = 0; i < n; i++) {
    // Create a temporary table to allow openmp
    wig_temp_init(2 * int(size));

    auto shape_parameters = band.ShapeParameters(i, T, 1, vmrs);
    W(i, i) = shape_parameters.G0;

    const Numeric& popi = population[i];
    for (Index j = 0; j < n; j++) {
      const Numeric& popj = population[j];

      if (i not_eq j) {
        OffDiagonalElementOutput X;
        switch (spec) {
          case Species::CO2:
            X = OffDiagonalElement::CO2_IR(band.UpperQuantumNumber(i, QuantumNumberType::J),
                                           band.UpperQuantumNumber(j, QuantumNumberType::J),
                                           band.LowerQuantumNumber(i, QuantumNumberType::J),
                                           band.LowerQuantumNumber(j, QuantumNumberType::J),
                                           band.UpperQuantumNumber(i, QuantumNumberType::l2),
                                           band.UpperQuantumNumber(j, QuantumNumberType::l2),
                                           band.LowerQuantumNumber(i, QuantumNumberType::l2),
                                           band.LowerQuantumNumber(j, QuantumNumberType::l2),
                                           popi,
                                           popj,
                                           br,
                                           af,
                                           T,
                                           B0,
                                           main.SpeciesMass(),
                                           collider.SpeciesMass());
            break;
          case Species::O2_66:
            X = OffDiagonalElement::O2_66_MW(band.UpperQuantumNumber(i, QuantumNumberType::J),
                                             band.UpperQuantumNumber(i, QuantumNumberType::N),
                                             band.LowerQuantumNumber(i, QuantumNumberType::J),
                                             band.LowerQuantumNumber(i, QuantumNumberType::N),
                                             band.UpperQuantumNumber(j, QuantumNumberType::J),
                                             band.UpperQuantumNumber(j, QuantumNumberType::N),
                                             band.LowerQuantumNumber(j, QuantumNumberType::J),
                                             band.LowerQuantumNumber(j, QuantumNumberType::N),
                                             popi, popj, T, collider.SpeciesMass());
            break;
          default:
            throw "DEVELOPER BUG:  Add species here as well, the other one is to check if the computations are valid at all...";
        }

        W(i, j) = X.ij;
        W(j, i) = X.ji;
      }
    }

    // Remove the temporary table
    wig_temp_free();
  }

  // rescale by the VMR
  W *= collider_vmr;
  return W;
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *relaxation_matrix_calculations*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *relaxation_matrix_calculations*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

void normalize_relaxation_matrix(Matrix& W,
                                 const Vector& population,
                                 const Vector& /* d0 */,
                                 const AbsorptionLines& band,
                                 const SpeciesAuxData& partition_functions,
                                 const Numeric& T) try {
  const Index n = band.NumLines();

  Vector d = reduced_dipole_vector(band, RedPoleType::ElectricRoVibDipole);

  // Index list showing sorting of W
  ArrayOfIndex sorted(n, 0);
  Vector test(n);
  if (n) {
    const Numeric QT =
        single_partition_function(T,
                                  partition_functions.getParamType(band.QuantumIdentity()),
                                  partition_functions.getParam(band.QuantumIdentity()));
    for (Index i = 0; i < n; i++) {
      const Numeric QT0 =
          single_partition_function(band.T0(),
                                    partition_functions.getParamType(band.QuantumIdentity()),
                                    partition_functions.getParam(band.QuantumIdentity()));
      test[i] = Linefunctions::lte_linestrength(band.I0(i),
                                                band.E0(i),
                                                band.F0(i),
                                                QT0,
                                                band.T0(),
                                                QT,
                                                T);
    }
  }

  get_sorted_indexes(sorted, test);
  for (Index i = 0; i < n - i - 1; i++) std::swap(sorted[i], sorted[n - i - 1]);

  // Sorted matrix
  Matrix Wr(n, n);
  for (Index i = 0; i < n; i++) {
    Wr(i, i) = W(sorted[i], sorted[i]);
    for (Index j = 0; j < n; j++) {
      if (i not_eq j) {
        Wr(i, j) = -std::abs(W(sorted[i], sorted[j]));
      }
    }
  }

  // Renormalization procedure
  for (Index i = 0; i < n; i++) {
    // Sum up upper and lower contributions
    Numeric Sup = 0, Slo = 0;
    for (Index j = 0; j < n; j++) {
      if (j <= i)
        Sup += std::abs(d[sorted[j]]) * Wr(i, j);
      else
        Slo += std::abs(d[sorted[j]]) * Wr(i, j);
    }

    Numeric UL = Sup / Slo;
    if (not std::isnormal(UL)) UL = 1.0;

    // Rescale to fulfill sum-rule, note how the loop for the upper triangle so the last row cannot be renormalized properly
    for (Index j = i; j < n; j++) {
      const Numeric r = population[sorted[i]] / population[sorted[j]];
      Wr(j, i) = r * Wr(i, j);
      if (j not_eq i) Wr(i, j) *= -UL;
    }
  }

  for (Index i = 0; i < n - 1; i++) Wr(n - 1, i) = 0;  // TEST!

  // Backsort matrix before returning
  for (Index i = 0; i < n; i++)
    for (Index j = 0; j < n; j++) W(sorted[i], sorted[j]) = Wr(i, j);
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *normalize_relaxation_matrix*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *normalize_relaxation_matrix*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

inline Matrix hartmann_relaxation_matrix(
    const AbsorptionLines& band,
    const Vector& population,
    const Vector& d0,
    const ArrayOfSpeciesTag& colliders,
    const Vector& colliders_vmr,
    const SpeciesAuxData& partition_functions,
    const Numeric& T,
    const Index& size) try {
  // Size of problem
  const Index n = band.NumLines();
  const Index c = colliders.nelem();

  // Create and normalize the matrix
  Matrix W(n, n, 0);
  if (n) {
    for (Index ic = 0; ic < c; ic++)
      W += relaxation_matrix_calculations(band,
                                          population,
                                          colliders[ic],
                                          colliders_vmr[ic],
                                          T,
                                          size);
    normalize_relaxation_matrix(
        W, population, d0, band, partition_functions, T);
  }

  return W;
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *hartmann_relaxation_matrix*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *hartmann_relaxation_matrix*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

inline Numeric population_density(Numeric T,
                                  Numeric E0,
                                  Numeric F0,
                                  Numeric QT) {
  return (1.0 - stimulated_emission(T, F0)) * boltzman_factor(T, E0) / QT;
}

inline Numeric dpopulation_densitydT(
    Numeric T, Numeric E0, Numeric F0, Numeric QT, Numeric dQTdT) {
  using namespace Constant;
  return (1.0 - dstimulated_emissiondT(T, F0)) * boltzman_factor(T, E0) / QT +
         (1.0 - stimulated_emission(T, F0)) * dboltzman_factordT(T, E0) / QT +
         (1.0 - stimulated_emission(T, F0)) * boltzman_factor(T, E0) *
             (-dQTdT) / pow2(QT);
}

Vector population_density_vector(const AbsorptionLines& band,
                                 const SpeciesAuxData& partition_functions,
                                 const Numeric& T) try {
  auto n = band.NumLines();
  Vector p(n);

  if (n) {
    const Numeric QT = single_partition_function(T,
        partition_functions.getParamType(band.QuantumIdentity()),
        partition_functions.getParam(band.QuantumIdentity()));

    for (auto k = 0; k < n; k++) {
      p[k] = population_density(T, band.E0(k), band.F0(k), QT);
    }
  }

  return p;
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *population_density_vector*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *population_density_vector*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

Vector dipole_vector(const AbsorptionLines& band,
                     const SpeciesAuxData& partition_functions) try {
  const Index n = band.NumLines();
  Vector d(n, 0);
  if (not n) return d;

  const Numeric T0 = band.T0();
  const Vector p0 =
  population_density_vector(band, partition_functions, T0);
  for (Index k = 0; k < n; k++) {
    d[k] = std::sqrt(band.I0(k) / p0[k]);
  }
  
  return d;
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *dipole_vector*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *dipole_vector*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

/*! Computes the dipole moment
 * \param abs_lines: all absorption lines of the band
 * \return dipole moment of the absorption lines
 */
Vector reduced_dipole_vector(const AbsorptionLines& band,
                             const RedPoleType type) try {
  const Index n = band.NumLines();
  Vector d(n, 0);
  if (not n) return d;

  switch (type) {
    case RedPoleType::ElectricRoVibDipole:
      for (Index i = 0; i < n; i++)
        d[i] = Absorption::reduced_rovibrational_dipole(band.LowerQuantumNumber(i, QuantumNumberType::J),
 band.UpperQuantumNumber(i, QuantumNumberType::J), band.LowerQuantumNumber(i, QuantumNumberType::l2), band.UpperQuantumNumber(i, QuantumNumberType::l2));
      break;
    case RedPoleType::MagneticQuadrapole:
      for (Index i = 0; i < n; i++)
        d[i] = Absorption::reduced_magnetic_quadrapole(band.LowerQuantumNumber(i, QuantumNumberType::J), band.UpperQuantumNumber(i, QuantumNumberType::J), band.UpperQuantumNumber(i, QuantumNumberType::N));
      break;
  }
  return d;
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *reduced_dipole_vector*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *reduced_dipole_vector*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

Vector rosenkranz_first_order(const AbsorptionLines& band,
                              const Matrix& W,
                              const Vector& d0) {
  const Index n = band.NumLines();
  Vector Y(n, 0);

  for (Index i = 0; i < n; i++) {
    Numeric sum = 0;
    for (Index j = 0; j < n; j++) {
      const Numeric df = band.F0(i) - band.F0(j);
      if (std::isnormal(df)) sum += d0[j] / d0[i] * W(i, j) / df;
    }
    Y[i] = -2 * sum;  // Sign change because ARTS uses (1-iY)
  }

  return Y;
}

Vector rosenkranz_shifting_second_order(const AbsorptionLines& band,
                                        const Matrix& W) {
  const Index n = band.NumLines();
  Vector DV(n, 0);

  for (Index i = 0; i < n; i++) {
    Numeric sum = 0;
    for (Index j = 0; j < n; j++) {
      const Numeric df = band.F0(j) - band.F0(i);
      if (std::isnormal(df)) sum += W(i, j) * W(j, i) / df;
    }
    DV[i] = sum;  // Does this require a sign change????
  }

  return DV;
}

Vector rosenkranz_scaling_second_order(const AbsorptionLines& band,
                                       const Matrix& W,
                                       const Vector& d0) {
  const Index n = band.NumLines();
  Vector G(n, 0);

  for (Index i = 0; i < n; i++) {
    const Numeric& di = d0[i];
    Numeric sum1 = 0;
    Numeric sum2 = 0;
    Numeric sum3 = 0;
    Numeric sum4 = 0;
    for (Index j = 0; j < n; j++) {
      const Numeric df = band.F0(j) - band.F0(i);
      if (std::isnormal(df)) {
        const Numeric& dj = d0[j];
        const Numeric r = dj / di;

        sum1 += W(i, j) * W(j, i) / Constant::pow2(df);
        sum2 += r * W(i, j) / df;
        sum3 += r * W(i, j) * W(i, i) / Constant::pow2(df);

        Numeric sum_tmp = 0;
        for (Index k = 0; k < n; k++) {
          const Numeric dfk = band.F0(k) - band.F0(i);
          if (std::isnormal(dfk)) sum_tmp += W(k, j) * W(i, k) / (dfk * df);
        }
        sum4 += r * sum_tmp;
      }
    }
    G[i] = sum1 - Constant::pow2(sum2) + 2 * sum3 - 2 * sum4;
  }
  return G;
}

Matrix hartmann_ecs_interface(const AbsorptionLines& band,
                              const ArrayOfSpeciesTag& collider_species,
                              const Vector& collider_species_vmr,
                              const SpeciesAuxData& partition_functions,
                              const Numeric& T,
                              const Index& size) try {
  const Vector population =
  population_density_vector(band, partition_functions, T);
  const Vector dipole = dipole_vector(band, partition_functions);
  const Matrix W = hartmann_relaxation_matrix(band,
                                              population,
                                              dipole,
                                              collider_species,
                                              collider_species_vmr,
                                              partition_functions,
                                              T,
                                              size);
  return W;
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *hartmann_ecs_interface*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

OffDiagonalElementOutput OffDiagonalElement::CO2_IR(
    const Rational& Jku,
    const Rational& Jju,
    const Rational& Jkl,
    const Rational& Jjl,
    const Rational& l2ku,
    const Rational& l2ju,
    const Rational& l2kl,
    const Rational& l2jl,
    const Numeric& j_rho,
    const Numeric& k_rho,
    const BasisRate& br,
    const AdiabaticFactor& af,
    const Numeric& T,
    const Numeric& B0,
    const Numeric& main_mass,
    const Numeric& collider_mass) {
  const bool jbig = Jjl >= Jkl;

  // Prepare for the wigner calculations --- NOTE: twice the values of J and l2 required for fast Wigner-solver
  // Also, prepare for initial and final phase to go from j to k if Jj >= Jk and vice verse
  const int Ji = (2 * (jbig ? Jju : Jku)).toInt();
  const int Jf = (2 * (jbig ? Jjl : Jkl)).toInt();
  const int Ji_p = (2 * (jbig ? Jku : Jju)).toInt();
  const int Jf_p = (2 * (jbig ? Jkl : Jjl)).toInt();
  const int li = (2 * (jbig ? l2ju : l2ku)).toInt();
  const int lf = (2 * (jbig ? l2jl : l2kl)).toInt();

  // Find best start and end-point of the summing loop
  //   const int st = std::max(Ji - Ji_p, Jf - Jf_p);
  const int en = std::min(Ji + Ji_p, Jf + Jf_p);

  // Adiabatic factor for Ji
  const Numeric AF1 = af.get(Ji / 2, B0, T, main_mass, collider_mass);

  // Scale constant for final state NOTE: "%4" and lack of 2*J because Fast library already doubles these numbers
  const Numeric K1 = ((li + lf) % 4 ? 1.0 : -1.0) * Numeric(Ji_p + 1) *
                     sqrt(Numeric((Jf + 1) * (Jf_p + 1))) * AF1;

  Numeric sum = 0;
  for (int L = 4; L <= en; L += 4) {
    //   for(int L = st>4?st:4; L <= en; L+=4) {
    // Basis rate for L
    const Numeric QL = br.get(L / 2, B0, T);

    // Adiabatic factor for L
    const Numeric AF2 = af.get(L / 2, B0, T, main_mass, collider_mass);

    // The wigner-symbol following Niro etal 2004
    const Numeric y = co2_ecs_wigner_symbol(Ji, Jf, Ji_p, Jf_p, L, li, lf);

    // Sum to the total
    sum += QL * y / AF2;
  }

  sum *= K1;

  const Numeric r = k_rho / j_rho;
  return {jbig ? sum : sum * r, jbig ? sum / r : sum};
}

constexpr auto params = 10;

void linearized_relaxation_matrix(Tensor3& M,
                                  const ArrayOfRational& Ji,
                                  const ArrayOfRational& Jf,
                                  const ArrayOfRational& l2i,
                                  const ArrayOfRational& l2f,
                                  [[maybe_unused]] const Vector& F0,
                                  const Vector& d,
                                  const Vector& rho,
                                  const Numeric& T,
                                  const Vector& a = Vector(params, 1)) {
  const Index n = Ji.nelem();
  Numeric c[params];

  Rational rmax = Ji[0];
  for (auto& r : Ji) rmax = max(r, rmax);
  for (auto& r : Jf) rmax = max(r, rmax);

  M = 0;

#pragma omp parallel for schedule( \
    guided, 1) if (DO_FAST_WIGNER && !arts_omp_in_parallel())
  for (Index i = 0; i < n; i++) {
    wig_temp_init(int(2 * rmax));
    for (Index j = 0; j < n; j++) {
      if (i == j) continue;         // to next loop
      if (Ji[i] < Ji[j]) continue;  // because we renormalize this

      const Numeric K1 = ((l2i[i] + l2f[i]) % 2 ? 1 : -1) *
                         Numeric(2 * Ji[j] + 1) *
                         sqrt((2 * Jf[i] + 1) * (2 * Jf[j] + 1));
      const Numeric d_ratio = d[j] / d[i];
      [[maybe_unused]] const Numeric exp =
          (Ji[i] == 0) ? Numeric(Ji[j] / Ji[i]) : Numeric(Ji[j]) / 0.5;

      const int en = std::min(int(Ji[i] * 2) + int(Ji[j] * 2),
                              int(Jf[i] * 2) + int(Jf[j] * 2));
      for (int L = 4; L <= en; L += 4) {
        const Numeric K2 = co2_ecs_wigner_symbol(int(2 * Ji[i]),
                                                 int(2 * Jf[i]),
                                                 int(2 * Ji[j]),
                                                 int(2 * Jf[j]),
                                                 L,
                                                 int(2 * l2i[i]),
                                                 int(2 * l2f[i]));

        const Numeric dE = L / 2 * (L / 2 + 1);
        c[0] = 1;
        c[1] = 1 / T;
        c[2] = dE;
        c[3] = dE / T;
        c[4] = dE * dE;
        c[5] = std::log(dE);
        c[6] = std::log(dE) / T;
        c[7] = std::log(dE) * std::log(dE);
        c[8] = 1 / dE;
        c[9] = 1 / dE / dE;

        for (auto kk = 0; kk < params; kk++)
          M(i, j, kk) += a[kk] * K2 * (std::isnormal(c[kk]) ? c[kk] : 0);
      }

      if (std::isnormal(d_ratio))
        M(i, j, joker) *= d_ratio * K1;
      else
        M(i, j, joker) *= K1;
      M(j, i, joker) = M(i, j, joker);

      M(i, j, joker) *= rho[j] / rho[i];
      M(j, i, joker) *= rho[i] / rho[j];
    }

    // Remove the temporary table
    wig_temp_free();
  }
}

Matrix CO2_ir_training(const ArrayOfRational& Ji,
                       const ArrayOfRational& Jf,
                       const ArrayOfRational& l2i,
                       const ArrayOfRational& l2f,
                       const Vector& F0,
                       const Vector& d,
                       const Vector& rho,
                       const Vector& gamma,
                       const Numeric& T) {
  const auto n = Ji.nelem();

  Vector a(params, 0);
  Matrix A(n, params, 0);
  Matrix W(n, n, 0);
  Tensor3 M(n, n, params, 0);

  // Create the training-tensor
  linearized_relaxation_matrix(M, Ji, Jf, l2i, l2f, F0, d, rho, T);

  // Set the over-determined matrix for least-square-fit
  for (Index i = 0; i < n; i++) {
    for (Index j = 0; j < n; j++) {
      A(i, joker) += M(i, j, joker);
    }
  }

  // Perform least-square-fit
  [[maybe_unused]] const Numeric R2 = lsf(a, A, gamma);

  // Create the proper tensor by recomputing the training tensor with new a-ratios
  linearized_relaxation_matrix(M, Ji, Jf, l2i, l2f, F0, d, rho, T, a);

  // Fill the relaxation matrix
  for (Index i = 0; i < n; i++) {
    for (Index j = 0; j < n; j++) {
      if (i == j)
        W(i, j) = gamma[i];
      else
        W(i, j) = M(i, j, joker).sum();
    }
  }

  return W;
}

Numeric o2_66_inelastic_cross_section_makarov(int L, Numeric T) {
  using namespace Constant;
  using namespace Molecule::O2_66;

  Numeric const1 = 0.086 + 8154e-7 * T;
  static constexpr Numeric const2 = 0.5805;

  return (2 * L + 1) / std::pow(L * L + L, const1) *
         std::exp(-const2 * h * hamiltonian_freq(L) / (k * T));
}

Numeric o2_66_adiabatic_factor_makarov(int L,
                                       Numeric T,
                                       Numeric collider_mass) {
  using namespace Constant;
  using namespace Conversion;
  using namespace Molecule::O2_66;

  static constexpr Numeric const1 = 0.545e-10;
  static constexpr Numeric constant = 2000 * R * inv_pi * pow2(inv_ln_2);

  // Mean speed of collisions
  const Numeric invmu = (1 / mass + 1 / collider_mass);
  const Numeric vm2 = constant * T * invmu;

  return 1 / pow2(1.0 + pow2(freq2angfreq(hamiltonian_freq(L) -
                                          hamiltonian_freq(L - 2)) *
                             const1) /
                            vm2 / 24.0);
}

OffDiagonalElementOutput OffDiagonalElement::O2_66_MW(
    const Rational& J1u,
    const Rational& N1u,
    const Rational& J1l,
    const Rational& N1l,
    const Rational& J2u,
    const Rational& N2u,
    const Rational& J2l,
    const Rational& N2l,
    const Numeric& rho1,
    const Numeric& rho2,
    const Numeric& T,
    const Numeric& collider_mass) {

  // Find which is the 'upper' transition
  const bool onebig =
      Molecule::O2_66::hamiltonian_freq(
          J1u.toNumeric(), (J1u - N1u).toInt(), (J1u - N1u).toInt()) >
      Molecule::O2_66::hamiltonian_freq(
          J2u.toNumeric(), (J2u - N2u).toInt(), (J2u - N2u).toInt());

  // Define the transitions as in Makarov etal 2013, double the value for wiglib
  const int Nk = (2 * (onebig ? N1u : N2u).toInt());
  const int Nkp = (2 * (onebig ? N1l : N2l).toInt());
  const int Jk = (2 * (onebig ? J1u : J2u).toInt());
  const int Jkp = (2 * (onebig ? J1l : J2l).toInt());
  const int Nl = (2 * (onebig ? N2u : N1u).toInt());
  const int Nlp = (2 * (onebig ? N2l : N1l).toInt());
  const int Jl = (2 * (onebig ? J2u : J1u).toInt());
  const int Jlp = (2 * (onebig ? J2l : J1l).toInt());

  if (Nl not_eq Nlp or Nk not_eq Nkp) throw "ERRROR, bad N-values";

  // 'length' of numbers
  const Numeric lNk = std::sqrt(Nk + 1.0);
  const Numeric lNl = std::sqrt(Nl + 1.0);
  const Numeric lJk = std::sqrt(Jk + 1.0);
  const Numeric lJl = std::sqrt(Jl + 1.0);
  const Numeric lJkp = std::sqrt(Jkp + 1.0);
  const Numeric lJlp = std::sqrt(Jlp + 1.0);

  // Constant independenf of L
  const Numeric const1 = lNk * lNl * std::sqrt(lJk * lJl * lJkp * lJlp) *
                         o2_66_inelastic_cross_section_makarov(Nk / 2, T);

  Numeric sum = 0;
  for (int L = 4; L < 400; L += 4) {
    const int sgn = ((Jk + Jl + L + 2) % 4) ? 1 : -1;

    const Numeric const2 =
        sgn * const1 * o2_66_adiabatic_factor_makarov(L / 2, T, collider_mass) /
        o2_66_inelastic_cross_section_makarov(L / 2, T);

    sum += o2_ecs_wigner_symbol(Nl, Nk, Jl, Jk, Jlp, Jkp, L) * const2;
  }

  return {onebig ? sum : sum * rho2 / rho1, onebig ? sum * rho1 / rho2 : sum};
}

Numeric AdiabaticFactor::mol_X(const Numeric& L,
                               const Numeric& B0,
                               const Numeric& T,
                               const Numeric& main_mass,
                               const Numeric& collider_mass) const {
  using namespace Constant;

  using namespace Conversion;

  if (L < 1) return 0.;

  static constexpr Numeric constant = 2000 * R * inv_pi * pow2(inv_ln_2);

  const Numeric invmu = (1 / main_mass + 1 / collider_mass);

  const Numeric vm2 = constant * T * invmu;

  return 1 / pow2(1.0 + pow2(freq2angfreq(B0) * (4 * L - 2) *
                             mdata[Index(HartmannPos::dc)]) /
                            vm2 / 24.0);
}

Numeric BasisRate::mol_X(const Numeric& L,
                         const Numeric& B0,
                         const Numeric& T) const {
  using namespace Constant;

  const Numeric& a1 = mdata[Index(HartmannPos::a1)];
  const Numeric& a2 = mdata[Index(HartmannPos::a2)];
  const Numeric& a3 = mdata[Index(HartmannPos::a3)];

  const Numeric EL = L * L + L;
  return a1 / std::pow(EL, a2) * std::exp(-a3 * h * B0 * EL / (k * T));
}

SecondOrderLineMixingCoeffs compute_2nd_order_lm_coeff(ConstVectorView y,
                                                       ConstVectorView x,
                                                       const Numeric exp,
                                                       const Numeric x0) {
  const auto n = y.nelem();

  if (n not_eq x.nelem())
    throw std::runtime_error("Bad input to compute_2nd_order_lm_coeff");

  Index best_x = 0;
  for (auto i = 0; i < n - 1; i++) {
    if (x[i] >= x[i + 1])
      throw std::runtime_error(
          "Bad structure on x-input; must be constantly increasing");
    if (x[i] > x0) best_x++;
  }

  Matrix A(n, 2);
  Vector ans(2, y[best_x] * 0.5), dans(2, 0), delta(n, 0);

  Numeric rel_res = 1e99;
  Numeric rel_res_limit = 1e-16;
  auto loop_count = 0;
  constexpr auto max_loop_count = 20;

  do {
    for (auto i = 0; i < n; i++) {
      const Numeric theta = x0 / x[i];
      const Numeric TP = pow(theta, exp);
      A(i, 0) = TP;
      A(i, 1) = (theta - 1.0) * TP;
      const Numeric f = ans[0] * TP + ans[1] * (theta - 1.0) * TP;
      delta[i] = y[i] - f;
    }

    // Least square fit
    const Numeric res = lsf(dans, A, delta);

    if (not std::isnormal(res)) break;

    rel_res = std::abs(dans[0] / ans[0]) + std::abs(dans[1] / ans[1]);
    ans += dans;
    loop_count += 1;

  } while (rel_res > rel_res_limit and loop_count < max_loop_count);

  return {ans[0], ans[1]};
}

ComplexVector equivalent_linestrengths(
    const Vector& population,
    const Vector& dipole,
    const Eigen::ComplexEigenSolver<Eigen::MatrixXcd>& M) {
  const auto n = population.nelem();

  auto& V = M.eigenvectors();
  auto Vinv = V.inverse().eval();

  ComplexVector B(n);

  for (auto i = 0; i < n; i++) {
    for (auto j1 = 0; j1 < n; j1++) {
      for (auto j2 = 0; j2 < n; j2++) {
        B[i] +=
            population[j1] * dipole[j1] * dipole[j2] * V(i, j2) * Vinv(j1, i);
      }
    }
  }
  return B;
}

Numeric total_linestrengths(const Vector& population, const Vector& dipole) {
  using namespace Constant;
  const auto n = population.nelem();

  Numeric I0 = 0;
  for (auto i = 0; i < n; i++) I0 += pi * population[i] * pow2(dipole[i]);
  return I0;
}

void relmatInAir(Matrix& relmat,
                 const AbsorptionLines& band,
                 const SpeciesAuxData& partition_functions,
                 const Index& wigner_initialized,
                 const Numeric& temperature) try {
  // Only for Earth's atmosphere
  const ArrayOfSpeciesTag collider_species = {SpeciesTag("O2-66"),
                                              SpeciesTag("N2-44")};
  const Vector collider_species_vmr = {0.21, 0.79};

  
  if (band.Species() == SpeciesTag("CO2-626").Species() and
      band.Isotopologue() == SpeciesTag("CO2-626").Isotopologue()) {
  } else if (band.Species() == SpeciesTag("O2-66").Species() and
             band.Isotopologue() == SpeciesTag("O2-66").Isotopologue()) {
  } else
    throw "Limit in functionality encountered.  We only support CO2-626 and O2-66 for now.";
  
  if (band.BroadeningSpecies().nelem() not_eq 2)
    throw "Bad species length, only works for self+air broadening for now";
  
  if (band.Population() == Absorption::PopulationType::ByRelmatHartmannLTE)
    relmat = hartmann_ecs_interface(band,
                                    collider_species,
                                    collider_species_vmr,
                                    partition_functions,
                                    temperature,
                                    wigner_initialized);
  else
    throw "Population type has not been implemented";
} catch (const char* e) {
  std::ostringstream os;
  os << "Errors raised by *relmatInAir*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
} catch (const std::exception& e) {
  std::ostringstream os;
  os << "Errors in calls by *relmatInAir*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}

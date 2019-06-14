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

// NOTE: This file contains only experimental code to test the F90-routines.
// It might evolve to replace it at some point but that is beyond present intent

#include "abs_species_tags.h"
#include "complex.h"
#include "linemixing.h"
#include "wigner_functions.h"
#include "sorting.h"
#include "linescaling.h"
#include "lin_alg.h"
#include "linefunctions.h"

//! Rotational constants of supported scenarios
/*!
 \param main: species of interest
 \return The rotational constant in Hertz
 */
inline Numeric getB0(const SpeciesTag& main)
{
  if(main.IsSpecies("CO2"))
    return  Conversion::kaycm2freq(0.39021);  // Herzberg 1966
  else if(main.IsSpecies("O2")) {
    if(main.IsIsotopologue("66"))
      return 43100.44276e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
    else if(main.IsIsotopologue("67"))
      return Conversion::kaycm2freq(1.35);
    else if(main.IsIsotopologue("68"))
      return 40707.38657e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
    else if(main.IsIsotopologue("88"))
      return 38313.72938e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
  }
  else if(main.IsSpecies("O2")) {
    return -1e99;
  }
  else if(main.IsSpecies("CH4"))
    return Conversion::kaycm2freq(5.2);
  
  throw std::runtime_error("Error in getB0: Unsupported main species/isotopologue.  Check your broadening data");
}

//! Adiabatic function format for supported scenarios
/*!
 \param main: species of interest
 \param collider: species broadening main
 \return An adiabatic factor function
 */
inline AdiabaticFactor adiabatic_factor(const SpeciesTag& main, const SpeciesTag& collider)
{
  if(main.IsSpecies("CO2") and collider.IsSpecies("N2"))
    return AdiabaticFactor({2.2e-10}, AdiabaticFactor::Type::Hartmann);
  else if(main.IsSpecies("CO2") and collider.IsSpecies("O2"))
    return AdiabaticFactor({2.4e-10}, AdiabaticFactor::Type::Hartmann);
  else if(main.IsSpecies("O2")) {
    return AdiabaticFactor({0.545e-10}, AdiabaticFactor::Type::Hartmann);
  }
      
  throw std::runtime_error("Error in adiabatic_factor: Unsupported species pairs, main-collider.  Check your broadening data");
}


//! BasisRate function format for supported scenarios
/*!
 \param main: species of interest
 \param collider: species broadening main
 \param T: Atmospheric temperature
 \param T0: Reference temperature
 \return An BasisRate function
 */
inline BasisRate basis_rate(const SpeciesTag& main, const SpeciesTag& collider, const Numeric& T, const Numeric& T0)
{
  if(main.IsSpecies("CO2") and collider.IsSpecies("N2"))
    return BasisRate({
      Conversion:: hitran2arts_broadening(0.0180) * std::pow(T0/T, 0.85),  // Hz/Pa
      0.81 * std::pow(T0/T, 0.0152),  // unitless
      0.008},  // unitless
      BasisRate::Type::Hartmann);
  else if(main.IsSpecies("CO2") and collider.IsSpecies("O2"))
    return BasisRate({
      Conversion:: hitran2arts_broadening(0.0168) * std::pow(T0/T, 0.50),  // Hz/Pa
      0.82 * std::pow(T0/T, -0.091),  // unitless 
      0.007},  // unitless
      BasisRate::Type::Hartmann);
  else if(main.IsSpecies("O2"))
    return BasisRate({-1e99, -1e99, -1e99}, BasisRate::Type::Hartmann);
  
  throw std::runtime_error("Error in basis_rate: Unsupported species pairs, main-collider.  Check your broadening data");
}


// A helper for insider relaxation_matrix_calculations
enum class Species {CO2, O2_66};


/*! Computes an individual W
 \param lines: all absorption lines
 \param population: population distribution of the absorption lines
 \param main: species of interest
 \param collider: species broadening main
 \param collider_vmr: vmr of species broadening main
 \param T: Atmospheric temperature
 \return An individual W
*/
Matrix relaxation_matrix_calculations(const ArrayOfLineRecord& lines,
                                      const Vector& population,
                                      const SpeciesTag& main,
                                      const SpeciesTag& collider,
                                      const Numeric& collider_vmr,
                                      const Numeric& T,
                                      const Index& size)
try
{
  const Index n = lines.nelem();
  Matrix W(n, n);
  
  const BasisRate br = basis_rate(main, collider, T, lines[0].Ti0());
  const AdiabaticFactor af = adiabatic_factor(main, collider);
  const Numeric B0   = getB0(main);
  const ArrayOfArrayOfSpeciesTag pseudo_species({ArrayOfSpeciesTag(1, collider), ArrayOfSpeciesTag(1, main)});
  const Vector pseudo_vmrs({1, 0});
  
  Species spec;
  if(main.IsSpecies("CO2"))
    spec = Species::CO2;
  else if(main.IsSpecies("O2") and main.IsIsotopologue("66"))
    spec = Species::O2_66;
  else throw "Unsupported species";
  
  #pragma omp parallel for schedule(guided, 1) if(DO_FAST_WIGNER && !arts_omp_in_parallel())
  for(Index i=0; i<n; i++) {
    // Create a temporary table to allow openmp
    wig_temp_init(2*int(size));
    
    const auto shape_parameters = lines[i].GetShapeParams(T, 1, 1, pseudo_vmrs, pseudo_species);
    W(i, i) = shape_parameters.G0;
    
    const Numeric& popi = population[i];
    for(Index j=0; j<n; j++) {
      const Numeric& popj = population[j];
      
      if(i not_eq j) {
      OffDiagonalElementOutput X;
        switch(spec) {
          case Species::CO2:
            X = OffDiagonalElement::CO2_IR(lines[i], lines[j], popi, popj, br, af, T, B0, main.SpeciesMass(), collider.SpeciesMass());
            break;
          case Species::O2_66:
            X = OffDiagonalElement::O2_66_MW(lines[i], lines[j], popi, popj, T, collider.SpeciesMass());
            break;
          default: throw "DEVELOPER BUG:  Add species here as well, the other one is to check if the computations are valid at all...";
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
}
catch(const char * e)
{
  std::ostringstream os;
  os << "Errors raised by *relaxation_matrix_calculations*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *relaxation_matrix_calculations*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


/*! Renormalize W
 \param W:  A full relaxation matrix for all absorption lines to be modified
 \param population: population distribution of the absorption lines
 \param d0: reduced dipole moment of the absorption lines
 */
void normalize_relaxation_matrix(Matrix& W,
                                 const Vector& population,
                                 const Vector& d0,
                                 const ArrayOfLineRecord& lines,
                                 const SpeciesAuxData& partition_functions,
                                 const Numeric& T)
try
{
  const Index n=lines.nelem();
  
  Vector d = reduced_dipole_vector(lines, RedPoleType::ElectricRoVibDipole);
  
  // Index list showing sorting of W
  ArrayOfIndex sorted(n, 0);
  Vector test(n);
  if(n) {
    const Numeric QT = single_partition_function(T, partition_functions.getParamType(lines[0]), partition_functions.getParam(lines[0]));
    for(Index i=0; i<n; i++){
      const Numeric QT0 = single_partition_function(lines[i].Ti0(), partition_functions.getParamType(lines[i]), partition_functions.getParam(lines[i]));
      test[i] = Linefunctions::lte_linestrength(lines[i].I0(), lines[i].Elow(), lines[i].F(), QT0, lines[i].Ti0(), QT, T);
    }
  }
  
  get_sorted_indexes(sorted, test);
  for(Index i=0; i<n-i-1; i++) std::swap(sorted[i], sorted[n-i-1]);
  
  // Sorted matrix
  Matrix Wr(n, n);
  for(Index i=0; i<n; i++) {
    Wr(i, i) = W(sorted[i], sorted[i]);
    for(Index j=0; j<n; j++) {
      if(i not_eq j) {
        Wr(i, j) = -std::abs(W(sorted[i], sorted[j]));
      }
    }
  }
  
  // Renormalization procedure
  for(Index i=0; i<n; i++) {
    
    // Sum up upper and lower contributions
    Numeric Sup=0, Slo=0;
    for(Index j=0; j<n; j++) {
      if(j <= i)
        Sup += std::abs(d[sorted[j]]) * Wr(i, j);
      else
        Slo += std::abs(d[sorted[j]]) * Wr(i, j);
    }
    
    Numeric UL = Sup / Slo;
    if(not std::isnormal(UL)) UL = 1.0;
    
    // Rescale to fulfill sum-rule, note how the loop for the upper triangle so the last row cannot be renormalized properly
    for(Index j=i; j<n; j++) {
      const Numeric r = population[sorted[i]] / population[sorted[j]];
      Wr(j, i) = r * Wr(i, j);
      if(j not_eq i)
        Wr(i, j) *= -UL;
    }
  }
  
  for(Index i=0; i<n-1; i++)
    Wr(n-1, i) = 0;  // TEST!
  
  // Backsort matrix before returning
  for(Index i=0; i<n; i++)
    for(Index j=0; j<n; j++)
      W(sorted[i], sorted[j]) = Wr(i, j);
}
catch(const char * e)
{
  std::ostringstream os;
  os << "Errors raised by *normalize_relaxation_matrix*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *normalize_relaxation_matrix*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


/*! Computes an individual W
 \param abs_lines: all absorption lines
 \param main_species: species of interest
 \param population: population distribution of the absorption lines
 \param d0: reduced dipole moment of the absorption lines
 \param colliders: species broadening main
 \param colliders_vmr: vmr of species broadening main
 \param T: Atmospheric temperature
 \return A normalized relaxation matrix
*/
inline Matrix hartmann_relaxation_matrix(const ArrayOfLineRecord& abs_lines,
                                         const SpeciesTag& main_species,
                                         const Vector& population,
                                         const Vector& d0,
                                         const ArrayOfSpeciesTag& colliders,
                                         const Vector& colliders_vmr,
                                         const SpeciesAuxData& partition_functions,
                                         const Numeric& T,
                                         const Index& size)
try
{
  // Size of problem
  const Index n=abs_lines.nelem();
  const Index c=colliders.nelem();
  
  // Create and normalize the matrix
  Matrix W(n, n, 0);
  if(n) {
    for(Index ic=0; ic<c; ic++)
      W += relaxation_matrix_calculations(abs_lines, population, main_species, colliders[ic], colliders_vmr[ic], T, size);
    normalize_relaxation_matrix(W, population, d0, abs_lines, partition_functions, T);
  }
  
  return W;
}
catch(const char * e)
{
  std::ostringstream os;
  os << "Errors raised by *hartmann_relaxation_matrix*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *hartmann_relaxation_matrix*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


inline Numeric population_density(Numeric T, Numeric E0, Numeric F0, Numeric QT) {
  return (1.0 - stimulated_emission(T, F0)) * boltzman_factor(T, E0) / QT;
}


inline Numeric dpopulation_densitydT(Numeric T, Numeric E0, Numeric F0, Numeric QT, Numeric dQTdT) {
  using namespace Constant;
  return (1.0 - dstimulated_emissiondT(T, F0)) *  boltzman_factor(  T, E0) / QT +
         (1.0 -  stimulated_emission(  T, F0)) * dboltzman_factordT(T, E0) / QT +
         (1.0 -  stimulated_emission(  T, F0)) *  boltzman_factor(  T, E0) * (-dQTdT) / pow2(QT);
}


/*! Computes the population distribution
 \param abs_lines: all absorption lines
 \param partition_functions: the partition functions
 \param T: Atmospheric temperature
 \return population distribution of the absorption lines
*/
Vector population_density_vector(const ArrayOfLineRecord& abs_lines,
                                 const SpeciesAuxData& partition_functions,
                                 const Numeric& T)
try
{
  auto n=abs_lines.nelem();
  Vector p(n);
  
  if(n) {
    const Numeric QT = single_partition_function(T, partition_functions.getParamType(abs_lines[0]), partition_functions.getParam(abs_lines[0]));
    
    for(auto i=0; i<n; i++) {
      if(abs_lines[i].IsNotSameSpecIso(abs_lines[0]))
        throw "Not the same species-isotopologue in a set of lines means they are not part of the same band.";
      p[i] = population_density(T, abs_lines[i].Elow(), abs_lines[i].F(), QT);
    }
  }
  
  return p;
}
catch(const char * e)
{
  std::ostringstream os;
  os << "Errors raised by *population_density_vector*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *population_density_vector*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


/*! Computes the dipole moment
 \param abs_lines: all absorption lines of the band
 \return dipole moment of the absorption lines
*/
Vector dipole_vector(const ArrayOfLineRecord& abs_lines,
                     const SpeciesAuxData& partition_functions)
try
{
  const Index n=abs_lines.nelem();
  Vector d(n, 0);
  if(not n) return d;
  
  const Numeric T0 = abs_lines[0].Ti0();
  const Vector p0 = population_density_vector(abs_lines, partition_functions, T0);
  for(Index i=0; i<n; i++) {
    if(T0 not_eq abs_lines[i].Ti0()) throw "Bad T0";
      d[i] = std::sqrt(abs_lines[i].I0()/p0[i]);
  }
  return d;
}
catch(const char * e)
{
  std::ostringstream os;
  os << "Errors raised by *dipole_vector*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *dipole_vector*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


/*! Computes the dipole moment
 * \param abs_lines: all absorption lines of the band
 * \return dipole moment of the absorption lines
 */
Vector reduced_dipole_vector(const ArrayOfLineRecord& abs_lines,
                             const RedPoleType type)
try
{
  const Index n=abs_lines.nelem();
  Vector d(n, 0);
  if(not n) return d;
  
  switch(type) {
    case RedPoleType::ElectricRoVibDipole:
      for(Index i=0; i<n; i++)
        d[i] = abs_lines[i].reduced_rovibrational_dipole();
      break;
    case RedPoleType::MagneticQuadrapole:
      for(Index i=0; i<n; i++)
        d[i] = abs_lines[i].reduced_magnetic_quadrapole();
      break;
  }
  return d;
}
catch(const char * e)
{
  std::ostringstream os;
  os << "Errors raised by *reduced_dipole_vector*:\n";
  os << "\tError: " << e << '\n';
  throw std::runtime_error(os.str());
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *reduced_dipole_vector*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


Vector rosenkranz_first_order(const ArrayOfLineRecord& abs_lines,
                              const Matrix& W,
                              const Vector& d0)
{
  const Index n=abs_lines.nelem();
  Vector Y(n, 0);
  
  for(Index i=0; i<n; i++) {
    Numeric sum=0;
    for(Index j=0; j<n; j++) {
      const Numeric df = abs_lines[i].F() - abs_lines[j].F();
      if(std::isnormal(df))
        sum += d0[j] / d0[i] * W(i, j) / df;
    }
    Y[i] = -2*sum;  // Sign change because ARTS uses (1-iY)
  }
  
  return Y;
}

Vector rosenkranz_shifting_second_order(const ArrayOfLineRecord& abs_lines,
                                        const Matrix& W)
{
  const Index n=abs_lines.nelem();
  Vector DV(n, 0);
  
  for(Index i=0; i<n; i++) {
    Numeric sum=0;
    for(Index j=0; j<n; j++) {
      const Numeric df = abs_lines[j].F() - abs_lines[i].F();
      if(std::isnormal(df))
        sum += W(i, j) * W(j, i) / df;
    }
    DV[i] = sum;  // Does this require a sign change????
  }
  
  return DV;
}


Vector rosenkranz_scaling_second_order(const ArrayOfLineRecord& abs_lines,
                                       const Matrix& W,
                                       const Vector& d0)
{
  const Index n=abs_lines.nelem();
  Vector G(n, 0);
  
  for(Index i=0; i<n; i++) {
    const Numeric& di = d0[i];
    Numeric sum1=0;
    Numeric sum2=0;
    Numeric sum3=0;
    Numeric sum4=0;
    for(Index j=0; j<n; j++) {
      const Numeric df = abs_lines[j].F() - abs_lines[i].F();
      if(std::isnormal(df)) {
        const Numeric& dj = d0[j];
        const Numeric r = dj / di;
        
        sum1 += W(i, j) * W(j, i) / Constant::pow2(df);
        sum2 += r * W(i, j) / df;
        sum3 += r * W(i, j) * W(i, i) / Constant::pow2(df);
        
        Numeric sum_tmp=0;
        for(Index k=0; k<n; k++) {
          const Numeric dfk = abs_lines[k].F() - abs_lines[i].F();
          if(std::isnormal(dfk))
            sum_tmp += W(k, j) * W(i, k) / (dfk * df);
        }
        sum4 += r * sum_tmp;
      }
    }
    G[i] = sum1 - Constant::pow2(sum2) + 2*sum3 - 2*sum4;
  }
  return G;
}


Matrix hartmann_ecs_interface(const ArrayOfLineRecord& abs_lines,
                              const ArrayOfSpeciesTag& main_species,
                              const ArrayOfSpeciesTag& collider_species,
                              const Vector& collider_species_vmr,
                              const SpeciesAuxData& partition_functions,
                              const Numeric& T,
                              const Index& size)
try
{
  const Vector population = population_density_vector(abs_lines, partition_functions, T);
  const Vector dipole = dipole_vector(abs_lines, partition_functions);
  const Matrix W = hartmann_relaxation_matrix(abs_lines, main_species[0], population, dipole, collider_species, collider_species_vmr, partition_functions, T, size);
  return W;
}
catch(const std::exception& e)
{
  std::ostringstream os;
  os << "Errors in calls by *hartmann_ecs_interface*:\n";
  os << e.what();
  throw std::runtime_error(os.str());
}


OffDiagonalElementOutput OffDiagonalElement::CO2_IR(const LineRecord& j_line,
                                                    const LineRecord& k_line,
                                                    const Numeric& j_rho,
                                                    const Numeric& k_rho,
                                                    const BasisRate& br,
                                                    const AdiabaticFactor& af,
                                                    const Numeric& T,
                                                    const Numeric& B0,
                                                    const Numeric& main_mass,
                                                    const Numeric& collider_mass)
{
  // Point at quantum numbers
  const Rational Jku  = k_line.UpperQuantumNumber(QuantumNumberType::J );
  const Rational Jkl  = k_line.LowerQuantumNumber(QuantumNumberType::J );
  const Rational Jju  = j_line.UpperQuantumNumber(QuantumNumberType::J );
  const Rational Jjl  = j_line.LowerQuantumNumber(QuantumNumberType::J );
  const Rational l2ju = j_line.UpperQuantumNumber(QuantumNumberType::l2);
  const Rational l2ku = k_line.UpperQuantumNumber(QuantumNumberType::l2);
  const Rational l2jl = j_line.LowerQuantumNumber(QuantumNumberType::l2);
  const Rational l2kl = k_line.LowerQuantumNumber(QuantumNumberType::l2);
  
  const bool jbig = Jjl >= Jkl;
  
  // Prepare for the wigner calculations --- NOTE: twice the values of J and l2 required for fast Wigner-solver
  // Also, prepare for initial and final phase to go from j to k if Jj >= Jk and vice verse
  const int Ji   = (2*(jbig ? Jju : Jku)).toInt();
  const int Jf   = (2*(jbig ? Jjl : Jkl)).toInt();
  const int Ji_p = (2*(jbig ? Jku : Jju)).toInt();
  const int Jf_p = (2*(jbig ? Jkl : Jjl)).toInt();
  const int li   = (2*(jbig ? l2ju : l2ku)).toInt();
  const int lf   = (2*(jbig ? l2jl : l2kl)).toInt();

  // Find best start and end-point of the summing loop
//   const int st = std::max(Ji - Ji_p, Jf - Jf_p);
  const int en = std::min(Ji + Ji_p, Jf + Jf_p);
  
  // Adiabatic factor for Ji
  const Numeric AF1 = af.get(Ji/2, B0, T, main_mass, collider_mass);
  
  // Scale constant for final state NOTE: "%4" and lack of 2*J because Fast library already doubles these numbers
  const Numeric K1  = ((li+lf)%4 ? 1.0:-1.0) *
        Numeric(Ji_p+1) * sqrt(Numeric((Jf+1)*(Jf_p+1))) * AF1;
  
  Numeric sum=0;
  for(int L = 4; L <= en; L+=4) { 
//   for(int L = st>4?st:4; L <= en; L+=4) { 
    // Basis rate for L
    const Numeric QL = br.get(L/2, B0, T);
    
    // Adiabatic factor for L
    const Numeric AF2 = af.get(L/2, B0, T, main_mass, collider_mass);
    
    // The wigner-symbol following Niro etal 2004
    const Numeric y = co2_ecs_wigner_symbol(Ji, Jf, Ji_p, Jf_p, L, li, lf);
    
    // Sum to the total
    sum += QL * y / AF2;
  }
  
  sum *= K1;
  
  const Numeric r = k_rho / j_rho;
  return {jbig ? sum : sum * r , jbig ? sum / r : sum};
}


Numeric o2_66_inelastic_cross_section_makarov(int L, Numeric T) {
  using namespace Constant;
  using namespace Molecule::O2_66;
  
  Numeric const1 = 0.086 + 8154e-7 * T;
  static constexpr Numeric const2 = 0.5805;
  
  return (2*L + 1) / std::pow(L*L + L, const1) * 
  std::exp(-const2 * h * hamiltonian_freq(L) / (k * T));
}


Numeric o2_66_adiabatic_factor_makarov(int L, Numeric T, Numeric collider_mass) {
  using namespace Constant;
  using namespace Conversion;
  using namespace Molecule::O2_66;
  
  static constexpr Numeric const1 = 0.545e-10;
  static constexpr Numeric constant = 2000 * R * inv_pi * pow2(inv_ln_2);
  
  // Mean speed of collisions
  const Numeric invmu = (1/mass + 1/collider_mass);
  const Numeric vm2 = constant * T * invmu;
  
  return inv_pow2(1.0 + pow2(freq2angfreq(hamiltonian_freq(L)-hamiltonian_freq(L-2)) * const1) / vm2 / 24.0);\
}


OffDiagonalElementOutput OffDiagonalElement::O2_66_MW(const LineRecord& line1,
                                                      const LineRecord& line2,
                                                      const Numeric& rho1,
                                                      const Numeric& rho2,
                                                      const Numeric& T,
                                                      const Numeric& collider_mass)
{
  const Rational J1u = line1.UpperQuantumNumber(QuantumNumberType::J);
  const Rational N1u = line1.UpperQuantumNumber(QuantumNumberType::N);
  const Rational J1l = line1.LowerQuantumNumber(QuantumNumberType::J);
  const Rational N1l = line1.LowerQuantumNumber(QuantumNumberType::N);
  const Rational J2u = line2.UpperQuantumNumber(QuantumNumberType::J);
  const Rational N2u = line2.UpperQuantumNumber(QuantumNumberType::N);
  const Rational J2l = line2.LowerQuantumNumber(QuantumNumberType::J);
  const Rational N2l = line2.LowerQuantumNumber(QuantumNumberType::N);
  
  // Find which is the 'upper' transition
  const bool onebig = Molecule::O2_66::hamiltonian_freq(J1u.toNumeric(), (J1u-N1u).toInt(), (J1u-N1u).toInt())
                    > Molecule::O2_66::hamiltonian_freq(J2u.toNumeric(), (J2u-N2u).toInt(), (J2u-N2u).toInt());
  
  // Define the transitions as in Makarov etal 2013, double the value for wiglib
  const int Nk  = (2 * (onebig ? N1u : N2u).toInt());
  const int Nkp = (2 * (onebig ? N1l : N2l).toInt());
  const int Jk  = (2 * (onebig ? J1u : J2u).toInt());
  const int Jkp = (2 * (onebig ? J1l : J2l).toInt());
  const int Nl  = (2 * (onebig ? N2u : N1u).toInt());
  const int Nlp = (2 * (onebig ? N2l : N1l).toInt());
  const int Jl  = (2 * (onebig ? J2u : J1u).toInt());
  const int Jlp = (2 * (onebig ? J2l : J1l).toInt());
  
  if(Nl not_eq Nlp or Nk not_eq Nkp)
    throw "ERRROR, bad N-values";
  
  // 'length' of numbers
  const Numeric lNk  = std::sqrt(Nk  + 1.0);
  const Numeric lNl  = std::sqrt(Nl  + 1.0);
  const Numeric lJk  = std::sqrt(Jk  + 1.0);
  const Numeric lJl  = std::sqrt(Jl  + 1.0);
  const Numeric lJkp = std::sqrt(Jkp + 1.0);
  const Numeric lJlp = std::sqrt(Jlp + 1.0);
  
  // Constant independenf of L
  const Numeric const1 = lNk * lNl * std::sqrt(lJk*lJl*lJkp*lJlp) * o2_66_inelastic_cross_section_makarov(Nk/2, T);
  
  Numeric sum=0;
  for(int L=4; L<400; L+=4) {
    const int sgn = ((Jk+Jl+L+2) % 4) ? 1 : -1;
    
    const Numeric const2 = sgn * const1 * o2_66_adiabatic_factor_makarov(L/2, T, collider_mass) / o2_66_inelastic_cross_section_makarov(L/2, T);
    
    sum += o2_ecs_wigner_symbol(Nl, Nk, Jl, Jk, Jlp, Jkp, L) * const2;
  }
  
  return {onebig ? sum : sum * rho2 / rho1 , onebig ? sum * rho1 / rho2 : sum};
}


/*! Computes adiabatic factor by
 * 
 * AF = (1 + 1/24 * (w(J) dc / v(T) )^2)^-2
 * 
 */
Numeric AdiabaticFactor::mol_X(const int L,
                               const Numeric& B0,
                               const Numeric& T,
                               const Numeric& main_mass,
                               const Numeric& collider_mass) const
{ 
  using namespace Constant;
  using namespace Conversion;
  
  if(L < 1) return 0.;
  
  static constexpr Numeric constant = 2000 * R * inv_pi * pow2(inv_ln_2);
  
  // Mean speed of collisions
  const Numeric invmu = (1/main_mass + 1/collider_mass);
  const Numeric vm2 = constant * T * invmu;
  
  // FIXME: With or without freq2angfreq???
  return inv_pow2(1.0 + pow2(freq2angfreq(B0) * (4*L-2) * mdata[Index(HartmannPos::dc)]) / vm2 / 24.0);
}


/*! Computes adiabatic factor by
 * 
 * BR = a1 (L*(L+1))^-a2 * exp(- a3 B0 L*(L+1) hc/kT)
 * 
 */
Numeric BasisRate::mol_X(const int L, const Numeric& B0, const Numeric& T) const
{
  using namespace Constant;
  
  const Numeric& a1 = mdata[Index(HartmannPos::a1)];
  const Numeric& a2 = mdata[Index(HartmannPos::a2)];
  const Numeric& a3 = mdata[Index(HartmannPos::a3)];
  
  const Numeric EL = Numeric(L*L + L);
  return a1 / std::pow( EL, a2 ) *  std::exp(- a3 * h * B0 * EL / (k * T));
}


SecondOrderLineMixingCoeffs compute_2nd_order_lm_coeff(ConstVectorView y, ConstVectorView x, const Numeric exp, const Numeric x0)
{
  const auto n=y.nelem();
  
  if(n not_eq x.nelem()) throw std::runtime_error("Bad input to compute_2nd_order_lm_coeff");
  
  Index best_x=0;
  for(auto i=0; i<n-1; i++) {
    if(x[i] >= x[i+1])
      throw std::runtime_error("Bad structure on x-input; must be constantly increasing");
    if(x[i] > x0)
      best_x++;
  }
  
  Matrix A(n, 2);
  Vector ans(2, y[best_x]*0.5), dans(2, 0), delta(n, 0);
  
  Numeric rel_res=1e99;
  Numeric rel_res_limit=1e-16;
  auto loop_count=0;
  constexpr auto max_loop_count=20;
  
  do {
    for(auto i=0; i<n; i++) {
      const Numeric theta = x0 / x[i];
      const Numeric TP = pow(theta, exp);
      A(i, 0) = TP;
      A(i, 1) = (theta - 1.0) * TP;
      const Numeric f = ans[0] * TP + ans[1] * (theta - 1.0) * TP;
      delta[i] = y[i] - f;
    }
    
    // Least square fit
    const Numeric res = lsf(dans, A, delta);
    
    if(not std::isnormal(res)) break;
    
    rel_res = std::abs(dans[0]/ans[0]) + std::abs(dans[1]/ans[1]);
    ans += dans;
    loop_count+=1;
    
  } while(rel_res > rel_res_limit and loop_count < max_loop_count);
  
  return {ans[0], ans[1]};
}


ComplexVector equivalent_linestrengths(const Vector& population,
                                       const Vector& dipole,
                                       const Eigen::ComplexEigenSolver<Eigen::MatrixXcd>& M)
{
  const auto n=population.nelem();
  
  auto& V = M.eigenvectors();
  auto Vinv = V.inverse().eval();
  
  ComplexVector B(n);
  
  for(auto i=0; i<n; i++) {
    for(auto j1=0; j1<n; j1++) {
      for(auto j2=0; j2<n; j2++) {
        B[i] += population[j1]*dipole[j1]*dipole[j2]*V(i, j2)*Vinv(j1, i);
      }
    }
  }
  return B;
}


Numeric total_linestrengths(const Vector& population,
                            const Vector& dipole)
{
  using namespace Constant;
  const auto n=population.nelem();
  
  Numeric I0=0;
  for(auto i=0; i<n; i++)
    I0 += pi * population[i] * pow2(dipole[i]);
  return I0;
}

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
#include "constants.h"

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
  else if(main.IsSpecies("CH4"))
    return Conversion::kaycm2freq(5.2);
  
  throw std::runtime_error("Unsupported main species/isotopologue.  Check your broadening data");
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
      
  throw std::runtime_error("Unsupported species pairs, main-collider.  Check your broadening data");
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
  
  throw std::runtime_error("Unsupported species pairs, main-collider.  Check your broadening data");
}


// A helper for insider relaxation_matrix_calculations
enum class Species {CO2};


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
                                      const Index& size
                                     )
{
  const Index n = lines.nelem();
  
  assert(n);
  assert(n == population.nelem());
  Matrix W(n, n);
  
  const BasisRate br = basis_rate(main, collider, T, lines[0].Ti0());
  AdiabaticFactor af = adiabatic_factor(main, collider);
  const Numeric B0   = getB0(main);
  
  Species spec;
  if(main.IsSpecies("CO2"))
    spec = Species::CO2;
  else throw std::runtime_error("Unsupported species");
  
  #pragma omp parallel for if(DO_FAST_WIGNER && !arts_omp_in_parallel())
  for(Index i=0; i<n; i++) {
    // Create a temporary table to allow openmp
    wig_temp_init(2*int(size));
    
    const Numeric& popi = population[i];
    for(Index j=0; j<n; j++) {
      const Numeric& popj = population[j];
      
      if(i == j) {
        W(i, j) = lines[i].Agam() * pow(lines[i].Ti0()/T, lines[i].Nair());  // FIXME: ADOPT THIS TO BE FOR COLLIDER SPECIES AND NOT JUST AIR
      }
      else {
      OffDiagonalElementOutput X;
        switch(spec) {
          case Species::CO2:
            X = OffDiagonalElement::CO2_IR(lines[i], lines[j], popi, popj, br, af, T, B0, main.SpeciesMass(), collider.SpeciesMass());
            break;
          default: throw std::runtime_error("DEVELOPER BUG:  Add species here as well, the other one is to check if the computations are valid at all...");
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


/*! Renormalize W
 \param W:  A full relaxation matrix for all absorption lines to be modified
 \param population: population distribution of the absorption lines
 \param d0: reduced dipole moment of the absorption lines
 */
void normalize_relaxation_matrix(Matrix& W,
                                 const Vector& population,
                                 const Vector& d0)
{
  const Index n=d0.nelem();
  
  assert(W.ncols() == n);
  assert(W.nrows() == n);
  assert(population.nelem() == n);
  
  // Index list showing sorting of W
  ArrayOfIndex back_sorted(n), sorted(n);
  get_sorted_indexes(back_sorted, population);
  for(Index i=0; i<n; i++) sorted[i] = back_sorted[n-i-1];
  
  // Sorted matrix
  Matrix Wr(n, n);
  for(Index i=0; i<n; i++) {
    for(Index j=0; j<n; j++) {
      if(j==i)
        Wr(i, j) = W(sorted[i], sorted[j]);
      else
        Wr(i, j) = -std::abs(W(sorted[i], sorted[j]));
    }
  }
  
  // Renormalization procedure
  for(Index i=0; i<n; i++) {
    
    // Sum up upper and lower contributions
    Numeric Sup=0, Slo=0;
    for(Index j=0; j<n; j++) {
      if(j <= i)
        Sup += std::abs(d0[sorted[j]]) * Wr(i, j);
      else
        Slo += std::abs(d0[sorted[j]]) * Wr(i, j);
    }
    
    // The ratio between upper and lower contributions
    const Numeric UL = Sup / Slo;
    
    // Rescale to fulfill sum-rule, note how the loop for the upper triangle so the last row cannot be renormalized properly
    for(Index j=i; j<n; j++) {
      const Numeric r = population[sorted[i]] / population[sorted[j]];
      Wr(j, i) = r * Wr(i, j);
      if(j not_eq i)
        Wr(i, j) *= -UL;
    }
  }
  
  // Backsort matrix before returning
  for(Index i=0; i<n; i++)
    for(Index j=0; j<n; j++)
      W(sorted[i], sorted[j]) = Wr(i, j);
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
                                         const Numeric& T,
                                         const Index& size)
{
  // Size of problem
  const Index n=abs_lines.nelem();
  const Index c=colliders.nelem();
  
  // Create and normalize the matrix
  Matrix W(n, n, 0);
  for(Index ic=0; ic<c; ic++)
    W += relaxation_matrix_calculations(abs_lines, population, main_species, colliders[ic], colliders_vmr[ic], T, size);
  normalize_relaxation_matrix(W, population, d0);
  
  return W;
}


/*! Computes the population distribution
 \param abs_lines: all absorption lines
 \param partition_type: the partition function type
 \param partition_data: the partition function data
 \param T: Atmospheric temperature
 \return population distribution of the absorption lines
*/
inline Vector population_density_vector(const ArrayOfLineRecord& abs_lines,
                                        const SpeciesAuxData::AuxType& partition_type,
                                        const ArrayOfGriddedField1& partition_data,
                                        const Numeric& T)
{
  const Index n=abs_lines.nelem();
  Vector a(n);
  
  const Numeric QT = single_partition_function(T, partition_type, partition_data);
  
  for(Index i=0; i<n; i++) {
    const Numeric bT = boltzman_factor(T, abs_lines[i].Elow());
    a[i] = abs_lines[i].G_lower() * bT / QT;
  }
  
  return a;
}


/*! Computes the dipole moment
 \param abs_lines: all absorption lines
 \param partition_type: the partition function type
 \param partition_data: the partition function data
 \return dipole moment of the absorption lines
*/
inline Vector dipole_vector(const ArrayOfLineRecord& abs_lines)
{
  const Index n=abs_lines.nelem();
  Vector d(n);
  
  for(Index i=0; i<n; i++)
    d[i] = std::sqrt(abs_lines[i].electric_dipole_moment_squared());
  return d;
}


inline Vector rosenkranz_first_order(const ArrayOfLineRecord& abs_lines,
                                     const Matrix& W,
                                     const Vector& d0)
{
  const Index n=abs_lines.nelem();
  Vector Y(n);
  
  for(Index i=0; i<n; i++) {
    Numeric sum=0;
    for(Index j=0; j<n; j++) {
      const Numeric df = abs_lines[i].F() - abs_lines[j].F();
      if(std::isnormal(df))
        sum += d0[j] / d0[i] * W(i, j) / df;
    }
    Y[i] = 2*sum;
  }
  
  return Y;
}

inline Vector rosenkranz_shifting_second_order(const ArrayOfLineRecord& abs_lines,
                                               const Matrix& W)
{
  const Index n=abs_lines.nelem();
  Vector DV(n);
  
  for(Index i=0; i<n; i++) {
    Numeric sum=0;
    for(Index j=0; j<n; j++) {
      const Numeric df = abs_lines[j].F() - abs_lines[i].F();
      if(std::isnormal(df))
        sum += W(i, j) * W(j, i) / df;
    }
    DV[i] = sum;
  }
  
  return DV;
}


inline Vector rosenkranz_scaling_second_order(const ArrayOfLineRecord& abs_lines,
                                              const Matrix& W,
                                              const Vector& d0)
{
  const Index n=abs_lines.nelem();
  Vector G(n);
  
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


inline Vector pressure_broadening_from_diagonal(const Matrix& W)
{
  const Index n=W.ncols();
  Vector g(n);
  for(Index i=0; i<n; i++) g[i] = W(i, i);
  return g;
}

Matrix hartmann_ecs_interface(const ArrayOfLineRecord& abs_lines,
                              const SpeciesTag& main_species,
                              const ArrayOfSpeciesTag& collider_species,
                              const Vector& collider_species_vmr,
                              const SpeciesAuxData::AuxType& partition_type,
                              const ArrayOfGriddedField1& partition_data,
                              const Numeric& T,
                              const Index& size,
                              const RelmatType type)
{
  switch(type) {
    case RelmatType::Population: {
        const Vector X = population_density_vector(abs_lines, partition_type, partition_data, T);
        return Matrix(X);
      }
      break;
    case RelmatType::Dipole: {
        const Vector X = dipole_vector(abs_lines);
        return Matrix(X);
      }
      break;
    default:
      break;
  }
  
  const Vector population = population_density_vector(abs_lines, partition_type, partition_data, T);
  const Vector dipole = dipole_vector(abs_lines);
  const Matrix W = hartmann_relaxation_matrix(abs_lines, main_species, population, dipole, collider_species, collider_species_vmr, T, size);
  
  if(type == RelmatType::FullMatrix) return W;
  
  Matrix X(type == RelmatType::SecondOrderRosenkranz ? 4 : 2, abs_lines.nelem());
  switch(type) {
    case RelmatType::SecondOrderRosenkranz:
      X(2, joker) = rosenkranz_scaling_second_order(abs_lines, W, dipole); 
      X(3, joker) = rosenkranz_shifting_second_order(abs_lines, W); /* fallthrough */
    case RelmatType::FirstOrderRosenkranz:
      X(0, joker) = pressure_broadening_from_diagonal(W);
      X(1, joker) = rosenkranz_first_order(abs_lines, W, dipole);
      break;
    default: throw std::runtime_error("Developer error: Invalid type-statement.\n");
  }
  return X;
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
  const Rational& Jku  = k_line.UpperQuantumNumbers()[QuantumNumberType::J ];
  const Rational& Jkl  = k_line.LowerQuantumNumbers()[QuantumNumberType::J ];
  const Rational& Jju  = j_line.UpperQuantumNumbers()[QuantumNumberType::J ];
  const Rational& Jjl  = j_line.LowerQuantumNumbers()[QuantumNumberType::J ];
  const Rational& l2ju = j_line.UpperQuantumNumbers()[QuantumNumberType::l2];
  const Rational& l2ku = k_line.UpperQuantumNumbers()[QuantumNumberType::l2];
  const Rational& l2jl = j_line.LowerQuantumNumbers()[QuantumNumberType::l2];
  const Rational& l2kl = k_line.LowerQuantumNumbers()[QuantumNumberType::l2];
  
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
  const int st = std::max(Ji - Ji_p, Jf - Jf_p);
  const int en = std::min(Ji + Ji_p, Jf + Jf_p);
  
  // Adiabatic factor for Ji
  const Numeric AF1 = af.get(Ji/2, B0, T, main_mass, collider_mass);
  
  // Scale constant for final state NOTE: "%4" and lack of 2*J because Fast library already doubles these numbers
  const Numeric K1  = ((li+lf)%4 ? 1.0:-1.0) *
        Numeric(Ji_p+1) * sqrt(Numeric((Jf+1)*(Jf_p+1))) * AF1;
  
  Numeric sum=0;
  
  for(int L = st?st:4; L <= en; L+=4) { 
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
  
  constexpr Numeric constant = 2000 * R * inv_pi * pow2(inv_ln_2);
  
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
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOLTZMAN_CONST;
  
  const Numeric& a1 = mdata[Index(HartmannPos::a1)];
  const Numeric& a2 = mdata[Index(HartmannPos::a2)];
  const Numeric& a3 = mdata[Index(HartmannPos::a3)];
  
  const Numeric EL = Numeric(L*L + L);
  
  return a1 * std::pow( EL, -a2 ) *  std::exp(- a3 * PLANCK_CONST * B0 * EL / (BOLTZMAN_CONST * T));
}

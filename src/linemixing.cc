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

//! Rotational constants of supported scenarios
/*!
 \param main: species of interest
 \return The rotational constant in Hertz
 */
inline Numeric getB0(const SpeciesTag& main)
{
  extern const Numeric SPEED_OF_LIGHT;
  const Numeric hi2arts = 1e-2 * SPEED_OF_LIGHT;
  
  if(main.IsSpecies("CO2"))
    return 0.39021 * hi2arts;  // Herzberg 1966
  else if(main.IsSpecies("O2")) {
    if(main.IsIsotopologue("66"))
      return 43100.44276e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
    else if(main.IsIsotopologue("67"))
      return 1.35 * hi2arts;
    else if(main.IsIsotopologue("68"))
      return 40707.38657e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
    else if(main.IsIsotopologue("88"))
      return 38313.72938e6;  // B.J. Drouin et al. Journal of Quantitative Spectroscopy & Radiative Transfer 111 (2010) 1167–1173
  }
  else if(main.IsSpecies("CH4"))
    return 5.2 * hi2arts;
  
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
  extern const Numeric HZ2CM;
  extern const Numeric ATM2PA;
  
  if(main.IsSpecies("CO2") and collider.IsSpecies("N2"))
    return BasisRate({0.0180 * pow(T0/T, 0.85) / (HZ2CM * ATM2PA),
      0.81 * pow(T0/T, 0.0152), 0.008}, BasisRate::Type::Hartmann);
  else if(main.IsSpecies("CO2") and collider.IsSpecies("O2"))
    return BasisRate({0.0168 * pow(T0/T, 0.50) / (HZ2CM * ATM2PA),
      0.82 * pow(T0/T, -0.091 ), 0.007}, BasisRate::Type::Hartmann);
  
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
        W(i, j) = 2* lines[i].Agam() * pow(lines[i].Ti0()/T, lines[i].Nair());  // FIXME: ADOPT THIS TO BE FOR COLLIDER SPECIES AND NOT JUST AIR
      }
      else {
      tuple<Numeric, Numeric> X;
        switch(spec) {
          case Species::CO2:
            X = OffDiagonalElement::CO2_IR(lines[i], lines[j], popi, popj, br, af, T, B0, main.SpeciesMass(), collider.SpeciesMass());
            break;
          default: throw std::runtime_error("DEVELOPER BUG:  Add species here as well, the other one is to check if the computations are valid at all...");
        }
        
        W(i, j) = std::move(std::get<0>(X));
        W(j, i) = std::move(std::get<1>(X));
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
  for(Index i=0; i<n; i++)
    for(Index j=0; j<n; j++)
      Wr(i, j) = (j == i) ?     W(sorted[i], sorted[j]) : 
                           -abs(W(sorted[i], sorted[j]));
  
  // Renormalization procedure
  for(Index i=0; i<n; i++) {
    
    // Sum up upper and lower contributions
    Numeric Sup=0, Slo=0;
    for(Index j=0; j<n; j++) {
      if(j <= i)
        Sup += abs(d0[sorted[j]]) * Wr(i, j);
      else
        Slo += abs(d0[sorted[j]]) * Wr(i, j);
    }
    
    // The ratio between upper and lower contributions
    const Numeric UL = Sup / Slo;
    
    // Rescale to fulfill sum-rule, note how the loop for the upper triangle so the last row cannot be renormalized properly
    const Numeric& rho_i = population[sorted[i]];
    for(Index j=i; j<n; j++) {
      if(j not_eq i) {
        Wr(i, j) *= -UL;
        Wr(j, i) = rho_i / population[sorted[j]] * Wr(i, j);
      }
    }
  }
  
  // Test the sum-rule
  for(Index i=0; i<n; i++) {
    Numeric sum=0;
    for(Index j=0; j<n; j++)
      sum += Wr(i, j) * d0[sorted[j]] / d0[sorted[i]];
  }
  
  // Resort matrix before returning
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
inline Vector dipole_vector(const ArrayOfLineRecord& abs_lines,
                            const SpeciesAuxData::AuxType& partition_type,
                            const ArrayOfGriddedField1& partition_data)
{
  const Index n=abs_lines.nelem();
  Vector a(n);
  
  for(Index i=0; i<n; i++) {
    const Numeric QT0 = single_partition_function(abs_lines[i].Ti0(), partition_type, partition_data);
    const Numeric gT0 = stimulated_emission(abs_lines[i].Ti0(), abs_lines[i].F());
    const Numeric bT0 = boltzman_factor(abs_lines[i].Ti0(), abs_lines[i].Elow());
    a[i] = sqrt(abs_lines[i].I0() * QT0 / (abs_lines[i].F() * bT0 * (1 - gT0)));
  }
  
  return a;
}


/*! Computes the reduced dipole moment
 \param abs_lines: all absorption lines
 \param main: the main species 
 \return reduced dipole moment of the absorption lines
 */
inline Vector reduced_dipole_vector(const ArrayOfLineRecord& abs_lines,
                                    const SpeciesTag& main,
                                    const Index& size)
{
  // Create a temporary table
  wig_temp_init(2*int(size));
  
  const Index n=abs_lines.nelem();
  Vector a(n);
  
  if(main.IsSpecies("CO2")) {
    for(Index i=0; i<n; i++) {
      const QuantumIdentifier& qi = abs_lines[i].QuantumIdentity();
      
      const Rational& J1 = qi.LowerQuantumNumbers()[QuantumNumberType::J];
      const Rational& l1 = qi.LowerQuantumNumbers()[QuantumNumberType::l2];
      const Rational& J2 = qi.UpperQuantumNumbers()[QuantumNumberType::J];
      const Rational& l2 = qi.UpperQuantumNumbers()[QuantumNumberType::l2];
      
      const int sqn = ((J2.toIndex() + l2.toIndex()) % 2 ? -1 : +1);
      const Numeric w3j = wigner3j(J1, 1, J2, l1, l2-l1, -l2);
      a[i] = sqn * sqrt((2*J2).toInt() + 1) * w3j;
    }
  }
  else {
    throw std::runtime_error("Cannot support species at this point");
  }
  
  // Remove the temporary table
  wig_temp_free();
  return a;
}


inline Vector rosenkranz_first_order(const ArrayOfLineRecord& abs_lines,
                                     const Matrix& W,
                                     const Vector& d0)
{
  const Index n=abs_lines.nelem();
  Vector Y(n);
  
  for(Index i=0; i<n; i++) {
    const Numeric& di = d0[i];
    Numeric sum=0;
    for(Index j=0; j<n; j++) {
      if(i not_eq j) {
        const Numeric& dj = d0[j];
        const Numeric r = dj / di;
        const Numeric df = abs_lines[i].F() - abs_lines[j].F();
        sum += r * W(i, j) / df;
      }
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
      if(i not_eq j) {
        const Numeric df = abs_lines[j].F() - abs_lines[i].F();
        sum += W(i, j) * W(j, i) / df;
      }
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
    Numeric sum_tmp=0;
    for(Index j=0; j<n; j++) {
      if(i not_eq j) {
        const Numeric& dj = d0[j];
        const Numeric r = dj / di;
        const Numeric df = abs_lines[j].F() - abs_lines[i].F();
        
        sum1 += W(i, j) * W(j, i) / pow(df, 2);
        sum2 += r * W(i, j) / df;
        sum3 += r * W(i, j) * W(i, i) / pow(df, 2);
        for(Index k=0; k<n; k++) {
          if(k not_eq i) {
            const Numeric dfk = abs_lines[k].F() - abs_lines[i].F();
            sum_tmp += W(k, j) * W(i, k) / (dfk * df);
          }
        }
        sum4 += r * sum_tmp;
      }
    }
    G[i] = sum1 - pow(sum2, 2) + 2*sum3 - 2*sum4;
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
                              const Index type)
{
  const Vector population = population_density_vector(abs_lines, partition_type, partition_data, T);
  if(type == - 1) return Matrix(population);
  
  const Vector dipole = dipole_vector(abs_lines, partition_type, partition_data);
  if(type == - 2) return Matrix(dipole);
  
  const Vector d0 = reduced_dipole_vector(abs_lines, main_species, size);
  if(type == - 3) return Matrix(d0);
  
  if(type < 0) throw std::runtime_error("Developer error: Invalid type-statement.\n");
  
  const Matrix W = hartmann_relaxation_matrix(abs_lines, main_species, population, d0, collider_species, collider_species_vmr, T, size);
  if(not type) return W;
  
  Matrix X(type*2, abs_lines.nelem());
  switch(type) {
    case 2:
      X(2, joker) = rosenkranz_scaling_second_order(abs_lines, W, d0); 
      X(3, joker) = rosenkranz_shifting_second_order(abs_lines, W); /* fallthrough */
    case 1:
      X(0, joker) = pressure_broadening_from_diagonal(W);
      X(1, joker) = rosenkranz_first_order(abs_lines, W, d0);
      break;
    default: throw std::runtime_error("Developer error: Invalid type-statement.\n");
  }
  return X;
}


tuple<Numeric, Numeric> OffDiagonalElement::CO2_IR(const LineRecord& j_line,
                                                   const LineRecord& k_line,
                                                   const Numeric& j_rho,
                                                   const Numeric& k_rho,
                                                   const BasisRate& br,
                                                   AdiabaticFactor& af, // non-const to save parts of calcs
                                                   const Numeric& T,
                                                   const Numeric& B0,
                                                   const Numeric& main_mass,
                                                   const Numeric& collider_mass)
{
  const bool jbig = j_line.LowerQuantumNumbers()[QuantumNumberType::J] >= k_line.LowerQuantumNumbers()[QuantumNumberType::J];
  
  // Prepare for the wigner calculations --- NOTE: twice the values of J and l required for Wigner-solver
  const int Ji   = jbig ? (2*j_line.UpperQuantumNumbers()[  QuantumNumberType::J]).toInt():
                          (2*k_line.UpperQuantumNumbers()[  QuantumNumberType::J]).toInt();
  const int Jf   = jbig ? (2*j_line.LowerQuantumNumbers()[  QuantumNumberType::J]).toInt():
                          (2*k_line.LowerQuantumNumbers()[  QuantumNumberType::J]).toInt();
  const int Ji_p = jbig ? (2*k_line.UpperQuantumNumbers()[  QuantumNumberType::J]).toInt():
                          (2*j_line.UpperQuantumNumbers()[  QuantumNumberType::J]).toInt();
  const int Jf_p = jbig ? (2*k_line.LowerQuantumNumbers()[  QuantumNumberType::J]).toInt():
                          (2*j_line.LowerQuantumNumbers()[  QuantumNumberType::J]).toInt();
  const int li   = jbig ? (2*j_line.UpperQuantumNumbers()[  QuantumNumberType::l2]).toInt():
                          (2*k_line.UpperQuantumNumbers()[  QuantumNumberType::l2]).toInt();
  const int lf   = jbig ? (2*j_line.LowerQuantumNumbers()[  QuantumNumberType::l2]).toInt():
                          (2*k_line.LowerQuantumNumbers()[  QuantumNumberType::l2]).toInt();

  // Find best start and end-point of the summing loop
  const int st = std::max(Ji - Ji_p, Jf - Jf_p);
  const int en = std::min(Ji + Ji_p, Jf + Jf_p);
  
  // Adiabatic factor for Ji
  const Numeric AF1 = af.get(Ji/2, 2, B0, T, main_mass, collider_mass, true);
  
  // Scale constant for final state NOTE: "%4" and lack of 2*J because wigner double these numbers
  const Numeric K1  = ((li+lf)%4 ? 1.0:-1.0) *
        Numeric(Ji_p+1) * sqrt(Numeric((Jf+1)*(Jf_p+1))) * AF1;
  
  Numeric sum=0;
  
  for(int L = st?st:4; L <= en; L+=4)
  { 
    // Basis rate for L
    const Numeric QL = br.get(L/2, B0, T);
    
    // Adiabatic factor for L
    const Numeric AF2 = af.get(L/2, 2, B0, T, main_mass, collider_mass, false);
    
    // The wigner-symbol following Niro etal 2004
    const Numeric y = co2_ecs_wigner_symbol(Ji, Jf, Ji_p, Jf_p, L, li, lf);
    
    // Sum to the total
    sum += QL * y / AF2;
  }
  
  sum *= K1;
  return tuple<Numeric, Numeric>(jbig ? sum : sum * k_rho / j_rho , jbig ? sum * j_rho / k_rho : sum);
}


/*! Computes adiabatic factor by
 * 
 * AF = (1 + 1/24 * (w(J) dc / v(T) )^2)^-2
 * 
 */
Numeric AdiabaticFactor::mol_X(const int L,
                               const int s,
                               const Numeric& B0,
                               const Numeric& T,
                               const Numeric& main_mass,
                               const Numeric& collider_mass,
                               const bool init)
{ 
  if(L < 1) return 0.;
  
  if(init) {
    extern const Numeric AVOGADROS_NUMB;
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric PI;
    
    // Mean speed of collisions
    const Numeric invmu = AVOGADROS_NUMB * (1/main_mass + 1/collider_mass);
    const Numeric vm = 2 * sqrt(2.0 * BOLTZMAN_CONST * T * invmu / PI);
    
    // Save the values that are constant
    msave = pow(B0 * mdata[Index(HartmannPos::dc)] / vm, 2) / 24.0;
  }
  
  // Energy of rotational states of the main molecule squared
  const Index wj2 = (s*(-2*L + s - 1)) * (s*(-2*L + s - 1));
  
  return pow(1.0 + msave * Numeric(wj2), -2);
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
  
  return a1 * pow( EL, -a2 ) *  exp(- a3 * PLANCK_CONST * B0 * EL / (BOLTZMAN_CONST * T));
}

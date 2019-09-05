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
 * @file linemixing.h
 * @author Richard Larsson
 * @date 2018-08-06
 * 
 * @brief Line mixing calculation implementation
 * 
 * This file contains only experimental code to test the F90-routines.
 * It might evolve to replace it at some point but that is beyond present ability
 */

#ifndef linemixing_h
#define linemixing_h

#include "absorption.h"
#include "complex.h"
#include "constants.h"
#include "linerecord.h"
#include "rational.h"

/** Compute the rotational energy of a Hund b case molecule
 * 
 * @param[in] J Rotational angular momentum w/ spin
 * @param[in] N Rotational angular momentum w/0 spin
 * @param[in] J2 Rotational angular momentum w/ spin for other level
 * @param[in] B Main rotation constant
 * @param[in] D Main rotation constant squared
 * @param[in] H Main rotation constant cubic
 * @param[in] gamma Second rotation constant offset
 * @param[in] gamma_D Second rotation constant
 * @param[in] gamma_H Second rotation constant squared
 * @param[in] lambda Third energy constant offset
 * @param[in] lambda_D Third energy constant 
 * @param[in] lambda_H Third energy constant squared
 * @return Numeric Energy of level
 */
template <class T>
inline Numeric rotational_energy_hund_b_molecule(T J,
                                                 T N,
                                                 T J2,
                                                 Numeric B,
                                                 Numeric D,
                                                 Numeric H,
                                                 Numeric gamma,
                                                 Numeric gamma_D,
                                                 Numeric gamma_H,
                                                 Numeric lambda,
                                                 Numeric lambda_D,
                                                 Numeric lambda_H) {
  if (N == J and N == J2) {
    const Numeric JJ1 = J * (J + 1);
    const Numeric JJ2 = JJ1 * JJ1;
    const Numeric JJ3 = JJ1 * JJ1 * JJ1;
    return +(B * JJ1 - D * JJ2 + H * JJ3) -
           (gamma + gamma_D * JJ1 + gamma_H * JJ2) +
           2 / 3 * (lambda + lambda_D * JJ1 + lambda_H * JJ2);
  } else if (N == (J - 1) and N == (J2 - 1)) {
    const Numeric JJ1 = J * (J - 1);
    const Numeric JJ2 = JJ1 * JJ1;
    const Numeric JJ3 = JJ1 * JJ1 * JJ1;
    return +(B * JJ1 - D * JJ2 + H * JJ3) +
           (gamma + gamma_D * JJ1 + gamma_H * JJ2) * (J - 1) +
           (lambda + lambda_D * JJ1 + lambda_H * JJ2) *
               (2 / 3 - 2 * J / (2 * J + 1));
  } else if (N == (J + 1) and N == (J2 + 1)) {
    const Numeric JJ1 = (J + 2) * (J + 1);
    const Numeric JJ2 = JJ1 * JJ1;
    const Numeric JJ3 = JJ1 * JJ1 * JJ1;
    return +(B * JJ1 - D * JJ2 + H * JJ3) -
           (gamma + gamma_D * JJ1 + gamma_H * JJ2) * (J + 2) +
           (lambda + lambda_D * JJ1 + lambda_H * JJ2) *
               (2 / 3 - 2 * (J + 1) / (2 * J + 1));

  } else if ((N == (J + 1) and N == (J2 - 1)) or
             (N == (J - 1) and N == (J2 + 1))) {
    const Numeric JJ1 = J * J + J + 1;
    const Numeric JJ2 = JJ1 * JJ1;
    return (lambda + lambda_D * JJ1 + lambda_H * JJ2) * 2 *
           std::sqrt(J * J + J) / (2 * J + 1);
  } else {
    return 0;
  }
}

/** Molecular constant */
namespace Molecule {

/**  Name of molecules */
enum class Name { O2_66, CO2_626 };

/** O2-66 constants */
namespace O2_66 {
constexpr Numeric g_S = 2.002084;
constexpr Numeric ge_l = 2.77e-3;
constexpr Numeric g_r = -1.16e-4;

constexpr Numeric B = 43100.44276e6;
constexpr Numeric D = 145.1271e3;
constexpr Numeric H = 49e-3;

constexpr Numeric lambda = 59501.3438e6;
constexpr Numeric lambda_D = 58.3680e3;
constexpr Numeric lambda_H = 290.8e-3;

constexpr Numeric gamma = -252.58634e6;
constexpr Numeric gamma_D = -243.42;
constexpr Numeric gamma_H = -1.46e-3;

constexpr Numeric mass = 31.989830;

/** Hamiltonian frequency
 * 
 * @param[in] J Rotational constant
 * @param[in] dcol +1, 0, -1
 * @param[in] drow +1, 0, -1
 * @return Numeric Frequency
 */
template <class T>
Numeric hamiltonian_freq(T J, int dcol = 0, int drow = 0) {
  /*
      Matrix:
        (-1,-1), (0,-1), (1,-1)
        (-1, 0), (0, 0), (1, 0)
        (-1, 1), (0, 1), (1, 1)
      */
  return rotational_energy_hund_b_molecule(J + dcol,
                                           J,
                                           J + drow,
                                           B,
                                           D,
                                           H,
                                           gamma,
                                           gamma_D,
                                           gamma_H,
                                           lambda,
                                           lambda_D,
                                           lambda_H);
}
};  // namespace O2_66

/** CO-626 constants */
namespace CO2_626 {
constexpr Numeric B = Conversion::kaycm2freq(0.39021);
template <class T>

/** Hamiltonian frequency
 * 
 * @param[in] J Rotational constant
 * @return constexpr Numeric Frequency
 */
constexpr Numeric hamiltonian_freq(T J) {
  return B * J * (J + 1);
}
}  // namespace CO2_626
};  // namespace Molecule

/** Adiabatic factor computations */
class AdiabaticFactor {
 public:

  /** Methods to compute factor */
  enum class Type { Hartmann };

  /** Description of data in Hartmann Type */
  enum class HartmannPos : Index { dc, SIZE };

  /** Construct a new Adiabatic Factor object
   * 
   * Checks that the type and the input agrees
   * 
   * @param[in] v data of how to compute the adiabatic factor
   * @param[in] t type of adiabatic factor
   */
  AdiabaticFactor(const Vector& v, Type t) : mtype(t), mdata(v) {
    bool error = false;
    switch (mtype) {
      case Type::Hartmann:
        if (Index(HartmannPos::SIZE) not_eq mdata.nelem()) error = true;
        break;
    }

    if (error)
      throw std::runtime_error(
          "Bad initializaton of BasisRate, type and size disagree...");
  }

  /** Hartmann AF
   * 
   * @param[in] L Rotational angular momentum
   * @param[in] B0 Rotational constant
   * @param[in] T Temperature
   * @param[in] main_mass Mass of main molecule
   * @param[in] collider_mass Mass of colliding molecule
   * @return Numeric Adiabatic factor
   */
  Numeric mol_X(const Numeric& L,
                const Numeric& B0,
                const Numeric& T,
                const Numeric& main_mass,
                const Numeric& collider_mass) const;

  /** Get AF
   * 
   * @param[in] L Rotational angular momentum
   * @param[in] B0 Rotational constant
   * @param[in] T Temperature
   * @param[in] main_mass Mass of main molecule
   * @param[in] collider_mass Mass of colliding molecule
   * @return Numeric Adiabatic factor
   */
  Numeric get(const Numeric& L,
              const Numeric& B0,
              const Numeric& T,
              const Numeric& main_mass,
              const Numeric& collider_mass) const {
    switch (mtype) {
      case Type::Hartmann:
        return mol_X(L, B0, T, main_mass, collider_mass);
    }
    return 0;
  }

 private:
  Type mtype;
  Vector mdata;
};

/** Basis rate of transitions */
class BasisRate {
 public:

  /** Type of basis rate */
  enum class Type { Hartmann };

  /** Type of basis rate */
  enum class HartmannPos : Index { a1, a2, a3, SIZE };

  /** Construct a new Basis Rate object
   * 
   * Checks that the type and the input agrees
   * 
   * @param[in] v data for the basis rate
   * @param[in] t type describing position of data
   */
  BasisRate(const Vector& v, Type t) : mtype(t), mdata(v) {
    bool error = false;
    switch (mtype) {
      case Type::Hartmann:
        if (Index(HartmannPos::SIZE) not_eq mdata.nelem()) error = true;
        break;
    }

    if (error)
      throw std::runtime_error(
          "Bad initializaton of BasisRate, type and size disagree...");
  }

  /** Computes the basis rate using Hartman method
   * 
   * @param[in] L Rotational angular momentum
   * @param[in] B0 Rotational angular momentum constant
   * @param[in] T Temperautre
   * @return Numeric Basis rate
   */
  Numeric mol_X(const Numeric& L, const Numeric& B0, const Numeric& T) const;

  /** Get the basis rate
   * 
   * @param[in] L Rotational angular momentum
   * @param[in] B0 Rotational angular momentum constant
   * @param[in] T Temperautre
   * @return Numeric Basis rate
   */
  Numeric get(const Numeric& L, const Numeric& B0, const Numeric& T) const {
    switch (mtype) {
      case Type::Hartmann:
        return mol_X(L, B0, T);
    }
    return 0;
  }

 private:
  Type mtype;
  Vector mdata;
};

/** Struct to help keep position of the two matched outputs clear */
struct OffDiagonalElementOutput {
  Numeric ij, ji;
};

/** Methods to compute off diagonal elements */
namespace OffDiagonalElement {

/** Type of off diagonal element computations */
enum class Type { CO2_IR, O2_66_MW };

/** CO2 IR off diagonal element computer
 * 
 * @param[in] j_line Line at pos j
 * @param[in] k_line Line at pos k
 * @param[in] j_rho Population density at pos j
 * @param[in] k_rho Population density at pos k
 * @param[in] br Basis rate computer
 * @param[in] af Adiabatic factor computer
 * @param[in] T Temperature
 * @param[in] B0 Rotational constant
 * @param[in] main_mass Mass of main molecule
 * @param[in] collider_mass Mass of collider molecule
 * @return OffDiagonalElementOutput for j and k
 */
OffDiagonalElementOutput CO2_IR(const LineRecord& j_line,
                                const LineRecord& k_line,
                                const Numeric& j_rho,
                                const Numeric& k_rho,
                                const BasisRate& br,
                                const AdiabaticFactor& af,
                                const Numeric& T,
                                const Numeric& B0,
                                const Numeric& main_mass,
                                const Numeric& collider_mass);

/** O2-66 MW off diagonal element computer 
 * 
 * @param[in] line1 Line at pos 1
 * @param[in] line2 Line at pos 2
 * @param[in] rho1 Population density at pos 1
 * @param[in] rho2 Population density at pos 2
 * @param[in] T Temperature
 * @param[in] collider_mass Mass of collider
 * @return OffDiagonalElementOutput for 1 and 2
 */
OffDiagonalElementOutput O2_66_MW(const LineRecord& line1,
                                  const LineRecord& line2,
                                  const Numeric& rho1,
                                  const Numeric& rho2,
                                  const Numeric& T,
                                  const Numeric& collider_mass);
};  // namespace OffDiagonalElement

/** Energy corrected sudden relaxation matrix using Hartmann's method
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] main_species Species tag of these lines
 * @param[in] collider_species Species tag of collider
 * @param[in] collider_species_vmr VMR of collider
 * @param[in] partition_functions Method to compute the partition function
 * @param[in] T Temperature
 * @param[in] size Number of elements
 * @return Relaxation Matrix
 */
Matrix hartmann_ecs_interface(const ArrayOfLineRecord& abs_lines,
                              const ArrayOfSpeciesTag& main_species,
                              const ArrayOfSpeciesTag& collider_species,
                              const Vector& collider_species_vmr,
                              const SpeciesAuxData& partition_functions,
                              const Numeric& T,
                              const Index& size);

/** Compute the population density
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] partition_functions Method to compute the partition function
 * @param[in] T Temperature
 * @return Vector The population density for each line
 */
Vector population_density_vector(const ArrayOfLineRecord& abs_lines,
                                 const SpeciesAuxData& partition_functions,
                                 const Numeric& T);

/** Dipole vector
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] partition_functions Method to compute the partition function
 * @return Vector Dipole for each line
 */
Vector dipole_vector(const ArrayOfLineRecord& abs_lines,
                     const SpeciesAuxData& partition_functions);

/** Type of reduced dipole */
enum class RedPoleType { ElectricRoVibDipole, MagneticQuadrapole };

/** Reduced dipole vector
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] type Type of reduced dipole
 * @return Vector Reduced dipole for each line
 */
Vector reduced_dipole_vector(const ArrayOfLineRecord& abs_lines,
                             const RedPoleType type);

/** Computes G for Rosenkranz's line mixing coefficients
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] W Relaxation Matrix
 * @param[in] d0 Dipole vector
 * @return Vector G for each line
 */
Vector rosenkranz_scaling_second_order(const ArrayOfLineRecord& abs_lines,
                                       const Matrix& W,
                                       const Vector& d0);

/** Computes DV for Rosenkranz's line mixing coefficients
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] W Relaxation Matrix
 * @return Vector DV for each line
 */
Vector rosenkranz_shifting_second_order(const ArrayOfLineRecord& abs_lines,
                                        const Matrix& W);

/** Computes Y for Rosenkranz's line mixing coefficients
 * 
 * @param[in] abs_lines One band of lines
 * @param[in] W Relaxation Matrix
 * @param[in] d0 Dipole vector
 * @return Vector Y for each line
 */
Vector rosenkranz_first_order(const ArrayOfLineRecord& abs_lines,
                              const Matrix& W,
                              const Vector& d0);

/** To keep track of (y0 + y1 * (T0 / T - 1.)) * pow(T0 / T, n) */
struct SecondOrderLineMixingCoeffs {
  Numeric y0, y1;
};

/** Fit to a second order line mixing formula the input
 * 
 * Finds best fit [c0, c1] of
 * y(x) = (c0 + c1 * (x0 / x - 1.)) * pow(x0 / x, exp)
 * 
 * @param[in] y Y-axis
 * @param[in] x X-axis
 * @param[in] exp Exponent
 * @param[in] x0 Zero-value of x0
 * @return SecondOrderLineMixingCoeffs 
 */
SecondOrderLineMixingCoeffs compute_2nd_order_lm_coeff(ConstVectorView y,
                                                       ConstVectorView x,
                                                       const Numeric exp,
                                                       const Numeric x0);

/** Equivalent line strengths
 * 
 * @param[in] population The population density for each line
 * @param[in] dipole Dipole for each line
 * @param[in] M Solver
 * @return ComplexVector Equivalent line strengths of each line
 */
ComplexVector equivalent_linestrengths(
    const Vector& population,
    const Vector& dipole,
    const Eigen::ComplexEigenSolver<Eigen::MatrixXcd>& M);

/** Sum of line strengths
 * 
 * @param[in] population The population density for each line
 * @param[in] dipole Dipole for each line
 * @return Numeric The sum
 */
Numeric total_linestrengths(const Vector& population, const Vector& dipole);

/** CO2 IR training algorithm for linearization
 * 
 * @param[in] Ji J init for all lines
 * @param[in] Jf J final for all lines
 * @param[in] l2i l2 init for all lines
 * @param[in] l2f l2 final for all lines
 * @param[in] F0 Central frequency
 * @param[in] d Dipole for each line
 * @param[in] rho The population density for each line
 * @param[in] gamma Pressure broadening for each line
 * @param[in] T Temperature
 * @return Matrix Relaxation Matrix
 */
Matrix CO2_ir_training(const ArrayOfRational& Ji,
                       const ArrayOfRational& Jf,
                       const ArrayOfRational& l2i,
                       const ArrayOfRational& l2f,
                       const Vector& F0,
                       const Vector& d,
                       const Vector& rho,
                       const Vector& gamma,
                       const Numeric& T);

/** Compute the relaxation matrix in air mixture
 * 
 * @param[out] relmat Relaxation matrix
 * @param[in] abs_lines Absorption band
 * @param[in] abs_species Atmospheric species
 * @param[in] partition_functions Partition functions
 * @param[in] wigner_initialized Indication of the Wigner state
 * @param[in] temperature Atmospheric temperature
 * @param[in] species This species index
 */
void relmatInAir(Matrix& relmat,
                 const ArrayOfLineRecord& abs_lines,
                 const ArrayOfArrayOfSpeciesTag& abs_species,
                 const SpeciesAuxData& partition_functions,
                 const Index& wigner_initialized,
                 const Numeric& temperature,
                 const Index& species);

/** Xsec from full relaxation matrix
 * 
 * @param[in,out] xsec Cross-section per level
 * @param[in,out] dxsec_dx Cross-section per species derivatives
 * @param[in] lines Absorption band
 * @param[in] derivatives_data The derivatives in dF
 * @param[in] derivatives_data_position The derivatives positions in dF
 * @param[in] Wmat Relaxation matrix
 * @param[in] Wmat_perturbedT Relaxation matrix perturbed by temperature
 * @param[in] f0 Line frequencies
 * @param[in] f_grid As WSV
 * @param[in] d0 Dipole moments
 * @param[in] rhoT Population density at temperature
 * @param[in] rhoT_perturbedT Population density at temperature perturbed by temperature
 * @param[in] psf Line shifting
 * @param[in] psf_perturbedT Line shifting perturbed by temperature
 * @param[in] T Atmospheric temperature
 * @param[in] isotopologue_ratio Ratio of isotopologue in atmosphere
 * @param[in] this_species Index pointing at this species
 * @param[in] this_level Index pointing at this species
 * @param[in] n Number of lines
 */
void calculate_xsec_from_full_relmat(
  ArrayOfMatrix& xsec,
  ArrayOfArrayOfMatrix& dxsec_dx,
  const ArrayOfLineRecord& lines,
  const ArrayOfRetrievalQuantity& derivatives_data,
  const ArrayOfIndex& derivatives_data_position,
  const ConstMatrixView Wmat,
  const ConstMatrixView Wmat_perturbedT,
  const ConstVectorView f0,
  const ConstVectorView f_grid,
  const ConstVectorView d0,
  const ConstVectorView rhoT,
  const ConstVectorView rhoT_perturbedT,
  const ConstVectorView psf,
  const ConstVectorView psf_perturbedT,
  const Numeric& T,
  const Numeric& isotopologue_ratio,
  const Index& this_species,
  const Index& this_level,
  const Index& n);


/** Xsec from full relaxation matrix
 * 
 * @param[in,out] xsec Cross-section per level
 * @param[in,out] dxsec_dx Cross-section per species derivatives
 * @param[in] lines Absorption band
 * @param[in] derivatives_data The derivatives in dF
 * @param[in] derivatives_data_position The derivatives positions in dF
 * @param[in] Line broadening
 * @param[in] Line broadening temperature derivative
 * @param[in] f0 Line frequencies
 * @param[in] f_grid As WSV
 * @param[in] d0 Dipole moments
 * @param[in] rhoT Population density at temperature
 * @param[in] rhoT_perturbedT Population density at temperature temperature derivative
 * @param[in] psf Line shifting
 * @param[in] psf_perturbedT Line shifting temperature derivative
 * @param[in] Y First order line mixing coefficient
 * @param[in] dY_dT First order line mixing coefficient temperature derivative
 * @param[in] G Second order line mixing coefficient
 * @param[in] dG_dT Second order line mixing coefficient temperature derivative
 * @param[in] DV Second order line mixing shifting coefficient
 * @param[in] dDV_dT Second order line mixing shifting coefficient temperature derivative
 * @param[in] T Atmospheric temperature
 * @param[in] isotopologue_mass Mass of isotopologue
 * @param[in] isotopologue_ratio Ratio of isotopologue in atmosphere
 * @param[in] this_species Index pointing at this species
 * @param[in] this_level Index pointing at this species
 * @param[in] n Number of lines
 */
void calculate_xsec_from_relmat_coefficients(
  ArrayOfMatrix& xsec,
  ArrayOfArrayOfMatrix& dxsec_dx,
  const ArrayOfRetrievalQuantity& derivatives_data,
  const ArrayOfIndex& derivatives_data_position,
  const ConstVectorView pressure_broadening,
  const ConstVectorView dpressure_broadening_dT,
  const ConstVectorView f0,
  const ConstVectorView f_grid,
  const ConstVectorView d0,
  const ConstVectorView rhoT,
  const ConstVectorView drhoT_dT,
  const ConstVectorView psf,
  const ConstVectorView dpsf_dT,
  const ConstVectorView Y,
  const ConstVectorView dY_dT,
  const ConstVectorView G,
  const ConstVectorView dG_dT,
  const ConstVectorView DV,
  const ConstVectorView dDV_dT,
  const Numeric& T,
  const Numeric& isotopologue_mass,
  const Numeric& isotopologue_ratio,
  const Index& this_species,
  const Index& this_level,
  const Index& n);

#ifdef ENABLE_RELMAT
extern "C" {
// This is the interfaces between the Fortran code that calculates W and ARTS
extern void arts_relmat_interface__hartmann_and_niro_type(
    long* nlines,
    double* fmin,
    double* fmax,
    long* M,
    long* I,
    double* v,
    double* S,
    double* gamma_air,
    double* E_double_prime,
    double* n_air,
    long* upper,
    long* lower,
    long* g_prime,
    long* g_double_prime,
    double* temperature,
    double* pressure,
    double* partition_function_t,
    double* partition_function_t0,
    double* isotopologue_mass,
    long* number_of_perturbers,
    long* molecule_code_perturber,
    long* iso_code_perturber,
    double* perturber_mass,
    double* vmr,
    //output+input
    long* debug_in__error_out,
    long* ordered,
    double* tolerance_in_rule_nr2,
    bool* use_adiabatic_factor,
    //outputs
    double* W,
    double* dipole,
    double* rhoT,
    double* Y,
    double* G,
    double* DV);

extern void arts_relmat_interface__linear_type(long* nlines,
                                               double* fmin,
                                               double* fmax,
                                               long* M,
                                               long* I,
                                               double* v,
                                               double* S,
                                               double* gamma_air,
                                               double* E_double_prime,
                                               double* n_air,
                                               long* upper,
                                               long* lower,
                                               long* g_prime,
                                               long* g_double_prime,
                                               double* temperature,
                                               double* pressure,
                                               double* partition_function_t,
                                               double* partition_function_t0,
                                               double* isotopologue_mass,
                                               long* number_of_perturbers,
                                               long* molecule_code_perturber,
                                               long* iso_code_perturber,
                                               double* perturber_mass,
                                               double* vmr,
                                               //output+input
                                               long* debug_in__error_out,
                                               long* ordered,
                                               double* tolerance_in_rule_nr2,
                                               bool* use_adiabatic_factor,
                                               //outputs
                                               double* W,
                                               double* dipole,
                                               double* rhoT,
                                               double* Y,
                                               double* G,
                                               double* DV);

extern double* wigner3j_(double*, double*, double*, double*, double*, double*);
extern double* wigner6j_(double*, double*, double*, double*, double*, double*);
}
#endif  //ENABLE_RELMAT

#endif  // linemixing_h

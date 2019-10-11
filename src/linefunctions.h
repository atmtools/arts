/* Copyright (C) 2017
 * Richard Larsson <ric.larsson@gmail.com>
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
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/*!
 * @file   linefunctions.h
 * @author Richard Larsson
 * @date   2017-05-16
 * 
 * @brief  Stuff related to lineshape functions.
 * 
 * This file should contain complete handling of individual lines.
 * The reason is that the old methods are cumbersome to adapt and need redesigning
 */

#ifndef linefunctions_h
#define linefunctions_h

#include "complex.h"
#include "jacobian.h"
#include "linerecord.h"
#include "absorptionlines.h"

/** Line functions related to line shapes and line strength */
namespace Linefunctions {

/** Size required for data buffer */
constexpr Index ExpectedDataSize() { return 2; }

/** Sets the lineshape normalized to unity.
 * 
 * No line mixing or linestrength is computed.
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     line The absortion line
 * @param[in]     vmrs Atmospheric absorption species VMRs
 * @param[in]     temperature Atmospheric temperature
 * @param[in]     pressure Atmospheric pressure
 * @param[in]     zeeman_df Zeeman splitting coefficient
 * @param[in]     magnetic_magnitude Strength of local magnetic field
 * @param[in]     abs_species Atmospheric absorption species
 */
void set_lineshape(Eigen::Ref<Eigen::VectorXcd> F,
                   const Eigen::Ref<const Eigen::VectorXd> f_grid,
                   const LineRecord& line,
                   const ConstVectorView vmrs,
                   const Numeric& temperature,
                   const Numeric& pressure,
                   const Numeric& zeeman_df,
                   const Numeric& magnetic_magnitude,
                   const ArrayOfArrayOfSpeciesTag& abs_species);

/** Sets the Lorentz line shape. Normalization is unity.
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] data Block of allocated memory.  Output nonsensical
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     zeeman_df Zeeman shift parameter for the line
 * @param[in]     magnetic_magnitude Absolute strength of the magnetic field
 * @param[in]     F0_noshift Central frequency without any shifts
 * @param[in]     X Line shape parameters
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dT Temperature derivatives of line shape parameters
 * @param[in]     dVMR VMR derivatives of line shape parameters
 */
void set_lorentz(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& F0_noshift,
    const LineShape::Output& X,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const LineShape::Output& dT = {0, 0, 0, 0, 0, 0, 0, 0, 0},
    const LineShape::Output& dxdVMR = {0, 0, 0, 0, 0, 0, 0, 0, 0});

/** Sets the HTP line shape. Normalization is unity.
 * 
 * Note:  We are not experienced with this line shape and cannot tell what
 * parameters depends on what input.  It is therefore likely that this function
 * will have to be adapted in the future
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     zeeman_df Zeeman shift parameter for the line
 * @param[in]     magnetic_magnitude Absolute strength of the magnetic field
 * @param[in]     F0_noshift Central frequency without any shifts
 * @param[in]     GD_div_F0 Frequency-independent part of the Doppler broadening
 * @param[in]     X Line shape parameters
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dGD_div_F0_dT Temperature derivative of GD_div_F0
 * @param[in]     dT Temperature derivatives of line shape parameters
 * @param[in]     dVMR VMR derivatives of line shape parameters
 */
void set_htp(Eigen::Ref<Eigen::VectorXcd> F,
             Eigen::Ref<Eigen::MatrixXcd> dF,
             const Eigen::Ref<const Eigen::VectorXd> f_grid,
             const Numeric& zeeman_df,
             const Numeric& magnetic_magnitude,
             const Numeric& F0_noshift,
             const Numeric& GD_div_F0,
             const LineShape::Output& X,
             const ArrayOfRetrievalQuantity& derivatives_data =
                 ArrayOfRetrievalQuantity(),
             const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
             const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
             const Numeric& dGD_div_F0_dT = 0.0,
             const LineShape::Output& dT = {0, 0, 0, 0, 0, 0, 0, 0, 0},
             const LineShape::Output& dVMR = {0, 0, 0, 0, 0, 0, 0, 0, 0});

/** Sets the Voigt line shape. Normalization is unity.
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] data Block of allocated memory.  Output nonsensical
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     zeeman_df Zeeman shift parameter for the line
 * @param[in]     magnetic_magnitude Absolute strength of the magnetic field
 * @param[in]     F0_noshift Central frequency without any shifts
 * @param[in]     GD_div_F0 Frequency-independent part of the Doppler broadening
 * @param[in]     X Line shape parameters
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dGD_div_F0_dT Temperature derivative of GD_div_F0
 * @param[in]     dT Temperature derivatives of line shape parameters
 * @param[in]     dVMR VMR derivatives of line shape parameters
 */
void set_voigt(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& F0_noshift,
    const Numeric& GD_div_F0,
    const LineShape::Output& X,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& dGD_div_F0_dT = 0.0,
    const LineShape::Output& dT = {0, 0, 0, 0, 0, 0, 0, 0, 0},
    const LineShape::Output& dVMR = {0, 0, 0, 0, 0, 0, 0, 0, 0});

/** Sets the Doppler line shape. Normalization is unity.
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] data Block of allocated memory.  Output nonsensical
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     zeeman_df Zeeman shift parameter for the line
 * @param[in]     magnetic_magnitude Absolute strength of the magnetic field
 * @param[in]     F0_noshift Central frequency without any shifts
 * @param[in]     GD_div_F0 Frequency-independent part of the Doppler broadening
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dGD_div_F0_dT Temperature derivative of GD_div_F0
 */
void set_doppler(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& F0_noshift,
    const Numeric& GD_div_F0,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& dGD_div_F0_dT = 0.0);

/** Applies line mixing scaling to already set lineshape and line mirror
 * 
 * Equation: 
 *   with_mirroring:
 *     F := (1+G-iY) * F + (1+G+iY) * Fm
 *   else:
 *     F := (1+G-iY) * F
 * 
 * and appropriate derivatives
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     Fm Mirrored lineshape.  Must be right size
 * @param[in]     dFm Mirrored lineshape derivative.  Must be right size
 * @param[in]     X Line shape parameters
 * @param[in]     with_mirroring Mirror lineshape check
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dT Temperature derivatives of line shape parameters
 * @param[in]     dVMR VMR derivatives of line shape parameters
 */
void apply_linemixing_scaling_and_mirroring(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<Eigen::VectorXcd> Fm,
    const Eigen::Ref<Eigen::MatrixXcd> dFm,
    const LineShape::Output& X,
    const bool with_mirroring,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const LineShape::Output& dT = {0, 0, 0, 0, 0, 0, 0, 0, 0},
    const LineShape::Output& dVMR = {0, 0, 0, 0, 0, 0, 0, 0, 0});

/** Applies Rosenkranz quadratic normalization to already set line shape
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     T Atmospheric temperature at level
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 */
void apply_rosenkranz_quadratic_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const Numeric& T,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier());

/** Applies Van Vleck and Huber normalization to already set line shape
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] data Block of allocated memory.  Output nonsensical
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     T Atmospheric temperature at level
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 */
void apply_VVH_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const Numeric& T,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier());

/** Applies Van Vleck and Weiskopf normalization to already set line shape
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 */
void apply_VVW_scaling(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const Numeric& F0,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier());

/** Gets the local thermodynamic equilibrium line strength
 * 
 * @param[in] S0 Reference linestrength
 * @param[in] E0 Reference line lower energy
 * @param[in] F0 Central frequency without any shifts
 * @param[in] QT0 Partition function at reference temperature
 * @param[in] T0 Reference temperature
 * @param[in] QT Partition function at atmospheric temperature
 * @param[in] T Atmospheric temperature at level
 * 
 * @return The local thermodynamic equilibrium line strength
 */
Numeric lte_linestrength(Numeric S0,
                         Numeric E0,
                         Numeric F0,
                         Numeric QT0,
                         Numeric T0,
                         Numeric QT,
                         Numeric T);

/** Applies linestrength to already set line shape by LTE population type
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in]     line The absortion line
 * @param[in]     T The atmospheric temperature
 * @param[in]     isotopic_ratio The ratio of the isotopologue in the atmosphere
 * @param[in]     QT Partition function at atmospheric temperature of level
 * @param[in]     QT0 Partition function at reference temperature
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dQT_dT Temperature derivative of QT
 */
void apply_linestrength_scaling_by_lte(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const LineRecord& line,
    const Numeric& T,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& dQT_dT = 0.0);

/** Applies linestrength to already set line shape by LTE population type
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in]     band The absorption lines
 * @param[in]     T The atmospheric temperature
 * @param[in]     isotopic_ratio The ratio of the isotopologue in the atmosphere
 * @param[in]     QT Partition function at atmospheric temperature of level
 * @param[in]     QT0 Partition function at reference temperature
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dQT_dT Temperature derivative of QT
 * @param[in]     line_ind The current line's ID
 */
void apply_linestrength_scaling_by_lte(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const AbsorptionLines& band,
    const Numeric& T,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& dQT_dT = 0.0,
    size_t line_ind = 0);

/** Applies linestrength to already set line shape by vibrational level temperatures
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in]     line The absortion line
 * @param[in]     T The atmospheric temperature
 * @param[in]     Tu The upper state vibrational temperature; must be T if level is LTE
 * @param[in]     Tl The lower state vibrational temperature; must be T if level is LTE
 * @param[in]     Evu The upper state vibrational energy; if set funny, yields funny results
 * @param[in]     Evl The lower state vibrational energy; if set funny, yields funny results
 * @param[in]     isotopic_ratio The ratio of the isotopologue in the atmosphere
 * @param[in]     QT Partition function at atmospheric temperature of level
 * @param[in]     QT0 Partition function at reference temperature
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dQT_dT Temperature derivative of QT
 */
void apply_linestrength_scaling_by_vibrational_nlte(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const LineRecord& line,
    const Numeric& T,
    const Numeric& Tu,
    const Numeric& Tl,
    const Numeric& Evu,
    const Numeric& Evl,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& dQT_dT = 0.0);

/** Applies linestrength to already set line shape by vibrational level temperatures
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in]     band The absorption lines
 * @param[in]     T The atmospheric temperature
 * @param[in]     Tu The upper state vibrational temperature; must be T if level is LTE
 * @param[in]     Tl The lower state vibrational temperature; must be T if level is LTE
 * @param[in]     Evu The upper state vibrational energy; if set funny, yields funny results
 * @param[in]     Evl The lower state vibrational energy; if set funny, yields funny results
 * @param[in]     isotopic_ratio The ratio of the isotopologue in the atmosphere
 * @param[in]     QT Partition function at atmospheric temperature of level
 * @param[in]     QT0 Partition function at reference temperature
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dQT_dT Temperature derivative of QT
 * @param[in]     line_ind The current line's ID
 */
void apply_linestrength_scaling_by_vibrational_nlte(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const AbsorptionLines& band,
    const Numeric& T,
    const Numeric& Tu,
    const Numeric& Tl,
    const Numeric& Evu,
    const Numeric& Evl,
    const Numeric& isotopic_ratio,
    const Numeric& QT,
    const Numeric& QT0,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& dQT_dT = 0.0,
    size_t line_ind = 0);

/** Applies the line strength to the line shape using the complex line strength of line mixing.
 * Lineshape is already set.
 * 
 * Note:  This is still under construction (until there is a publication, then this comment should be changed)
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     T Atmospheric temperature at level
 * @param[in]     S_LM Complex linestrength
 * @param[in]     isotopic_ratio The ratio of the isotopologue in the atmosphere
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     dS_LM_dT Temperature derivative of S_LM
 */
void apply_linestrength_from_full_linemixing(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Numeric& F0,
    const Numeric& T,
    const Complex& S_LM,
    const Numeric& isotopic_ratio,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Complex& dS_LM_dT = 0.0);

/** Applies the line strength to the line shape using the dipole information
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     T Atmospheric temperature at level
 * @param[in]     d0 Dipole of the absorption line
 * @param[in]     rho Density (of molecules at level
 * @param[in]     isotopic_ratio The ratio of the isotopologue in the atmosphere
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     drho_dT Temperature derivative of rho
 */
void apply_dipole(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const Numeric& F0,
    const Numeric& T,
    const Numeric& d0,
    const Numeric& rho,
    const Numeric& isotopic_ratio,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
    const QuantumIdentifier& quantum_identity = QuantumIdentifier(),
    const Numeric& drho_dT = 0.0);

/** Applies the line-by-line pressure broadening jacobian for the matching lines
 * 
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     line The absorption line
 * @param[in]     T Atmospheric temperature
 * @param[in]     P Atmospheric pressure
 * @param[in]     vmrs VMRs for line shape broadeners
 */
void apply_lineshapemodel_jacobian_scaling(
    Eigen::Ref<Eigen::MatrixXcd> dF,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const QuantumIdentifier& quantum_identity,
    const LineRecord& line,
    const Numeric& T,
    const Numeric& P,
    const Vector& vmrs);

/** Applies the line-by-line pressure broadening jacobian for the matching lines
 * 
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 * @param[in]     band The absorption lines
 * @param[in]     T Atmospheric temperature
 * @param[in]     P Atmospheric pressure
 * @param[in]     vmrs VMRs for line shape broadeners
 * @param[in]     line_ind The current line's ID
 */
void apply_lineshapemodel_jacobian_scaling(
  Eigen::Ref<Eigen::MatrixXcd> dF,
  const ArrayOfRetrievalQuantity& derivatives_data,
  const ArrayOfIndex& derivatives_data_position,
  const QuantumIdentifier& quantum_identity,
  const AbsorptionLines& band,
  const Numeric& T,
  const Numeric& P,
  const Vector& vmrs,
  size_t line_ind);

/** Returns the frequency-independent part of the Doppler broadening
 * 
 * @param[in] T Atmospheric temperature at level
 * @param[in] mass Mass of molecule under consideration
 * 
 * @return Doppler broadening constant
 */
Numeric DopplerConstant(Numeric T, Numeric mass);

/** Returns the temperature derivative of the frequency-independent part of the Doppler broadening
 * 
 * @param[in] T Atmospheric temperature at level
 * @param[in] dc Output of Linefunctions::DopplerConstant(T, mass)
 * 
 * @return Doppler broadening constant temperature derivative
 */
Numeric dDopplerConstant_dT(const Numeric& T, const Numeric& dc);

/** Cross section and derivatives
 * 
 * Combination function using standard setup to compute
 * line strength and lineshape of a single line.
 * Computes in order the lineshape, the linemirroring,
 * the linenormalization, the linemixing, the 
 * linestrength, the non-lte, and the cutoff frequency.
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in,out] data Block of allocated memory.  Output nonsensical
 * @param[out]    start_cutoff Start pos of cutoff frequency
 * @param[out]    end_cutoff End pos of cutoff frequency
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     line The absortion line
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     volume_mixing_ratio_of_lineshape As name suggests
 * @param[in]     nlte_distribution As name suggests
 * @param[in]     pressure As name suggests
 * @param[in]     temperature As name suggests
 * @param[in]     doppler_constant Frequency-independent part of the Doppler broadening
 * @param[in]     isotopologue_ratio The ratio of the isotopologue in the atmosphere at this level
 * @param[in]     zeeman_df Zeeman shift parameter for the line
 * @param[in]     magnetic_magnitude Absolute strength of the magnetic field
 * @param[in]     ddoppler_constant_dT Temperature derivative of doppler_constant
 * @param[in]     partition_function_at_temperature As name suggests
 * @param[in]     dpartition_function_at_temperature_dT Temeperature derivative of partition_function_at_temperature
 * @param[in]     partition_function_at_line_temperature As name suggests
 * @param[in]     cutoff_call Flag to ignore some functions inside if this call is from the cutoff-computations
 */
void set_cross_section_for_single_line(
    Eigen::Ref<Eigen::VectorXcd> F_full,
    Eigen::Ref<Eigen::MatrixXcd> dF_full,
    Eigen::Ref<Eigen::VectorXcd> N_full,
    Eigen::Ref<Eigen::MatrixXcd> dN_full,
    Eigen::Ref<Eigen::MatrixXcd> data_block,
    Index& start_cutoff,
    Index& nelem_cutoff,
    const Eigen::Ref<const Eigen::VectorXd> f_grid,
    const LineRecord& line,
    const ArrayOfRetrievalQuantity& derivatives_data,
    const ArrayOfIndex& derivatives_data_position,
    const Vector& volume_mixing_ratio_of_lineshape,
    const ConstVectorView nlte_distribution,
    const Numeric& pressure,
    const Numeric& temperature,
    const Numeric& doppler_constant,
    const Numeric& isotopologue_ratio,
    const Numeric& zeeman_df,
    const Numeric& magnetic_magnitude,
    const Numeric& ddoppler_constant_dT,
    const Numeric& partition_function_at_temperature,
    const Numeric& dpartition_function_at_temperature_dT,
    const Numeric& partition_function_at_line_temperature,
    const bool cutoff_call = false);

/** Perform cutoff calculations
 * 
 * Combination function using standard setup to compute line strength and lineshape of a single line.
 * Computes in order the lineshape, the linemirroring, the linenormalization, the linemixing, the 
 * linestrength, the non-lte, and the cutoff frequency.
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     line The absortion line
 * @param[in]     volume_mixing_ratio_of_lineshape As name suggests
 * @param[in]     nlte_distribution As name suggests
 * @param[in]     pressure As name suggests
 * @param[in]     temperature As name suggests
 * @param[in]     doppler_constant Frequency-independent part of the Doppler broadening
 * @param[in]     isotopologue_ratio The ratio of the isotopologue in the atmosphere at this level
 * @param[in]     zeeman_df Zeeman shift parameter for the line
 * @param[in]     magnetic_magnitude Absolute strength of the magnetic field
 * @param[in]     ddoppler_constant_dT Temperature derivative of doppler_constant
 * @param[in]     partition_function_at_temperature As name suggests
 * @param[in]     dpartition_function_at_temperature_dT Temeperature derivative of partition_function_at_temperature
 * @param[in]     partition_function_at_line_temperature As name suggests
 */
void apply_cutoff(Eigen::Ref<Eigen::VectorXcd> F,
                  Eigen::Ref<Eigen::MatrixXcd> dF,
                  Eigen::Ref<Eigen::VectorXcd> N,
                  Eigen::Ref<Eigen::MatrixXcd> dN,
                  const ArrayOfRetrievalQuantity& derivatives_data,
                  const ArrayOfIndex& derivatives_data_position,
                  const LineRecord& line,
                  const Vector& volume_mixing_ratio_of_lineshape,
                  const ConstVectorView nlte_distribution,
                  const Numeric& pressure,
                  const Numeric& temperature,
                  const Numeric& doppler_constant,
                  const Numeric& isotopologue_ratio,
                  const Numeric& zeeman_df,
                  const Numeric& magnetic_magnitude,
                  const Numeric& ddoppler_constant_dT,
                  const Numeric& partition_function_at_temperature,
                  const Numeric& dpartition_function_at_temperature_dT,
                  const Numeric& partition_function_at_line_temperature);

/** Sets cutoff frequency indices
 * 
 * @param[out]    start_cutoff Start pos of cutoff frequency
 * @param[out]    end_cutoff End pos of cutoff frequency
 * @param[in]     f_grid Frequency grid of computations
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     cutoff Cutoff frequency
 */
void find_cutoff_ranges(Index& start_cutoff,
                        Index& nelem_cutoff,
                        const Eigen::Ref<const Eigen::VectorXd> f_grid,
                        const Numeric& F0,
                        const Numeric& cutoff);

/** Applies non-lte linestrength to already set line shape
 * 
 * Works on ratio-inputs, meaning that the total distribution does not have to be known
 * 
 * Cannot support partial derivatives at this point due to ARTS not possessing its own
 * NLTE ratio calculation agenda
 * 
 * @param[in,out] F Lineshape.  Must be right size
 * @param[in,out] dF Lineshape derivative.  Must be right size
 * @param[in,out] N Source lineshape
 * @param[in,out] dN Source lineshape derivative
 * @param[in]     r1 Ratio of molecules at energy level 1
 * @param[in]     r2 Ratio of molecules at energy level 2 
 * @param[in]     g1 Statistical weight of energy level 1
 * @param[in]     g2 Statistical weight of energy level 2
 * @param[in]     A21 Einstein coefficient for the transition from energy level 2 to energy level 1
 * @param[in]     F0 Central frequency without any shifts
 * @param[in]     T Atmospheric temperature
 * @param[in]     derivatives_data The derivatives in dF
 * @param[in]     derivatives_data_position The derivatives positions in dF
 * @param[in]     quantum_identity ID of the absorption line
 */
void apply_linestrength_from_nlte_level_distributions(
    Eigen::Ref<Eigen::VectorXcd> F,
    Eigen::Ref<Eigen::MatrixXcd> dF,
    Eigen::Ref<Eigen::VectorXcd> N,
    Eigen::Ref<Eigen::MatrixXcd> dN,
    const Numeric& r1,
    const Numeric& r2,
    const Numeric& g1,
    const Numeric& g2,
    const Numeric& A21,
    const Numeric& F0,
    const Numeric& T,
    const ArrayOfRetrievalQuantity& derivatives_data =
        ArrayOfRetrievalQuantity(),
    const ArrayOfIndex& derivatives_data_position = ArrayOfIndex(),
                                                      const QuantumIdentifier& quantum_identity = QuantumIdentifier());

class InternalData {
public:
  Eigen::VectorXcd F;
  Eigen::VectorXcd N;
  Eigen::MatrixXcd dF;
  Eigen::MatrixXcd dN;
  
  Eigen::Matrix<Complex, 1, 1> Fc;
  Eigen::Matrix<Complex, 1, 1> Nc;
  Eigen::Matrix<Complex, 1, Eigen::Dynamic> dFc;
  Eigen::Matrix<Complex, 1, Eigen::Dynamic> dNc;
  
  Eigen::Matrix<Complex, Eigen::Dynamic, Linefunctions::ExpectedDataSize()>
  data;
  
  InternalData(Index nf, Index nj) noexcept :
  F(nf), N(nf),
  dF(nf, nj), dN(nf, nj),
  dFc(1, nj), dNc(1, nj),
  data(nf, Linefunctions::ExpectedDataSize()) {}
  
  void SetZero() {
    F.setZero();
    N.setZero();
    dF.setZero();
    dN.setZero();
  }
  
  InternalData& operator+=(const InternalData& other) {
    F.noalias() += other.F;
    N.noalias() += other.N;
    dF.noalias() += other.dF;
    dN.noalias() += other.dN;
    return *this;
  }
};  // InternalData

void set_cross_section_of_band(
  InternalData& scratch,
  InternalData& sum,
  const ConstVectorView f_grid,
  const AbsorptionLines& band,
  const ArrayOfRetrievalQuantity& derivatives_data,
  const ArrayOfIndex& derivatives_data_active,
  const Vector& vmrs,
  const ConstVectorView nlte,  // This must be turned into a map of some kind...
  const Numeric& P,
  const Numeric& T,
  const Numeric& isot_ratio,
  const Numeric& H,
  const Numeric& DC,
  const Numeric& dDCdT,
  const Numeric& QT,
  const Numeric& dQTdT,
  const Numeric& QT0,
  const bool no_negatives=false,
  const bool zeeman=false,
  const Zeeman::Polarization zeeman_polarization=Zeeman::Polarization::Pi);
};  // namespace Linefunctions

#endif  //linefunctions_h



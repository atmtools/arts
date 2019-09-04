/* Copyright (C) 2015
 Richard Larsson <ric.larsson@gmail.com>
 
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
 * @file partial_derivatives.h
 * @author Richard Larsson
 * @date 2015-10-07
 * 
 * @brief Computes partial derivatives (old method)
 */

#ifndef propmatpartials_h
#define propmatpartials_h

#include "jacobian.h"

/** Computes all partial derivatives in old method
 * 
 * @param[in,out] partials_attenuation Partial derivatives of the attenualtion
 * @param[in,out] partials_phase Partial derivatives of the phase
 * @param[in,out] partials_src Partial derivatives of the source
 * @param[in] flag_partials List of derivatives
 * @param[in] flag_partials_position list of relevant derivatives
 * @param[in] CF_A Line shape attenuation and strength without line mixing
 * @param[in] CF_B Line shape phase and strength without line mixing
 * @param[in] C Renormalization factor for the line shape
 * @param[in] dFa_dF Derivative of attenuation line shape wrt frequency factor (strength and line mixing scaled)
 * @param[in] dFb_dF Derivative of phase line shape wrt frequency factor (strength and line mixing scaled)
 * @param[in] dFa_dF Derivative of attenuation line shape wrt pressure factor (strength and line mixing scaled)
 * @param[in] dFb_dF Derivative of phase line shape wrt pressure factor (strength and line mixing scaled)
 * @param[in] f_grid Frequency grid
 * @param[in] this_f_grid Frequency grid range
 * @param[in] temperature Temperature
 * @param[in] sigma Doppler broadening
 * @param[in] K2 Stimulated emission
 * @param[in] dK2_dT Stimulated emission temperature derivative
 * @param[in] K3 Absorption NLTE ratio
 * @param[in] dK3_dT Absorption NLTE ratio temperature derivative
 * @param[in] K4 Source NLTE ratio
 * @param[in] line_frequency Line frequency
 * @param[in] line_strength Line reference strength
 * @param[in] line_temperature Line reference temperature
 * @param[in] line_E_low Line reference lower state energy
 * @param[in] line_E_v_low Line reference lower state vibrational energy
 * @param[in] line_E_v_upp Line reference upper state vibrational energy
 * @param[in] line_T_v_low Line reference lower state vibrational temperature
 * @param[in] line_T_v_upp Line reference upper state vibrational temperature
 * @param[in] Y_LM First order line mixing
 * @param[in] dY_LM_dT First order line mixing temperature derivative
 * @param[in] G_LM Second order line mixing
 * @param[in] dG_LM_dT Second order line mixing temperature derivative
 * @param[in] DF_LM Second order shifting line mixing
 * @param[in] dDF_LM_dT Second order shifting line mixing temperature derivative
 * @param[in] qi Line identifier
 * @param[in] ind_ls Line shape number
 * @param[in] ind_lsn Line shape normalization number
 * @param[in] df_0 Speed-independent line center shift by pressure
 * @param[in] ddf_0_dT Speed-independent line center shift by pressure temperature derivative
 * @param[in] gamma Speed-independent line broadening by pressure
 * @param[in] dgamma_dT Speed-independent line broadening by pressure temperature derivative
 * @param[in] dQ_dT Partition function temperature derivative
 * @param[in] DF_Zeeman Zeeman shifting coefficient
 * @param[in] H_mag_Zeeman Strength of the magnetic field
 * @param[in] do_zeeman Do Zeeman calculations
 * @param[in] pressure_level_index Index pointing at this pressure level
 * @param[in] do_partials_phase Use phase calculations
 * @param[in] do_src Use source calculations
 */
void partial_derivatives_lineshape_dependency(
    ArrayOfMatrix& partials_attenuation,
    ArrayOfMatrix& partials_phase,
    ArrayOfMatrix& partials_src,
    const ArrayOfRetrievalQuantity& flag_partials,
    const ArrayOfIndex& flag_partials_position,
    ConstVectorView CF_A,
    ConstVectorView CF_B,
    ConstVectorView C,
    ConstVectorView dFa_dF,
    ConstVectorView dFb_dF,
    ConstVectorView dFa_dP,
    ConstVectorView dFb_dP,
    ConstVectorView f_grid,
    const Range& this_f_grid,
    const Numeric& temperature,
    const Numeric& sigma,
    const Numeric& K2,
    const Numeric& dK2_dT,
    const Numeric& K3,
    const Numeric& dK3_dT,
    const Numeric& K4,
    const Numeric& line_frequency,
    const Numeric& line_strength,
    const Numeric& line_temperature,
    const Numeric& line_E_low,
    const Numeric& line_E_v_low,
    const Numeric& line_E_v_upp,
    const Numeric& line_T_v_low,
    const Numeric& line_T_v_upp,
    const Numeric& Y_LM,
    const Numeric& dY_LM_dT,
    const Numeric& G_LM,
    const Numeric& dG_LM_dT,
    const Numeric& DF_LM,
    const Numeric& dDF_LM_dT,
    const QuantumIdentifier& qi,
    const Index& ind_ls,
    const Index& ind_lsn,
    const Numeric& df_0,
    const Numeric& ddf_0_dT,
    const Numeric& gamma,
    const Numeric& dgamma_dT,
    const Numeric& dQ_dT,
    const Numeric& DF_Zeeman,
    const Numeric& H_mag_Zeeman,
    const bool do_zeeman,
    const Index& pressure_level_index,
    const bool do_partials_phase,
    const bool do_src);

/** Does this line match the identifier
 * 
 * @param[in] from_jac A line identifier
 * @param[in] species A species-mapped index
 * @param[in] isotopologue A isotopologue-mapped index
 * @param[in] lower_qn Lower state quantum numbers
 * @param[in] upper_qn Upper state quantum numbers
 * @return true if lower_qn and upper_qn perfectly match from_jac
 * @return false Otherwise
 */
bool line_match_line(const QuantumIdentifier& from_jac,
                     const Index& species,
                     const Index& isotopologue,
                     const QuantumNumbers& lower_qn,
                     const QuantumNumbers& upper_qn);

/** Match energy level to the identifier
 * 
 * @param[out] lower_energy_level Does this match the lower energy level?
 * @param[out] upper_energy_level Does this match the upper energy level?
 * @param[in] from_jac A line identifier
 * @param[in] species A species-mapped index
 * @param[in] isotopologue A isotopologue-mapped index
 * @param[in] lower_qn Lower state quantum numbers
 * @param[in] upper_qn Upper state quantum numbers
 */
void line_match_level(bool& lower_energy_level,
                      bool& upper_energy_level,
                      const QuantumIdentifier& from_jac,
                      const Index& species,
                      const Index& isotopologue,
                      const QuantumNumbers& lower_qn,
                      const QuantumNumbers& upper_qn);

#endif

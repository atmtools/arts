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


#ifndef propmatpartials_h
#define propmatpartials_h

#include "jacobian.h"


// Generic function that will calculate all partial derivatives of the line-shape that depends on input per line
void partial_derivatives_lineshape_dependency(ArrayOfMatrix& partials_attenuation,
                                              ArrayOfMatrix& partials_phase, 
                                              ArrayOfMatrix&  partials_src,  
                                              const ArrayOfRetrievalQuantity&  flag_partials, 
                                              const ArrayOfIndex&  flag_partials_position, 
                                              ConstVectorView CF_A,
                                              ConstVectorView CF_B,
                                              ConstVectorView C,
                                              ConstVectorView dFa_dF,
                                              ConstVectorView dFb_dF, 
                                              ConstVectorView dFa_dP, 
                                              ConstVectorView dFb_dP, 
                                              ConstVectorView f_grid, 
                                              const Range&    this_f_grid,
                                              const Numeric&  temperature,
                                              const Numeric&  sigma,
                                              const Numeric&  K2,
                                              const Numeric&  dK2_dT,
                                              const Numeric&  K3,
                                              const Numeric&  dK3_dT,
                                              const Numeric&  K4,
                                              // Line parameters
                                              const Numeric&  line_frequency,
                                              const Numeric&  line_strength,
                                              const Numeric&  line_temperature,
                                              const Numeric&  line_E_low,
                                              const Numeric&  line_E_v_low,
                                              const Numeric&  line_E_v_upp,
                                              const Numeric&  line_T_v_low,
                                              const Numeric&  line_T_v_upp,
                                              const Numeric&  Y_LM,
                                              const Numeric&  dY_LM_dT,
                                              const Numeric&  dY_LM0,
                                              const Numeric&  dY_LM1,
                                              const Numeric&  dY_LMexp,
                                              const Numeric&  G_LM,
                                              const Numeric&  dG_LM_dT,
                                              const Numeric&  dG_LM0,
                                              const Numeric&  dG_LM1,
                                              const Numeric&  dG_LMexp,
                                              const Numeric&  DF_LM,
                                              const Numeric&  dDF_LM_dT,
                                              const Numeric&  dDF_LM0,
                                              const Numeric&  dDF_LM1,
                                              const Numeric&  dDF_LMexp,
                                              const QuantumIdentifier&  qi,
                                              // LINE SHAPE
                                              const Index& ind_ls,
                                              const Index& ind_lsn,
                                              const Numeric& df_0,
                                              const Numeric& ddf_0_dT,
                                              const Numeric& gamma,
                                              const Numeric& dgamma_dT,
                                              const Numeric& dgamma_dSelf,
                                              const Numeric& dgamma_dForeign,
                                              const Numeric& dgamma_dWater,
                                              const Numeric& dpsf_dSelf,
                                              const Numeric& dpsf_dForeign,
                                              const Numeric& dpsf_dWater,
                                              const Numeric& dgamma_dSelfExponent,
                                              const Numeric& dgamma_dForeignExponent,
                                              const Numeric& dgamma_dWaterExponent,
                                              const Numeric& dpsf_dSelfExponent,
                                              const Numeric& dpsf_dForeignExponent,
                                              const Numeric& dpsf_dWaterExponent,
                                              // Partition data parameters
                                              const Numeric& dQ_dT,
                                              // Magnetic variables
                                              const Numeric&  DF_Zeeman,
                                              const Numeric&  H_mag_Zeeman,
                                              const bool      do_zeeman,
                                              // Programming variables
                                              const Index&    pressure_level_index,
                                              const bool      do_partials_phase,
                                              const bool      do_src);

void partial_derivatives_of_stokes_along_path(ArrayOfMatrix& dIn_dt,  //Size is [n_retrieval_points][f_grid,stokes_dim]
                                              const ArrayOfRetrievalQuantity& jacobian_quantities,
                                              ConstTensor4View                extinction_matrices,
                                              const ArrayOfTensor4&           partial_propmat_clearsky,
                                              const Numeric&                  path_length,
                                              const Numeric&                  atmospheric_temperature);

bool line_match_line(const QuantumIdentifier& from_jac,
                     const Index& species,
                     const Index& isotopologue,
                     const QuantumNumbers& lower_qn, 
                     const QuantumNumbers& upper_qn);

void line_match_level(bool& lower_energy_level,
                      bool& upper_energy_level,
                      const QuantumIdentifier& from_jac,
                      const Index& species,
                      const Index& isotopologue,
                      const QuantumNumbers& lower_qn, 
                      const QuantumNumbers& upper_qn);

#endif

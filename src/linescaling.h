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

/** Constains various line scaling functions
 * @file linescaling.h
 * @author Richard Larsson
 * @date 2015-06-25
 * 
 * @brief Headers for line scaling functions
 * 
 * Helps to solve line scaling and derivatives for lbl calculations
 */

#ifndef linescaling_h
#define linescaling_h

#include "absorption.h"

/**  Calculates the line strength scaling parameters for cross section calculations.
 * 
 * @param[out] q_t The partition function at atmospheric temperature (LTE)
 * @param[out] q_ref The partition function at line temperature (LTE)
 * @param[out] partition_ratio The partition ratio (LTE)
 * @param[out] K1 The boltzmann ratio to atmospheric temperature (LTE)
 * @param[out] K2 The relative extra absorption due to NLTE effects
 * @param[out] abs_nlte_ratio The relative extra absorption due to NLTE effects
 * @param[out] src_nlte_ratio The relative extra source due to NLTE effects
 * @param[in]  partition_type Switch for partition type of line
 * @param[in]  partition_data Partition data of line
 * @param[in]  atm_t The path point atmospheric temperature
 * @param[in]  line_t The line reference temperature
 * @param[in]  line_f The line central frequency
 * @param[in]  line_elow The line lower energy level
 * @param[in]  do_nlte Bool for "We need to to NLTE calculations"
 * @param[in]  line_evlow The line lower vibrational energy level
 * @param[in]  line_evupp The line upper vibrational energy level
 * @param[in]  line_evlow_index The line lower vibrational energy level index
 * @param[in]  line_evupp_index The line upper vibrational energy level index
 * @param[in]  atm_t_nlte Vector of NLTE temperatures.  The line knows which ones belong to it.
 */
void GetLineScalingData(Numeric& q_t,
                        Numeric& q_ref,
                        Numeric& partition_ratio,
                        Numeric& K1,
                        Numeric& K2,
                        Numeric& abs_nlte_ratio,
                        Numeric& src_nlte_ratio,
                        const SpeciesAuxData::AuxType& partition_type,
                        const ArrayOfGriddedField1& partition_data,
                        const Numeric& atm_t,
                        const Numeric& line_t,
                        const Numeric& line_f,
                        const Numeric& line_elow,
                        const bool& do_nlte,
                        const Numeric& line_evlow,
                        const Numeric& line_evupp,
                        const Index& line_evlow_index,
                        const Index& line_evupp_index,
                        ConstVectorView atm_t_nlte);

/**  Temperature derivatives of GetLineScalingData(...)
 * 
 * @param[out] dq_t_dT The partition function at atmospheric temperature (LTE)
 * @param[out] dK2_dT The relative extra absorption due to NLTE effects
 * @param[out] dpartition_ratio_dT The partition ratio (LTE)
 * @param[out] dabs_nlte_ratio_dT The relative extra absorption due to NLTE effects
 * @param[out] atm_tv_low NLTE temperature output lower level
 * @param[out] atm_tv_upp NLTE temperature output upper level
 * @param[in]  q_t The partition function at atmospheric temperature (LTE)
 * @param[in]  abs_nlte_ratio The relative extra absorption due to NLTE effects
 * @param[in]  partition_type Switch for partition type of line
 * @param[in]  partition_data Partition data of line
 * @param[in]  atm_t The path point atmospheric temperature
 * @param[in]  line_t The line reference temperature
 * @param[in]  dt The temperature perturbance for perturbed derivatives
 * @param[in]  line_f The line central frequency
 * @param[in]  do_nlte Bool for "We need to to NLTE calculations"
 * @param[in]  line_evlow The line lower vibrational energy level
 * @param[in]  line_evupp The line upper vibrational energy level
 * @param[in]  line_evlow_index The line lower vibrational energy level index
 * @param[in]  line_evupp_index The line upper vibrational energy level index
 * @param[in]  atm_t_nlte Vector of NLTE temperatures.  The line knows which ones belong to it.
 */
void GetLineScalingData_dT(Numeric& dq_t_dT,
                           Numeric& dK2_dT,
                           Numeric& dpartition_ratio_dT,
                           Numeric& dabs_nlte_ratio_dT,
                           Numeric& atm_tv_low,
                           Numeric& atm_tv_upp,
                           const Numeric& q_t,
                           const Numeric& abs_nlte_ratio,
                           const SpeciesAuxData::AuxType& partition_type,
                           const ArrayOfGriddedField1& partition_data,
                           const Numeric& atm_t,
                           const Numeric& line_t,
                           const Numeric& dt,
                           const Numeric& line_f,
                           const bool& do_nlte,
                           const Numeric& line_evlow,
                           const Numeric& line_evupp,
                           const Index& line_evlow_index,
                           const Index& line_evupp_index,
                           ConstVectorView atm_t_nlte);

/**  Central frequency derivatives of GetLineScalingData(...)
 * 
 * @param[out] dK2_dF0 The relative extra absorption due to NLTE effects
 * @param[out] dabs_nlte_ratio_dF0 The relative extra absorption due to NLTE effects
 * @param[in]  atm_t The path point atmospheric temperature
 * @param[in]  line_t The line reference temperature
 * @param[in]  atm_tv_low NLTE temperature output lower level
 * @param[in]  atm_tv_upp NLTE temperature output upper level
 * @param[in]  line_evlow The line lower vibrational energy level
 * @param[in]  line_evupp The line upper vibrational energy level
 * @param[in]  line_f The line central frequency
 */
void GetLineScalingData_dF0(Numeric& dK2_dF0,
                            Numeric& dabs_nlte_ratio_dF0,
                            const Numeric& atm_t,
                            const Numeric& line_t,
                            const Numeric& atm_tv_low,
                            const Numeric& atm_tv_upp,
                            const Numeric& line_evlow,
                            const Numeric& line_evupp,
                            const Numeric& line_f);

/**  Calculates the partition function
 * 
 * @param[out] q_ref The partition function at line temperature (LTE)
 * @param[out] q_t The partition function at atmospheric temperature (LTE)
 * @param[in]  line_t The line reference temperature
 * @param[in]  atm_t The path point atmospheric temperature
 * @param[in]  partition_type Switch for partition type of line
 * @param[in]  partition_data Partition data of line
 */
void partition_function(Numeric& q_ref,
                        Numeric& q_t,
                        const Numeric& line_t,
                        const Numeric& atm_t,
                        const SpeciesAuxData::AuxType& partition_type,
                        const ArrayOfGriddedField1& partition_data);

/**  Calculates the partition function temperature derivative
 * 
 * @param[out] dq_t_dT The partition function at atmospheric temperature (LTE)
 * @param[in]  q_t The partition function at atmospheric temperature (LTE)
 * @param[in]  atm_t The path point atmospheric temperature
 * @param[in]  dt The temperature perturbance for perturbed derivatives
 * @param[in]  partition_type Switch for partition type of line
 * @param[in]  partition_data Partition data of line
 */
void dpartition_function_dT(Numeric& dq_t_dT,
                            const Numeric& q_t,
                            const Numeric& atm_t,
                            const Numeric& dt,
                            const SpeciesAuxData::AuxType& partition_type,
                            const ArrayOfGriddedField1& partition_data);

/** Interpolate the partition functions to temperatures
 * 
 * @param[out] q_ref The partition function at line temperature (LTE)
 * @param[out] q_t The partition function at atmospheric temperature (LTE)
 * @param[in]  ref The line reference temperature
 * @param[in]  t The path point atmospheric temperature
 * @param[in]  t_grid The temperature grid
 * @param[in]  q_grid The value grid to interpolate from
 * @param[in]  interp_order Order of interpolation
 */
void CalculatePartitionFctFromData(Numeric& q_ref,
                                   Numeric& q_t,
                                   const Numeric& ref,
                                   const Numeric& t,
                                   ConstVectorView t_grid,
                                   ConstVectorView q_grid,
                                   const Index& interp_order);

/** Temperature derivative of CalculatePartitionFctFromData
 * 
 * @param[out] dQ_dT The partition function at line temperature (LTE)
 * @param[in]  t The path point atmospheric temperature
 * @param[in]  dT The temperature perturbance
 * @param[in]  t_grid The temperature grid
 * @param[in]  q_grid The value grid to interpolate from
 * @param[in]  interp_order Order of interpolation
 */
void CalculatePartitionFctFromData_perturbed(Numeric& dQ_dT,
                                             const Numeric& t,
                                             const Numeric& dT,
                                             const Numeric& q_t,
                                             ConstVectorView t_grid,
                                             ConstVectorView q_grid,
                                             const Index& interp_order);

/** Get the partition functions at the temperatures by polynominal fit
 * 
 * @param[out] q_ref The partition function at line temperature (LTE)
 * @param[out] q_t The partition function at atmospheric temperature (LTE)
 * @param[in]  ref The line reference temperature
 * @param[in]  t The path point atmospheric temperature
 * @param[in]  q_grid The polynominal coefficients
 */
void CalculatePartitionFctFromCoeff(Numeric& q_ref,
                                    Numeric& q_t,
                                    const Numeric& ref,
                                    const Numeric& t,
                                    ConstVectorView q_grid);

/** Temperature derivative of CalculatePartitionFctFromCoeff
 * 
 * @param[out] dQ_dT The partition function at line temperature (LTE)
 * @param[in]  t The path point atmospheric temperature
 * @param[in]  q_grid The polynominal coefficients
 */
void CalculatePartitionFctFromCoeff_dT(Numeric& dQ_dT,
                                       const Numeric& t,
                                       ConstVectorView q_grid);

/** Computes the partition function at one temperature
 * 
 * @param[in] T Temperature
 * @param[in] partition_type Switch for partition type of line
 * @param[in] partition_data Partition data of line
 * 
 * @return partition function
 */
Numeric single_partition_function(const Numeric& T,
                                  const SpeciesAuxData::AuxType& partition_type,
                                  const ArrayOfGriddedField1& partition_data);

/** Computes the partition function temperature derivative
 * 
 * @param[in] QT partition function
 * @param[in] T Temperature
 * @param[in] dT Temperature perturbance
 * @param[in] partition_type Switch for partition type of line
 * @param[in] partition_data Partition data of line
 * 
 * @return partition function derivative wrt temperature
 */
Numeric dsingle_partition_function_dT(
    const Numeric& QT,
    const Numeric& T,
    const Numeric& dT,
    const SpeciesAuxData::AuxType& partition_type,
    const ArrayOfGriddedField1& partition_data);

Numeric stimulated_emission(Numeric T, Numeric F0);

Numeric dstimulated_emissiondT(Numeric, Numeric);

Numeric dstimulated_emissiondF0(Numeric, Numeric);

Numeric stimulated_relative_emission(const Numeric& gamma,
                                     const Numeric& gamma_ref);

Numeric dstimulated_relative_emission_dT(const Numeric& gamma,
                                         const Numeric& gamma_ref,
                                         const Numeric& F0,
                                         const Numeric& T);

Numeric dstimulated_relative_emission_dF0(const Numeric& gamma,
                                          const Numeric& gamma_ref,
                                          const Numeric& T,
                                          const Numeric& T0);

Numeric boltzman_ratio(const Numeric& T, const Numeric& T0, const Numeric& E0);

Numeric dboltzman_ratio_dT(const Numeric& boltzmann_ratio,
                           const Numeric& T,
                           const Numeric& E0);

Numeric boltzman_factor(Numeric, Numeric);

Numeric dboltzman_factordT(Numeric, Numeric);

Numeric dboltzman_factordE0(Numeric, Numeric);

Numeric absorption_nlte_ratio(const Numeric& gamma,
                              const Numeric& r_upp = 1.0,
                              const Numeric& r_low = 1.0);

Numeric dabsorption_nlte_rate_dT(const Numeric& gamma,
                                 const Numeric& T,
                                 const Numeric& F0,
                                 const Numeric& El,
                                 const Numeric& Eu,
                                 const Numeric& r_upp = 1.0,
                                 const Numeric& r_low = 1.0);

Numeric dabsorption_nlte_rate_dF0(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& r_upp = 1.0,
                                  const Numeric& r_low = 1.0);

Numeric dabsorption_nlte_rate_dTl(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tl,
                                  const Numeric& El,
                                  const Numeric& r_low = 1.0);

Numeric dabsorption_nlte_rate_dTu(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tu,
                                  const Numeric& Eu,
                                  const Numeric& r_upp = 1.0);

#endif  // linescaling_h

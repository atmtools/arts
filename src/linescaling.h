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
#include "constants.h"

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
                                  const ArrayOfGriddedField1& partition_data) noexcept;

/** Computes the partition function temperature derivative
 * 
 * @param[in] T Temperature
 * @param[in] partition_type Switch for partition type of line
 * @param[in] partition_data Partition data of line
 * 
 * @return partition function derivative wrt temperature
 */
Numeric dsingle_partition_function_dT(
    const Numeric& T,
    const SpeciesAuxData::AuxType& partition_type,
    const ArrayOfGriddedField1& partition_data) noexcept;

/** Computes exp(-hf/kT)
 * 
 * @param[in] T Temperatures
 * @param[in] F0 Frequency
 * 
 * @return exp(-hf/kT)
 */
Numeric stimulated_emission(Numeric T, Numeric F0);

/** Computes temperature derivative of exp(-hf/kT)
 * 
 * @param[in] T Temperatures
 * @param[in] F0 Frequency
 * 
 * @return dexp(-hf/kT) / dT
 */
Numeric dstimulated_emissiondT(Numeric T, Numeric F0);

/** Computes frequency derivative of exp(-hf/kT)
 * 
 * @param[in] T Temperatures
 * @param[in] F0 Frequency
 * 
 * @return dexp(-hf/kT) / dF0
 */
Numeric dstimulated_emissiondF0(Numeric T, Numeric F0);

/** Computes
 * 
 * \f[ \frac{1 - e^{h f_0/k T}}{1 - e^{h f_0 / k T_0}} \f]
 * 
 * using std::expm1 for increased low number accuracies
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] gamma_ref Stimulated emission at reference temperature
 * 
 * @return \f$ \frac{1 - e^{h f_0/k T}}{1 - e^{h f_0 / k T_0}} \f$
 */
Numeric stimulated_relative_emission(const Numeric F0, const Numeric T0, const Numeric T) noexcept ;

/** Computes
 * 
 * \f[ \frac{1}{T} \frac{h f_0}{k T} \frac{e^{h f_0/k T}}{1 - e^{h f_0 / k T_0}} \f]
 * 
 * using std::expm1 for increased low number accuracies
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] gamma_ref Stimulated emission at reference temperature
 * 
 * @return  \f$ \frac{1}{T} \frac{h f_0}{k T} \frac{e^{h f_0/k T}}{1 - e^{h f_0 / k T_0}} \f$
 */
Numeric dstimulated_relative_emission_dT(const Numeric F0, const Numeric T0, const Numeric T) noexcept;

/** Computes
 * 
 * \f[ \frac{h\left(T\left[1 -e^{h f_0/k T}\right] e^{h f_0 / k T_0} - T_0 \left[1 -e^{h f_0/k T_0}\right] e^{h f_0 / k T} \right)}{kTT_0 \left(1 -e^{h f_0/k T_0}\right)^2} \f]
 * 
 * using std::expm1 for increased low number accuracies
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] gamma_ref Stimulated emission at reference temperature
 * 
 * @return  \f$ \frac{h\left(T\left[1 -e^{h f_0/k T}\right] e^{h f_0 / k T_0} - T_0 \left[1 -e^{h f_0/k T_0}\right] e^{h f_0 / k T} \right)}{kTT_0 \left(1 -e^{h f_0/k T_0}\right)^2} \f$
 */
Numeric dstimulated_relative_emission_dF0(const Numeric F0, const Numeric T0, const Numeric T) noexcept;

/** Computes (1 - gamma) / (1 - gamma_ref)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] gamma_ref Stimulated emission at reference temperature
 * 
 * @return (1 - gamma) / (1 - gamma_ref)
 */
Numeric stimulated_relative_emission(const Numeric& gamma,
                                     const Numeric& gamma_ref);

/** Computes temperature derivative of (1 - gamma) / (1 - gamma_ref)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] gamma_ref Stimulated emission at reference temperature
 * @param[in] F0 Frequency
 * @param[in] T Temperature
 * 
 * @return d[(1 - gamma) / (1 - gamma_ref)] / dT
 */
Numeric dstimulated_relative_emission_dT(const Numeric& gamma,
                                         const Numeric& gamma_ref,
                                         const Numeric& F0,
                                         const Numeric& T);

/** Computes frequency derivative of (1 - gamma) / (1 - gamma_ref)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] gamma_ref Stimulated emission at reference temperature
 * @param[in] T Temperature
 * @param[in] T0 Reference temperature
 * 
 * @return d[(1 - gamma) / (1 - gamma_ref)] / dF0
 */
Numeric dstimulated_relative_emission_dF0(const Numeric& gamma,
                                          const Numeric& gamma_ref,
                                          const Numeric& T,
                                          const Numeric& T0);

/** Computes exp(E0/c  (T - T0) / (T * T0))
 * 
 * @param[in] T Temperature
 * @param[in] T0 Reference temperature
 * @param[in] E0 Lower state energy
 * 
 * @return exp(E0/c  (T - T0) / (T * T0))
 */
Numeric boltzman_ratio(const Numeric& T, const Numeric& T0, const Numeric& E0);

/** Computes temperature derivatives exp(E0/k  (T - T0) / (T * T0))
 * 
 * @param[in] boltzmann_ratio Output of boltzmann_ratio(...)
 * @param[in] T Temperature
 * @param[in] E0 Lower state energy
 * 
 * @return exp(E0/k  (T - T0) / (T * T0))
 */
Numeric dboltzman_ratio_dT(const Numeric& boltzmann_ratio,
                           const Numeric& T,
                           const Numeric& E0);

/** Computes temperature derivatives exp(E0/k  (T - T0) / (T * T0)) / 
 * exp(E0/c  (T - T0) / (T * T0))
 * 
 * @param[in] T Temperature
 * @param[in] E0 Lower state energy
 * 
 * @return E0 k / T^2
 */
constexpr Numeric dboltzman_ratio_dT_div_boltzmann_ratio(Numeric T,
                                                         Numeric E0)
{
  return E0 / (Constant::k * T * T);
}

/** Computes exp(- E0/kT)
 * 
 * @param[in] T Temperature
 * @param[in] E0 Lower state energy
 * 
 * @return exp(- E0/kT)
 */
Numeric boltzman_factor(Numeric T, Numeric E0);

/** Computes temperature derivatives exp(- E0/kT)
 * 
 * @param[in] T Temperature
 * @param[in] E0 Lower state energy
 * 
 * @return dexp(- E0/kT) / dT
 */
Numeric dboltzman_factordT(Numeric T, Numeric E0);

/** Computes lower state energy derivatives exp(- E0/kT)
 * 
 * @param[in] T Temperature
 * @param[in] E0 Lower state energy
 * 
 * @return dexp(- E0/kT) / dE0
 */
Numeric dboltzman_factordE0(Numeric T, Numeric E0);

/** Computes (r_low - r_upp * gamma) / (1 - gamma)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] r_upp Relative ratio NLTE/LTE in upper ratio
 * @param[in] r_low Relative ratio NLTE/LTE in lower ratio
 * 
 * @return (r_low - r_upp * gamma) / (1 - gamma)
 */
Numeric absorption_nlte_ratio(const Numeric& gamma,
                              const Numeric& r_upp = 1.0,
                              const Numeric& r_low = 1.0) noexcept;

/** Computes temperature derivatives of (r_low - r_upp * gamma) / (1 - gamma)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] T Temperature
 * @param[in] F0 Central frequency
 * @param[in] El Lower state energy
 * @param[in] Eu Upper state energy
 * @param[in] r_upp Relative ratio NLTE/LTE in upper ratio
 * @param[in] r_low Relative ratio NLTE/LTE in lower ratio
 * 
 * @return d[(r_low - r_upp * gamma) / (1 - gamma)] / dT
 */
Numeric dabsorption_nlte_rate_dT(const Numeric& gamma,
                                 const Numeric& T,
                                 const Numeric& F0,
                                 const Numeric& El,
                                 const Numeric& Eu,
                                 const Numeric& r_upp = 1.0,
                                 const Numeric& r_low = 1.0);

/** Computes  frequency derivative of (r_low - r_upp * gamma) / (1 - gamma)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] T Temperature
 * @param[in] r_upp Relative ratio NLTE/LTE in upper ratio
 * @param[in] r_low Relative ratio NLTE/LTE in lower ratio
 * 
 * @return d[(r_low - r_upp * gamma) / (1 - gamma)] / dF0
 */
Numeric dabsorption_nlte_rate_dF0(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& r_upp = 1.0,
                                  const Numeric& r_low = 1.0);

/** Computes lower state temperature derivative of (r_low - r_upp * gamma) / (1 - gamma)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] T Temperature
 * @param[in] Tl Temperature of lower level
 * @param[in] El Lower state energy
 * @param[in] r_low Relative ratio NLTE/LTE in lower ratio
 * 
 * @return d[(r_low - r_upp * gamma) / (1 - gamma)] / dTl
 */
Numeric dabsorption_nlte_rate_dTl(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tl,
                                  const Numeric& El,
                                  const Numeric& r_low = 1.0);

/** Computes upper state temperature derivative of (r_low - r_upp * gamma) / (1 - gamma)
 * 
 * @param[in] gamma Stimulated emission at temperature
 * @param[in] T Temperature
 * @param[in] Tu Temperature of lower level
 * @param[in] Eu Lower state energy
 * @param[in] r_upp Relative ratio NLTE/LTE in upper ratio
 * 
 * @return d[(r_low - r_upp * gamma) / (1 - gamma)] / dTu
 */
Numeric dabsorption_nlte_rate_dTu(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tu,
                                  const Numeric& Eu,
                                  const Numeric& r_upp = 1.0);

#endif  // linescaling_h

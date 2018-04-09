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
 * \file   lineshapesdata.h
 * \brief  Stuff related to lineshape functions.
 * 
 * This file should contain complete handling of individual lines.
 * The reason is that the old methods are cumbersome to adapt and need redesigning
 * 
 * \author Richard Larsson
 * \date   2017-05-16
 */

#ifndef linefunctions_h
#define linefunctions_h

#include "complex.h"
#include "partial_derivatives.h"
#include "linerecord.h"


/*
 * Class to solve the problem
 * 
 *      cross-section of a line equals line strength times line shape times rescaling and normalizations
 * 
 * which means
 *      
 *      dsigma = dS x F + S x dF
 * 
 * TODO: Find work-around for incomplete line-shapes like "Voigt Kuntz"
 */


namespace Linefunctions
{
  void set_lineshape(ComplexVectorView F, 
                     const LineRecord& line, 
                     ConstVectorView f_grid, 
                     ConstVectorView vmrs, 
                     const Numeric& temperature, 
                     const Numeric& pressure, 
                     const Numeric& pressure_limit_for_linemixing, 
                     const Numeric& zeeman_df, 
                     const Numeric& magnetic_magnitude,
                     const ArrayOfIndex& broad_spec_locations,
                     const Index& this_species,
                     const Index& water_species,
                     const Verbosity& verbosity);
  
  void set_lorentz(ComplexVectorView F, 
                   ArrayOfComplexVector& dF, 
                   ConstVectorView f_grid, 
                   const Numeric& zeeman_df, 
                   const Numeric& magnetic_magnitude, 
                   const Numeric& F0_noshift, 
                   const Numeric& G0, 
                   const Numeric& L0, 
                   const Numeric& dF0,
                   const PropmatPartialsData& derivatives_data=PropmatPartialsData(), 
                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(), 
                   const Numeric& dG0_dT=0.0, 
                   const Numeric& dL0_dT=0.0,
                   const Numeric& ddF0_dT=0.0,
                   const Range& df_range=Range(joker));
  
  void set_htp(ComplexVectorView F,
               ArrayOfComplexVector& dF,
               ConstVectorView f_grid,
               const Numeric& zeeman_df,
               const Numeric& magnetic_magnitude,
               const Numeric& F0_noshift,
               const Numeric& GD_div_F0,
               const Numeric& G0,
               const Numeric& L0,
               const Numeric& G2,
               const Numeric& L2,
               const Numeric& eta,
               const Numeric& FVC,
               const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
               const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
               const Numeric& dGD_div_F0_dT=0.0,
               const Numeric& dG0_dT=0.0,
               const Numeric& dL0_dT=0.0,
               const Numeric& dG2_dT=0.0,
               const Numeric& dL2_dT=0.0,
               const Numeric& deta_dT=0.0,
               const Numeric& dFVC_dT=0.0,
               const Range& df_range=Range(joker));
  
  void set_faddeeva_algorithm916(ComplexVectorView F,
                                 ArrayOfComplexVector& dF,
                                 ConstVectorView f_grid,
                                 const Numeric& zeeman_df,
                                 const Numeric& magnetic_magnitude,
                                 const Numeric& F0_noshift,
                                 const Numeric& GD_div_F0,
                                 const Numeric& G0,
                                 const Numeric& L0,
                                 const Numeric& dF0,
                                 const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                 const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                 const Numeric& dGD_div_F0_dT=0.0,
                                 const Numeric& dG0_dT=0.0,
                                 const Numeric& dL0_dT=0.0,
                                 const Numeric& ddF0_dT=0.0,
                                 const Range& df_range=Range(joker));
  
  void set_doppler(ComplexVectorView F,
                   ArrayOfComplexVector& dF,
                   ConstVectorView f_grid,
                   const Numeric& zeeman_df,
                   const Numeric& magnetic_magnitude,
                   const Numeric& F0_noshift,
                   const Numeric& GD_div_F0,
                   const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                   const Numeric& dGD_div_F0_dT=0.0,
                   const Range& df_range=Range(joker));
  
  void set_faddeeva_from_full_linemixing(ComplexVectorView F,
                                         ArrayOfComplexVector& dF,
                                         ConstVectorView f_grid,
                                         const Complex& eigenvalue_no_shift,
                                         const Numeric& GD_div_F0,
                                         const Numeric& L0,
                                         const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                         const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                         const Numeric& dGD_div_F0_dT=0.0,
                                         const Complex& deigenvalue_dT=0.0,
                                         const Numeric& dL0_dT=0.0);
  
  void apply_linemixing_scaling(ComplexVectorView F,
                                ArrayOfComplexVector& dF,
                                const Numeric& Y,
                                const Numeric& G,
                                const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                const Numeric& dY_dT=0.0,
                                const Numeric& dG_dT=0.0,
                                const Range& df_range=Range(joker));
  
  void apply_rosenkranz_quadratic_scaling(ComplexVectorView F,
                                          ArrayOfComplexVector& dF,
                                          ConstVectorView f_grid,
                                          const Numeric& F0,
                                          const Numeric& T,
                                          const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                          const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                          const Range& df_range=Range(joker));
  
  void apply_VVH_scaling(ComplexVectorView F,
                         ArrayOfComplexVector& dF,
                         ConstVectorView f_grid,
                         const Numeric& F0,
                         const Numeric& T,
                         const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                         const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                         const Range& df_range=Range(joker));
  
  void apply_VVW_scaling(ComplexVectorView F,
                         ArrayOfComplexVector& dF,
                         ConstVectorView f_grid,
                         const Numeric& F0,
                         const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                         const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                         const Range& df_range=Range(joker));
  
  void apply_linestrength_scaling(ComplexVectorView F,
                                  ArrayOfComplexVector& dF,
                                  const Numeric& S0,
                                  const Numeric& isotopic_ratio,
                                  const Numeric& QT,
                                  const Numeric& QT0,
                                  const Numeric& K1,
                                  const Numeric& K2,
                                  const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                  const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                  const Numeric& dQT_dT=0.0,
                                  const Numeric& dK1_dT=0.0,
                                  const Numeric& dK2_dT=0.0,
                                  const Numeric& dK2_dF0=0.0,
                                  const Range& df_range=Range(joker));
  
  void set_nonlte_source_and_apply_absorption_scaling(ComplexVectorView F,
                                                      ArrayOfComplexVector& dF,
                                                      ComplexVectorView N,
                                                      ArrayOfComplexVector& dN,
                                                      const Numeric& K3=1.0,
                                                      const Numeric& K4=1.0,
                                                      const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                                      const QuantumIdentifier& quantum_identity=QuantumIdentifier(), 
                                                      const Numeric& dK3_dT=0.0, 
                                                      const Numeric& dK4_dT=0.0,
                                                      const Numeric& dK3_dF0=0.0, 
                                                      const Numeric& dK3_dTl=0.0, 
                                                      const Numeric& dK3_dTu=0.0, 
                                                      const Numeric& dK4_dTu=0.0,
                                                      const Range& df_range=Range(joker));
  
  void apply_linestrength_from_full_linemixing(ComplexVectorView F,
                                               ArrayOfComplexVector& dF,
                                               const Numeric& F0,
                                               const Numeric& T,
                                               const Complex& S_LM,
                                               const Numeric& isotopic_ratio,
                                               const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                               const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                               const Complex& dS_LM_dT=0.0);
  
  void apply_dipole(ComplexVectorView F,
                    ArrayOfComplexVector& dF,
                    const Numeric& F0,
                    const Numeric& T,
                    const Numeric& d0,
                    const Numeric& rho,
                    const Numeric& isotopic_ratio,
                    const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                    const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                    const Numeric& drho_dT=0.0);
  
  void apply_pressurebroadening_jacobian_scaling(ArrayOfComplexVector& dF,
                                                 const PropmatPartialsData& derivatives_data,
                                                 const QuantumIdentifier& quantum_identity,
                                                 const ComplexVector& dgamma,
                                                 const Range& df_range=Range(joker));
  
  void apply_linemixing_jacobian_scaling(ArrayOfComplexVector& dF,
                                         const PropmatPartialsData& derivatives_data,
                                         const QuantumIdentifier& quantum_identity,
                                         const ComplexVector& dlm,
                                         const Range& df_range=Range(joker));
  
  Numeric DopplerConstant(const Numeric T, const Numeric mass);
  
  Numeric dDopplerConstant_dT(const Numeric T, const Numeric mass);
  
  void set_cross_section_for_single_line(ComplexVectorView F,
                                         ArrayOfComplexVector& dF,
                                         ComplexVectorView N,
                                         ArrayOfComplexVector& dN,
                                         Range& this_xsec_range,
                                         const PropmatPartialsData& derivatives_data,
                                         const LineRecord& line,
                                         ConstVectorView f_grid,
                                         ConstVectorView volume_mixing_ratio_of_all_species,
                                         ConstVectorView nlte_distribution,
                                         const Numeric& pressure,
                                         const Numeric& temperature,
                                         const Numeric& doppler_constant,
                                         const Numeric& partial_pressure,
                                         const Numeric& isotopologue_ratio,
                                         const Numeric& magnetic_magnitude,
                                         const Numeric& ddoppler_constant_dT,
                                         const Numeric& pressure_limit_for_linemixing,
                                         const Numeric& partition_function_at_temperature,
                                         const Numeric& dpartition_function_at_temperature_dT,
                                         const Numeric& partition_function_at_line_temperature,
                                         const ArrayOfIndex& broad_spec_locations,
                                         const Index& this_species_location_in_tags,
                                         const Index& water_index_location_in_tags,
                                         const Verbosity& verbosity,
                                         const bool cutoff_call=false,
                                         const Index binary_speedup=false);
  
  void apply_cutoff(ComplexVectorView F,
                    ArrayOfComplexVector& dF,
                    ComplexVectorView N,
                    ArrayOfComplexVector& dN,
                    const PropmatPartialsData& derivatives_data,
                    const LineRecord& line,
                    ConstVectorView volume_mixing_ratio_of_all_species,
                    ConstVectorView nlte_distribution,
                    const Numeric& pressure,
                    const Numeric& temperature,
                    const Numeric& doppler_constant,
                    const Numeric& partial_pressure,
                    const Numeric& isotopologue_ratio,
                    const Numeric& magnetic_magnitude,
                    const Numeric& ddoppler_constant_dT,
                    const Numeric& pressure_limit_for_linemixing,
                    const Numeric& partition_function_at_temperature,
                    const Numeric& dpartition_function_at_temperature_dT,
                    const Numeric& partition_function_at_line_temperature,
                    const ArrayOfIndex& broad_spec_locations,
                    const Index& this_species_location_in_tags,
                    const Index& water_index_location_in_tags,
                    const Range& df_range,
                    const Verbosity& verbosity);
  
  bool find_cutoff_ranges(Range& range,
                          ConstVectorView f_grid,
                          const Numeric& F0,
                          const Numeric& cutoff);
  
  void apply_linestrength_from_nlte_level_distributions(ComplexVectorView F, 
                                                        ArrayOfComplexVector& dF, 
                                                        ComplexVectorView N, 
                                                        ArrayOfComplexVector& dN, 
                                                        const Numeric& r1,
                                                        const Numeric& r2,
                                                        const Numeric& g1,
                                                        const Numeric& g2,
                                                        const Numeric& A21,
                                                        const Numeric& F0,
                                                        const Numeric& T,
                                                        const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                                        const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                                        const Range& df_range=Range(joker));
  
  Range binary_range(const ConstVectorView f, const Index& i, const bool full=false);
  
  Range speedup_binary_range(const ConstVectorView f, const Numeric& u, const Numeric& l);
  
  void interp_up_inside_binary_range(ComplexVectorView f, const Range& l);
  
  void interp_to_boundary_of_binary_range(ComplexVectorView f, const Index upper_boundary, const Index lower_boundary);
  
  Numeric speedup_distance_binary_range(const Index, const Numeric, const Numeric, const Numeric, const Numeric);
  
  void find_boundary_of_binary_range(Index& u, Index& l, ConstVectorView f, const Numeric C, 
                                     const Numeric G0, const Numeric F0, const Numeric GD_div_F0,
                                     const Index binary_speedup);
};

#endif //linefunctions_h

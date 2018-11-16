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
#include "jacobian.h"
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
                     const Numeric& magnetic_magnitude,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const Index& this_species,
                     const Index& zeeman_index);
  
  void set_lorentz(ComplexVectorView F, 
                   ComplexMatrixView dF, 
                   ConstVectorView f_grid, 
                   const Numeric& zeeman_df, 
                   const Numeric& magnetic_magnitude, 
                   const Numeric& F0_noshift, 
                   const Numeric& G0, 
                   const Numeric& L0, 
                   const Numeric& dF0,
                   const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                   const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(), 
                   const Numeric& dG0_dT=0.0, 
                   const Numeric& dL0_dT=0.0,
                   const Numeric& ddF0_dT=0.0,
                   const LineFunctionData::Output& dVMR=NoLineFunctionDataOutput());
  
  void set_htp(ComplexVectorView F,
               ComplexMatrixView dF,
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
               const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
               const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
               const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
               const Numeric& dGD_div_F0_dT=0.0,
               const Numeric& dG0_dT=0.0,
               const Numeric& dL0_dT=0.0,
               const Numeric& dG2_dT=0.0,
               const Numeric& dL2_dT=0.0,
               const Numeric& deta_dT=0.0,
               const Numeric& dFVC_dT=0.0);
  
  void set_voigt(ComplexVectorView F,
                 ComplexMatrixView dF,
                 ConstVectorView f_grid,
                 const Numeric& zeeman_df,
                 const Numeric& magnetic_magnitude,
                 const Numeric& F0_noshift,
                 const Numeric& GD_div_F0,
                 const Numeric& G0,
                 const Numeric& L0,
                 const Numeric& dF0,
                 const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                 const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                 const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                 const Numeric& dGD_div_F0_dT=0.0,
                 const Numeric& dG0_dT=0.0,
                 const Numeric& dL0_dT=0.0,
                 const Numeric& ddF0_dT=0.0,
                 const LineFunctionData::Output& dVMR=NoLineFunctionDataOutput());
  
  void set_doppler(ComplexVectorView F,
                   ComplexMatrixView dF,
                   ConstVectorView f_grid,
                   const Numeric& zeeman_df,
                   const Numeric& magnetic_magnitude,
                   const Numeric& F0_noshift,
                   const Numeric& GD_div_F0,
                   const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                   const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                   const Numeric& dGD_div_F0_dT=0.0);
  
  void set_voigt_from_full_linemixing(ComplexVectorView F,
                                      ComplexMatrixView dF,
                                      ConstVectorView f_grid,
                                      const Complex& eigenvalue_no_shift,
                                      const Numeric& GD_div_F0,
                                      const Numeric& L0,
                                      const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                      const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                      const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                      const Numeric& dGD_div_F0_dT=0.0,
                                      const Complex& deigenvalue_dT=0.0,
                                      const Numeric& dL0_dT=0.0);
  
  void apply_linemixing_scaling(ComplexVectorView F,
                                ComplexMatrixView dF,
                                const Numeric& Y,
                                const Numeric& G,
                                const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                const Numeric& dY_dT=0.0,
                                const Numeric& dG_dT=0.0,
                                const LineFunctionData::Output& dVMR=NoLineFunctionDataOutput());
  
  void apply_rosenkranz_quadratic_scaling(ComplexVectorView F,
                                          ComplexMatrixView dF,
                                          ConstVectorView f_grid,
                                          const Numeric& F0,
                                          const Numeric& T,
                                          const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                          const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                          const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  void apply_VVH_scaling(ComplexVectorView F,
                         ComplexMatrixView dF,
                         ConstVectorView f_grid,
                         const Numeric& F0,
                         const Numeric& T,
                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                         const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  void apply_VVW_scaling(ComplexVectorView F,
                         ComplexMatrixView dF,
                         ConstVectorView f_grid,
                         const Numeric& F0,
                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                         const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  void apply_linestrength_scaling(ComplexVectorView F,
                                  ComplexMatrixView dF,
                                  const Numeric& S0,
                                  const Numeric& isotopic_ratio,
                                  const Numeric& QT,
                                  const Numeric& QT0,
                                  const Numeric& K1,
                                  const Numeric& K2,
                                  const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                  const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                  const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                  const Numeric& dQT_dT=0.0,
                                  const Numeric& dK1_dT=0.0,
                                  const Numeric& dK2_dT=0.0,
                                  const Numeric& dK2_dF0=0.0);
  
  void apply_linestrength_scaling_vibrational_nlte(ComplexVectorView F,
                                                   ComplexMatrixView dF,
                                                   ComplexVectorView N,
                                                   ComplexMatrixView dN,
                                                   const Numeric& S0,
                                                   const Numeric& isotopic_ratio,
                                                   const Numeric& QT,
                                                   const Numeric& QT0,
                                                   const Numeric& K1,
                                                   const Numeric& K2,
                                                   const Numeric& K3,
                                                   const Numeric& K4,
                                                   const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                   const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                                   const Numeric& dQT_dT=0.0,
                                                   const Numeric& dK1_dT=0.0,
                                                   const Numeric& dK2_dT=0.0,
                                                   const Numeric& dK2_dF0=0.0,
                                                   const Numeric& dK3_dT=0.0,
                                                   const Numeric& dK3_dF0=0.0,
                                                   const Numeric& dK3_dTl=0.0,
                                                   const Numeric& dK3_dTu=0.0,
                                                   const Numeric& dK4_dT=0.0,
                                                   const Numeric& dK4_dTu=0.0);
  
  void apply_linestrength_from_full_linemixing(ComplexVectorView F,
                                               ComplexMatrixView dF,
                                               const Numeric& F0,
                                               const Numeric& T,
                                               const Complex& S_LM,
                                               const Numeric& isotopic_ratio,
                                               const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                               const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                               const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                               const Complex& dS_LM_dT=0.0);
  
  void apply_dipole(ComplexVectorView F,
                    ComplexMatrixView dF,
                    const Numeric& F0,
                    const Numeric& T,
                    const Numeric& d0,
                    const Numeric& rho,
                    const Numeric& isotopic_ratio,
                    const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                    const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                    const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                    const Numeric& drho_dT=0.0);
  
  void apply_linefunctiondata_jacobian_scaling(ComplexMatrixView dF,
                                              const ArrayOfRetrievalQuantity& derivatives_data,
                                              const ArrayOfIndex& derivatives_data_position,
                                              const QuantumIdentifier& quantum_identity,
                                              const LineRecord& line,
                                              const Numeric& T,
                                              const Numeric& P,
                                              const Index& this_species,
                                              const ConstVectorView& vmrs,
                                              const ArrayOfArrayOfSpeciesTag& species);
  
  Numeric DopplerConstant(const Numeric& T, const Numeric& mass);
  
  Numeric dDopplerConstant_dT(const Numeric& T, const Numeric& mass);
  
  void set_cross_section_for_single_line(ComplexVectorView F,
                                         ComplexMatrixView dF,
                                         ComplexVectorView N,
                                         ComplexMatrixView dN,
                                         Range& this_xsec_range,
                                         const ArrayOfRetrievalQuantity& derivatives_data,
                                         const ArrayOfIndex& derivatives_data_position,
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
                                         const Numeric& partition_function_at_temperature,
                                         const Numeric& dpartition_function_at_temperature_dT,
                                         const Numeric& partition_function_at_line_temperature,
                                         const ArrayOfArrayOfSpeciesTag& abs_species,
                                         const Index& this_species_location_in_tags,
                                         const Index& zeeman_index,
                                         const Verbosity& verbosity,
                                         const bool cutoff_call=false);
  
  void apply_cutoff(ComplexVectorView F,
                    ComplexMatrixView dF,
                    ComplexVectorView N,
                    ComplexMatrixView dN,
                    const ArrayOfRetrievalQuantity& derivatives_data,
                    const ArrayOfIndex& derivatives_data_position,
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
                    const Numeric& partition_function_at_temperature,
                    const Numeric& dpartition_function_at_temperature_dT,
                    const Numeric& partition_function_at_line_temperature,
                    const ArrayOfArrayOfSpeciesTag& abs_species,
                    const Index& this_species_location_in_tags,
                    const Index& zeeman_index,
                    const Verbosity& verbosity);
  
  bool find_cutoff_ranges(Range& range,
                          const ConstVectorView& f_grid,
                          const Numeric& F0,
                          const Numeric& cutoff);
  
  void apply_linestrength_from_nlte_level_distributions(ComplexVectorView F, 
                                                        ComplexMatrixView dF, 
                                                        ComplexVectorView N, 
                                                        ComplexMatrixView dN, 
                                                        const Numeric& r1,
                                                        const Numeric& r2,
                                                        const Numeric& g1,
                                                        const Numeric& g2,
                                                        const Numeric& A21,
                                                        const Numeric& F0,
                                                        const Numeric& T,
                                                        const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                        const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                                        const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  Index first_binary_level(const Numeric& dg, const Numeric& df, const Index N);
  
  Numeric binary_central_step_size(const Numeric& ds, const Index i, const Index i0);
  
  Index binary_frequency_position(const Numeric& df, const Numeric& dg, const Numeric& f, const Numeric& F0, const Index n, const Index dn);
  
  ArrayOfArrayOfIndex binary_boundaries(const Numeric& F0, const ConstVectorView& f_grid, const Numeric& G0,
                                        const Numeric& GD_div_F0, const Numeric& binary_speedup_coef,
                                        const Index binary_speedup, const Index binary_frequency_count,
                                        const Index steps=3);
  
  void binary_interpolation(ComplexVectorView f, const ArrayOfArrayOfIndex& binary_bounds);
  
  Range binary_level_range(const ArrayOfArrayOfIndex& binary_bounds, const Index nf, const Index i, const bool lower_range);
  
  };

#endif //linefunctions_h

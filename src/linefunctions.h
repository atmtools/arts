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
  constexpr Index ExpectedDataSize(){return 2;}
  
  void set_lineshape(Eigen::Ref<Eigen::VectorXcd> F, 
                     const Eigen::Ref<const Eigen::VectorXd> f_grid, 
                     const LineRecord& line, 
                     const ConstVectorView vmrs, 
                     const Numeric& temperature, 
                     const Numeric& pressure, 
                     const Numeric& magnetic_magnitude,
                     const ArrayOfArrayOfSpeciesTag& abs_species,
                     const Index& this_species,
                     const Index& zeeman_index);
  
  void set_lorentz(Eigen::Ref<Eigen::VectorXcd> F, 
                   Eigen::Ref<Eigen::MatrixXcd> dF, 
                   Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                   const Eigen::Ref<const Eigen::VectorXd> f_grid, 
                   const Numeric& zeeman_df, 
                   const Numeric& magnetic_magnitude, 
                   const Numeric& F0_noshift, 
                   const LineFunctionDataOutput& X,
                   const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                   const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(), 
                   const LineFunctionDataOutput& dT=NoLineFunctionDataOutput(),
                   const LineFunctionDataOutput& dxdVMR=NoLineFunctionDataOutput());
  
  void set_htp(Eigen::Ref<Eigen::VectorXcd> F,
               Eigen::Ref<Eigen::MatrixXcd> dF,
               const Eigen::Ref<const Eigen::VectorXd> f_grid,
               const Numeric& zeeman_df,
               const Numeric& magnetic_magnitude,
               const Numeric& F0_noshift,
               const Numeric& GD_div_F0,
               const LineFunctionDataOutput& X,
               const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
               const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
               const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
               const Numeric& dGD_div_F0_dT=0.0,
               const LineFunctionDataOutput& dT=NoLineFunctionDataOutput(),
               const LineFunctionDataOutput& dVMR=NoLineFunctionDataOutput());
  
  void set_voigt(Eigen::Ref<Eigen::VectorXcd> F,
                 Eigen::Ref<Eigen::MatrixXcd> dF,
                 Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                 const Eigen::Ref<const Eigen::VectorXd> f_grid,
                 const Numeric& zeeman_df,
                 const Numeric& magnetic_magnitude,
                 const Numeric& F0_noshift,
                 const Numeric& GD_div_F0,
                 const LineFunctionDataOutput& X,
                 const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                 const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                 const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                 const Numeric& dGD_div_F0_dT=0.0,
                 const LineFunctionDataOutput& dT=NoLineFunctionDataOutput(),
                 const LineFunctionDataOutput& dVMR=NoLineFunctionDataOutput());
  
  void set_doppler(Eigen::Ref<Eigen::VectorXcd> F,
                   Eigen::Ref<Eigen::MatrixXcd> dF,
                   Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                   const Eigen::Ref<const Eigen::VectorXd> f_grid,
                   const Numeric& zeeman_df,
                   const Numeric& magnetic_magnitude,
                   const Numeric& F0_noshift,
                   const Numeric& GD_div_F0,
                   const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                   const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                   const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                   const Numeric& dGD_div_F0_dT=0.0);
  
  void set_voigt_from_full_linemixing(Eigen::Ref<Eigen::VectorXcd> F,
                                      Eigen::Ref<Eigen::MatrixXcd> dF,
                                      const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                      const Complex& eigenvalue_no_shift,
                                      const Numeric& GD_div_F0,
                                      const Numeric& L0,
                                      const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                      const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                      const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                      const Numeric& dGD_div_F0_dT=0.0,
                                      const Complex& deigenvalue_dT=0.0,
                                      const Numeric& dL0_dT=0.0);
  
  void apply_linemixing_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                                Eigen::Ref<Eigen::MatrixXcd> dF,
                                const LineFunctionDataOutput& X,
                                const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                const LineFunctionDataOutput& dT=NoLineFunctionDataOutput(),
                                const LineFunctionDataOutput& dVMR=NoLineFunctionDataOutput());
  
  void apply_rosenkranz_quadratic_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                                          Eigen::Ref<Eigen::MatrixXcd> dF,
                                          const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                          const Numeric& F0,
                                          const Numeric& T,
                                          const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                          const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                          const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  void apply_VVH_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                         Eigen::Ref<Eigen::MatrixXcd> dF,
                         Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                         const Eigen::Ref<const Eigen::VectorXd> f_grid,
                         const Numeric& F0,
                         const Numeric& T,
                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                         const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  void apply_VVW_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                         Eigen::Ref<Eigen::MatrixXcd> dF,
                         const Eigen::Ref<const Eigen::VectorXd> f_grid,
                         const Numeric& F0,
                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                         const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  
  void apply_linestrength_scaling_by_lte(Eigen::Ref<Eigen::VectorXcd> F,
                                         Eigen::Ref<Eigen::MatrixXcd> dF,
                                         Eigen::Ref<Eigen::VectorXcd> N,
                                         Eigen::Ref<Eigen::MatrixXcd> dN,
                                         const LineRecord& line,
                                         const Numeric& T,
                                         const Numeric& isotopic_ratio,
                                         const Numeric& zeeman_scaling,
                                         const Numeric& QT,
                                         const Numeric& QT0,
                                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                         const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                         const Numeric& dQT_dT=0.0);
  
  void apply_linestrength_scaling_by_vibrational_nlte(Eigen::Ref<Eigen::VectorXcd> F,
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
                                                      const Numeric& zeeman_scaling,
                                                      const Numeric& QT,
                                                      const Numeric& QT0,
                                                      const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                      const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                                      const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                                      const Numeric& dQT_dT=0.0);
  
  void apply_linestrength_from_full_linemixing(Eigen::Ref<Eigen::VectorXcd> F,
                                               Eigen::Ref<Eigen::MatrixXcd> dF,
                                               const Numeric& F0,
                                               const Numeric& T,
                                               const Complex& S_LM,
                                               const Numeric& isotopic_ratio,
                                               const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                               const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                               const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                                               const Complex& dS_LM_dT=0.0);
  
  void apply_dipole(Eigen::Ref<Eigen::VectorXcd> F,
                    Eigen::Ref<Eigen::MatrixXcd> dF,
                    const Numeric& F0,
                    const Numeric& T,
                    const Numeric& d0,
                    const Numeric& rho,
                    const Numeric& isotopic_ratio,
                    const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                    const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                    const QuantumIdentifier& quantum_identity=QuantumIdentifier(),
                    const Numeric& drho_dT=0.0);
  
  void apply_linefunctiondata_jacobian_scaling(Eigen::Ref<Eigen::MatrixXcd> dF,
                                              const ArrayOfRetrievalQuantity& derivatives_data,
                                              const ArrayOfIndex& derivatives_data_position,
                                              const QuantumIdentifier& quantum_identity,
                                              const LineRecord& line,
                                              const Numeric& T,
                                              const Numeric& P,
                                              const Index& this_species,
                                              const ConstVectorView& vmrs,
                                              const ArrayOfArrayOfSpeciesTag& species);
  
  Numeric DopplerConstant(Numeric T, Numeric mass);
  
  Numeric dDopplerConstant_dT(const Numeric& T, const Numeric& dc);
  
  void set_cross_section_for_single_line(Eigen::Ref<Eigen::VectorXcd>  F_full,
                                         Eigen::Ref<Eigen::MatrixXcd> dF_full,
                                         Eigen::Ref<Eigen::VectorXcd>  N_full,
                                         Eigen::Ref<Eigen::MatrixXcd> dN_full,
                                         Eigen::Ref<Eigen::MatrixXcd> data_block,
                                         Index& start_cutoff,
                                         Index& nelem_cutoff,
                                         const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                         const LineRecord& line,
                                         const ArrayOfRetrievalQuantity& derivatives_data,
                                         const ArrayOfIndex& derivatives_data_position,
                                         const ConstVectorView volume_mixing_ratio_of_all_species,
                                         const ConstVectorView nlte_distribution,
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
                                         const bool cutoff_call=false);
  
  void apply_cutoff(Eigen::Ref<Eigen::VectorXcd> F,
                    Eigen::Ref<Eigen::MatrixXcd> dF,
                    Eigen::Ref<Eigen::VectorXcd> N,
                    Eigen::Ref<Eigen::MatrixXcd> dN,
                    const ArrayOfRetrievalQuantity& derivatives_data,
                    const ArrayOfIndex& derivatives_data_position,
                    const LineRecord& line,
                    const ConstVectorView volume_mixing_ratio_of_all_species,
                    const ConstVectorView nlte_distribution,
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
                    const Index& zeeman_index);
  
  void find_cutoff_ranges(Index& start_cutoff,
                          Index& nelem_cutoff,
                          const Eigen::Ref<const Eigen::VectorXd> f_grid,
                          const Numeric& F0,
                          const Numeric& cutoff);
  
  void apply_linestrength_from_nlte_level_distributions(Eigen::Ref<Eigen::VectorXcd> F, 
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
                                                        const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                        const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                                        const QuantumIdentifier& quantum_identity=QuantumIdentifier());
  };

#endif //linefunctions_h

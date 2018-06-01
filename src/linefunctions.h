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
  class SingleLevelLineData {
  public:
    SingleLevelLineData(const LineRecord& line,
                        const ConstVectorView vmrs,
                        const ConstVectorView nlte_distribution,
                        const Numeric& T,
                        const Numeric& P,
                        const Numeric& H,
                        const Numeric& lm_p_lim,
                        const Numeric& QT,
                        const Numeric& dQTdT,
                        const Numeric& QT0,
                        const Numeric& isotopic_ratio,
                        const ArrayOfIndex& broadening_species,
                        const Index this_species,
                        const Index water_species,
                        const Index zeeman_index,
                        const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                        const ArrayOfIndex& derivatives_data_position=ArrayOfIndex());
    
    enum class SpectroscopyDerivivatives : Index {Upper, Lower, Both, FullTransition, None};
    
    inline Complex doppler_z(const Numeric& f0, const Numeric& f) const noexcept {return Complex(minvGD * (f - f0 - mZ), 0.0);}
    
    inline Complex lorentz_z(const Numeric& f0, const Numeric& f) const noexcept {return Complex(f - f0 - mC0.imag() - mZ - mDV, -mC0.real());}
    inline Complex lorentz_dT() const noexcept {return Complex(mdC0dT.imag() + mdDVdT, -mdC0dT.imag());;}
    
    inline Complex voigt_z(const Numeric& f0, const Numeric& f) const noexcept {return Complex(f - f0 - mC0.imag() - mZ - mDV, mC0.real()) * minvGD;}
    inline Complex voigt_dT(const Complex& z) const noexcept {return (Complex(-mdC0dT.imag() - mdDVdT, mdC0dT.real()) - z*mdGDdT) * minvGD;}
    inline Complex voigt_dF0(const Complex& z) const noexcept {return -z * mGD_div_F0 - 1.0;}
    
    inline Complex htp_x(const Numeric& f0, const Numeric& f) const noexcept;
    inline Complex htp_dxdT(const Complex& x) const noexcept;
    inline Complex htp_dxdf() const noexcept;
    inline Complex htp_dxdf0() const noexcept;
    inline Complex htp_dxdmag(const Numeric& zeeman_df) const noexcept;
    inline Complex htp_dxdC0() const noexcept;
    inline Complex htp_dxdC2(const Complex& x) const noexcept;
    inline Complex htp_dxdFVC() const noexcept;
    inline Complex htp_dxdeta(const Complex& x) const noexcept;
    
    inline Numeric xy_ratio(const Complex& x) const noexcept {return abs(x) / abs(my);}
    inline Complex htp_z1(const Complex& x, const Numeric& ratio_xy) const noexcept;
    inline Complex htp_dz1dt(const Complex& z1, const Complex& dx, const Numeric& ratio_xy, const bool for_temperature=false) const noexcept;
    inline Complex htp_z2(const Complex& x, const Complex& z1, const Numeric& ratio_xy) const noexcept;
    inline Complex htp_dz2dt(const Complex& z2, const Complex& dz1, const Complex& dx, const Numeric& ratio_xy, const bool for_temperature=false) const noexcept;
    
    inline Complex htp_w1(const Complex& z1) const noexcept;
    inline Complex htp_dw1_over_dz1(const Complex& z1, const Complex& w1) const noexcept;
    inline Complex htp_w2(const Complex& z2) const noexcept;
    inline Complex htp_dw2_over_dz2(const Complex& z2, const Complex& w2) const noexcept;
    
    inline Complex htp_A(const Complex& w1, const Complex& w2, const Complex& z1, const Numeric& ratio_xy) const noexcept;
    inline Complex htp_dAdt(const Complex& A,  const Complex& w1, const Complex& dw1, const Complex& dw2, const Complex& z1, const Complex& dz1, const Numeric& ratio_xy, const bool for_temperature=false) const noexcept;
    inline Complex htp_B(const Complex& w1, const Complex& w2, const Complex& z1, const Complex& z2, const Numeric& ratio_xy) const noexcept;
    inline Complex htp_dBdt(const Complex& w1, const Complex& dw1, const Complex& w2, const Complex& dw2, const Complex& z1, const Complex& dz1, const Complex& z2, const Complex& dz2, const Numeric& ratio_xy, const bool for_temperature=false) const noexcept;
    inline Complex htp_G(const Complex& A, const Complex& B) const noexcept;
    inline Complex htp_dGdt(const Complex& A, const Complex& dA, const Complex& dB, const bool for_temperature=false) const noexcept;
    
    inline Complex dw_over_dz(const Complex& z, const Complex& w) const noexcept;
    inline Complex scale_w(const Complex& w) const noexcept;
    
    inline bool no_more_pressure_jacs(const Index i) const noexcept {return i >= mpressure_derivatives.nelem();}
    inline Complex dgamma(const Index i) const noexcept {return mpressure_derivatives[i];}
    
    inline bool no_more_linemixing_jacs(const Index i) const noexcept {return i >= mlinemixing_derivatives.nelem();}
    inline Complex dlm(const Index i) const noexcept {return mlinemixing_derivatives[i];}
    
    const Numeric& invGD() const {return minvGD;}
    const Numeric& dGDdT() const {return mdGDdT;}
    const Numeric& GD_div_F0() const {return mGD_div_F0;}
    const Numeric& norm() const {return mnorm;}
    const Numeric& dnormdT() const {return mdnormdT;}
    const Numeric& dnormdf0() const {return mdnormdf0;}
    const Complex& LM() const {return mLM;}
    const Complex& dLMdT() const {return mdLMdT;}
    const Numeric& S() const {return mS;}
    const Numeric& dSdT() const {return mdSdT;}
    const Numeric& dSdf0() const {return mdSdf0;}
    const Numeric& nlte_source_factor() const {return mnlte_src;}
    const Numeric& dnlte_source_factordT() const {return mdnlte_srcdT;}
    const Numeric& dnlte_source_factordf0() const {return mdnlte_srcdf0;}
    const Numeric& dnlte_source_factordlow() const {return mdnlte_srcdlow;}
    const Numeric& dnlte_source_factordupp() const {return mdnlte_srcdupp;}
    const Numeric& nlte_absorption_factor() const {return mnlte_abs;}
    const Numeric& dnlte_absorption_factordT() const {return mdnlte_absdT;}
    const Numeric& dnlte_absorption_factordf0() const {return mdnlte_absdf0;}
    const Numeric& dnlte_absorption_factordupp() const {return mdnlte_absdupp;}
    const Numeric& dnlte_absorption_factordlow() const {return mdnlte_absdlow;}
    Numeric G0() const {return mC0.real();}
    Numeric L0() const {return mC0.imag();}
    const Numeric& LM_DF() const {return mDV;}
    const SpectroscopyDerivivatives& operator()(const Index& id) const {return mspectroscopy_derivatives[id];}
    
    friend std::ostream& operator<<(std::ostream& os, const SingleLevelLineData& slld);
  private:
    Complex mC0, mdC0dT, mC2, mdC2dT, mLM, mdLMdT;  // Pressure broadening and line mixing in complex terms
    Numeric mFVC, mdFVCdT, meta;  // Pressure broadening terms for HTP
    Numeric mDV, mdDVdT;  // Line mixing shifts
    ComplexVector mpressure_derivatives, mlinemixing_derivatives;  // values for the pressure derivatives
    Array<SpectroscopyDerivivatives> mspectroscopy_derivatives;  // flag for spectroscopic parameters of given line
    Numeric mGD_div_F0, minvGD, mdGDdT;  // Doppler factors (used a lot for scaling)
    Numeric mZ;  // Zeeman shift
    Complex mC0t, mdC0tdT, mC2t, mdC2tdT, my, msqrty, mdydT;
    bool mC2t_is_zero;  // HTP test
    Numeric mnorm, mdnormdT, mdnormdf0; // Normalization factors
    Numeric mS, mdSdT, mdSdf0;
    Numeric mnlte_abs, mdnlte_absdT, mdnlte_absdf0, mdnlte_absdlow, mdnlte_absdupp;
    Numeric mnlte_src, mdnlte_srcdT, mdnlte_srcdf0, mdnlte_srcdlow, mdnlte_srcdupp;
  };
  
  std::ostream& operator<<(std::ostream& os, const SingleLevelLineData& slld);
  
  void set_lineshape_from_level_line_data(Complex& F,
                                          Complex& N,
                                          ComplexVectorView dF,
                                          ComplexVectorView dN,
                                          const Numeric& f,
                                          const Numeric& T,
                                          const SingleLevelLineData& level_line_data,
                                          const LineRecord& line,
                                          const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                          const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void set_doppler_from_level_line_data(Complex& F,
                                        ComplexVectorView dF,
                                        const Numeric& f,
                                        const Numeric& f0,
                                        const Numeric& dZdH,
                                        const SingleLevelLineData& level_line_data,
                                        const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                        const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void set_lorentz_from_level_line_data(Complex& F,
                                        ComplexVectorView dF,
                                        const Numeric& f,
                                        const Numeric& f0,
                                        const Numeric& dZdH,
                                        const SingleLevelLineData& level_line_data,
                                        const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                        const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void set_voigt_from_level_line_data(Complex& F,
                                      ComplexVectorView dF,
                                      const Numeric& f,
                                      const Numeric& f0,
                                      const Numeric& dZdH,
                                      const SingleLevelLineData& level_line_data,
                                      const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                      const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void set_htp_from_level_line_data(Complex& F,
                                    ComplexVectorView dF,
                                    const Numeric& f,
                                    const Numeric& f0,
                                    const Numeric& dZdH,
                                    const SingleLevelLineData& level_line_data,
                                    const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                    const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_rosenkranz_quadratic_scaling_from_level_data(Complex& F,
                                                          ComplexVectorView dF,
                                                          const Numeric& f,
                                                          const SingleLevelLineData& level_line_data,
                                                          const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                          const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_VVH_scaling_from_level_data(Complex& F, ComplexVectorView dF,
                                         const Numeric& f, const Numeric& T, const SingleLevelLineData& level_line_data,
                                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_VVW_scaling_from_level_data(Complex& F, ComplexVectorView dF,
                                         const Numeric& f, const SingleLevelLineData& level_line_data,
                                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_linemixing_from_level_data(Complex& F, ComplexVectorView dF,
                                        const SingleLevelLineData& level_line_data,
                                        const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                        const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;

  void apply_pressurebroadening_jacobian_scaling_from_level_data(ComplexVectorView dF,
                                                                 const SingleLevelLineData& level_line_data,
                                                                 const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                                 const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;

  void apply_linemixing_jacobian_scaling_from_level_data(ComplexVectorView dF,
                                                         const SingleLevelLineData& level_line_data,
                                                         const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                         const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_LTE_linestrength_from_level_data(Complex& F, ComplexVectorView dF,
                                              const SingleLevelLineData& level_line_data,
                                              const LineRecord& line,
                                              const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                              const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_NLTE_vibrational_temperature_linestrength_from_level_data(Complex& F, Complex& N,
                                                                       ComplexVectorView dF, ComplexVectorView dN,
                                                                       const SingleLevelLineData& level_line_data,
                                                                       const LineRecord& line,
                                                                       const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                                       const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void apply_NLTE_population_distribution_linestrength_from_level_data(Complex& F, Complex& N,
                                                                       ComplexVectorView dF, ComplexVectorView dN,
                                                                       const SingleLevelLineData& level_line_data,
                                                                       const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                                       const ArrayOfIndex& derivatives_data_position=ArrayOfIndex()) noexcept;
  
  void set_lineshape(ComplexVectorView F, 
                     const LineRecord& line, 
                     ConstVectorView f_grid, 
                     ConstVectorView vmrs, 
                     const Numeric& temperature, 
                     const Numeric& pressure, 
                     const Numeric& pressure_limit_for_linemixing, 
                     const Numeric& magnetic_magnitude,
                     const ArrayOfIndex& broad_spec_locations,
                     const Index& this_species,
                     const Index& water_species,
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
                   const Numeric& ddF0_dT=0.0);
  
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
                 const Numeric& ddF0_dT=0.0);
  
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
                                const Numeric& dG_dT=0.0);
  
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
  
  void set_nonlte_source_and_apply_absorption_scaling(ComplexVectorView F,
                                                      ComplexMatrixView dF,
                                                      ComplexVectorView N,
                                                      ComplexMatrixView dN,
                                                      const Numeric& K3=1.0,
                                                      const Numeric& K4=1.0,
                                                      const ArrayOfRetrievalQuantity& derivatives_data=ArrayOfRetrievalQuantity(),
                                                      const ArrayOfIndex& derivatives_data_position=ArrayOfIndex(),
                                                      const QuantumIdentifier& quantum_identity=QuantumIdentifier(), 
                                                      const Numeric& dK3_dT=0.0, 
                                                      const Numeric& dK4_dT=0.0,
                                                      const Numeric& dK3_dF0=0.0, 
                                                      const Numeric& dK3_dTl=0.0, 
                                                      const Numeric& dK3_dTu=0.0, 
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
  
  void apply_pressurebroadening_jacobian_scaling(ComplexMatrixView dF,
                                                 const ArrayOfRetrievalQuantity& derivatives_data,
                                                 const ArrayOfIndex& derivatives_data_pos,
                                                 const QuantumIdentifier& quantum_identity,
                                                 const ComplexVector& dgamma);
  
  void apply_linemixing_jacobian_scaling(ComplexMatrixView dF,
                                         const ArrayOfRetrievalQuantity& derivatives_data,
                                         const ArrayOfIndex& derivatives_data_pos,
                                         const QuantumIdentifier& quantum_identity,
                                         const ComplexVector& dlm);
  
  Numeric DopplerConstant(const Numeric T, const Numeric mass);
  
  Numeric dDopplerConstant_dT(const Numeric T, const Numeric mass);
  
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
                                         const Numeric& pressure_limit_for_linemixing,
                                         const Numeric& partition_function_at_temperature,
                                         const Numeric& dpartition_function_at_temperature_dT,
                                         const Numeric& partition_function_at_line_temperature,
                                         const ArrayOfIndex& broad_spec_locations,
                                         const Index& this_species_location_in_tags,
                                         const Index& water_index_location_in_tags,
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
                    const Numeric& pressure_limit_for_linemixing,
                    const Numeric& partition_function_at_temperature,
                    const Numeric& dpartition_function_at_temperature_dT,
                    const Numeric& partition_function_at_line_temperature,
                    const ArrayOfIndex& broad_spec_locations,
                    const Index& this_species_location_in_tags,
                    const Index& water_index_location_in_tags,
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

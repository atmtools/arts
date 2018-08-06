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
 * Example usage for simple Lorentz line shape with line strength scaled to the correct integration
 * set_lorentz(...)
 * apply_linestrength_scaling(...)
 * 
 * \author Richard Larsson
 * \date   2017-05-16
 */

#include "Faddeeva.hh"
#include "linefunctions.h"
#include "linescaling.h"

// Basic constants 
extern const Numeric PI;
extern const Numeric PLANCK_CONST;
extern const Numeric BOLTZMAN_CONST;
extern const Numeric AVOGADROS_NUMB;
extern const Numeric SPEED_OF_LIGHT;

// Derived constants 
static const Numeric invPI = 1.0 / PI;
static const Numeric sqrtInvPI =  sqrt(invPI);
static const Numeric sqrtPI = sqrt(PI);
static const Numeric C1 = - PLANCK_CONST / BOLTZMAN_CONST;
static const Numeric doppler_const = sqrt(2.0 * BOLTZMAN_CONST * AVOGADROS_NUMB ) / SPEED_OF_LIGHT; 


void Linefunctions::set_lineshape(ComplexVectorView F, 
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
                                  const Index& zeeman_index)
{
  // Pressure broadening terms
  const Numeric partial_pressure = pressure * vmrs[this_species];
  Numeric G0, G2, e, L0, L2, FVC;
  line.PressureBroadening().GetPressureBroadeningParams(
    G0, G2, e, L0, L2, FVC, temperature, line.Ti0(), pressure, partial_pressure, 
    this_species, water_species, broad_spec_locations, vmrs);
  
  // Line mixing terms
  Numeric Y=0, G=0, DV=0;
  line.LineMixing().GetLineMixingParams(Y, G, DV, temperature, pressure,  pressure_limit_for_linemixing);
  
  // Line shape usage remembering variable
  LineShapeType lst = LineShapeType::End;
  
  ComplexMatrix dF(0, 0);
  const Numeric doppler_constant = DopplerConstant(temperature, line.IsotopologueData().Mass());
  
  switch(line.GetLineShapeType())
  {
    case LineShapeType::ByPressureBroadeningData:
      switch(line.PressureBroadening().Type())
      {
        // Use data as per speed dependent air
        case PressureBroadeningData::PB_SD_AIR_VOLUME:
        case PressureBroadeningData::PB_PURELY_FOR_TESTING:
        case PressureBroadeningData::PB_HTP_AIR_VOLUME:
        case PressureBroadeningData::PB_SD_TEST_WATER:
          lst = LineShapeType::HTP;
          set_htp(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, L0, G2, L2, e, FVC);
          break;
          // Use for data that requires air and water Voigt broadening
        case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
          // Use for data that requires planetary Voigt broadening
        case PressureBroadeningData::PB_PLANETARY_BROADENING:
          // Use for data that requires air Voigt broadening
        case PressureBroadeningData::PB_AIR_BROADENING:
        case PressureBroadeningData::PB_VOIGT_TEST_WATER:
          // Above should be all methods of pressure broadening requiring Voigt in ARTS by default
          lst = LineShapeType::Voigt;
          set_voigt(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, L0, DV);
          break;
        case PressureBroadeningData::PB_NONE:
          throw std::runtime_error("Cannot understand the pressure broadening scheme.");
      }
      break;
        case LineShapeType::Doppler:
          lst = LineShapeType::Doppler;
          set_doppler(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant);
          break;
          // This line only needs Hartmann-Tran
        case LineShapeType::HTP:
          set_htp(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, L0, G2, L2, e, FVC);
          lst = LineShapeType::HTP;
          break;
          // This line only needs Lorentz
        case LineShapeType::Lorentz:
          lst = LineShapeType::Lorentz;
          set_lorentz(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), G0, L0, DV);
          break;
          // This line only needs Voigt
        case LineShapeType::Voigt:
          lst = LineShapeType::Voigt;
          set_voigt(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, L0, DV);
          break;
        case LineShapeType::End:
          throw std::runtime_error("Cannot understand the requested line shape type.");
  }
  
  switch(line.GetMirroringType())
  {
    // No mirroring
    case MirroringType::None:
      break;
      // Lorentz mirroring
    case MirroringType::Lorentz:
    {
      // Set the mirroring computational vectors and size them as needed
      ComplexVector Fm(F.nelem());
      ComplexMatrix dFm(0, 0);
      
      set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), G0, -L0, -DV);
      
      // Apply mirroring
      F -= Fm;
    }
    break;
    // Same type of mirroring as before
    case MirroringType::SameAsLineShape:
    {
      // Set the mirroring computational vectors and size them as needed
      ComplexVector Fm(F.nelem());
      ComplexMatrix dFm(0, 0);
      
      switch(lst)
      {
        case LineShapeType::Doppler:
          set_doppler(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant);
          break;
        case LineShapeType::Lorentz:
          set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), G0, -L0, -DV);
          break;
        case LineShapeType::Voigt:
          set_voigt(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant, G0, -L0, -DV);
          break;
        case LineShapeType::HTP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant, G0, -L0, G2, -L2, e, FVC);
          break;
        case LineShapeType::ByPressureBroadeningData:
        case LineShapeType::End:
          throw std::runtime_error("Cannot understand the requested line shape type for mirroring.");
      }
      F -= Fm;
      break;
    }
        case MirroringType::End:
          throw std::runtime_error("Cannot understand the requested mirroring type for mirroring.");
  }
  
  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by and index
  // that is internally interpreted to mean some kind of lineshape normalization
  switch(line.GetLineNormalizationType())
  {
    // No normalization
    case LineNormalizationType::None:
      break;
      // van Vleck and Huber normalization
    case LineNormalizationType::VVH:
      apply_VVH_scaling(F, dF, f_grid, line.F(), temperature);
      break;
      // van Vleck and Weiskopf normalization
    case LineNormalizationType::VVW:
      apply_VVW_scaling(F, dF, f_grid, line.F());
      break;
      // Rosenkranz's Quadratic normalization
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F, dF, f_grid, line.F(), temperature);
      break;
    case LineNormalizationType::End:
      throw std::runtime_error("Cannot understand the requested line normalization type.");
  }
}



/*!
 * Sets the line shape to Lorentz line shape. Normalization is unity.
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param G0 Speed-independent pressure broadening term
 * \param L0 Speed-independent pressure shift term
 * \param dF0 Second order line mixing shift
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
 * \param ddF0_dT Temperature derivative of dF0
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_lorentz(ComplexVectorView F,
                                ComplexMatrixView dF,
                                ConstVectorView f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& G0,
                                const Numeric& L0, 
                                const Numeric& dF0,
                                const ArrayOfRetrievalQuantity& derivatives_data,
                                const ArrayOfIndex& derivatives_data_position,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dG0_dT,
                                const Numeric& dL0_dT,
                                const Numeric& ddF0_dT)
{ 
  // Size of the problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  // The central frequency
  const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude + dF0;
  
  // Constant part of the denominator
  const Complex denom0 = Complex(G0, F0);
  
  Complex d, denom;
  
  #pragma omp simd
  for(Index iv = 0; iv < nf; iv++)
  {
    denom = 1.0 / ((denom0 - Complex(0.0, f_grid[iv])));
    
    F[iv] = invPI * denom;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq == 0)
        d = - F[iv] * denom;
      
      if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
      {
        // Temperature derivative only depends on how pressure shift and broadening change
        dF(iq, iv) = d * Complex(dG0_dT, dL0_dT + ddF0_dT);
      }
      else if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        // Frequency scale 1 to -1 linearly
        dF(iq, iv) = d * Complex(0.0, -1.0);
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
        {
          // Line center scales 1 to 1 linearly
          dF(iq, iv) = d * Complex(0.0, 1.0);
        }
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR)
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
        {
          // Pressure broadening will be dealt with in another function, though the partial derivative
          dF(iq, iv) = d * Complex(0.0, -1.0);
        }
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        // Magnetic magnitude changes like line center in part
        // FIXME: Add magnetic components here
        dF(iq, iv) = d * Complex(0.0, zeeman_df);
      }
    }
  }
}


/*!
 * Sets the line shape to Hartmann-Tran line shape. Normalization is unity.
 * 
 * Note:  We are not experienced with this line shape and cannot tell what
 * parameters depends on what input.  It is therefore likely that this function
 * will have to be adapted in the future
 * 
 * Takes only one asymptotic limit into account...
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param G0 Speed-independent pressure broadening term
 * \param L0 Speed-independent pressure shift term plus second order line mixing shift
 * \param G2 Speed-dependent pressure broadening term
 * \param L2 Speed-dependent pressure shift term
 * \param eta Correlation between velocity and rotational changes due to collisions
 * \param FVC Velocity--changing collisional parameter
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
 * \param dG2_dT Temperature derivative of G2
 * \param dL2_dT Temperature derivative of L2
 * \param deta_dT Temperature derivative of eta
 * \param dFVC_dT Temperature derivative of FVC
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_htp(ComplexVectorView F, // Sets the full complex line shape without line mixing
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
                            const ArrayOfRetrievalQuantity& derivatives_data,
                            const ArrayOfIndex& derivatives_data_position,
                            const QuantumIdentifier& quantum_identity,
                            const Numeric& dGD_div_F0_dT,
                            const Numeric& dG0_dT,
                            const Numeric& dL0_dT,
                            const Numeric& dG2_dT,
                            const Numeric& dL2_dT,
                            const Numeric& deta_dT,
                            const Numeric& dFVC_dT)
{
  // Size of the problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  static const Complex i(0.0, 1.0), one_plus_one_i(1.0, 1.0);
  
  // Main lineshape parameters
  Complex A, Zp, Zm, Zp2, Zm2, wiZp, wiZm, X, sqrtXY, G, invG;
  
  // Derivatives parameters
  Complex dA, dZp, dZm, dwiZp, dwiZm, dG;
  
  // Test parameters for deciding how to compute certain asymptotes
  Numeric absX, ratioXY;
  
  // The line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude;
  
  // Doppler terms
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // Index of derivative for first occurrence of similarly natured partial derivatives
  const Index first_frequency = get_first_frequency_index(derivatives_data, derivatives_data_position); 
  
  // Pressure broadening terms
  const Complex C0 = G0 + i*L0;
  const Complex C2 = G2 + i*L2;
  
  // Correlation term
  const Numeric one_minus_eta = 1.0 - eta;
  
  // Pressure terms adjusted by collisions and correlations
  const Complex C0_m1p5_C2 = (C0 - 1.5 * C2);
  const Complex C0t = one_minus_eta * C0_m1p5_C2 + FVC;
  const Complex C2t = one_minus_eta * C2;
  
  // Limits and computational tracks
  const bool eta_zero_limit = eta == 0;
  const bool eta_one_limit = eta == 1;
  const bool C2t_zero_limit = C2t == Complex(0.0, 0.0) or eta_one_limit;
  bool ratioXY_low_limit, ratioXY_high_limit;
  
  // Practical term
  const Complex invC2t = C2t_zero_limit ? -1.0 : (1.0 / C2t);
  
  // Relative pressure broadening terms to be used in the shape function
  const Complex sqrtY = 0.5 * invC2t * GD;
  const Complex Y = sqrtY * sqrtY;
  const Numeric invAbsY = 1.0 / abs(Y);
  
  // Temperature derivatives precomputed
  Complex dC0_dT, dC2_dT, dC0t_dT, dC2t_dT, dC0_m1p5_C2_dT, dY_dT;
  if(do_temperature_jacobian(derivatives_data)) {
    dC0_dT = dG0_dT + i*dL0_dT;
    dC2_dT = dG2_dT + i*dL2_dT;
    dC0t_dT = one_minus_eta * (dC0_dT - 1.5 * dC2_dT) - deta_dT * C0_m1p5_C2 + dFVC_dT;
    dC2t_dT = one_minus_eta * dC2_dT - deta_dT * C2;
    dC0_m1p5_C2_dT = dC0_dT - 1.5 * dC2_dT;
    dY_dT = C2t_zero_limit ? -1.0 : (GD / 2.0*invC2t*invC2t * (dGD_dT - GD * invC2t * dC2t_dT));
  }
  
  // Scale factor (normalizes to PI)
  const Numeric fac = sqrtPI * invGD;
  
  #pragma omp simd
  for(Index iv = 0; iv < nf; iv++) {
    // Relative frequency
    X = (C0t + i*(F0 - f_grid[iv])) * invC2t;
    absX = abs(X);
    sqrtXY = sqrt(X + Y);
    
    // Limit tester.  Not valid if C2t_zero_limit is true;
    ratioXY = absX * invAbsY;
    
    // Limits are hard-coded from original paper by Tran etal
    ratioXY_high_limit = ratioXY >= 10e15;
    ratioXY_low_limit  = ratioXY <= 3.0*10e-8;
    
    // Compute Zm
    if(C2t_zero_limit or ratioXY_low_limit)
      Zm = (C0t + i*(F0 - f_grid[iv])) * invGD;
    else if(ratioXY_high_limit)
      Zm = sqrt(X);
    else
      Zm = sqrtXY - sqrtY;
    
    // Compute Zp
    if(C2t_zero_limit)
    {/* Do nothing since Zp will be infinity large */}
    else if(ratioXY_high_limit)
      Zp = sqrtXY;
    else
      Zp = sqrtXY + sqrtY;
      
    //Helpers
    Zm2 = Zm * Zm;
    if(not C2t_zero_limit)
      Zp2 = Zp * Zp;
    
    // Compute W-minus
    wiZm = Faddeeva::w(i*Zm);
    
    // Compute W-plus
    if(C2t_zero_limit)
      wiZp = Complex(0.0, 0.0);  // sqrtInvPI / Zp**2
    else
      wiZp = Faddeeva::w(i*Zp);
    
    // Compute A
    if(ratioXY_high_limit)
      A = fac * (sqrtInvPI - Zm * wiZm);
    else
      A = fac * (wiZm - wiZp);
    
    // Compute G
    if(eta_zero_limit)
      G = 1.0 - FVC * A;
    else if(C2t_zero_limit)
      G = 1.0 - (FVC - eta * C0_m1p5_C2) * A + C2 * fac * ((1.0 - Zm2)*wiZm + Zm * sqrtInvPI);
    else if (ratioXY_high_limit)
      G = 1.0 - (FVC - eta * C0_m1p5_C2) * A + eta / one_minus_eta * (-1.0 + 2.0 * sqrtPI * (1.0-X-2.0*Y) * (sqrtInvPI - Zm * wiZm) + 2.0*sqrtPI*Zp*wiZp);
    else
      G = 1.0 - (FVC - eta * C0_m1p5_C2) * A + eta / one_minus_eta * (-1.0 + sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*wiZm - (1.0-Zp2)*wiZp));
   
    // Compute denominator
    invG = 1.0 / G;
    
    // Compute line shape
    F[iv] = A * invG * invPI;
    
    for(Index iq = 0; iq < nppd; iq++) {
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]])) {
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq) {
          // Compute dZm
          if(C2t_zero_limit or ratioXY_low_limit)
            dZm = -invGD * i;
          else if(ratioXY_high_limit)
            dZm = -0.5 * invC2t * i / Zm;
          else
            dZm = -0.5 * invC2t * i / sqrtXY;
          
          // Compute dZp
          if(not C2t_zero_limit)
            dZp = 0.5 * invC2t * i / sqrtXY;
          
          // Compute dW-minus 
          dwiZm = 2.0 * (Zm * wiZm - sqrtInvPI) * dZm; // FIXME: test this for different asymptotes
          
          // Compute dW-plus
          if(C2t_zero_limit)
            dwiZp = Complex(0.0, 0.0);  // sqrtInvPI / Zp
          else
            dwiZp = 2.0 * (Zp * wiZp - sqrtInvPI) * dZp; // FIXME: test this for different asymptotes
          
          // Compute dA
          if(ratioXY_high_limit)
            dA = fac * (- dZm * wiZm - Zm * dwiZm);
          else
            dA = fac * (dwiZm - dwiZp);
          
          // Compute dG
          if(eta_zero_limit)
            dG = - FVC * dA;
          else if(C2t_zero_limit)
            dG = - (FVC - eta * C0_m1p5_C2) * dA + C2 * fac * ((1.0 - Zm2)*dwiZm + (-2.0*Zm)*wiZm + dZm * sqrtInvPI);
          else if (ratioXY_high_limit)
            dG = - (FVC - eta*C0_m1p5_C2) * dA + eta / one_minus_eta * (2.0 * sqrtPI * (-invC2t * i) * (sqrtInvPI - Zm * wiZm) + 2.0 * sqrtPI * (1.0-X-2.0*Y) * (- dZm * wiZm - Zm * dwiZm)+ 2.0*sqrtPI*(dZp*wiZp+Zp*dwiZp));
          else
            dG = - (FVC - eta*C0_m1p5_C2) * dA + eta / one_minus_eta * (sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*dwiZm + (-2.0*Zm)*wiZm - (1.0-Zp2)*dwiZp)- (-2.0*Zp)*wiZp);
          
          dF(iq, iv) = invG * (invPI * dA - F[iv] * dG); 
        }
        else  // copy for repeated occurrences
          dF(iq, iv) = dF(first_frequency, iv); 
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature) {
        // Compute dZm
        if(C2t_zero_limit or ratioXY_low_limit)
          dZm = (dC0t_dT - Zm * dGD_dT) * invGD;
        else if(ratioXY_high_limit)
          dZm = 0.5 / Zm  * (dC0t_dT - X * dC2t_dT) * invC2t;
        else
          dZm = 0.5 / sqrtXY * ((dC0t_dT - X * dC2t_dT) * invC2t + dY_dT) - 0.5 / sqrtY * dY_dT;
        
        // Compute Zp
        if(C2t_zero_limit)
        {/* Do nothing since Zp will be infinity large */}
        else if(ratioXY_high_limit)
          dZp = 0.5 / sqrtXY * ((dC0t_dT - X * dC2t_dT) * invC2t + dY_dT);
        else
          dZp = 0.5 / sqrtXY * ((dC0t_dT - X * dC2t_dT) * invC2t + dY_dT) + 0.5 / sqrtY * dY_dT;
          
        // Compute dW-minus 
        dwiZm = 2.0 * (Zm * wiZm - sqrtInvPI) * dZm; // FIXME: test this for different asymptotes
        
        // Compute dW-plus
        if(C2t_zero_limit)
          dwiZp = Complex(0.0, 0.0);  // sqrtInvPI / Zp
        else
          dwiZp = 2.0 * (Zp * wiZp - sqrtInvPI) * dZp; // FIXME: test this for different asymptotes
          
        // Compute dA
        if(ratioXY_high_limit)
          dA = fac * (- dZm * wiZm - Zm * dwiZm) - A * invGD * dGD_dT;
        else
          dA = fac * (dwiZm - dwiZp) - A * invGD * dGD_dT;
        
        // Compute dG
        if(eta_zero_limit)
          dG = - dFVC_dT * A - FVC * dA;
        else if(C2t_zero_limit)
          dG = - (dFVC_dT - deta_dT * C0_m1p5_C2 - eta * dC0_m1p5_C2_dT) * A - (FVC - eta*C0_m1p5_C2) * dA 
          + dC2_dT * fac * ((1.0 - Zm2)*wiZm + Zm * sqrtInvPI)
          + C2 * fac * ((1.0 - Zm2)*dwiZm - 2.0 * Zm * dwiZm + dZm * sqrtInvPI)
          - C2 * fac * ((1.0 - Zm2)*wiZm + Zm * sqrtInvPI) * invGD * dGD_dT;
        else if (ratioXY_high_limit)
          dG = - (dFVC_dT - deta_dT * C0_m1p5_C2 - eta * dC0_m1p5_C2_dT) * A - (FVC - eta*C0_m1p5_C2) * dA 
          +
          deta_dT / one_minus_eta * (-1.0 + 2.0 * sqrtPI * (1.0-X-2.0*Y) * (sqrtInvPI - Zm * wiZm) + 2.0*sqrtPI*Zp*wiZp)
          -
          deta_dT * eta / (one_minus_eta*one_minus_eta) * (-1.0 + 2.0 * sqrtPI * (1.0-X-2.0*Y) * (sqrtInvPI - Zm * wiZm) + 2.0*sqrtPI*Zp*wiZp)
          +
          eta / one_minus_eta * (
            2.0 * sqrtPI * (-((dC0t_dT - X * dC2t_dT) * invC2t)-2.0*dY_dT) * (sqrtInvPI - Zm * wiZm)  +
            2.0 * sqrtPI * (1.0-X-2.0*Y) * (- dZm * wiZm - Zm * dwiZm) 
          + 2.0 * sqrtPI*(dZp*wiZp + Zp*dwiZp))
          ;
        else
          dG = - (dFVC_dT - deta_dT * C0_m1p5_C2 - eta * dC0_m1p5_C2_dT) * A - (FVC - eta*C0_m1p5_C2) * dA 
          + 
          deta_dT / one_minus_eta * (-1.0 + sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*wiZm - (1.0-Zp2)*wiZp))
          -
          deta_dT * eta / (one_minus_eta*one_minus_eta) * (-1.0 + sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*wiZm - (1.0-Zp2)*wiZp))
          +
          eta / one_minus_eta * 
          (
            -0.5 * sqrtPI/(2.0*sqrtY)/sqrtY * dY_dT * ((1.0-Zm2)*wiZm - (1.0-Zp2)*wiZp)
            +
            (sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*dwiZm - (1.0-Zp2)*dwiZp)+ (-2.0*Zm)*wiZm - (-2.0*Zp)*wiZp)
          );
          
        dF(iq, iv) = invG * (invPI * dA - F[iv] * dG); 
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter) {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
          // Compute dZm
          if(C2t_zero_limit or ratioXY_low_limit)
            dZm = (i - Zm) * invGD;
          else if(ratioXY_high_limit)
            dZm = 0.5/Zm * i * invC2t;
          else
            dZm = 0.5/sqrtXY * i * invC2t;
          
          // Compute dZp
          if(C2t_zero_limit)
          {/* Do nothing since Zp will be infinity large */}
          else if(ratioXY_high_limit or ratioXY_low_limit)
            dZp = 0.5/sqrtXY * i * invC2t;
          else
            dZp = dZm;
          
          // Compute dW-minus 
          dwiZm = 2.0 * (Zm * wiZm - sqrtInvPI) * dZm; // FIXME: test this for different asymptotes
          
          // Compute dW-plus
          if(C2t_zero_limit)
            dwiZp = Complex(0.0, 0.0);  // sqrtInvPI / Zp
          else
            dwiZp = 2.0 * (Zp * wiZp - sqrtInvPI) * dZp; // FIXME: test this for different asymptotes
            
          // Compute dA
          if(ratioXY_high_limit)
            dA = fac * (- dZm * wiZm - Zm * dwiZm) - A * invGD;
          else
            dA = fac * (dwiZm - dwiZp) - A * invGD;
          
          // Compute G
          if(eta_zero_limit)
            dG = - FVC * dA;
          else if(C2t_zero_limit)
            dG = - (FVC - eta * C0_m1p5_C2) * dA + C2 * fac * (((1.0 - Zm2)*dwiZm + dZm * sqrtInvPI - 2.0*Zm*dZm*wiZm) - ((1.0 - Zm2)*wiZm + Zm * sqrtInvPI) * invGD);
          else if (ratioXY_high_limit)
            dG = - (FVC - eta*C0_m1p5_C2) * dA + eta / one_minus_eta * (2.0 * sqrtPI * (-i*invC2t * (sqrtInvPI - Zm * wiZm) + (1.0-X-2.0*Y) * (- dZm * wiZm - Zm * dwiZm)) +             2.0*sqrtPI*(dZp*wiZp + Zp*dwiZp));
          else
            dG = - (FVC - eta*C0_m1p5_C2) * dA + eta / one_minus_eta * (sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*dwiZm - (1.0-Zp2)*dwiZp - 2.0 * Zm * dZm * dwiZm + 2.0 * Zp * dZp * wiZp));
          
          dF(iq, iv) = invG * (invPI * dA - F[iv] * dG); 
        }
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR)  {
        // NOTE:  These are first order Voigt-like.  
        // The variables that are not Voigt-like must be dealt with separately
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
          // Compute dZm
          if(C2t_zero_limit or ratioXY_low_limit)
            dZm = one_minus_eta * invGD;
          else if(ratioXY_high_limit)
            dZm = 0.5 / Zm * one_minus_eta * invC2t;
          else
            dZm = 0.5 * sqrtXY * one_minus_eta * invC2t;
          
          // Compute dZp
          if(C2t_zero_limit)
          {/* Do nothing since Zp will be infinity large */}
          else if(ratioXY_high_limit or ratioXY_low_limit)
            dZp = 0.5 * sqrtXY * one_minus_eta * invC2t;
          else
            dZp = dZm;
          
          // Compute dW-minus 
          dwiZm = 2.0 * (Zm * wiZm - sqrtInvPI) * dZm; // FIXME: test this for different asymptotes
          
          // Compute dW-plus
          if(C2t_zero_limit)
            dwiZp = Complex(0.0, 0.0);  // sqrtInvPI / Zp
          else
            dwiZp = 2.0 * (Zp * wiZp - sqrtInvPI) * dZp; // FIXME: test this for different asymptotes
            
          // Compute dA
          if(ratioXY_high_limit)
            dA = fac * (- dZm * wiZm - Zm * dwiZm);
          else
            dA = fac * (dwiZm - dwiZp);
          
          // Compute dG
          if(eta_zero_limit)
            dG = - FVC * dA;
          else if(C2t_zero_limit)
            dG = - (FVC - eta * C0_m1p5_C2) * dA + eta * A + C2 * fac * ((1.0 - Zm2)*dwiZm + dZm * sqrtInvPI - 2.0*Zm*dZm*wiZm);
          else if (ratioXY_high_limit)
            dG = - (FVC - eta*C0_m1p5_C2) * dA + eta * A + eta / one_minus_eta * (2.0 * sqrtPI * (1.0-X-2.0*Y) * (- dZm * wiZm - Zm * dwiZm) + 2.0 * sqrtPI * (-one_minus_eta * invC2t) * (sqrtInvPI - Zm * wiZm) + 2.0*sqrtPI*(dZp*wiZp+Zp*dwiZp));
          else
            dG = - (FVC - eta*C0_m1p5_C2) * dA + eta * A + eta / one_minus_eta * (sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*dwiZm - (1.0-Zp2)*dwiZp - 2.0*Zm*dZm*wiZm + 2.0*Zp*dZp*wiZp));
          
          dF(iq, iv) = invG * (invPI * dA - F[iv] * dG) * Complex(0.0, -1.0); 
        }
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[iq]])) {
        // Compute dZm
        if(C2t_zero_limit or ratioXY_low_limit)
          dZm = i * invGD;
        else if(ratioXY_high_limit)
          dZm = 0.5/Zm * i * invC2t;
        else
          dZm = 0.5/sqrtXY * i * invC2t;
        
        // Compute dZp
        if(C2t_zero_limit)
        {/* Do nothing since Zp will be infinity large */}
        else if(ratioXY_high_limit or ratioXY_low_limit)
          dZp = 0.5/sqrtXY * i * invC2t;
        else
          dZp = dZm;
        
        // Compute dW-minus 
        dwiZm = 2.0 * (Zm * wiZm - sqrtInvPI) * dZm; // FIXME: test this for different asymptotes
        
        // Compute dW-plus
        if(C2t_zero_limit)
          dwiZp = Complex(0.0, 0.0);  // sqrtInvPI / Zp
        else
          dwiZp = 2.0 * (Zp * wiZp - sqrtInvPI) * dZp; // FIXME: test this for different asymptotes
          
        // Compute dA
        if(ratioXY_high_limit)
          dA = fac * (- dZm * wiZm - Zm * dwiZm);
        else
          dA = fac * (dwiZm - dwiZp);
        
        // Compute G
        if(eta_zero_limit)
          dG = - FVC * dA;
        else if(C2t_zero_limit)
          dG = - (FVC - eta * C0_m1p5_C2) * dA + C2 * fac * ((1.0 - Zm2)*dwiZm + dZm * sqrtInvPI - 2.0*Zm*dZm*wiZm);
        else if (ratioXY_high_limit)
          dG = - (FVC - eta*C0_m1p5_C2) * dA + eta / one_minus_eta * (2.0 * sqrtPI * (-i*invC2t * (sqrtInvPI - Zm * wiZm) + (1.0-X-2.0*Y) * (- dZm * wiZm - Zm * dwiZm)) + 2.0*sqrtPI*(dZp*wiZp + Zp*dwiZp));
        else
          dG = - (FVC - eta*C0_m1p5_C2) * dA + eta / one_minus_eta * (sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*dwiZm - (1.0-Zp2)*dwiZp - 2.0 * Zm * dZm * dwiZm + 2.0 * Zp * dZp * wiZp));
        
        dF(iq, iv) = invG * (invPI * dA - F[iv] * dG); 
      }
    }
  }
}


/*!
 * Sets the line shape to Voigt line shape. Normalization is unity.
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param G0 Speed-independent pressure broadening term
 * \param L0 Speed-independent pressure shift term
 * \param dF0 Second order line mixing shift
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
 * \param dF0_dT Temperature derivative of dF0_dT
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_voigt(ComplexVectorView F, 
                              ComplexMatrixView dF, 
                              ConstVectorView f_grid, 
                              const Numeric& zeeman_df, 
                              const Numeric& magnetic_magnitude,
                              const Numeric& F0_noshift, 
                              const Numeric& GD_div_F0,
                              const Numeric& G0, 
                              const Numeric& L0,
                              const Numeric& dF0,
                              const ArrayOfRetrievalQuantity& derivatives_data,
                              const ArrayOfIndex& derivatives_data_position,
                              const QuantumIdentifier& quantum_identity,
                              const Numeric& dGD_div_F0_dT,
                              const Numeric& dG0_dT,
                              const Numeric& dL0_dT,
                              const Numeric& dF0_dT)
{
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  // For calculations
  Complex w, z, dw_over_dz, dz;
  
  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + L0 + dF0;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0 + GD_div_F0 * (-dL0_dT - dF0_dT);
  
  // constant normalization factor for Voigt
  const Numeric fac = sqrtInvPI * invGD;
  
  // Ratio of the Lorentz halfwidth to the Doppler halfwidth
  const Complex z0 = Complex(-F0, G0) * invGD;
  
  // frequency in units of Doppler
  #pragma omp simd
  for (Index iv=0; iv<nf; iv++) {
    z = z0 + f_grid[iv] * invGD;
    w = Faddeeva::w(z);
    F[iv] = fac * w;
    
    for(Index iq = 0; iq < nppd; iq++) {
      if(iq==0)  // Standard basic form for all transitions
        dw_over_dz = 2.0 * fac *  (Complex(0, sqrtInvPI) - z * w);
      
      // switch-like statement for all relevant partials in the xsec calculations
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]]))
        dF(iq, iv) = dw_over_dz * invGD;
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
        dF(iq, iv) = -F[iv] * dGD_dT * invGD + dw_over_dz * ((Complex(-dL0_dT - dF0_dT, dG0_dT) - z * dGD_dT) * invGD);
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]])) {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
          dF(iq, iv) = -F[iv]/F0 + dw_over_dz * (-invGD - z/F0);
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR) {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
          dF(iq, iv) = dw_over_dz * invGD;
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[iq]]))
        dF(iq, iv) = dw_over_dz * (- zeeman_df * invGD); //* dz;
    }
  }
}


/*!
 * Sets the line shape to Doppler line shape. Normalization is unity.
 * 
 * Note:  Uses the Voigt function with special parameters to get imaginary
 * part of spectra.  This should be fixed to gain speed...
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_doppler(ComplexVectorView F, // Sets the full complex line shape without line mixing
                                ComplexMatrixView dF,
                                ConstVectorView f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& GD_div_F0,
                                const ArrayOfRetrievalQuantity& derivatives_data,
                                const ArrayOfIndex& derivatives_data_position,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dGD_div_F0_dT)
{
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // Computational variables
  //Numeric x, w_im, w_re;
  Complex w, dw_over_dx;
  Complex x;
  
  const Numeric fac = invGD * sqrtInvPI;
  
  // Computational speed-up
  const Index first_frequency = get_first_frequency_index(derivatives_data, derivatives_data_position);
  
  #pragma omp simd
  for(Index iv = 0; iv < nf; iv++)
  {
    x = (f_grid[iv] - F0) * invGD;
    w = Faddeeva::w(x);
    
    F[iv] = fac * w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq == 0)
        dw_over_dx = 2.0 * fac * (Complex(0.0, sqrtInvPI) - x * w);
      
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]] ))
      {
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          dF(iq, iv) = dw_over_dx * invGD;
        }
        else  // copy for repeated occurrences
        {
          dF(iq, iv) = dF(first_frequency, iv);
        }
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
      {
        dF(iq, iv) = F[iv] * dGD_dT + x * dGD_dT * dw_over_dx;
        dF(iq, iv) *= -invGD ;
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
        {
          dF(iq, iv) = -F[iv] * GD_div_F0;
          dF(iq, iv) += dw_over_dx * (-x * GD_div_F0 - 1.0);
          dF(iq, iv) *= invGD;
        }
      }
    }
  }
}


/*!
 * Sets the line shape using eigenvalue decomposition for line mixing. Normalization is unity.
 * 
 * Note:  This is still under construction (until there is a publication, then this comment should be changed)
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param eigenvalue_no_shift Output from eigenvalue decomposition of band relaxation at some frequency
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param L0 Speed-independent pressure shift term
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dL0_dT Temperature derivative of L0
 * 
 */
void Linefunctions::set_voigt_from_full_linemixing(ComplexVectorView F, 
                                                   ComplexMatrixView dF,
                                                   ConstVectorView f_grid,
                                                   const Complex& eigenvalue_no_shift,
                                                   const Numeric& GD_div_F0,
                                                   const Numeric& L0,
                                                   const ArrayOfRetrievalQuantity& derivatives_data,
                                                   const ArrayOfIndex& derivatives_data_position,
                                                   const QuantumIdentifier& quantum_identity,
                                                   const Numeric& dGD_div_F0_dT,
                                                   const Complex& deigenvalue_dT,
                                                   const Numeric& dL0_dT)
{
  // For calculations
  Numeric dx;
  Complex w, z, dw_over_dz, dz;
  
  const Index nf = f_grid.nelem(), nppd = derivatives_data_position.nelem();
  
  // Doppler broadening and line center
  const Complex eigenvalue = eigenvalue_no_shift + L0;
  const Numeric F0 = eigenvalue.real() + L0;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // constant normalization factor for voigt
  const Numeric fac = sqrtInvPI * invGD;
  
  // Ratio of the Lorentz halfwidth to the Doppler halfwidth
  const Complex z0 = -eigenvalue * invGD;
  
  // frequency in units of Doppler
  #pragma omp simd
  for (Index iv=0; iv<nf; iv++)
  {
    dx = f_grid[iv] * invGD;
    z = z0 + dx;
    w = Faddeeva::w(z);
    
    F[iv] = fac * w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq==0)
        dw_over_dz = 2.0 * (z * w - sqrtInvPI);
      
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]] ))
        dF(iq, iv) = fac * dw_over_dz * invGD;
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
      {
        dz = (deigenvalue_dT - dL0_dT) - z * dGD_dT;
        
        dF(iq, iv) = -F[iv] * dGD_dT;
        dF(iq, iv) += fac * dw_over_dz * dz;
        dF(iq, iv) *= invGD;
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_parameter(derivatives_data[derivatives_data_position[iq]])) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF(iq, iv) = -F[iv] * GD_div_F0;
          dF(iq, iv) += dw_over_dz * dz;
          dF(iq, iv) *= fac * invGD;
        }
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR)
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
          dF(iq, iv) = fac * dw_over_dz * ( Complex(-1.0, 1.0) * invGD);
      }
    }
    
  }
}


/*!
 * Applies line mixing scaling to already set lineshape
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param Y First-order line-mixing coefficient
 * \param G Second order line-mixing coefficient
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dY_dT Temperature derivative of Y
 * \param dG_dT Temperature derivative of G
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linemixing_scaling(ComplexVectorView F,
                                             ComplexMatrixView dF,
                                             const Numeric& Y,
                                             const Numeric& G,
                                             const ArrayOfRetrievalQuantity& derivatives_data,
                                             const ArrayOfIndex& derivatives_data_position,
                                             const QuantumIdentifier& quantum_identity,
                                             const Numeric& dY_dT,
                                             const Numeric& dG_dT)
{
  const Index nf = F.nelem(), nppd = derivatives_data_position.nelem();
  
  const Complex LM = Complex(1.0 + G, -Y);
  const Complex dLM_dT = Complex(dG_dT, -dY_dT);
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
    {
      dF(iq, joker) *= LM;
      #pragma omp simd
      for(Index iv = 0; iv < nf; iv++)
        dF(iq, iv) += F[iv] * dLM_dT;
    }
    else if(is_line_mixing_line_strength_parameter(derivatives_data[derivatives_data_position[iq]]))
    {
       if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
         dF(iq, joker) = F;
    }
    else
      dF(iq, joker) *= LM;
  }
  
  F *= LM;
}

/*!
 * Applies Rosenkranz quadratic normalization to already set line shape
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param F0 Central frequency without any shifts
 * \param T Atmospheric temperature at level
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_rosenkranz_quadratic_scaling(ComplexVectorView F,
                                                       ComplexMatrixView dF,
                                                       ConstVectorView f_grid,
                                                       const Numeric& F0,
                                                       const Numeric& T,
                                                       const ArrayOfRetrievalQuantity& derivatives_data,
                                                       const ArrayOfIndex& derivatives_data_position,
                                                       const QuantumIdentifier& quantum_identity)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data_position.nelem();
  
  const Numeric invF0 = 1.0/F0;
  const Numeric mafac = (PLANCK_CONST) / (2.0 * BOLTZMAN_CONST * T) /
  sinh((PLANCK_CONST * F0) / (2.0 * BOLTZMAN_CONST * T)) * invF0;
  
  Numeric dmafac_dT_div_fun = 0;
  if(do_temperature_jacobian(derivatives_data))
  {
    dmafac_dT_div_fun = -(BOLTZMAN_CONST*T - F0*PLANCK_CONST/
    (2.0*tanh(F0*PLANCK_CONST/(2.0*BOLTZMAN_CONST*T))))/(BOLTZMAN_CONST*T*T);
  }
  
  Numeric fun;
  
  #pragma omp simd
  for (Index iv=0; iv < nf; iv++)
  {
    fun = mafac * (f_grid[iv] * f_grid[iv]);
    F[iv] *= fun;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      dF(iq, iv) *= fun;
      if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
      {
        dF(iq, iv) += dmafac_dT_div_fun * F[iv];
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
          const Numeric dmafac_dF0_div_fun = -invF0 - PLANCK_CONST/(2.0*BOLTZMAN_CONST*T*tanh(F0*PLANCK_CONST/(2.0*BOLTZMAN_CONST*T)));
          dF(iq, iv) += dmafac_dF0_div_fun * F[iv];
        }
      }
      else if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        dF(iq, iv) += (2.0 / f_grid[iv]) * F[iv];
      }
    }
  }
}


/*!
 * Applies Van Vleck and Huber normalization to already set line shape
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param F0 Central frequency without any shifts
 * \param T Atmospheric temperature at level
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_VVH_scaling(ComplexVectorView F,
                                      ComplexMatrixView dF,
                                      ConstVectorView f_grid,
                                      const Numeric& F0,
                                      const Numeric& T,
                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                      const ArrayOfIndex& derivatives_data_position,
                                      const QuantumIdentifier& quantum_identity)
{ 
  const Index nf = f_grid.nelem(), nppd = derivatives_data_position.nelem();
  
  // 2kT is constant for the loop
  const Numeric kT = 2.0 * BOLTZMAN_CONST * T;
  
  // denominator is constant for the loop
  const Numeric tanh_f0part = tanh(PLANCK_CONST * F0 / kT);
  const Numeric denom = F0 * tanh_f0part;
  
  #pragma omp simd
  for(Index iv=0; iv < nf; iv++)
  {
    const Numeric tanh_fpart = tanh( PLANCK_CONST * f_grid[iv] / kT );
    const Numeric fun = f_grid[iv] * tanh_fpart / denom;
    F[iv] *= fun;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      dF(iq, iv) *= fun;
      if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
      {
        dF(iq, iv) += (-PLANCK_CONST*(denom - F0/tanh_f0part - 
        f_grid[iv]*tanh_fpart + f_grid[iv]/tanh_fpart)/(kT*T)) * F[iv];
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
          const Numeric fac_df0 = (-1.0/F0 + PLANCK_CONST*tanh_f0part/(kT) - PLANCK_CONST/(kT*tanh_f0part)) * F0/F0;
          dF(iq, iv) += fac_df0 * F[iv];
        }
      }
      else if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        dF(iq, iv) += (1.0/f_grid[iv] -PLANCK_CONST*tanh_fpart/kT + PLANCK_CONST/(kT*tanh_fpart)) * F[iv];
      }
    }
  }
}


/*!
 * Applies Van Vleck and Weiskopf normalization to already set line shape
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param f_grid Frequency grid of computations
 * \param F0 Central frequency without any shifts
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_VVW_scaling(ComplexVectorView F,
                                      ComplexMatrixView dF,
                                      ConstVectorView f_grid,
                                      const Numeric& F0,
                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                      const ArrayOfIndex& derivatives_data_position,
                                      const QuantumIdentifier& quantum_identity)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data_position.nelem();
  
  // denominator is constant for the loop
  const Numeric invF0 = 1.0 / F0;
  const Numeric invF02 = invF0 * invF0;
  
  #pragma omp simd
  for(Index iv = 0; iv < nf; iv++)
  {
    // Set the factor
    const Numeric fac = f_grid[iv] * f_grid[iv] * invF02;
    
    // Set the line shape
    F[iv] *= fac;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      // The factor is applied to all partial derivatives
      dF(iq, iv) *= fac;
      
      // These partial derivatives are special
      if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
          dF(iq, iv) -= 2.0 * invF0 * F[iv] ;
      }
      else if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]]))
      {
        dF(iq, iv) += 2.0 / f_grid[iv] * F[iv];
      }
    }
  }
}


/*!
 * Applies linestrength to already set line shape
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param S0 Strength of line at reference temperature
 * \param isotopic_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param QT Partition function at atmospheric temperature of level
 * \param QT0 Partition function at reference temperature
 * \param K1 Boltzmann ratio between reference and atmospheric temperatures
 * \param K2 Stimulated emission ratio between reference and atmospheric temperatures
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dQT_dT Temperature derivative of QT
 * \param dK1_dT Temperature derivative of K1
 * \param dK2_dT Temperature derivative of K2
 * \param dK2_dF0 Central frequency derivative of K2
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linestrength_scaling(ComplexVectorView F,
                                               ComplexMatrixView dF,
                                               const Numeric& S0,
                                               const Numeric& isotopic_ratio,
                                               const Numeric& QT,
                                               const Numeric& QT0,
                                               const Numeric& K1,
                                               const Numeric& K2,
                                               const ArrayOfRetrievalQuantity& derivatives_data,
                                               const ArrayOfIndex& derivatives_data_position,
                                               const QuantumIdentifier& quantum_identity,
                                               const Numeric& dQT_dT,
                                               const Numeric& dK1_dT,
                                               const Numeric& dK2_dT,
                                               const Numeric& dK2_dF0)
{
  const Index nf = F.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  const Numeric invQT = 1.0/QT;
  const Numeric QT_ratio = QT0 * invQT;
  
  const Numeric dS_dS0 = isotopic_ratio * QT_ratio * K1 * K2;
  const Numeric S = S0 * dS_dS0;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
    {
      const Numeric dS_dT = S * (dK2_dT / K2 + dK1_dT / K1 - invQT * dQT_dT);
      
      dF(iq, joker) *= S;
      
      #pragma omp simd
      for(Index iv = 0; iv < nf; iv++)
        dF(iq, iv) += F[iv] * dS_dT;
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineStrength)
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        dF(iq, joker) = F;
        dF(iq, joker) *= dS_dS0;
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter)
    {
      
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        const Numeric dS_dF0 = S * dK2_dF0 / K2;
        
        dF(iq, joker) *= S;
        
        #pragma omp simd
        for(Index iv = 0; iv < nf; iv++)
          dF(iq, iv) += F[iv] * dS_dF0;
      }
    }
    else
    {
      dF(iq, joker) *= S;
    }
  }
  
  // Set lineshape at the end
  F *= S;
}


/*! Applies the line strength to the line shape using the complex line strength of line mixing.
 * Lineshape is already set.
 * 
 * Note:  This is still under construction (until there is a publication, then this comment should be changed)
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param F0 Central frequency without any shifts
 * \param T Atmospheric temperature at level
 * \param S_LM Complex linestrength
 * \param isotopic_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dS_LM_dT Temperature derivative of S_LM
 * 
 */
void Linefunctions::apply_linestrength_from_full_linemixing(ComplexVectorView F,
                                                            ComplexMatrixView dF,
                                                            const Numeric& F0,
                                                            const Numeric& T,
                                                            const Complex& S_LM,
                                                            const Numeric& isotopic_ratio,
                                                            const ArrayOfRetrievalQuantity& derivatives_data,
                                                            const ArrayOfIndex& derivatives_data_position,
                                                            const QuantumIdentifier& quantum_identity,
                                                            const Complex& dS_LM_dT)
{
  const Index nppd = derivatives_data_position.nelem();
  
  const Numeric invT = 1.0 / T, 
  F0_invT = F0 * invT,
  exp_factor = exp(C1 * F0_invT), 
  f0_factor = F0 * (1.0 - exp_factor); 
  
  const Complex S = S_LM * f0_factor * isotopic_ratio;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
    {
      Eigen::VectorXcd eig_dF = MapToEigen(dF(iq, joker));
      
      eig_dF *= S_LM;
      eig_dF += MapToEigen(F) * (dS_LM_dT * f0_factor + 
      S_LM * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio;
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineStrength)
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else
    {
      dF(iq, joker) *= S;
    }
  }
  
  F *= S;
}


/*! Applies the line strength to the line shape using the dipole information
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param F0 Central frequency without any shifts
 * \param T Atmospheric temperature at level
 * \param d0 Dipole of the absorption line
 * \param rho Density (of molecules at level
 * \param isotopic_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param drho_dT Temperature derivative of rho
 * 
 */
void Linefunctions::apply_dipole(ComplexVectorView F,
                                 ComplexMatrixView,
                                 const Numeric& F0,
                                 const Numeric& T,
                                 const Numeric& d0,
                                 const Numeric& rho,
                                 const Numeric& isotopic_ratio,
                                 const ArrayOfRetrievalQuantity&,
                                 const ArrayOfIndex& derivatives_data_position,
                                 const QuantumIdentifier&,
                                 const Numeric&)
{
  // Output is d0^2 * rho * F * isotopic_ratio * F0 * (1-e^(hF0/kT))
  
  const Index nppd = derivatives_data_position.nelem();
  
  const Numeric S = d0 * d0 * rho * isotopic_ratio, 
  invT = 1.0 / T, 
  F0_invT = F0 * invT,
  exp_factor = exp(C1 * F0_invT), 
  f0_factor = F0 * (1.0 - exp_factor);
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    throw std::runtime_error("Cannot support Jacobian from dipole calculations yet");
  }
  
  F *= S * f0_factor;
}


/*! Applies the line-by-line pressure broadening jacobian for the matching lines
 * 
 * \retval dF Lineshape derivative
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dgamma Derivatives in order as they appear that are related to pressure broadening coefficients
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_pressurebroadening_jacobian_scaling(ComplexMatrixView dF,
                                                              const ArrayOfRetrievalQuantity& derivatives_data,
                                                              const ArrayOfIndex& derivatives_data_position,
                                                              const QuantumIdentifier& quantum_identity,
                                                              const ComplexVector& dgamma)
{
  const Index nppd = derivatives_data_position.nelem(), ng = dgamma.nelem();
  
  Index ipd = 0;
  
  // Length of dgamma must be the same as total number of instances of pressure broadening jacobians
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(ipd == ng) break;
    
    if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR)
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        dF(iq, joker) *= dgamma[ipd];
        ++ipd;
      }
    }
  }
}


/*! Applies the line-by-line pressure broadening jacobian for the matching lines
 * 
 * \retval dF Lineshape derivative
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dlm Derivatives in order as they appear that are related to line mixing
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linemixing_jacobian_scaling(ComplexMatrixView dF,
                                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                                      const ArrayOfIndex& derivatives_data_position,
                                                      const QuantumIdentifier& quantum_identity,
                                                      const ComplexVector& dlm)
{
  const Index nppd = derivatives_data_position.nelem(), ng = dlm.nelem();
  
  Index ipd = 0;
  
  // Length of dlm must be the same as total number of instances of line mixing jacobians
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(ipd == ng) break;
    
    if(is_line_mixing_parameter(derivatives_data[derivatives_data_position[iq]]))
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        dF(iq, joker) *= dlm[ipd];
        ++ipd;
      }
    }
  }
}


/*! Returns the frequency-independent part of the Doppler broadening
 * 
 * \param T Atmospheric temperature at level
 * \param mass Mass of molecule under consideration
 * 
 */
Numeric Linefunctions::DopplerConstant(const Numeric& T, const Numeric& mass)
{
  return doppler_const * sqrt(T / mass);
}


/*! Returns the temperature derivative of the frequency-independent part of the Doppler broadening
 * 
 * \param T Atmospheric temperature at level
 * \param mass Mass of molecule under consideration
 * 
 */
Numeric Linefunctions::dDopplerConstant_dT(const Numeric& T, const Numeric& mass)
{
  return doppler_const * 0.5 * sqrt(1.0 / mass / T);
}


/*!
 * Applies non-lte linestrength to already set line shape
 * 
 * Only works for scaling-problems (i.e., when total number density distribution is already known)
 * 
 * \retval F Lineshape (absorption)
 * \retval dF Lineshape derivative (absorption)
 * \retval N Non-lte lineshape (source)
 * \retval dN Non-lte lineshape derivative (source)
 * 
 * \param K3 Ratio of absorption due to non-lte effects
 * \param K4 Ratio of emission due to non-lte effects
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dK3_dT Temperature derivative of K3
 * \param dK4_dT Temperature derivative of K4
 * \param dK3_dF0 Central frequency derivative of K3
 * \param dK3_dTl Lower energy state equivalent non-lte temperature derivative of K3
 * \param dK3_dTu Upper energy state equivalent non-lte temperature derivative of K3
 * \param dK4_dTu Upper energy state equivalent non-lte temperature derivative of K4
 * \param df_range Frequency range to use inside dF and dN
 * 
 */
void Linefunctions::set_nonlte_source_and_apply_absorption_scaling(ComplexVectorView F, 
                                                                   ComplexMatrixView dF,
                                                                   ComplexVectorView N,
                                                                   ComplexMatrixView dN,
                                                                   const Numeric& K3, 
                                                                   const Numeric& K4,
                                                                   const ArrayOfRetrievalQuantity& derivatives_data,
                                                                   const ArrayOfIndex& derivatives_data_position,
                                                                   const QuantumIdentifier& quantum_identity,
                                                                   const Numeric& dK3_dT, 
                                                                   const Numeric& dK4_dT,
                                                                   const Numeric& dK3_dF0, 
                                                                   const Numeric& dK3_dTl, 
                                                                   const Numeric& dK3_dTu, 
                                                                   const Numeric& dK4_dTu)
{
  const Index nppd = derivatives_data_position.nelem(), nf = F.nelem();
  
  const Numeric scaled_ratio = K4/K3 - 1.0;
  
  // Set the non-lte source factors
  N = F;
  N *= scaled_ratio;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    dN(iq, joker) = dF(iq, joker);
    dN(iq, joker) *= scaled_ratio;
    dF(iq, joker) *= K3;
    
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
    {
      const Numeric dscaled_ratio_dT = (dK4_dT - dK3_dT / K3) / K3;
      
      #pragma omp simd
      for(Index iv = 0; iv < nf; iv++)
      {
        dF(iq, iv) += F[iv] * dK3_dT;
        dN(iq, iv) += F[iv] * dscaled_ratio_dT;
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        const Numeric dscaled_ratio_dF0 = - dK3_dF0 / K3 / K3;
        
        #pragma omp simd
        for(Index iv = 0; iv < nf; iv++)
        {
          dF(iq, iv) += F[iv] * dK3_dF0;
          dN(iq, iv) += F[iv] * dscaled_ratio_dF0;
        }
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::NLTE)
    {
      if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().Species() not_eq quantum_identity.Species() or
        derivatives_data[derivatives_data_position[iq]].QuantumIdentity().Isotopologue() not_eq quantum_identity.Isotopologue())
        continue;  // Wrong species or wrong isotopologue
      
      if(quantum_identity.QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX] >
        derivatives_data[derivatives_data_position[iq]].QuantumIdentity().QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL] or
        derivatives_data[derivatives_data_position[iq]].QuantumIdentity().Type() == QuantumIdentifier::ALL)
      {
        const Numeric dscaled_ratio_dTl = - dK3_dTl / K3 / K3;
        
        #pragma omp simd
        for(Index iv = 0; iv < nf; iv++)
        {
          dF(iq, iv) += F[iv] * dK3_dTl;
          dN(iq, iv) += F[iv] * dscaled_ratio_dTl;
        }
      }
      
      if(quantum_identity.QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX] >
        derivatives_data[derivatives_data_position[iq]].QuantumIdentity().QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL] or
        derivatives_data[derivatives_data_position[iq]].QuantumIdentity().Type() == QuantumIdentifier::ALL)
      {
        const Numeric dscaled_ratio_dTu = (dK4_dTu - dK3_dTu / K3) / K3;
        
        #pragma omp simd
        for(Index iv = 0; iv < nf; iv++)
        {
          dF(iq, iv) += F[iv] * dK3_dTu;
          dN(iq, iv) += F[iv] * dscaled_ratio_dTu;
        }
      }
    }
  }
  
  // Finish by scaling F to the true value
  F *= K3;
}


/*!
 * Combination function using standard setup to compute line strength and lineshape of a single line.
 * Computes in order the lineshape, the linemirroring, the linenormalization, the linemixing, the 
 * linestrength, the non-lte, and the cutoff frequency.
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * \retval N Non-lte lineshape
 * \retval dN Non-lte lineshape derivative
 * \retval this_xsec_range Range indicating which frequency grids have been altered in F, dF, N, and dN
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param line line-record containing most line parameters
 * \param volume_mixing_ratio_of_all_species As name suggests
 * \param nlte_distribution As name suggests
 * \param pressure As name suggests
 * \param temperature As name suggests
 * \param doppler_constant Frequency-independent part of the Doppler broadening
 * \param partial_pressure Pressure of species that line belongs to at this level
 * \param isotopologue_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param ddoppler_constant_dT Temperature derivative of doppler_constant
 * \param pressure_limit_for_linemixing As WSV lm_p_lim
 * \param partition_function_at_temperature As name suggests
 * \param dpartition_function_at_temperature_dT Temeperature derivative of partition_function_at_temperature
 * \param partition_function_at_line_temperature As name suggests
 * \param broad_spec_locations Locations of broadening species using all-planetary broadening inside volume_mixing_ratio_of_all_species
 * \param this_species_location_in_tags Location of species of line in volume_mixing_ratio_of_all_species
 * \param water_index_location_in_tags Location of water in volume_mixing_ratio_of_all_species
 * \param verbosity Verbosity level
 * \param cutoff_call Flag to ignore some functions inside if this call is from the cutoff-computations
 * 
 */
void Linefunctions::set_cross_section_for_single_line(ComplexVectorView F_full,
                                                      ComplexMatrixView dF_full,
                                                      ComplexVectorView N_full, 
                                                      ComplexMatrixView dN_full,
                                                      Range& this_f_range,
                                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                                      const ArrayOfIndex& derivatives_data_position,
                                                      const LineRecord& line, 
                                                      ConstVectorView f_grid_full, 
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
                                                      const bool cutoff_call)
{  
  /* Single line shape solver
     
     The equation being solved here:
     F = LM NORM K1 K2 K3 QT/QT0 S0 f(...) - CUT(F),
     N = LM NORM K1 K2 (K4/K3 - 1) QT/QT0 S0 f(...) - CUT(N),
     dF = d(F)/dx
     dN = d(N)/dx
     
     or if using external population distribution calculations
     
     F = A f(...) - CUT(F)
     N = r f(...) = CUT(N)
     
     Where 
     LM:  Line mixing
     NORM: Normalization
     K1: Boltzmann ratio
     K2: Stimulated emission ratio
     K3: NLTE absorption ratio
     K4: NLTE emission ratio
     QT: Partition function
     QT0: Partition function at reference temperature
     S0: Line strength at reference temperature
     f(...): Line shape, including mirroring
     CUT(X): Computation of F and N at some provided cutoff frequency
     A is the absorption from external level distributions
     r is the emission ratio addition to Planck from external level distributions
  */ 
  
  /* Get the cutoff range if applicable
   * The line will have had a number set either by default or by finding the LSM CUT flag
   * followed by the numeric that is used to set the cutoff frequency.  By default, the 
   * cutoff frequency is set to a negative number, and it will only be used if there is
   * a non-negative number in the frequency range.
   */
  
  // Size and type of problem
  const Numeric& cutoff = line.CutOff();
  const bool need_cutoff = cutoff_call ? false : find_cutoff_ranges(this_f_range, f_grid_full, line.F(), cutoff);
  const bool do_temperature = do_temperature_jacobian(derivatives_data);
  const bool do_line_center = do_line_center_jacobian(derivatives_data);
                           
  // Leave this function if there is nothing to compute
  if(this_f_range.get_extent() == 0)
    return;
  
  // Extract the quantum identify of the line to be used in-case there are derivatives
  const QuantumIdentifier& QI = line.QuantumIdentity();
  
  /* Pressure broadening terms
   * These are set by the line catalog.  There are no defaults.
   */
  Numeric G0, G2, e, L0, L2, FVC;
  line.PressureBroadening().GetPressureBroadeningParams(
    G0, G2, e, L0, L2, FVC,
    temperature, line.Ti0(), pressure, partial_pressure, 
    this_species_location_in_tags, water_index_location_in_tags,
    broad_spec_locations, volume_mixing_ratio_of_all_species);
  
  // Line mixing terms
  Numeric Y, G, DV;
  line.LineMixing().GetLineMixingParams(Y, G, DV, temperature, pressure, 
                                        pressure_limit_for_linemixing);
  
  // Partial derivatives for temperature
  Numeric dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT, dY_dT, dG_dT, dDV_dT, dK1_dT, dK2_dT;
  if(do_temperature)
  {
    // Pressure broadening partial derivatives
    line.PressureBroadening().GetPressureBroadeningParams_dT(
      dG0_dT, dG2_dT, de_dT, dL0_dT, dL2_dT, dFVC_dT,
      temperature, line.Ti0(), pressure, partial_pressure,
      this_species_location_in_tags, water_index_location_in_tags,
      broad_spec_locations, volume_mixing_ratio_of_all_species);
    
    // Line mixing partial derivatives
    line.LineMixing().GetLineMixingParams_dT(
      dY_dT, dG_dT, dDV_dT, temperature,
      temperature_perturbation(derivatives_data),
      pressure, pressure_limit_for_linemixing);
  }
  
  /* Partial derivatives due to pressure
   * The vector below will be rescaled by the set internal derivatives function such that
   * the order of their occurrences are the same as in the partial derivative output
   */
  ComplexVector pressure_derivatives;
  if(do_pressure_jacobian(derivatives_data))
    line.PressureBroadening().SetInternalDerivatives(pressure_derivatives, derivatives_data, QI, 
                                                     line.Ti0()/temperature, pressure, partial_pressure, 
                                                     this_species_location_in_tags, water_index_location_in_tags,
                                                     volume_mixing_ratio_of_all_species);
    
  /* Partial derivatives due to line mixing
    * The vector below will be rescaled by the set internal derivatives function such that
    * the order of their occurrences are the same as in the partial derivative output
    */
  ComplexVector linemixing_derivatives;
  line.LineMixing().SetInternalDerivatives(linemixing_derivatives, derivatives_data, QI, 
                                            temperature, pressure, pressure_limit_for_linemixing);
  
  // Line shape usage remembering variable. 
  // Is only used if the user has set to use mirroring 
  // type to the same line shape as the main line.
  LineShapeType lst = LineShapeType::End;

  /*! Set the line shape normalized to unity integration
  * The user can set this by LSM LST followed by an index that 
  * is interpreted internally as a line shape.
  * The main point is not that the user should use such functions 
  * but that support functions can set the catalog, and that once
  * stored the catalog will use that line shape.  If no line shape 
  * tag is given, the line shape will be set by the type of pressure
  * broadening data that has been provided.
  */
  
  // Compute values
  ComplexVectorView F = F_full[this_f_range];
  ComplexMatrixView dF = dF_full(joker, this_f_range);
  ConstVectorView f_grid = f_grid_full[this_f_range];
  ComplexVectorView N = N_full[N_full.nelem() ? this_f_range : joker];
  ComplexMatrixView dN = dN_full(joker, N_full.nelem() ? this_f_range : joker);
  
  switch(line.GetLineShapeType())
  {
    // Use data as provided by the pressure broadening scheme
    case LineShapeType::ByPressureBroadeningData:
      switch(line.PressureBroadening().Type())
      {
        // Use data as per speed dependent air
        case PressureBroadeningData::PB_SD_AIR_VOLUME:
        case PressureBroadeningData::PB_HTP_AIR_VOLUME:
        case PressureBroadeningData::PB_PURELY_FOR_TESTING:
        case PressureBroadeningData::PB_SD_TEST_WATER:
          lst = LineShapeType::HTP;
          set_htp(F, dF, 
                  f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), doppler_constant, 
                  G0, L0, G2, L2, e, FVC,
                  derivatives_data, derivatives_data_position, QI,
                  ddoppler_constant_dT, 
                  dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT);
          break;
        // Use for data that requires air and water Voigt broadening
        case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
        // Use for data that requires planetary Voigt broadening
        case PressureBroadeningData::PB_PLANETARY_BROADENING:
          // Use for data that requires air Voigt broadening
        case PressureBroadeningData::PB_AIR_BROADENING:
        case PressureBroadeningData::PB_VOIGT_TEST_WATER:
          // Above should be all methods of pressure broadening requiring Voigt in ARTS by default
          lst = LineShapeType::Voigt;
          set_voigt(F, dF, f_grid, 
                    line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    line.F(), doppler_constant, 
                    G0, L0, DV, derivatives_data, derivatives_data_position, QI,
                    ddoppler_constant_dT, dG0_dT, dL0_dT, dDV_dT);
          break;
        case PressureBroadeningData::PB_NONE:
          throw std::runtime_error("Cannot understand the pressure broadening scheme");
      }
      break;
    // This line only needs the Doppler effect
    case LineShapeType::Doppler:
      lst = LineShapeType::Doppler;
      set_doppler(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), doppler_constant, derivatives_data, derivatives_data_position, QI, ddoppler_constant_dT);
      break;
    // This line only needs Hartmann-Tran
    case LineShapeType::HTP:
      set_htp(F, dF, 
              f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
              line.F(), doppler_constant, 
              G0, L0, G2, L2, e, FVC,
              derivatives_data, derivatives_data_position, QI,
              ddoppler_constant_dT, 
              dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT);
      lst = LineShapeType::HTP;
      break;
    // This line only needs Lorentz
    case LineShapeType::Lorentz:
      lst = LineShapeType::Lorentz;
      set_lorentz(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), G0, L0, DV, derivatives_data, derivatives_data_position, QI, dG0_dT, dL0_dT, dDV_dT);
      break;
    // This line only needs Voigt
    case LineShapeType::Voigt:
      lst = LineShapeType::Voigt;
      set_voigt(F, dF, f_grid, 
                line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                line.F(), doppler_constant, 
                G0, L0, DV, derivatives_data, derivatives_data_position, QI,
                ddoppler_constant_dT, dG0_dT, dL0_dT, dDV_dT);
      break;
    case LineShapeType::End:
      throw std::runtime_error("Cannot understand the requested line shape type.");
  }
  
  // Set the mirroring by repeating computations above using 
  // negative numbers for frequency of line related terms
  // The user sets if they want mirroring by LSM MTM followed by an index
  // that is interpreted as either mirroring by the same line shape or as 
  // mirroring by Lorentz lineshape
  switch(line.GetMirroringType())
  {
    // No mirroring
    case MirroringType::None:
      break;
    // Lorentz mirroring
    case MirroringType::Lorentz:
      {
        // Set the mirroring computational vectors and size them as needed
        ComplexVector Fm(F.nelem());
        ComplexMatrix dFm(dF.nrows(), dF.ncols());
        
        set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    -line.F(), G0, -L0, -DV, derivatives_data, derivatives_data_position, QI, dG0_dT, -dL0_dT, -dDV_dT);
        
        // Apply mirroring
        F -= Fm;
        dF -= dFm;
      }
      break;
    // Same type of mirroring as before
    case MirroringType::SameAsLineShape:
    {
      // Set the mirroring computational vectors and size them as needed
      ComplexVector Fm(F.nelem());
      ComplexMatrix dFm(dF.nrows(), dF.ncols());
      
      switch(lst)
      {
        case LineShapeType::Doppler:
          set_doppler(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      -line.F(), -doppler_constant, derivatives_data, derivatives_data_position, QI, -ddoppler_constant_dT);
          break;
        case LineShapeType::Lorentz:
          set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      -line.F(), G0, -L0, -DV, derivatives_data, derivatives_data_position, QI, dG0_dT, -dL0_dT, -dDV_dT);
          break;
        case LineShapeType::Voigt:
          set_voigt(Fm, dFm, f_grid, 
                    -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    -line.F(), -doppler_constant, 
                    G0, -L0, -DV, derivatives_data, derivatives_data_position, QI,
                    -ddoppler_constant_dT, dG0_dT, -dL0_dT, -dDV_dT);
          break;
        case LineShapeType::HTP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dFm, f_grid, 
                  -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  -line.F(), -doppler_constant, 
                  G0, -L0, G2, -L2, e, FVC,
                  derivatives_data, derivatives_data_position, QI,
                  -ddoppler_constant_dT, 
                  dG0_dT, -dL0_dT, dG2_dT, -dL2_dT, de_dT, dFVC_dT);
          break;
        case LineShapeType::ByPressureBroadeningData:
        case LineShapeType::End:
          throw std::runtime_error("Cannot understand the requested line shape type for mirroring.");
      }
      F -= Fm;
      dF -= dFm;
      break;
    }
    case MirroringType::End:
      throw std::runtime_error("Cannot understand the requested mirroring type for mirroring.");
  }
  
  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by and index
  // that is internally interpreted to mean some kind of lineshape normalization
  switch(line.GetLineNormalizationType())
  {
    // No normalization
    case LineNormalizationType::None:
      break;
      // van Vleck and Huber normalization
    case LineNormalizationType::VVH:
      apply_VVH_scaling(F, dF, f_grid, line.F(), temperature, derivatives_data, derivatives_data_position, QI);
      break;
      // van Vleck and Weiskopf normalization
    case LineNormalizationType::VVW:
      apply_VVW_scaling(F, dF, f_grid, line.F(), derivatives_data, derivatives_data_position, QI);
      break;
      // Rosenkranz's Quadratic normalization
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F, dF, f_grid, line.F(), temperature, derivatives_data, derivatives_data_position, QI);
      break;
    case LineNormalizationType::End:
      throw std::runtime_error("Cannot understand the requested line normalization type.");
  }
  
  // Apply line mixing if relevant
  if(Y not_eq 0 or G not_eq 0)
    apply_linemixing_scaling(F, dF, Y, G, derivatives_data, derivatives_data_position, QI, dY_dT, dG_dT);
  
  // Apply pressure broadening partial derivative vector if necessary
  if(pressure_derivatives.nelem()>0)
    apply_pressurebroadening_jacobian_scaling(dF, derivatives_data, derivatives_data_position, QI, pressure_derivatives);
  
  // Apply line mixing partial derivative vector if necessary
  if(linemixing_derivatives.nelem()>0)
    apply_linemixing_jacobian_scaling(dF, derivatives_data, derivatives_data_position, QI, linemixing_derivatives);
  
  // Apply line strength by whatever method is necessary
  switch(line.GetLinePopulationType())
  {
    case LinePopulationType::ByLTE:
    case LinePopulationType::ByVibrationalTemperatures:
      {
        // Line strength scaling that are line-dependent ---
        // partition functions are species dependent and computed at a higher level
        const Numeric gamma = stimulated_emission(temperature, line.F());
        const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
        const Numeric K1 = boltzman_ratio(temperature, line.Ti0(), line.Elow());
        const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
        
        // Line strength partial derivatives
        
        if(do_temperature)
        {
          dK1_dT = dboltzman_ratio_dT(K1, temperature, line.Elow());
          dK2_dT = dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), temperature);
        }
        
        // Partial derivatives due to central frequency of the stimulated emission
        Numeric dK2_dF0;
        if(do_line_center)
          dK2_dF0 = dstimulated_relative_emission_dF0(gamma, gamma_ref, temperature, line.Ti0());
        
        // Multiply the line strength by the line shape
        apply_linestrength_scaling(F, dF,  line.I0() * line.ZeemanEffect().StrengthScaling(zeeman_index), isotopologue_ratio,
                                   partition_function_at_temperature, partition_function_at_line_temperature, K1, K2,
                                   derivatives_data, derivatives_data_position, QI, dpartition_function_at_temperature_dT, dK1_dT, dK2_dT, dK2_dF0);
        
        if(line.GetLinePopulationType() == LinePopulationType::ByLTE)
          break;
        else if(nlte_distribution.nelem())
        {
          // Internal parameters
          Numeric Tu, Tl, K4, r_low, dK3_dF0, dK3_dT, dK3_dTl, dK4_dT, dK3_dTu, dK4_dTu;
          
          // These four are set by user on controlfile level
          // They are indexes to find the energy level in the nlte-temperature 
          // vector and the energy level of the states
          const Index evlow_index = line.NLTELowerIndex();
          const Index evupp_index = line.NLTEUpperIndex();
          const Numeric El = line.Evlow();
          const Numeric Eu = line.Evupp();
          
          // If the user set this parameters, another set of calculations are needed
          if(evupp_index > -1)
          {
            Tu = nlte_distribution[evupp_index];
            
            // Additional emission is from upper state
            K4 = boltzman_ratio(Tu, temperature, Eu);
          }
          // Otherwise the ratios are unity and nothing needs be done
          else
          {
            Tu = temperature;
            K4 = 1.0;
          }
          
          // The same as above but for the lower state level
          if(evlow_index > -1)
          {
            Tl = nlte_distribution[evlow_index];
            r_low = boltzman_ratio(Tl, temperature, El);
          }
          else
          {
            Tl = temperature;
            r_low = 1.0;
          }
          
          // Any additional absorption requires the ratio between upper and lower state number distributions
          const Numeric K3 = absorption_nlte_ratio(gamma, K4, r_low);
          
          // Are we computing the line center derivatives?
          if(do_line_center)
            dK3_dF0 = dabsorption_nlte_rate_dF0(gamma, temperature, K4, r_low);
          
          // Are we computing the temperature derivatives?
          // NOTE:  Having NLTE active AT ALL will change the jacobian because of this part of the code,
          // though this requires setting El and Eu for all lines, though this is not yet default...
          // So if you see this part of the code after having a runtime_error, 
          // you will need to write those functions yourself...
          if(do_temperature)
            dK3_dT = dabsorption_nlte_rate_dT(gamma, temperature, line.F(), El, Eu, K4, r_low);
          
          // Does the lower state level energy exist?
          if(El > 0)
            dK3_dTl = dabsorption_nlte_rate_dTl(gamma, temperature, Tl, El, r_low);
          
          // Does the upper state level energy exist?
          if(Eu > 0)
          {
            dK3_dTu = dabsorption_nlte_rate_dTu(gamma, temperature, Tu, Eu, K4);
            dK4_dTu = dboltzman_ratio_dT(K4, Tu, Eu);
          }
          
          // Apply this knowledge to set N and dN
          set_nonlte_source_and_apply_absorption_scaling(F, dF, N, dN, K3, K4, 
                                                         derivatives_data, derivatives_data_position,
                                                         QI,  dK3_dT, dK4_dT, dK3_dF0, dK3_dTl, 
                                                         dK3_dTu, dK4_dTu);
        }
      }
      break;
    case LinePopulationType::ByPopulationDistribution:
    {
        if(derivatives_data.nelem())
        {
          ostringstream os;
          os << "There will be no support for partial derivatives in population distribution\n" <<
                "mode before we can compute our own distributions.\n" << 
                "Use LTE if you know you have LTE distribution or ByVibrationalTemperatures if\n" <<
                "you can define vibrational temperatures for your problem.\n";
          throw std::runtime_error(os.str());
        }
        
        const Index nlte_low_index = line.NLTELowerIndex();
        const Index nlte_upp_index = line.NLTEUpperIndex();
        
        if(nlte_low_index < 0)
          throw std::runtime_error("No lower level distribution number in population distribution mode");
        if(nlte_upp_index < 0)
          throw std::runtime_error("No upper level distribution number in population distribution mode");
        
        apply_linestrength_from_nlte_level_distributions(F, dF, N, dN,
                                                         nlte_distribution[nlte_low_index],
                                                         nlte_distribution[nlte_upp_index],
                                                         line.G_lower(), line.G_upper(),
                                                         line.A(), line.F(),
                                                         temperature, derivatives_data, 
                                                         derivatives_data_position, QI);
      }
      break;
    case LinePopulationType::End: 
      throw std::runtime_error("Cannot understand the line strength computations");
  }
  
  // Cutoff frequency is applied at the end because 
  // the entire process above is complicated and applying
  // cutoff last means that the code is kept cleaner
  if(need_cutoff)
  {
    apply_cutoff(F, dF, N, dN,
                 derivatives_data, derivatives_data_position, line,
                 volume_mixing_ratio_of_all_species,
                 nlte_distribution, pressure, temperature,
                 doppler_constant, partial_pressure,
                 isotopologue_ratio, magnetic_magnitude,
                 ddoppler_constant_dT, pressure_limit_for_linemixing,
                 partition_function_at_temperature,
                 dpartition_function_at_temperature_dT,
                 partition_function_at_line_temperature,
                 broad_spec_locations, this_species_location_in_tags,
                 water_index_location_in_tags, zeeman_index, verbosity);
  }
}


/*!
 * Combination function using standard setup to compute line strength and lineshape of a single line.
 * Computes in order the lineshape, the linemirroring, the linenormalization, the linemixing, the 
 * linestrength, the non-lte, and the cutoff frequency.
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * \retval N Non-lte lineshape
 * \retval dN Non-lte lineshape derivative
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param line line-record containing most line parameters
 * \param volume_mixing_ratio_of_all_species As name suggests
 * \param nlte_distribution As name suggests
 * \param pressure As name suggests
 * \param temperature As name suggests
 * \param doppler_constant Frequency-independent part of the Doppler broadening
 * \param partial_pressure Pressure of species that line belongs to at this level
 * \param isotopologue_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param ddoppler_constant_dT Temperature derivative of doppler_constant
 * \param pressure_limit_for_linemixing As WSV lm_p_lim
 * \param partition_function_at_temperature As name suggests
 * \param dpartition_function_at_temperature_dT Temeperature derivative of partition_function_at_temperature
 * \param partition_function_at_line_temperature As name suggests
 * \param broad_spec_locations Locations of broadening species using all-planetary broadening inside volume_mixing_ratio_of_all_species
 * \param this_species_location_in_tags Location of species of line in volume_mixing_ratio_of_all_species
 * \param water_index_location_in_tags Location of water in volume_mixing_ratio_of_all_species
 * \param df_range Frequency range to use inside dF and dN
 * \param verbosity Verbosity level
 * 
 */
void Linefunctions::apply_cutoff(ComplexVectorView F,
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
                                 const Verbosity& verbosity)
{ 
  // Size of derivatives
  const Index nj = dF.nrows(); 
  const Index nn = dN.nrows(); 
  
  // Setup compute variables
  Vector f_grid_cutoff(1, line.F() + line.CutOff());
  ComplexVector Fc(1), Nc(1);
  ComplexMatrix dFc(nj, 1), dNc(nn, 1);
  Range tmp(joker);
  
  // Recompute the line for a single frequency
  set_cross_section_for_single_line(Fc, dFc, Nc, dNc, tmp,
                                    derivatives_data, derivatives_data_position, line, f_grid_cutoff,
                                    volume_mixing_ratio_of_all_species,
                                    nlte_distribution, pressure, temperature,
                                    doppler_constant, partial_pressure, isotopologue_ratio,
                                    magnetic_magnitude, ddoppler_constant_dT,
                                    pressure_limit_for_linemixing,
                                    partition_function_at_temperature,
                                    dpartition_function_at_temperature_dT,
                                    partition_function_at_line_temperature,
                                    broad_spec_locations, this_species_location_in_tags,
                                    water_index_location_in_tags, zeeman_index, verbosity, true);
  
  // Apply cutoff values
  F -= Fc[0];
  if(N.nelem())
    N -= Nc[0];
  #pragma omp simd
  for(Index i = 0; i < nj; i++)
    dF(i, joker) -= dFc(i, 0);
  #pragma omp simd
  for(Index i = 0; i < nn; i++)
    dN(i, joker) -= dNc(i, 0);
}


bool Linefunctions::find_cutoff_ranges(Range& range,
                                       const ConstVectorView& f_grid,
                                       const Numeric& F0,
                                       const Numeric& cutoff)
{
  const Index nf = f_grid.nelem();
  
  const bool need_cutoff = (cutoff > 0);
  if(need_cutoff)
  {
    // Find range of simulations
    Index i_f_min = 0;
    Index i_f_max = nf-1;
    
    // Loop over positions to compute the line shape cutoff point
    while(i_f_min < nf and (F0 - cutoff) > f_grid[i_f_min])
      ++i_f_min;
    while(i_f_max >= i_f_min and (F0 + cutoff) < f_grid[i_f_max])
      --i_f_max;
    
    //  The extent is one more than the difference between the indices of interest
    const Index extent = i_f_max - i_f_min + 1; // min is 0, max is nf
    
    range = Range(i_f_min, extent);
  }
  else
  {
    range = Range(joker);
  } 
  return need_cutoff;
}


/*!
 * Applies non-lte linestrength to already set line shape
 * 
 * Works on ratio-inputs, meaning that the total distribution does not have to be known
 * 
 * Cannot support partial derivatives at this point due to ARTS not possessing its own
 * NLTE ratio calculation agenda
 * 
 * \retval F Lineshape (absorption)
 * \retval N Non-lte lineshape (source)
 * 
 * \param r1 Ratio of molecules at energy level 1
 * \param r2 Ratio of molecules at energy level 2 
 * \param g1 Statistical weight of energy level 1
 * \param g2 Statistical weight of energy level 2
 * \param A21 Einstein coefficient for the transition from energy level 2 to energy level 1
 * \param F0 Central frequency without any shifts
 * \param T Atmospheric temperature at level
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF and dN
 * 
 */
void Linefunctions::apply_linestrength_from_nlte_level_distributions(ComplexVectorView F, 
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
                                                                     const ArrayOfRetrievalQuantity& derivatives_data,
                                                                     const ArrayOfIndex& derivatives_data_position,
                                                                     const QuantumIdentifier& quantum_identity)
{
  // Size of the problem
  const Index nf = F.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  // Physical constants
  const static Numeric c0 = 2.0 * PLANCK_CONST / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
  const static Numeric c1 = PLANCK_CONST / 4 / PI;
  
  // Constants based on input
  const Numeric c2 = c0 * F0 * F0 * F0;
  const Numeric c3 = c1 * F0;
  const Numeric x = g2 / g1;
  
  /*
    Einstein 'active' coefficients are as follows:
    const Numeric B21 =  A21 / c2;
    const Numeric B12 = x * B21;
  */
  
  // Planck function of this line
  const Numeric exp_T = exp(PLANCK_CONST * F0 / BOLTZMAN_CONST / T),  b = c2/(exp_T - 1);
  
  // Absorption strength
  const Numeric k = c3 * (r1*x - r2) * (A21 / c2);
  
  // Emission strength
  const Numeric e = c3 * r2 * A21;
  
  const Numeric ratio = (e/b-k);
  
  // Partial derivatives
  for(Index iq=0; iq<nppd; iq++)
  {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)  // nb.  Only here tou counter-act the latter on in 
    { 
      Numeric done_over_b_dT = PLANCK_CONST*F0*exp_T/(c2*BOLTZMAN_CONST*T*T);
      
      // dN is unset so far.  It should return to just be lineshape later...
      #pragma omp simd
      for(Index iv=0; iv<nf; iv++)
        dN(iq, iv) = F[iv]*e*done_over_b_dT + dF(iq, iv)*ratio;
      
      // dk_dT = 0...
      dF(iq, joker) *= k;
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]]))
    {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        Numeric done_over_b_df0 = PLANCK_CONST*exp_T/(c2*BOLTZMAN_CONST*T) - 3.0*b/F0;
        Numeric de_df0 = c1 * r2 * A21;
        Numeric dk_df0 = c1 * (r1*x - r2) * (A21 / c2) - 3.0*k/F0;
        
        #pragma omp simd
        for(Index iv=0; iv<nf; iv++)
        {
          dN(iq, iv) = F[iv]*(e*done_over_b_df0 + de_df0/b - dk_df0) + dF(iq, iv)*ratio;
          dF(iq, iv) = dF(iq, iv)*k + F[iv]*dk_df0;
        }
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::NLTE)
    {
      if(quantum_identity.LowerQuantumId() > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        const Numeric dk_dr2 = - c3 * A21 / c2, de_dr2 = c3 * A21, dratio_dr2 = de_dr2/b - dk_dr2;
        
        #pragma omp simd
        for(Index iv=0; iv<nf; iv++)
        {
          dN(iq, iv) = F[iv]*dratio_dr2;
          dF(iq, iv) = F[iv]*dk_dr2 + dF(iq, iv)*k;
        } 
      }
      else if(quantum_identity.UpperQuantumId() > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
      {
        const Numeric dk_dr1 = c3 * x * A21 / c2;
        
        #pragma omp simd
        for(Index iv=0; iv<nf; iv++)
          dF(iq, iv) = F[iv]*dk_dr1 + dF(iq, iv)*k;
      }
    }
  }
  
  // Set source function to be the relative amount emitted by the line divided by the Planck function
  N = F;
  N *= ratio;
  
  // Set absorption
  F *= k;
}


Index Linefunctions::first_binary_level(const Numeric& dg, const Numeric& df, const Index N)
{
  const Numeric p = dg / df;
  Index i = 0;
  while(p > ((1 << i) * N)) i++;  // nb, N is the minimum number of points within range of p
  return i;
}


Numeric Linefunctions::binary_central_step_size(const Numeric& ds, const Index i, const Index i0)
{
  return ds * (1 << (i - i0));
}


Index Linefunctions::binary_frequency_position(const Numeric& df, 
                                               const Numeric& dg,
                                               const Numeric& f,
                                               const Numeric& F0,
                                               const Index n,
                                               const Index dn)
{
  // Frequency equivalen Index 
  const Numeric x = (F0 + dg - f) / df;
  Index i;
  if(x < 0)
    i = 0;
  else if(x > Numeric(n))
    i = n-1;
  else
    i = Index(x);
  
  // Adjust position by grid
  const Index DN = abs(dn);
  const Index sgn = dn / DN;
  while(i < n-1 and i > 0 and (i % 2 or i % DN)) i += sgn;
  return i;
}


/*! Returns binary_speedup x 6 array of array of indexes
 * 
 * These indexes are defined such that:
 *       A[i][0] is the lower bound of the lower range
 *       A[i][1] is the upper bound of the lower range
 *       A[i][2] is the stepsize of lower range (used as a check)
 *       A[i][3] is the lower bound of the upper range
 *       A[i][4] is the upper bound of the upper range
 *       A[i][5] is the stepsize of upper range (used as a check)
 * 
 * Rules:
 *       0)     -1 is used to indicate that values do not exist except...
 * 
 *       1)     If A[i][3] =:= A[i][4] =:= A[i][5] =:= -1, 
 *       1)     but A[i][0] =!= A[i][1] =!= A[i][2] =!= -1, 
 *       1)     then the otherwise lower range values represents
 *       1)     a central range that is extrapolated, if 
 *       1)     possible, at both ends.  Note that these should also have 
 *       1)     A[i][0] % A[i][2] =:= A[i][1] % A[i][2] =:= 0
 * 
 *       2)     If A[i][3] =:= A[i][4] =:= A[i][5] =:= N, 
 *       2)     then there should be a linear interpolation from N 
 *       2)     to the upper edge.  Note that N == 1 + 2^binary_speedup
 *       2)     requires no interpolation!
 * 
 *       3)     If A[i][0] =:= A[i][1] =:= A[i][2] =:= N,
 *       3)     then there should be a linear interpolation from 0 to N.
 *       3)     Note that N == 0 requires no interpolation!
 * 
 *       4)     These are expected to be true in all other cases
 *              a) A[i][j] >= 0 and A[i][j] <= 2^binary_speedup
 *              b) A[i][3] := A[i-1][4] + A[i-1][5]; 
 *              c) A[i][0] := A[i-1][1] + A[i-1][2]; 
 *              d) A[i][5] := 2 * A[i-1][5]; 
 *              e) A[i][0] := 2 * A[i-1][2]; 
 *   
 */
ArrayOfArrayOfIndex Linefunctions::binary_boundaries(const Numeric& F0, 
                                                     const ConstVectorView& f_grid,
                                                     const Numeric& G0,
                                                     const Numeric& GD_div_F0,
                                                     const Numeric& binary_speedup_coef,
                                                     const Index binary_speedup,
                                                     const Index binary_frequency_count,
                                                     const Index steps)
{
  assert(steps > 1 and steps % 2);
  
  const Index& nf  = f_grid.nelem();
  const Numeric ds = binary_speedup_coef * sqrt(GD_div_F0*GD_div_F0*F0*F0 + G0*G0);
  
  const Numeric df = f_grid[1] - f_grid[0];
  const Index i0   = first_binary_level(binary_central_step_size(ds, 0, 0), df, binary_frequency_count);
  
  bool any_count=false;
  Index l=-1, u=-1, ol=-1, ou=-1, stepsize=1<<i0;
  ArrayOfArrayOfIndex A(binary_speedup, ArrayOfIndex(6, -1));
  for(Index i=i0; i<binary_speedup; i++) {
    if(i == binary_speedup - 1 and not any_count) {
      A[i][0] = 0;
      A[i][1] = nf-1;
      A[i][2] = nf/2;
    }
    else if(i == i0) {
      // Distance of a signle broadening
      const Numeric dg = binary_central_step_size(ds, i, i0);
      
      // l is nf-1, 0, or a number inbetween
      A[i][0] = l = ol = binary_frequency_position(df, -dg, f_grid[0], F0, nf, -stepsize);
      
      // u is nf-1, 0, or a number inbetween
      A[i][1] = u = ou = binary_frequency_position(df, dg, f_grid[0], F0, nf, stepsize);
      
      // original stepsize
      A[i][2] = stepsize;
      
      // Flag the others as -1 to indicate that this represents a central range
      A[i][3] = A[i][4] = A[i][5] = -1;
    }
    else {
      // The next steps ideally begins 1 old stepszie away
      u = ou + stepsize;
      l = ol - stepsize;
      // Now we can double the stepsize
      stepsize <<= 1;
      
      // Fix upper positions
      if(u >= nf-1) {
        // The upper range is already outside the range of importance so we set it to be zero-long value at it
        A[i][3] = A[i][4] = ou;
        ou = u = nf-1;
        A[i][5] = 0;
      }
      else {  // There are relevant steps
        if(ou == 0)
          // the central range was outside the lower edge so we initialize calculations at position 0
          // FIXME:  Should we use a distance-based beginning instead?
          A[i][3] = u = 0;
        else
          // the range will be start at what is currently a half-step away from the previous range
          A[i][3] = u;
        
        // The upper range should contain steps indexes above the beginning if possible
        u = A[i][3] + (steps-1)*stepsize;
        
        // If this puts the range outside the reach of the grid, we move back a step at a time
        while(u >= nf) u -= stepsize;
        
        if(A[i][3] == u)
          // If these numbers are the same, then this is the point that will be used to interpolate to the upper edge
          A[i][4] = A[i][5] = ou = u = A[i][3];
        else {
          // Otherwise we have at least 2 points and we move on
          A[i][4] = ou = u;
          A[i][5] = stepsize;
        }
      }

      // Fix lower positions
      if(l <= 0) {
        // The lower range is already at the limit and will not do anything more
        A[i][0] = A[i][1] = ol;
        A[i][2] = 0;
        l = ol = 0;
      }
      else {
        if(ol >= nf-1)
          // the central range was outside the upper edge so we initialize calculations at position nf-1
          // FIXME:  Should we use a distance-based beginning instead?
          A[i][1] = l = nf-1;
        else
          // the range will be start at what is currently a half-step away from the previous range
          A[i][1] = l;
        
        // The lower range should contain steps indexes below the current one if possible
        l = A[i][1] - (steps-1)*stepsize;
        
        // Increment by stepsize until we are within the range of f_grid
        while(l < 0) l += stepsize;
        
        if(A[i][1] == l) {
          // We have reached the point where we need a linear interpolation to the lower edge
          A[i][0] = A[i][2] = ol = l = A[i][1];
        }
        else{
          // Otherwise we have at least 2 points and we move on
          A[i][0] = ol = l;
          A[i][2] = stepsize;
        }
      }
    }
    
    if(l not_eq u)
      any_count = true;
  }
  
  return A;
}


/*! Lagrange functions for even-spaced functions
 * 
 * [y1 ··· yN] are assumed found for values of x between [0 ··· N-1], where x1 gets y1 and so on.
 * 
 * Code copy-pasted from sympy so errors are possible
 */ 
inline Complex evenspaced_interp_9point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3, const Complex& y4, const Complex& y5,
                                                 const Complex& y6, const Complex& y7, const Complex& y8,
                                                 const Complex& y9) noexcept
{
  return 
  + y1 * (  (x - 8)*(x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 40320)
  - y2 * (x*(x - 8)*(x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)         / 5040)
  + y3 * (x*(x - 8)*(x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*        (x - 1) / 1440)
  - y4 * (x*(x - 8)*(x - 7)*(x - 6)*(x - 5)*(x - 4)*        (x - 2)*(x - 1) / 720)
  + y5 * (x*(x - 8)*(x - 7)*(x - 6)*(x - 5)*        (x - 3)*(x - 2)*(x - 1) / 576)
  - y6 * (x*(x - 8)*(x - 7)*(x - 6)*        (x - 4)*(x - 3)*(x - 2)*(x - 1) / 720)
  + y7 * (x*(x - 8)*(x - 7)*        (x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 1440)
  - y8 * (x*(x - 8)*        (x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 5040)
  + y9 * (x*        (x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 40320);
}
inline Complex evenspaced_interp_8point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3, const Complex& y4, const Complex& y5,
                                                 const Complex& y6, const Complex& y7, const Complex& y8) noexcept
{
  return
  - y1 * (  (x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 5040)
  + y2 * (x*(x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)         / 720)
  - y3 * (x*(x - 7)*(x - 6)*(x - 5)*(x - 4)*(x - 3)*        (x - 1) / 240)
  + y4 * (x*(x - 7)*(x - 6)*(x - 5)*(x - 4)*        (x - 2)*(x - 1) / 144)
  - y5 * (x*(x - 7)*(x - 6)*(x - 5)*        (x - 3)*(x - 2)*(x - 1) / 144)
  + y6 * (x*(x - 7)*(x - 6)*        (x - 4)*(x - 3)*(x - 2)*(x - 1) / 240)
  - y7 * (x*(x - 7)*        (x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 720)
  + y8 * (x*        (x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 5040);
}
inline Complex evenspaced_interp_7point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3, const Complex& y4, const Complex& y5,
                                                 const Complex& y6, const Complex& y7) noexcept
{
  return
  + y1 * (  (x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 720)
  - y2 * (x*(x - 6)*(x - 5)*(x - 4)*(x - 3)*(x - 2)         / 120) 
  + y3 * (x*(x - 6)*(x - 5)*(x - 4)*(x - 3)*        (x - 1) / 48)
  - y4 * (x*(x - 6)*(x - 5)*(x - 4)*        (x - 2)*(x - 1) / 36)
  + y5 * (x*(x - 6)*(x - 5)*        (x - 3)*(x - 2)*(x - 1) / 48)
  - y6 * (x*(x - 6)*        (x - 4)*(x - 3)*(x - 2)*(x - 1) / 120)
  + y7 * (x*        (x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 720);
}
inline Complex evenspaced_interp_6point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3, const Complex& y4, const Complex& y5,
                                                 const Complex& y6) noexcept
{
  return
  - y1 * (  (x - 5)*(x - 4)*(x - 3)*(x - 2)*(x - 1) / 120)
  + y2 * (x*(x - 5)*(x - 4)*(x - 3)*(x - 2)         / 24)
  - y3 * (x*(x - 5)*(x - 4)*(x - 3)*        (x - 1) / 12)
  + y4 * (x*(x - 5)*(x - 4)*        (x - 2)*(x - 1) / 12)
  - y5 * (x*(x - 5)*        (x - 3)*(x - 2)*(x - 1) / 24)
  + y6 * (x*        (x - 4)*(x - 3)*(x - 2)*(x - 1) / 120);
}
inline Complex evenspaced_interp_5point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3, const Complex& y4, const Complex& y5) noexcept
{
  return
  + y1 * (  (x - 4)*(x - 3)*(x - 2)*(x - 1) / 24)
  - y2 * (x*(x - 4)*(x - 3)*(x - 2)         / 6)
  + y3 * (x*(x - 4)*(x - 3)*        (x - 1) / 4)
  - y4 * (x*(x - 4)*        (x - 2)*(x - 1) / 6)
  + y5 * (x*        (x - 3)*(x - 2)*(x - 1) / 24);
}
inline Complex evenspaced_interp_4point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3, const Complex& y4) noexcept
{
  return
  - y1 * (  (x - 3)*(x - 2)*(x - 1) / 6)
  + y2 * (x*(x - 3)*(x - 2)         / 2)
  - y3 * (x*(x - 3)*        (x - 1) / 2)
  + y4 * (x*        (x - 2)*(x - 1) / 6);
}
inline Complex evenspaced_interp_3point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2,
                                                 const Complex& y3) noexcept
{
  return
  + y1 * (  (x - 2)*(x - 1) / 2)
  - y2 * (x*(x - 2))
  + y3 * (x*        (x - 1) / 2);
}
inline Complex evenspaced_interp_2point_lagrange(const Numeric& x, const Complex& y1, const Complex& y2) noexcept
{
  return y1 + (y2-y1)*x;
}


inline void central_binary_interpolation_step(ComplexVectorView f, const Index& l, const Index& u, const Index& s) noexcept
{
  const Index nf = f.nelem();
  
  // Stepsize
  const Numeric dx=1.0/Numeric(s);
  Numeric x;  // stepper
  
  // Are there values to extrapolate for "anchoring" the interpolation
  const bool upper_extrapolation = u + s <= nf-1;
  const bool lower_extrapolation = l - s >= 0;
  const Index nelem = 1 + (u - l) / s;
  
  if(u < 1 or l > nf-2 or s == 1) {
    // Do nothing if the range is far away or is the main range
  }
  else if(nelem == 1) {
    // These functions must have extrapolation on either side to do anything in the end
    
    x = dx;
    
    // When both upper and lower extrapolation points are included, there are 3 points and the interpolation should be between these three
    if(upper_extrapolation and lower_extrapolation) {
      const Complex &y1=f[l-s], &y2=f[l], &y3=f[u+s];
      for(Index i = l - s + 1; i < u + s; i++) {
        if(i not_eq l)
          f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    
    // When only the upper or lower extrapolation point is available, then we can only perform linera interpolation
    else if(upper_extrapolation) {
      const Complex &y1=f[u], &y2=f[u+s];
      for(Index i = l + 1; i < u + s; i++) {
        f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
        x += dx;
      }
    }
    else if(lower_extrapolation) {
      const Complex &y1=f[l-s], &y2=f[l];
      for(Index i = l - s + 1; i < u; i++) {
        f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
        x += dx;
      }
    }
  }
  else if(nelem == 2) {
    // These functions must be interpolated between eachother and also to the available extrapolation points
    
    x = dx;
    
    // When both extrapolation points are available, there are four points to interpolate between
    if(upper_extrapolation and lower_extrapolation) {
      const Complex &y1=f[l-s], &y2=f[l], &y3=f[u], &y4=f[u+s];
      for(Index i = l - s + 1; i < u + s; i++) {
        if(i not_eq l and i not_eq u)
          f[i] = evenspaced_interp_4point_lagrange(x, y1, y2, y3, y4);
        x += dx;
      }
    }
    
    // When only upper or lower exists outside the range, then there are only three points to interpolate between
    else if(upper_extrapolation) {
      const Complex &y1=f[l], &y2=f[u], &y3=f[u+s];
      for(Index i = l + 1; i < u + s; i++) {
        if(i not_eq u)
          f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    else if(lower_extrapolation) {
      const Complex &y1=f[l-s], &y2=f[l], &y3=f[u];
      for(Index i = l - s + 1; i < u; i++) {
        if(i not_eq l)
          f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    
    // Otherwise we have only linear interpolation between the available points
    else {
      const Complex &y1=f[l], &y2=f[u];
      for(Index i = l + 1; i < u; i++) {
        f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
        x += dx;
      }
    }
  }
  else {
    // With more than three points, we can use the interpolation regime at the center of four points
    // for the internals, and only extrapolate to the outer ranges from the closes two internals
    // We already know that there are enough points internally if the extrapolation points are available
    
    // Interpolate at X in "l|X|X| |" or "|X| |", where l| is the lower extrapolation point
    if(lower_extrapolation) {
      x = dx;
      const Complex &y1=f[l-s], &y2=f[l], &y3=f[l+s], &y4=f[l+2*s];
      for(Index i = l - s + 1; i < l+s; i++) {
        if(i not_eq l)
          f[i] = evenspaced_interp_4point_lagrange(x, y1, y2, y3, y4);
        x += dx;
      }
    }
    else {
      x = dx;
      const Complex &y1=f[l], &y2=f[l+s], &y3=f[l+2*s];
      for(Index i = l + 1; i < l + s; i++) {
        f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    
    // Interpolate at X in "| |X|X|u" or "| |X|", where |u is the upper extrapolation point
    if(upper_extrapolation) {
      x = 1.0 + dx;
      const Complex &y1=f[u-2*s], &y2=f[u-s], &y3=f[u], &y4=f[u+s];
      for(Index i = u - s + 1; i < u + s; i++) {
        if(i not_eq u)
          f[i] = evenspaced_interp_4point_lagrange(x, y1, y2, y3, y4);
        x += dx;
      }
    }
    else {
      x = 1.0 + dx;
      const Complex &y1=f[u-2*s], &y2=f[u-s], &y3=f[u];
      for(Index i = u - s + 1; i < u; i++) {
        if(i not_eq u)
          f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    
    // Interpolation between all internal points happen at X in "| |X| |"
    for(Index j=1; j<nelem-1; j++) {
      const Index p = l + j*s;
      const Complex &y1=f[p-s], &y2=f[p], &y3=f[p+s], &y4=f[p+2*s];
      x = 1.0 + dx;
      for(Index i = p + 1; i < p + s; i++) {
        f[i] = evenspaced_interp_4point_lagrange(x, y1, y2, y3, y4);
        x += dx;
      }
    }
  }
}


inline void lower_binary_interpolation_step(ComplexVectorView f, const Index& l, const Index& u, const Index& s) noexcept
{
  if(l == u and s == 0 and l) {
    // The stepsize in relative Numerical values
    const Numeric dx = 1.0 / Numeric(l);
    Numeric x = dx;  // interpolation coefficient
    
    // The lower range is to be extrapolated from 0 to l
    Complex &y1=f[0], &y2=f[l];
    for(Index i=1; i<l; i++) {
      f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
      x += dx;
    }
  }
  else {
    // This is a normal range that has to be interpolated from start to finish
    
    // The stepsize in relative Numerical values
    const Numeric dx = 1.0 / Numeric(s);
    Numeric x = dx;  // interpolation coefficient
    
    // Is the extrapolation point below available?
    const bool lower_extrapolation = l-s >= 0;
    
    // There is one extra point if there is a lower extrapolation point
    const Index nelem = (u-l) / s + 1 + (lower_extrapolation?1:0);
    
    // We need the start to be a stepsize below if the lower extrapolation point is there
    const Index start = l - (lower_extrapolation?s:0);
    
    // We sadly need special rules for the number of points here so that the references can be set ahead of the loop...
    if(nelem == 2) {
      const Complex &y1=f[start], &y2=f[start+s];
      for(Index i=start+1; i<u; i++) {
        f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
        x += dx;
      }
    }
    else if(nelem == 3) {
      const Complex &y1=f[start], &y2=f[start+s], &y3=f[start+s+s];
      for(Index i=start+1; i<u; i++) {
        f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    else if(nelem == 4) {
      const Complex &y1=f[start], &y2=f[start+s], &y3=f[start+s+s], &y4=f[start+s+s+s];
      for(Index i=start+1; i<u; i++) {
        f[i] = evenspaced_interp_4point_lagrange(x, y1, y2, y3, y4);
        x += dx;
      }
    }
    /* else if(nelem == 5) add these if the code changes later */
  }
}


inline void upper_binary_interpolation_step(ComplexVectorView f, const Index& l, const Index& u, const Index& s) noexcept
{
  const Index nf=f.nelem();
  
  if(l == u and s == 0 and l) {
    // The stepsize in relative Numerical values
    const Numeric dx = 1.0 / Numeric(nf-l);
    Numeric x = dx;  // interpolation coefficient
    
    // The upper range is to be extrapolated from u to nf
    Complex &y1=f[l], &y2=f[nf-1];
    for(Index i=l+1; i<nf; i++) {
      f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
      x += dx;
    }
  }
  else {
    // This is a normal range that has to be interpolated from start to finish
    
    // The stepsize in relative Numerical values
    const Numeric dx = 1.0 / Numeric(s);
    Numeric x = dx;  // interpolation coefficient
    
    // Is the extrapolation point above available?
    const bool upper_extrapolation = u+s < nf;
    
    // There is one extra point if there is a upper extrapolation point
    const Index nelem = (u-l) / s + 1 + (upper_extrapolation?1:0);
    
    // We need the end to be a stepsize above if the upper extrapolation point is there
    const Index end = u + (upper_extrapolation?s:0);
    
    // We sadly need special rules for the number of points here so that the references can be set ahead of the loop...
    if(nelem == 2) {
      const Complex &y1=f[l], &y2=f[l+s];
      for(Index i=l+1; i<end; i++) {
        f[i] = evenspaced_interp_2point_lagrange(x, y1, y2);
        x += dx;
      }
    }
    else if(nelem == 3) {
      const Complex &y1=f[l], &y2=f[l+s], &y3=f[l+s+s];
      for(Index i=l+1; i<end; i++) {
        f[i] = evenspaced_interp_3point_lagrange(x, y1, y2, y3);
        x += dx;
      }
    }
    else if(nelem == 4) {
      const Complex &y1=f[l], &y2=f[l+s], &y3=f[l+s+s], &y4=f[l+s+s+s];
      for(Index i=l+1; i<end; i++) {
        f[i] = evenspaced_interp_4point_lagrange(x, y1, y2, y3, y4);
        x += dx;
      }
    }
    /* else if(nelem == 5) add these if the code changes later */
  }
}


inline void binary_interpolation_step(ComplexVectorView f, const ArrayOfIndex& bounds) noexcept
{
  const Index lupp=bounds[3], uupp=bounds[4], supp=bounds[5];  // Upper range
  const Index llow=bounds[0], ulow=bounds[1], slow=bounds[2];  // Lower range
  const Index nf = f.nelem();
  
  if(llow == -1 and ulow == -1 and slow == -1) {
    /* Nothing happens here because the range is considered too dense for the given grid */ 
  }
  else if(lupp == -1 and uupp == -1 and supp == -1) {
    central_binary_interpolation_step(f, llow, ulow, slow);
  }
  else {
    if(llow < nf and llow >= 0 and ulow)
      lower_binary_interpolation_step(f, llow, ulow, slow);
    if(uupp >= 0 and uupp < nf)
      upper_binary_interpolation_step(f, lupp, uupp, supp);
  }
}


void Linefunctions::binary_interpolation(ComplexVectorView f, const ArrayOfArrayOfIndex& binary_bounds)
{
  for(const ArrayOfIndex& bounds : binary_bounds) {
    assert(bounds.nelem() == 6);
    assert(f.nelem() > bounds[0]); assert(f.nelem() > bounds[1]); assert(f.nelem() > bounds[2]);
    assert(f.nelem() > bounds[3]); assert(f.nelem() > bounds[4]); assert(f.nelem() > bounds[5]);
    binary_interpolation_step(f, bounds);
  }
}


Range Linefunctions::binary_level_range(const ArrayOfArrayOfIndex& binary_bounds, const Index nf, const Index i, const bool lower_range)
{
  assert(i < binary_bounds.nelem());
  Index l, u, s;
  
  if(lower_range) {
    l = binary_bounds[i][0];
    u = binary_bounds[i][1];
    s = binary_bounds[i][2];
  }
  else {
    l = binary_bounds[i][3];
    u = binary_bounds[i][4];
    s = binary_bounds[i][5];
  }
  
  if(l == -1 and u == -1 and s == -1)
    return Range(0, 0);
  else if(lower_range and l == u and s == 0 and l not_eq 0)
    return Range(0, 1);
  else if(lower_range and l == nf-1)
    return Range(nf-1, 0);
  else if(not lower_range and l == nf-1)
    return Range(nf-1, 0);
  else if(not lower_range and l == u and s == 0)
    return Range(nf-1, 1);
  else if(s==0)
    return Range(0, 0);
  else
    return Range(l, 1 + (u-l)/s, s);
}


void Linefunctions::set_lineshape_from_level_line_data(Complex& F,
                                                       Complex& N,
                                                       ComplexVectorView dF,
                                                       ComplexVectorView dN,
                                                       const Numeric& f,
                                                       const Numeric& T,
                                                       const SingleLevelLineData& level_line_data,
                                                       const LineRecord& line,
                                                       const ArrayOfRetrievalQuantity& derivatives_data,
                                                       const ArrayOfIndex& derivatives_data_position) noexcept
{
  // FIXME:  Not compatible with Zeeman effect yet (perhaps as easy as moving a loop to a lower level?)
  
  switch(line.GetLineShapeType()) {
    case LineShapeType::Doppler:
      set_doppler_from_level_line_data(F, dF, f, line.F(), 0.0, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineShapeType::Lorentz:
      set_lorentz_from_level_line_data(F, dF, f, line.F(), 0.0, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineShapeType::Voigt:
      set_voigt_from_level_line_data(F, dF, f, line.F(), 0.0, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineShapeType::HTP:
      set_htp_from_level_line_data(F, dF, f, line.F(), 0.0, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineShapeType::ByPressureBroadeningData:
      switch(line.PressureBroadening().Type()) {
        case PressureBroadeningData::PB_SD_AIR_VOLUME:
        case PressureBroadeningData::PB_HTP_AIR_VOLUME:
        case PressureBroadeningData::PB_PURELY_FOR_TESTING:
        case PressureBroadeningData::PB_SD_TEST_WATER:
          set_htp_from_level_line_data(F, dF, f, line.F(), 0.0, level_line_data, derivatives_data, derivatives_data_position);
          break;
        case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
        case PressureBroadeningData::PB_PLANETARY_BROADENING:
        case PressureBroadeningData::PB_AIR_BROADENING:
        case PressureBroadeningData::PB_VOIGT_TEST_WATER:
          set_voigt_from_level_line_data(F, dF, f, line.F(), 0.0, level_line_data, derivatives_data, derivatives_data_position);
          break;
        case PressureBroadeningData::PB_NONE:
          F = Complex(1.0, 0.0);
          dF = Complex(0.0, 0.0);
      }
      break;
    case LineShapeType::End:
      break;
  }
  
  switch(line.GetMirroringType()) {
    case MirroringType::Lorentz:
      {
        Complex Fm;
        ComplexVector dFm(dF.nelem());
        set_lorentz_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
        F -= Fm;
        dF -= dFm;
      }
      break;
    case MirroringType::SameAsLineShape:
      {
        Complex Fm;
        ComplexVector dFm(dF.nelem());
        switch(line.GetLineShapeType()) {
          case LineShapeType::Doppler:
            set_doppler_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
            break;
          case LineShapeType::Lorentz:
            set_lorentz_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
            break;
          case LineShapeType::Voigt:
            set_voigt_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
            break;
          case LineShapeType::HTP:
            set_htp_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
            break;
          case LineShapeType::ByPressureBroadeningData:
            switch(line.PressureBroadening().Type()) {
              case PressureBroadeningData::PB_SD_AIR_VOLUME:
              case PressureBroadeningData::PB_HTP_AIR_VOLUME:
              case PressureBroadeningData::PB_PURELY_FOR_TESTING:
              case PressureBroadeningData::PB_SD_TEST_WATER:
                set_htp_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
                break;
              case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
              case PressureBroadeningData::PB_PLANETARY_BROADENING:
              case PressureBroadeningData::PB_AIR_BROADENING:
              case PressureBroadeningData::PB_VOIGT_TEST_WATER:
                set_voigt_from_level_line_data(Fm, dFm, f, -line.F(), -0.0, level_line_data, derivatives_data, derivatives_data_position);
                break;
              case PressureBroadeningData::PB_NONE:
                Fm = Complex(0.0, 0.0);
                dFm = Complex(0.0, 0.0);
            }
            break;
          case LineShapeType::End: 
            break;
        }
        F -= Fm;
        dF -= dFm;
      }
      break;
    case MirroringType::None:
    case MirroringType::End:
      break;
  }
  
  switch(line.GetLineNormalizationType()) {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling_from_level_data(F, dF, f, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineNormalizationType::VVH:
      apply_VVH_scaling_from_level_data(F, dF, f, T, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineNormalizationType::VVW:
      apply_VVW_scaling_from_level_data(F, dF, f, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LineNormalizationType::End:
      break;
  }
  
  apply_linemixing_from_level_data(F, dF, level_line_data, derivatives_data, derivatives_data_position);
  
  apply_pressurebroadening_jacobian_scaling_from_level_data(dF, level_line_data, derivatives_data, derivatives_data_position);
  
  apply_linemixing_jacobian_scaling_from_level_data(dF, level_line_data, derivatives_data, derivatives_data_position);
  
  // Line population is last since it generates source variables with NLTE on
  switch(line.GetLinePopulationType()) {
    case LinePopulationType::ByLTE:
      apply_LTE_linestrength_from_level_data(F, dF, level_line_data, line, derivatives_data, derivatives_data_position);
      break;
    case LinePopulationType::ByVibrationalTemperatures:
      apply_NLTE_vibrational_temperature_linestrength_from_level_data(F, N, dF, dN, level_line_data, line, derivatives_data, derivatives_data_position);
      break;
    case LinePopulationType::ByPopulationDistribution:
      apply_NLTE_population_distribution_linestrength_from_level_data(F, N, dF, dN, level_line_data, derivatives_data, derivatives_data_position);
      break;
    case LinePopulationType::End:
      break;
  }
}


void Linefunctions::set_doppler_from_level_line_data(Complex& F,
                                                     ComplexVectorView dF,
                                                     const Numeric& f,
                                                     const Numeric& f0,
                                                     const Numeric& dZdH,
                                                     const SingleLevelLineData& level_line_data,
                                                     const ArrayOfRetrievalQuantity& derivatives_data,
                                                     const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  const Complex x = level_line_data.doppler_z(f0, f);
  const Complex w = Faddeeva::w(x);
  F = level_line_data.scale_w(w);
  
  if(nd) {
    const Complex dw_over_dx = level_line_data.dw_over_dz(x, w);
    for(Index id=0; id<nd; id++) {
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[id]]))
        dF[id] = level_line_data.invGD() * dw_over_dx;
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature)
        dF[id] = level_line_data.dGDdT() * (F + x * dw_over_dx) * (-level_line_data.invGD());
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
          dF[id] = (- F * level_line_data.GD_div_F0() 
                    + dw_over_dx * (x * level_line_data.GD_div_F0() - 1.0)) * level_line_data.invGD();
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[id]]))
        dF[id] = dw_over_dx * (-dZdH);
    }
  }
}


void Linefunctions::set_lorentz_from_level_line_data(Complex& F,
                                                     ComplexVectorView dF,
                                                     const Numeric& f,
                                                     const Numeric& f0,
                                                     const Numeric& dZdH,
                                                     const SingleLevelLineData& level_line_data,
                                                     const ArrayOfRetrievalQuantity& derivatives_data,
                                                     const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  const Complex z = level_line_data.lorentz_z(f0, f);
  const Complex L = invPI / z;
  
  // nb.  ARTS believes imag == refraction and real == absorption.
  // We must therefore always change order of arguments because that's not how Lorentz thinks
  F.real(L.imag());
  F.imag(L.real());
  
  if(nd) {
    const Complex dL2_over_dz = - L * L * PI;
    for(Index id=0; id<nd; id++) {
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[id]])) {
        dF[id].real(dL2_over_dz.imag());
        dF[id].imag(dL2_over_dz.real());
      }
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature) {
        const Complex d = dL2_over_dz * level_line_data.lorentz_dT();
        dF[id].real(d.imag());
        dF[id].imag(d.real());
      }
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
          dF[id].real(-dL2_over_dz.imag());
          dF[id].imag(-dL2_over_dz.real());
        }
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
          dF[id].real(dL2_over_dz.imag());
          dF[id].imag(dL2_over_dz.real());
        }
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[id]])) {
        dF[id].real(dL2_over_dz.imag() * (-dZdH));
        dF[id].imag(dL2_over_dz.real() * (-dZdH));
      }
      else if(is_line_mixing_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
          dF[id].real(dL2_over_dz.imag());
          dF[id].imag(dL2_over_dz.real());
        }
      }
    }
  }
}


void Linefunctions::set_voigt_from_level_line_data(Complex& F,
                                                   ComplexVectorView dF,
                                                   const Numeric& f,
                                                   const Numeric& f0,
                                                   const Numeric& dZdH,
                                                   const SingleLevelLineData& level_line_data,
                                                   const ArrayOfRetrievalQuantity& derivatives_data,
                                                   const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
 
  const Complex z = level_line_data.voigt_z(f0, f);
  const Complex w = Faddeeva::w(z);
  F = level_line_data.scale_w(w);
  
  if(nd) {
    const Complex dw_over_dz = level_line_data.dw_over_dz(z, w);
    for(Index id=0; id<nd; id++) {
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[id]]))
        dF[id] = dw_over_dz * level_line_data.invGD();
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature)
        dF[id] = - F * level_line_data.dGDdT() * level_line_data.invGD() + dw_over_dz * level_line_data.voigt_dT(z);
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
          dF[id] = (- F * level_line_data.dGDdT() * level_line_data.invGD()
                    + dw_over_dz * level_line_data.voigt_dF0(z)) * level_line_data.invGD();
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
          dF[id] = dw_over_dz * level_line_data.invGD();
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[id]]))
        dF[id] = - dw_over_dz * dZdH * level_line_data.invGD();
      else if(is_line_mixing_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
          dF[id] = dw_over_dz * level_line_data.invGD();
        }
      }
    }
  }
}


void Linefunctions::set_htp_from_level_line_data(Complex& F,
                                                 ComplexVectorView dF,
                                                 const Numeric& f,
                                                 const Numeric& f0,
                                                 const Numeric& dZdH,
                                                 const SingleLevelLineData& level_line_data,
                                                 const ArrayOfRetrievalQuantity& derivatives_data,
                                                 const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  const Complex x = level_line_data.htp_x(f0, f);
  const Numeric ratio_xy = level_line_data.xy_ratio(x);
  const Complex z1 = level_line_data.htp_z1(x, ratio_xy);
  const Complex z2 = level_line_data.htp_z2(x, z1, ratio_xy);
  const Complex w1 = level_line_data.htp_w1(z1);
  const Complex w2 = level_line_data.htp_w1(z2);
  const Complex A = level_line_data.htp_A(w1, w2, z1, ratio_xy);
  const Complex B = level_line_data.htp_B(w1, w2, z1, z2, ratio_xy);
  const Complex G = level_line_data.htp_G(A, B);
  F = A/G * invPI;
  
  if(nd) {
    const Complex dw1_over_dz1 = level_line_data.htp_dw1_over_dz1(z1, w1);
    const Complex dw2_over_dz2 = level_line_data.htp_dw2_over_dz2(z2, w2);
    
    for(Index id=0; id<nd; id++) {
      Complex dx;
      bool for_temperature = false;
      
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[id]]))
        dx = level_line_data.htp_dxdf();
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature) {
        dx = level_line_data.htp_dxdT(x);
        for_temperature = true;
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
          dx = level_line_data.htp_dxdC0();
        else 
          continue;
      }
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
          dx = level_line_data.htp_dxdf0();
        else
          continue;
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[id]]))
        dx = level_line_data.htp_dxdmag(dZdH);
      else
        continue; // to not repeat the code below!  FIXME:  must change this once pressure derivatives are possible
      
      const Complex dz1 = level_line_data.htp_dz1dt(z1, dx, ratio_xy, for_temperature);
      const Complex dz2 = level_line_data.htp_dz2dt(z2, dz1, dx, ratio_xy, for_temperature);
      const Complex dw1 = dz1 * dw1_over_dz1;
      const Complex dw2 = dz2 * dw2_over_dz2;
      const Complex dA = level_line_data.htp_dAdt(A, w1, dw1, dw2, z1, dz1, ratio_xy, for_temperature);
      const Complex dB = level_line_data.htp_dBdt(w1, dw1, w2, dw2, z1, dz1, z2, dz2, ratio_xy, for_temperature);
      const Complex dG = level_line_data.htp_dGdt(A, dA, dB, for_temperature);
      dF[id] = (invPI * dA - F * dG) / G;
    }
  }
}


void Linefunctions::apply_rosenkranz_quadratic_scaling_from_level_data(Complex& F,
                                                                       ComplexVectorView dF,
                                                                       const Numeric& f,
                                                                       const SingleLevelLineData& level_line_data,
                                                                       const ArrayOfRetrievalQuantity& derivatives_data,
                                                                       const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  const Numeric f2 = f * f;
  const Numeric scaling = f2 * level_line_data.norm();
  
  for(Index id=0; id<nd; id++) {
    if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature) 
      dF[id] = f2 * level_line_data.dnormdT() * F + dF[id] * scaling;
    else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
        dF[id] = f2 * level_line_data.dnormdf0() * F + dF[id] * scaling;
    }
    else if(is_frequency_parameter(derivatives_data[derivatives_data_position[id]]))
      dF[id] = 2.0 * f * level_line_data.norm() * F + dF[id] * scaling;
    else
      dF[id] *= scaling;
  }
  
  F *= scaling;
}


void Linefunctions::apply_VVH_scaling_from_level_data(Complex& F,
                                                      ComplexVectorView dF,
                                                      const Numeric& f,
                                                      const Numeric& T,
                                                      const SingleLevelLineData& level_line_data,
                                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                                      const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  const static Numeric c1 = 0.5 * PLANCK_CONST / BOLTZMAN_CONST;
  const Numeric tanh_part = tanh(c1 * f / T);
  const Numeric scaling = f * tanh_part * level_line_data.norm();
  
  if(nd) {
    const Numeric cosh_part = cosh(c1 * f / T);
    for(Index id=0; id<nd; id++) {
      if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature)
        dF[id] = (f * tanh_part * level_line_data.dnormdT() - c1 * f * f / (T * T * cosh_part * cosh_part) * level_line_data.norm()) * F + scaling * dF[id];
      else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
        if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
          dF[id] = f * tanh_part * level_line_data.dnormdf0() * F + scaling * dF[id];
      }
      else if(is_frequency_parameter(derivatives_data[derivatives_data_position[id]]))
        dF[id] = (f * c1 / T / (cosh_part * cosh_part) + tanh_part) * level_line_data.norm() * F + scaling * dF[id];
      else
        dF[id] *= scaling;
    }
  }
  
  F *= scaling;
}


void Linefunctions::apply_VVW_scaling_from_level_data(Complex& F,
                                                      ComplexVectorView dF,
                                                      const Numeric& f,
                                                      const SingleLevelLineData& level_line_data,
                                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                                      const ArrayOfIndex& derivatives_data_position) noexcept
{
  // Same code executed but the constant factor is different!  This is set elsewhere so call the same function
  apply_rosenkranz_quadratic_scaling_from_level_data(F, dF, f, level_line_data, derivatives_data, derivatives_data_position);
}


void Linefunctions::apply_linemixing_from_level_data(Complex& F,
                                                     ComplexVectorView dF,
                                                     const SingleLevelLineData& level_line_data,
                                                     const ArrayOfRetrievalQuantity& derivatives_data,
                                                     const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  for(Index id=0; id<nd; id++) {
    if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature)
      dF[id] = level_line_data.LM() * dF[id] + level_line_data.dLMdT() * F;
    else if(is_line_mixing_parameter(derivatives_data[derivatives_data_position[id]]))
      dF[id] = F;
    else
      dF[id] *= level_line_data.LM();
  }
  
  F *= level_line_data.LM();
}


void Linefunctions::apply_pressurebroadening_jacobian_scaling_from_level_data(ComplexVectorView dF,
                                                                              const SingleLevelLineData& level_line_data,
                                                                              const ArrayOfRetrievalQuantity& derivatives_data,
                                                                              const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  Index ig = 0;
  
  for(Index id=0; id<nd; id++) {
    if(level_line_data.no_more_pressure_jacs(ig)) break;
    if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[id]])) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
        dF[id] *= level_line_data.dgamma(ig);
        ++ig;
      }
    }
  }
}


void Linefunctions::apply_linemixing_jacobian_scaling_from_level_data(ComplexVectorView dF,
                                                                      const SingleLevelLineData& level_line_data,
                                                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                                                      const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  Index ilm = 0;
  
  for(Index id=0; id<nd; id++) {
    if(level_line_data.no_more_linemixing_jacs(ilm)) break;
    if(is_line_mixing_parameter(derivatives_data[derivatives_data_position[id]])) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
        dF[id] *= level_line_data.dlm(ilm);
        ++ilm;
      }
    }
  }
}


void Linefunctions::apply_LTE_linestrength_from_level_data(Complex& F, 
                                                           ComplexVectorView dF, 
                                                           const SingleLevelLineData& level_line_data, 
                                                           const LineRecord& line,
                                                           const ArrayOfRetrievalQuantity& derivatives_data, 
                                                           const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  for(Index id=0; id<nd; id++) {
    if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature)
      dF[id] = F * level_line_data.dSdT() + level_line_data.S() * dF[id];
    else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineStrength) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
        dF[id] = F * level_line_data.S() / line.I0();
    }
    else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition)
        dF[id] = F * level_line_data.dSdf0() + level_line_data.S() * dF[id];
    }
    else
      dF[id] *= level_line_data.S();
  }
  
  F *= level_line_data.S();
}


void Linefunctions::apply_NLTE_vibrational_temperature_linestrength_from_level_data(Complex& F, 
                                                                                    Complex& N,
                                                                                    ComplexVectorView dF, 
                                                                                    ComplexVectorView dN, 
                                                                                    const SingleLevelLineData& level_line_data,
                                                                                    const LineRecord& line,
                                                                                    const ArrayOfRetrievalQuantity& derivatives_data,
                                                                                    const ArrayOfIndex& derivatives_data_position) noexcept
{
  // First apply LTE variables for vib-temps
  apply_LTE_linestrength_from_level_data(F, dF, level_line_data, line, derivatives_data, derivatives_data_position);
  
  // Same code path but the constants are derived differently
  apply_NLTE_population_distribution_linestrength_from_level_data(F, N, dF, dN, level_line_data, derivatives_data, derivatives_data_position);
}


void Linefunctions::apply_NLTE_population_distribution_linestrength_from_level_data(Complex& F,
                                                                                    Complex& N,
                                                                                    ComplexVectorView dF,
                                                                                    ComplexVectorView dN,
                                                                                    const SingleLevelLineData& level_line_data,
                                                                                    const ArrayOfRetrievalQuantity& derivatives_data,
                                                                                    const ArrayOfIndex& derivatives_data_position) noexcept
{
  const Index nd = dF.nelem();
  
  for(Index id=0; id<nd; id++) {
    if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::Temperature) {
      dN[id] = F * level_line_data.dnlte_source_factordT() + dF[id] * level_line_data.nlte_source_factor();
      dF[id] = F * level_line_data.dnlte_absorption_factordT() + dF[id] * level_line_data.nlte_absorption_factor();
    }
    else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[id]])) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::FullTransition) {
        dN[id] = F * level_line_data.dnlte_source_factordf0() + dF[id] * level_line_data.nlte_source_factor();
        dF[id] = F * level_line_data.dnlte_absorption_factordf0() + dF[id] * level_line_data.nlte_absorption_factor();
      }
    }
    else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::NLTE) {
      if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::Lower or level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::Both) {
        dN[id] = F * level_line_data.dnlte_source_factordlow() + dF[id] * level_line_data.nlte_source_factor();
        dF[id] = F * level_line_data.dnlte_absorption_factordlow() + dF[id] * level_line_data.nlte_absorption_factor();
      }
      else if(level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::Upper or level_line_data(id) == SingleLevelLineData::SpectroscopyDerivivatives::Both) {
        dN[id] = F * level_line_data.dnlte_source_factordupp() + dF[id] * level_line_data.nlte_source_factor();
        dF[id] = F * level_line_data.dnlte_absorption_factordupp() + dF[id] * level_line_data.nlte_absorption_factor();
      }
    }
    else {
      dN[id] = dF[id] * level_line_data.nlte_source_factor();
      dF[id] *= level_line_data.nlte_absorption_factor();
    }
  }
  
  N = F * level_line_data.nlte_source_factor();
  F *= level_line_data.nlte_absorption_factor();
}



inline Complex Linefunctions::SingleLevelLineData::scale_w(const Complex& w) const noexcept
{
  return w * minvGD * sqrtInvPI;
}


inline Complex Linefunctions::SingleLevelLineData::dw_over_dz(const Complex& z, const Complex& w) const noexcept
{
  const static Complex s = Complex(0.0, sqrtInvPI);
  return 2.0 * minvGD * sqrtInvPI * (s - z * w);
}


inline Complex Linefunctions::SingleLevelLineData::htp_x(const Numeric& f0, const Numeric& f) const noexcept
{
  if(mC2t_is_zero)
    return (Complex(0.0, 0.0 + f0 + mZ - f) + mC0t) * minvGD;
  else 
    return (Complex(0.0, 0.0 + f0 + mZ - f) + mC0t) / mC2t;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dxdT(const Complex& x) const noexcept
{
  if(mC2t_is_zero)
    return (mdC0tdT - x * mdGDdT) * minvGD;
  else 
    return (mdC0tdT - x * mdC2tdT) / mC2t;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dxdf() const noexcept
{
  if(mC2t_is_zero)
    return Complex(0.0, -1.0) * minvGD;
  else 
    return Complex(0.0, - 1.0) / mC2t;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dxdf0() const noexcept
{
  return -htp_dxdf();
}



inline Complex Linefunctions::SingleLevelLineData::htp_dxdmag(const Numeric& zeeman_df) const noexcept
{
  return htp_dxdf0() * zeeman_df;
}



inline Complex Linefunctions::SingleLevelLineData::htp_dxdC0() const noexcept
{ 
  if(mC2t_is_zero)
    return (1 - meta) * minvGD;
  else 
    return (1 - meta) / mC2t;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dxdC2(const Complex& x) const noexcept
{
  if(mC2t_is_zero)
    return Complex(0.0, 0.0);// FIXME:  Find analytical solution for this...  maybe "- x * minvGD * (1 - meta)"?
  else 
    return - x / mC2t * (1 - meta);
}


inline Complex Linefunctions::SingleLevelLineData::htp_dxdFVC() const noexcept
{
  if(mC2t_is_zero)
    return 1.0 * minvGD;
  else 
    return 1.0 / mC2t;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dxdeta(const Complex& x) const noexcept
{
  if(mC2t_is_zero)
    return - 1.0 * minvGD;
  else 
    return (x - mC0 + 1.5*mC2) / mC2t;
}


static const Numeric MAXIMUM_XY_RATIO = 10e15;
static const Numeric MINIMUM_XY_RATIO = 3e-8;


inline Complex Linefunctions::SingleLevelLineData::htp_z1(const Complex& x, const Numeric& ratio_xy) const noexcept
{
  if(mC2t_is_zero)
    return x;
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot
    return sqrt(x);
  else if(ratio_xy < MINIMUM_XY_RATIO)  //  Y more than X by a lot
    return x * mC2t * minvGD;
  else
    return sqrt(x + my) - msqrty;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dz1dt(const Complex& z1, const Complex& dx, const Numeric& ratio_xy, const bool for_temperature) const noexcept
{
  if(mC2t_is_zero)
    return dx;
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot
    return dx / (2. * z1);
  else if(ratio_xy < MINIMUM_XY_RATIO and for_temperature)  //  Y more than X by a lot
    return dx * mC2t * minvGD + z1 * (mdC2dT / mC2t - minvGD * mdGDdT);
  else if(ratio_xy < MINIMUM_XY_RATIO)  //  Y more than X by a lot
    return dx * mC2t * minvGD;
  else if(for_temperature)
    return  (dx + mdydT) / (2.0 * (z1 + msqrty)) - mdydT / (2.0 * msqrty);
  else
    return  dx / (2.0 * (z1 + msqrty));
}


inline Complex Linefunctions::SingleLevelLineData::htp_z2(const Complex& x, const Complex& z1, const Numeric& ratio_xy) const noexcept
{
  if(mC2t_is_zero)  // z2 infinite, nb sqrt(z → inf) * w(z → inf) → 1 / sqrt(pi) & w(z → inf) → 0.
    return Complex(1e99, 1e99);
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot
    return sqrt(x + my);
  else if(ratio_xy < MINIMUM_XY_RATIO)  //  Y more than X by a lot
    return sqrt(x + my) + msqrty;
  else
    return z1 + 2.0 * msqrty;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dz2dt(const Complex& z2, const Complex& dz1, const Complex& dx, const Numeric& ratio_xy, const bool for_temperature) const noexcept
{
  if(mC2t_is_zero)  // z2 infinite, nb sqrt(z → inf) * w(z → inf) → 1 / sqrt(pi) & w(z → inf) → 0.
    return Complex(0.0, 0.0);
  else if(ratio_xy > MAXIMUM_XY_RATIO and for_temperature)  //  Y less than X by a lot
    return (dx + mdydT) / (2.0 * z2);
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot
    return dx / (2.0 * z2);
  else if(ratio_xy < MINIMUM_XY_RATIO and for_temperature)  //  Y more than X by a lot
    return  (dx + mdydT) / (2.0 * (z2 - msqrty)) + mdydT / (2.0 * msqrty);
  else if(ratio_xy < MINIMUM_XY_RATIO)  //  Y more than X by a lot
    return dx / (2.0 * (z2 - msqrty));
  else if(for_temperature)
    return dz1 + mdydT / msqrty;
  else
    return dz1;
}


inline Complex Linefunctions::SingleLevelLineData::htp_w1(const Complex& z1) const noexcept
{
  return Faddeeva::w(Complex(0.0, 1.0) * z1);  // z1 is set to always give results
}


inline Complex Linefunctions::SingleLevelLineData::htp_dw1_over_dz1(const Complex& z1, const Complex& w1) const noexcept
{
  return 2.0 * (z1 * w1 - sqrtInvPI);
}


inline Complex Linefunctions::SingleLevelLineData::htp_w2(const Complex& z2) const noexcept
{
  if(mC2t_is_zero)  // w(z → inf) → 0
    return Complex(0.0, 0.0);
  else
    return Faddeeva::w(Complex(0.0, 1.0) * z2);
}


inline Complex Linefunctions::SingleLevelLineData::htp_dw2_over_dz2(const Complex& z2, const Complex& w2) const noexcept
{
  if(mC2t_is_zero)  // w(z → inf) → 0?
    return Complex(0.0, 0.0);
  else
    return 2.0 * (z2 * w2 - sqrtInvPI);
}


inline Complex Linefunctions::SingleLevelLineData::htp_A(const Complex& w1, const Complex& w2, const Complex& z1, const Numeric& ratio_xy) const noexcept
{
  if(mC2t_is_zero)  // Standard Voigt function!
    return sqrtPI * minvGD * w1;
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot;  z1 is already sqrt(x);  Note problem for large x!
    return 2 * sqrtPI / mC2t * (sqrtInvPI - z1 * w1);
  else  // Diff of two Voigt functions
    return sqrtPI * minvGD * (w1 - w2);
}


inline Complex Linefunctions::SingleLevelLineData::htp_dAdt(const Complex& A,  const Complex& w1, const Complex& dw1, const Complex& dw2,
                                                            const Complex& z1, const Complex& dz1, const Numeric& ratio_xy, const bool for_temperature) const noexcept
{
  // FIXME:  Add derivatives for C2 and for eta here in future!
  
  if(mC2t_is_zero and for_temperature)  // Standard Voigt function!
    return minvGD * (sqrtPI * dw1 - A * mdGDdT);
  if(mC2t_is_zero)  // Standard Voigt function!
    return sqrtPI * minvGD * dw1;
  else if(ratio_xy > MAXIMUM_XY_RATIO and for_temperature)  //  Y less than X by a lot;  z1 is already sqrt(x);  Note problem for large x!
    return (- A * mdC2tdT + 2 * sqrtPI * (- dz1 * w1 - z1 * dw1)) / mC2t;
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot;  z1 is already sqrt(x);  Note problem for large x!
    return 2 * sqrtPI / mC2t * (- dz1 * w1 - z1 * dw1);
  else if(for_temperature)
    return (- A * mdGDdT + sqrtPI * (dw1 - dw2)) * minvGD;
  else  
    return sqrtPI * minvGD * (dw1 - dw2);
}


inline Complex Linefunctions::SingleLevelLineData::htp_B(const Complex& w1, const Complex& w2, const Complex& z1, const Complex& z2, const Numeric& ratio_xy) const noexcept
{
  if(mC2t_is_zero)  // In original equation, this term is unimportant because eta * C2 is 0.
    return Complex(0.0, 0.0);
  else if(ratio_xy > MAXIMUM_XY_RATIO)  //  Y less than X by a lot;  z1 is sqrt(x), z2 is sqrt(x + y);  Note problem for large z1/z2!
    return meta / (1 - meta) * (-1. + 2 * sqrtPI * ((1. - z2*z2 - my) * (sqrtInvPI - z1 * w1) +  z2 * w2));
  else
    return meta / (1 - meta) * (-1. + sqrtPI / (2. * msqrty) * ((1. - z1*z1) * w1 - (1. - z2*z2) * w2));
}


inline Complex Linefunctions::SingleLevelLineData::htp_dBdt(const Complex& w1, const Complex& dw1,
                                                            const Complex& w2, const Complex& dw2,
                                                            const Complex& z1, const Complex& dz1,
                                                            const Complex& z2, const Complex& dz2,
                                                            const Numeric& ratio_xy, const bool for_temperature) const noexcept
{
  if(mC2t_is_zero)  // In original equation, this term is unimportant because eta * C2 is 0
    return Complex(0.0, 0.0);
  else if(ratio_xy > MAXIMUM_XY_RATIO and for_temperature)  //  Y less than X by a lot;  z1 is sqrt(x), z2 is sqrt(x + y);  Note problem for large z1/z2!
    return meta / (1 - meta) * (2 * sqrtPI * ((- 2.0*dz2*z2 - mdydT) * (sqrtInvPI - z1 * w1) + (1. - z2*z2 - my) * (- dz1 * w1 - z1 * dw1) + dz2 * w2 + z2 * dw2));
  else if(ratio_xy > MAXIMUM_XY_RATIO and for_temperature)  //  Y less than X by a lot;  z1 is sqrt(x), z2 is sqrt(x + y);  Note problem for large z1/z2!
    return meta / (1 - meta) * (2 * sqrtPI * ((- 2.0*dz2*z2) * (sqrtInvPI - z1 * w1) + (1. - z2*z2 - my) * (- dz1 * w1 - z1 * dw1) + dz2 * w2 + z2 * dw2));
  else if(for_temperature)
    return meta / (1 - meta) * (sqrtPI / (2. * msqrty) * ((- 2.0*dz1*z1) * w1 + (1. - z1*z1) * dw1 - (- 2.0*dz2*z2) * w2 - (1. - z2*z2) * dw2) - mdydT * sqrtPI / (4. * my * msqrty) * ((1. - z1*z1) * w1 - (1. - z2*z2) * w2));
  else
    return meta / (1 - meta) * (sqrtPI / (2. * msqrty) * ((- 2.0*dz1*z1) * w1 + (1. - z1*z1) * dw1 - (- 2.0*dz2*z2) * w2 - (1. - z2*z2) * dw2));
    
}


inline Complex Linefunctions::SingleLevelLineData::htp_G(const Complex& A, const Complex& B) const noexcept
{
  return 1.0 - (mFVC - meta * (mC0 - 1.5 * mC2)) * A + B;
}


inline Complex Linefunctions::SingleLevelLineData::htp_dGdt(const Complex& A, const Complex& dA, const Complex& dB, const bool for_temperature) const noexcept
{
  // FIXME:  Add derivatives for C2 and for eta here in future!
  if(for_temperature)
    return - (mdFVCdT - meta * (mdC0dT - 1.5 * mdC2dT)) * A - (mFVC - meta * (mC0 - 1.5 * mC2)) * dA + dB;
  else
    return - (mFVC - meta * (mC0 - 1.5 * mC2)) * dA + dB;
}


Linefunctions::SingleLevelLineData::SingleLevelLineData(const LineRecord& line,
                                                        const ConstVectorView vmrs,
                                                        const ConstVectorView nlte_distribution,
                                                        const Numeric& T,
                                                        const Numeric& P,
                                                        const Numeric&,
                                                        const Numeric& lm_p_lim,
                                                        const Numeric& isotopic_ratio,
                                                        const Numeric& QT,
                                                        const Numeric& dQTdT,
                                                        const Numeric& QT0,
                                                        const ArrayOfIndex& broadening_species,
                                                        const Index this_species,
                                                        const Index water_species,
                                                        const Index,
                                                        const ArrayOfRetrievalQuantity& derivatives_data,
                                                        const ArrayOfIndex& derivatives_data_position)
{
  // Extract the quantum identify of the line to be used in-case there are derivatives
  const QuantumIdentifier& QI = line.QuantumIdentity();
  
  Numeric A, B, C, D, E;
  line.PressureBroadening().GetPressureBroadeningParams(
     A, B, meta, C, D, mFVC, T, line.Ti0(), P, P*vmrs[this_species], 
     this_species, water_species, broadening_species, vmrs);
  mC0 = Complex(A, C); mC2 = Complex(B, D);
  
  line.LineMixing().GetLineMixingParams(A, B, mDV, T, P, lm_p_lim);
  mLM = Complex(1.0 + B, A);
  
  if(do_temperature_jacobian(derivatives_data)) {
    // Pressure broadening partial derivatives
    line.PressureBroadening().GetPressureBroadeningParams_dT(
      A, B, E, C, D, mdFVCdT, T, line.Ti0(), P, P*vmrs[this_species],
      this_species, water_species, broadening_species, vmrs);
    mdC0dT.imag(C); mdC0dT.real(A); mdC2dT.imag(D); mdC2dT.real(B);
    
    line.LineMixing().GetLineMixingParams_dT(A, B, mdDVdT, T, temperature_perturbation(derivatives_data), P, lm_p_lim);
    mdLMdT = Complex(B, A);
  }
  
  const Index nd = derivatives_data_position.nelem();
  mspectroscopy_derivatives = Array<SpectroscopyDerivivatives>(nd, SpectroscopyDerivivatives::None);
  for(Index id=0; id<nd; id++) {
    if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[id]]) or
       is_line_mixing_parameter(derivatives_data[derivatives_data_position[id]])         or
       derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineCenter     or
       derivatives_data[derivatives_data_position[id]] == JacPropMatType::LineStrength) {
      if(QI > derivatives_data[derivatives_data_position[id]].QuantumIdentity()) {
        mspectroscopy_derivatives[id] = SpectroscopyDerivivatives::FullTransition;
      }
    }
    else if(derivatives_data[derivatives_data_position[id]] == JacPropMatType::NLTE) {
      if(QI.LowerQuantumId() > derivatives_data[derivatives_data_position[id]].QuantumIdentity()) {
        mspectroscopy_derivatives[id] = SpectroscopyDerivivatives::Lower;
      }
      if(QI.UpperQuantumId() > derivatives_data[derivatives_data_position[id]].QuantumIdentity()) {
        if(mspectroscopy_derivatives[id] == SpectroscopyDerivivatives::Lower)
          mspectroscopy_derivatives[id] = SpectroscopyDerivivatives::Both;
        else
          mspectroscopy_derivatives[id] = SpectroscopyDerivivatives::Upper;
      }
    }
  }
  
  line.PressureBroadening().SetInternalDerivatives(
    mpressure_derivatives, derivatives_data, QI, line.Ti0()/T, P, P*vmrs[this_species], 
    this_species, water_species, vmrs);
  
  line.LineMixing().SetInternalDerivatives(mlinemixing_derivatives, derivatives_data, QI, T, P, lm_p_lim);
  
  mZ = 0.0;
  
  switch(line.GetLineNormalizationType()) {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::VVH:
      A = 2*BOLTZMAN_CONST*T;
      B = tanh(PLANCK_CONST*line.F()/A);
      mnorm = 1 / (line.F()*B);
      mdnormdT = PLANCK_CONST*(1-B*B)/(A*A*B*B);
      mdnormdf0 = - (mnorm/line.F() + PLANCK_CONST*(1-B*B)/(A*line.F()*B*B));
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      A = 2*BOLTZMAN_CONST*T;
      B = sinh(PLANCK_CONST*line.F()/A);
      C = cosh(PLANCK_CONST*line.F()/A);
      mnorm = PLANCK_CONST/(A*line.F()*B);
      mdnormdT = -mnorm/T + PLANCK_CONST*PLANCK_CONST*C/(A*A*T*B*B);
      mdnormdf0 = -mnorm/line.F() - PLANCK_CONST*PLANCK_CONST*C/(A*A*line.F()*B*B);
      break;
    case LineNormalizationType::VVW:
      A = line.F() * line.F();
      mnorm = 1 / A;
      mdnormdT = 0.0;
      mdnormdf0 = - 2 * mnorm / line.F();
      break;
    case LineNormalizationType::End:
      throw std::runtime_error("Bad line normalization");
  }
  
  switch(line.GetLinePopulationType()) {
    case LinePopulationType::ByLTE:
    case LinePopulationType::ByVibrationalTemperatures:
    {
      // Line strength scaling that are line-dependent ---
      const Numeric gamma = stimulated_emission(T, line.F());
      const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
      const Numeric K1 = boltzman_ratio(T, line.Ti0(), line.Elow());
      const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
      mS = K1 * K2 * QT0/QT * isotopic_ratio * line.I0();
      if(do_temperature_jacobian(derivatives_data)) {
        const Numeric dK1_dT = dboltzman_ratio_dT(K1, T, line.Elow());
        const Numeric dK2_dT = dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T);
        mdSdT = line.I0() * isotopic_ratio * QT0/QT  * (K1 * dK2_dT + dK1_dT * K2) - mS / QT * dQTdT;
      }
      mdSdf0 = line.I0() * isotopic_ratio * QT0/QT * K1 * dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0());
      if(line.GetLinePopulationType() != LinePopulationType::ByLTE) {
        const Numeric K4 = boltzman_ratio(nlte_distribution[line.NLTEUpperIndex()], T, line.Evupp());
        const Numeric rlow = boltzman_ratio(nlte_distribution[line.NLTELowerIndex()], T, line.Evlow());
        const Numeric K3 = absorption_nlte_ratio(gamma, K4, rlow);
        const Numeric dK3_dF0 = dabsorption_nlte_rate_dF0(gamma, T, K4, rlow);
        const Numeric dK3_dT = dabsorption_nlte_rate_dT(gamma, T, line.F(), line.Evlow(), line.Evupp(), K4, rlow);
        const Numeric dK3_dTl = dabsorption_nlte_rate_dTl(gamma, T, nlte_distribution[line.NLTELowerIndex()], line.Evlow(), rlow);
        const Numeric dK3_dTu = dabsorption_nlte_rate_dTu(gamma, T, nlte_distribution[line.NLTEUpperIndex()], line.Evupp(), K4);
        const Numeric dK4_dTu = dboltzman_ratio_dT(K4, nlte_distribution[line.NLTEUpperIndex()], line.Evupp());
        
        // NOTE:  This code might have bit-rot in it...  This is the third copy...
        mnlte_abs = K3;
        mdnlte_absdT = dK3_dT;
        mdnlte_absdf0 = dK3_dF0;
        mdnlte_absdlow = dK3_dTl;
        mdnlte_absdupp = dK3_dTu;
        mnlte_src = K4/K3 - 1.0;
        mdnlte_absdf0 = - dK3_dF0 / K3 / K3;
        mdnlte_absdT = - dK3_dT / K3 / K3;
        mdnlte_srcdlow = - dK3_dTl / K3 / K3;
        mdnlte_srcdupp = (dK4_dTu - dK3_dTu / K3) / K3;
      }
      break;
    }
    case LinePopulationType::ByPopulationDistribution:
    {
      // Physical constants
      const static Numeric c0 = 2.0 * PLANCK_CONST / SPEED_OF_LIGHT / SPEED_OF_LIGHT;
      const static Numeric c1 = PLANCK_CONST / 4 / PI;
      const static Numeric fac = c1 / c0;
      const Numeric x = line.G_upper() / line.G_lower();
      const Numeric f = fac * line.A() / (line.F() * line.F());
      const Numeric exp_T = exp(PLANCK_CONST * line.F() / BOLTZMAN_CONST / T);
      
      const Numeric& r1 = nlte_distribution[line.NLTELowerIndex()];
      const Numeric& r2 = nlte_distribution[line.NLTEUpperIndex()];
      
      mdnlte_absdlow = -f;
      mdnlte_absdupp = x * f;
      mnlte_abs = r1 * mdnlte_absdlow + r2 * mdnlte_absdupp;
      mdnlte_absdT = 0.0;
      mdnlte_absdf0 = -2.0*mnlte_abs/line.F();
      
      
      mdnlte_srcdupp = (exp_T - 1) * f - mdnlte_absdupp;
      mdnlte_srcdlow = -mdnlte_absdlow;
      mnlte_src = r1 * mdnlte_srcdlow + r2 * mdnlte_srcdupp;
      mdnlte_absdf0 = -2.0*mnlte_src/line.F() + PLANCK_CONST / BOLTZMAN_CONST / T * f;
      mdnlte_srcdT = - r2 * PLANCK_CONST * line.F() / BOLTZMAN_CONST * exp_T * f / T / T;
      break;
    }
    case LinePopulationType::End:
      throw std::runtime_error("Error in line populations");
  }
  
  mGD_div_F0 = DopplerConstant(T, line.IsotopologueData().Mass());
  minvGD = 1 / (mGD_div_F0 * line.F());
  mdGDdT = dDopplerConstant_dT(T, line.IsotopologueData().Mass()) * line.F();
}


ostream& Linefunctions::operator<<(ostream& os, const Linefunctions::SingleLevelLineData& slld)
{
  // For internal review of values only
  
  os << "Pressure Data:\n\t";
  os << "C0: " << slld.mC0<<", ";
  os << "C2: " << slld.mC2<<", ";
  os << "FVC: " << slld.mFVC<<", ";
  os << "eta: " << slld.meta<<", ";
  os << "\n\t";
  os << "dC0/dT: " << slld.mdC0dT<<", ";
  os << "dC2/dT: " << slld.mdC2dT<<", ";
  os << "dFVC/dT: " << slld.mdFVCdT<<", ";
  os << "\n";
  
  os << "Pressure Data Spectroscopic Derivivatives:\n\t";
  for(const auto& c : slld.mpressure_derivatives)
    os << "Value: " << c << ", ";
  os << "\n";
  
  os << "Line Mixing:\n\t";
  os << "LM: " << slld.mLM<<", ";
  os << "DV: " << slld.mDV<<", ";
  os << "dLM/dT: " << slld.mdLMdT<<", ";
  os << "dDV/dT: " << slld.mdDVdT<<", ";
  os << "\n";
  
  os << "Line Mixing Spectroscopic Derivivatives:\n\t";
  for(const auto& c : slld.mlinemixing_derivatives)
    os << "Value: " << c << ", ";
  os << "\n";
  
  os << "Doppler Effects:\n\t";
  os << "GD/F0: " << slld.mGD_div_F0 << ", ";
  os << "1/GD: " << slld.minvGD<< ", ";
  os << "dGD/dT: " << slld.mdGDdT<< ", ";
  os << "\n";
  
  os << "Zeeman shift:\n\t";
  os << slld.mZ << ", ";
  os << "\n";
  
  os << "HTP Derived Parameters:\n\t";
  os << "C0-t: " << slld.mC0t << ", ";
  os << "C2-t: " << slld.mC2t << ", ";
  os << "Y: " << slld.my << ", ";
  os << "sqrt(Y): " << slld.msqrty << ", ";
  os << "\n\t";
  os << "dC0-t/dT: " << slld.mdC0tdT << ", ";
  os << "dC2-t/dT: " << slld.mdC2tdT << ", ";
  os << "dY/dT: " << slld.mdydT << ", ";
  os << "\n";
  
  os << "Line Shape Normalization Factors:\n\t";
  os << "Norm: " << slld.mnorm << ", ";
  os << "dNorm/dT: " << slld.mdnormdT << ", ";
  os << "dNorm/df0: " << slld.mdnormdf0 << ", ";
  os << "\n";
  
  os << "Line Strength  LTE Factors:\n\t";
  os << "S: " << slld.mS << ", ";
  os << "dS/dT: " << slld.mdSdT<< ", ";
  os << "dS/df0: " << slld.mdSdf0 << ", ";
  os << "\n";
  
  os << "Line Strength NLTE Factors:\n\t";
  os << "a: " << slld.mnlte_abs << ", ";
  os << "da/dT: " << slld.mdnlte_absdT << ", ";
  os << "da/df0: " << slld.mdnlte_absdf0 << ", ";
  os << "da/dn2: " << slld.mdnlte_absdupp << ", ";
  os << "da/dn1: " << slld.mdnlte_absdlow << ", ";
  os << "\n\t";
  os << "s: " << slld.mnlte_src << ", ";
  os << "ds/dT: " << slld.mdnlte_srcdT << ", ";
  os << "ds/df0: " << slld.mdnlte_srcdf0 << ", ";
  os << "ds/dn2: " << slld.mdnlte_srcdupp << ", ";
  os << "ds/dn1: " << slld.mdnlte_srcdlow << ", ";
  os << "\n";
  
  return os;
}

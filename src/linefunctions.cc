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
                                  const Numeric& zeeman_df, 
                                  const Numeric& magnetic_magnitude,
                                  const ArrayOfIndex& broad_spec_locations,
                                  const Index& this_species,
                                  const Index& water_species,
                                  const Verbosity& verbosity)
{
  // Pressure broadening terms
  const Numeric partial_pressure = pressure * vmrs[this_species];
  Numeric G0, G2, e, L0, L2, FVC;
  line.PressureBroadening().GetPressureBroadeningParams(
    G0, G2, e, L0, L2, FVC, line.Ti0()/temperature, pressure, partial_pressure, 
    this_species, water_species, broad_spec_locations, vmrs, verbosity);
  
  // Line mixing terms
  Numeric Y=0, G=0, DV=0;
  line.LineMixing().GetLineMixingParams(Y, G, DV, temperature, pressure,  pressure_limit_for_linemixing);
  
  // Line shape usage remembering variable
  LineShapeType lst = LineShapeType::End;
  
  ArrayOfComplexVector dF(0);
  const Numeric doppler_constant = DopplerConstant(temperature, line.IsotopologueData().Mass());
  
  switch(line.GetLineShapeType())
  {
    case LineShapeType::ByPressureBroadeningData:
      switch(line.PressureBroadening().Type())
      {
        // Use data as per speed dependent air
        case PressureBroadeningData::PB_SD_AIR_VOLUME:
        case PressureBroadeningData::PB_PURELY_FOR_TESTING:
          lst = LineShapeType::HTP;
          set_htp(F, dF, f_grid, zeeman_df, magnetic_magnitude, line.F(), doppler_constant, G0, L0, G2, L2, e, FVC);
          break;
          // Use for data that requires air and water Voigt broadening
        case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
          // Use for data that requires planetary Voigt broadening
        case PressureBroadeningData::PB_PLANETARY_BROADENING:
          // Use for data that requires air Voigt broadening
        case PressureBroadeningData::PB_AIR_BROADENING:
          // Above should be all methods of pressure broadening requiring Voigt in ARTS by default
          lst = LineShapeType::Voigt;
          set_faddeeva_algorithm916(F, dF, f_grid, zeeman_df, magnetic_magnitude, line.F(), doppler_constant, G0, L0, DV);
          break;
        default:
          throw std::runtime_error("Developer has messed up and needs to add the key to the code above this error");
      }
      break;
        case LineShapeType::Doppler:
          lst = LineShapeType::Doppler;
          set_doppler(F, dF, f_grid, zeeman_df, magnetic_magnitude, line.F(), doppler_constant);
          break;
          // This line only needs Hartmann-Tran
        case LineShapeType::HTP:
          set_htp(F, dF, f_grid, zeeman_df, magnetic_magnitude, line.F(), doppler_constant, G0, L0, G2, L2, e, FVC);
          lst = LineShapeType::HTP;
          break;
          // This line only needs Lorentz
        case LineShapeType::Lorentz:
          lst = LineShapeType::Lorentz;
          set_lorentz(F, dF, f_grid, zeeman_df, magnetic_magnitude, line.F(), G0, L0, DV);
          break;
          // This line only needs Voigt
        case LineShapeType::Voigt:
          lst = LineShapeType::Voigt;
          set_faddeeva_algorithm916(F, dF, f_grid, zeeman_df, magnetic_magnitude, line.F(), doppler_constant, G0, L0, DV);
          break;
        case LineShapeType::End:
        default:
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
      ArrayOfComplexVector dFm(0);
      
      set_lorentz(Fm, dFm, f_grid, -zeeman_df, magnetic_magnitude, -line.F(), G0, -L0, -DV);
      
      // Apply mirroring
      F -= Fm;
    }
    break;
    // Same type of mirroring as before
    case MirroringType::SameAsLineShape:
    {
      // Set the mirroring computational vectors and size them as needed
      ComplexVector Fm(F.nelem());
      ArrayOfComplexVector dFm(0);
      
      switch(lst)
      {
        case LineShapeType::Doppler:
          set_doppler(Fm, dFm, f_grid, -zeeman_df, magnetic_magnitude, -line.F(), -doppler_constant);
          break;
        case LineShapeType::Lorentz:
          set_lorentz(Fm, dFm, f_grid, -zeeman_df, magnetic_magnitude, -line.F(), G0, -L0, -DV);
          break;
        case LineShapeType::Voigt:
          set_faddeeva_algorithm916(Fm, dFm, f_grid, -zeeman_df, magnetic_magnitude, -line.F(), -doppler_constant, G0, -L0, -DV);
          break;
        case LineShapeType::HTP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dFm, f_grid, -zeeman_df, magnetic_magnitude, -line.F(), -doppler_constant, G0, -L0, G2, -L2, e, FVC);
          break;
        case LineShapeType::ByPressureBroadeningData:
        case LineShapeType::End:
        default:
          throw std::runtime_error("Cannot understand the requested line shape type for mirroring.");
      }
      F -= Fm;
      break;
    }
        case MirroringType::End:
        default:
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
    default:
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
 * \param quantum_identity ID of the absorption line
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
 * \param ddF0_dT Temperature derivative of dF0
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_lorentz(ComplexVectorView F,
                                ArrayOfComplexVector& dF,
                                ConstVectorView f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& G0,
                                const Numeric& L0, 
                                const Numeric& dF0,
                                const PropmatPartialsData& derivatives_data,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dG0_dT,
                                const Numeric& dL0_dT,
                                const Numeric& ddF0_dT,
                                const Range& df_range)
{ 
  // Size of the problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data.nelem();
  
  // The central frequency
  const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude + dF0;
  
  // Constant part of the denominator
  const Complex denom0 = Complex(G0, F0);
  
  Complex d, denom;
  
  for(Index iv = 0; iv < nf; iv++)
  {
    denom = 1.0 / ((denom0 - Complex(0.0, f_grid[iv])));
    
    F[iv] = invPI * denom;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq == 0)
        d = - F[iv] * denom;
      
      if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
        // Temperature derivative only depends on how pressure shift and broadening change
        dF[iq][df_range][iv] = d * Complex(dG0_dT, dL0_dT + ddF0_dT);
      }
      else if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      {
        // Frequency scale 1 to -1 linearly
        dF[iq][df_range][iv] = d * Complex(0.0, -1.0);
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          // Line center scales 1 to 1 linearly
          dF[iq][df_range][iv] = d * Complex(0.0, 1.0);
        }
      }
      else if(derivatives_data.IsPressureBroadeningParameter(derivatives_data(iq)))
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          // Pressure broadening will be dealt with in another function, though the partial derivative
          dF[iq][df_range][iv] = d;
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::MagneticMagnitude)
      {
        // Magnetic magnitude changes like line center in part
        // FIXME: Add magnetic components here
        dF[iq][df_range][iv] = d * Complex(0.0, zeeman_df);
      }
      else if(derivatives_data.IsLineMixingDFParameter(derivatives_data(iq)))
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][df_range][iv] = d * Complex(0.0, -1.0);
        }
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
                            const PropmatPartialsData& derivatives_data,
                            const QuantumIdentifier& quantum_identity,
                            const Numeric& dGD_div_F0_dT,
                            const Numeric& dG0_dT,
                            const Numeric& dL0_dT,
                            const Numeric& dG2_dT,
                            const Numeric& dL2_dT,
                            const Numeric& deta_dT,
                            const Numeric& dFVC_dT,
                            const Range& df_range)
{
  // Size of the problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data.nelem();
  
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
  const Index first_frequency = derivatives_data.get_first_frequency(), 
  first_pressure_broadening = derivatives_data.get_first_pressure_term();
  
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
  const Complex invC2t = C2t_zero_limit ? Complex(0.0, 0.0) : (1.0 / C2t);
  
  // Relative pressure broadening terms to be used in the shape function
  const Complex sqrtY =  C2t_zero_limit ? -1 : (0.5 * invC2t * GD);
  const Complex Y =      C2t_zero_limit ? -1 : (sqrtY * sqrtY);
  const Numeric invAbsY =   C2t_zero_limit ? -1.0 : 1.0 / abs(Y);
  
  // Temperature derivatives precomputed
  Complex dC0_dT, dC2_dT, dC0t_dT, dC2t_dT, dC0_m1p5_C2_dT, dY_dT;
  if(derivatives_data.do_temperature())
  {
    dC0_dT = dG0_dT + i*dL0_dT;
    dC2_dT = dG2_dT + i*dL2_dT;
    dC0t_dT = one_minus_eta * (dC0_dT - 1.5 * dC2_dT) - deta_dT * C0_m1p5_C2 + dFVC_dT;
    dC2t_dT = one_minus_eta * dC2_dT - deta_dT * C2;
    dC0_m1p5_C2_dT = dC0_dT - 1.5 * dC2_dT;
    dY_dT = C2t_zero_limit ? -1.0 : (GD / 2.0*invC2t*invC2t * (dGD_dT - GD * invC2t * dC2t_dT));
  }
  
  // Scale factor (normalizes to PI)
  const Numeric fac = sqrtPI * invGD;
  
  for(Index iv = 0; iv < nf; iv++)
  {
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
    if(C2t_zero_limit)
    {
      Zm2 = Zm * Zm;
    }
    else if(not ratioXY_high_limit)
    {
      Zm2 = Zm * Zm;
      Zp2 = Zp * Zp;
    }
    
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
      G = 1.0 - (FVC - eta*C0_m1p5_C2) * A + eta / one_minus_eta * (-1.0 + 2.0 * sqrtPI * (1.0-X-2.0*Y) * (sqrtInvPI - Zm * wiZm) + 2.0*sqrtPI*Zp*wiZp);
    else
      G = 1.0 - (FVC - eta*C0_m1p5_C2) * A + eta / one_minus_eta * (-1.0 + sqrtPI/(2.0*sqrtY) * ((1.0-Zm2)*wiZm - (1.0-Zp2)*wiZp));
   
    // Compute denominator
    invG = 1.0 / G;
    
    // Compute line shape
    F[iv] = A * invG * invPI;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      {
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
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
          
          dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dG); 
        }
        else  // copy for repeated occurrences
        {
          dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
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
          
        dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dG); 
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
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
          
          dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dG); 
        }
      }
      else if(derivatives_data.IsPressureBroadeningParameter(derivatives_data(iq))) 
      {
        // NOTE:  These are first order Voigt-like.  
        // The variables that are not Voigt-like must be dealt with separately
        
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          if(first_pressure_broadening == iq)
          {
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
            
            dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dG); 
          }
          else  // copy for repeated occurrences
          {
            dF[iq][df_range][iv] = dF[first_pressure_broadening][df_range][iv]; 
          }
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::MagneticMagnitude)
      {
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
        
        dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dG); 
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
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
 * \param dF0_dT Temperature derivative of dF0_dT
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_faddeeva_algorithm916(ComplexVectorView F, 
                                              ArrayOfComplexVector& dF, 
                                              ConstVectorView f_grid, 
                                              const Numeric& zeeman_df, 
                                              const Numeric& magnetic_magnitude,
                                              const Numeric& F0_noshift, 
                                              const Numeric& GD_div_F0,
                                              const Numeric& G0, 
                                              const Numeric& L0,
                                              const Numeric& dF0,
                                              const PropmatPartialsData& derivatives_data,
                                              const QuantumIdentifier& quantum_identity,
                                              const Numeric& dGD_div_F0_dT,
                                              const Numeric& dG0_dT,
                                              const Numeric& dL0_dT,
                                              const Numeric& dF0_dT,
                                              const Range& df_range)
{
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data.nelem();
  
  // For calculations
  Numeric dx;
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
  
  const Index first_pressure_broadening = derivatives_data.get_first_pressure_term(),
  first_frequency = derivatives_data.get_first_frequency();
  
  // frequency in units of Doppler
  for (Index iv=0; iv<nf; iv++)
  {
    dx = f_grid[iv] * invGD;
    z = z0 + dx;
    w = Faddeeva::w(z);
    
    F[iv] = fac * w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq==0)
        dw_over_dz = 2.0 * fac *  (Complex(0, sqrtInvPI) - z * w);
      
      if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      { 
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          //dz = Complex(invGD, 0.0);
          
          dF[iq][df_range][iv] = dw_over_dz * invGD; //dz; 
        }
        else  // copy for repeated occurrences
        {
          dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
        dz = (Complex(-dL0_dT - dF0_dT, dG0_dT) - z * dGD_dT) * invGD;
        
        dF[iq][df_range][iv] = -F[iv] * dGD_dT * invGD;
        dF[iq][df_range][iv] += dw_over_dz * dz;
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[iq][df_range][iv] = -F[iv] * GD_div_F0;
          dF[iq][df_range][iv] += dw_over_dz * dz;
          dF[iq][df_range][iv] *= invGD;
        }
      }
      else if(derivatives_data.IsPressureBroadeningParameter(derivatives_data(iq))) // Only the zeroth order terms --- the derivative with respect to these have to happen later
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == iq)
          {
            dz = Complex(0.0, 1.0) * invGD;
            dF[iq][df_range][iv] = dw_over_dz * dz; 
          }
          else
          {
            dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
          }
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::MagneticMagnitude)// No external inputs --- errors because of frequency shift when Zeeman is used?
      {
        // dz = Complex(- zeeman_df * invGD, 0.0);
        
        dF[iq][df_range][iv] = dw_over_dz * (- zeeman_df * invGD); //* dz; 
      }
      else if(derivatives_data.IsLineMixingDFParameter(derivatives_data(iq)))
      {
        // dz = Complex(-invGD, 0.0);
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][df_range][iv] = - (dw_over_dz * invGD); //* dz;
        }
      }
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
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_doppler(ComplexVectorView F, // Sets the full complex line shape without line mixing
                                ArrayOfComplexVector& dF,
                                ConstVectorView f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& GD_div_F0,
                                const PropmatPartialsData& derivatives_data,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dGD_div_F0_dT,
                                const Range& df_range)
{
  const Index nf = f_grid.nelem();
  const Index nppd = dF.nelem();
  
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
  const Index first_frequency = derivatives_data.get_first_frequency();
  
  for(Index iv = 0; iv < nf; iv++)
  {
    x = (f_grid[iv] - F0) * invGD;
    w = Faddeeva::w(x);
    
    F[iv] = fac * w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq == 0)
        dw_over_dx = 2.0 * fac * (Complex(0.0, sqrtInvPI) - x * w);
      
      if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      {
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          dF[iq][df_range][iv] = dw_over_dx * invGD;
        }
        else  // copy for repeated occurrences
        {
          dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
        dF[iq][df_range][iv] = F[iv] * dGD_dT + x * dGD_dT * dw_over_dx;
        dF[iq][df_range][iv] *= -invGD ;
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][df_range][iv] = -F[iv] * GD_div_F0;
          dF[iq][df_range][iv] += dw_over_dx * (-x * GD_div_F0 - 1.0);
          dF[iq][df_range][iv] *= invGD;
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
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dL0_dT Temperature derivative of L0
 * 
 */
void Linefunctions::set_faddeeva_from_full_linemixing(ComplexVectorView F, 
                                                      ArrayOfComplexVector& dF,
                                                      ConstVectorView f_grid,
                                                      const Complex& eigenvalue_no_shift,
                                                      const Numeric& GD_div_F0,
                                                      const Numeric& L0,
                                                      const PropmatPartialsData& derivatives_data,
                                                      const QuantumIdentifier& quantum_identity,
                                                      const Numeric& dGD_div_F0_dT,
                                                      const Complex& deigenvalue_dT,
                                                      const Numeric& dL0_dT)
{
  // For calculations
  Numeric dx;
  Complex w, z, dw_over_dz, dz;
  
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
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
  
  const Index first_pressure_broadening = derivatives_data.get_first_pressure_term();
  const Index first_frequency = derivatives_data.get_first_frequency();
  
  // frequency in units of Doppler
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
      
      if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      { 
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          //dz = Complex(invGD, 0.0);
          
          dF[iq][iv] = fac * dw_over_dz * invGD; //dz; 
        }
        else  // copy for repeated occurrences
        {
          dF[iq][iv] = dF[first_frequency][iv]; 
        }
      }
      else if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
        dz = (deigenvalue_dT - dL0_dT) - z * dGD_dT;
        
        dF[iq][iv] = -F[iv] * dGD_dT;
        dF[iq][iv] += fac * dw_over_dz * dz;
        dF[iq][iv] *= invGD;
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[iq][iv] = -F[iv] * GD_div_F0;
          dF[iq][iv] += dw_over_dz * dz;
          dF[iq][iv] *= fac * invGD;
        }
      }
      else if(derivatives_data.IsPressureBroadeningParameter(derivatives_data(iq))) // Only the zeroth order terms --- the derivative with respect to these have to happen later
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == iq)
          {
            dz = Complex(-1.0, 1.0) * invGD;
            dF[iq][iv] = fac * dw_over_dz * dz; 
          }
          else
          {
            dF[iq][iv] = dF[first_frequency][iv]; 
          }
        }
      }
      else if(derivatives_data.IsLineMixingDFParameter(derivatives_data(iq)))
      {
        // dz = Complex(-invGD, 0.0);
        
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][iv] = fac * dw_over_dz * (-invGD); //* dz;
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
 * \param quantum_identity ID of the absorption line
 * \param dY_dT Temperature derivative of Y
 * \param dG_dT Temperature derivative of G
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linemixing_scaling(ComplexVectorView F,
                                             ArrayOfComplexVector& dF,
                                             const Numeric& Y,
                                             const Numeric& G,
                                             const PropmatPartialsData& derivatives_data,
                                             const QuantumIdentifier& quantum_identity,
                                             const Numeric& dY_dT,
                                             const Numeric& dG_dT,
                                             const Range& df_range)
{
  const Index nf = F.nelem(), nppd = derivatives_data.nelem();
  
  const Complex LM = Complex(1.0 + G, -Y);
  const Complex dLM_dT = Complex(dG_dT, -dY_dT);
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JacPropMatType::Temperature)
    {
      dF[iq][df_range] *= LM;
      for(Index iv = 0; iv < nf; iv++)
        dF[iq][df_range][iv] += F[iv] * dLM_dT;
    }
    else if(derivatives_data.IsLineMixingParameter(derivatives_data(iq)))
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] = F;
    }
    else
      dF[iq][df_range] *= LM;
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
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_rosenkranz_quadratic_scaling(ComplexVectorView F,
                                                       ArrayOfComplexVector& dF,
                                                       ConstVectorView f_grid,
                                                       const Numeric& F0,
                                                       const Numeric& T,
                                                       const PropmatPartialsData& derivatives_data,
                                                       const QuantumIdentifier& quantum_identity,
                                                       const Range& df_range)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  const Numeric invF0 = 1.0/F0;
  const Numeric mafac = (PLANCK_CONST) / (2.0 * BOLTZMAN_CONST * T) /
  sinh((PLANCK_CONST * F0) / (2.0 * BOLTZMAN_CONST * T)) * invF0;
  
  Numeric dmafac_dF0_div_fun = 0, dmafac_dT_div_fun = 0;
  if(derivatives_data.do_line_center())
  {
    dmafac_dF0_div_fun = -invF0 - 
    PLANCK_CONST/(2.0*BOLTZMAN_CONST*T*tanh(F0*PLANCK_CONST/(2.0*BOLTZMAN_CONST*T)));
  }
  if(derivatives_data.do_temperature())
  {
    dmafac_dT_div_fun = -(BOLTZMAN_CONST*T - F0*PLANCK_CONST/
    (2.0*tanh(F0*PLANCK_CONST/(2.0*BOLTZMAN_CONST*T))))/(BOLTZMAN_CONST*T*T);
  }
  
  Numeric fun;
  
  for (Index iv=0; iv < nf; iv++)
  {
    fun = mafac * (f_grid[iv] * f_grid[iv]);
    F[iv] *= fun;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      dF[iq][df_range][iv] *= fun;
      if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
        dF[iq][df_range][iv] += dmafac_dT_div_fun * F[iv];
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][df_range][iv] += dmafac_dF0_div_fun * F[iv];
      }
      else if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      {
        dF[iq][df_range][iv] += (2.0 / f_grid[iv]) * F[iv];
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
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_VVH_scaling(ComplexVectorView F,
                                      ArrayOfComplexVector& dF,
                                      ConstVectorView f_grid,
                                      const Numeric& F0,
                                      const Numeric& T,
                                      const PropmatPartialsData& derivatives_data,
                                      const QuantumIdentifier& quantum_identity,
                                      const Range& df_range)
{ 
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // 2kT is constant for the loop
  const Numeric kT = 2.0 * BOLTZMAN_CONST * T;
  
  // denominator is constant for the loop
  const Numeric tanh_f0part = tanh(PLANCK_CONST * F0 / kT);
  const Numeric denom = F0 * tanh_f0part;
  
  Numeric fac_df0 = 0; 
  if(derivatives_data.do_line_center())
    fac_df0 = (-1.0/F0 + PLANCK_CONST*tanh_f0part/(kT) - PLANCK_CONST/(kT*tanh_f0part)) * F0/F0;
  
  for(Index iv=0; iv < nf; iv++)
  {
    const Numeric tanh_fpart = tanh( PLANCK_CONST * f_grid[iv] / kT );
    const Numeric fun = f_grid[iv] * tanh_fpart / denom;
    F[iv] *= fun;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      dF[iq][df_range][iv] *= fun;
      if(derivatives_data(iq) == JacPropMatType::Temperature)
      {
        dF[iq][df_range][iv] += (-PLANCK_CONST*(denom - F0/tanh_f0part - 
        f_grid[iv]*tanh_fpart + f_grid[iv]/tanh_fpart)/(kT*T)) * F[iv];
      }
      else if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][df_range][iv] += fac_df0 * F[iv];
      }
      else if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      {
        dF[iq][df_range][iv] += (1.0/f_grid[iv] -PLANCK_CONST*tanh_fpart/kT + PLANCK_CONST/(kT*tanh_fpart)) * F[iv];
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
 * \param quantum_identity ID of the absorption line
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_VVW_scaling(ComplexVectorView F,
                                      ArrayOfComplexVector& dF,
                                      ConstVectorView f_grid,
                                      const Numeric& F0,
                                      const PropmatPartialsData& derivatives_data,
                                      const QuantumIdentifier& quantum_identity,
                                      const Range& df_range)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // denominator is constant for the loop
  const Numeric invF0 = 1.0 / F0;
  const Numeric invF02 = invF0 * invF0;
  
  for(Index iv = 0; iv < nf; iv++)
  {
    // Set the factor
    const Numeric fac = f_grid[iv] * invF02;
    
    // Set the line shape
    F[iv] *= fac;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      // The factor is applied to all partial derivatives
      dF[iq][df_range][iv] *= fac;
      
      // These partial derivatives are special
      if(derivatives_data(iq) == JacPropMatType::LineCenter)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][df_range][iv] -= 2.0 * invF0 * F[iv] ;
      }
      else if(derivatives_data.IsFrequencyParameter(derivatives_data(iq)))
      {
        dF[iq][df_range][iv] += 2.0 / f_grid[iv] * F[iv];
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
 * \param quantum_identity ID of the absorption line
 * \param dQT_dT Temperature derivative of QT
 * \param dK1_dT Temperature derivative of K1
 * \param dK2_dT Temperature derivative of K2
 * \param dK2_dF0 Central frequency derivative of K2
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linestrength_scaling(ComplexVectorView F,
                                               ArrayOfComplexVector& dF,
                                               const Numeric& S0,
                                               const Numeric& isotopic_ratio,
                                               const Numeric& QT,
                                               const Numeric& QT0,
                                               const Numeric& K1,
                                               const Numeric& K2,
                                               const PropmatPartialsData& derivatives_data,
                                               const QuantumIdentifier& quantum_identity,
                                               const Numeric& dQT_dT,
                                               const Numeric& dK1_dT,
                                               const Numeric& dK2_dT,
                                               const Numeric& dK2_dF0,
                                               const Range& df_range)
{
  const Index nf = F.nelem();
  const Index nppd = derivatives_data.nelem();
  
  const Numeric invQT = 1.0/QT;
  const Numeric QT_ratio = QT0 * invQT;
  
  const Numeric dS_dS0 = isotopic_ratio * QT_ratio * K1 * K2;
  const Numeric S = S0 * dS_dS0;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JacPropMatType::Temperature)
    {
      const Numeric dS_dT = S0 * isotopic_ratio * QT_ratio  * (K1 * dK2_dT + dK1_dT * K2) - S * invQT * dQT_dT;
      
      dF[iq][df_range] *= S;
      
      for(Index iv = 0; iv < nf; iv++)
        dF[iq][df_range][iv] += F[iv] * dS_dT;
    }
    else if(derivatives_data(iq) == JacPropMatType::LineStrength)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq][df_range] = F;
        dF[iq][df_range] *= dS_dS0;
      }
    }
    else if(derivatives_data(iq) == JacPropMatType::LineCenter)
    {
      
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        const Numeric dS_dF0 = (S0 * isotopic_ratio * QT_ratio * K1 * dK2_dF0);
        
        dF[iq][df_range] *= S;
        
        for(Index iv = 0; iv < nf; iv++)
          dF[iq][df_range][iv] += F[iv] * dS_dF0;
      }
    }
    else
    {
      dF[iq][df_range] *= S;
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
 * \param quantum_identity ID of the absorption line
 * \param dS_LM_dT Temperature derivative of S_LM
 * 
 */
void Linefunctions::apply_linestrength_from_full_linemixing(ComplexVectorView F,
                                                            ArrayOfComplexVector& dF,
                                                            const Numeric& F0,
                                                            const Numeric& T,
                                                            const Complex& S_LM,
                                                            const Numeric& isotopic_ratio,
                                                            const PropmatPartialsData& derivatives_data,
                                                            const QuantumIdentifier& quantum_identity,
                                                            const Complex& dS_LM_dT)
{
  const Index nppd = derivatives_data.nelem();
  
  const Numeric invT = 1.0 / T, 
  F0_invT = F0 * invT,
  exp_factor = exp(C1 * F0_invT), 
  f0_factor = F0 * (1.0 - exp_factor); 
  
  const Complex S = S_LM * f0_factor * isotopic_ratio;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JacPropMatType::Temperature)
    {
      Eigen::VectorXcd eig_dF = MapToEigen(dF[iq]);
      
      eig_dF *= S_LM;
      eig_dF += MapToEigen(F) * (dS_LM_dT * f0_factor + 
      S_LM * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio;
    }
    else if(derivatives_data(iq) == JacPropMatType::LineStrength)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else if(derivatives_data(iq) == JacPropMatType::LineCenter)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else
    {
      dF[iq] *= S;
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
 * \param quantum_identity ID of the absorption line
 * \param drho_dT Temperature derivative of rho
 * 
 */
void Linefunctions::apply_dipole(ComplexVectorView F,
                                 ArrayOfComplexVector& dF,
                                 const Numeric& F0,
                                 const Numeric& T,
                                 const Numeric& d0,
                                 const Numeric& rho,
                                 const Numeric& isotopic_ratio,
                                 const PropmatPartialsData& derivatives_data,
                                 const QuantumIdentifier& quantum_identity,
                                 const Numeric& drho_dT)
{
  // Output is d0^2 * rho * F * isotopic_ratio * F0 * (1-e^(hF0/kT))
  
  const Index nppd = derivatives_data.nelem();
  
  const Numeric S = d0 * d0 * rho * isotopic_ratio, 
  invT = 1.0 / T, 
  F0_invT = F0 * invT,
  exp_factor = exp(C1 * F0_invT), 
  f0_factor = F0 * (1.0 - exp_factor);
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JacPropMatType::Temperature)
    {
      ComplexMatrixViewMap eig_dF = MapToEigen(dF[iq]);
      
      eig_dF *= S * f0_factor;
      eig_dF += MapToEigen(F) * (d0 * d0 * (drho_dT * f0_factor + 
      rho * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio);
    }
    else if(derivatives_data(iq) == JacPropMatType::LineCenter)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else if(derivatives_data(iq) == JacPropMatType::LineStrength)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else
    {
      dF[iq] *= S * f0_factor;
    }
  }
  
  F *= S * f0_factor;
}


/*! Applies the line-by-line pressure broadening jacobian for the matching lines
 * 
 * \retval dF Lineshape derivative
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param quantum_identity ID of the absorption line
 * \param dgamma Derivatives in order as they appear that are related to pressure broadening coefficients
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_pressurebroadening_jacobian_scaling(ArrayOfComplexVector& dF,
                                                              const PropmatPartialsData& derivatives_data,
                                                              const QuantumIdentifier& quantum_identity,
                                                              const ComplexVector& dgamma,
                                                              const Range& df_range)
{
  const Index nppd = derivatives_data.nelem(), ng = dgamma.nelem();
  
  Index ipd = 0;
  
  // Length of dgamma must be the same as total number of instances of pressure broadening jacobians
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(ipd == ng) break;
    
    if(derivatives_data.IsPressureBroadeningParameter(derivatives_data(iq)))
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq][df_range] *= dgamma[ipd];
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
 * \param quantum_identity ID of the absorption line
 * \param dlm Derivatives in order as they appear that are related to line mixing
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linemixing_jacobian_scaling(ArrayOfComplexVector& dF,
                                                      const PropmatPartialsData& derivatives_data,
                                                      const QuantumIdentifier& quantum_identity,
                                                      const ComplexVector& dlm,
                                                      const Range& df_range)
{
  const Index nppd = derivatives_data.nelem(), ng = dlm.nelem();
  
  Index ipd = 0;
  
  // Length of dlm must be the same as total number of instances of line mixing jacobians
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(ipd == ng) break;
    
    if(derivatives_data.IsLineMixingParameter(derivatives_data(iq)))
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq][df_range] *= dlm[ipd];
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
Numeric Linefunctions::DopplerConstant(const Numeric T, const Numeric mass)
{
  return doppler_const * sqrt(T / mass);
}


/*! Returns the temperature derivative of the frequency-independent part of the Doppler broadening
 * 
 * \param T Atmospheric temperature at level
 * \param mass Mass of molecule under consideration
 * 
 */
Numeric Linefunctions::dDopplerConstant_dT(const Numeric T, const Numeric mass)
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
                                                                   ArrayOfComplexVector& dF, 
                                                                   ComplexVectorView N, 
                                                                   ArrayOfComplexVector& dN, 
                                                                   const Numeric& K3, 
                                                                   const Numeric& K4,
                                                                   const PropmatPartialsData& derivatives_data,
                                                                   const QuantumIdentifier& quantum_identity,
                                                                   const Numeric& dK3_dT, 
                                                                   const Numeric& dK4_dT,
                                                                   const Numeric& dK3_dF0, 
                                                                   const Numeric& dK3_dTl, 
                                                                   const Numeric& dK3_dTu, 
                                                                   const Numeric& dK4_dTu,
                                                                   const Range& df_range)
{
  const Index nppd = derivatives_data.nelem(), nf = F.nelem();
  
  const Numeric scaled_ratio = K4/K3 - 1.0;
  
  // Set the non-lte source factors
  N = F;
  N *= scaled_ratio;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    dN[iq][df_range] = dF[iq][df_range];
    dN[iq][df_range] *= scaled_ratio;
    dF[iq][df_range] *= K3;
    
    if(derivatives_data(iq) == JacPropMatType::Temperature)
    {
      const Numeric dscaled_ratio_dT = (dK4_dT - dK3_dT / K3) / K3;
      
      for(Index iv = 0; iv < nf; iv++)
      {
        dF[iq][df_range][iv] += F[iv] * dK3_dT;
        dN[iq][df_range][iv] += F[iv] * dscaled_ratio_dT;
      }
    }
    else if(derivatives_data(iq) == JacPropMatType::LineCenter)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        const Numeric dscaled_ratio_dF0 = - dK3_dF0 / K3 / K3;
        
        for(Index iv = 0; iv < nf; iv++)
        {
          dF[iq][df_range][iv] += F[iv] * dK3_dF0;
          dN[iq][df_range][iv] += F[iv] * dscaled_ratio_dF0;
        }
      }
    }
    else if(derivatives_data(iq) == JacPropMatType::VibrationalTemperature)
    {
      if(derivatives_data.jac(iq).QuantumIdentity().Species() not_eq quantum_identity.Species() or
        derivatives_data.jac(iq).QuantumIdentity().Isotopologue() not_eq quantum_identity.Isotopologue())
        continue;  // Wrong species or wrong isotopologue
      
      if(quantum_identity.QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX] >
        derivatives_data.jac(iq).QuantumIdentity().QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL] or
        derivatives_data.jac(iq).QuantumIdentity().Type() == QuantumIdentifier::ALL)
      {
        const Numeric dscaled_ratio_dTl = - dK3_dTl / K3 / K3;
        
        for(Index iv = 0; iv < nf; iv++)
        {
          dF[iq][df_range][iv] += F[iv] * dK3_dTl;
          dN[iq][df_range][iv] += F[iv] * dscaled_ratio_dTl;
        }
      }
      
      if(quantum_identity.QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX] >
        derivatives_data.jac(iq).QuantumIdentity().QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL] or
        derivatives_data.jac(iq).QuantumIdentity().Type() == QuantumIdentifier::ALL)
      {
        const Numeric dscaled_ratio_dTu = (dK4_dTu - dK3_dTu / K3) / K3;
        
        for(Index iv = 0; iv < nf; iv++)
        {
          dF[iq][df_range][iv] += F[iv] * dK3_dTu;
          dN[iq][df_range][iv] += F[iv] * dscaled_ratio_dTu;
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
void Linefunctions::set_cross_section_for_single_line(ComplexVectorView F,
                                                      ArrayOfComplexVector& dF,
                                                      ComplexVectorView N, 
                                                      ArrayOfComplexVector& dN,
                                                      Range& this_f_range,
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
                                                      const bool cutoff_call,
                                                      const Index binary_speedup)
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
  const bool need_cutoff = cutoff_call ? false : find_cutoff_ranges(this_f_range, f_grid, line.F(), cutoff);
  const bool need_binary_speedup = binary_speedup and line.SpeedUpCoeff() > 0.0;
  if(need_cutoff and need_binary_speedup)
    throw std::runtime_error("Cannot have booth cutoff and binary speedup at same time...");
                           
  // Leave this function if there is nothing to compute
  if(this_f_range.get_extent() == 0)
    return;
  
  // Extract the quantum identify of the line to be used in-case there are derivatives
  const QuantumIdentifier QI = line.QuantumIdentity();
  
  /* Pressure broadening terms
   * These are set by the line catalog.  There are no defaults.
   */
  Numeric G0, G2, e, L0, L2, FVC;
  line.PressureBroadening().GetPressureBroadeningParams(
    G0, G2, e, L0, L2, FVC,
    line.Ti0()/temperature, pressure, partial_pressure, 
    this_species_location_in_tags, water_index_location_in_tags,
    broad_spec_locations, volume_mixing_ratio_of_all_species,
    verbosity);
  
  // Line mixing terms
  Numeric Y=0, G=0, DV=0;
  line.LineMixing().GetLineMixingParams(Y, G, DV, temperature, pressure, 
                                        pressure_limit_for_linemixing);
  
  // Partial derivatives for temperature
  Numeric dG0_dT, dL0_dT, dG2_dT=0.0, dL2_dT=0.0, de_dT=0.0, dFVC_dT=0.0, 
          dY_dT=0, dG_dT=0, dDV_dT=0, dK1_dT, dK2_dT;
  if(derivatives_data.do_temperature())
  {
    // Pressure broadening partial derivatives
    line.PressureBroadening().GetPressureBroadeningParams_dT(
      dG0_dT, dG2_dT, de_dT, dL0_dT, dL2_dT, dFVC_dT,
      temperature, line.Ti0(), pressure, partial_pressure,
      this_species_location_in_tags, water_index_location_in_tags,
      broad_spec_locations, volume_mixing_ratio_of_all_species,
      verbosity);
    
    // Line mixing partial derivatives
    line.LineMixing().GetLineMixingParams_dT(
      dY_dT, dG_dT, dDV_dT, temperature, 
      derivatives_data.Temperature_Perturbation(),
      pressure, pressure_limit_for_linemixing);
  }
  
  /* Partial derivatives due to pressure
   * The vector below will be rescaled by the set internal derivatives function such that
   * the order of their occurrences are the same as in the partial derivative output
   */
  ComplexVector pressure_derivatives;
  if(derivatives_data.get_first_pressure_term() > -1)
    line.PressureBroadening().SetInternalDerivatives(pressure_derivatives, derivatives_data, QI, 
                                                     line.Ti0()/temperature, pressure, partial_pressure, 
                                                     this_species_location_in_tags, water_index_location_in_tags,
                                                     volume_mixing_ratio_of_all_species, verbosity);
    
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
  
  // Speedup variables
  Index speedup_index=0, computational_points=0, lower_boundary, upper_boundary;
  Range speedup_range(joker);
  
  if(need_binary_speedup)
    find_boundary_of_binary_range(upper_boundary, lower_boundary, f_grid, line.SpeedUpCoeff(),
                                  G0, line.F(), doppler_constant, binary_speedup);
    
  do
  { 
    // Setup the speedup range if necessary
    if(need_binary_speedup)
    {
      // binary_range sets the grids to [1, n-1], [(n-1)/2], [(n-1)/4, 3*(n-1)/4], [n/8, 3*n/8, 5*n/8, 7*n/8], ...
      speedup_range = binary_range(f_grid, speedup_index, false);
      
      // Distance from the line center where the range is computed
      const Numeric d = speedup_distance_binary_range(binary_speedup-speedup_index, line.SpeedUpCoeff(), G0, line.F(), doppler_constant);
      
      // Only select a part of the range
      if(speedup_index)
        speedup_range = speedup_range(speedup_binary_range(f_grid[speedup_range], line.F() + d, line.F() - d));
      
      /* 
       * If the first range to be computed after the boundary stretches far outside
       * the range of the speedup-region, the interpolation to an upper level will
       * contain zeroes.  To combat this, the first range to be computed is computed
       * fully, so that the interpolation to a denser range is made easy.
       * 
       * These cases are either the first time a speedup-range has an extent or if the 
       * middle point has beed computed
       */
      if((computational_points == 2 and speedup_range.get_extent()) or computational_points == 3)
        speedup_range = Range(speedup_range.get_start(), 2*speedup_range.get_extent(), speedup_range.get_stride()/2);
      
      // Add up the number of points being computed
      computational_points += speedup_range.get_extent();
    }
    
    if(speedup_range.get_extent())
    {
      /*! Set the line shape normalized to unity integration
      * The user can set this by LSM LST followed by an index that 
      * is interpreted internally as a line shape.
      * The main point is not that the user should use such functions 
      * but that support functions can set the catalog, and that once
      * stored the catalog will use that line shape.  If no line shape 
      * tag is given, the line shape will be set by the type of pressure
      * broadening data that has been provided.
      */
      switch(line.GetLineShapeType())
      {
        // Use data as provided by the pressure broadening scheme
        case LineShapeType::ByPressureBroadeningData:
          switch(line.PressureBroadening().Type())
          {
            // Use data as per speed dependent air
            case PressureBroadeningData::PB_SD_AIR_VOLUME:
            case PressureBroadeningData::PB_PURELY_FOR_TESTING:
              lst = LineShapeType::HTP;
              set_htp(F[this_f_range(speedup_range)], dF, 
                      f_grid[this_f_range(speedup_range)], line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                      line.F(), doppler_constant, 
                      G0, L0, G2, L2, e, FVC,
                      derivatives_data, QI,
                      ddoppler_constant_dT, 
                      dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT,
                      this_f_range(speedup_range));
              break;
            // Use for data that requires air and water Voigt broadening
            case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
            // Use for data that requires planetary Voigt broadening
            case PressureBroadeningData::PB_PLANETARY_BROADENING:
              // Use for data that requires air Voigt broadening
            case PressureBroadeningData::PB_AIR_BROADENING:
              // Above should be all methods of pressure broadening requiring Voigt in ARTS by default
              lst = LineShapeType::Voigt;
              set_faddeeva_algorithm916(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], 
                                        line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                                        line.F(), doppler_constant, 
                                        G0, L0, DV, derivatives_data, QI,
                                        ddoppler_constant_dT, dG0_dT, dL0_dT, dDV_dT,
                                        this_f_range(speedup_range));
              break;
            default:
              throw std::runtime_error("Developer has messed up and needs to add the key to the code above this error");
          }
          break;
        // This line only needs the Doppler effect
        case LineShapeType::Doppler:
          lst = LineShapeType::Doppler;
          set_doppler(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                      line.F(), doppler_constant, derivatives_data, QI, ddoppler_constant_dT, this_f_range(speedup_range));
          break;
        // This line only needs Hartmann-Tran
        case LineShapeType::HTP:
          set_htp(F[this_f_range(speedup_range)], dF, 
                  f_grid[this_f_range(speedup_range)], line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                  line.F(), doppler_constant, 
                  G0, L0, G2, L2, e, FVC,
                  derivatives_data, QI,
                  ddoppler_constant_dT, 
                  dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT,
                  this_f_range(speedup_range));
          lst = LineShapeType::HTP;
          break;
        // This line only needs Lorentz
        case LineShapeType::Lorentz:
          lst = LineShapeType::Lorentz;
          set_lorentz(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                      line.F(), G0, L0, DV, derivatives_data, QI, dG0_dT, dL0_dT, dDV_dT, this_f_range(speedup_range));
          break;
        // This line only needs Voigt
        case LineShapeType::Voigt:
          lst = LineShapeType::Voigt;
          set_faddeeva_algorithm916(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], 
                                    line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                                    line.F(), doppler_constant, 
                                    G0, L0, DV, derivatives_data, QI,
                                    ddoppler_constant_dT, dG0_dT, dL0_dT, dDV_dT,
                                    this_f_range(speedup_range));
          break;
        case LineShapeType::End:
        default:
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
            ComplexVector Fm(F[this_f_range(speedup_range)].nelem());
            ArrayOfComplexVector dFm(dF.nelem());
            for(auto& aocv : dFm) aocv.resize(F[this_f_range(speedup_range)].nelem());
            
            set_lorentz(Fm, dFm, f_grid[this_f_range(speedup_range)], -line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                        -line.F(), G0, -L0, -DV, derivatives_data, QI, dG0_dT, -dL0_dT, -dDV_dT, this_f_range(speedup_range));
            
            // Apply mirroring
            F[this_f_range(speedup_range)] -= Fm;
            for(Index i = 0; i < dF.nelem(); i++) dF[i][this_f_range(speedup_range)] -= dFm[i];
          }
          break;
        // Same type of mirroring as before
        case MirroringType::SameAsLineShape:
        {
          // Set the mirroring computational vectors and size them as needed
          ComplexVector Fm(F[this_f_range(speedup_range)].nelem());
          ArrayOfComplexVector dFm(dF.nelem());
          for(auto& aocv : dFm) aocv.resize(F[this_f_range(speedup_range)].nelem());
          
          switch(lst)
          {
            case LineShapeType::Doppler:
              set_doppler(Fm, dFm, f_grid[this_f_range(speedup_range)], -line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                          -line.F(), -doppler_constant, derivatives_data, QI, -ddoppler_constant_dT, this_f_range(speedup_range));
              break;
            case LineShapeType::Lorentz:
              set_lorentz(Fm, dFm, f_grid[this_f_range(speedup_range)], -line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                          -line.F(), G0, -L0, -DV, derivatives_data, QI, dG0_dT, -dL0_dT, -dDV_dT, this_f_range(speedup_range));
              break;
            case LineShapeType::Voigt:
              set_faddeeva_algorithm916(Fm, dFm, f_grid[this_f_range(speedup_range)], 
                                        -line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                                        -line.F(), -doppler_constant, 
                                        G0, -L0, -DV, derivatives_data, QI,
                                        -ddoppler_constant_dT, dG0_dT, -dL0_dT, -dDV_dT, this_f_range(speedup_range));
              break;
            case LineShapeType::HTP:
              // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
              set_htp(Fm, dFm, f_grid[this_f_range(speedup_range)], 
                      -line.ZeemanEffect().frequency_shift_per_tesla(line.QuantumNumbers(), line.Species()), magnetic_magnitude, 
                      -line.F(), -doppler_constant, 
                      G0, -L0, G2, -L2, e, FVC,
                      derivatives_data, QI,
                      -ddoppler_constant_dT, 
                      dG0_dT, -dL0_dT, dG2_dT, -dL2_dT, de_dT, dFVC_dT, this_f_range(speedup_range));
              break;
            case LineShapeType::ByPressureBroadeningData:
            case LineShapeType::End:
            default:
              throw std::runtime_error("Cannot understand the requested line shape type for mirroring.");
          }
          F[this_f_range(speedup_range)] -= Fm;
          for(Index i = 0; i < dF.nelem(); i++) dF[i][this_f_range(speedup_range)] -= dFm[i];
          break;
        }
        case MirroringType::End:
        default:
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
          apply_VVH_scaling(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], line.F(), temperature, derivatives_data, QI, this_f_range(speedup_range));
          break;
          // van Vleck and Weiskopf normalization
        case LineNormalizationType::VVW:
          apply_VVW_scaling(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], line.F(), derivatives_data, QI, this_f_range(speedup_range));
          break;
          // Rosenkranz's Quadratic normalization
        case LineNormalizationType::RosenkranzQuadratic:
          apply_rosenkranz_quadratic_scaling(F[this_f_range(speedup_range)], dF, f_grid[this_f_range(speedup_range)], line.F(), temperature, derivatives_data, QI, this_f_range(speedup_range));
          break;
        case LineNormalizationType::End:
        default:
          throw std::runtime_error("Cannot understand the requested line normalization type.");
      }
    }
    
    // Interpolate to a higher level if using speedup
    if(need_binary_speedup)
    {
      // Interpolation up needs to happen if we are at a stage of having more than 2 points
      if((speedup_index != binary_speedup and computational_points > 2) or computational_points == 3)
      {
        interp_up_inside_binary_range(F, binary_range(f_grid, speedup_index, true));
        for(auto& cv : dF)
          interp_up_inside_binary_range(cv, binary_range(f_grid, speedup_index, true));
      }
      else if(speedup_index == binary_speedup)
      {
        if(upper_boundary > lower_boundary) 
        {
          interp_to_boundary_of_binary_range(F, upper_boundary, lower_boundary); 
          for(auto& cv : dF)
            interp_to_boundary_of_binary_range(cv, upper_boundary, lower_boundary); 
        } 
        else 
        {
          interp_to_boundary_of_binary_range(F, 0, 0); 
          for(auto& cv : dF)
            interp_to_boundary_of_binary_range(cv, 0, 0); 
        }
      }
      speedup_index++;
    }
  } while(speedup_index <= binary_speedup and need_binary_speedup); // End of frequency loop
  
  // Apply line mixing if relevant
  if(Y not_eq 0 or G not_eq 0)
    apply_linemixing_scaling(F[this_f_range], dF, Y, G, derivatives_data, QI, dY_dT, dG_dT, this_f_range);
  
  // Apply pressure broadening partial derivative vector if necessary
  if(pressure_derivatives.nelem()>0)
    apply_pressurebroadening_jacobian_scaling(dF, derivatives_data, QI, pressure_derivatives, this_f_range);
  
  // Apply line mixing partial derivative vector if necessary
  if(linemixing_derivatives.nelem()>0)
    apply_linemixing_jacobian_scaling(dF, derivatives_data, QI, linemixing_derivatives, this_f_range);
  
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
        
        if(derivatives_data.do_temperature())
        {
          dK1_dT = dboltzman_ratio_dT(K1, temperature, line.Elow());
          dK2_dT = dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), temperature);
        }
        
        // Partial derivatives due to central frequency of the stimulated emission
        Numeric dK2_dF0;
        if(derivatives_data.do_line_center())
          dK2_dF0 = dstimulated_relative_emission_dF0(gamma, gamma_ref, temperature);
        
        // Multiply the line strength by the line shape
        apply_linestrength_scaling(F[this_f_range], dF,  line.I0(), isotopologue_ratio,
                                   partition_function_at_temperature, partition_function_at_line_temperature, K1, K2,
                                   derivatives_data, QI, dpartition_function_at_temperature_dT, dK1_dT, dK2_dT, dK2_dF0, this_f_range);
        
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
          if(derivatives_data.do_line_center())
            dK3_dF0 = dabsorption_nlte_rate_dF0(gamma, temperature, K4, r_low);
          
          // Are we computing the temperature derivatives?
          // NOTE:  Having NLTE active AT ALL will change the jacobian because of this part of the code,
          // though this requires setting El and Eu for all lines, though this is not yet default...
          // So if you see this part of the code after having a runtime_error, 
          // you will need to write those functions yourself...
          if(derivatives_data.do_temperature())
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
          set_nonlte_source_and_apply_absorption_scaling(F[this_f_range], dF, N[this_f_range], dN, K3, K4, 
                                                         derivatives_data, QI,  dK3_dT, dK4_dT, dK3_dF0, dK3_dTl, 
                                                         dK3_dTu, dK4_dTu, this_f_range);
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
        
        apply_linestrength_from_nlte_level_distributions(F[this_f_range], 
                                                         dF,
                                                         N[this_f_range], 
                                                         dN,
                                                         nlte_distribution[nlte_low_index],
                                                         nlte_distribution[nlte_upp_index],
                                                         line.G_lower(),
                                                         line.G_upper(),
                                                         line.A(),
                                                         line.F(),
                                                         temperature, 
                                                         derivatives_data, 
                                                         QI,
                                                         this_f_range);
      }
      break;
    default:
      throw std::runtime_error("Unknown population distribution type for this line");
  }
  
  // Cutoff frequency is applied at the end because 
  // the entire process above is complicated and applying
  // cutoff last means that the code is kept cleaner
  if(need_cutoff)
  {
    apply_cutoff(F[this_f_range], dF,
                 nlte_distribution.nelem()?N[this_f_range]:N, dN,
                 derivatives_data, line,
                 volume_mixing_ratio_of_all_species,
                 nlte_distribution, pressure, temperature,
                 doppler_constant, partial_pressure,
                 isotopologue_ratio, magnetic_magnitude,
                 ddoppler_constant_dT, pressure_limit_for_linemixing,
                 partition_function_at_temperature,
                 dpartition_function_at_temperature_dT,
                 partition_function_at_line_temperature,
                 broad_spec_locations, this_species_location_in_tags,
                 water_index_location_in_tags, this_f_range, verbosity);
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
                                 const Verbosity& verbosity)
{ 
  // Size of derivatives
  const Index nj = dF.nelem(); 
  const Index nn = dN.nelem(); 
  
  // Setup compute variables
  Vector f_grid_cutoff(1, line.F() + line.CutOff());
  ComplexVector Fc(1), Nc(1);
  ArrayOfComplexVector dFc(nj), dNc(nn);
  for(auto& aovc : dFc) aovc.resize(1);
  for(auto& aovc : dNc) aovc.resize(1);
  Range tmp(joker);
  
  // Recompute the line for a single frequency
  set_cross_section_for_single_line(Fc, dFc, Nc, dNc, tmp,
                                    derivatives_data, line, f_grid_cutoff,
                                    volume_mixing_ratio_of_all_species,
                                    nlte_distribution, pressure, temperature,
                                    doppler_constant, partial_pressure, isotopologue_ratio,
                                    magnetic_magnitude, ddoppler_constant_dT,
                                    pressure_limit_for_linemixing,
                                    partition_function_at_temperature,
                                    dpartition_function_at_temperature_dT,
                                    partition_function_at_line_temperature,
                                    broad_spec_locations, this_species_location_in_tags,
                                    water_index_location_in_tags, verbosity, true, false);
  
  // Apply cutoff values
  F -= Fc[0];
  if(N.nelem())
    N -= Nc[0];
  for(Index i = 0; i < nj; i++)
    dF[i][df_range] -= dFc[i][0];
  for(Index i = 0; i < nn; i++)
    dN[i][df_range] -= dNc[i][0];
}


bool Linefunctions::find_cutoff_ranges(Range& range,
                                       ConstVectorView f_grid,
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
 * 
 */
void Linefunctions::apply_linestrength_from_nlte_level_distributions(ComplexVectorView F, 
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
                                                                     const PropmatPartialsData& derivatives_data,
                                                                     const QuantumIdentifier& quantum_identity,
                                                                     const Range& df_range)
{
  // Size of the problem
  const Index nf = F.nelem();
  const Index nppd = derivatives_data.nelem();
  
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
    if(derivatives_data(iq) == JacPropMatType::Temperature)  // nb.  Only here tou counter-act the latter on in 
    { 
      Numeric done_over_b_dT = PLANCK_CONST*F0*exp_T/(c2*BOLTZMAN_CONST*T*T);
      
      // dN is unset so far.  It should return to just be lineshape later...
      for(Index iv=0; iv<nf; iv++)
        dN[iq][df_range][iv] = F[iv]*e*done_over_b_dT + dF[iq][df_range][iv]*ratio;
      
      // dk_dT = 0...
      dF[iq][df_range] *= k;
    }
    else if(derivatives_data(iq) == JacPropMatType::LineCenter)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        Numeric done_over_b_df0 = PLANCK_CONST*exp_T/(c2*BOLTZMAN_CONST*T) - 3.0*b/F0;
        Numeric de_df0 = c1 * r2 * A21;
        Numeric dk_df0 = c1 * (r1*x - r2) * (A21 / c2) - 3.0*k/F0;
        
        for(Index iv=0; iv<nf; iv++)
        {
          dN[iq][df_range][iv] = F[iv]*(e*done_over_b_df0 + de_df0/b - dk_df0) + dF[iq][df_range][iv]*ratio;
          dF[iq][df_range][iv] = dF[iq][df_range][iv]*k + F[iv]*dk_df0;
        }
      }
    }
    else if(derivatives_data(iq) == JacPropMatType::PopulationRatio)
    {
      if(quantum_identity.Lower() > derivatives_data.jac(iq).QuantumIdentity())
      {
        const Numeric dk_dr2 = - c3 * A21 / c2, de_dr2 = c3 * A21, dratio_dr2 = de_dr2/b - dk_dr2;
        
        for(Index iv=0; iv<nf; iv++)
        {
          dN[iq][df_range][iv] = F[iv]*dratio_dr2;
          dF[iq][df_range][iv] = F[iv]*dk_dr2 + dF[iq][df_range][iv]*k;
        } 
      }
      else if(quantum_identity.Upper() > derivatives_data.jac(iq).QuantumIdentity())
      {
        const Numeric dk_dr1 = c3 * x * A21 / c2;
        
        for(Index iv=0; iv<nf; iv++)
          dF[iq][df_range][iv] = F[iv]*dk_dr1 + dF[iq][df_range][iv]*k;
      }
    }
  }
  
  // Set source function to be the relative amount emitted by the line divided by the Planck function
  N = F;
  N *= ratio;
  
  // Set absorption
  F *= k;
}


/*!
 * Returns a range on a binary scale
 * 
 * \param f A vector that has the length 2^N + 1
 * \param i An index less or equal to N
 * \param full_range If false, output range does not contain values of lower i
 * 
 */
Range Linefunctions::binary_range(const ConstVectorView f, const Index& i, const bool full_range)
{
  // Size of problem
  const Index n = f.nelem();
  
  // The first range is always complete and contains end-points
  if(not i)
    return Range(0, 2, n-1);
  
  // The full range is 2^i + 1, the partial range is half minus 1 long
  const Index d = 1 << i;
  
  // The full range stepsize is thus simply
  const Index s = (n - 1) / d;
  
  // And we now have the range, either the full or the partial range...
  if(full_range)
    return Range(0, d + 1, s);
  else
    return Range(s, d >> 1, s << 1);
}


/*!
 * Returns the range of f where l < f[i] < = u
 * 
 * Assumes f is evenly spaced
 * 
 * \param f A vector that has the length 2^N + 1
 * \param u An upper boundary
 * \param l A lower boundary
 * 
 */
Range Linefunctions::speedup_binary_range(const ConstVectorView f, 
                                          const Numeric& u, 
                                          const Numeric& l)
{
  // Size of problem
  const Index n = f.nelem();
  
  // Treat the special cases first
  if(not n)
    return Range(0, 0);
  else if(n == 1)
  {
    if(f[0] < u and f[0] > l)
      return Range(0, 1);
    else
      return Range(0, 0);
  }
  else if(n == 2)
  {
    if((f[0] > u) or (f[1] < l) or (f[0] < l and f[1] > u))
      return Range(0, 0);
    else if(f[1] < u and f[0] > l)
      return Range(0, 2);
    else if(f[1] < u)
      return Range(1, 1);
    else
      return Range(0, 1);
  }
  
  // Equidistance in f required
  const Numeric d = f[1] - f[0];
  
  Index nl = (Index)((l-f[0])/d);
  Index nu = (Index)((u-f[0])/d)+2;
  
  // Lower start-position adjustments
  if(nl < 0)
    nl = 0;
  else if(nl > n)
    nl = n;
  
  // Upper end-position adjustment
  if(nu < 0)
    nu = 0;
  else if(nu > n)
    nu = n;
  
  return Range(nl, nu - nl);
}


/*!
 * Interpolates to a denser grid 
 * 
 * \param f A vector that is to become denser by linear interpolation
 * \param l A range of the vector representing a lower density grid by a factor 2
 * 
 */
void Linefunctions::interp_up_inside_binary_range(ComplexVectorView f, const Range& l)
{
  Index II = l.get_start();
  const Index dII = l.get_stride();
  const Index dII2 = dII / 2;
  
  for(Index i = 0; i < (l.get_extent()-1); i++)
  {
    f[II + dII2] = 0.5 * (f[II] + f[II + dII]);
    II += dII;
  }
}


/*!
 * Interpolates to the boundary
 * 
 * \param f A vector that has values at f[0], f[l], f[u], and f[n], 
 * where 0 < l < u < n, and n is the length of f minus one
 * \param u Upper boundary index
 * \param l Lower boundary index
 * 
 */
void Linefunctions::interp_to_boundary_of_binary_range(ComplexVectorView f, const Index u, const Index l)
{
  const Index n = f.nelem();
  Index i;
  Complex cl, cu, df;
  
  if(l > 1)
  {
    cl = f[0];
    cu = f[l];
    df = (cu - cl) / Numeric(l);
    i = 1;
    while(i < l)
    {
      f[i] = cl + df * Numeric(i);
      i++;
    }
  }
  
  if(u < n-2)
  {
    cl = f[u];
    cu = f[n-1];
    const Index N = n - u - 1;
    df = (cu - cl) / Numeric(N);
    i = 1;
    while(i < N)
    {
      f[u + i] = cl + df * Numeric(i);
      i++;
    }
  }
}


/*!
 * The distance of the speedup
 * 
 * \param s An index of the speedup so that s>0 and s less than the speedup order
 * \param C A constant describing a scaling factor for the distance and the line
 * \param G0 Speed-independent pressure broadening term
 * \param F0 Central frequency
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * 
 */
Numeric Linefunctions::speedup_distance_binary_range(const Index s, 
                                                     const Numeric C, 
                                                     const Numeric G0, 
                                                     const Numeric F0,
                                                     const Numeric GD_div_F0)
{
  const Numeric GD = GD_div_F0*F0;
  
  // 2^(s+1) * C * L
  return Numeric(2 << s) * C * sqrt(G0*G0 + GD*GD);
}


/*!
 * Find out what the lower and upper boundary of the binary range will be...
 * 
 * \param l An index indicating the lower boundary
 * \param u An index indicating the upper boundary
 * \param f A vector that has the length 2^N + 1
 * \param C A constant describing a scaling factor for the distance and the line
 * \param G0 Speed-independent pressure broadening term
 * \param F0 Central frequency
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param binary_speedup N in the description of f...
 */
void Linefunctions::find_boundary_of_binary_range(Index& u, 
                                                  Index& l,
                                                  ConstVectorView f, 
                                                  const Numeric C, 
                                                  const Numeric G0,
                                                  const Numeric F0,
                                                  const Numeric GD_div_F0,
                                                  const Index binary_speedup)
{
  const Index nf = f.nelem();
  Index t;
  u = 0;
  l = nf - 1;
  
  for(Index i = 1; i <= binary_speedup; i++)
  {
    const Numeric d = speedup_distance_binary_range(binary_speedup - i, C, G0, F0, GD_div_F0);
    
    Range r = binary_range(f, i, false);
    r = r(speedup_binary_range(f[r], F0 + d, F0 - d));
    
    if(r.get_extent())
    {
      t = r.get_start();
      if(t < l)
        l = t;
      
      t = (r.get_extent() - 1) * r.get_stride() + r.get_start();
      if(t > u)
        u = t;
    }
  }
}

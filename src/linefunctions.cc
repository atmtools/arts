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
                                  const Numeric& magnetic_magnitude,
                                  const ArrayOfArrayOfSpeciesTag& abs_species,
                                  const Index& this_species,
                                  const Index& zeeman_index)
{
  // Pressure broadening and line mixing terms
  const auto X = line.GetShapeParams(temperature, pressure, this_species, vmrs, abs_species);
  const auto& G0  = std::get<Index(LineFunctionData::TuplePos::G0 )>(X);
  const auto& D0  = std::get<Index(LineFunctionData::TuplePos::D0 )>(X);
  const auto& G2  = std::get<Index(LineFunctionData::TuplePos::G2 )>(X);
  const auto& D2  = std::get<Index(LineFunctionData::TuplePos::D2 )>(X);
  const auto& FVC = std::get<Index(LineFunctionData::TuplePos::FVC)>(X);
  const auto& ETA = std::get<Index(LineFunctionData::TuplePos::ETA)>(X);
//   const auto& Y   = std::get<Index(LineFunctionData::TuplePos::Y  )>(X);
//   const auto& G   = std::get<Index(LineFunctionData::TuplePos::G  )>(X);
  const auto& DV  = std::get<Index(LineFunctionData::TuplePos::DV )>(X);
  
  // Line shape usage remembering variable
  LineShapeType lst = LineShapeType::End;
  
  ComplexMatrix dF(0, 0);
  const Numeric doppler_constant = DopplerConstant(temperature, line.IsotopologueData().Mass());
  
  switch(line.GetExternalLineShapeType())
  {
    case LineShapeType::ByPressureBroadeningData:
      switch(line.GetInternalLineShapeType())
      {
        case LineFunctionData::LineShapeType::HTP:
        case LineFunctionData::LineShapeType::SDVP:
          lst = LineShapeType::HTP;
          set_htp(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, D0, G2, D2, ETA, FVC);
          break;
        case LineFunctionData::LineShapeType::VP:
          lst = LineShapeType::Voigt;
          set_voigt(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, D0, DV);
          break;
        case LineFunctionData::LineShapeType::DP:
          lst = LineShapeType::Doppler;
          set_doppler(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant);
          break;
        case LineFunctionData::LineShapeType::LP:
          lst = LineShapeType::Lorentz;
          set_lorentz(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), G0, D0, DV);
          break;
      }
      break;
        case LineShapeType::Doppler:
          lst = LineShapeType::Doppler;
          set_doppler(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant);
          break;
          // This line only needs Hartmann-Tran
        case LineShapeType::HTP:
          set_htp(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, D0, G2, D2, ETA, FVC);
          lst = LineShapeType::HTP;
          break;
          // This line only needs Lorentz
        case LineShapeType::Lorentz:
          lst = LineShapeType::Lorentz;
          set_lorentz(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), G0, D0, DV);
          break;
          // This line only needs Voigt
        case LineShapeType::Voigt:
          lst = LineShapeType::Voigt;
          set_voigt(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, G0, D0, DV);
          break;
        case LineShapeType::End:
          throw std::runtime_error("Cannot understand the requested line shape type.");
  }
  
  switch(line.GetMirroringType())
  {
    // No mirroring
    case MirroringType::None:
    case MirroringType::Manual:
      break;
      // Lorentz mirroring
    case MirroringType::Lorentz:
    {
      // Set the mirroring computational vectors and size them as needed
      ComplexVector Fm(F.nelem());
      ComplexMatrix dFm(0, 0);
      
      set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), G0, -D0, -DV);
      
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
          set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), G0, -D0, -DV);
          break;
        case LineShapeType::Voigt:
          set_voigt(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant, G0, -D0, -DV);
          break;
        case LineShapeType::HTP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant, G0, -D0, G2, -D2, ETA, FVC);
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
                                const Numeric& ddF0_dT,
                                const LineFunctionData::Output& dVMR)
{ 
  // Size of the problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
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
      
      if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
        dF(iq, iv) = d * Complex(dG0_dT, dL0_dT + ddF0_dT);  // Temperature derivative only depends on how pressure shift and broadening change
      else if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]]))
        dF(iq, iv) = d * Complex(0.0, -1.0);  // Frequency scale 1 to -1 linearly
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter or 
              is_line_mixing_DF_parameter(derivatives_data[derivatives_data_position[iq]])) {
        if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity))
          dF(iq, iv) = d * Complex(0.0, 1.0);  // Line center scales 1 to 1 linearly
      }
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]])) {
        if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity))
          dF(iq, iv) = d * Complex(0.0, -1.0);  // Pressure broadening will be dealt with in another function, though the partial derivative
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR) {
        if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity))
          dF(iq, iv) = d * Complex(std::get<Index(LineFunctionData::TuplePos::G0)>(dVMR), 
                                   std::get<Index(LineFunctionData::TuplePos::D0)>(dVMR) +
                                   std::get<Index(LineFunctionData::TuplePos::DV)>(dVMR));
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[iq]]))
        dF(iq, iv) = d * Complex(0.0, zeeman_df);  // Magnetic magnitude changes like line center in part
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
  
  static const Complex i(0.0, 1.0);
  
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
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]))  {
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
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR) {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
          throw std::runtime_error("No support for internal VMR partial derivatives for this line, use derivatives at a higher level");
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
                              const Numeric& dF0_dT,
                              const LineFunctionData::Output& dVMR)
{
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  // For calculations
  Complex w, z, dw_over_dz;
  
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
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) /*or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR*/) {
        if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
          dF(iq, iv) = dw_over_dz * invGD;
      }
      else if(is_magnetic_magnitude_parameter(derivatives_data[derivatives_data_position[iq]]))
        dF(iq, iv) = dw_over_dz * (- zeeman_df * invGD); //* dz;
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR) {
        if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity))
          dF(iq, iv) = dw_over_dz * Complex(- std::get<Index(LineFunctionData::TuplePos::D0)>(dVMR)
                                            - std::get<Index(LineFunctionData::TuplePos::DV)>(dVMR),
                                              std::get<Index(LineFunctionData::TuplePos::G0)>(dVMR)) * invGD;
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
  const Numeric fac = invGD * sqrtInvPI;
  
  // Computational speed-up
  const Index first_frequency = get_first_frequency_index(derivatives_data, derivatives_data_position);
  
  for(Index iv = 0; iv < nf; iv++) {
    const Numeric x = (f_grid[iv] - F0) * invGD;
    F[iv] = fac * exp(- x*x);
    
    for(Index iq = 0; iq < nppd; iq++) { 
      if(is_frequency_parameter(derivatives_data[derivatives_data_position[iq]] )) {
        if(first_frequency == iq)
          dF(iq, iv) = -2 * F[iv] * x * invGD;  // If this is the first time it is calculated this frequency bin, do the full calculation
        else // copy for repeated occurrences
          dF(iq, iv) = dF(first_frequency, iv);
      }
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature)
        dF(iq, iv) = F[iv] * (2*x*x - 1) * dGD_dT * invGD;
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter) {
        if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity)) 
          dF(iq, iv) = F[iv] * 2*x*((x-1)/F0 + invGD);
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
      else if(is_pressure_broadening_parameter(derivatives_data[derivatives_data_position[iq]]) /*or derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR*/)
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
                                             const Numeric& dG_dT,
                                             const LineFunctionData::Output& dVMR)
{
  const Index nf = F.nelem(), nppd = derivatives_data_position.nelem();
  
  const Complex LM = Complex(1.0 + G, -Y);
  const Complex dLM_dT = Complex(dG_dT, -dY_dT);
  
  for(Index iq = 0; iq < nppd; iq++) {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature) {
      dF(iq, joker) *= LM;
      for(Index iv = 0; iv < nf; iv++)
        dF(iq, iv) += F[iv] * dLM_dT;
    }
    else if(is_line_mixing_line_strength_parameter(derivatives_data[derivatives_data_position[iq]])) {
       if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity())
         dF(iq, joker) = F;
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::VMR) {
      if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity))
        for(Index iv = 0; iv < nf; iv++)
          dF(iq, iv) = dF(iq, iv) * LM +  F[iv] * Complex(  std::get<Index(LineFunctionData::TuplePos::G)>(dVMR),
                                                          - std::get<Index(LineFunctionData::TuplePos::Y)>(dVMR));
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
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter)
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
      else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter)
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
      if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter)
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
  
  for(Index iq = 0; iq < nppd; iq++) {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature) {
      const Numeric dS_dT = S * (dK2_dT / K2 + dK1_dT / K1 - invQT * dQT_dT);
      
      dF(iq, joker) *= S;
      for(Index iv = 0; iv < nf; iv++)
        dF(iq, iv) += F[iv] * dS_dT;
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineStrength) {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        dF(iq, joker) = F;
        dF(iq, joker) *= dS_dS0;
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter) {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        const Numeric dS_dF0 = S * dK2_dF0 / K2;
        
        dF(iq, joker) *= S;
        for(Index iv = 0; iv < nf; iv++)
          dF(iq, iv) += F[iv] * dS_dF0;
      }
    }
    else
      dF(iq, joker) *= S;
  }
  
  // Set lineshape at the end
  F *= S;
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
void Linefunctions::apply_linestrength_scaling_vibrational_nlte(ComplexVectorView F,
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
                                                                const ArrayOfRetrievalQuantity& derivatives_data,
                                                                const ArrayOfIndex& derivatives_data_position,
                                                                const QuantumIdentifier& quantum_identity,
                                                                const Numeric& dQT_dT,
                                                                const Numeric& dK1_dT,
                                                                const Numeric& dK2_dT,
                                                                const Numeric& dK2_dF0,
                                                                const Numeric& dK3_dT,
                                                                const Numeric& dK3_dF0,
                                                                const Numeric& dK3_dTl,
                                                                const Numeric& dK3_dTu,
                                                                const Numeric& dK4_dT,
                                                                const Numeric& dK4_dTu)
{
  const Index nf = F.nelem();
  const Index nppd = derivatives_data_position.nelem();
  
  const Numeric invQT = 1.0 / QT;
  const Numeric QT_ratio = QT0 * invQT;
  
  const Numeric dS_dS0_abs = isotopic_ratio * QT_ratio * K1 * K2 * K3;
  const Numeric S_abs = S0 * dS_dS0_abs;
  const Numeric dS_dS0_src = isotopic_ratio * QT_ratio * K1 * K2 * K4;
  const Numeric S_src = S0 * dS_dS0_src;
  
  for(Index iq = 0; iq < nppd; iq++) {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature) {
      const Numeric dS_dT_abs = S_abs * (dK2_dT / K2 + dK1_dT / K1 + dK3_dT / K3 - invQT * dQT_dT);
      const Numeric dS_dT_src = S_src * (dK2_dT / K2 + dK1_dT / K1 + dK4_dT / K4 - invQT * dQT_dT);
      
      dN(iq, joker)  = dF(iq, joker);
      dN(iq, joker) *= S_src - S_abs;
      dF(iq, joker) *=         S_abs;
      for(Index iv = 0; iv < nf; iv++) {
        dN.get(iq, iv) += F.get(iv) * (dS_dT_src - dS_dT_abs);
        dF.get(iq, iv) += F.get(iv) *              dS_dT_abs;
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineStrength) {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        dF(iq, joker) = F;
        dN(iq, joker) = dF(iq, joker);
        dN(iq, joker) *= dS_dS0_src - dS_dS0_abs;
        dF(iq, joker) *=              dS_dS0_abs;
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter) {
      if(quantum_identity > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        const Numeric dS_dF0_abs = S_abs * (dK2_dF0 / K2 + dK3_dF0 / K3);
        const Numeric dS_dF0_src = S_src * (dK2_dF0 / K2);
        
        dN(iq, joker) = dF(iq, joker);
        dN(iq, joker) *= S_src - S_abs;
        dF(iq, joker) *=         S_abs;
        for(Index iv = 0; iv < nf; iv++) {
          dN.get(iq, iv) += F.get(iv) * (dS_dF0_src - dS_dF0_abs);
          dF.get(iq, iv) += F.get(iv) *               dS_dF0_abs;
        }
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::NLTE) {
      if(quantum_identity.LowerQuantumId() > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        const Numeric dS_dTl_abs = S_abs * dK3_dTl / K3;
        const Numeric dS_dTl_src = 0;
        for(Index iv = 0; iv < nf; iv++) {
          dN.get(iq, iv) = F.get(iv) * (dS_dTl_src - dS_dTl_abs);
          dF.get(iq, iv) = F.get(iv) *               dS_dTl_abs;
        }
      }
      else if(quantum_identity.UpperQuantumId() > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        const Numeric dS_dTu_abs = S_abs * dK3_dTu / K3;
        const Numeric dS_dTu_src = S_src * dK4_dTu / K4;
        for(Index iv = 0; iv < nf; iv++) {
          dN.get(iq, iv) = F.get(iv) * (dS_dTu_src - dS_dTu_abs);
          dF.get(iq, iv) = F.get(iv) *               dS_dTu_abs;
        }
      }
    }
    else {
      dN(iq, joker)  = dF(iq, joker);
      dN(iq, joker) *= S_src - S_abs;
      dF(iq, joker) *=         S_abs;
    }
  }
  
  // Set lineshape at the end
  N = F;
  N *= S_src - S_abs;
  F *=         S_abs;
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
 * \param dlfd Derivatives in order as they appear that are related to internal line parameters
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_linefunctiondata_jacobian_scaling(ComplexMatrixView dF,
                                                            const ArrayOfRetrievalQuantity& derivatives_data,
                                                            const ArrayOfIndex& derivatives_data_position,
                                                            const QuantumIdentifier& quantum_identity,
                                                            const LineRecord& line,
                                                            const Numeric& T,
                                                            const Numeric& P,
                                                            const Index& this_species,
                                                            const ConstVectorView& vmrs,
                                                            const ArrayOfArrayOfSpeciesTag& species)
{
  const Index nppd = derivatives_data_position.nelem();
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    const Index& pos = derivatives_data_position[iq];
    const RetrievalQuantity& rt = derivatives_data[pos];
    
    if(rt.QuantumIdentity().In(quantum_identity)) {
      switch(rt.PropMatType()) {
        case JacPropMatType::LineFunctionDataG0X0:
        case JacPropMatType::LineFunctionDataG0X1:
        case JacPropMatType::LineFunctionDataG0X2:
        case JacPropMatType::LineFunctionDataG2X0:
        case JacPropMatType::LineFunctionDataG2X1:
        case JacPropMatType::LineFunctionDataG2X2:
          dF(iq, joker) *= Complex(0, line.GetInternalDerivative(T,P,this_species,vmrs,species, rt));
          break;
        case JacPropMatType::LineFunctionDataD0X0:
        case JacPropMatType::LineFunctionDataD0X1:
        case JacPropMatType::LineFunctionDataD0X2:
        case JacPropMatType::LineFunctionDataD2X0:
        case JacPropMatType::LineFunctionDataD2X1:
        case JacPropMatType::LineFunctionDataD2X2:
          dF(iq, joker) *= -line.GetInternalDerivative(T,P,this_species,vmrs,species, rt);
          break;
        case JacPropMatType::LineFunctionDataFVCX0:
        case JacPropMatType::LineFunctionDataFVCX1:
        case JacPropMatType::LineFunctionDataFVCX2:
        case JacPropMatType::LineFunctionDataETAX0:
        case JacPropMatType::LineFunctionDataETAX1:
        case JacPropMatType::LineFunctionDataETAX2:
        case JacPropMatType::LineFunctionDataGX0:
        case JacPropMatType::LineFunctionDataGX1:
        case JacPropMatType::LineFunctionDataGX2:
        case JacPropMatType::LineFunctionDataDVX0:
        case JacPropMatType::LineFunctionDataDVX1:
        case JacPropMatType::LineFunctionDataDVX2:
          dF(iq, joker) *= line.GetInternalDerivative(T,P,this_species,vmrs,species, rt);
          break;
        case JacPropMatType::LineFunctionDataYX0:
        case JacPropMatType::LineFunctionDataYX1:
        case JacPropMatType::LineFunctionDataYX2:
          dF(iq, joker) *= Complex(0, -line.GetInternalDerivative(T,P,this_species,vmrs,species, rt));
          break;
        default:
        {/*pass*/}
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
                                                      const Numeric& partition_function_at_temperature,
                                                      const Numeric& dpartition_function_at_temperature_dT,
                                                      const Numeric& partition_function_at_line_temperature,
                                                      const ArrayOfArrayOfSpeciesTag& abs_species,
                                                      const Index& this_species_location_in_tags,
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
  
  // Pressure broadening and line mixing terms
  const auto X = line.GetShapeParams(temperature, pressure, this_species_location_in_tags, volume_mixing_ratio_of_all_species, abs_species);
  const auto& G0  = std::get<Index(LineFunctionData::TuplePos::G0 )>(X);
  const auto& D0  = std::get<Index(LineFunctionData::TuplePos::D0 )>(X);
  const auto& G2  = std::get<Index(LineFunctionData::TuplePos::G2 )>(X);
  const auto& D2  = std::get<Index(LineFunctionData::TuplePos::D2 )>(X);
  const auto& FVC = std::get<Index(LineFunctionData::TuplePos::FVC)>(X);
  const auto& ETA = std::get<Index(LineFunctionData::TuplePos::ETA)>(X);
  const auto& Y   = std::get<Index(LineFunctionData::TuplePos::Y  )>(X);
  const auto& G   = std::get<Index(LineFunctionData::TuplePos::G  )>(X);
  const auto& DV  = std::get<Index(LineFunctionData::TuplePos::DV )>(X);
  
  // Partial derivatives for temperature
  const LineFunctionData::Output dXdT = do_temperature ?  line.GetShapeParams_dT(temperature, temperature_perturbation(derivatives_data), 
                                                             pressure, this_species_location_in_tags, volume_mixing_ratio_of_all_species, 
                                                             abs_species) : NoLineFunctionDataOutput();
  const auto& dG0_dT  = std::get<Index(LineFunctionData::TuplePos::G0 )>(dXdT);
  const auto& dD0_dT  = std::get<Index(LineFunctionData::TuplePos::D0 )>(dXdT);
  const auto& dG2_dT  = std::get<Index(LineFunctionData::TuplePos::G2 )>(dXdT);
  const auto& dD2_dT  = std::get<Index(LineFunctionData::TuplePos::D2 )>(dXdT);
  const auto& dFVC_dT = std::get<Index(LineFunctionData::TuplePos::FVC)>(dXdT);
  const auto& dETA_dT = std::get<Index(LineFunctionData::TuplePos::ETA)>(dXdT);
  const auto& dY_dT   = std::get<Index(LineFunctionData::TuplePos::Y  )>(dXdT);
  const auto& dG_dT   = std::get<Index(LineFunctionData::TuplePos::G  )>(dXdT);
  const auto& dDV_dT  = std::get<Index(LineFunctionData::TuplePos::DV )>(dXdT);
  
  // Partial derivatives for VMR
  const std::tuple<bool, const QuantumIdentifier&> do_vmr = do_vmr_jacobian(derivatives_data, line.QuantumIdentity());
  const LineFunctionData::Output dXdVMR = std::get<0>(do_vmr) ?  line.GetShapeParams_dVMR(temperature, pressure, this_species_location_in_tags,
                                                                      volume_mixing_ratio_of_all_species, abs_species, std::get<1>(do_vmr)) : NoLineFunctionDataOutput();
                                                                      
  // Partial derivatives for temperature
  Numeric dK1_dT, dK2_dT;
  
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
  
  switch(line.GetExternalLineShapeType()) {
    // Use data as provided by the pressure broadening scheme
    case LineShapeType::ByPressureBroadeningData:
      switch(line.GetInternalLineShapeType()) {
        // Use data as per speed dependent air
        case LineFunctionData::LineShapeType::HTP:
        case LineFunctionData::LineShapeType::SDVP:
          lst = LineShapeType::HTP;
          set_htp(F, dF, 
                  f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), doppler_constant, 
                  G0, D0, G2, D2, ETA, FVC,
                  derivatives_data, derivatives_data_position, QI,
                  ddoppler_constant_dT, 
                  dG0_dT, dD0_dT, dG2_dT, dD2_dT, dETA_dT, dFVC_dT);
          break;
        case LineFunctionData::LineShapeType::VP:
          lst = LineShapeType::Voigt;
          set_voigt(F, dF, f_grid, 
                    line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    line.F(), doppler_constant, 
                    G0, D0, DV, derivatives_data, derivatives_data_position, QI,
                    ddoppler_constant_dT, dG0_dT, dD0_dT, dDV_dT, dXdVMR);
          break;
        case LineFunctionData::LineShapeType::LP:
          lst = LineShapeType::Lorentz;
          set_lorentz(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      line.F(), G0, D0, DV, derivatives_data, derivatives_data_position, QI, dG0_dT, dD0_dT, dDV_dT, dXdVMR);
          break;
        case LineFunctionData::LineShapeType::DP:
          lst = LineShapeType::Doppler;
          set_doppler(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      line.F(), doppler_constant, derivatives_data, derivatives_data_position, QI, ddoppler_constant_dT);
          break;
          
      } break;
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
              G0, D0, G2, D2, ETA, FVC,
              derivatives_data, derivatives_data_position, QI,
              ddoppler_constant_dT, 
              dG0_dT, dD0_dT, dG2_dT, dD2_dT, dETA_dT, dFVC_dT);
      lst = LineShapeType::HTP;
      break;
    // This line only needs Lorentz
    case LineShapeType::Lorentz:
      lst = LineShapeType::Lorentz;
      set_lorentz(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), G0, D0, DV, derivatives_data, derivatives_data_position, QI, dG0_dT, dD0_dT, dDV_dT, dXdVMR);
      break;
    // This line only needs Voigt
    case LineShapeType::Voigt:
      lst = LineShapeType::Voigt;
      set_voigt(F, dF, f_grid, 
                line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                line.F(), doppler_constant, 
                G0, D0, DV, derivatives_data, derivatives_data_position, QI,
                ddoppler_constant_dT, dG0_dT, dD0_dT, dDV_dT, dXdVMR);
      break;
    case LineShapeType::End:
      throw std::runtime_error("Cannot understand the requested line shape type.");
  }
  
  // Set the mirroring by repeating computations above using 
  // negative numbers for frequency of line related terms
  // The user sets if they want mirroring by LSM MTM followed by an index
  // that is interpreted as either mirroring by the same line shape or as 
  // mirroring by Lorentz lineshape
  switch(line.GetMirroringType()) {
    // No mirroring
    case MirroringType::None:
    case MirroringType::Manual:
      break;
    // Lorentz mirroring
    case MirroringType::Lorentz: {
        if(lst == LineShapeType::Doppler) throw std::runtime_error("Cannot apply Lorentz mirroring for Doppler line shape");
        // Set the mirroring computational vectors and size them as needed
        ComplexVector Fm(F.nelem());
        ComplexMatrix dFm(dF.nrows(), dF.ncols());
        
        set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    -line.F(), G0, -D0, -DV, derivatives_data, derivatives_data_position, QI, dG0_dT, -dD0_dT, -dDV_dT, dXdVMR);
        
        // Apply mirroring
        F -= Fm;
        dF -= dFm;
    } break;
    // Same type of mirroring as before
    case MirroringType::SameAsLineShape: {
      // Set the mirroring computational vectors and size them as needed
      ComplexVector Fm(F.nelem());
      ComplexMatrix dFm(dF.nrows(), dF.ncols());
      
      switch(lst) {
        case LineShapeType::Doppler:
          set_doppler(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      -line.F(), -doppler_constant, derivatives_data, derivatives_data_position, QI, -ddoppler_constant_dT);
          break;
        case LineShapeType::Lorentz:
          set_lorentz(Fm, dFm, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      -line.F(), G0, -D0, -DV, derivatives_data, derivatives_data_position, QI, dG0_dT, -dD0_dT, -dDV_dT, dXdVMR);
          break;
        case LineShapeType::Voigt:
          set_voigt(Fm, dFm, f_grid, 
                    -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    -line.F(), -doppler_constant, 
                    G0, -D0, -DV, derivatives_data, derivatives_data_position, QI,
                    -ddoppler_constant_dT, dG0_dT, -dD0_dT, -dDV_dT, dXdVMR);
          break;
        case LineShapeType::HTP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dFm, f_grid, 
                  -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  -line.F(), -doppler_constant, 
                  G0, -D0, G2, -D2, ETA, FVC,
                  derivatives_data, derivatives_data_position, QI,
                  -ddoppler_constant_dT, 
                  dG0_dT, -dD0_dT, dG2_dT, -dD2_dT, dETA_dT, dFVC_dT);
          break;
        case LineShapeType::ByPressureBroadeningData:
        case LineShapeType::End:
          throw std::runtime_error("Cannot understand the requested line shape type for mirroring.");
      }
      F -= Fm;
      dF -= dFm;
      break;
    } break;
    case MirroringType::End:
      throw std::runtime_error("Cannot understand the requested mirroring type for mirroring.");
  }
  
  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by and index
  // that is internally interpreted to mean some kind of lineshape normalization
  switch(line.GetLineNormalizationType()) {
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
  
  // Only for non-Doppler shapes
  if(lst not_eq LineShapeType::Doppler) {
    apply_linemixing_scaling(F, dF, Y, G, derivatives_data, derivatives_data_position, QI, dY_dT, dG_dT, dXdVMR);
  
    // Apply line mixing and pressure broadening partial derivatives        
    apply_linefunctiondata_jacobian_scaling(dF, derivatives_data, derivatives_data_position, QI, line,
                                            temperature, pressure,  this_species_location_in_tags,
                                            volume_mixing_ratio_of_all_species, abs_species);
  }
  
  // Apply line strength by whatever method is necessary
  switch(line.GetLinePopulationType()) {
    case LinePopulationType::ByLTE:
    case LinePopulationType::ByVibrationalTemperatures: {
      // Line strength scaling that are line-dependent ---
      // partition functions are species dependent and computed at a higher level
      const Numeric gamma = stimulated_emission(temperature, line.F());
      const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
      const Numeric K1 = boltzman_ratio(temperature, line.Ti0(), line.Elow());
      const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
      
      // Line strength partial derivatives
      
      if(do_temperature) {
        dK1_dT = dboltzman_ratio_dT(K1, temperature, line.Elow());
        dK2_dT = dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), temperature);
      }
      
      // Partial derivatives due to central frequency of the stimulated emission
      Numeric dK2_dF0;
      if(do_line_center)
        dK2_dF0 = dstimulated_relative_emission_dF0(gamma, gamma_ref, temperature, line.Ti0());
      
      // Multiply the line strength by the line shape
      if(line.GetLinePopulationType() == LinePopulationType::ByLTE)
        apply_linestrength_scaling(F, dF,  line.I0() * line.ZeemanEffect().StrengthScaling(zeeman_index), isotopologue_ratio,
                                   partition_function_at_temperature, partition_function_at_line_temperature, K1, K2,
                                   derivatives_data, derivatives_data_position, QI,
                                   dpartition_function_at_temperature_dT, dK1_dT, dK2_dT, dK2_dF0);
      else if (line.GetLinePopulationType() == LinePopulationType::ByVibrationalTemperatures) {
        // NLTE parameters
        Numeric Tu, Tl, K4, r_low, dK3_dF0, dK3_dT, dK3_dTl, dK4_dT, dK3_dTu, dK4_dTu;
        
        // These four are set by user on controlfile level
        // They are indexes to find the energy level in the nlte-temperature 
        // vector and the energy level of the states
        const Index evlow_index = line.NLTELowerIndex();
        const Index evupp_index = line.NLTEUpperIndex();
        const Numeric El = line.Evlow();
        const Numeric Eu = line.Evupp();
        
        // If the user set this parameters, another set of calculations are needed
        if(evupp_index > -1 and nlte_distribution.nelem() > evupp_index)
        {
          Tu = nlte_distribution[evupp_index];
          
          // Additional emission is from upper state
          K4 = boltzman_ratio(Tu, temperature, Eu);
        }
        else if(evupp_index > -1)
          throw std::runtime_error("Bad size nlte_distribution for declared NLTE calculations");
        // Otherwise the ratios are unity and nothing needs be done
        else {
          Tu = temperature;
          K4 = 1.0;
        }
        
        // The same as above but for the lower state level
        if(evlow_index > -1 and nlte_distribution.nelem() > evlow_index) {
          Tl = nlte_distribution[evlow_index];
          r_low = boltzman_ratio(Tl, temperature, El);
        }
        else if(evlow_index > -1)
          throw std::runtime_error("Bad size nlte_distribution for declared NLTE calculations");
        else {
          Tl = temperature;
          r_low = 1.0;
        }
        
        // Any additional absorption requires the ratio between upper and lower state number distributions
        const Numeric K3 = absorption_nlte_ratio(gamma, K4, r_low);
        
        // Are we computing the line center derivatives?
        if(do_line_center)
          dK3_dF0 = dabsorption_nlte_rate_dF0(gamma, temperature, K4, r_low);
        
        // Are we computing the temperature derivatives?
        // NOTE:  Having vibrational NLTE active AT ALL will change the jacobian because of this part of the code,
        // though this requires setting El and Eu for all lines, though this is not yet default...
        // So if you see this part of the code after having a runtime_error, 
        // you will need to write those functions yourself...
        if(do_temperature)
          dK3_dT = dabsorption_nlte_rate_dT(gamma, temperature, line.F(), El, Eu, K4, r_low);
        
        // Does the lower state level energy exist?
        if(El > 0)
          dK3_dTl = dabsorption_nlte_rate_dTl(gamma, temperature, Tl, El, r_low);
        
        // Does the upper state level energy exist?
        if(Eu > 0) {
          dK3_dTu = dabsorption_nlte_rate_dTu(gamma, temperature, Tu, Eu, K4);
          dK4_dTu = dboltzman_ratio_dT(K4, Tu, Eu);
        }
        
        // Apply this knowledge to set N and dN
        apply_linestrength_scaling_vibrational_nlte(F, dF, N, dN, line.I0() * line.ZeemanEffect().StrengthScaling(zeeman_index),
                                                    isotopologue_ratio, partition_function_at_temperature, 
                                                    partition_function_at_line_temperature, K1, K2, K3, K4,
                                                    derivatives_data, derivatives_data_position, QI,
                                                    dpartition_function_at_temperature_dT, dK1_dT, dK2_dT, dK2_dF0, dK3_dT, dK3_dF0, dK3_dTl,
                                                    dK3_dTu, dK4_dT, dK4_dTu);
      }
    } break;
    case LinePopulationType::ByPopulationDistribution: {
      const Index nlte_low_index = line.NLTELowerIndex();
      const Index nlte_upp_index = line.NLTEUpperIndex();
      
      if(nlte_low_index < 0 or nlte_distribution.nelem() <= nlte_low_index)
        throw std::runtime_error("No lower level distribution number in population distribution mode");
      if(nlte_upp_index < 0 or nlte_distribution.nelem() <= nlte_upp_index)
        throw std::runtime_error("No upper level distribution number in population distribution mode");
      
      apply_linestrength_from_nlte_level_distributions(F, dF, N, dN,
                                                       nlte_distribution[nlte_low_index],
                                                       nlte_distribution[nlte_upp_index],
                                                       line.G_lower(), line.G_upper(),
                                                       line.A(), line.F(),
                                                       temperature, derivatives_data, 
                                                       derivatives_data_position, QI);
      
    } break;
    case LinePopulationType::End: 
      throw std::runtime_error("Cannot understand the line strength computations");
  }
  
  // Cutoff frequency is applied at the end because 
  // the entire process above is complicated and applying
  // cutoff last means that the code is kept cleaner
  if(need_cutoff) {
    apply_cutoff(F, dF, N, dN,
                 derivatives_data, derivatives_data_position, line,
                 volume_mixing_ratio_of_all_species,
                 nlte_distribution, pressure, temperature,
                 doppler_constant, partial_pressure,
                 isotopologue_ratio, magnetic_magnitude,
                 ddoppler_constant_dT,
                 partition_function_at_temperature,
                 dpartition_function_at_temperature_dT,
                 partition_function_at_line_temperature,
                 abs_species, this_species_location_in_tags,
                 zeeman_index, verbosity);
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
                                 const Numeric& partition_function_at_temperature,
                                 const Numeric& dpartition_function_at_temperature_dT,
                                 const Numeric& partition_function_at_line_temperature,
                                 const ArrayOfArrayOfSpeciesTag& abs_species,
                                 const Index& this_species_location_in_tags,
                                 const Index& zeeman_index,
                                 const Verbosity& verbosity)
{ 
  // Size of derivatives
  const Index nj = dF.nrows(); 
  const Index nn = dN.nrows(); 
  
  // Setup compute variables
  
  Vector f_grid_cutoff(1);
  if(line.F() > 0) f_grid_cutoff[0] = line.F() + line.CutOff();
  else             f_grid_cutoff[0] = line.F() - line.CutOff();
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
                                    partition_function_at_temperature,
                                    dpartition_function_at_temperature_dT,
                                    partition_function_at_line_temperature,
                                    abs_species, this_species_location_in_tags,
                                    zeeman_index, verbosity, true);
  
  // Apply cutoff values
  F -= Fc.get(0);
  
  if(N.nelem())
    N -= Nc.get(0);
  for(Index i = 0; i < nj; i++)
    dF(i, joker) -= dFc.get(i, 0);
  for(Index i = 0; i < nn; i++)
    dN(i, joker) -= dNc.get(i, 0);
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
  
  const Numeric ratio = e/b - k;
  
  // Partial derivatives
  for(Index iq=0; iq<nppd; iq++) {
    if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::Temperature) { // nb.  Only here tou counter-act the latter on in
      Numeric done_over_b_dT = PLANCK_CONST*F0*exp_T/(c2*BOLTZMAN_CONST*T*T);
      
      // dN is unset so far.  It should return to just be lineshape later...
      for(Index iv=0; iv<nf; iv++)
        dN(iq, iv) = F[iv]*e*done_over_b_dT + dF(iq, iv)*ratio;
      
      // dk_dT = 0...
      dF(iq, joker) *= k;
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::LineCenter) {
      if(derivatives_data[derivatives_data_position[iq]].QuantumIdentity().In(quantum_identity)) {
        Numeric done_over_b_df0 = PLANCK_CONST*exp_T/(c2*BOLTZMAN_CONST*T) - 3.0*b/F0;
        Numeric de_df0 = c1 * r2 * A21;
        Numeric dk_df0 = c1 * (r1*x - r2) * (A21 / c2) - 3.0*k/F0;
        
        for(Index iv=0; iv<nf; iv++) {
          dN(iq, iv) = F[iv]*(e*done_over_b_df0 + de_df0/b - dk_df0) + dF(iq, iv)*ratio;
          dF(iq, iv) = dF(iq, iv)*k + F[iv]*dk_df0;
        }
      }
    }
    else if(derivatives_data[derivatives_data_position[iq]] == JacPropMatType::NLTE) {
      if(quantum_identity.LowerQuantumId() > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        const Numeric dk_dr2 = - c3 * A21 / c2, de_dr2 = c3 * A21, dratio_dr2 = de_dr2/b - dk_dr2;
        
        for(Index iv=0; iv<nf; iv++) {
          dN(iq, iv) = F[iv]*dratio_dr2;
          dF(iq, iv) = F[iv]*dk_dr2 + dF(iq, iv)*k;
        } 
      }
      else if(quantum_identity.UpperQuantumId() > derivatives_data[derivatives_data_position[iq]].QuantumIdentity()) {
        const Numeric dk_dr1 = c3 * x * A21 / c2;
        for(Index iv=0; iv<nf; iv++)
          dF(iq, iv) = F[iv]*dk_dr1 + dF(iq, iv)*k;
      }
    }
    else {
      dN(iq, joker)  = dF(iq, joker);
      dN(iq, joker) *= ratio;
      dF(iq, joker) *= k;
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
  while(p > Numeric((1 << i) * N)) i++;  // nb, N is the minimum number of points within range of p
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

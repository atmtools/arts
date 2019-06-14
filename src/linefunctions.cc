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
#include "constants.h"
#include <Eigen/Core>

// The Faddeeva function
inline Complex w(Complex z) noexcept {return Faddeeva::w(z);} 

void Linefunctions::set_lineshape(Eigen::Ref<Eigen::VectorXcd> F, 
                                  const Eigen::Ref<const Eigen::VectorXd> f_grid, 
                                  const LineRecord& line, 
                                  const ConstVectorView vmrs, 
                                  const Numeric& temperature, 
                                  const Numeric& pressure, 
                                  const Numeric& magnetic_magnitude,
                                  const ArrayOfArrayOfSpeciesTag& abs_species,
                                  const Index& this_species,
                                  const Index& zeeman_index)
{
  // Pressure broadening and line mixing terms
  const auto X = line.GetShapeParams(temperature, pressure, this_species, vmrs, abs_species);
  
  Eigen::MatrixXcd dF(0, 0), data(F.size(), Linefunctions::ExpectedDataSize());
  const Numeric doppler_constant = DopplerConstant(temperature, line.IsotopologueData().Mass());
  
  switch(line.GetLineShapeType()) {
    case LineFunctionData::LineShapeType::HTP:
    case LineFunctionData::LineShapeType::SDVP:
      set_htp(F, dF, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, X);
      break;
    case LineFunctionData::LineShapeType::VP:
      set_voigt(F, dF, data, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant, X);
      break;
    case LineFunctionData::LineShapeType::DP:
      set_doppler(F, dF, data, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), doppler_constant);
      break;
    case LineFunctionData::LineShapeType::LP:
      set_lorentz(F, dF, data, f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, line.F(), X);
      break;
  }
  
  switch(line.GetMirroringType()) {
    case MirroringType::None:
    case MirroringType::Manual:
      break;
    case MirroringType::Lorentz: {
      // Set the mirroring computational vectors and size them as needed
      Eigen::VectorXcd Fm(F.size());
      
      set_lorentz(Fm, dF, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), mirroredOutput(X));
      
      // Apply mirroring;  FIXME: Add conjugate?
      F.noalias() += Fm;
    }
    break;
    case MirroringType::SameAsLineShape: {
      // Set the mirroring computational vectors and size them as needed
      Eigen::VectorXcd Fm(F.size());
      
      switch(line.GetLineShapeType())
      {
        case LineFunctionData::LineShapeType::DP:
          set_doppler(Fm, dF, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant);
          break;
        case LineFunctionData::LineShapeType::LP:
          set_lorentz(Fm, dF, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), mirroredOutput(X));
          break;
        case LineFunctionData::LineShapeType::VP:
          set_voigt(Fm, dF, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant, mirroredOutput(X));
          break;
        case LineFunctionData::LineShapeType::HTP:
        case LineFunctionData::LineShapeType::SDVP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dF, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, -line.F(), -doppler_constant, mirroredOutput(X));
          break;
      }
      
      F.noalias() += Fm;
      break;
    }
  }
  
  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by a type
  // that is internally interpreted to mean some kind of lineshape normalization
  switch(line.GetLineNormalizationType()) {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::VVH:
      apply_VVH_scaling(F, dF, data, f_grid, line.F(), temperature);
      break;
    case LineNormalizationType::VVW:
      apply_VVW_scaling(F, dF, f_grid, line.F());
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F, dF, f_grid, line.F(), temperature);
      break;
  }
}



/*!
 * Sets the line shape to Lorentz line shape. Normalization is unity.
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param data Block of allocated memory
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param x Line shape parameters
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dxdT Temperature derivatives of line shape parameters
 * \param dxdVMR VMR derivatives of line shape parameters
 * 
 */
void Linefunctions::set_lorentz(Eigen::Ref<Eigen::VectorXcd> F,
                                Eigen::Ref<Eigen::MatrixXcd> dF,
                                Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                                const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const LineFunctionDataOutput& x,
                                const ArrayOfRetrievalQuantity& derivatives_data,
                                const ArrayOfIndex& derivatives_data_position,
                                const QuantumIdentifier& quantum_identity,
                                const LineFunctionDataOutput& dxdT,
                                const LineFunctionDataOutput& dxdVMR)
{
  // Size of the problem
  auto nppd = derivatives_data_position.nelem();
  
  // The central frequency
  const Numeric F0 = F0_noshift + x.D0 + zeeman_df * magnetic_magnitude + x.DV;
  
  // Naming req data blocks
  auto z  = data.col(0);
  auto dw = data.col(1);
  
  z.noalias() = (Constant::pi * Complex(x.G0, F0) - Complex(0, Constant::pi) * f_grid.array()).matrix();
  F.noalias() = z.cwiseInverse();
  
  if(nppd) {
    dw.noalias() = - Constant::pi * F.array().square().matrix();
    
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      if(deriv == JacPropMatType::Temperature)
        dF.col(iq).noalias() = Complex(dxdT.G0, dxdT.D0 + dxdT.DV) * dw;
      else if(is_frequency_parameter(deriv))
        dF.col(iq).noalias() = Complex(0.0, -1.0) * dw;
      else if((deriv == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(deriv)) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = Complex(0.0, 1.0) * dw;
      else if(is_pressure_broadening_parameter(deriv) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() =  Complex(0.0, -1.0) * dw;
      else if(is_magnetic_parameter(deriv))
        dF.col(iq).noalias() = Complex(0.0, zeeman_df) * dw;
      else if(deriv == JacPropMatType::VMR) {
        if(deriv.QuantumIdentity().In(quantum_identity))
          dF.col(iq).noalias() = Complex(dxdVMR.G0, dxdVMR.D0 + dxdVMR.DV) * dw;
        else
          dF.col(iq).setZero();
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
 * \param x Line shape parameters
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dxdT Temperature derivatives of line shape parameters
 * \param dxdVMR VMR derivatives of line shape parameters
 * 
 */
void Linefunctions::set_htp(Eigen::Ref<Eigen::VectorXcd> F,
                            Eigen::Ref<Eigen::MatrixXcd> dF,
                            const Eigen::Ref<const Eigen::VectorXd> f_grid, 
                            const Numeric& zeeman_df, 
                            const Numeric& magnetic_magnitude,
                            const Numeric& F0_noshift, 
                            const Numeric& GD_div_F0,
                            const LineFunctionDataOutput& x,
                            const ArrayOfRetrievalQuantity& derivatives_data,
                            const ArrayOfIndex& derivatives_data_position,
                            const QuantumIdentifier& quantum_identity,
                            const Numeric& dGD_div_F0_dT,
                            const LineFunctionDataOutput& dxdT,
                            const LineFunctionDataOutput& dxdVMR)
{
  auto nq = derivatives_data_position.nelem();
  auto nf = f_grid.size();
  
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric fac = Constant::sqrt_pi / GD;
  
  const Complex C0(x.G0, x.D0);
  const Complex C2(x.G2, x.D2);
  const Complex C0t = (1 - x.ETA) * (C0 - 1.5 * C2) + x.FVC;
  const Complex C2t = (1 - x.ETA) * C2;
  
  if(C2t == Complex(0, 0)) {
    const Complex Z0 = (Complex(0, F0) + C0t) / GD;
    for(auto i=0; i<nf; i++) {
      const Complex Zm = Complex(0, -f_grid[i] / GD) + Z0;
      const Complex wm = Faddeeva::w(Complex(0, 1) * Zm);
      const Complex A = fac * wm;
      const Complex G = 1.0 - (x.FVC - x.ETA * (C0 - 1.5 * C2)) * A;
      
      F[i] = Constant::inv_pi * A / G;
      
      if(nq) {
        const Complex dwm_divdZ = 2.0 * (Zm * wm - Constant::inv_sqrt_pi);
        for(auto iq=0; iq<nq; iq++) {
          const RetrievalQuantity& rt = derivatives_data[derivatives_data_position[iq]];
          if(rt == JacPropMatType::Temperature) {
            const Numeric dGD = dGD_div_F0_dT * F0;
            const Complex dC0 =  Complex(dxdT.G0, dxdT.D0);
            const Complex dC2 =  Complex(dxdT.G2, dxdT.D2);
            const Complex dC0t =  -(C0 - 1.5*C2)*dxdT.ETA - (x.ETA - 1.0)*(dC0 - 1.5*dC2) + dxdT.FVC;
            const Complex dZm =  ((Complex(0, f_grid[i]-F0) - C0t)*dGD + GD*dC0t)/(GD*GD);
            const Complex dwm = dwm_divdZ * dZm;
            const Complex dfac =  -fac*dGD/GD;
            const Complex dA =  fac*dwm + wm*dfac;
            const Complex dG =  ((C0 - 1.5*C2)*x.ETA - x.FVC)*dA + ((C0 - 1.5*C2)*dxdT.ETA + (dC0 - 1.5*dC2)*x.ETA - dxdT.FVC)*A;
            
            dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
          }
          else if(is_frequency_parameter(rt)) {
            const Complex dZm = Complex(0, -1 / GD);
            const Complex dwm = dwm_divdZ * dZm;
            const Complex dA = fac * dwm;
            const Complex dG = (x.ETA * (C0 - 1.5 * C2) - x.FVC) * dA;
            
            dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
          }
          else if(is_magnetic_parameter(rt)) {
            const Complex dZm =  zeeman_df*(-C0t + Complex(0, f_grid[i]))/(GD_div_F0*F0);
            const Complex dfac =  -Constant::sqrt_pi*zeeman_df/(GD_div_F0*(F0*F0));
            const Complex dwm = dwm_divdZ * dZm;
            const Complex dA =  fac*dwm + wm*dfac;
            const Complex dG =  (x.ETA*(C0 - 1.5*C2) - x.FVC)*dA;
            
            dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
          }
          else if(rt == JacPropMatType::VMR) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC0 =  Complex(dxdVMR.G0, dxdVMR.D0);
              const Complex dC2 =  Complex(dxdVMR.G2, dxdVMR.D2);
              const Complex dC0t =  -(C0 - 1.5*C2)*dxdVMR.ETA - (x.ETA - 1.0)*(dC0 - 1.5*dC2) + dxdVMR.FVC;
              const Complex dZm =  dC0t/GD;
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dA =  fac*dwm;
              const Complex dG =  ((C0 - 1.5*C2)*x.ETA - x.FVC)*dA + ((C0 - 1.5*C2)*dxdVMR.ETA + (dC0 - 1.5*dC2)*x.ETA - dxdVMR.FVC)*A;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
	    else if(not i)
	      dF.col(iq).setZero();
          }
          else if(rt == JacPropMatType::LineCenter) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dZm =  (-C0t + Complex(0, f_grid[i]))/((F0*F0)*GD_div_F0);
              const Complex dfac =  -Constant::sqrt_pi/((F0*F0)*GD_div_F0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dA =  fac*dwm + wm*dfac;
              const Complex dG =  (x.ETA*(C0 - 1.5*C2) - x.FVC)*dA;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_velocity_changing_collision_frequency(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Numeric dZm = 1 / GD;
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dA = fac * dwm;
              const Complex dG = - A - (x.FVC - x.ETA * (C0 - 1.5 * C2)) * dA;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_correlation(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC0t = - (C0 - 1.5 * C2);
              const Complex dZm = dC0t / GD;
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dA = fac * dwm;
              const Complex dG = (C0 - 1.5 * C2) * A - (x.FVC - x.ETA * (C0 - 1.5 * C2)) * dA;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_speed_dependent(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC2(0, -1);  // Note change of sign to fit with Voigt/Lorentz
              const Complex dC0t = (1 - x.ETA) * (- 1.5 * dC2);
              const Complex dZm = dC0t / GD;
              const Complex dwm = dwm_divdZ  * dZm;
              const Complex dA = fac * dwm;
              const Complex dG = - (x.FVC - x.ETA * (C0 - 1.5 * C2)) * dA - 1.5 * dC2 * A;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_speed_independent(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC0(0, -1);  // Note change of sign to fit with Voigt/Lorentz
              const Complex dC0t = (1 - x.ETA) * dC0;
              const Complex dZm = dC0t / GD;
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dA = fac * dwm;
              const Complex dG = - (x.FVC - x.ETA * (C0 - 1.5 * C2)) * dA + x.ETA * dC0 * A;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
        }
      }
    }
  }
  else {
    const Complex X0 = (Complex(0, -F0) + C0t) / C2t;
    const Complex sqrtY = 0.5 * GD / C2t;
    const Complex Y = sqrtY * sqrtY;
    for(auto i=0; i<nf; i++) {
      const Complex X = Complex(0, f_grid[i]) / C2t + X0;
      const Complex Z0 = std::sqrt(X + Y);
      const Complex Zm = Z0 - sqrtY;
      const Complex Zp = Z0 + sqrtY;
      const Complex wm = Faddeeva::w(Complex(0, 1) * Zm);
      const Complex wp = Faddeeva::w(Complex(0, 1) * Zp);
      const Complex A = fac * (wm - wp);
      const Complex B = x.ETA / (1 - x.ETA) * (-1.0 + Constant::sqrt_pi/(2.0*sqrtY) * ((1.0-Zm*Zm)*wm - (1.0-Zp*Zp)*wp));
      const Complex G = 1.0 - (x.FVC - x.ETA * (C0 - 1.5 * C2)) * A + B;
      
      F[i] = Constant::inv_pi * A / G;
      if(nq) {
        const Complex dwm_divdZ = 2.0 * (Zm * wm - Constant::inv_sqrt_pi);
        const Complex dwp_divdZ = 2.0 * (Zp * wp - Constant::inv_sqrt_pi);
        for(auto iq=0; iq<nq; iq++) {
          const RetrievalQuantity& rt = derivatives_data[derivatives_data_position[iq]];
          if(rt == JacPropMatType::Temperature) {
            const Numeric dGD = dGD_div_F0_dT * F0;
            const Complex dC0 =  Complex(dxdT.G0, dxdT.D0);
            const Complex dC2 =  Complex(dxdT.G2, dxdT.D2);
            const Complex dfac =  -Constant::sqrt_pi*dGD/GD/GD;
            const Complex dC0t =  -(C0 - 1.5*C2)*dxdT.ETA - (x.ETA - 1)*(dC0 - 1.5*dC2) + dxdT.FVC;
            const Complex dC2t =  -(x.ETA - 1)*dC2 - C2*dxdT.ETA;
            const Complex dY =  (C2t*dGD - GD*dC2t)*GD/(2*C2t*C2t*C2t);
            const Complex dX =  ((Complex(0, F0-f_grid[i]) - C0t)*dC2t + C2t*dC0t)/C2t/C2t;
            const Complex dZm =  -dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
            const Complex dZp =  dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
            const Complex dwm = dwm_divdZ * dZm;
            const Complex dwp = dwp_divdZ * dZp;
            const Complex dA =  (wm - wp)*dfac + (dwm - dwp)*fac;
            const Complex dB =  (2*(Constant::sqrt_pi*((Zm*Zm - 1.)*wm - (Zp*Zp - 1.)*wp) + 2*sqrtY)*(x.ETA - 1)*Y*dxdT.ETA - 2*(Constant::sqrt_pi*((Zm*Zm - 1.)*wm - (Zp*Zp - 1.)*wp) + 2*sqrtY)*x.ETA*Y*dxdT.ETA - Constant::sqrt_pi*(((Zm*Zm - 1.)*wm - (Zp*Zp - 1.)*wp)*dY + 2*(-(Zm*Zm - 1.)*dwm + (Zp*Zp - 1.)*dwp - 2*Zm*wm*dZm + 2*Zp*wp*dZp)*Y)*(x.ETA - 1)*x.ETA)/(4*(x.ETA - 1)*(x.ETA - 1)*Y*sqrtY);
            const Complex dG =  ((2*C0 - 3*C2)*x.ETA - 2*x.FVC)*dA/2. + ((2*C0 - 3*C2)*dxdT.ETA + (2*dC0 - 3*dC2)*x.ETA - 2*dxdT.FVC)*A/2. + dB;
            
            dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
          }
          else if(is_frequency_parameter(rt)) {
            const Complex dX =  Complex(0, 1)/C2t;
            const Complex dZm =  dX/(2*Z0);
            const Complex dZp =  dX/(2*Z0);
            const Complex dwm = dwm_divdZ * dZm;
            const Complex dwp = dwp_divdZ * dZp;
            const Complex dA =  (dwm - dwp)*fac;
            const Complex dB =  Constant::sqrt_pi*x.ETA*((Zm*Zm - 1.)*dwm - (Zp*Zp - 1.)*dwp + 2*Zm*wm*dZm - 2*Zp*wp*dZp)/(2*sqrtY*(x.ETA - 1));
            const Complex dG =  (x.ETA*(2*C0 - 3*C2) - 2*x.FVC)*dA/2. + dB;
            
            dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
          }
          else if(is_magnetic_parameter(rt)) {
            const Complex dfac =  -Constant::sqrt_pi*zeeman_df/(GD_div_F0*F0*F0);
            const Complex dY =  GD_div_F0*GD_div_F0*zeeman_df*F0/(2*C2t*C2t);
            const Complex dX =  -Complex(0, zeeman_df)/C2t;
            const Complex dZm =  -dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
            const Complex dZp =  dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
            const Complex dwm = dwm_divdZ * dZm;
            const Complex dwp = dwp_divdZ * dZp;
            const Complex dA =  (wm - wp)*dfac + (dwm - dwp)*fac;
            const Complex dB =  Constant::sqrt_pi*x.ETA*((-(Zm*Zm - 1.)*wm + (Zp*Zp - 1.)*wp)*dY + 2*((Zm*Zm - 1.)*dwm - (Zp*Zp - 1.)*dwp + 2*Zm*wm*dZm - 2*Zp*wp*dZp)*Y)/(4*(x.ETA - 1)*Y*sqrtY);
            const Complex dG =  (x.ETA*(2*C0 - 3*C2) - 2*x.FVC)*dA/2. + dB;
            
            dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
          }
          else if(rt == JacPropMatType::VMR) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC0 =  Complex(dxdVMR.G0, dxdVMR.D0);
              const Complex dC2 =  Complex(dxdVMR.G2, dxdVMR.D2);
              const Complex dC0t =  -(C0 - 1.5*C2)*dxdVMR.ETA - (x.ETA - 1)*(dC0 - 1.5*dC2) + dxdVMR.FVC;
              const Complex dC2t =  -(x.ETA - 1)*dC2 - C2*dxdVMR.ETA;
              const Complex dY =  -2*dC2t/C2t * Y;
              const Complex dX =  ((Complex(0, F0-f_grid[i]) - C0t)*dC2t + C2t*dC0t)/C2t/C2t;
              const Complex dZm =  -dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dZp =  dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dwp = dwp_divdZ * dZp;
              const Complex dA =  fac*(dwm - dwp);
              const Complex dB =  (2*(Constant::sqrt_pi*((Zm*Zm - 1.)*wm - (Zp*Zp - 1.)*wp) + 2*sqrtY)*(x.ETA - 1)*Y*dxdVMR.ETA - 2*(Constant::sqrt_pi*((Zm*Zm - 1.)*wm - (Zp*Zp - 1.)*wp) + 2*sqrtY)*x.ETA*Y*dxdVMR.ETA - Constant::sqrt_pi*(((Zm*Zm - 1.)*wm - (Zp*Zp - 1.)*wp)*dY + 2*(-(Zm*Zm - 1.)*dwm + (Zp*Zp - 1.)*dwp - 2*Zm*wm*dZm + 2*Zp*wp*dZp)*Y)*(x.ETA - 1)*x.ETA)/(4*(x.ETA - 1)*(x.ETA - 1)*Y*sqrtY);
              const Complex dG =  ((2*C0 - 3*C2)*x.ETA - 2*x.FVC)*dA/2. + ((2*C0 - 3*C2)*dxdVMR.ETA + (2*dC0 - 3*dC2)*x.ETA - 2*dxdVMR.FVC)*A/2. + dB;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
	    else if(not i)
	      dF.col(iq).setZero();
          }
          else if(rt == JacPropMatType::LineCenter) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dfac =  -Constant::sqrt_pi/(F0*F0*GD_div_F0);
              const Complex dY =  F0*GD_div_F0*GD_div_F0/(2*C2t*C2t);
              const Complex dX =  -Complex(0, 1)/C2t;
              const Complex dZm =  -dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dZp =  dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dwp = dwp_divdZ * dZp;
              const Complex dA =  (wm - wp)*dfac + (dwm - dwp)*fac;
              const Complex dB =  Constant::sqrt_pi*x.ETA*((-(Zm*Zm - 1.)*wm + (Zp*Zp - 1.)*wp)*dY + 2*((Zm*Zm - 1.)*dwm - (Zp*Zp - 1.)*dwp + 2*Zm*wm*dZm - 2*Zp*wp*dZp)*Y)/(4*(x.ETA - 1)*Y*sqrtY);
              const Complex dG =  (x.ETA*(2*C0 - 3*C2) - 2*x.FVC)*dA/2. + dB;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_velocity_changing_collision_frequency(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dX =  -1./(C2*(x.ETA - 1));
              const Complex dZm =  dX/(2*Z0);
              const Complex dZp =  dX/(2*Z0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dwp = dwp_divdZ * dZp;
              const Complex dA =  fac*(dwm - dwp);
              const Complex dB =  Constant::sqrt_pi*x.ETA*((Zm*Zm - 1.)*dwm - (Zp*Zp - 1.)*dwp + 2*Zm*wm*dZm - 2*Zp*wp*dZp)/(2*sqrtY*(x.ETA - 1));
              const Complex dG =  (x.ETA*(2*C0 - 3*C2) - 2*x.FVC)*dA/2. - A + dB;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_correlation(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Numeric dETA = 1;
              const Complex dY =  -(GD*GD)*dETA/(2*C2*C2*(x.ETA - 1)*(x.ETA - 1)*(x.ETA - 1));
              const Complex dX =  (x.FVC - Complex(0, F0 - f_grid[i]))*dETA/(C2*(x.ETA - 1)*(x.ETA - 1));
              const Complex dZm =  -dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dZp =  dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dwp = dwp_divdZ * dZp;
              const Complex dA =  fac*(dwm - dwp);
              const Complex dB =  (2*(Constant::sqrt_pi*(((Zm*Zm) - 1.0)*wm - ((Zp*Zp) - 1.0)*wp) + 2*sqrtY)*(x.ETA - 1)*Y*dETA - 2*(Constant::sqrt_pi*(((Zm*Zm) - 1.0)*wm - ((Zp*Zp) - 1.0)*wp) + 2*sqrtY)*x.ETA*Y*dETA - Constant::sqrt_pi*((((Zm*Zm) - 1.0)*wm - ((Zp*Zp) - 1.0)*wp)*dY + 2*(-((Zm*Zm) - 1.0)*dwm + ((Zp*Zp) - 1.0)*dwp - 2*Zm*wm*dZm + 2*Zp*wp*dZp)*Y)*(x.ETA - 1)*x.ETA)/(4*(x.ETA - 1)*(x.ETA - 1)*Y*sqrtY);
              const Complex dG =  (C0 - 1.5*C2)*A*dETA - (x.FVC - (C0 - 1.5*C2)*x.ETA)*dA + dB;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
              
            }
          }
          else if(is_pressure_broadening_speed_dependent(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC2(0, -1);  // Note change of sign to fit with Voigt/Lorentz (for speed independent variables)
              const Complex dY =  -(GD*GD)*dC2/(2*(x.ETA - 1)*(x.ETA - 1)*C2*C2*C2);
              const Complex dX =  (-C0*x.ETA + C0 - Complex(0, F0-f_grid[i]) + x.FVC)*dC2/((x.ETA - 1)*C2*C2);
              const Complex dZm =  -dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dZp =  dY/(2*sqrtY) + dX/(2*Z0) + dY/(2*Z0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dwp = dwp_divdZ * dZp;
              const Complex dA =  fac*(dwm - dwp);
              const Complex dB =  Constant::sqrt_pi*x.ETA*((-((Zm*Zm) - 1.)*wm + ((Zp*Zp) - 1.)*wp)*dY + 2*(((Zm*Zm) - 1.)*dwm - ((Zp*Zp) - 1.)*dwp + 2*Zm*wm*dZm - 2*Zp*wp*dZp)*Y)/(4*(x.ETA - 1)*Y*sqrtY);
              const Complex dG =  -3*x.ETA*A*dC2/2. + (x.ETA*(2*C0 - 3*C2) - 2*x.FVC)*dA/2. + dB;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
          else if(is_pressure_broadening_speed_independent(rt)) {
            if(rt.QuantumIdentity().In(quantum_identity)) {
              const Complex dC0(0, -1);  // Note change of sign to fit with Voigt/Lorentz
              const Complex dX =  dC0/C2;
              const Complex dZm =  dX/(2*Z0);
              const Complex dZp =  dX/(2*Z0);
              const Complex dwm = dwm_divdZ * dZm;
              const Complex dwp = dwp_divdZ * dZp;
              const Complex dA =  fac*(dwm - dwp);
              const Complex dB =  Constant::sqrt_pi*x.ETA*(((Zm*Zm) - 1.)*dwm - ((Zp*Zp) - 1.)*dwp + 2*Zm*wm*dZm - 2*Zp*wp*dZp)/(2*sqrtY*(x.ETA - 1));
              const Complex dG =  x.ETA*A*dC0 - (x.ETA*(3*C2 - 2*C0) + 2*x.FVC)*dA/2. + dB;
              
              dF(i, iq) = Constant::inv_pi*(-A*dG + G*dA)/(G*G);
            }
          }
        }
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
 * \param data Block of allocated memory
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param x Line shape parameters
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dxdT Temperature derivatives of line shape parameters
 * \param dxdVMR VMR derivatives of line shape parameters
 * 
 */
void Linefunctions::set_voigt(Eigen::Ref<Eigen::VectorXcd> F,
                              Eigen::Ref<Eigen::MatrixXcd> dF,
                              Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                              const Eigen::Ref<const Eigen::VectorXd> f_grid, 
                              const Numeric& zeeman_df, 
                              const Numeric& magnetic_magnitude,
                              const Numeric& F0_noshift, 
                              const Numeric& GD_div_F0,
                              const LineFunctionDataOutput& x,
                              const ArrayOfRetrievalQuantity& derivatives_data,
                              const ArrayOfIndex& derivatives_data_position,
                              const QuantumIdentifier& quantum_identity,
                              const Numeric& dGD_div_F0_dT,
                              const LineFunctionDataOutput& dxdT,
                              const LineFunctionDataOutput& dxdVMR)
{ 
  // Size of problem
  auto nppd = derivatives_data_position.nelem();
  
  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + x.D0 + x.DV;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0 - GD_div_F0 * (dxdT.D0 + dxdT.DV);
  
  // constant normalization factor for Voigt
  const Numeric fac = Constant::inv_sqrt_pi * invGD;
  
  // Naming req data blocks
  auto z  = data.col(0);
  auto dw = data.col(1);
  
  // Frequency grid
  z.noalias() = invGD * (Complex(-F0, x.G0) + f_grid.array()).matrix();
  
  // Line shape
  F.noalias() = fac * z.unaryExpr(&w);
  
  if(nppd) {
    dw.noalias() = 2 * (Complex(0, fac * Constant::inv_sqrt_pi) - z.cwiseProduct(F).array()).matrix();
    
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      if(is_frequency_parameter(deriv))
        dF.col(iq).noalias() = invGD * dw;
      else if(deriv == JacPropMatType::Temperature)
        dF.col(iq).noalias() = dw * Complex(-dxdT.D0 - dxdT.DV, dxdT.G0) * invGD - F * dGD_dT * invGD - dw.cwiseProduct(z) * dGD_dT * invGD;
      else if((deriv == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(deriv)) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = -F/F0 - dw * invGD - dw.cwiseProduct(z)/F0;
      else if(is_pressure_broadening_parameter(deriv) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = dw * invGD;
      else if(is_magnetic_parameter(deriv))
        dF.col(iq).noalias() = dw * (- zeeman_df * invGD);
      else if(deriv == JacPropMatType::VMR and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = dw * Complex(- dxdVMR.D0 - dxdVMR.DV, dxdVMR.G0) * invGD;
      else
        dF.col(iq).setZero();
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
 * \param data Block of allocated memory
 * \param f_grid Frequency grid of computations
 * \param zeeman_df Zeeman shift parameter for the line
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param F0_noshift Central frequency without any shifts
 * \param GD_div_F0 Frequency-independent part of the Doppler broadening
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 */
void Linefunctions::set_doppler(Eigen::Ref<Eigen::VectorXcd> F,
                                Eigen::Ref<Eigen::MatrixXcd> dF,
                                Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                                const Eigen::Ref<const Eigen::VectorXd> f_grid, 
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& GD_div_F0,
                                const ArrayOfRetrievalQuantity& derivatives_data,
                                const ArrayOfIndex& derivatives_data_position,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dGD_div_F0_dT)
{
  auto nppd = derivatives_data_position.nelem();
  
  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude;
  const Numeric invGD = 1.0 / (GD_div_F0 * F0);
  
  // Naming data blocks
  auto x = data.col(0);
  auto mx2 = data.col(1);
  
  x.noalias() = (f_grid.array() - F0).matrix() * invGD;
  mx2.noalias() = -data.col(0).cwiseAbs2();
  F.noalias() = invGD * Constant::inv_sqrt_pi * mx2.array().exp().matrix();
  
  for(auto iq=0; iq<nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];
    
    if(is_frequency_parameter(deriv))
      dF.col(iq).noalias() = - 2 * invGD * F.cwiseProduct(x);
    else if(deriv == JacPropMatType::Temperature)
      dF.col(iq).noalias() = - dGD_div_F0_dT * F0 * invGD * (2.0 * F.cwiseProduct(mx2) + F);
    else if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() = - F.cwiseProduct(mx2) * (2/F0) + F.cwiseProduct(x) * 2 * (invGD - 1/F0);
    else if(deriv == JacPropMatType::VMR)
      dF.col(iq).setZero();  // Must reset incase other lineshapes are mixed in
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
void Linefunctions::set_voigt_from_full_linemixing(Eigen::Ref<Eigen::VectorXcd> F,
                                                   Eigen::Ref<Eigen::MatrixXcd> dF,
                                                   const Eigen::Ref<const Eigen::VectorXd> f_grid, 
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
  
  auto nf = f_grid.size(), nppd = derivatives_data_position.nelem();
  
  // Doppler broadening and line center
  const Complex eigenvalue = eigenvalue_no_shift + L0;
  const Numeric F0 = eigenvalue.real() + L0;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // constant normalization factor for voigt
  const Numeric fac = Constant::inv_sqrt_pi * invGD;
  
  // Ratio of the Lorentz halfwidth to the Doppler halfwidth
  const Complex z0 = -eigenvalue * invGD;
  
  // frequency in units of Doppler
  for (auto iv=0; iv<nf; iv++) {
    dx = f_grid[iv] * invGD;
    z = z0 + dx;
    w = Faddeeva::w(z);
    
    F[iv] = fac * w;
    
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      if(iq==0)
        dw_over_dz = 2.0 * (z * w - Constant::inv_sqrt_pi);
      
      if(is_frequency_parameter(deriv))
        dF(iv, iq) = fac * dw_over_dz * invGD;
      else if(deriv == JacPropMatType::Temperature) {
        dz = (deigenvalue_dT - dL0_dT) - z * dGD_dT;
        
        dF(iv, iq) = -F[iv] * dGD_dT;
        dF(iv, iq) += fac * dw_over_dz * dz;
        dF(iv, iq) *= invGD;
      }
      else if((deriv == JacPropMatType::LineCenter or is_line_mixing_parameter(deriv)) and deriv.QuantumIdentity().In(quantum_identity))
      {
        dz = -z * GD_div_F0 - 1.0;
        
        dF(iv, iq) = -F[iv] * GD_div_F0;
        dF(iv, iq) += dw_over_dz * dz;
        dF(iv, iq) *= fac * invGD;
      }
      else if(is_pressure_broadening_parameter(deriv))
      {
        if(deriv.QuantumIdentity().In(quantum_identity))
          dF(iv, iq) = fac * dw_over_dz * ( Complex(-1.0, 1.0) * invGD);
        else if(not iv)
          dF.col(iq).setZero();
      }
    }
  }
}


/*!
 * Applies line mixing scaling to already set lineshape and line mirror
 * 
 * Equation: 
 *   with_mirroring:
 *     F := (1+G-iY) * F + (1+G+iY) * Fm
 *   else:
 *     F := (1+G-iY) * F
 * 
 * and appropriate derivatives
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * 
 * \param Fm Mirrored lineshape
 * \param dFm Mirrored lineshape derivative
 * \param X Variable with line mixing information
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dT Variable with line mixing temperature derivative information
 * \param dVMR Variable with line mixing VMR derivative information
 * 
 */
void Linefunctions::apply_linemixing_scaling_and_mirroring(Eigen::Ref<Eigen::VectorXcd> F,
                                                           Eigen::Ref<Eigen::MatrixXcd> dF,
                                                           const Eigen::Ref<Eigen::VectorXcd> Fm,
                                                           const Eigen::Ref<Eigen::MatrixXcd> dFm,
                                                           const LineFunctionDataOutput& X,
                                                           const bool with_mirroring,
                                                           const ArrayOfRetrievalQuantity& derivatives_data,
                                                           const ArrayOfIndex& derivatives_data_position,
                                                           const QuantumIdentifier& quantum_identity,
                                                           const LineFunctionDataOutput& dT,
                                                           const LineFunctionDataOutput& dVMR)
{
  auto nppd = derivatives_data_position.nelem();
  
  const Complex LM = Complex(1.0 + X.G, -X.Y);
  
  dF *= LM;
  if(with_mirroring) dF.noalias() += dFm * std::conj(LM);

  if(with_mirroring) {
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      if(deriv == JacPropMatType::Temperature) {
        const auto c = Complex(dT.G, -dT.Y);
        dF.col(iq).noalias() += F * c + Fm * std::conj(c);
      }
      else if(is_real_line_mixing_strength_parameter(deriv) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = F + Fm;
      else if(is_imag_line_mixing_strength_parameter(deriv) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = F - Fm;
      else if(deriv == JacPropMatType::VMR and deriv.QuantumIdentity().In(quantum_identity)) {
        const auto c = Complex(dVMR.G, -dVMR.Y);
        dF.col(iq).noalias() += F * c + Fm * std::conj(c);
      }
    }
  }
  else {
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      if(deriv == JacPropMatType::Temperature)
        dF.col(iq).noalias() += F * Complex(dT.G, -dT.Y);
      else if(is_line_mixing_line_strength_parameter(deriv) and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() = F;
      else if(deriv == JacPropMatType::VMR and deriv.QuantumIdentity().In(quantum_identity))
        dF.col(iq).noalias() += F * Complex(dVMR.G, -dVMR.Y);
    }
  }
  
  F *= LM;
  if(with_mirroring) F.noalias() += Fm * std::conj(LM);
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
void Linefunctions::apply_rosenkranz_quadratic_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                                                       Eigen::Ref<Eigen::MatrixXcd> dF,
                                                       const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                                       const Numeric& F0,
                                                       const Numeric& T,
                                                       const ArrayOfRetrievalQuantity& derivatives_data,
                                                       const ArrayOfIndex& derivatives_data_position,
                                                       const QuantumIdentifier& quantum_identity)
{
  auto nf = f_grid.size(), nppd = derivatives_data_position.nelem();
  
  const Numeric invF0 = 1.0/F0;
  const Numeric mafac = (Constant::h) / (2.0 * Constant::k * T) /
  std::sinh((Constant::h * F0) / (2.0 * Constant::k * T)) * invF0;
  
  Numeric dmafac_dT_div_fun = 0;
  if(do_temperature_jacobian(derivatives_data))
    dmafac_dT_div_fun = -(Constant::k*T - F0*Constant::h/
    (2.0*std::tanh(F0*Constant::h/(2.0*Constant::k*T))))/(Constant::k*T*T);
  
  Numeric fun;
  
  for (auto iv=0; iv<nf; iv++) {
    fun = mafac * (f_grid[iv] * f_grid[iv]);
    F[iv] *= fun;
    
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      dF(iv, iq) *= fun;
      if(deriv == JacPropMatType::Temperature)
        dF(iv, iq) += dmafac_dT_div_fun * F[iv];
      else if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity))
        dF(iv, iq) += (-invF0 - Constant::h/(2.0*Constant::k*T*std::tanh(F0*Constant::h/(2.0*Constant::k*T)))) * F[iv];
      else if(is_frequency_parameter(deriv))
        dF(iv, iq) += (2.0 / f_grid[iv]) * F[iv];
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
void Linefunctions::apply_VVH_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                                      Eigen::Ref<Eigen::MatrixXcd> dF,
                                      Eigen::Ref<Eigen::Matrix<Complex, Eigen::Dynamic, ExpectedDataSize()>> data,
                                      const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                      const Numeric& F0,
                                      const Numeric& T,
                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                      const ArrayOfIndex& derivatives_data_position,
                                      const QuantumIdentifier& quantum_identity)
{ 
  auto nppd = derivatives_data_position.nelem();
  
  // 2kT is constant for the loop
  const Numeric kT = 2.0 * Constant::k * T;
  const Numeric c1 = Constant::h / kT;
  
  // denominator is constant for the loop
  const Numeric tanh_f0part = std::tanh(c1 * F0);
  const Numeric denom = F0 * tanh_f0part;
  
  // Name the data
  auto tanh_f = data.col(0);
  auto ftanh_f = data.col(1);
  
  tanh_f.noalias() = (c1 * f_grid.array()).tanh().matrix();
  ftanh_f.noalias() = f_grid.cwiseProduct(tanh_f) / denom;
  F.array() *= ftanh_f.array();
  
  dF.array().colwise() *= ftanh_f.array();
  for(auto iq=0; iq<nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];
    
    if(deriv == JacPropMatType::Temperature)
      dF.col(iq).noalias() += -c1/T * ((denom - F0/tanh_f0part) * F - denom * ftanh_f.cwiseProduct(F) + f_grid.cwiseProduct(F).cwiseQuotient(tanh_f));
    else if(is_frequency_parameter(deriv))
      dF.col(iq).noalias() += F.cwiseQuotient(f_grid) + c1 * (F.cwiseQuotient(tanh_f) - F.cwiseProduct(tanh_f));
    else if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() += (-1.0/F0 + c1*tanh_f0part - c1/tanh_f0part) * F;
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
void Linefunctions::apply_VVW_scaling(Eigen::Ref<Eigen::VectorXcd> F,
                                      Eigen::Ref<Eigen::MatrixXcd> dF,
                                      const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                      const Numeric& F0,
                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                      const ArrayOfIndex& derivatives_data_position,
                                      const QuantumIdentifier& quantum_identity)
{
  auto nf = f_grid.size(), nppd = derivatives_data_position.nelem();
  
  // denominator is constant for the loop
  const Numeric invF0 = 1.0 / F0;
  const Numeric invF02 = invF0 * invF0;
  
  for(auto iv=0; iv<nf; iv++) {
    // Set the factor
    const Numeric fac = f_grid[iv] * f_grid[iv] * invF02;
    
    // Set the line shape
    F[iv] *= fac;
    
    for(auto iq=0; iq<nppd; iq++) {
      const auto& deriv = derivatives_data[derivatives_data_position[iq]];
      
      // The factor is applied to all partial derivatives
      dF(iv, iq) *= fac;
      
      // These partial derivatives are special
      if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity))
        dF(iv, iq) -= 2.0 * invF0 * F[iv];
      else if(is_frequency_parameter(deriv))
        dF(iv, iq) += 2.0 / f_grid[iv] * F[iv];
    }
  }
}


Numeric Linefunctions::lte_linestrength(Numeric S0, Numeric E0, Numeric F0, Numeric QT0, Numeric T0, Numeric QT, Numeric T)
{
  const Numeric gamma = stimulated_emission(T, F0);
  const Numeric gamma_ref = stimulated_emission(T0, F0);
  const Numeric K1 = boltzman_ratio(T, T0, E0);
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
  return S0 * K1 * K2 * QT0/QT;
}


/*!
 * Applies linestrength to already set line shape by LTE population type
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * \retval N Source lineshape
 * \retval dN Source lineshape derivative
 * 
 * \param line The line data
 * \param T The atmospheric temperature
 * \param isotopic_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param zeeman_scaling A line strength scaling for the Zeeman effect
 * \param QT Partition function at atmospheric temperature of level
 * \param QT0 Partition function at reference temperature
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dQT_dT Temperature derivative of QT
 * 
 */
void Linefunctions::apply_linestrength_scaling_by_lte(Eigen::Ref<Eigen::VectorXcd> F,
                                                      Eigen::Ref<Eigen::MatrixXcd> dF,
                                                      Eigen::Ref<Eigen::VectorXcd> N,
                                                      Eigen::Ref<Eigen::MatrixXcd> dN,
                                                      const LineRecord& line,
                                                      const Numeric& T,
                                                      const Numeric& isotopic_ratio,
                                                      const Numeric& zeeman_scaling,
                                                      const Numeric& QT,
                                                      const Numeric& QT0,
                                                      const ArrayOfRetrievalQuantity& derivatives_data,
                                                      const ArrayOfIndex& derivatives_data_position,
                                                      const QuantumIdentifier& quantum_identity,
                                                      const Numeric& dQT_dT)
{
  auto nppd = derivatives_data_position.nelem();
  
  const Numeric gamma = stimulated_emission(T, line.F());
  const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
  const Numeric K1 = boltzman_ratio(T, line.Ti0(), line.Elow());
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
  
  const Numeric invQT = 1.0/QT;
  const Numeric S = zeeman_scaling * line.I0() * isotopic_ratio * QT0 * invQT * K1 * K2;
  
  F *= S;
  dF *= S;
  for(auto iq=0; iq<nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];
    
    if(deriv == JacPropMatType::Temperature)
      dF.col(iq).noalias() += F * (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T) / K2
                                 + dboltzman_ratio_dT(K1, T, line.Elow()) / K1
                                 - invQT * dQT_dT);
    else if(deriv == JacPropMatType::LineStrength and deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() = F / line.I0();  //nb. overwrite
    else if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity))
      dF.col(iq).noalias() += F * dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0()) / K2;
  }
  
  // No NLTE variables
  N.setZero();
  dN.setZero();
}


/*!
 * Applies linestrength to already set line shape by vibrational level temperatures
 * 
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * \retval N Source lineshape
 * \retval dN Source Lineshape derivative
 * 
 * \param line The absorption line
 * \param T The atmospheric temperature
 * \param Tu The upper state vibrational temperature; must be T if level is LTE
 * \param Tl The lower state vibrational temperature; must be T if level is LTE
 * \param Evu The upper state vibrational energy; if set funny, yields funny results
 * \param Evl The lower state vibrational energy; if set funny, yields funny results
 * \param isotopic_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param zeeman_scaling A line strength scaling for the Zeeman effect
 * \param QT Partition function at atmospheric temperature of level
 * \param QT0 Partition function at reference temperature
 * \param derivatives_data Information about the derivatives in dF
 * \param derivatives_data_position Information about the derivatives positions in dF
 * \param quantum_identity ID of the absorption line
 * \param dQT_dT Temperature derivative of QT
 * 
 */
void Linefunctions::apply_linestrength_scaling_by_vibrational_nlte(Eigen::Ref<Eigen::VectorXcd> F,
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
                                                                   const ArrayOfRetrievalQuantity& derivatives_data,
                                                                   const ArrayOfIndex& derivatives_data_position,
                                                                   const QuantumIdentifier& quantum_identity,
                                                                   const Numeric& dQT_dT)
{
  auto nppd = derivatives_data_position.nelem();
  
  const Numeric gamma = stimulated_emission(T, line.F());
  const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
  const Numeric r_low = boltzman_ratio(Tl, T, Evl);
  const Numeric r_upp = boltzman_ratio(Tu, T, Evu);
  
  const Numeric K1 = boltzman_ratio(T, line.Ti0(), line.Elow());
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
  const Numeric K3 = absorption_nlte_ratio(gamma, r_upp, r_low);
  const Numeric K4 = r_upp;
  
  const Numeric invQT = 1.0 / QT;
  const Numeric QT_ratio = QT0 * invQT;
  
  const Numeric dS_dS0_abs = isotopic_ratio * QT_ratio * K1 * K2 * K3;
  const Numeric S_abs = zeeman_scaling * line.I0() * dS_dS0_abs;
  const Numeric dS_dS0_src = isotopic_ratio * QT_ratio * K1 * K2 * K4;
  const Numeric S_src = zeeman_scaling * line.I0() * dS_dS0_src;
  
  dN.noalias() = dF * (S_src - S_abs);
  dF *= S_abs;
  for(auto iq=0; iq<nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];
    
    if(deriv == JacPropMatType::Temperature) {
      const Numeric dS_dT_abs = S_abs * (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T) / K2
                                       + dboltzman_ratio_dT(K1, T, line.Elow()) / K1
                                       + dabsorption_nlte_rate_dT(gamma, T, line.F(), Evl, Evu, K4, r_low) / K3
                                       - invQT * dQT_dT);
      const Numeric dS_dT_src = S_src * (dstimulated_relative_emission_dT(gamma, gamma_ref, line.F(), T) / K2
                                       + dboltzman_ratio_dT(K1, T, line.Elow()) / K1
                                       - dboltzman_ratio_dT(K4, T, Evu) / K4
                                       - invQT * dQT_dT);
      
      dN.col(iq).noalias() += F * (dS_dT_src - dS_dT_abs);
      dF.col(iq).noalias() += F * dS_dT_abs;
    }
    else if(deriv == JacPropMatType::LineStrength and deriv.QuantumIdentity().In(quantum_identity)) {
      dF.col(iq).noalias() = F * dS_dS0_abs;
      dN.col(iq).noalias() = F * (dS_dS0_src - dS_dS0_abs);
    }
    else if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity)) {
      const Numeric dS_dF0_abs = S_abs * (dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0()) / K2
                                        + dabsorption_nlte_rate_dF0(gamma, T, K4, r_low) / K3);
      const Numeric dS_dF0_src = S_src * dstimulated_relative_emission_dF0(gamma, gamma_ref, T, line.Ti0()) / K2;
      
      dN.col(iq).noalias() += F * (dS_dF0_src - dS_dF0_abs);
      dF.col(iq).noalias() += F * dS_dF0_abs;
    }
    else if(deriv == JacPropMatType::NLTE) {
      if(deriv.QuantumIdentity().InLower(quantum_identity)) {
        const Numeric dS_dTl_abs = S_abs * dabsorption_nlte_rate_dTl(gamma, T, Tl, Evl, r_low) / K3;
        
        dN.col(iq).noalias() = - F * dS_dTl_abs;
        dF.col(iq).noalias() = - dN.col(iq);
      }
      else if(deriv.QuantumIdentity().InUpper(quantum_identity)) {
        const Numeric dS_dTu_abs = S_abs * dabsorption_nlte_rate_dTu(gamma, T, Tu, Evu, K4) / K3;
        const Numeric dS_dTu_src = S_src * dboltzman_ratio_dT(K4, Tu, Evu) / K4;
        
        dN.col(iq).noalias() = F * (dS_dTu_src - dS_dTu_abs);
        dF.col(iq).noalias() = F * dS_dTu_abs;
      }
    }
  }
  
  N.noalias() = F * (S_src - S_abs);
  F *= S_abs;
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
void Linefunctions::apply_linestrength_from_full_linemixing(Eigen::Ref<Eigen::VectorXcd> F,
                                                            Eigen::Ref<Eigen::MatrixXcd> dF,
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
  
  constexpr Numeric C1 = - Constant::h / Constant::k;
  
  const Numeric invT = 1.0 / T, 
  F0_invT = F0 * invT,
  exp_factor = std::exp(C1 * F0_invT), 
  f0_factor = F0 * (1.0 - exp_factor); 
  
  const Complex S = S_LM * f0_factor * isotopic_ratio;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];
    
    if(deriv == JacPropMatType::Temperature)
    {
      dF.col(iq) *= S_LM;
      dF.col(iq).noalias() += F * (dS_LM_dT * f0_factor + 
      S_LM * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio;
    }
    else if(deriv == JacPropMatType::LineStrength)
    {
      if(deriv.QuantumIdentity().In(quantum_identity))
        throw std::runtime_error("Not working yet");
    }
    else if(deriv == JacPropMatType::LineCenter or is_line_mixing_DF_parameter(deriv))
    {
      if(deriv.QuantumIdentity().In(quantum_identity))
        throw std::runtime_error("Not working yet");
    }
    else
    {
      dF.col(iq) *= S;
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
void Linefunctions::apply_dipole(Eigen::Ref<Eigen::VectorXcd> F,
                                 Eigen::Ref<Eigen::MatrixXcd> /*dF*/,
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
  
  constexpr Numeric C1 = - Constant::h / Constant::k;
  
  const Index nppd = derivatives_data_position.nelem();
  
  const Numeric S = d0 * d0 * rho * isotopic_ratio, 
  invT = 1.0 / T, 
  F0_invT = F0 * invT,
  exp_factor = exp(C1 * F0_invT), 
  f0_factor = F0 * (1.0 - exp_factor);
  
  if(nppd)
    throw std::runtime_error("Cannot support Jacobian from dipole calculations yet");
  
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
void Linefunctions::apply_linefunctiondata_jacobian_scaling(Eigen::Ref<Eigen::MatrixXcd> dF,
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
  auto nppd = derivatives_data_position.nelem();
  
  for(auto iq=0; iq<nppd; iq++) {
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
          dF.col(iq) *= Complex(0, line.GetInternalDerivative(T,P,this_species,vmrs,species, rt));
          break;
        case JacPropMatType::LineFunctionDataD0X0:
        case JacPropMatType::LineFunctionDataD0X1:
        case JacPropMatType::LineFunctionDataD0X2:
        case JacPropMatType::LineFunctionDataD2X0:
        case JacPropMatType::LineFunctionDataD2X1:
        case JacPropMatType::LineFunctionDataD2X2:
          dF.col(iq) *= -line.GetInternalDerivative(T,P,this_species,vmrs,species, rt);
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
          dF.col(iq) *= line.GetInternalDerivative(T,P,this_species,vmrs,species, rt);
          break;
        case JacPropMatType::LineFunctionDataYX0:
        case JacPropMatType::LineFunctionDataYX1:
        case JacPropMatType::LineFunctionDataYX2:
          dF.col(iq) *= Complex(0, -line.GetInternalDerivative(T,P,this_species,vmrs,species, rt));
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
Numeric Linefunctions::DopplerConstant(Numeric T, Numeric mass)
{
  return std::sqrt(Constant::doppler_broadening_const_squared * T / mass);
}


/*! Returns the temperature derivative of the frequency-independent part of the Doppler broadening
 * 
 * \param T Atmospheric temperature at level
 * \param dc Output of Linefunctions::DopplerConstant(T, mass)
 * 
 */
Numeric Linefunctions::dDopplerConstant_dT(const Numeric& T, const Numeric& dc)
{
  return dc / (2 * T);
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
 * \param partition_function_at_temperature As name suggests
 * \param dpartition_function_at_temperature_dT Temeperature derivative of partition_function_at_temperature
 * \param partition_function_at_line_temperature As name suggests
 * \param broad_spec_locations Locations of broadening species using all-planetary broadening inside volume_mixing_ratio_of_all_species
 * \param this_species_location_in_tags Location of species of line in volume_mixing_ratio_of_all_species
 * \param water_index_location_in_tags Location of water in volume_mixing_ratio_of_all_species
 * \param verbosity Verbosity level
 * \param cutoff_call Flag to ignore some functions inside if this call is from the cutoff-computations
 */
void Linefunctions::set_cross_section_for_single_line(Eigen::Ref<Eigen::VectorXcd>  F_full,
                                                      Eigen::Ref<Eigen::MatrixXcd> dF_full,
                                                      Eigen::Ref<Eigen::VectorXcd>  N_full,
                                                      Eigen::Ref<Eigen::MatrixXcd> dN_full,
                                                      Eigen::Ref<Eigen::MatrixXcd> data_block_full,
                                                      Index& start_cutoff,
                                                      Index& nelem_cutoff,
                                                      const Eigen::Ref<const Eigen::VectorXd> f_grid_full,
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
  if(not cutoff_call)
    find_cutoff_ranges(start_cutoff, nelem_cutoff, f_grid_full, line.F(), cutoff);
  else {
    start_cutoff = 0;
    nelem_cutoff = 1;
  }

  const bool need_cutoff = cutoff > 0 and not cutoff_call;
  const bool do_temperature = do_temperature_jacobian(derivatives_data);
                           
  // Leave this function if there is nothing to compute
  if(nelem_cutoff == 0)
    return;
  
  // Extract the quantum identify of the line to be used in-case there are derivatives
  const QuantumIdentifier& QI = line.QuantumIdentity();
  
  // Pressure broadening and line mixing terms
  const auto X = line.GetShapeParams(temperature, pressure, this_species_location_in_tags, volume_mixing_ratio_of_all_species, abs_species);
  
  // Partial derivatives for temperature
  const auto dXdT = do_temperature ? line.GetShapeParams_dT(temperature, temperature_perturbation(derivatives_data), 
                                                            pressure, this_species_location_in_tags, 
                                                            volume_mixing_ratio_of_all_species,  abs_species)
                                   : NoLineFunctionDataOutput();
  
  // Partial derivatives for VMR... the first function 
  auto do_vmr = do_vmr_jacobian(derivatives_data, line.QuantumIdentity());  // At all, Species
  const auto dXdVMR = do_vmr.test ? line.GetShapeParams_dVMR(temperature, pressure, this_species_location_in_tags,
                                                             volume_mixing_ratio_of_all_species, abs_species, do_vmr.qid)
                                  : NoLineFunctionDataOutput();
  
  // Arrays on which all computations happen are segments of the full input, and the segmenting is part of the output
  auto  F = F_full.segment(start_cutoff, nelem_cutoff);
  auto  N = N_full.segment(start_cutoff, nelem_cutoff);
  auto dF = dF_full.middleRows(start_cutoff, nelem_cutoff);
  auto dN = dN_full.middleRows(start_cutoff, nelem_cutoff);
  auto data = data_block_full.middleRows(start_cutoff, nelem_cutoff);
  auto f_grid = f_grid_full.middleRows(start_cutoff, nelem_cutoff);
  
  switch(line.GetLineShapeType()) {
    case LineFunctionData::LineShapeType::DP:
      set_doppler(F, dF, data, 
                  f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), doppler_constant, derivatives_data, derivatives_data_position, QI, ddoppler_constant_dT);
      break;
    case LineFunctionData::LineShapeType::HTP:
    case LineFunctionData::LineShapeType::SDVP:
      set_htp(F, dF, 
              f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
              line.F(), doppler_constant, X, derivatives_data, derivatives_data_position, QI,
              ddoppler_constant_dT, dXdT, dXdVMR);
      break;
    case LineFunctionData::LineShapeType::LP:
      set_lorentz(F, dF, data,
                  f_grid, line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  line.F(), X, derivatives_data, derivatives_data_position, QI, dXdT, dXdVMR);
      break;
    case LineFunctionData::LineShapeType::VP:
      set_voigt(F, dF, data, f_grid, 
                line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                line.F(), doppler_constant, X, derivatives_data, derivatives_data_position, QI,
                ddoppler_constant_dT, dXdT, dXdVMR);
      break;
  }
  
  // Set the mirroring by repeating computations above using 
  // negative numbers for frequency of line related terms
  // The user sets if they want mirroring by LSM MTM followed by an index
  // that is interpreted as either mirroring by the same line shape or as 
  // mirroring by Lorentz lineshape
  const bool with_mirroring = line.GetMirroringType() not_eq MirroringType::None and 
                              line.GetMirroringType() not_eq MirroringType::Manual;
  switch(line.GetMirroringType()) {
    case MirroringType::None:
    case MirroringType::Manual:
      break;
    case MirroringType::Lorentz:
      set_lorentz(N, dN, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  -line.F(), mirroredOutput(X), derivatives_data, derivatives_data_position, QI, mirroredOutput(dXdT), mirroredOutput(dXdVMR));
      break;
    case MirroringType::SameAsLineShape:
      switch(line.GetLineShapeType()) {
        case LineFunctionData::LineShapeType::DP:
          set_doppler(N, dN, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      -line.F(), -doppler_constant, derivatives_data, derivatives_data_position, QI, -ddoppler_constant_dT);
          break;
        case LineFunctionData::LineShapeType::LP:
          set_lorentz(N, dN, data, f_grid, -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                      -line.F(), mirroredOutput(X), derivatives_data, derivatives_data_position, QI, mirroredOutput(dXdT), mirroredOutput(dXdVMR));
          break;
        case LineFunctionData::LineShapeType::VP:
          set_voigt(N, dN, data, f_grid, 
                    -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                    -line.F(), -doppler_constant, mirroredOutput(X), derivatives_data, derivatives_data_position, QI,
                    -ddoppler_constant_dT, mirroredOutput(dXdT), mirroredOutput(dXdVMR));
          break;
        case LineFunctionData::LineShapeType::HTP:
        case LineFunctionData::LineShapeType::SDVP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(N, dN, f_grid, 
                  -line.ZeemanEffect().SplittingConstant(zeeman_index), magnetic_magnitude, 
                  -line.F(), -doppler_constant, mirroredOutput(X),
                  derivatives_data, derivatives_data_position, QI,
                  -ddoppler_constant_dT, mirroredOutput(dXdT), mirroredOutput(dXdVMR));
          break;
      }
      break;
  }
  
  // Mixing and mirroring can only apply to non-Doppler shapes
  if(line.GetLineShapeType() not_eq LineFunctionData::LineShapeType::DP) {
    apply_linemixing_scaling_and_mirroring(F, dF, N, dN, X, with_mirroring,
                                           derivatives_data, derivatives_data_position, QI, dXdT, dXdVMR);
    
    // Apply line mixing and pressure broadening partial derivatives        
    apply_linefunctiondata_jacobian_scaling(dF, derivatives_data, derivatives_data_position, QI, line,
                                            temperature, pressure,  this_species_location_in_tags,
                                            volume_mixing_ratio_of_all_species, abs_species);
  }
  
  // Line normalization if necessary
  // The user sets this by setting LSM LNT followed by an index
  // that is internally interpreted to mean some kind of lineshape normalization
  switch(line.GetLineNormalizationType()) {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::VVH:
      apply_VVH_scaling(F, dF, data, f_grid, line.F(), temperature, derivatives_data, derivatives_data_position, QI);
      break;
    case LineNormalizationType::VVW:
      apply_VVW_scaling(F, dF, f_grid, line.F(), derivatives_data, derivatives_data_position, QI);
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic_scaling(F, dF, f_grid, line.F(), temperature, derivatives_data, derivatives_data_position, QI);
      break;
  }
  
  // Apply line strength by whatever method is necessary
  switch(line.GetLinePopulationType()) {
    case LinePopulationType::ByLTE:
      apply_linestrength_scaling_by_lte(F, dF, N, dN, line, temperature, isotopologue_ratio, line.ZeemanEffect().StrengthScaling(zeeman_index),
                                        partition_function_at_temperature, partition_function_at_line_temperature,
                                        derivatives_data, derivatives_data_position, QI, dpartition_function_at_temperature_dT);
      break;
    case LinePopulationType::ByVibrationalTemperatures:
      if(line.NLTELowerIndex() >= 0 and nlte_distribution.nelem() <= line.NLTELowerIndex() )
        throw std::runtime_error("Bad lower vibrational temperature record.  It contains fewer temperatures than LineRecord require");
      
      if(line.NLTEUpperIndex() >= 0 and nlte_distribution.nelem() <= line.NLTEUpperIndex() )
        throw std::runtime_error("Bad upper vibrational temperature record.  It contains fewer temperatures than LineRecord require");
    
      apply_linestrength_scaling_by_vibrational_nlte(F, dF, N, dN, line, temperature,
                                                     line.NLTEUpperIndex() < 0 ? temperature : nlte_distribution[line.NLTEUpperIndex()],
                                                     line.NLTELowerIndex() < 0 ? temperature : nlte_distribution[line.NLTELowerIndex()],
                                                     line.Evupp() < 0 ? 0 : line.Evupp(),
                                                     line.Evlow() < 0 ? 0 : line.Evlow(),
                                                     isotopologue_ratio, line.ZeemanEffect().StrengthScaling(zeeman_index),
                                                     partition_function_at_temperature, partition_function_at_line_temperature,
                                                     derivatives_data, derivatives_data_position, QI, dpartition_function_at_temperature_dT);
      break;
    case LinePopulationType::ByPopulationDistribution:
      if(line.NLTELowerIndex() < 0)
        throw std::runtime_error("No lower NLTE distribution data for line marked to require it");
      if(nlte_distribution.nelem() <= line.NLTELowerIndex() )
        throw std::runtime_error("Bad lower NLTE distribution record.  It contains fewer values than LineRecord require");
      if(nlte_distribution[line.NLTELowerIndex()] <= 0)
        throw std::runtime_error("Bad lower NLTE distribution number.  It should be strictly above 0.");
      
      if(line.NLTEUpperIndex() < 0)
        throw std::runtime_error("No upper NLTE distribution data for line marked to require it");
      if(nlte_distribution.nelem() <= line.NLTEUpperIndex() )
        throw std::runtime_error("Bad upper NLTE distribution record.  It contains fewer values than LineRecord require");
      if(nlte_distribution[line.NLTEUpperIndex()] <= 0)
        throw std::runtime_error("Bad upper NLTE distribution number.  It should be strictly above 0.");
      
      apply_linestrength_from_nlte_level_distributions(F, dF, N, dN,
                                                       nlte_distribution[line.NLTELowerIndex()], nlte_distribution[line.NLTEUpperIndex()],
                                                       line.G_lower(), line.G_upper(), line.A(), line.F(), temperature, derivatives_data, 
                                                       derivatives_data_position, QI);
      break;
  }
  
  // Cutoff frequency is applied at the end because 
  // the entire process above is complicated and applying
  // cutoff last means that the code is kept cleaner

  if(need_cutoff)
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
                 zeeman_index);
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
void Linefunctions::apply_cutoff(Eigen::Ref<Eigen::VectorXcd> F,
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
                                 const Index& zeeman_index)
{ 
  // Size of derivatives
  auto nj = dF.cols(); 
  
  // Setup compute variables

  const auto v = Vector(1, (line.F() > 0) ? line.F() + line.CutOff() : line.F() - line.CutOff());
  const auto f_grid_cutoff = MapToEigen(v);
  Eigen::VectorXcd Fc(1), Nc(1);
  Eigen::RowVectorXcd dFc(nj), dNc(nj), data(Linefunctions::ExpectedDataSize());
  Index _tmp1, _tmp2;
  
  // Recompute the line for a single frequency
  set_cross_section_for_single_line(Fc, dFc, Nc, dNc, data, _tmp1, _tmp2, f_grid_cutoff, line,
                                    derivatives_data, derivatives_data_position,
                                    volume_mixing_ratio_of_all_species,
                                    nlte_distribution, pressure, temperature,
                                    doppler_constant, partial_pressure, isotopologue_ratio,
                                    magnetic_magnitude, ddoppler_constant_dT,
                                    partition_function_at_temperature,
                                    dpartition_function_at_temperature_dT,
                                    partition_function_at_line_temperature,
                                    abs_species, this_species_location_in_tags,
                                    zeeman_index, true);
  
  // Apply cutoff values
  F.array() -= Fc[0];
  N.array() -= Nc[0];
  for(Index i = 0; i < nj; i++) dF.col(i).array() -= dFc[i];
  for(Index i = 0; i < nj; i++) dN.col(i).array() -= dNc[i];
}


void Linefunctions::find_cutoff_ranges(Index& start_cutoff,
                                       Index& nelem_cutoff,
                                       const Eigen::Ref<const Eigen::VectorXd> f_grid,
                                       const Numeric& F0,
                                       const Numeric& cutoff)
{
  auto nf = f_grid.size();
  
  const bool need_cutoff = (cutoff > 0);
  if(need_cutoff) {
    // Find range of simulations
    start_cutoff = 0;
    Index i_f_max = nf-1;
    
    // Loop over positions to compute the line shape cutoff point
    while(start_cutoff < nf and (F0 - cutoff) > f_grid[start_cutoff])  ++start_cutoff;
    while(i_f_max >= start_cutoff and (F0 + cutoff) < f_grid[i_f_max]) --i_f_max;
    
    //  The extent is one more than the difference between the indices of interest
    nelem_cutoff = i_f_max - start_cutoff + 1; // min is 0, max is nf
  }
  else {
    start_cutoff = 0;
    nelem_cutoff = nf;
  }
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
void Linefunctions::apply_linestrength_from_nlte_level_distributions(Eigen::Ref<Eigen::VectorXcd> F,
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
                                                                     const ArrayOfRetrievalQuantity& derivatives_data,
                                                                     const ArrayOfIndex& derivatives_data_position,
                                                                     const QuantumIdentifier& quantum_identity)
{
  // Size of the problem
  auto nppd = derivatives_data_position.nelem();
  
  // Physical constants
  constexpr Numeric c0 = 2.0 * Constant::h / Constant::pow2(Constant::c);
  constexpr Numeric c1 = Constant::h / (4 * Constant::pi);
  
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
  const Numeric exp_T = std::exp(Constant::h / Constant::k * F0 / T),  b = c2/(exp_T - 1);
  
  // Absorption strength
  const Numeric k = c3 * (r1*x - r2) * (A21 / c2);
  
  // Emission strength
  const Numeric e = c3 * r2 * A21;
  
  // Ratio between emission and absorption constant  
  const Numeric ratio = e/b - k;
  
  // Constants ALMOST everywhere inside these loops
  dN.noalias() = dF * ratio;
  dF *= k;

  // Partial derivatives
  for(auto iq=0; iq<nppd; iq++) {
    const auto& deriv = derivatives_data[derivatives_data_position[iq]];
    
    if(deriv == JacPropMatType::Temperature)
      dN.col(iq).noalias() += F*e*Constant::h*F0*exp_T/(c2*Constant::k*T*T);
    else if(deriv == JacPropMatType::LineCenter and deriv.QuantumIdentity().In(quantum_identity)) {
      Numeric done_over_b_df0 = Constant::h*exp_T/(c2*Constant::k*T) - 3.0*b/F0;
      Numeric de_df0 = c1 * r2 * A21;
      Numeric dk_df0 = c1 * (r1*x - r2) * (A21 / c2) - 3.0*k/F0;
      
      dN.col(iq).noalias() += F*(e*done_over_b_df0 + de_df0/b - dk_df0);
      dF.col(iq).noalias() += F*dk_df0;
    }
    else if(deriv == JacPropMatType::NLTE) {
      if(deriv.QuantumIdentity().InLower(quantum_identity)) {
        const Numeric dk_dr2 = - c3 * A21 / c2, de_dr2 = c3 * A21, dratio_dr2 = de_dr2/b - dk_dr2;
        dN.col(iq).noalias() = F*dratio_dr2;
        dF.col(iq).noalias() += F*dk_dr2;
      }
      else if(deriv.QuantumIdentity().InUpper(quantum_identity)) {
        const Numeric dk_dr1 = c3 * x * A21 / c2;
        dF.col(iq).noalias() += F*dk_dr1;
      }
    }
  }
  
  // Set source function to be the relative amount emitted by the line divided by the Planck function
  N.noalias() = F * ratio;
  
  // Set absorption
  F *= k;
}


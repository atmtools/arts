/* Copyright (C) 2000-2012
 * Axel von Engeln <ric.larsson@gmail.com>
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
 * This file contains both the lineshape functions themselves and the
 * function define_lineshape_data which sets the lineshape lookup
 * data. 
 * 
 * This is the file from arts-1-0, back-ported to arts-1-1.
 * 
 * \author Richard Larsson
 * \date   2017-05-16
 */

#include "Faddeeva.hh"
#include "lineshapesdata.h"

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

void Linefunctions::set_lorentz(ComplexVector& F, // Sets the full complex line shape without line mixing
                                ArrayOfComplexVector& dF,
                                const Vector& f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& G0,
                                const Numeric& L0,
                                const PropmatPartialsData& derivatives_data,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dG0_dT,
                                const Numeric& dL0_dT)
{ 
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude;
  
  // Signa change of F and F0?
  const Complex denom0 = Complex(G0, F0);
  
  F = invPI;
  
  for(Index iv = 0; iv < nf; iv++)
  {
    F[iv] /= (denom0 - Complex(0.0, f_grid[iv]));
  }
  
  if(nppd > 0)
  {
    dF[nppd-1] = F;
    dF[nppd-1] *= F;
    dF[nppd-1] *= PI;
  }
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    
    if(derivatives_data(iq) == JQT_temperature)
    {
      dF[iq] = dF[nppd-1];
      dF[iq] *= Complex(dG0_dT, dL0_dT);
    }
    else if(derivatives_data(iq) == JQT_frequency or
      derivatives_data(iq) == JQT_wind_magnitude or
      derivatives_data(iq) == JQT_wind_u or
      derivatives_data(iq) == JQT_wind_v or
      derivatives_data(iq) == JQT_wind_w)
    {
      dF[iq] = dF[nppd-1];
      dF[iq] *= Complex(0.0, -1.0);
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq] = dF[nppd-1];
        dF[iq] *= Complex(0.0, 1.0);
      }
    }
    else if(derivatives_data(iq) == JQT_line_gamma_self or 
      derivatives_data(iq) == JQT_line_gamma_selfexponent or
      derivatives_data(iq) == JQT_line_pressureshift_self or
      derivatives_data(iq) == JQT_line_gamma_foreign or
      derivatives_data(iq) == JQT_line_gamma_foreignexponent or
      derivatives_data(iq) == JQT_line_pressureshift_foreign or
      derivatives_data(iq) == JQT_line_gamma_water or
      derivatives_data(iq) == JQT_line_gamma_waterexponent or 
      derivatives_data(iq) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq] = dF[nppd-1];
        dF[iq] *= Complex(1.0, 1.0);
      }
    }
    else if(derivatives_data(iq) == JQT_magnetic_magntitude)// No //external inputs --- errors because of frequency shift when Zeeman is used?
    {
      dF[iq] = dF[nppd-1];
      dF[iq] *= Complex(0.0, zeeman_df);
    }
  }
}


void Linefunctions::set_mirrored_lorentz(ComplexVector& F, // Sets the full complex line shape without line mixing
                                         ArrayOfComplexVector& dF, 
                                         const Vector& f_grid,
                                         const Numeric& zeeman_df,
                                         const Numeric& magnetic_magnitude,
                                         const Numeric& F0_noshift,
                                         const Numeric& G0,
                                         const Numeric& L0,
                                         const PropmatPartialsData& derivatives_data,
                                         const QuantumIdentifier& quantum_identity,
                                         const Numeric& dG0_dT,
                                         const Numeric& dL0_dT)
{ 
  Complex Fplus, Fminus;
  
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude;
  
  // Signa change of F and F0?
  const Complex denom0 = Complex(G0, F0), denom1 = Complex(G0, -F0);
  
  for(Index iv = 0; iv < nf; iv++)
  {
    Fplus = invPI / (denom0 - Complex(0.0, f_grid[iv]));
    Fminus = invPI / (denom1 - Complex(0.0, f_grid[iv]));
    F[iv] = Fplus + Fminus;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(derivatives_data(iq) == JQT_temperature)
      {
        dF[iq][iv] = (Fplus * Fplus * Complex(dG0_dT, dL0_dT) + Fminus * Fminus * Complex(dG0_dT, -dL0_dT)) * PI;
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
      {
        dF[iq][iv] = (Fplus * Fplus + Fminus * Fminus) * PI;
        dF[iq][iv] *= Complex(0.0, -1.0);
      }
      else if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][iv] = (Fplus * Fplus * Complex(0.0, 1.0) + Fminus * Fminus * Complex(0.0, -1.0)) * PI;
        }
      }
      else if(derivatives_data(iq) == JQT_line_gamma_self or 
        derivatives_data(iq) == JQT_line_gamma_selfexponent or
        derivatives_data(iq) == JQT_line_pressureshift_self or
        derivatives_data(iq) == JQT_line_gamma_foreign or
        derivatives_data(iq) == JQT_line_gamma_foreignexponent or
        derivatives_data(iq) == JQT_line_pressureshift_foreign or
        derivatives_data(iq) == JQT_line_gamma_water or
        derivatives_data(iq) == JQT_line_gamma_waterexponent or 
        derivatives_data(iq) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][iv] = (Fplus * Fplus * Complex(1.0, 1.0) + Fminus * Fminus * Complex(1.0, -1.0)) * PI;
        }
      }
      else if(derivatives_data(iq) == JQT_magnetic_magntitude)// No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        dF[iq][iv] = (Fplus * Fplus * Complex(0.0, zeeman_df) + Fminus * Fminus * Complex(0.0, -zeeman_df)) * PI;
      }
    }
  }
}


void Linefunctions::set_htp(ComplexVector& F, // Sets the full complex line shape without line mixing
                            ArrayOfComplexVector& dF,
                            const Vector& f_grid,
                            const Numeric& zeeman_df,
                            const Numeric& magnetic_magnitude,
                            const Numeric& F0_noshift,
                            const Numeric& GD_div_F0,
                            const Numeric& G0,
                            const Numeric& L0,
                            const Numeric& L2,
                            const Numeric& G2,
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
                            const Numeric& dFVC_dT)
{
  /*
   * 
   * Function is meant to compute the Hartmann-Tran lineshape
   * 
   * The assumptions on input are described in accompanying pdf-documents
   * 
   * Note that Hartmann-Tran lineshape does not support line mixing in this version
   * 
   */
  //extern const Numeric PI;
  static const Complex i(0.0, 1.0), one_plus_one_i(1.0, 1.0);
  static const Numeric ln2 = log(2.0), sqrtLN2 = sqrt(ln2); 
  
  // Main lineshape parameters
  Complex A, B, Zp, Zm, Zp2, Zm2, wiZp, wiZm, X, sqrtXY, invG;
  
  // Derivatives parameters
  Complex dA, dB, dC0, dC0t, dC2, dC2t, dZp, dZm, dwiZp, dwiZm, dX, dY, dg;
  
  // Number of frequencies
  const Index nf = f_grid.nelem();
  
  // Number of derivatives
  const Index nppd = derivatives_data.nelem();
  
  // Zeeman effect and other line center shifts means that the true line center is shifted
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude;
  const Numeric F0_ratio = F0_noshift / F0;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric one_minus_eta = 1.0 - eta;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // Index of derivative for first occurence of similarly natured partial derivatives
  const Index first_frequency = derivatives_data.get_first_frequency(), 
  first_pressure_broadening = derivatives_data.get_first_pressure_term();
  
  // Pressure broadening terms
  const Complex C0 = G0 + i*L0;
  const Complex C2 = G2 + i*L2;
  
  // Pressure terms adjusted by collisions and correlations
  const Complex C0_m1p5_C2 = (C0 - 1.5 * C2);
  const Complex C0t = one_minus_eta * C0_m1p5_C2 + FVC;
  const Complex invC2t = 1.0 / (one_minus_eta * C2);
  
  // Relative pressure broadening terms to be used in the shape function
  const Complex sqrtY = 0.5 * invC2t * GD / sqrtLN2;
  const Complex Y = sqrtY * sqrtY;
  
  // Scale factor
  const Numeric const1 = ((sqrtLN2 * sqrtPI) / GD);
  
  // Loop over frequencies that cannot be parallelized
  for(Index iv = 0; iv < nf; iv++)
  {
    // Relative frequency
    X = (C0t - i*(F0 - f_grid[iv])) * invC2t;
    
    // Setting up the Z terms
    sqrtXY = sqrt(X+Y);
    Zm = Zp = sqrtXY;
    Zm -= sqrtY;
    Zp += sqrtY;
    
    // The shape functions
    wiZm = Faddeeva::w(i*Zm);
    wiZp = Faddeeva::w(i*Zp);
    
    // Main line shape is computed here --- WARNING when wiZM and wiZp are too close, this could lead to numerical errors --- should be caught before frequency loop
    A = const1 * (wiZm - wiZp);
    
    // If there are correlations or if there is a temperature dependency
    if(eta not_eq 0.0 or deta_dT not_eq 0.0)
    {
      Zm2 = Zm * Zm;
      Zp2 = Zp * Zp;
      B = (1.0 / one_minus_eta) * (-1.0 + 0.5 * sqrtPI / sqrtY * ((1.0 - Zm2)*wiZm - (1.0 - Zp2)*wiZp));
      invG = 1.0 / (1.0 - (FVC - eta * C0_m1p5_C2)*A + eta * B);
    }
    else
    {
      invG = 1.0 / (1.0 - FVC*A);  // WARNING if FVC x A == 1 then this fail --- no way around it, though A should be very small and FVC should not be too large
    }
    F[iv] = A * invG * invPI;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w) // No //external inputs
      { 
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          dX = invC2t * i;
          dZp = dZm = 0.5 * dX / sqrtXY;
          
          dwiZm = dwiZp = i * sqrtInvPI;
          dwiZm -= Zm * wiZm;
          dwiZp -= Zp * wiZp;
          dwiZm *= 2.0 * dZm;
          dwiZp *= 2.0 * dZp;
          
          dA = const1 * (dwiZm - dwiZp);
          
          if(eta not_eq 0.0)
          {
            dB = 0.5 * sqrtPI / (sqrtY * one_minus_eta) * ((1.0 - Zm2) * dwiZm - 2.0 * Zm * dZm * wiZm + (1.0 - Zp2) * dwiZp + 2.0 * Zp * dZp * wiZp);
            dg = eta * (C0_m1p5_C2 * dA - dB) - FVC * dA;
          }
          else
          {
            dg = - FVC * dA;
          }
          
          dF[iq][iv] = invG * (invPI * dA - F[iv] * dg); 
        }
        else  // copy for repeated occurences
        {
          dF[iq][iv] = dF[first_frequency][iv]; 
        }
      }
      else if(derivatives_data(iq) == JQT_temperature)
      {
        dC0 = dG0_dT + i*dL0_dT;
        dC2 = dG2_dT + i*dL2_dT;
        
        dC0t = one_minus_eta * (dC0 - 1.5 * dC2) - deta_dT * C0_m1p5_C2 + dFVC_dT;
        dC2t = one_minus_eta * dC2 - deta_dT * C2;
        
        dY = GD / (2.0 * ln2) * invC2t*invC2t * (dGD_dT - GD * invC2t * dC2t);
        dX = invC2t * (dC0t - X*dC2t);
        
        dZm = dZp = 0.5 * (dX + dY) / sqrtXY;
        dZm -= 0.5 / sqrtY * dY;
        dZp += 0.5 / sqrtY * dY;
        
        dwiZm = dwiZp = i * sqrtInvPI;
        dwiZm -= Zm * wiZm;
        dwiZp -= Zp * wiZp;
        dwiZm *= 2.0 * dZm;
        dwiZp *= 2.0 * dZp;
        
        dA = const1 * (dwiZm - dwiZp - A * GD_div_F0);
        
        if(eta not_eq 0)
        {
          dB = 1.0 / one_minus_eta * (1.0 / one_minus_eta * deta_dT + 0.5 * sqrtPI / sqrtY * ( - 0.5 * ((1.0 - Zm2) * wiZm - (1.0 - Zp2) * wiZp) / Y * dY + (1.0 - Zm2) * dwiZm - 2.0 * Zm * dZm * wiZm + (1.0 - Zp2) * dwiZp + 2.0 * Zp * dZp * wiZp));
          dg = (dFVC_dT - deta_dT * C0_m1p5_C2 - eta * (dC0 - 1.5 * dC2)) * A - (FVC - eta * C0_m1p5_C2) * dA + deta_dT * B + eta * dB;
        }
        else
        {
          dg = (dFVC_dT - deta_dT * C0_m1p5_C2) * A - FVC * dA + deta_dT * B;
        }
        
        dF[iq][iv] = invG * (invPI * dA - F[iv] * dg); 
      }
      else if(derivatives_data(iq) == JQT_line_center) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dY = GD / (2.0 * ln2) * invC2t*invC2t * GD_div_F0 * F0_ratio;
          dX = - i * invC2t * F0_ratio;
          
          dZm = dZp = 0.5 * (dX + dY) / sqrtXY;
          dZm -= 0.5 / sqrtY * dY;
          dZp += 0.5 / sqrtY * dY;
          
          dwiZm = dwiZp = i * sqrtInvPI;
          dwiZm -= Zm * wiZm;
          dwiZp -= Zp * wiZp;
          dwiZm *= 2.0 * dZm;
          dwiZp *= 2.0 * dZp;
          
          dA = const1 * (dwiZm - dwiZp - A * GD_div_F0);
          
          if(eta not_eq 0.0)
          {
            dB = 0.5 * sqrtPI / (sqrtY * one_minus_eta) * (- 0.5 * ((1.0 - Zm2) * wiZm - (1.0 - Zp2) * wiZp) / Y * dY + (1.0 - Zm2) * dwiZm - 2.0 * Zm * dZm * wiZm + (1.0 - Zp2) * dwiZp + 2.0 * Zp * dZp * wiZp);
            dg = eta * (C0_m1p5_C2 * dA - dB) - FVC * dA;
          }
          else
          {
            dg = - FVC * dA;
          }
          
          dF[iq][iv] = invG * (invPI * dA - F[iv] * dg); 
        }
      }
      else if(derivatives_data(iq) == JQT_line_gamma_self or 
        derivatives_data(iq) == JQT_line_gamma_selfexponent or
        derivatives_data(iq) == JQT_line_pressureshift_self or
        derivatives_data(iq) == JQT_line_gamma_foreign or
        derivatives_data(iq) == JQT_line_gamma_foreignexponent or
        derivatives_data(iq) == JQT_line_pressureshift_foreign or
        derivatives_data(iq) == JQT_line_gamma_water or
        derivatives_data(iq) == JQT_line_gamma_waterexponent or 
        derivatives_data(iq) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == iq)
          {
            dC0t = one_minus_eta * one_plus_one_i;
            dX = invC2t * dC0t;
            dZm = dZp = 0.5 * dX / sqrtXY;
            
            dwiZm = dwiZp = i * sqrtInvPI;
            dwiZm -= Zm * wiZm;
            dwiZp -= Zp * wiZp;
            dwiZm *= 2.0 * dZm;
            dwiZp *= 2.0 * dZp;
            
            dA = const1 * (dwiZm - dwiZp);
            
            if(eta not_eq 0.0)
            {
              dB = 0.5 * sqrtPI / (sqrtY * one_minus_eta) * ((1.0 - Zm2) * dwiZm - 2.0 * Zm * dZm * wiZm + (1.0 - Zp2) * dwiZp + 2.0 * Zp * dZp * wiZp);
              dg = eta * (C0_m1p5_C2 * dA - dB - dC0 * A) - FVC * dA;
            }
            else
            {
              dg = - FVC * dA;
            }
            
            dF[iq][iv] = invG * (invPI * dA - F[iv] * dg); 
          }
          else  // copy for repeated occurences
          {
            dF[iq][iv] = dF[first_frequency][iv]; 
          }
        }
      }
      else if(derivatives_data(iq) == JQT_magnetic_magntitude)// No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        dY = GD / (2.0 * ln2) * invC2t*invC2t * GD_div_F0 * (1.0 - F0_ratio) / magnetic_magnitude;
        dX = - i * invC2t * (1.0 - F0_ratio) / magnetic_magnitude;
        
        dZm = dZp = 0.5 * (dX + dY) / sqrtXY;
        dZm -= 0.5 / sqrtY * dY;
        dZp += 0.5 / sqrtY * dY;
        
        dwiZm = dwiZp = i * sqrtInvPI;
        dwiZm -= Zm * wiZm;
        dwiZp -= Zp * wiZp;
        dwiZm *= 2.0 * dZm;
        dwiZp *= 2.0 * dZp;
        
        dA = const1 * (dwiZm - dwiZp - A * GD_div_F0);
        
        if(eta not_eq 0.0)
        {
          dB = 0.5 * sqrtPI / (sqrtY * one_minus_eta) * 
          (- 0.5 * ((1.0 - Zm2) * wiZm - (1.0 - Zp2) * wiZp) / Y * dY +
          (1.0 - Zm2) * dwiZm - 2.0 * Zm * dZm * wiZm +
          (1.0 - Zp2) * dwiZp + 2.0 * Zp * dZp * wiZp);
          dg = eta * (C0_m1p5_C2 * dA - dB) - FVC * dA;
        }
        else
        {
          dg = - FVC * dA;
        }
        
        dF[iq][iv] = invG * (invPI * dA - F[iv] * dg); 
      }
    }
  }
}


void Linefunctions::set_faddeeva_algorithm916(ComplexVector& F, 
                                              ArrayOfComplexVector& dF, 
                                              const Vector& f_grid, 
                                              const Numeric& zeeman_df, 
                                              const Numeric& magnetic_magnitude,
                                              const Numeric& F0_noshift, 
                                              const Numeric& GD_div_F0,
                                              const Numeric& G0, 
                                              const Numeric& L0,
                                              const PropmatPartialsData& derivatives_data,
                                              const QuantumIdentifier& quantum_identity,
                                              const Numeric& dGD_div_F0_dT,
                                              const Numeric& dG0_dT,
                                              const Numeric& dL0_dT)
{
  // For calculations
  Numeric dx;
  Complex w, z, dw_over_dz, dz;
  
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + L0;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // constant normalization factor for voigt
  const Numeric fac = sqrtInvPI * invGD;
  F = fac;
  
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
    
    F[iv] *= w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq==0)
        dw_over_dz = 2.0 * (z * w - sqrtInvPI);
      
      if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w) // No //external inputs
      { 
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          //dz = Complex(invGD, 0.0);
          
          dF[iq][iv] = fac * dw_over_dz * invGD; //dz; 
        }
        else  // copy for repeated occurences
        {
          dF[iq][iv] = dF[first_frequency][iv]; 
        }
      }
      else if(derivatives_data(iq) == JQT_temperature)
      {
        dz = Complex(-dL0_dT, dG0_dT) - z * dGD_dT;
        
        dF[iq][iv] = -F[iv] * dGD_dT;
        dF[iq][iv] += fac * dw_over_dz * dz;
        dF[iq][iv] *= invGD;
      }
      else if(derivatives_data(iq) == JQT_line_center) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[iq][iv] = -F[iv] * GD_div_F0;
          dF[iq][iv] += dw_over_dz * dz;
          dF[iq][iv] *= fac * invGD;
        }
      }
      else if(derivatives_data(iq) == JQT_line_gamma_self or 
        derivatives_data(iq) == JQT_line_gamma_selfexponent or
        derivatives_data(iq) == JQT_line_pressureshift_self or
        derivatives_data(iq) == JQT_line_gamma_foreign or
        derivatives_data(iq) == JQT_line_gamma_foreignexponent or
        derivatives_data(iq) == JQT_line_pressureshift_foreign or
        derivatives_data(iq) == JQT_line_gamma_water or
        derivatives_data(iq) == JQT_line_gamma_waterexponent or 
        derivatives_data(iq) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
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
      else if(derivatives_data(iq) == JQT_magnetic_magntitude)// No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        // dz = Complex(- zeeman_df * invGD, 0.0);
        
        dF[iq][iv] = fac * dw_over_dz * (- zeeman_df * invGD); //* dz; 
      }
      else if(derivatives_data(iq) == JQT_line_mixing_DF or
        derivatives_data(iq) == JQT_line_mixing_DF0 or
        derivatives_data(iq) == JQT_line_mixing_DF1 or
        derivatives_data(iq) == JQT_line_mixing_DFexp)
      {
        // dz = Complex(-invGD, 0.0);
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][iv] = fac * dw_over_dz * (-invGD); //* dz;
        }
      }
    }
  }
}


void Linefunctions::set_doppler(ComplexVector& F, // Sets the full complex line shape without line mixing
                                ArrayOfComplexVector& dF,
                                const Vector& f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& GD_div_F0,
                                const PropmatPartialsData& derivatives_data,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dGD_div_F0_dT)
{
  set_faddeeva_algorithm916(F,
                            dF, 
                            f_grid, 
                            zeeman_df, 
                            magnetic_magnitude, 
                            F0_noshift, 
                            GD_div_F0, 0.0, 0.0,
                            derivatives_data,
                            quantum_identity,
                            dGD_div_F0_dT);
}

void Linefunctions::set_faddeeva_from_full_linemixing(ComplexVector& F, 
                                                      ArrayOfComplexVector& dF,
                                                      const Vector& f_grid,
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
  F = fac;
  
  // Ratio of the Lorentz halfwidth to the Doppler halfwidth
  const Complex z0 = -eigenvalue * invGD;
  
  const Index first_pressure_broadening = derivatives_data.get_first_pressure_term(),
  first_frequency = derivatives_data.get_first_frequency();
  
  // frequency in units of Doppler
  for (Index iv=0; iv<nf; iv++)
  {
    dx = f_grid[iv] * invGD;
    z = z0 + dx;
    w = Faddeeva::w(z);
    
    F[iv] *= w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq==0)
        dw_over_dz = 2.0 * (z * w - sqrtInvPI);
      
      if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w) // No //external inputs
      { 
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          //dz = Complex(invGD, 0.0);
          
          dF[iq][iv] = fac * dw_over_dz * invGD; //dz; 
        }
        else  // copy for repeated occurences
        {
          dF[iq][iv] = dF[first_frequency][iv]; 
        }
      }
      else if(derivatives_data(iq) == JQT_temperature)
      {
        dz = (deigenvalue_dT - dL0_dT) - z * dGD_dT;
        
        dF[iq][iv] = -F[iv] * dGD_dT;
        dF[iq][iv] += fac * dw_over_dz * dz;
        dF[iq][iv] *= invGD;
      }
      else if(derivatives_data(iq) == JQT_line_center) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[iq][iv] = -F[iv] * GD_div_F0;
          dF[iq][iv] += dw_over_dz * dz;
          dF[iq][iv] *= fac * invGD;
        }
      }
      else if(derivatives_data(iq) == JQT_line_gamma_self or 
        derivatives_data(iq) == JQT_line_gamma_selfexponent or
        derivatives_data(iq) == JQT_line_pressureshift_self or
        derivatives_data(iq) == JQT_line_gamma_foreign or
        derivatives_data(iq) == JQT_line_gamma_foreignexponent or
        derivatives_data(iq) == JQT_line_pressureshift_foreign or
        derivatives_data(iq) == JQT_line_gamma_water or
        derivatives_data(iq) == JQT_line_gamma_waterexponent or 
        derivatives_data(iq) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
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
      else if(derivatives_data(iq) == JQT_line_mixing_DF or
        derivatives_data(iq) == JQT_line_mixing_DF0 or
        derivatives_data(iq) == JQT_line_mixing_DF1 or
        derivatives_data(iq) == JQT_line_mixing_DFexp)
      {
        // dz = Complex(-invGD, 0.0);
        
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][iv] = fac * dw_over_dz * (-invGD); //* dz;
      }
    }
    
  }
}


void Linefunctions::set_hui_etal_1978(ComplexVector& F,
                                      ArrayOfComplexVector& dF,
                                      const Vector& f_grid,
                                      const Numeric& zeeman_df,
                                      const Numeric& magnetic_magnitude,
                                      const Numeric& F0_noshift,
                                      const Numeric& GD_div_F0,
                                      const Numeric& G0,
                                      const Numeric& L0,
                                      const PropmatPartialsData& derivatives_data,
                                      const QuantumIdentifier& quantum_identity,
                                      const Numeric& dGD_div_F0_dT,
                                      const Numeric& dG0_dT,
                                      const Numeric& dL0_dT)
{
  /*
   * Exact copy of get_faddeeva_algorithm916 but with the change of A and B calculations from empirical numbers
   * 
   * For future extensions and corrections, make sure get_faddeeva_algorithm916 first works and then update this
   */
  
  // For calculations
  Numeric dx;
  Complex w, z, dw_over_dz, dz, A, B;
  
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // Doppler broadening and line center
  const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + L0;
  const Numeric GD = GD_div_F0 * F0;
  const Numeric invGD = 1.0 / GD;
  const Numeric dGD_dT = dGD_div_F0_dT * F0;
  
  // constant normalization factor for voigt
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
    //  Since this is ported from old FORTRAN code, intention must be that A and B coeffs are floats?
    A = (((((.5641896*z+5.912626)*z+30.18014)*z+
    93.15558)*z+181.9285)*z+214.3824)*z+122.6079;
    B = ((((((z+10.47986)*z+53.99291)*z+170.3540)*z+
    348.7039)*z+457.3345)*z+352.7306)*z+122.6079;
    
    F[iv] = fac * w;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(iq==0)
        dw_over_dz = 2.0 * (z * w - sqrtInvPI);
      
      if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w) // No //external inputs
      { 
        // If this is the first time it is calculated this frequency bin, do the full calculation
        if(first_frequency == iq)
        {
          //dz = Complex(invGD, 0.0);
          
          dF[iq][iv] = fac * dw_over_dz * invGD; //dz; 
        }
        else  // copy for repeated occurences
        {
          dF[iq][iv] = dF[first_frequency][iv]; 
        }
      }
      else if(derivatives_data(iq) == JQT_temperature)
      {
        dz = Complex(-dL0_dT, dG0_dT) - z * dGD_dT;
        
        dF[iq][iv] = -F[iv] * dGD_dT;
        dF[iq][iv] += fac * dw_over_dz * dz;
        dF[iq][iv] *= invGD;
      }
      else if(derivatives_data(iq) == JQT_line_center) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[iq][iv] = -F[iv] * GD_div_F0;
          dF[iq][iv] += dw_over_dz * dz;
          dF[iq][iv] *= fac * invGD;
        }
      }
      else if(derivatives_data(iq) == JQT_line_gamma_self or 
        derivatives_data(iq) == JQT_line_gamma_selfexponent or
        derivatives_data(iq) == JQT_line_pressureshift_self or
        derivatives_data(iq) == JQT_line_gamma_foreign or
        derivatives_data(iq) == JQT_line_gamma_foreignexponent or
        derivatives_data(iq) == JQT_line_pressureshift_foreign or
        derivatives_data(iq) == JQT_line_gamma_water or
        derivatives_data(iq) == JQT_line_gamma_waterexponent or 
        derivatives_data(iq) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
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
      else if(derivatives_data(iq) == JQT_magnetic_magntitude)// No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        // dz = Complex(- zeeman_df * invGD, 0.0);
        
        dF[iq][iv] = fac * dw_over_dz * (- zeeman_df * invGD); //* dz; 
      }
      else if(derivatives_data(iq) == JQT_line_mixing_DF or
        derivatives_data(iq) == JQT_line_mixing_DF0 or
        derivatives_data(iq) == JQT_line_mixing_DF1 or
        derivatives_data(iq) == JQT_line_mixing_DFexp)
      {
        // dz = Complex(-invGD, 0.0);
        
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][iv] = fac * dw_over_dz * (-invGD); //* dz;
      }
    } 
  }
}


void Linefunctions::set_o2_non_resonant(ComplexVector& F,
                                        ArrayOfComplexVector& dF,
                                        const Vector& f_grid,
                                        const Numeric& F0,
                                        const Numeric& G0,
                                        const PropmatPartialsData& derivatives_data,
                                        const QuantumIdentifier& quantum_identity,
                                        const Numeric& dG0_dT)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  const Numeric fac = G0 * invPI, invG0 = 1.0/G0;
  
  for(Index iv = 0; iv < nf; iv++)
  {
    F[iv] =  fac / ((f_grid[iv]-F0) * (f_grid[iv]-F0));
  }
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JQT_temperature)
    {
      dF[iq] = F;
      dF[iq] *= dG0_dT * invG0;
    }
    else if(derivatives_data(iq) == JQT_frequency or
      derivatives_data(iq) == JQT_wind_magnitude or
      derivatives_data(iq) == JQT_wind_u or
      derivatives_data(iq) == JQT_wind_v or
      derivatives_data(iq) == JQT_wind_w)
    {
      dF[iq] = F;
      for(Index iv = 0; iv < nf; iv++)
      {
        dF[iq][iv] /= -2.0 / (f_grid[iv] - F0);
      }
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq] = F;
        for(Index iv = 0; iv < nf; iv++)
        {
          dF[iq][iv] *= 2.0 / (f_grid[iv] - F0);
        }
      }
    }
    else if(derivatives_data(iq) == JQT_line_gamma_self or 
      derivatives_data(iq) == JQT_line_gamma_selfexponent or
      derivatives_data(iq) == JQT_line_gamma_foreign or
      derivatives_data(iq) == JQT_line_gamma_foreignexponent or
      derivatives_data(iq) == JQT_line_gamma_water or
      derivatives_data(iq) == JQT_line_gamma_waterexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq] = F;
        dF[iq] *= invG0;
      }
    }
  }
}


void Linefunctions::apply_linemixing(ComplexVector& F,
                                     ArrayOfComplexVector& dF,
                                     const Numeric& Y,
                                     const Numeric& G,
                                     const PropmatPartialsData& derivatives_data,
                                     const QuantumIdentifier& quantum_identity,
                                     const Numeric& dY_dT,
                                     const Numeric& dG_dT)
{
  const Index nf = F.nelem(), nppd = derivatives_data.nelem();
  
  const Complex LM = Complex(1.0 + G, Y);
  const Complex dLM_dT = Complex(dG_dT, dY_dT);
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JQT_temperature)
    {
      dF[iq] *= LM;
      for(Index iv = 0; iv < nf; iv++)
      {
        dF[iq][iv] += F[iv] * dLM_dT;
      }
    }
    else if(derivatives_data(iq) == JQT_line_mixing_Y or
      derivatives_data(iq) == JQT_line_mixing_Y0 or 
      derivatives_data(iq) == JQT_line_mixing_Y1 or 
      derivatives_data(iq) == JQT_line_mixing_Yexp)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq] = F;
        dF[iq] *= Complex(0.0, 1.0);
      }
    }
    else if(derivatives_data(iq) == JQT_line_mixing_G or
      derivatives_data(iq) == JQT_line_mixing_G0 or 
      derivatives_data(iq) == JQT_line_mixing_G1 or 
      derivatives_data(iq) == JQT_line_mixing_Gexp)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] = F;
    }
    else
    {
      dF[iq] *= LM;
    }
  }
  
  F *= LM;
}


void Linefunctions::apply_rosenkranz_quadratic(ComplexVector& F,
                                               ArrayOfComplexVector& dF,
                                               const Vector& f_grid,
                                               const Numeric& F0,
                                               const Numeric& T,
                                               const PropmatPartialsData& derivatives_data,
                                               const QuantumIdentifier& quantum_identity)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  const Numeric invF0 = 1.0/F0;
  const Numeric mafac = (PLANCK_CONST) / (2 * BOLTZMAN_CONST * T) /
  sinh((PLANCK_CONST * F0) / (2 * BOLTZMAN_CONST * T)) * invF0;
  
  Numeric dmafac_dF0_div_fun = 0, dmafac_dT_div_fun = 0;
  if(derivatives_data.do_line_center())
  {
    dmafac_dF0_div_fun = -invF0 - 
    PLANCK_CONST/(2*BOLTZMAN_CONST*T*tanh(F0*PLANCK_CONST/(2*BOLTZMAN_CONST*T)));
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
      dF[iq][iv] *= fun;
      if(derivatives_data(iq) == JQT_temperature)
      {
        dF[iq][iv] += dmafac_dT_div_fun * F[iv];
      }
      else if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][iv] += dmafac_dF0_div_fun * F[iv];
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
      {
        dF[iq][iv] += (2.0 / f_grid[iv]) * F[iv];
      }
    }
  }
}


void Linefunctions::apply_VVH(ComplexVector& F,
                              ArrayOfComplexVector& dF,
                              const Vector& f_grid,
                              const Numeric& F0,
                              const Numeric& T,
                              const PropmatPartialsData& derivatives_data,
                              const QuantumIdentifier& quantum_identity)
{ 
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // 2kT is constant for the loop
  const Numeric kT = 2.0 * BOLTZMAN_CONST * T;
  
  // denominator is constant for the loop
  const Numeric absF0 = abs(F0);
  const Numeric tanh_f0part = tanh(PLANCK_CONST * absF0 / kT);
  const Numeric denom = absF0 * tanh_f0part;
  
  Numeric fac_df0 = 0; 
  if(derivatives_data.do_line_center())
    fac_df0 = (-1.0/absF0 + PLANCK_CONST*tanh_f0part/(kT) - PLANCK_CONST/(kT*tanh_f0part)) * F0/absF0;
  
  for(Index iv=0; iv < nf; iv++)
  {
    const Numeric tanh_fpart = tanh( PLANCK_CONST * f_grid[iv] / kT );
    const Numeric fun = f_grid[iv] * tanh_fpart / denom;
    F[iv] *= fun;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      dF[iq][iv] *= fun;
      if(derivatives_data(iq) == JQT_temperature)
      {
        dF[iq][iv] += (-PLANCK_CONST*(denom - absF0/tanh_f0part - 
        f_grid[iv]*tanh_fpart + f_grid[iv]/tanh_fpart)/(kT*T)) * F[iv];
      }
      else if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][iv] += fac_df0 * F[iv];
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
      {
        dF[iq][iv] += (1.0/f_grid[iv] -PLANCK_CONST*tanh_fpart/kT + PLANCK_CONST/(kT*tanh_fpart)) * F[iv];
      }
    }
  }
}


void Linefunctions::apply_VVW(ComplexVector& F,
                              ArrayOfComplexVector& dF,
                              const Vector& f_grid,
                              const Numeric& F0,
                              const PropmatPartialsData& derivatives_data,
                              const QuantumIdentifier& quantum_identity)
{
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // denominator is constant for the loop
  const Numeric invF0 = 1.0 / F0;
  const Numeric invF02 = invF0 * invF0;
  
  for(Index iv = 0; iv < nf; iv++)
  {
    const Numeric fun = f_grid[iv] * f_grid[iv] * invF02;
    F[iv] *= fun;
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      dF[iq][iv] *= fun;
      if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][iv] -= 2.0 * invF0 * F[iv];
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
      {
        dF[iq][iv] += 2.0 * f_grid[iv] * F[iv];
      }
    }
  }
}


void Linefunctions::apply_linestrength(ComplexVector& F,
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
                                       const Numeric& dK2_dF0)
{
  const Index nppd = derivatives_data.nelem();
  
  const Numeric invQT = 1.0/QT;
  const Numeric QT_ratio = QT0 * invQT;
  const Numeric S = S0 * isotopic_ratio * QT_ratio * K1 * K2;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JQT_temperature)
    {
      Eigen::VectorXcd eig_dF = MapToEigen(dF[iq]);
      
      eig_dF *= S;
      eig_dF += MapToEigen(F) * (S0 * isotopic_ratio * QT_ratio * 
      (K1 * dK2_dT + 
      dK1_dT * K2 - 
      invQT * dQT_dT * K1 * K2));
    }
    else if(derivatives_data(iq) == JQT_line_strength)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq] = F;
        dF[iq] *= isotopic_ratio;
      }
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        Eigen::VectorXcd eig_dF = MapToEigen(dF[iq]);
        
        eig_dF *= S;
        eig_dF += MapToEigen(F) * (S0 * isotopic_ratio * QT_ratio * K1 * dK2_dF0);
      }
    }
    else
    {
      dF[iq] *= S;
    }
  }
  F *= S;
}


/*! Applies the line strength to the line shape using the complex line strenth of line mixing
 * 
 * \param F                        Lineshape vector (already set)
 * \param dF                       Lineshape derivative vector (already set)
 * \param F0                       Line center frequency[Hz]
 * \param T                        Atmospheric temperature at level
 * \param S_LM                     Complex linestrength
 * \param isotopic_ratio           Ratio of isotopologue in atmosphere at level
 * \param derivatives_data         Information on the partial derivative
 * \param dS_LM_dT                 Derivative of S_LM with respect to T
 * 
 */
void Linefunctions::apply_linestrength_from_full_linemixing(ComplexVector& F,
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
    if(derivatives_data(iq) == JQT_temperature)
    {
      Eigen::VectorXcd eig_dF = MapToEigen(dF[iq]);
      
      eig_dF *= S_LM;
      eig_dF += MapToEigen(F) * (dS_LM_dT * f0_factor + 
      S_LM * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio;
    }
    else if(derivatives_data(iq) == JQT_line_strength)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else if(derivatives_data(iq) == JQT_line_center)
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
 * \param F                        Lineshape vector (already set)
 * \param dF                       Lineshape derivative vector (already set)
 * \param F0                       Line center frequency[Hz]
 * \param T                        Atmospheric temperature at level
 * \param rho                      Density (of molecules at level)
 * \param isotopic_ratio           Ratio of isotopologue in atmosphere at level
 * \param derivatives_data         Information on the partial derivative
 * \param drho_dT                  Derivative of rho with respect to T
 * 
 */
void Linefunctions::apply_dipole(ComplexVector& F,
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
    if(derivatives_data(iq) == JQT_temperature)
    {
      ComplexMatrixViewMap eig_dF = MapToEigen(dF[iq]);
      
      eig_dF *= S * f0_factor;
      eig_dF += MapToEigen(F) * (d0 * d0 * (drho_dT * f0_factor + 
      rho * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio);
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        throw std::runtime_error("Not working yet");
    }
    else if(derivatives_data(iq) == JQT_line_strength)
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
 * \param dF                       array of already applied lineshape jacobians that need 
 *                                 to have the derivatives data applied
 * \param derivatives_data         information on where the derivatives need to be applied
 * \param dgamma                   vector the length of the number of pressure broadening
 *                                 parameters that need to be applied sorted by their
 *                                 apprearance in derivatives_data
 * 
 */
void Linefunctions::apply_pressurebroadening_jacobian(ArrayOfComplexVector& dF,
                                                      const PropmatPartialsData& derivatives_data,
                                                      const QuantumIdentifier& quantum_identity,
                                                      const ComplexVector& dgamma)
{
  const Index nppd = derivatives_data.nelem(), ng = dgamma.nelem();
  
  Index ipd = 0;
  
  // Length of dgamma must be the same as total number of instances of pressure broadening jacobians
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(ipd == ng) break;
    
    if(derivatives_data(iq) == JQT_line_gamma_self)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_foreign)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_water)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_pressureshift_self)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_pressureshift_foreign)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_pressureshift_water)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= 0.0;
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_selfexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_foreignexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_waterexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq] *= dgamma[ipd];
      else
        continue;
    }
    else
      continue;
    
    // Only activate this when something hit the target
    ++ipd;
  }
}


Numeric Linefunctions::DopplerConstant(const Numeric T, const Numeric mass)
{
  return doppler_const * sqrt(T / mass);
}


Numeric Linefunctions::dDopplerConstant_dT(const Numeric T, const Numeric mass)
{
  return doppler_const * 0.5 * sqrt(1.0 / mass / T);
}


void Linefunctions::apply_nonlte(ComplexVector& F, 
                                 ArrayOfComplexVector& dF, 
                                 ComplexVector& N, 
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
                                 const Numeric& dK4_dTu)
{
  const Index nppd = derivatives_data.nelem(), nf = F.nelem();
  
  const Numeric scaled_ratio = K4/K3 - 1.0;
  
  // Set the non-lte source factors
  N = F;
  N *= scaled_ratio;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    dN[iq] = dF[iq];
    dN[iq] *= scaled_ratio;
    dF[iq] *= K3;
    
    if(derivatives_data(iq) == JQT_temperature)
    {
      const Numeric dscaled_ratio_dT = (dK4_dT - dK3_dT / K3) / K3;
      
      for(Index iv = 0; iv < nf; iv++)
      {
        dF[iq][iv] += F[iv] * dK3_dT;
        dN[iq][iv] += F[iv] * dscaled_ratio_dT;
      }
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        const Numeric dscaled_ratio_dF0 = - dK3_dF0 / K3 / K3;
        
        for(Index iv = 0; iv < nf; iv++)
        {
          dF[iq][iv] += F[iv] * dK3_dF0;
          dN[iq][iv] += F[iv] * dscaled_ratio_dF0;
        }
      }
    }
    else if(derivatives_data(iq) == JQT_nlte_temperature)
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
          dF[iq][iv] += F[iv] * dK3_dTl;
          dN[iq][iv] += F[iv] * dscaled_ratio_dTl;
        }
      }
      
      if(quantum_identity.QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX] >
        derivatives_data.jac(iq).QuantumIdentity().QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL] or
        derivatives_data.jac(iq).QuantumIdentity().Type() == QuantumIdentifier::ALL)
      {
        const Numeric dscaled_ratio_dTu = (dK4_dTu - dK3_dTu / K3) / K3;
        
        for(Index iv = 0; iv < nf; iv++)
        {
          dF[iq][iv] += F[iv] * dK3_dTu;
          dN[iq][iv] += F[iv] * dscaled_ratio_dTu;
        }
      }
    }
  }
  
  // Finish by scaling F to the true value
  F *= K3;
}


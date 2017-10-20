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
 * \param L0 Speed-independent pressure shift term plus second order line mixing shift
 * \param derivatives_data Information about the derivatives in dF
 * \param quantum_identity ID of the absorption line
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::set_lorentz(ComplexVectorView F, // Sets the full complex line shape without line mixing
                                ArrayOfComplexVector& dF,
                                ConstVectorView f_grid,
                                const Numeric& zeeman_df,
                                const Numeric& magnetic_magnitude,
                                const Numeric& F0_noshift,
                                const Numeric& G0,
                                const Numeric& L0, 
                                const PropmatPartialsData& derivatives_data,
                                const QuantumIdentifier& quantum_identity,
                                const Numeric& dG0_dT,
                                const Numeric& dL0_dT,
                                const ComplexRange& df_range)
{ 
  const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
  
  // Compute the true central frequency
  const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude;
  
  // Sigma change of F and F0?
  const Complex denom0 = Complex(G0, F0);
  
  for(Index iv = 0; iv < nf; iv++)
  {
    F[iv] = invPI / (denom0 - Complex(0.0, f_grid[iv]));
  }
  
  // If there are partial derivatives, then the last 
  // one will now contain a temporary speed increasing 
  // variable since d (1/F(x)) / dx = - (dF(x)/dx) / F(x)^2
  // hopefully, C++ knows that A[i] = A[i] does nothing....
  if(nppd > 0)
  {
    dF[nppd-1][df_range] = F;
    dF[nppd-1][df_range] *= F;
    dF[nppd-1][df_range] *= PI;
  }
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    
    if(derivatives_data(iq) == JQT_temperature)
    {
      // Temperature derivative only depends on how pressure shift and broadening change
      dF[iq][df_range] = dF[nppd-1][df_range];
      dF[iq][df_range] *= Complex(dG0_dT, dL0_dT);
    }
    else if(derivatives_data(iq) == JQT_frequency or
      derivatives_data(iq) == JQT_wind_magnitude or
      derivatives_data(iq) == JQT_wind_u or
      derivatives_data(iq) == JQT_wind_v or
      derivatives_data(iq) == JQT_wind_w)
    {
      // Frequency scale 1 to -1 linearly
      dF[iq][df_range] = dF[nppd-1][df_range];
      dF[iq][df_range] *= Complex(0.0, -1.0);
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        // Line center scales 1 to 1 linearly
        dF[iq][df_range] = dF[nppd-1][df_range];
        dF[iq][df_range] *= Complex(0.0, 1.0);
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
      derivatives_data(iq) == JQT_line_pressureshift_water) 
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        // Pressure broadening will be dealt with in another function, though the partial derivative
        // messes some things up (should this 1+1i be removed?)
        dF[iq][df_range] = dF[nppd-1][df_range];
        dF[iq][df_range] *= Complex(1.0, 1.0);
      }
    }
    else if(derivatives_data(iq) == JQT_magnetic_magntitude)
    {
      // Magnetic magnitude changes like line center in part
      // FIXME: Add magnetic components here
      dF[iq][df_range] = dF[nppd-1][df_range];
      dF[iq][df_range] *= Complex(0.0, zeeman_df);
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
                            const ComplexRange& df_range)
{
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
      invG = 1.0 / (1.0 - (FVC - eta * C0_m1p5_C2)*A + eta * B); // WARNING What to do at denom 0?  Catch earlier?
    }
    else
    {
      invG = 1.0 / (1.0 - FVC*A);  // WARNING What to do at denom 0?  Catch earlier?
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
          
          dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dg); 
        }
        else  // copy for repeated occurences
        {
          dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
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
        
        dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dg); 
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
          
          dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dg); 
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
            
            dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dg); 
          }
          else  // copy for repeated occurences
          {
            dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
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
        
        dF[iq][df_range][iv] = invG * (invPI * dA - F[iv] * dg); 
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
 * \param L0 Speed-independent pressure shift term plus second order line mixing shift
 * \param derivatives_data Information about the derivatives in dF
 * \param quantum_identity ID of the absorption line
 * \param dGD_div_F0_dT Temperature derivative of GD_div_F0
 * \param dG0_dT Temperature derivative of G0
 * \param dL0_dT Temperature derivative of L0
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
                                              const PropmatPartialsData& derivatives_data,
                                              const QuantumIdentifier& quantum_identity,
                                              const Numeric& dGD_div_F0_dT,
                                              const Numeric& dG0_dT,
                                              const Numeric& dL0_dT,
                                              const ComplexRange& df_range)
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
          
          dF[iq][df_range][iv] = fac * dw_over_dz * invGD; //dz; 
        }
        else  // copy for repeated occurences
        {
          dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
        }
      }
      else if(derivatives_data(iq) == JQT_temperature)
      {
        dz = Complex(-dL0_dT, dG0_dT) - z * dGD_dT;
        
        dF[iq][df_range][iv] = -F[iv] * dGD_dT;
        dF[iq][df_range][iv] += fac * dw_over_dz * dz;
        dF[iq][df_range][iv] *= invGD;
      }
      else if(derivatives_data(iq) == JQT_line_center) // No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[iq][df_range][iv] = -F[iv] * GD_div_F0;
          dF[iq][df_range][iv] += dw_over_dz * dz;
          dF[iq][df_range][iv] *= fac * invGD;
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
            dF[iq][df_range][iv] = fac * dw_over_dz * dz; 
          }
          else
          {
            dF[iq][df_range][iv] = dF[first_frequency][df_range][iv]; 
          }
        }
      }
      else if(derivatives_data(iq) == JQT_magnetic_magntitude)// No //external inputs --- errors because of frequency shift when Zeeman is used?
      {
        // dz = Complex(- zeeman_df * invGD, 0.0);
        
        dF[iq][df_range][iv] = fac * dw_over_dz * (- zeeman_df * invGD); //* dz; 
      }
      else if(derivatives_data(iq) == JQT_line_mixing_DF or
        derivatives_data(iq) == JQT_line_mixing_DF0 or
        derivatives_data(iq) == JQT_line_mixing_DF1 or
        derivatives_data(iq) == JQT_line_mixing_DFexp)
      {
        // dz = Complex(-invGD, 0.0);
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        {
          dF[iq][df_range][iv] = fac * dw_over_dz * (-invGD); //* dz;
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
                                const ComplexRange& df_range)
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
                            dGD_div_F0_dT, 0.0, 0.0,
                            df_range);
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
 * \param L0 Speed-independent pressure shift term plus second order line mixing shift
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


/*!
 * Applies line mixing to already set lineshape
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
void Linefunctions::apply_linemixing(ComplexVectorView F,
                                     ArrayOfComplexVector& dF,
                                     const Numeric& Y,
                                     const Numeric& G,
                                     const PropmatPartialsData& derivatives_data,
                                     const QuantumIdentifier& quantum_identity,
                                     const Numeric& dY_dT,
                                     const Numeric& dG_dT,
                                     const ComplexRange& df_range)
{
  const Index nf = F.nelem(), nppd = derivatives_data.nelem();
  
  const Complex LM = Complex(1.0 + G, -Y);
  const Complex dLM_dT = Complex(dG_dT, -dY_dT);
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JQT_temperature)
    {
      dF[iq][df_range] *= LM;
      for(Index iv = 0; iv < nf; iv++)
      {
        dF[iq][df_range][iv] += F[iv] * dLM_dT;
      }
    }
    else if(derivatives_data(iq) == JQT_line_mixing_Y or
      derivatives_data(iq) == JQT_line_mixing_Y0 or 
      derivatives_data(iq) == JQT_line_mixing_Y1 or 
      derivatives_data(iq) == JQT_line_mixing_Yexp)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        dF[iq][df_range] = F;
        dF[iq][df_range] *= Complex(0.0, -1.0);
      }
    }
    else if(derivatives_data(iq) == JQT_line_mixing_G or
      derivatives_data(iq) == JQT_line_mixing_G0 or 
      derivatives_data(iq) == JQT_line_mixing_G1 or 
      derivatives_data(iq) == JQT_line_mixing_Gexp)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] = F;
    }
    else
    {
      dF[iq][df_range] *= LM;
    }
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
void Linefunctions::apply_rosenkranz_quadratic(ComplexVectorView F,
                                               ArrayOfComplexVector& dF,
                                               ConstVectorView f_grid,
                                               const Numeric& F0,
                                               const Numeric& T,
                                               const PropmatPartialsData& derivatives_data,
                                               const QuantumIdentifier& quantum_identity,
                                               const ComplexRange& df_range)
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
      dF[iq][df_range][iv] *= fun;
      if(derivatives_data(iq) == JQT_temperature)
      {
        dF[iq][df_range][iv] += dmafac_dT_div_fun * F[iv];
      }
      else if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][df_range][iv] += dmafac_dF0_div_fun * F[iv];
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
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
void Linefunctions::apply_VVH(ComplexVectorView F,
                              ArrayOfComplexVector& dF,
                              ConstVectorView f_grid,
                              const Numeric& F0,
                              const Numeric& T,
                              const PropmatPartialsData& derivatives_data,
                              const QuantumIdentifier& quantum_identity,
                              const ComplexRange& df_range)
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
      if(derivatives_data(iq) == JQT_temperature)
      {
        dF[iq][df_range][iv] += (-PLANCK_CONST*(denom - F0/tanh_f0part - 
        f_grid[iv]*tanh_fpart + f_grid[iv]/tanh_fpart)/(kT*T)) * F[iv];
      }
      else if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][df_range][iv] += fac_df0 * F[iv];
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
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
void Linefunctions::apply_VVW(ComplexVectorView F,
                              ArrayOfComplexVector& dF,
                              ConstVectorView f_grid,
                              const Numeric& F0,
                              const PropmatPartialsData& derivatives_data,
                              const QuantumIdentifier& quantum_identity,
                              const ComplexRange& df_range)
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
      dF[iq][df_range][iv] *= fun;
      if(derivatives_data(iq) == JQT_line_center)
      {
        if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
          dF[iq][df_range][iv] -= 2.0 * invF0 * F[iv];
      }
      else if(derivatives_data(iq) == JQT_frequency or
        derivatives_data(iq) == JQT_wind_magnitude or
        derivatives_data(iq) == JQT_wind_u or
        derivatives_data(iq) == JQT_wind_v or
        derivatives_data(iq) == JQT_wind_w)
      {
        dF[iq][df_range][iv] += 2.0 * f_grid[iv] * F[iv];
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
void Linefunctions::apply_linestrength(ComplexVectorView F,
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
                                       const ComplexRange& df_range)
{
  const Index nppd = derivatives_data.nelem();
  
  const Numeric invQT = 1.0/QT;
  const Numeric QT_ratio = QT0 * invQT;
  const Numeric S = S0 * isotopic_ratio * QT_ratio * K1 * K2;
  
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(derivatives_data(iq) == JQT_temperature)
    {
      Eigen::VectorXcd eig_dF = MapToEigen(dF[iq][df_range]);
      
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
        dF[iq][df_range] = F;
        dF[iq][df_range] *= isotopic_ratio;
      }
    }
    else if(derivatives_data(iq) == JQT_line_center)
    {
      
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
      {
        Eigen::VectorXcd eig_dF = MapToEigen(dF[iq][df_range]);
        
        eig_dF *= S;
        eig_dF += MapToEigen(F) * (S0 * isotopic_ratio * QT_ratio * K1 * dK2_dF0);
      }
    }
    else
    {
      dF[iq][df_range] *= S;
    }
  }
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
 * \retval dF Lineshape derivative
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param quantum_identity ID of the absorption line
 * \param dgamma Derivatives in order as they appear that are related to pressure broadening coefficients
 * \param df_range Frequency range to use inside dF
 * 
 */
void Linefunctions::apply_pressurebroadening_jacobian(ArrayOfComplexVector& dF,
                                                      const PropmatPartialsData& derivatives_data,
                                                      const QuantumIdentifier& quantum_identity,
                                                      const ComplexVector& dgamma,
                                                      const ComplexRange& df_range)
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
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_foreign)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_water)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_pressureshift_self)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_pressureshift_foreign)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_pressureshift_water)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_selfexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_foreignexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else if(derivatives_data(iq) == JQT_line_gamma_waterexponent)
    {
      if(quantum_identity > derivatives_data.jac(iq).QuantumIdentity())
        dF[iq][df_range] *= dgamma[ipd];
      else
        continue;
    }
    else
      continue;
    
    // Only activate this when something hit the target
    ++ipd;
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
 * \retval F Lineshape
 * \retval dF Lineshape derivative
 * \retval N Non-lte lineshape
 * \retval dN Non-lte lineshape derivative
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
void Linefunctions::apply_nonlte(ComplexVectorView F, 
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
                                 const ComplexRange& df_range)
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
    
    if(derivatives_data(iq) == JQT_temperature)
    {
      const Numeric dscaled_ratio_dT = (dK4_dT - dK3_dT / K3) / K3;
      
      for(Index iv = 0; iv < nf; iv++)
      {
        dF[iq][df_range][iv] += F[iv] * dK3_dT;
        dN[iq][df_range][iv] += F[iv] * dscaled_ratio_dT;
      }
    }
    else if(derivatives_data(iq) == JQT_line_center)
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
 * 
 * \param derivatives_data Information about the derivatives in dF
 * \param line line-record containing most line parameters
 * \param volume_mixing_ratio_of_all_species As name suggests
 * \param nlte_temperatures As name suggests
 * \param pressure As name suggests
 * \param temperature As name suggests
 * \param doppler_constant Frequency-independent part of the Doppler broadening
 * \param partial_pressure Pressure of species that line belongs to at this level
 * \param isotopologue_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param ddoppler_constant_dT Temperature derivative of doppler_constant
 * \param pressure_limit_for_linemixing As WSV lm_p_lim
 * \param zeeman_frequency_shift_constant Zeeman shift parameter for the line
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
                                                      ComplexRange& this_xsec_range,
                                                      const PropmatPartialsData& derivatives_data, 
                                                      const LineRecord& line, 
                                                      ConstVectorView f_grid, 
                                                      ConstVectorView volume_mixing_ratio_of_all_species, 
                                                      ConstVectorView nlte_temperatures, 
                                                      const Numeric& pressure, 
                                                      const Numeric& temperature, 
                                                      const Numeric& doppler_constant, 
                                                      const Numeric& partial_pressure, 
                                                      const Numeric& isotopologue_ratio,
                                                      const Numeric& magnetic_magnitude,
                                                      const Numeric& ddoppler_constant_dT,
                                                      const Numeric& pressure_limit_for_linemixing,
                                                      const Numeric& zeeman_frequency_shift_constant,
                                                      const Numeric& partition_function_at_temperature,
                                                      const Numeric& dpartition_function_at_temperature_dT,
                                                      const Numeric& partition_function_at_line_temperature,
                                                      const ArrayOfIndex& broad_spec_locations,
                                                      const Index& this_species_location_in_tags,
                                                      const Index& water_index_location_in_tags,
                                                      const Verbosity& verbosity,
                                                      const bool cutoff_call)
{
  // This call assumes the size of all variables are as 
  // they should be and makes no tests in this is the case
  
  // The equation being solved here:
  //   F = LM K1 K2 K3 QT/QT0 S0 f(...),
  //   N = LM K1 K2 (K4/K3 - 1) QT/QT0 S0 f(...),
  //   dF = d(F)/dx
  //   dN = d(N)/dx
  
  // Get the cutoff range if applicable
  const Index nf = f_grid.nelem();
  Range this_f_range(joker);
  const Numeric cutoff = line.CutOff();
  const bool need_cutoff = (cutoff > 0) and not cutoff_call;
  
  if(need_cutoff)
  {
    // Find range of simulations
    Index i_f_min = 0;
    Index i_f_max = nf-1;
    
    while(i_f_min < nf and (line.F() - cutoff) > f_grid[i_f_min])
      ++i_f_min;
    while(i_f_max >= i_f_min and (line.F() + cutoff) < f_grid[i_f_max])
      --i_f_max;
    ++i_f_max;
    
    const Index extent = i_f_max - i_f_min;
    
    this_f_range = Range(i_f_min, extent);
    this_xsec_range = ComplexRange(i_f_min, extent);
    
    if(extent < 1)
      return;
  }
  else
  {
    this_f_range = Range(joker);
    this_xsec_range = ComplexRange(joker);
  }
  
  // Calculation vectors
  ComplexVectorView F_calc = F[this_xsec_range];
  ConstVectorView f_grid_calc = f_grid[this_f_range];
  
  // Extract the quantum identify of the line
  const QuantumIdentifier QI = line.QuantumIdentity();
  
  // Line strength scaling that are line-dependent
  const Numeric gamma = stimulated_emission(temperature, line.F());
  const Numeric gamma_ref = stimulated_emission(line.Ti0(), line.F());
  const Numeric K1 = boltzman_ratio(temperature, line.Ti0(), line.Elow());
  const Numeric K2 = stimulated_relative_emission(gamma, gamma_ref);
  
  // Pressure broadening terms
  Numeric G0, G2, e, L0, L2, FVC;
  line.PressureBroadening().GetPressureBroadeningParams(G0, G2, e, L0, L2, FVC,
    line.Ti0()/temperature, pressure, partial_pressure, 
    this_species_location_in_tags, water_index_location_in_tags,
    broad_spec_locations, volume_mixing_ratio_of_all_species,
    verbosity);
  
  // Line mixing terms
  Numeric Y=0, G=0, DV=0;
  line.LineMixing().GetLineMixingParams(Y, G, DV, temperature, pressure, 
                                        pressure_limit_for_linemixing);
  
  // Treat line mixing as a shift on first order of pressure shift
  // Note that HTP-like line shapes might only use first order line 
  // mixing because of this assumption, so if anything changes on
  // this from this part has to be updated
  if(DV not_eq 0)
    L0 += DV;
  
  // Partial derivatives for temperature
  Numeric dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT, dY_dT, dG_dT, dDV_dT, dK1_dT, dK2_dT;
  if(derivatives_data.do_temperature())
  {
    // NOTE:  This will change if and when we have a good temperature understanding of HTP,
    // and uninitialized variables above will then be set...
    line.PressureBroadening().GetPressureBroadeningParams_dT(dG0_dT, dL0_dT,
      temperature, line.Ti0(), pressure, partial_pressure,
      this_species_location_in_tags, water_index_location_in_tags,
      broad_spec_locations, volume_mixing_ratio_of_all_species,
      verbosity);
    
    // Line mixing partial derivatives
    line.LineMixing().GetLineMixingParams_dT(dY_dT, dG_dT, dDV_dT, temperature, 
                                              derivatives_data.Temperature_Perturbation(),
                                              pressure, pressure_limit_for_linemixing);
   
    // See comment above L0 += DV;
    if(dDV_dT not_eq 0)
      dL0_dT += dDV_dT;
  
    // Line strength partial derivatives
    dK1_dT = dboltzman_ratio_dT(K1, temperature, line.Elow());
    dK2_dT = dstimulated_relative_emission_dT(gamma, gamma_ref, line.F());
  }
  
  // Partial derivatives due to pressure
  ComplexVector pressure_derivatives;
  if(derivatives_data.get_first_pressure_term() > -1)
    line.PressureBroadening().SetInternalDerivatives(pressure_derivatives, derivatives_data, QI, 
                                                     line.Ti0()/temperature, pressure, partial_pressure, 
                                                     this_species_location_in_tags, water_index_location_in_tags,
                                                     volume_mixing_ratio_of_all_species, verbosity);
  
  // Partial derivatives due to central frequency of the stimulated emission
  Numeric dK2_dF0;
  if(derivatives_data.do_line_center())
    dK2_dF0 = dstimulated_relative_emission_dF0(gamma, gamma_ref, temperature);
  
  // Set the line shape normalized to unity integration
  LineShapeType lst = LineShapeType::End;
  switch(line.GetLineShapeType())
  {
    case LineShapeType::ByPressureBroadeningData:
      switch(line.PressureBroadening().Type())
      {
        case PressureBroadeningData::PB_SD_AIR_VOLUME:
          // Above should be all methods of pressure broadening requiring HTP in ARTS to work by default
          lst = LineShapeType::HTP;
          set_htp(F_calc, dF, 
                  f_grid_calc, zeeman_frequency_shift_constant, magnetic_magnitude, 
                  line.F(), doppler_constant, 
                  G0, L0, G2, L2, e, FVC,
                  derivatives_data, QI,
                  ddoppler_constant_dT, 
                  dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT,
                  this_xsec_range);
          break;
        case PressureBroadeningData::PB_AIR_AND_WATER_BROADENING:
        case PressureBroadeningData::PB_PLANETARY_BROADENING:
        case PressureBroadeningData::PB_AIR_BROADENING:
          // Above should be all methods of pressure broadening requiring Voigt in ARTS by default
          lst = LineShapeType::Voigt;
          set_faddeeva_algorithm916(F_calc, dF, f_grid_calc, 
                                    zeeman_frequency_shift_constant, magnetic_magnitude, 
                                    line.F(), doppler_constant, 
                                    G0, L0, derivatives_data, QI,
                                    ddoppler_constant_dT, dG0_dT, dL0_dT,
                                    this_xsec_range);
          break;
        default:
          throw std::runtime_error("Developer has messed up and needs to add the key to the code above this error");
      }
      break;
    case LineShapeType::Doppler:
      lst = LineShapeType::Doppler;
      set_doppler(F_calc, dF, f_grid_calc, zeeman_frequency_shift_constant, magnetic_magnitude, 
                  line.F(), doppler_constant, derivatives_data, QI, ddoppler_constant_dT, this_xsec_range);
      break;
    case LineShapeType::HTP:
      set_htp(F_calc, dF, 
              f_grid_calc, zeeman_frequency_shift_constant, magnetic_magnitude, 
              line.F(), doppler_constant, 
              G0, L0, G2, L2, e, FVC,
              derivatives_data, QI,
              ddoppler_constant_dT, 
              dG0_dT, dL0_dT, dG2_dT, dL2_dT, de_dT, dFVC_dT,
              this_xsec_range);
      lst = LineShapeType::HTP;
      break;
    case LineShapeType::Lorentz:
      lst = LineShapeType::Lorentz;
      set_lorentz(F_calc, dF, f_grid_calc, zeeman_frequency_shift_constant, magnetic_magnitude, 
                  line.F(), G0, L0, derivatives_data, QI, dG0_dT, dL0_dT, this_xsec_range);
      break;
    case LineShapeType::Voigt:
      lst = LineShapeType::Voigt;
      set_faddeeva_algorithm916(F_calc, dF, f_grid_calc, 
                                zeeman_frequency_shift_constant, magnetic_magnitude, 
                                line.F(), doppler_constant, 
                                G0, L0, derivatives_data, QI,
                                ddoppler_constant_dT, dG0_dT, dL0_dT,
                                this_xsec_range);
      break;
    case LineShapeType::End:
    default:
      throw std::runtime_error("Cannot understand the requested line shape type.");
  }
  
  // Set the mirroring by repeating computations above using negative numbers for frequency of line related terms
  switch(line.GetMirroringType())
  {
    case MirroringType::None:
      break;
    case MirroringType::Lorentz:
    {
      ComplexVector Fm(F_calc.nelem());
      ArrayOfComplexVector dFm(dF.nelem());
      for(auto& aocv : dFm) aocv.resize(F_calc.nelem());
      set_lorentz(Fm, dFm, f_grid_calc, -zeeman_frequency_shift_constant, magnetic_magnitude, 
                  -line.F(), G0, -L0, derivatives_data, QI, dG0_dT, -dL0_dT);
      F_calc -= Fm;
      for(Index i = 0; i < dF.nelem(); i++) dF[i][this_xsec_range] -= dFm[i];
    }
    break;
    case MirroringType::SameAsLineShape:
    {
      ComplexVector Fm(F_calc.nelem());
      ArrayOfComplexVector dFm(dF.nelem());
      for(auto& aocv : dFm) aocv.resize(F_calc.nelem());
      switch(lst)
      {
        case LineShapeType::Doppler:
          set_doppler(Fm, dFm, f_grid_calc, -zeeman_frequency_shift_constant, magnetic_magnitude, 
                      -line.F(), -doppler_constant, derivatives_data, QI, -ddoppler_constant_dT);
          break;
        case LineShapeType::Lorentz:
          set_lorentz(Fm, dFm, f_grid_calc, -zeeman_frequency_shift_constant, magnetic_magnitude, 
                      -line.F(), G0, -L0, derivatives_data, QI, dG0_dT, -dL0_dT);
          break;
        case LineShapeType::Voigt:
          set_faddeeva_algorithm916(Fm, dFm, f_grid_calc, 
                                    -zeeman_frequency_shift_constant, magnetic_magnitude, 
                                    -line.F(), -doppler_constant, 
                                    G0, -L0, derivatives_data, QI,
                                    -ddoppler_constant_dT, dG0_dT, -dL0_dT);
          break;
        case LineShapeType::HTP:
          // WARNING: This mirroring is not tested and it might require, e.g., FVC to be treated differently
          set_htp(Fm, dFm, f_grid_calc, 
                  -zeeman_frequency_shift_constant, magnetic_magnitude, 
                  -line.F(), -doppler_constant, 
                  G0, -L0, G2, -L2, e, FVC,
                  derivatives_data, QI,
                  -ddoppler_constant_dT, 
                  dG0_dT, -dL0_dT, dG2_dT, -dL2_dT, de_dT, dFVC_dT);
          break;
        case LineShapeType::ByPressureBroadeningData:
        case LineShapeType::End:
        default:
          throw std::runtime_error("Cannot understand the requested line shape type for mirroring by same line shape.");
      }
      F_calc -= Fm;
      for(Index i = 0; i < dF.nelem(); i++) dF[i][this_xsec_range] -= dFm[i];
      break;
    }
    case MirroringType::End:
    default:
      throw std::runtime_error("Cannot understand the requested mirroring type for mirroring.");
  }
  
  // Line normalization
  switch(line.GetLineNormalizationType())
  {
    case LineNormalizationType::None:
      break;
    case LineNormalizationType::VVH:
      apply_VVH(F_calc, dF, f_grid_calc, line.F(), temperature, derivatives_data, QI, this_xsec_range);
      break;
    case LineNormalizationType::VVW:
      apply_VVW(F_calc, dF, f_grid_calc, line.F(), derivatives_data, QI, this_xsec_range);
      break;
    case LineNormalizationType::RosenkranzQuadratic:
      apply_rosenkranz_quadratic(F_calc, dF, f_grid_calc, line.F(), temperature, derivatives_data, QI, this_xsec_range);
      break;
    case LineNormalizationType::End:
    default:
      throw std::runtime_error("Cannot understand the requested line normalization type.");
  }
  
  // Apply line mixing if relevant
  if(Y not_eq 0 or G not_eq 0)
    apply_linemixing(F_calc, dF, Y, G, 
                    derivatives_data, 
                    QI, dY_dT, dG_dT);
  
  // Change the line integration value
  apply_linestrength(F_calc, dF,  line.I0(), isotopologue_ratio,
                     partition_function_at_temperature, partition_function_at_line_temperature, K1, K2,
                     derivatives_data, QI, dpartition_function_at_temperature_dT, dK1_dT, dK2_dT, dK2_dF0);
  
  // Apply pressure broadening vector if there are matching cases
  if(derivatives_data.get_first_pressure_term() > -1)
    apply_pressurebroadening_jacobian(dF, derivatives_data, QI, pressure_derivatives);
  
  // Non-local thermodynamic equilibrium terms
  if(nlte_temperatures.nelem())
  {
    Numeric Tu, Tl, K4, r_low, dK3_dF0, dK3_dT, dK3_dTl, dK4_dT, dK3_dTu, dK4_dTu;
    
    // These four are set by user on controlfile level
    const Index evlow_index = line.EvlowIndex();
    const Index evupp_index = line.EvuppIndex();
    const Numeric El = line.Evlow();
    const Numeric Eu = line.Evupp();
    
    // If anything is said, the index is set and will match with the input vector
    if(evupp_index > -1)
    {
      Tu = nlte_temperatures[evupp_index];
      
      // Additional emission is from upper state
      K4 = boltzman_ratio(Tu, temperature, Eu);
    }
    else // vibrational temperature is atmospheric temperature and no additional numbers in upper state
    {
      Tu = temperature;
      K4 = 1.0;
    }
    
    // Same as upper state above
    if(evlow_index > -1)
    {
      Tl = nlte_temperatures[evlow_index];
      r_low = boltzman_ratio(Tl, temperature, El);
    }
    else // vibrational temperature is atmospheric temperature and no additional numbers in lower state
    {
      Tl = temperature;
      r_low = 1.0;
    }
    
    // The additional absorption requires the ratio between upper and lower state number distributions
    const Numeric K3 = absorption_nlte_ratio(gamma, K4, r_low);
    
    if(derivatives_data.do_frequency())
      dK3_dF0 = dabsorption_nlte_rate_dF0(gamma, temperature, K4, r_low);
    
    // NOTE:  Having NLTE active AT ALL will change the jacobian because of this part of the code,
    // though this requires setting El and Eu for all lines, though this is not yet default...
    // So if you see this part of the code after having a runtime_error, 
    // you will need to write those functions yourself...
    if(derivatives_data.do_temperature())
      dK3_dT = dabsorption_nlte_rate_dT(gamma, temperature, line.F(), El, Eu, K4, r_low);
    
    if(El > 0)
      dK3_dTl = dabsorption_nlte_rate_dTl(gamma, temperature, Tl, El, r_low);
    
    if(Eu > 0)
    {
      dK3_dTu = dabsorption_nlte_rate_dTu(gamma, temperature, Tu, Eu, K4);
      dK4_dTu = dboltzman_ratio_dT(K4, Tu, Eu);
    }
    
    apply_nonlte(F_calc, dF, N[this_xsec_range], dN, K3, K4, 
                 derivatives_data, QI, 
                 dK3_dT, dK4_dT, dK3_dF0, dK3_dTl, dK3_dTu, dK4_dTu, this_xsec_range);
  }
  
  // Cutoff frequency
  if(need_cutoff)
  {
    if(nlte_temperatures.nelem())
      apply_cutoff(F_calc, dF, N[this_xsec_range], dN,
                  derivatives_data, line,
                  volume_mixing_ratio_of_all_species,
                  nlte_temperatures, pressure, temperature,
                  doppler_constant, partial_pressure,
                  isotopologue_ratio, magnetic_magnitude,
                  ddoppler_constant_dT, pressure_limit_for_linemixing,
                  zeeman_frequency_shift_constant,
                  partition_function_at_temperature,
                  dpartition_function_at_temperature_dT,
                  partition_function_at_line_temperature,
                  broad_spec_locations, this_species_location_in_tags,
                  water_index_location_in_tags, this_xsec_range, verbosity);
    else
      apply_cutoff(F_calc, dF, N, dN,
                  derivatives_data, line,
                  volume_mixing_ratio_of_all_species,
                  nlte_temperatures, pressure, temperature,
                  doppler_constant, partial_pressure,
                  isotopologue_ratio, magnetic_magnitude,
                  ddoppler_constant_dT, pressure_limit_for_linemixing,
                  zeeman_frequency_shift_constant,
                  partition_function_at_temperature,
                  dpartition_function_at_temperature_dT,
                  partition_function_at_line_temperature,
                  broad_spec_locations, this_species_location_in_tags,
                  water_index_location_in_tags, this_xsec_range, verbosity);
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
 * \param nlte_temperatures As name suggests
 * \param pressure As name suggests
 * \param temperature As name suggests
 * \param doppler_constant Frequency-independent part of the Doppler broadening
 * \param partial_pressure Pressure of species that line belongs to at this level
 * \param isotopologue_ratio The ratio of the isotopologue in the atmosphere at this level
 * \param magnetic_magnitude Absolute strength of the magnetic field
 * \param ddoppler_constant_dT Temperature derivative of doppler_constant
 * \param pressure_limit_for_linemixing As WSV lm_p_lim
 * \param zeeman_frequency_shift_constant Zeeman shift parameter for the line
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
                                 ConstVectorView nlte_temperatures,
                                 const Numeric& pressure,
                                 const Numeric& temperature,
                                 const Numeric& doppler_constant,
                                 const Numeric& partial_pressure,
                                 const Numeric& isotopologue_ratio,
                                 const Numeric& magnetic_magnitude,
                                 const Numeric& ddoppler_constant_dT,
                                 const Numeric& pressure_limit_for_linemixing,
                                 const Numeric& zeeman_frequency_shift_constant,
                                 const Numeric& partition_function_at_temperature,
                                 const Numeric& dpartition_function_at_temperature_dT,
                                 const Numeric& partition_function_at_line_temperature,
                                 const ArrayOfIndex& broad_spec_locations,
                                 const Index& this_species_location_in_tags,
                                 const Index& water_index_location_in_tags,
                                 const ComplexRange& df_range,
                                 const Verbosity& verbosity)
{ 
  // Size of derivatives
  const Index nj = dF.nelem(); 
  const Index nn = dN.nelem(); 
  
  // Setup compute variables
  Vector f_grid_cutoff(1, line.F() +line.CutOff());
  ComplexVector Fc(1), Nc(1);
  ArrayOfComplexVector dFc(nj), dNc(nn);
  for(auto& aovc : dFc) aovc.resize(1);
  for(auto& aovc : dNc) aovc.resize(1);
  ComplexRange tmp(joker);
  
  // Recompute the line for a single frequency
  set_cross_section_for_single_line(Fc, dFc, Nc, dNc, tmp,
                                    derivatives_data, line, f_grid_cutoff,
                                    volume_mixing_ratio_of_all_species,
                                    nlte_temperatures, pressure, temperature,
                                    doppler_constant, partial_pressure, isotopologue_ratio,
                                    magnetic_magnitude, ddoppler_constant_dT,
                                    pressure_limit_for_linemixing,
                                    zeeman_frequency_shift_constant,
                                    partition_function_at_temperature,
                                    dpartition_function_at_temperature_dT,
                                    partition_function_at_line_temperature,
                                    broad_spec_locations, this_species_location_in_tags,
                                    water_index_location_in_tags, verbosity, true);
  
  // Apply cutoff values
  F -= Fc[0];
  if(N.nelem())
    N -= N[0];
  for(Index i = 0; i < nj; i++)
    dF[i][df_range] -= dFc[i][0];
  for(Index i = 0; i < nn; i++)
    dN[i][df_range] -= dNc[i][0];
}

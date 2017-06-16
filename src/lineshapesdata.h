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
 * This file contains all handling of individual lines.
 * The reason is that the old methods are cumbersome to adapt and need redesigning
 * 
 * \author Richard Larsson
 * \date   2017-05-16
 */

#ifndef lineshapedata_h
#define lineshapedata_h

#include "complex.h"
#include "linerecord.h"
#include "Faddeeva.hh"
#include "absorption.h"
#include "partial_derivatives.h"
#include <Eigen/Dense>

/*
 * Class to solve the problem
 * 
 *      cross-section of a line equals line strength times line shape
 *      
 *      or
 * 
 *      sigma = S x F,
 * 
 * which means
 *      
 *      dsigma = dS x F + S x dF
 * 
 * TODO: Add QuantumIdentifier input to test for line catalog-parameters derivatives in relevant places
 * 
 * TODO: Add NLTE line strength calculator
 * 
 * TODO: Better combination with Zeeman calculations
 * 
 * TODO: Find work-around for incomplete line-shapes like "Voigt Kuntz"
 */
class LineFunctions
{
public:
  
  enum LineShapeType {
    LST_NONE,                       // Holder that should always fail if evoked  X
    LST_O2_NON_RESONANT,            // Debye lineshape  X
    LST_DOPPLER,                    // Doppler lineshape  X
    LST_LORENTZ,                    // Lorentz lineshape  X
    LST_MIRRORED_LORENTZ,           // Mirrored Lorentz lineshape  X
    LST_FADDEEVA_ALGORITHM916,      // Faddeeva lineshape  X
    LST_HUI_ETAL_1978,              // Faddeeva lineshape  X
    LST_HTP                         // Hartmann-Tran lineshape  X
  };
  
  enum LineShapeNorm {
    LSN_NONE,                      // No normalization  X
    LSN_ROSENKRANZ_QUADRATIC,      // Quadratic normalization (f/f0)^2*h*f0/(2*k*T)/sinh(h*f0/(2*k*T))  X
    LSN_VVW,                       // Van Vleck Weiskopf normalization (f*f) / (f0*f0)  X
    LSN_VVH                        // Van Vleck Huber normalization f*tanh(h*f/(2*k*T))) / (f0*tanh(h*f0/(2*k*T)))  X
  };
  
  void set_lorentz(ComplexVector& F, // Sets the full complex line shape without line mixing
                   ArrayOfComplexVector& dF,
                   const Vector& f_grid,
                   const Numeric& zeeman_df,
                   const Numeric& magnetic_magnitude,
                   const Numeric& F0_noshift,
                   const Numeric& G0,
                   const Numeric& L0,
                   const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                   const Numeric& dG0_dT=0.0,
                   const Numeric& dL0_dT=0.0) const
  { 
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    extern const Numeric PI;
    
    const Numeric invPI = 1.0 / PI;
    
    const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude;
    
    // Signa change of F and F0?
    const Complex denom0 = Complex(G0, F0);
    
    F = invPI;
    
    for(Index ii = 0; ii < nf; ii++)
    {
      F[ii] /= (denom0 - Complex(0.0, f_grid[ii]));
    }
    
    if(nppd > 0)
    {
      dF[nppd-1] = F;
      dF[nppd-1] *= F;
      dF[nppd-1] *= PI;
    }
    
    for(Index jj = 0; jj < nppd; jj++)
    {
      
      if(derivatives_data(jj) == JQT_temperature)
      {
        dF[jj] = dF[nppd-1];
        dF[jj] *= Complex(dG0_dT, dL0_dT);
      }
      else if(derivatives_data(jj) == JQT_frequency or
        derivatives_data(jj) == JQT_wind_magnitude or
        derivatives_data(jj) == JQT_wind_u or
        derivatives_data(jj) == JQT_wind_v or
        derivatives_data(jj) == JQT_wind_w)
      {
        dF[jj] = dF[nppd-1];
        dF[jj] *= Complex(0.0, -1.0);
      }
      else if(derivatives_data(jj) == JQT_line_center)
      {
        dF[jj] = dF[nppd-1];
        dF[jj] *= Complex(0.0, 1.0);
      }
      else if(derivatives_data(jj) == JQT_line_gamma_self or 
        derivatives_data(jj) == JQT_line_gamma_selfexponent or
        derivatives_data(jj) == JQT_line_pressureshift_self or
        derivatives_data(jj) == JQT_line_gamma_foreign or
        derivatives_data(jj) == JQT_line_gamma_foreignexponent or
        derivatives_data(jj) == JQT_line_pressureshift_foreign or
        derivatives_data(jj) == JQT_line_gamma_water or
        derivatives_data(jj) == JQT_line_gamma_waterexponent or 
        derivatives_data(jj) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
      {
        dF[jj] = dF[nppd-1];
        dF[jj] *= Complex(1.0, 1.0);
      }
      else if(derivatives_data(jj) == JQT_magnetic_magntitude)// No external inputs --- errors because of frequency shift when Zeeman is used?
      {
        dF[jj] = dF[nppd-1];
        dF[jj] *= Complex(0.0, zeeman_df);
      }
    }
  }
  
  void set_mirrored_lorentz(ComplexVector& F, // Sets the full complex line shape without line mixing
                            ArrayOfComplexVector& dF, 
                            const Vector& f_grid,
                            const Numeric& zeeman_df,
                            const Numeric& magnetic_magnitude,
                            const Numeric& F0_noshift,
                            const Numeric& G0,
                            const Numeric& L0,
                            const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                            const Numeric& dG0_dT=0.0,
                            const Numeric& dL0_dT=0.0) const
  { 
    Complex Fplus, Fminus;
    
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    extern const Numeric PI;
    
    const Numeric invPI = 1.0 / PI;
    
    const Numeric F0 = F0_noshift + L0 + zeeman_df * magnetic_magnitude;
    
    // Signa change of F and F0?
    const Complex denom0 = Complex(G0, F0), denom1 = Complex(G0, -F0);
    
    for(Index ii = 0; ii < nf; ii++)
    {
      Fplus = invPI / (denom0 - Complex(0.0, f_grid[ii]));
      Fminus = invPI / (denom1 - Complex(0.0, f_grid[ii]));
      F[ii] = Fplus + Fminus;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        if(derivatives_data(jj) == JQT_temperature)
        {
          dF[jj][ii] = (Fplus * Fplus * Complex(dG0_dT, dL0_dT) + Fminus * Fminus * Complex(dG0_dT, -dL0_dT)) * PI;
        }
        else if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w)
        {
          dF[jj][ii] = (Fplus * Fplus + Fminus * Fminus) * PI;
          dF[jj][ii] *= Complex(0.0, -1.0);
        }
        else if(derivatives_data(jj) == JQT_line_center)
        {
          dF[jj][ii] = (Fplus * Fplus * Complex(0.0, 1.0) + Fminus * Fminus * Complex(0.0, -1.0)) * PI;
        }
        else if(derivatives_data(jj) == JQT_line_gamma_self or 
          derivatives_data(jj) == JQT_line_gamma_selfexponent or
          derivatives_data(jj) == JQT_line_pressureshift_self or
          derivatives_data(jj) == JQT_line_gamma_foreign or
          derivatives_data(jj) == JQT_line_gamma_foreignexponent or
          derivatives_data(jj) == JQT_line_pressureshift_foreign or
          derivatives_data(jj) == JQT_line_gamma_water or
          derivatives_data(jj) == JQT_line_gamma_waterexponent or 
          derivatives_data(jj) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
        {
          dF[jj][ii] = (Fplus * Fplus * Complex(1.0, 1.0) + Fminus * Fminus * Complex(1.0, -1.0)) * PI;
        }
        else if(derivatives_data(jj) == JQT_magnetic_magntitude)// No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          dF[jj][ii] = (Fplus * Fplus * Complex(0.0, zeeman_df) + Fminus * Fminus * Complex(0.0, -zeeman_df)) * PI;
        }
      }
    }
  }
  
  void set_doppler(ComplexVector& F, // Sets the full complex line shape without line mixing
                   ArrayOfComplexVector& dF,
                   const Vector& f_grid,
                   const Numeric& zeeman_df,
                   const Numeric& magnetic_magnitude,
                   const Numeric& F0_noshift,
                   const Numeric& GD_div_F0,
                   const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                   const Numeric& dGD_div_F0_dT=0.0) const
  {
    set_faddeeva_algorithm916(F, 
                              dF, 
                              f_grid, 
                              zeeman_df, 
                              magnetic_magnitude, 
                              F0_noshift, 
                              GD_div_F0, 0.0, 0.0,
                              derivatives_data, 
                              dGD_div_F0_dT);
  }
  
  void set_htp(ComplexVector& F, // Sets the full complex line shape without line mixing
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
               const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
               const Numeric& dGD_div_F0_dT=0.0,
               const Numeric& dG0_dT=0.0,
               const Numeric& dL0_dT=0.0,
               const Numeric& dG2_dT=0.0,
               const Numeric& dL2_dT=0.0,
               const Numeric& deta_dT=0.0,
               const Numeric& dFVC_dT=0.0) const
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
    extern const Numeric PI;
    static const Complex i(0.0, 1.0), one_plus_one_i(1.0, 1.0);
    static const Numeric ln2 = log(2.0), sqrtLN2 = sqrt(ln2); 
    static const Numeric sqrtPI = sqrt(PI), invPI = 1.0 / PI, invSqrtPI = sqrt(invPI);
    
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
    for(Index ii = 0; ii < nf; ii++)
    {
      // Relative frequency
      X = (C0t - i*(F0 - f_grid[ii])) * invC2t;
      
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
      F[ii] = A * invG * invPI;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w) // No external inputs
        { 
          // If this is the first time it is calculated this frequency bin, do the full calculation
          if(first_frequency == jj)
          {
            dX = invC2t * i;
            dZp = dZm = 0.5 * dX / sqrtXY;
            
            dwiZm = dwiZp = i * invSqrtPI;
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
            
            dF[jj][ii] = invG * (invPI * dA - F[ii] * dg); 
          }
          else  // copy for repeated occurences
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_temperature)
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
          
          dwiZm = dwiZp = i * invSqrtPI;
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
          
          dF[jj][ii] = invG * (invPI * dA - F[ii] * dg); 
        }
        else if(derivatives_data(jj) == JQT_line_center) // No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          dY = GD / (2.0 * ln2) * invC2t*invC2t * GD_div_F0 * F0_ratio;
          dX = - i * invC2t * F0_ratio;
          
          dZm = dZp = 0.5 * (dX + dY) / sqrtXY;
          dZm -= 0.5 / sqrtY * dY;
          dZp += 0.5 / sqrtY * dY;
          
          dwiZm = dwiZp = i * invSqrtPI;
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
          
          dF[jj][ii] = invG * (invPI * dA - F[ii] * dg); 
        }
        else if(derivatives_data(jj) == JQT_line_gamma_self or 
          derivatives_data(jj) == JQT_line_gamma_selfexponent or
          derivatives_data(jj) == JQT_line_pressureshift_self or
          derivatives_data(jj) == JQT_line_gamma_foreign or
          derivatives_data(jj) == JQT_line_gamma_foreignexponent or
          derivatives_data(jj) == JQT_line_pressureshift_foreign or
          derivatives_data(jj) == JQT_line_gamma_water or
          derivatives_data(jj) == JQT_line_gamma_waterexponent or 
          derivatives_data(jj) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == jj)
          {
            dC0t = one_minus_eta * one_plus_one_i;
            dX = invC2t * dC0t;
            dZm = dZp = 0.5 * dX / sqrtXY;
            
            dwiZm = dwiZp = i * invSqrtPI;
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
            
            dF[jj][ii] = invG * (invPI * dA - F[ii] * dg); 
          }
          else  // copy for repeated occurences
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_magnetic_magntitude)// No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          dY = GD / (2.0 * ln2) * invC2t*invC2t * GD_div_F0 * (1.0 - F0_ratio) / magnetic_magnitude;
          dX = - i * invC2t * (1.0 - F0_ratio) / magnetic_magnitude;
          
          dZm = dZp = 0.5 * (dX + dY) / sqrtXY;
          dZm -= 0.5 / sqrtY * dY;
          dZp += 0.5 / sqrtY * dY;
          
          dwiZm = dwiZp = i * invSqrtPI;
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
          
          dF[jj][ii] = invG * (invPI * dA - F[ii] * dg); 
        }
      }
    }
  }
  
  void set_faddeeva_algorithm916(ComplexVector& F, // Sets the full complex line shape without line mixing
                                 ArrayOfComplexVector& dF,
                                 const Vector& f_grid,
                                 const Numeric& zeeman_df,
                                 const Numeric& magnetic_magnitude,
                                 const Numeric& F0_noshift,
                                 const Numeric& GD_div_F0,
                                 const Numeric& G0,
                                 const Numeric& L0,
                                 const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                 const Numeric& dGD_div_F0_dT=0.0,
                                 const Numeric& dG0_dT=0.0,
                                 const Numeric& dL0_dT=0.0) const
  {
    // For calculations
    Numeric dx;
    Complex w, z, dw_over_dz, dz;
    
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    // PI
    extern const Numeric PI;
    
    // constant sqrt(1/pi)
    static const Numeric sqrt_invPI =  sqrt(1.0/PI);
    
    // Doppler broadening and line center
    const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + L0;
    const Numeric GD = GD_div_F0 * F0;
    const Numeric invGD = 1.0 / GD;
    const Numeric dGD_dT = dGD_div_F0_dT * F0;
    
    // constant normalization factor for voigt
    const Numeric fac = sqrt_invPI * invGD;
    F = fac;
    
    // Ratio of the Lorentz halfwidth to the Doppler halfwidth
    const Complex z0 = Complex(-F0, G0) * invGD;
    
    const Index first_pressure_broadening = derivatives_data.get_first_pressure_term(),
    first_frequency = derivatives_data.get_first_frequency();
    
    // frequency in units of Doppler
    for (Index ii=0; ii<nf; ii++)
    {
      dx = f_grid[ii] * invGD;
      z = z0 + dx;
      w = Faddeeva::w(z);
      
      F[ii] *= w;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        if(jj==0)
          dw_over_dz = 2.0 * (z * w - sqrt_invPI);
        
        if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w) // No external inputs
        { 
          // If this is the first time it is calculated this frequency bin, do the full calculation
          if(first_frequency == jj)
          {
            //dz = Complex(invGD, 0.0);
            
            dF[jj][ii] = fac * dw_over_dz * invGD; //dz; 
          }
          else  // copy for repeated occurences
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_temperature)
        {
          dz = Complex(-dL0_dT, dG0_dT) - z * dGD_dT;
          
          dF[jj][ii] = -F[ii] * dGD_dT;
          dF[jj][ii] += fac * dw_over_dz * dz;
          dF[jj][ii] *= invGD;
        }
        else if(derivatives_data(jj) == JQT_line_center) // No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[jj][ii] = -F[ii] * GD_div_F0;
          dF[jj][ii] += dw_over_dz * dz;
          dF[jj][ii] *= fac * invGD;
        }
        else if(derivatives_data(jj) == JQT_line_gamma_self or 
          derivatives_data(jj) == JQT_line_gamma_selfexponent or
          derivatives_data(jj) == JQT_line_pressureshift_self or
          derivatives_data(jj) == JQT_line_gamma_foreign or
          derivatives_data(jj) == JQT_line_gamma_foreignexponent or
          derivatives_data(jj) == JQT_line_pressureshift_foreign or
          derivatives_data(jj) == JQT_line_gamma_water or
          derivatives_data(jj) == JQT_line_gamma_waterexponent or 
          derivatives_data(jj) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == jj)
          {
            dz = Complex(-1.0, 1.0) * invGD;
            dF[jj][ii] = fac * dw_over_dz * dz; 
          }
          else
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_magnetic_magntitude)// No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          // dz = Complex(- zeeman_df * invGD, 0.0);
          
          dF[jj][ii] = fac * dw_over_dz * (- zeeman_df * invGD); //* dz; 
        }
        else if(derivatives_data(jj) == JQT_line_mixing_DF or
          derivatives_data(jj) == JQT_line_mixing_DF0 or
          derivatives_data(jj) == JQT_line_mixing_DF1 or
          derivatives_data(jj) == JQT_line_mixing_DFexp)
        {
          // dz = Complex(-invGD, 0.0);
          
          dF[jj][ii] = fac * dw_over_dz * (-invGD); //* dz;
        }
      }
      
    }
  }
  
  void set_faddeeva_from_full_linemixing(ComplexVector& F, // Sets the full complex line shape without line mixing
                                        ArrayOfComplexVector& dF,
                                        const Vector& f_grid,
                                        const Complex& eigenvalue_no_shift,
                                        const Numeric& GD_div_F0,
                                        const Numeric& L0,
                                        const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                        const Numeric& dGD_div_F0_dT=0.0,
                                        const Complex& deigenvalue_dT=0.0,
                                        const Numeric& dL0_dT=0.0) const
  {
    // For calculations
    Numeric dx;
    Complex w, z, dw_over_dz, dz;
    
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    // PI
    extern const Numeric PI;
    
    // constant sqrt(1/pi)
    static const Numeric sqrt_invPI =  sqrt(1.0/PI);
    
    // Doppler broadening and line center
    const Complex eigenvalue = eigenvalue_no_shift + L0;
    const Numeric F0 = eigenvalue.real() + L0;
    const Numeric GD = GD_div_F0 * F0;
    const Numeric invGD = 1.0 / GD;
    const Numeric dGD_dT = dGD_div_F0_dT * F0;
    
    // constant normalization factor for voigt
    const Numeric fac = sqrt_invPI * invGD;
    F = fac;
    
    // Ratio of the Lorentz halfwidth to the Doppler halfwidth
    const Complex z0 = -eigenvalue * invGD;
    
    const Index first_pressure_broadening = derivatives_data.get_first_pressure_term(),
    first_frequency = derivatives_data.get_first_frequency();
    
    // frequency in units of Doppler
    for (Index ii=0; ii<nf; ii++)
    {
      dx = f_grid[ii] * invGD;
      z = z0 + dx;
      w = Faddeeva::w(z);
      
      F[ii] *= w;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        if(jj==0)
          dw_over_dz = 2.0 * (z * w - sqrt_invPI);
        
        if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w) // No external inputs
        { 
          // If this is the first time it is calculated this frequency bin, do the full calculation
          if(first_frequency == jj)
          {
            //dz = Complex(invGD, 0.0);
            
            dF[jj][ii] = fac * dw_over_dz * invGD; //dz; 
          }
          else  // copy for repeated occurences
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_temperature)
        {
          dz = (deigenvalue_dT - dL0_dT) - z * dGD_dT;
          
          dF[jj][ii] = -F[ii] * dGD_dT;
          dF[jj][ii] += fac * dw_over_dz * dz;
          dF[jj][ii] *= invGD;
        }
        else if(derivatives_data(jj) == JQT_line_center) // No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[jj][ii] = -F[ii] * GD_div_F0;
          dF[jj][ii] += dw_over_dz * dz;
          dF[jj][ii] *= fac * invGD;
        }
        else if(derivatives_data(jj) == JQT_line_gamma_self or 
          derivatives_data(jj) == JQT_line_gamma_selfexponent or
          derivatives_data(jj) == JQT_line_pressureshift_self or
          derivatives_data(jj) == JQT_line_gamma_foreign or
          derivatives_data(jj) == JQT_line_gamma_foreignexponent or
          derivatives_data(jj) == JQT_line_pressureshift_foreign or
          derivatives_data(jj) == JQT_line_gamma_water or
          derivatives_data(jj) == JQT_line_gamma_waterexponent or 
          derivatives_data(jj) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == jj)
          {
            dz = Complex(-1.0, 1.0) * invGD;
            dF[jj][ii] = fac * dw_over_dz * dz; 
          }
          else
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_line_mixing_DF or
          derivatives_data(jj) == JQT_line_mixing_DF0 or
          derivatives_data(jj) == JQT_line_mixing_DF1 or
          derivatives_data(jj) == JQT_line_mixing_DFexp)
        {
          // dz = Complex(-invGD, 0.0);
          
          dF[jj][ii] = fac * dw_over_dz * (-invGD); //* dz;
        }
      }
      
    }
  }
  
  void set_hui_etal_1978(ComplexVector& F, // Sets the full complex line shape without line mixing
                         ArrayOfComplexVector& dF,
                         const Vector& f_grid,
                         const Numeric& zeeman_df,
                         const Numeric& magnetic_magnitude,
                         const Numeric& F0_noshift,
                         const Numeric& GD_div_F0,
                         const Numeric& G0,
                         const Numeric& L0,
                         const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                         const Numeric& dGD_div_F0_dT=0.0,
                         const Numeric& dG0_dT=0.0,
                         const Numeric& dL0_dT=0.0) const
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
    
    // PI
    extern const Numeric PI;
    
    // constant sqrt(1/pi)
    static const Numeric sqrt_invPI =  sqrt(1/PI);
    
    // Doppler broadening and line center
    const Numeric F0 = F0_noshift + zeeman_df * magnetic_magnitude + L0;
    const Numeric GD = GD_div_F0 * F0;
    const Numeric invGD = 1.0 / GD;
    const Numeric dGD_dT = dGD_div_F0_dT * F0;
    
    // constant normalization factor for voigt
    const Numeric fac = sqrt_invPI * invGD;
    
    // Ratio of the Lorentz halfwidth to the Doppler halfwidth
    const Complex z0 = Complex(-F0, G0) * invGD;
    
    const Index first_pressure_broadening = derivatives_data.get_first_pressure_term(),
    first_frequency = derivatives_data.get_first_frequency();
    
    // frequency in units of Doppler
    for (Index ii=0; ii<nf; ii++)
    {
      dx = f_grid[ii] * invGD;
      z = z0 + dx;
      /*  Since this is ported from old FORTRAN code, intention must be that A and B coeffs are floats?
      A = (((((.5641896*z+5.912626)*z+30.18014)*z+
      93.15558)*z+181.9285)*z+214.3824)*z+122.6079;
      B = ((((((z+10.47986)*z+53.99291)*z+170.3540)*z+
      348.7039)*z+457.3345)*z+352.7306)*z+122.6079;
       */
      A = (((((.5641896f*z+5.912626f)*z+30.18014f)*z+
      93.15558f)*z+181.9285f)*z+214.3824f)*z+122.6079f;
      B = ((((((z+10.47986f)*z+53.99291f)*z+170.3540f)*z+
      348.7039f)*z+457.3345f)*z+352.7306f)*z+122.6079f;
      w = A / B;
      
      F[ii] = fac * w;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        if(jj==0)
          dw_over_dz = 2.0 * (z * w - sqrt_invPI);
        
        if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w) // No external inputs
        { 
          // If this is the first time it is calculated this frequency bin, do the full calculation
          if(first_frequency == jj)
          {
            //dz = Complex(invGD, 0.0);
            
            dF[jj][ii] = fac * dw_over_dz * invGD; //dz; 
          }
          else  // copy for repeated occurences
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_temperature)
        {
          dz = Complex(-dL0_dT, dG0_dT) - z * dGD_dT;
          
          dF[jj][ii] = -F[ii] * dGD_dT;
          dF[jj][ii] += fac * dw_over_dz * dz;
          dF[jj][ii] *= invGD;
        }
        else if(derivatives_data(jj) == JQT_line_center) // No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          dz = -z * GD_div_F0 - 1.0;
          
          dF[jj][ii] = -F[ii] * GD_div_F0;
          dF[jj][ii] += dw_over_dz * dz;
          dF[jj][ii] *= fac * invGD;
        }
        else if(derivatives_data(jj) == JQT_line_gamma_self or 
          derivatives_data(jj) == JQT_line_gamma_selfexponent or
          derivatives_data(jj) == JQT_line_pressureshift_self or
          derivatives_data(jj) == JQT_line_gamma_foreign or
          derivatives_data(jj) == JQT_line_gamma_foreignexponent or
          derivatives_data(jj) == JQT_line_pressureshift_foreign or
          derivatives_data(jj) == JQT_line_gamma_water or
          derivatives_data(jj) == JQT_line_gamma_waterexponent or 
          derivatives_data(jj) == JQT_line_pressureshift_water) // Only the zeroth order terms --- the derivative with respect to these have to happen later
        {
          // Note that if the species vmr partial derivative is necessary here is where it goes
          if(first_pressure_broadening == jj)
          {
            dz = Complex(-1.0, 1.0) * invGD;
            dF[jj][ii] = fac * dw_over_dz * dz; 
          }
          else
          {
            dF[jj][ii] = dF[first_frequency][ii]; 
          }
        }
        else if(derivatives_data(jj) == JQT_magnetic_magntitude)// No external inputs --- errors because of frequency shift when Zeeman is used?
        {
          // dz = Complex(- zeeman_df * invGD, 0.0);
          
          dF[jj][ii] = fac * dw_over_dz * (- zeeman_df * invGD); //* dz; 
        }
        else if(derivatives_data(jj) == JQT_line_mixing_DF or
          derivatives_data(jj) == JQT_line_mixing_DF0 or
          derivatives_data(jj) == JQT_line_mixing_DF1 or
          derivatives_data(jj) == JQT_line_mixing_DFexp)
        {
          // dz = Complex(-invGD, 0.0);
          
          dF[jj][ii] = fac * dw_over_dz * (-invGD); //* dz;
        }
      } 
    }
  }
  
  void set_o2_non_resonant(ComplexVector& F, // Sets the full complex line shape without line mixing
                           ArrayOfComplexVector& dF,
                           const Vector& f_grid,
                           const Numeric& F0,
                           const Numeric& G0,
                           const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                           const Numeric& dG0_dT=0.0)
  {
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();

    extern const Numeric PI;
    
    const Numeric fac = G0/PI, invG0 = 1.0/G0;
    
    for(Index ii = 0; ii < nf; ii++)
    {
      F[ii] =  fac / ((f_grid[ii]-F0) * (f_grid[ii]-F0));
    }
    
    for(Index jj = 0; jj < nppd; jj++)
    {
      if(derivatives_data(jj) == JQT_temperature)
      {
        dF[jj] = F;
        dF[jj] *= dG0_dT * invG0;
      }
      else if(derivatives_data(jj) == JQT_frequency or
        derivatives_data(jj) == JQT_wind_magnitude or
        derivatives_data(jj) == JQT_wind_u or
        derivatives_data(jj) == JQT_wind_v or
        derivatives_data(jj) == JQT_wind_w)
      {
        dF[jj] = F;
        for(Index ii = 0; ii < nf; ii++)
        {
          dF[jj][ii] /= -2.0 / (f_grid[ii] - F0);
        }
      }
      else if(derivatives_data(jj) == JQT_line_center)
      {
        dF[jj] = F;
        for(Index ii = 0; ii < nf; ii++)
        {
          dF[jj][ii] *= 2.0 / (f_grid[ii] - F0);
        }
      }
      else if(derivatives_data(jj) == JQT_line_gamma_self or 
        derivatives_data(jj) == JQT_line_gamma_selfexponent or
        derivatives_data(jj) == JQT_line_gamma_foreign or
        derivatives_data(jj) == JQT_line_gamma_foreignexponent or
        derivatives_data(jj) == JQT_line_gamma_water or
        derivatives_data(jj) == JQT_line_gamma_waterexponent)
      {
        dF[jj] = F;
        dF[jj] *= invG0;
      }
    }
  }
  
  void apply_linemixing(ComplexVector& F, // Returns the full complex or normalized line shape with line mixing
                        ArrayOfComplexVector& dF,
                        const Numeric& Y,
                        const Numeric& G,
                        const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                        const Numeric& dY_dT=0.0,
                        const Numeric& dG_dT=0.0) const
  {
    const Index nf = F.nelem(), nppd = derivatives_data.nelem();
    
    const Complex LM = Complex(1.0 + G, Y);
    const Complex dLM_dT = Complex(dG_dT, dY_dT);
    
    for(Index iq = 0; iq < nppd; iq++)
    {
      if(derivatives_data(iq) == JQT_temperature)
      {
        dF[iq] *= LM;
        for(Index ii = 0; ii < nf; ii++)
        {
          dF[iq][ii] += F[ii] * dLM_dT;
        }
      }
      else if(derivatives_data(iq) == JQT_line_mixing_Y or
        derivatives_data(iq) == JQT_line_mixing_Y0 or 
        derivatives_data(iq) == JQT_line_mixing_Y1 or 
        derivatives_data(iq) == JQT_line_mixing_Yexp)
      {
        dF[iq] = F;
        dF[iq] *= Complex(0.0, 1.0);
      }
      else if(derivatives_data(iq) == JQT_line_mixing_G or
        derivatives_data(iq) == JQT_line_mixing_G0 or 
        derivatives_data(iq) == JQT_line_mixing_G1 or 
        derivatives_data(iq) == JQT_line_mixing_Gexp)
      {
        dF[iq] = F;
      }
      else
      {
        dF[iq] *= LM;
      }
    }
    
    F *= LM;
  }
  
  void apply_rosenkranz_quadratic(ComplexVector& F, // Returns as normalized complex line shape with or without line mixing
                                  ArrayOfComplexVector& dF,
                                  const Vector& f_grid,
                                  const Numeric& F0, // Only line center without any shifts
                                  const Numeric& T,
                                  const PropmatPartialsData& derivatives_data=PropmatPartialsData()) const
  {
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const Numeric invF0 = 1.0/F0;
    const Numeric mafac = (PLANCK_CONST) / (2.000e0 * BOLTZMAN_CONST * T) /
    sinh((PLANCK_CONST * F0) / (2.000e0 * BOLTZMAN_CONST * T)) * invF0;
    
    Numeric dmafac_dF0_div_fun, dmafac_dT_div_fun;
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
    
    for (Index ii=0; ii < nf; ii++)
    {
      fun = mafac * (f_grid[ii] * f_grid[ii]);
      F[ii] *= fun;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        dF[jj][ii] *= fun;
        if(derivatives_data(jj) == JQT_temperature)
        {
          dF[jj][ii] += dmafac_dT_div_fun * F[ii];
        }
        else if(derivatives_data(jj) == JQT_line_center)
        {
          dF[jj][ii] += dmafac_dF0_div_fun * F[ii];
        }
        else if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w)
        {
          dF[jj][ii] += (2.0 / f_grid[ii]) * F[ii];
        }
      }
    }
  }
  
  void apply_VVH(ComplexVector& F, // Returns as normalized complex line shape with or without line mixing
                 ArrayOfComplexVector& dF,
                 const Vector& f_grid,
                 const Numeric& F0, // Only line center without any shifts
                 const Numeric& T,
                 const PropmatPartialsData& derivatives_data=PropmatPartialsData()) const
  {
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    // 2kT is constant for the loop
    const Numeric kT = 2.0 * BOLTZMAN_CONST * T;
    
    // denominator is constant for the loop
    const Numeric absF0 = abs(F0);
    const Numeric tanh_f0part = tanh(PLANCK_CONST * absF0 / kT);
    const Numeric denom = absF0 * tanh_f0part;
    
    Numeric fac_df0;
    if(derivatives_data.do_line_center())
      fac_df0 = (-1.0/absF0 + PLANCK_CONST*tanh_f0part/(kT) - PLANCK_CONST/(kT*tanh_f0part)) * F0/absF0;
    
    for(Index ii=0; ii < nf; ii++)
    {
      const Numeric tanh_fpart = tanh( PLANCK_CONST * f_grid[ii] / kT );
      const Numeric fun = f_grid[ii] * tanh_fpart / denom;
      F[ii] *= fun;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        dF[jj][ii] *= fun;
        if(derivatives_data(jj) == JQT_temperature)
        {
          dF[jj][ii] += (-PLANCK_CONST*(denom - absF0/tanh_f0part - 
          f_grid[ii]*tanh_fpart + f_grid[ii]/tanh_fpart)/(kT*T)) * F[ii];
        }
        else if(derivatives_data(jj) == JQT_line_center)
        {
          dF[jj][ii] += fac_df0 * F[ii];
        }
        else if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w)
        {
          dF[jj][ii] += (1.0/f_grid[ii] -PLANCK_CONST*tanh_fpart/kT + PLANCK_CONST/(kT*tanh_fpart)) * F[ii];
        }
      }
    }
  }
  
  void apply_VVW(ComplexVector& F, // Returns as normalized line shape with or without line mixing
                 ArrayOfComplexVector& dF,
                 const Vector& f_grid,
                 const Numeric& F0, // Only line center without any shifts
                 const PropmatPartialsData& derivatives_data=PropmatPartialsData()) const
  {
    const Index nf = f_grid.nelem(), nppd = derivatives_data.nelem();
    
    // denominator is constant for the loop
    const Numeric invF0 = 1.0 / F0;
    const Numeric invF02 = invF0 * invF0;
    
    for(Index ii = 0; ii < nf; ii++)
    {
      const Numeric fun = f_grid[ii] * f_grid[ii] * invF02;
      F[ii] *= fun;
      
      for(Index jj = 0; jj < nppd; jj++)
      {
        dF[jj][ii] *= fun;
        if(derivatives_data(jj) == JQT_line_center)
        {
          dF[jj][ii] -= 2.0 * invF0 * F[ii];
        }
        else if(derivatives_data(jj) == JQT_frequency or
          derivatives_data(jj) == JQT_wind_magnitude or
          derivatives_data(jj) == JQT_wind_u or
          derivatives_data(jj) == JQT_wind_v or
          derivatives_data(jj) == JQT_wind_w)
        {
          dF[jj][ii] += 2.0 * f_grid[ii] * F[ii];
        }
      }
    }
  }
  
  void apply_linestrength(ComplexVector& F, // Returns the full complex line shape with or without line mixing
                          ArrayOfComplexVector& dF,
                          const Numeric& S0,
                          const Numeric& isotopic_ratio,
                          const Numeric& QT,
                          const Numeric& QT0,
                          const Numeric& K1,
                          const Numeric& K2,
                          //const Numeric& K3,  Add again when NLTE is considered
                          //const Numeric& K4,  Add again when NLTE is considered
                          const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                          const Numeric& dQT_dT=0.0,
                          const Numeric& dK1_dT=0.0,
                          const Numeric& dK2_dT=0.0,
                          const Numeric& dK2_dF0=0.0) const
  {
    const Index nppd = derivatives_data.nelem();
    
    const Numeric invQT = 1.0/QT;
    const Numeric QT_ratio = QT0 * invQT;
    const Numeric S = S0 * isotopic_ratio * QT_ratio * K1 * K2;
    
    for(Index jj = 0; jj < nppd; jj++)
    {
      if(derivatives_data(jj) == JQT_temperature)
      {
        Eigen::VectorXcd eig_dF = MapToEigen(dF[jj]);
        
        eig_dF *= S;
        eig_dF += MapToEigen(F) * (S0 * isotopic_ratio * QT_ratio * 
                                  (K1 * dK2_dT + 
                                   dK1_dT * K2 - 
                                   invQT * dQT_dT * K1 * K2));
      }
      else if(derivatives_data(jj) == JQT_line_strength)
      {
        dF[jj] = F;
        dF[jj] *= isotopic_ratio;
      }
      else if(derivatives_data(jj) == JQT_line_center)
      {
        Eigen::VectorXcd eig_dF = MapToEigen(dF[jj]);
        
        eig_dF *= S;
        eig_dF += MapToEigen(F) * (S0 * isotopic_ratio * QT_ratio * K1 * dK2_dF0);
      }
      else
      {
        dF[jj] *= S;
      }
    }
    F *= S;
  }
  
  void apply_linestrength_from_full_linemixing(ComplexVector& F, // Returns the full complex line shape with line mixing
                                               ArrayOfComplexVector& dF,
                                               const Numeric& F0,
                                               const Numeric& T,
                                               const Complex& S_LM,
                                               const Numeric& isotopic_ratio,
                                               const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                                               const Complex& dS_LM_dT=0.0) const
  {
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    static const Numeric C1 = - PLANCK_CONST / BOLTZMAN_CONST;
    
    const Index nppd = derivatives_data.nelem();
    
    const Numeric invT = 1.0 / T, 
      F0_invT = F0 * invT,
      exp_factor = exp(C1 * F0_invT), 
      f0_factor = F0 * (1.0 - exp_factor); 
      
    const Complex S = S_LM * f0_factor * isotopic_ratio;
      
    for(Index jj = 0; jj < nppd; jj++)
    {
      if(derivatives_data(jj) == JQT_temperature)
      {
        Eigen::VectorXcd eig_dF = MapToEigen(dF[jj]);
        
        eig_dF *= S_LM;
        eig_dF += MapToEigen(F) * (dS_LM_dT * f0_factor + 
                                   S_LM * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio;
      }
      else if(derivatives_data(jj) == JQT_line_strength)
      {
        throw std::runtime_error("Not working yet");
      }
      else if(derivatives_data(jj) == JQT_line_center)
      {
        throw std::runtime_error("Not working yet");
      }
      else
      {
        dF[jj] *= S;
      }
    }
    
    F *= S;
  }
  
  void apply_dipole(ComplexVector& F, // Returns the full complex line shape without line mixing
                    ArrayOfComplexVector& dF, // Returns the derivatives of the full line shape for list_of_derivatives
                    const Numeric& F0,
                    const Numeric& T,
                    const Numeric& d0,
                    const Numeric& rho,
                    const Numeric& isotopic_ratio,
                    const PropmatPartialsData& derivatives_data=PropmatPartialsData(),
                    const Numeric& drho_dT=0.0) const
  {
    // Output is d0^2 * rho * F * isotopic_ratio * F0 * (1-e^(hF0/kT))
    
    extern const Numeric PLANCK_CONST;
    extern const Numeric BOLTZMAN_CONST;
    static const Numeric C1 = - PLANCK_CONST / BOLTZMAN_CONST;
    
    const Index nppd = derivatives_data.nelem();
    
    const Numeric S = d0 * d0 * rho * isotopic_ratio, 
      invT = 1.0 / T, 
      F0_invT = F0 * invT,
      exp_factor = exp(C1 * F0_invT), 
      f0_factor = F0 * (1.0 - exp_factor);
    
    for(Index jj = 0; jj < nppd; jj++)
    {
      if(derivatives_data(jj) == JQT_temperature)
      {
        ComplexMatrixViewMap eig_dF = MapToEigen(dF[jj]);
        
        eig_dF *= S * f0_factor;
        eig_dF += MapToEigen(F) * (d0 * d0 * (drho_dT * f0_factor + 
          rho * C1 * F0_invT * F0_invT * exp_factor) * isotopic_ratio);
      }
      else if(derivatives_data(jj) == JQT_line_center)
      {
        throw std::runtime_error("Not working yet");
      }
      else
      {
        dF[jj] *= S * f0_factor;
      }
    }
    
    F *= S * f0_factor;
  }
                          
  void apply_pressurebroadening_jacobian(ArrayOfComplexVector& dF,
                                         const PropmatPartialsData& derivatives_data,
                                         const Numeric& dgamma_dline_gamma_self=0.0,
                                         const Numeric& dgamma_dline_gamma_foreign=0.0,
                                         const Numeric& dgamma_dline_gamma_water=0.0,
                                         const Numeric& dgamma_dline_pressureshift_self=0.0,
                                         const Numeric& dgamma_dline_pressureshift_foreign=0.0,
                                         const Numeric& dgamma_dline_pressureshift_water=0.0,
                                         const Complex& dgamma_dline_gamma_selfexponent=0.0,
                                         const Complex& dgamma_dline_gamma_foreignexponent=0.0,
                                         const Complex& dgamma_dline_gamma_waterexponent=0.0)
  {
    const Index nppd = derivatives_data.nelem();
    
    for(Index jj = 0; jj < nppd; jj++)
    {
      if(derivatives_data(jj) == JQT_line_gamma_self)
      {
        dF[jj] *= dgamma_dline_gamma_self;
      }
      else if(derivatives_data(jj) == JQT_line_gamma_foreign)
      {
        dF[jj] *= dgamma_dline_gamma_foreign;
      }
      else if(derivatives_data(jj) == JQT_line_gamma_water)
      {
        dF[jj] *= dgamma_dline_gamma_water;
      }
      else if(derivatives_data(jj) == JQT_line_pressureshift_self)
      {
        dF[jj] *= Complex(0.0, dgamma_dline_pressureshift_self);
      }
      else if(derivatives_data(jj) == JQT_line_pressureshift_foreign)
      {
        dF[jj] *= Complex(0.0, dgamma_dline_pressureshift_foreign);
      }
      else if(derivatives_data(jj) == JQT_line_pressureshift_water)
      {
        dF[jj] *= Complex(0.0, dgamma_dline_pressureshift_water);
      }
      else if(derivatives_data(jj) == JQT_line_gamma_selfexponent)
      {
        dF[jj] *= dgamma_dline_gamma_selfexponent;
      }
      else if(derivatives_data(jj) == JQT_line_gamma_foreignexponent)
      {
        dF[jj] *= dgamma_dline_gamma_foreignexponent;
      }
      else if(derivatives_data(jj) == JQT_line_gamma_waterexponent)
      {
        dF[jj] *= dgamma_dline_gamma_waterexponent;
      }
    }
  }
  
private:
  LineShapeType mtype;
  LineShapeNorm mnorm;
  Numeric cutoff;
};

typedef Array<LineFunctions> ArrayOfLineFunctions;

#endif //lineshapedata_h



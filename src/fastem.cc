/* Copyright (C) 2002,2003 
   Sreerekha Ravi<rekha@sat.physik.uni-bremen.de>
   Stefan Buehler <sbuehler@uni-bremen.de>                  

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. 
*/
  
/*!
  \file   fastem.cc
  \author Sreerekha Ravi <rekha@sat.physik.uni-bremen.de>
  \date   Tue Aug 10 15:16:31 2004
  
  \brief  This file contains functions that are adapted from FASTEM 
  code which is used to calculate surface emissivity.
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <cmath>
#include "matpackI.h"
#include "exceptions.h"
#include "complex.h"

extern const Numeric PI;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
//! Calculate the surface emissivity using Fastem model
/*! 
  Calculate surface emissivity using fastem  

  \param surface_emiss surface emissivity at one one point
  \param surface_temp surface temperature at one point
  \param surface_wind surface wind at one point
  \param f_grid frequency grid
  \param f_index index for the monochromatic frequency
  \param stokes_dim stokes dimension

  \author Sreerekha Ravi
  \date 2004-08-10
*/
void fastem(// Output:
            VectorView surface_emiss,
            // Input:
            const Numeric& surface_temp,
            ConstVectorView surface_wind,
            ConstVectorView surface_fastem_constants,
            const Numeric& freq
           )
{
  //  Calculate PIOM (Ellison et al.) xperm
  //Calculate xperm using the dclamkaouchi method
  //  Calculate PIOM (Ellison et al.) xperm

  // Convert the surface temperature in Kelvin to centigrade.
  Numeric temp_c  = surface_temp - 273.15;
  Numeric temp_cc = temp_c * temp_c;
  Numeric temp_ccc = temp_cc * temp_c;

  if (  (temp_c < -5.0)  ||  (temp_c > 100.0)  || 
        (freq < 10e+9) || (freq > 500e+9) )
    {
      
      ostringstream os;
      os << "Severe warning from dclamkaouchi: "
        << "The accepted temperature range in centigrade is "
        << "[-5,100],\nbut a value of " << temp_c 
        << "°C was found. Also the allowed frequency range is "
        << "[10 GHz,500 GHz],\nbut a value of " <<  freq
        << " was found.";
      
      throw runtime_error( os.str() );
    }
  else 
  if (	(freq < 20e+9) || (freq > 200e+9) )
    {
      
      ostringstream os;
      os << "Warning from dclamkaouchi: "
        << "The accepted temperature range in centigrade is "
        << "[-5,100],\nbut a value of " << temp_c 
        << "°C was found. Also the allowed frequency range is "
        << "[10 GHz,500 GHz],\nbut a value of " <<  freq
        << " was found."<< surface_wind; //remove surface_wind, it
                                         //was only to avoid
                                         //
                                         //compilation error due to unused variable  

      
      throw runtime_error( os.str() );
    }

  // define the two relaxation frequencies, tau1 and tau2
  Numeric tau1 = surface_fastem_constants[0] + surface_fastem_constants[1]* temp_c +  
    surface_fastem_constants[2] * temp_cc;

  Numeric tau2 = surface_fastem_constants[3] + surface_fastem_constants[4]* temp_c +  
    surface_fastem_constants[5] * temp_cc + surface_fastem_constants[6] * temp_ccc; 

  // define static xperm - FIXME TRS
  Numeric del1 = surface_fastem_constants[7] + surface_fastem_constants[8]* temp_c + 
    surface_fastem_constants[9] * temp_cc + surface_fastem_constants[10] * temp_ccc;

  Numeric del2 = surface_fastem_constants[11] + surface_fastem_constants[12]* temp_c + 
    surface_fastem_constants[13] * temp_cc + surface_fastem_constants[14] * temp_ccc;

  Numeric einf = surface_fastem_constants[17] + surface_fastem_constants[18] * temp_c;

  //calculate xperm using double debye formula
  Numeric fen = 2.0 * surface_fastem_constants[19] * freq/1e+9 * 0.001;
  Numeric fen2 = pow(fen,2);
  Numeric den1 = 1.0 + fen2 * tau1 * tau1;
  Numeric den2 = 1.0 + fen2 * tau2 * tau2;
  Numeric perm_real1 = del1/den1;
  Numeric perm_real2 = del2/den2;
  Numeric perm_imag1 = del1 * fen * tau1/den1;
  Numeric perm_imag2 = del2 * fen * tau2/den2;
  Numeric perm_real = perm_real1 + perm_real2 + einf;
  Numeric perm_imag = perm_imag1 + perm_imag2;
  
  Complex xperm (perm_real, perm_imag); //FIXME use complex here

  //Now the fresnel calculations
  // This is used to calculate vertical and horizontal polarised 
  //reflectivities given xperm at local incidence angle. I am not sure
  // how to include this theta now!!! FIXME
  Numeric theta = 55.0;
 
  Numeric cos_theta = cos(theta * DEG2RAD);
  Numeric sin_theta = sin(theta * DEG2RAD);
  
  //Numeric cos_2 = pow(cos_theta, 2);
  Numeric sin_2 = pow(sin_theta, 2);

  Complex perm1 = sqrt(xperm - sin_2);
  Complex perm2 = xperm * cos_theta;

  Complex rhth = (cos_theta - perm1)/(cos_theta + perm1);
  Complex rvth = (perm2 - perm1)/(perm2 + perm1);

  //Numeric rvertsr = real.rvth();
  //Numeric rvertsi = imag.rvth();
  Numeric rvertsr = real(rvth);
  Numeric rvertsi = imag(rvth);

  Numeric rverts = pow(rvertsr, 2) + pow(rvertsi, 2);

  //Numeric rhorzsr = real.rhth();
  //Numeric rhorzsi = imag.rhth();
  Numeric rhorzsr = real(rhth);
  Numeric rhorzsi = imag(rhth);
  Numeric rhorzs = pow(rhorzsr, 2) + pow(rhorzsi, 2);

  surface_emiss[0] = rverts;
  surface_emiss[1] = rhorzs;
  surface_emiss[2] = 0;
  surface_emiss[3] = 0;
  
    
}



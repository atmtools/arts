/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>

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
   USA. */




/*===========================================================================
  === File description 
  ===========================================================================*/

/*!
   \file   physics_funcs.cc
   \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   \date   2002-05-08 

   This file contains the code of functions of physical character.
   Modified by Claudia Emde (2002-05-28).
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "physics_funcs.h"
#include "messages.h"          
#include "mystring.h"
#include "physics_funcs.h"

extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric BOLTZMAN_CONST;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! invplanck
/*!
   Convert radiance to Plack brightness temperature.

    \return y Output:spectrum vector       
    \param  f       frequency
    \param  za      zenith angle

    \author Patrick Eriksson 
    \date   2000-09-28 
*/
Numeric invplanck (
		const Numeric&  f,
		const Numeric&  za )
{
  Numeric y;

  // Use always double to avoid numerical problem (see invrayjean)
  const double   a = PLANCK_CONST/BOLTZMAN_CONST;
  const double   b = 2*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
        double   c,d;

  // Check input
  if ( y > 1e-4 )  
    throw runtime_error("The spectrum cannot be in expected intensity unit "
                        "(impossible value detected).");
  
  c = a*f;
  d = b*pow(f,3);
  y = c / ( log(d/y+1) );
  return y;
}



//! invrayjean
/*! 
   Converts radiance to Rayleigh-Jean brightness temperature.

    \return y Output:spectrum vector       
    \param  f       frequency
    \param  za      zenith angle

    \author Patrick Eriksson 
    \date   2000-09-28 
*/
Numeric invrayjean (
		 const Numeric&  f,
		 const Numeric&  za )
{
  Numeric y;
  
 // The function returned NaNs when a and b were set to be Numeric (PE 010404)
  const double   a = SPEED_OF_LIGHT*SPEED_OF_LIGHT/(2*BOLTZMAN_CONST);
        double   b;

  // Check input
  if (y > 1e-4 )  
    throw runtime_error("The spectrum is not in expected intensity unit "
                                               "(impossible value detected).");
	
  b = a/(f*f);
  y = b * y;
  return y;
}


//! number_density
/*! 
   Calculates the number density.
   
   \return  nd Output: number density
   \param  p  Input: pressure
   \param  t  Input: temperature
   
   \author Patrick Eriksson 
   \date   2000-04-08 
*/
Numeric number_density (  
		     const Numeric& p,
		     const Numeric& t )
{
  Numeric nd;
  // Calculate p / (t*BOLTZMAN_CONST):
  assert( t > 0 );
  nd = p/(t*BOLTZMAN_CONST);
  return nd;
}



//! planck
/*! 
  Calculates the Planck function for a single temperature.
  
  \return B Output: the blackbody radiation
  \param  f Input: frequency value
  \param  t Input: temperature value
  
  \author Patrick Eriksson 
  \date   2000-04-08 
*/
Numeric planck (
	       const Numeric& f,
	       const Numeric& t )
{
  Numeric B;

  // Double must be used here (if not, a becomes 0 when using float)
  static const double  a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const double  b = PLANCK_CONST/BOLTZMAN_CONST;
  const double  c = b/t; 
  
  B = a * f*f*f / ( exp( f*c ) - 1.0 );
  return B;
}





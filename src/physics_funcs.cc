/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>

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
   \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
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

extern const Numeric BOLTZMAN_CONST;
extern const Numeric DEG2RAD;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! dinvplanckdI
/*!
   Calculates the derivative of inverse-Planck with respect to intensity.

    \return     The derivative
    \param  i   radiance
    \param  f   frequency

    \author Patrick Eriksson 
    \date   2010-10-26
*/
Numeric dinvplanckdI(
        const Numeric&  i,
        const Numeric&  f )
{
  assert( i >= 0 );
  assert( f > 0 );

  static const Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  static const Numeric b = 2 * PLANCK_CONST / ( SPEED_OF_LIGHT*SPEED_OF_LIGHT );

  const Numeric d    = b * pow(f,3.0) / i;
  const Numeric binv = a * f / log( d + 1 );

  return pow(binv,2.0) / ( a * f * i * (1/d+1) ); 
}



//! fresnel
/*!
    Calculates complex AMPLITUDE reflection coeffcients for a specular
    reflection

    The properties of the two involved media are given as the complex
    refractive index, n. A dielectric constant, eps, is converted as 
    n = sqrt( eps ). The power reflection coefficient, r, for one 
    polarisation is r = abs(R)^2.

    \param  Rv    Out: Reflection coefficient for vertical polarisation
    \param  Rh    Out: Reflection coefficient for vertical polarisation
    \param  n1    In: Refractive index of medium where radiation propagates
    \param  n2    In: Refractive index of reflecting medium 
    \param  theta In: Propagation angle from normal of radiation to be 
                      reflected

    \author Patrick Eriksson 
    \date   2004-09-21
*/
void fresnel(
             Complex&   Rv,
             Complex&   Rh,
       const Complex&   n1,
       const Complex&   n2,
       const Numeric&   theta )
{
  const Numeric theta1 = DEG2RAD * theta;
  const Numeric costheta1 = cos( theta1 );
  const Numeric costheta2 = cos( asin( n1.real() * sin(theta1) / n2.real() ) );

  Complex a, b;
  a  = n2 * costheta1;
  b  = n1 * costheta2;
  Rv = ( a - b ) / ( a + b );
  a  = n1 * costheta1;
  b  = n2 * costheta2;
  Rh = ( a - b ) / ( a + b );
}


//! invplanck
/*!
   Converts a radiance to Plack brightness temperature.

    \return     Planck brightness temperature
    \param  i   radiance
    \param  f   frequency

    \author Patrick Eriksson 
    \date   2002-08-11
*/
Numeric invplanck(
        const Numeric&  i,
        const Numeric&  f )
{
  assert( i >= 0 );
  assert( f > 0 );

  static const Numeric a = PLANCK_CONST / BOLTZMAN_CONST;
  static const Numeric b = 2 * PLANCK_CONST / ( SPEED_OF_LIGHT*SPEED_OF_LIGHT );

  return   ( a * f ) / log( (b*pow(f,3.0))/i + 1.0 );
}



//! invrayjean
/*! 
   Converts a radiance to Rayleigh-Jean brightness temperature.

    \return     RJ brightness temperature
    \param  i   radiance
    \param  f   frequency

    \author Patrick Eriksson 
    \date   2000-09-28 
*/
Numeric invrayjean(
        const Numeric&  i,
        const Numeric&  f )
{
  assert( f > 0 );

  static const Numeric   a = SPEED_OF_LIGHT*SPEED_OF_LIGHT/(2*BOLTZMAN_CONST);

  return  ( a * i ) / ( f * f );
}


//! number_density
/*! 
   Calculates the atmospheric number density.
   
   \return     number density
   \param  p   pressure
   \param  t   temperature
   
   \author Patrick Eriksson 
   \date   2000-04-08 
*/
Numeric number_density(  
        const Numeric&   p,
        const Numeric&   t )
{
  assert( p >= 0 );
  assert( t > 0 );
  return   p / ( t * BOLTZMAN_CONST );
}



//! planck
/*! 
  Calculates the Planck function for a single temperature.

  Note that this expression gives the intensity for both polarisations.
  
  \return     blackbody radiation
  \param  f   frequency
  \param  t   temperature
  
  \author Patrick Eriksson 
  \date   2000-04-08 
*/
Numeric planck( 
        const Numeric&   f, 
        const Numeric&   t )
{
  assert( t >= 0 );
  assert( f > 0 );

  if( t == 0 )
    { return 0; }

  else
    {
    static const Numeric a = 2 * PLANCK_CONST / (SPEED_OF_LIGHT*SPEED_OF_LIGHT);
    static const Numeric b = PLANCK_CONST / BOLTZMAN_CONST;
  
    return  ( a * pow(f,3.0) ) / ( exp( (b*f)/t ) - 1.0 );
    }
}



//! rayjean
/*! 
    Converts a Rayleigh-Jean brightness temperature to radiance

    \return     radiance
    \param  tb  RJ brightness temperature
    \param  f   frequency

    \author Patrick Eriksson 
    \date   2011-07-13 
*/
Numeric rayjean(
        const Numeric&  f,
        const Numeric&  t )
{
  assert( f > 0 );

  static const Numeric   a = SPEED_OF_LIGHT*SPEED_OF_LIGHT/(2*BOLTZMAN_CONST);

  return  ( f * f ) / ( a * t );
}

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




/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
   \file   physics_funcs.cc
   \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   \date   2002-05-08 

   This file contains the code of functions of physical character.
*/



/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include <math.h>
#include <stdexcept>
#include "physics_funcs.h"
#include "messages.h"          
#include "mystring.h"
#include "physics_funcs.h"

extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric BOLTZMAN_CONST;



/*****************************************************************************
 *** The functions (in alphabetical order)
 *****************************************************************************/

//! invplanck
/*!
   Converts a vector with radiances to Plack brightness temperatures.

    \param y Output:       spectrum vector       
    \param  f       frequencies
    \param  za      zenith angles

    \author Patrick Eriksson 
    \date   2000-09-28 
*/
void invplanck (
                   VectorView   y,
              ConstVectorView   f,
              ConstVectorView   za )
{
  const Index   nf  = f.nelem();
  const Index   nza = za.nelem();
  const Index   ny  = y.nelem();
        Index   i0;

  // Use always double to avoid numerical problem (see invrayjean)
  const double   a = PLANCK_CONST/BOLTZMAN_CONST;
  const double   b = 2*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
        double   c,d;

  // Check input
  if ( max(y) > 1e-4 )  
    throw runtime_error("The spectrum cannot be in expected intensity unit "
                        "(impossible value detected).");
  //
  if ( nf*nza != ny )  
  {
    ostringstream os;
    os << "The length of *y* does not match *f_mono* and *za_pencil*.\n"
       << "y.nelem():         " << y.nelem() << "\n"
       << "Should be f_mono.nelem()*za_pencil.nelem(): "
       << f.nelem() * za.nelem() << "\n"
       << "f_mono.nelem():  " << f.nelem() << "\n"
       << "za_pencil.nelem(): " << za.nelem();
    throw runtime_error(os.str());
  }

  for ( Index i=0; i<nf; i++ )
  {
    c = a*f[i];
    d = b*f[i]*f[i]*f[i];
    for ( Index j=0; j<nza; j++ )    
    {
      i0 = j*nf + i;
      y[i0] = c / ( log(d/y[i0]+1) );
    }
  }
}



//! invrayjean
/*! 
   Converts a vector with radiances to Rayleigh-Jean brightness temperatures.

    \param y Output:       spectrum vector       
    \param  f       frequencies
    \param  za      zenith angles

    \author Patrick Eriksson 
    \date   2000-09-28 
*/
void invrayjean (
                   VectorView   y,
              ConstVectorView   f,
              ConstVectorView   za )
{
  const Index   nf  = f.nelem();
  const Index   nza = za.nelem();
  const Index   ny  = y.nelem();
        Index   i0;

  // The function returned NaNs when a and b were set to be Numeric (PE 010404)
  const double   a = SPEED_OF_LIGHT*SPEED_OF_LIGHT/(2*BOLTZMAN_CONST);
        double   b;

  // Check input
  if ( max(y) > 1e-4 )  
    throw runtime_error("The spectrum is not in expected intensity unit "
                        "(impossible value detected).");
  //
  if ( nf*nza != ny )  
  {
    ostringstream os;
    os << "The length of *y* does not match *f_mono* and *za_pencil*.\n"
       << "y.nelem():         " << y.nelem() << "\n"
       << "Should be f_mono.nelem()*za_pencil.nelem(): "
       << f.nelem() * za.nelem() << "\n"
       << "f_mono.nelem():  " << f.nelem() << "\n"
       << "za_pencil.nelem(): " << za.nelem();
    throw runtime_error(os.str()); 
  }

  for ( Index i=0; i<nf; i++ )
  {
    b = a/(f[i]*f[i]);
    for ( Index j=0; j<nza; j++ )    
    {
      i0 = j*nf + i;
      y[i0] = b * y[i0];
    }
  }
}



//! number_density
/*! 
   Calculates the number density (scalar version).

    \return         number density
    \param  p       pressure
    \param  t       temperature

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
Numeric number_density (
			Numeric   p,
			Numeric   t )
{
  assert( 0!=t );
  return  p/t/BOLTZMAN_CONST;
}



//! number_density
/*! 
   Calculates the number density (vector version).

    \return number density
    \param  p       pressure
    \param  t       temperature

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
Vector number_density (
		       ConstVectorView    p,
		       ConstVectorView    t )
{
  assert( p.nelem()==t.nelem() );

  // Calculate p / (t*BOLTZMAN_CONST):

  Vector dummy(p);		// Matpack can initialize a
				// new Vector from another
				// Vector.
  dummy /= t;			// Element-vise divide by t.
  dummy /= BOLTZMAN_CONST;	// Divide all elements by BOLTZMAN_CONST.

  return dummy; 
}



//! planck 
/*! 
    Calculates a blackbody radiation (the Planck function) matrix.

    Each row of the returned matrix corresponds to a frequency, while each
    column corresponds to a temperature.

    \param B Output: the blackbody radiation
    \param  f       a frequency grid
    \param  t       a temperature profile

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void planck (
	     MatrixView      B, 
	     ConstVectorView f,
	     ConstVectorView t )
{
  // Double must be used here (if not, a becomes 0 when using float)
  static const double  a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const double  b = PLANCK_CONST/BOLTZMAN_CONST;

  const Index    n_f  = f.nelem();
  const Index    n_t  = t.nelem();
  Index    i_f, i_t;
  Numeric   c, d;

  assert( n_f==B.nrows() );
  assert( n_t==B.ncols() );

  for ( i_f=0; i_f<n_f; i_f++ )
  {
    c = a * f[i_f]*f[i_f]*f[i_f];
    d = b * f[i_f];
    for ( i_t=0; i_t<n_t; i_t++ )
      B(i_f,i_t) = c / (exp(d/t[i_t]) - 1.0);
  }
}



//! planck
/*! 
    Calculates the Planck function for a single temperature.

    \param B Output: the blackbody radiation
    \param  f       a frequency grid
    \param  t       a temperature value

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void planck (
             VectorView    B,
	     ConstVectorView    f,
	     Numeric   t )
{
  // Double must be used here (if not, a becomes 0 when using float)
  static const double  a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const double  b = PLANCK_CONST/BOLTZMAN_CONST;
         const double  c = b/t; 

  assert( B.nelem()==f.nelem() );

  for ( Index i=0; i<f.nelem(); i++ )
  {
    B[i] = a * f[i]*f[i]*f[i] / ( exp( f[i]*c ) - 1.0 );
  }
}




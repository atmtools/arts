/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>
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




////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file   atm_funcs.cc

   This file contains the code of functions releated to atmospheric 
   physics or geometry.

   \author Patrick Eriksson
   \date 2000-09-18 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdexcept>
#include "arts.h"
#include "matpackI.h"
#include "messages.h"          
#include "math_funcs.h"          
#include "make_vector.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric PLANCK_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric BOLTZMAN_CONST;



////////////////////////////////////////////////////////////////////////////
//   Physical functions
////////////////////////////////////////////////////////////////////////////

//// planck (matrix version) ///////////////////////////////////////////////
//
/** Calculates a blackbody radiation (the Planck function) matrix.

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



//// planck (vector version) ////////////////////////////////////////////////
//
/** Calculates the Planck function for a single temperature.

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



//// invplanck ////////////////////////////////////////////////////////////////
//
/** Converts a vector with radiances to Plack brightness temperatures.

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



//// invrayjean ///////////////////////////////////////////////////////////////
//
/** Converts a vector with radiances to Rayleigh-Jean brightness temperatures.

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



//// number_density (scalar version) ////////////////////////////////////////
//
/** Calculates the number density (scalar version).

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



//// number_density (vector version) ////////////////////////////////////////
//
/** Calculates the number density (vector version).

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



//// g_of_z ////////////////////////////////////////
//
/** Calculates the gravitational accelaration for a geometrical altitude.

    \return           the gravitational constant
    \param  r_geoid   radius of the geoid
    \param  g0        the gravitational constant at the geoid surface
    \param  z         geometrical altitude

    \author Patrick Eriksson 
    \date   2000-12-04
*/
Numeric g_of_z (
		Numeric   r_geoid,
		Numeric   g0,
		Numeric   z )
{
  return g0 * pow( r_geoid/(r_geoid+z), 2 );
}



////////////////////////////////////////////////////////////////////////////
//   Conversion and interpolation of pressure and altitude grids.
////////////////////////////////////////////////////////////////////////////

//// z2p ///////////////////////////////////////////////////////////////////
//
/** Converts an altitude vector to pressures.

    The log of the pressures are interpolated linearly.
    In Matlab notation:

      p = exp(interp1(z0,log(p0),z,'linear'))

    \param p Output: the pressures at z
    \param  z0      original altitude grid
    \param  p0      original pressure grid
    \param  z       new altitude grid

    \author Patrick Eriksson 
    \date   2000-04-08
*/
void z2p(
	 VectorView      p,
	 ConstVectorView z0,
	 ConstVectorView p0,
	 ConstVectorView z )
{
  assert( p.nelem()==z.nelem() );
  if ( z.nelem() > 0 )
  {
    // Local variable to store log pressure grid:
    Vector logp0(p0.nelem());
    transform( logp0, log, p0 );	// This calculates logp0 = log(p0).

    interp_lin_vector( p, z0, logp0, z );
    transform( p, exp, p );	        // This calculates p = exp(p).
  }
}



//// interpp (vector version) ///////////////////////////////////////////////
//
/** Interpolates a vertical profile at a new set of pressures.

    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      p = interp1(log(p0),x,log(p),'linear')

    \param x Output: the interpolated values at p
    \param  p0      original pressure grid
    \param  x0      the profile to be interpolated
    \param  p       new pressure grid

    \author Patrick Eriksson 
    \date   2000-04-08
*/
void interpp(
	     VectorView          x, 
	     ConstVectorView     p0,
	     ConstVectorView     x0,
	     ConstVectorView     p )
{
  assert( x.nelem()==p.nelem() );

  // Local variables to store log pressure grids:
  Vector logp0(p0.nelem());
  Vector logp(p.nelem());
  transform( logp0, log, p0 );	// This calculates logp0 = log(p0).
  transform( logp,  log, p  );	// This calculates logp  = log(p).

  interp_lin_vector( x, logp0, x0, logp );
}

void interpp_cloud(
		   VectorView     x, 
		   ConstVectorView     p0,
		   ConstVectorView     x0,
		   ConstVectorView     p )
{
  assert( x.nelem()==p.nelem() );

  interp_lin_vector( x, p0, x0, p );
}


//// interpp (matrix version) ///////////////////////////////////////////////
//
/** Interpolates a matrix, such as an absorption matrix, at a new 
    set of pressures.

    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      A = interp1(log(p0),A0,log(p),'linear')

    \param A Output: the interpolated values at p
    \param  p0      original pressure grid
    \param  A0      the matrix to be interpolated
    \param  p       new pressure grid

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void interpp(
	     MatrixView  A,
	     ConstVectorView  p0, 
	     ConstMatrixView  A0, 
	     ConstVectorView  p )
{
  assert( A.nrows()  == A0.nrows() );
  assert( A.ncols()  == p.nelem()  ); 

  // Local variables to store log pressure grids:
  Vector logp0(p0.nelem());
  Vector logp(p.nelem());
  transform( logp0, log, p0 );	// This calculates logp0 = log(p0).
  transform( logp,  log, p );	// This calculates logp0 = log(p0).

  interp_lin_matrix( A, logp0, A0, logp );
}



//// interpp (scalar return version) ////////////////////////////////////////
//
/** Interpolates a vertical profile at one pressure level.

    See the vector version.

    \param x Output: the interpolated values at p
    \param  p0      original pressure grid
    \param  x0      the profile to be interpolated
    \param  p       a pressure level

    \author Patrick Eriksson 
    \date   2000-12-04
*/
Numeric interpp(
		ConstVectorView     p0,
		ConstVectorView     x0,
		const Numeric       p )
{
  // Local variable to store log pressure grid:
  Vector logp0(p0.nelem());
  transform( logp0, log, p0 );	// This calculates logp0 = log(p0).

  return interp_lin( logp0, x0, log(p) );
}



//// interpz (vector version) ///////////////////////////////////////////////
//
/** Interpolates a vertical profile at a new set of vertical altitudes.

    NOTE!! Avoid to use this function, interpolation should mainly be done
    in pressure, that is, use interpp when possible.

    This function uses z2p and interpp to make an interpolation for vertical 
    altitudes. 

    Used mainly for LOS calculations with refraction.

    \param x Output: the interpolated values at z
    \param  p0      original pressure grid
    \param  z0      original vertical altitude grid
    \param  x0      the profile to be interpolated
    \param  z       new vertical altitude grid

    \author Patrick Eriksson 
    \date   2000-10-02 
*/
void interpz(
	     VectorView     x, 
	     ConstVectorView     p0,
	     ConstVectorView     z0,
	     ConstVectorView     x0,
	     ConstVectorView     z )
{
  assert( x.nelem()==z.nelem() ); 
  Vector p(z.nelem());
  z2p( p, z0, p0, z );
  interpp( x, p0, x0, p );
}



//// interpz (scalar version) ///////////////////////////////////////////////
//
/** Interpolates a vertical profile at a single vertical altitude.

    NOTE!! Avoid to use this function, interpolation should mainly be done
    in pressure, that is, use interpp when possible.

    This function uses z2p and interpp to make an interpolation for a vertical 
    altitude. 

    Used mainly for LOS calculations with refraction.

    \param x Output: the interpolated values at z
    \param  p0      original pressure grid
    \param  z0      original vertical altitude grid
    \param  x0      the profile to be interpolated
    \param  z       new vertical altitude grid

    \author Patrick Eriksson 
    \date   2000-10-02 
*/
Numeric interpz(
        ConstVectorView     p0,
        ConstVectorView     z0,
        ConstVectorView     x0,
        const Numeric    z )
{
  Vector x(1);
  MakeVector Z(z);
  interpz( x, p0, z0, x0, Z );
  return x[0];
}




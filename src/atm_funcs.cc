/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

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

#include "arts.h"
#include "vecmat.h"
#include "messages.h"          
#include "math_funcs.h"          
extern const Numeric EARTH_RADIUS;
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

    \retval B       output: the blackbody radiation
    \param  f       a frequency grid
    \param  t       a temperature profile

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void planck (
              MATRIX&     B, 
        const VECTOR&     f,
        const VECTOR&     t )
{
  static const Numeric a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST/BOLTZMAN_CONST;
  const size_t  n_f  = f.dim();
  const size_t  n_t  = t.dim();
  size_t        i_f, i_t;
  Numeric       c, d;
  B.newsize(n_f,n_t);
  
  for ( i_f=1; i_f<=n_f; i_f++ )
  {
    c = a * f(i_f)*f(i_f)*f(i_f);
    d = b * f(i_f);
    for ( i_t=1; i_t<=n_t; i_t++ )
      B(i_f,i_t) = c / (exp(d/t(i_t)) - 1.0);
  }
}



//// planck (vector version) ////////////////////////////////////////////////
//
/** Calculates the Planck function for a single temperature.

    \retval B       output: the blackbody radiation
    \param  f       a frequency grid
    \param  t       a temperature value

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void planck (
             VECTOR&    B,
       const VECTOR&    f,
       const Numeric&   t )
{
  static const Numeric a = 2.0*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
  static const Numeric b = PLANCK_CONST/BOLTZMAN_CONST;
  
  B = ediv( a*emult(f,emult(f,f)), exp(f*(b/t))-1.0 ) ;
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
       const Numeric&   p,
       const Numeric&   t )
{
  return  p/t/BOLTZMAN_CONST;
}



//// number_density (vector version) ////////////////////////////////////////
//
/** Calculates the number density (vector version).

    \return nd      number density
    \param  p       pressure
    \param  t       temperature

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
VECTOR number_density (
       const VECTOR&    p,
       const VECTOR&    t )
{
  return ediv(p,t)/BOLTZMAN_CONST;
}



////////////////////////////////////////////////////////////////////////////
//   Core functions for RTE and BL 
////////////////////////////////////////////////////////////////////////////

//// rte_iterate
//
/** Performs a single iteration for RTE calculations (one zenith angle).

    The vector Y is not initilised, the obtained values are added to Y.
    Note that only a single iteration is performed.

    This function can be used to calculate emission spectra for parts of
    the atmosphere.
        
    \retval y             the spectrum
    \param  start_index   start index for the integration
    \param  stop_index    stop index for the integration
    \param  Tr            transmission matrix
    \param  S             source function matrix
    \param  n_f           number of frequencies

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void rte_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   tr,
       const MATRIX&   s,
       const size_t    n_f )
{
        size_t   i_f;        // frequency index
           int   i_z;        // LOS index
           int   i_step;     // step order, -1 or 1

  if ( start_index >= stop_index )
    i_step = -1;
  else
    i_step = 1;

  for ( i_z=start_index; i_z!=(stop_index+i_step); i_z+=i_step ) 
  {
    for ( i_f=1; i_f<=n_f; i_f++ )    
      y(i_f) = y(i_f)*tr(i_f,i_z) + s(i_f,i_z) * ( 1.0-tr(i_f,i_z) );
  }
}



//// rte ////////////////////////////////////////////////////////////////////
//
/** Performs the RTE calculations for one zenith angle.

    This function allows calculation of emission spectra for single
    zenith angles in function beside yRteXx.
        
    \retval y             the spectrum
    \param  start_index   start index for the integration
    \param  stop_index    stop index for the integration
    \param  Tr            transmission matrix
    \param  S             source function matrix
    \param  y_space       intensity entering the atmosphre at start of LOS
    \param  ground        flag/index for ground intersection
    \param  e_ground      ground emissivity
    \param  y_ground      ground blackbody radiation 

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void rte (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   tr,
       const MATRIX&   s,
       const VECTOR&   y_space,
       const int&      ground,
       const VECTOR&   e_ground,
       const VECTOR&   y_ground )
{
  const int   n_f = tr.dim(1);               // number of frequencies
        int   i_f;                           // frequency index
        int   i_break;                       // break index for looping
        int   i_start;                       // variable for second loop

  // Init Y with Y_SPACE
  y = y_space;

  // Check if LOS inside the atmosphere (if START_INDEX=0, Y=Y_SPACE)
  if ( start_index > 0 )
  {
    // Determine break index for looping, either 1 or the ground
    if ( ground > 0 )
      i_break = ground;
    else
      i_break = 1;       

    // Make first loop
    rte_iterate( y, start_index-1, i_break, tr, s, n_f );

    // We are now at the sensor, the ground or the tangent point
    // We are ready only if we are at the sensor.
    // If at sensor, we have that STOP_INDEX=1 and not GROUND
    if ( !(stop_index==1 && !ground) )
    {
      // Set most common values for I_START and I_BREAK
      i_start = 1;
      i_break = stop_index - 1;
      
      // If at the ground, include ground reflection. 
      // The loop can continue both downwards or upwards
      if ( ground )
      {            
        for ( i_f=1; i_f<=n_f; i_f++ )    
          y(i_f) = y(i_f)*(1.0-e_ground(i_f)) + y_ground(i_f)*e_ground(i_f);
        
        if ( ground != 1 )  // 2D case, loop downwards
	{
         i_start = ground - 1;
         i_break = 1;
        }
      }

      // Make second loop
      rte_iterate( y, i_start, i_break, tr, s, n_f );

    } // second part
  } // if any values
}



//// bl_iterate /////////////////////////////////////////////////////////////
//
/** Performs a single iteration for BL calculations (one zenith angle).

    The vector Y is not initilised, Y is multiplied with the obtained values.
    Note that only a single iteration is performed.

    This function can be used to calculate transmissions for parts of
    the atmosphere.
        
    \retval y             the spectrum
    \param  start_index   start index for the integration
    \param  stop_index    stop index for the integration
    \param  Tr            transmission matrix
    \param  S             source function matrix
    \param  n_f           number of frequencies

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void bl_iterate (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   tr,
       const size_t    n_f )
{
        size_t   i_f;        // frequency index
           int   i_z;        // LOS index
           int   i_step;     // step order, -1 or 1

  if ( start_index >= stop_index )
    i_step = -1;
  else
    i_step = 1;

  for ( i_z=start_index; i_z!=(stop_index+i_step); i_z+=i_step ) 
  {
    for ( i_f=1; i_f<=n_f; i_f++ )    
      y(i_f) *= tr(i_f,i_z);
  }
}



//// bl //////////////////////////////////////////////////////////////////////
//
/** Performs the BL (transmission) calculations for one zenith angle.

    This function allows calculation of transmission spectra for single
    zenith angles in functions beside yBlXx.
        
    \retval y             the spectrum
    \param  start_index   start index for the integration
    \param  stop_index    stop index for the integration
    \param  Tr            transmission matrix
    \param  ground        flag/index for ground intersection
    \param  e_ground      ground emissivity

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void bl (
             VECTOR&   y,
       const int&      start_index,
       const int&      stop_index,
       const MATRIX&   tr,
       const int&      ground,
       const VECTOR&   e_ground )
{
  const int   nf = tr.dim(1);          // number of frequencies
        int   iy;                      // frequency index

  // Init Y
  y.newsize(nf);
  y = 1.0;

  // Loop steps passed twice
  if ( stop_index > 1 )
  {
    bl_iterate( y, 1, stop_index-1, tr, nf );
    y = emult( y, y );  
  }

  // Loop remaining steps
  if ( start_index != stop_index )
    bl_iterate( y, stop_index, start_index-1, tr, nf );

  // Include effect of ground reflection
  if ( ground )
  {
    for ( iy=1; iy<=nf; iy++ )    
      y(iy) *= ( 1.0 - e_ground(iy) );
  }
}



////////////////////////////////////////////////////////////////////////////
//   Conversion and interpolation of pressure and altitude grids.
////////////////////////////////////////////////////////////////////////////

//// z2p ///////////////////////////////////////////////////////////////////
//
/** Converts an altitude vector to pressures.

    The log. of the pressures are interpolated linearly.
    In Matlab notation:

      p = exp(interp1(z0,log(p0),z,'linear'))

    \retval p       output: the pressures at z
    \param  z0      original altitude grid
    \param  p0      original pressure grid
    \param  z       new altitude grid

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void z2p(
              VECTOR&     p,
        const VECTOR&     z0,
        const VECTOR&     p0,
        const VECTOR&     z )
{
  if ( z.dim() > 0 )
    p = exp( interp_lin( z0, log(p0), z ) );
}



//// interpp (vector version) ///////////////////////////////////////////////
//
/** Interpolates a vertical profile at a new set of pressures.

    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      p = interp1(log(p0),x,log(p),'linear')

    \retval x       output: the interpolated values at p
    \param  p0      original pressure grid
    \param  x0      the profile to be interpolated
    \param  p       new pressure grid

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void interpp(
              VECTOR&     x, 
        const VECTOR&     p0,
        const VECTOR&     x0,
        const VECTOR&     p )
{
  interp_lin( x, log(p0), x0, log(p) );
}



//// interpp (matrix version) ///////////////////////////////////////////////
//
/** Interpolates a matrix, such as an absorption matrix, at a new 
    set of pressures.

    A linear interpolation using log. pressure is applied.
    In Matlab notation, the following expression is used:

      A = interp1(log(p0),A0,log(p),'linear')

    \retval A       output: the interpolated values at p
    \param  p0      original pressure grid
    \param  A0      the matrix to be interpolated
    \param  p       new pressure grid

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
void interpp(
              MATRIX&  A,
        const VECTOR&  p0, 
        const MATRIX&  A0, 
        const VECTOR&  p )
{
  interp_lin_row( A, log(p0), A0, log(p) );
}



//// interpz (vector version) ///////////////////////////////////////////////
//
/** Interpolates a vertical profile at a new set of vertical altitudes.

    NOTE!! Avoid to use this function, interpolation should mainly be done
    in pressure, that is, use interpp when possible.

    This function uses z2p and interpp to make an interpolation for vertical 
    altitudes. 

    Used mainly for LOS calculations with refraction.

    \retval x       output: the interpolated values at z
    \param  p0      original pressure grid
    \param  z0      original vertical altitude grid
    \param  x0      the profile to be interpolated
    \param  z       new vertical altitude grid

    \author Patrick Eriksson 
    \date   2000-10-02 
*/
void interpz(
              VECTOR&     x, 
        const VECTOR&     p0,
        const VECTOR&     z0,
        const VECTOR&     x0,
        const VECTOR&     z )
{
  VECTOR p;
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

    \retval x       output: the interpolated values at z
    \param  p0      original pressure grid
    \param  z0      original vertical altitude grid
    \param  x0      the profile to be interpolated
    \param  z       new vertical altitude grid

    \author Patrick Eriksson 
    \date   2000-10-02 
*/
Numeric interpz(
        const VECTOR&     p0,
        const VECTOR&     z0,
        const VECTOR&     x0,
        const Numeric&    z )
{
  VECTOR x;
  interpz( x, p0, z0, x0, VECTOR(1,z) );
  return x(1);
}



/////////////////////////////////////////////////////////////////////////////
//   Tangent altitudes.
/////////////////////////////////////////////////////////////////////////////

//// ztan_geom //////////////////////////////////////////////////////////////
//
/** Calculates the geometrical tangent altitude (no refraction).

    \return        the tangent altitude
    \param za      the angle between zenith and the LOS
    \param z_plat  the platform altitude

    \author Patrick Eriksson 
    \date   2000-04-08 
*/
Numeric ztan_geom(
        const Numeric&     za,
        const Numeric&     z_plat )
{
  Numeric  z_tan;
  if ( za >= 90 )   
    z_tan = (EARTH_RADIUS+z_plat)*sin(DEG2RAD*za) - EARTH_RADIUS; 
  else
    z_tan = 9.9999e6;
  return z_tan;
}



//// ztan_refr //////////////////////////////////////////////////////////////
//
/** Calculates the tangent altitude with refraction.

    \return               the tangent altitude
    \param    c           LOS constant
    \param    za          the angle between zenith and the LOS
    \param    z_plat      the platform altitude
    \param    z_ground    the ground altitude
    \param    p_abs       absorption pressure grid
    \param    z_abs       absorption altitude grid
    \param    refr_index  refrective index corresponding to p_abs

    \author Patrick Eriksson 
    \date   2000-10-02
*/
Numeric ztan_refr(
        const Numeric&     c,
        const Numeric&     za,
        const Numeric&     z_plat,
        const Numeric&     z_ground,
        const VECTOR&      p_abs,
        const VECTOR&      z_abs,
        const VECTOR&      refr_index )
{
  if ( za < 90 )   //=== Upward ==========================================
    return ztan_geom( za, z_plat);
  else
  {
    const size_t  n = z_abs.dim();
          size_t  i;

    for ( i=n; (i>0) && (EARTH_RADIUS+z_abs(i))*refr_index(i)>c; i-- ) 
    {
      if ( z_abs(i) <= z_ground ) //=== Ground intersection ==============
      {
        Numeric n_ground = interpz( p_abs, z_abs, refr_index, z_ground );
        Numeric theta = RAD2DEG*asin(c/n_ground/(EARTH_RADIUS+z_ground));
        return ztan_geom( 180-theta, z_ground );
      }
    }
    if ( i == n )      //=== outside the atmosphere ======================
      return ztan_geom( za, z_plat);
    else               //=== z_tan inside the atmosphere =================
    {
      VECTOR zs(2), cs(2);
      zs(1) = z_abs(i);
      zs(2) = z_abs(i+1);
      cs(1) = (EARTH_RADIUS+z_abs(i))*refr_index(i);
      cs(2) = (EARTH_RADIUS+z_abs(i+1))*refr_index(i+1);  
      return interp_lin( cs, zs, c ); 
    }
  }
}




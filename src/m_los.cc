/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>

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
   \file   m_los.cc

   This file contains functions associated with 1D line-of-sight (LOS)
   calculations.

   Types of functions are:

     1. determination of LOS
     2. calculation of transmissions and source function along the LOS
     3. solving the radiative transfer equation
   
   Functions in this file assumes LTE and no scattering.
   The LOS is defined by a structure of type LOS, defined in los.h.

   \author Patrick Eriksson
   \date 2000-09-14 
*/



////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include "arts.h"
#include "atm_funcs.h"          
#include "matpackI.h"
#include "los.h"
#include "math_funcs.h"          
#include "messages.h"          
#include "auto_wsv.h"          
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric COSMIC_BG_TEMP;
extern const Numeric SUN_TEMP;
extern const Numeric PLANCK_CONST;
extern const Numeric BOLTZMAN_CONST;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric EARTH_GRAV_CONST;



////////////////////////////////////////////////////////////////////////////
//   LOS help functions 
////////////////////////////////////////////////////////////////////////////

//// any_ground /////////////////////////////////////////////////////////////
/**
   Checks if any of the pencil beam spectra corresponds to a ground reflection.

   This function is most likely called as any_ground( los.ground )

   \return           1 if any ground intersection, 0 otherwise
   \param    ground  array of ground index/flag values

   \author Patrick Eriksson
   \date   2000-12-12
*/
bool any_ground( const ArrayOfIndex& ground )  
{
  for ( Index i=0; i<ground.nelem(); i++ )
  {
    if ( ground[i] )
      return 1;
  }  
  return 0;
}




//// los_geometric /////////////////////////////////////////////////////////
/**
   Core function for geometric LOS calculations.

   All observations geometries are handled. The variables za and z_plat
   shall be treated as the "zenith angle" and the vertical altitude of
   the lowest point of the LOS.

   The variable l_step should be set to before calling this function.
   The value of this variable is normally untouched, but when the distance
   between the starting point and the atmsopheric limit is smaller than
   the given value, l_step is adjusted.

   \retval   z           vertical altitudes of the LOS points
   \retval   psi         the angle in the observation plane between the vectors
                         going from the sensor and the LOS point, respectively,
                         to the centre of the earth geoid.angle from the sensor
                         of the LOS points
   \retval   l_step      distance along the LOS between the points
   \param    z_plat      platform altitude
   \param    za          zentith angles
   \param    atm_limit   maximum altitude of the absorption grid
   \param    r_geoid     the local geoid radius

   \author Patrick Eriksson
   \date   2001-02-15
*/
void los_geometric(
		    Vector&     z,
                    Vector&     psi,
                    Numeric&    l_step,
              const Numeric&    z_plat,
              const Numeric&    za,
              const Numeric&    atm_limit,
	      const Numeric&    r_geoid )
{
  // A safety check
  assert( za <= 90 );

  Vector    l;      // Lengths along the LOS from the tangent point
  Index     nz;     // Length of z and psi

  // Some temporary values are always double to avoid numerical problems
  // (this is especially a problem where r_geoid is squared).
  double    a, b;   // Temporary values
  double    llim;   // distance to atmospheric limit

  // Distance from the lowest point of the LOS to the atmospheric limit
  a     = r_geoid + atm_limit;
  b     = (r_geoid+z_plat) * sin(DEG2RAD*za);
  llim  = sqrt( a*a - b*b ) - (r_geoid+z_plat)*cos(DEG2RAD*za) ;

  // Handle the rare case that llim < l_step
  if ( llim < l_step )         
    l_step = llim*0.9999;       // *0.9999 to avoid problem in interpolations

  // Create equally spaced points along the LOS
  linspace( l, 0, llim, l_step );

  nz = l.nelem();
  z.resize(   nz );
  psi.resize( nz );

  // Calculate vertical altitudes and angles
  b = r_geoid + z_plat;  
  a = b * b;
  //
  for ( Index i=0; i<nz; i++ )
  {
    z[i]   = sqrt( a + l[i]*l[i] + 2.0*b*l[i]*cos(DEG2RAD*za) );
    psi[i] = RAD2DEG * acos( (a+z[i]*z[i]-l[i]*l[i]) / (2.0*b*z[i]) ); 

    // Nan can in some cases be obtained for very small angles 
    if ( isnan(psi[i]) )
      psi[i] = 0;

    z[i]     = z[i] - r_geoid;
  }
}



//// los_refraction /////////////////////////////////////////////////////////
/**
   Core function for LOS calculations with refraction.

   See further los_geometric.

   \retval   z           vertical altitudes of the LOS points
   \retval   psi         the angle in the observation plane between the LOS
                         points and the lowest point of the LOS
   \retval   l_step      distance along the LOS between the points
   \param    z_plat      platform altitude
   \param    za          zentith angles
   \param    atm_limit   maximum altitude of the absorption grid
   \param    r_geoid     the local geoid radius
   \param    p_abs       absorption pressure grid
   \param    z_abs       absorption altitude grid
   \param    refr        refraction flag
   \param    refr_lfac   ray tracing length factor
   \param    refr_index  refrective index corresponding to p_refr
   \param    c           constant for the LOS

   \author Patrick Eriksson
   \date   2001-02-18
*/
void los_refraction(
		    Vector&     z,
                    Vector&     psi,
                    Numeric&    l_step,
              const Numeric&    z_plat,
              const Numeric&    za,
              const Numeric&    atm_limit,
	      const Numeric&    r_geoid,
	      const Vector&     p_abs,
	      const Vector&     z_abs,
              const Index&      refr,
              const Index&      refr_lfac,
              const Vector&     refr_index,
              const Numeric&    c )
{
  // A safety check
  assert( za <= 90 );

  // Allocate memory for temporary z and psi. To be safe, make vectors 50 %
  // as long than for the geometric case
  Index np;
  {
    // Distance from the lowest point of the LOS to the atmospheric limit
    double  a    = r_geoid + atm_limit;
    double  b    = (r_geoid+z_plat) * sin(DEG2RAD*za);
    double  llim = sqrt( a*a - b*b ) - (r_geoid+z_plat)*cos(DEG2RAD*za) ;

    // Handle the rare case that llim < l_step
    if ( llim < l_step )         
      l_step = llim*0.9999;       // *0.9999 to avoid problem in interpolations

    np = Index( ceil( 1.5 * ( llim/l_step + 1) ) );
  }
  Vector   zv(np), pv(np); 

  // Double is used here instead of Numeric to avoid nuerical problems
  const double l = l_step / refr_lfac;   // Step length of ray tracing
  Index    i = refr_lfac;                // Ray tracing step counter
  double   z1;                           // Old altitude of the LOS
  double   z2 = z_plat;                  // New altitude of the LOS
  double   rz1, rz2;                     // As z1 and z2 but + r_geoid
  double   psi1;                         // Old angle of the LOS
  double   psi2 = 0;                     // New angle of the LOS
  double   n1, n2;                       // Refractive index at z1 and z2
  double   n;                            // Either n1 or the mean of n1 and n2
  double   c2=c; c2 = c2 * c2;           // Square of the LOS constant
  Index    j;                            // See below
  double   d, e, f;                      // Some temporary values

  np = 0;

  // To save computational time, the interpolation is handled locally so
  // the indeces for the refr_index vector can be remembered from one 
  // interpolation to next.
  const Index   nz = z_abs.nelem();
        Index   iz; 
  for ( iz=0; (iz<nz) && (z_abs[iz]<=z2); iz++ ) {}
  if ( iz < nz )
    n2 = refr_index[iz-1] + (refr_index[iz]-refr_index[iz-1])*
                                      (z2-z_abs[iz-1])/(z_abs[iz]-z_abs[iz-1]);
  else
    n2 = 1;

  while ( z2 <= atm_limit )
  {

    z1   = z2;
    psi1 = psi2;
    n1   = n2;

    if ( i == Index(refr_lfac) )
    {    
      zv[np] = z2;
      pv[np] = RAD2DEG * psi2;
      i     = 1;
      np++;
      assert( np < zv.nelem() );
    }
    else
      i++; 

    // We repeat the calculation of z2 some times to get better estimates for
    // a mean of value of n between z1 and z2. For first iteration n is set 
    // to n1, and for later iteration as the mean of n1 and n2.
    //
    // A practical test showed a clear improvement when doing 2 iterations 
    // instead of a single iteration, but just a marginally improvement when
    // going to 3 iterations. So 2 iterations seem to be the best choice.
    //
    for ( j=1; j<=2; j++ )
    {

      if ( j == 1 )
        n = n1;
      else
        n = ( n1 + n2 ) / 2;

      rz1 = z1 + r_geoid;
      d   = rz1 * rz1;
      e   = c2/(n*n);
      f   = d - e;
 
      // When using float, there have been NaNs here (due to z1 < c/n).
      // So we must make a check to avoid these NaNs.
      // 
      if ( f <= 0 )
        rz2 = sqrt( l*l + e );
      else
        rz2 = sqrt( pow( l + sqrt(f), 2 ) + e );

      z2 = rz2 - r_geoid;

      // Determine n at z2
      for ( ; (iz<nz) && (z_abs[iz]<=z2); iz++ ) {}
      if ( iz < nz )
        n2 = refr_index[iz-1] + (refr_index[iz]-refr_index[iz-1])*
                                      (z2-z_abs[iz-1])/(z_abs[iz]-z_abs[iz-1]);
      else
        n2 = 1;
    }

    psi2 = psi1 + acos( (d+rz2*rz2-l*l) / (2*rz1*rz2) ); 
  }

  // Move values from temporary vectors
  z.resize(   np );
  psi.resize( np );
  z   = zv[Range(0,np)];	// Copy first np values of zv to z.
  psi = pv[Range(0,np)];	// Copy first np values of pv to psi.
}



////////////////////////////////////////////////////////////////////////////
//   The sub-function to losCalc
////////////////////////////////////////////////////////////////////////////

//// los_1za ///////////////////////////////////////////////////////////////
/**
   Performs the LOS calculations for one zenith angle.

   All observations geometries are handled.   

   \retval   z           vertical altitudes of the LOS points
   \retval   psi         the angle in the observation plane between the vectors
                         going from the sensor and the LOS point, respectively,
                         to the centre of the earth geoid.angle from the sensor
                         of the LOS points
   \retval   l_step      distance along the LOS between the points
   \retval   ground      ground flag (0 = no ground intersection)
   \retval   start       start index when solving the RTE
   \retval   stop        stop index when solving the RTE
   \param    z_plat      platform altitude
   \param    za          zentith angles
   \param    l_step_max  the user defined maximum step length along the LOS
   \param    atm_limit   maximum altitude of the absorption grid
   \param    z_ground    altitude of the ground (above the geoid)
   \param    r_geoid     the local geoid radius
   \param    p_abs       absorption pressure grid
   \param    z_abs       absorption altitude grid
   \param    refr        refraction flag
   \param    refr_lfac   ray tracing length factor
   \param    refr_index  refrective index corresponding to p_refr

   \author Patrick Eriksson
   \date   2001-02-18
*/
void los_1za(
		    Vector&     z,
                    Vector&     psi,
                    Numeric&    l_step,
                    Index&      ground,
                    Index&      start,
                    Index&      stop,
                    Numeric&    z_tan,
              const Numeric&    z_plat,
              const Numeric&    za,
              const Numeric&    l_step_max,
              const Numeric&    z_ground,
	      const Numeric&    r_geoid,
	      const Vector&     p_abs,
	      const Vector&     z_abs,
              const Index&        refr,
              const Index&        refr_lfac,
              const Vector&     refr_index )
{
  Numeric   c;        // LOS constant when considering refraction

  // Determine the upper limit of the atmosphere
  const Numeric atm_limit = last(z_abs);

  if ( refr )
  {
    c = refr_constant( r_geoid, za, z_plat, p_abs, z_abs, atm_limit, 
                                                                  refr_index );
    z_tan = ztan_refr( c, za, z_plat, z_ground, p_abs, z_abs, refr_index, 
                                                                     r_geoid );
  }
  else
    z_tan  = ztan_geom( za, z_plat, r_geoid );

  // Set l_step to its most probable value
  l_step = l_step_max;


  //=== Observation from space ================================================
  if ( z_plat >= atm_limit )
  {
    Index     nz;          // Length of z and psi
    Numeric   psi0 = 0;    // Correction value for psi

    out3 << " (z_tan = " << z_tan/1e3 << " km)";
    
    // If LOS outside the atmosphere, return empty vectors
    if ( z_tan >= atm_limit )
    {
      ground = 0;
      z.resize(   0 );
      psi.resize( 0 );
      nz     = 1;
    }
  
    // Only through the atmosphere
    else if ( z_tan >= z_ground )
    {
      if ( !refr )
      {
        los_geometric( z, psi, l_step, z_tan, 90.0, atm_limit, r_geoid );
        psi0 = za - 90.0;
      }
      else
      {
        los_refraction( z, psi, l_step, z_tan, 90.0, atm_limit, r_geoid, 
                                p_abs, z_abs, refr, refr_lfac, refr_index, c );

        // Determine the "zenith angle" of the LOS at the top of the atmosphere
        Numeric zmax = last( z );
        Numeric n = interpz( p_abs, z_abs, refr_index, zmax );
        Numeric theta = RAD2DEG * asin( c / ((r_geoid+zmax)*n)  );

        psi0 = theta + za - 180.0 + last(psi);
      }

      ground = 0;
      nz     = z.nelem();
    }   
  
    // Intersection with the ground
    else
    {
      // The "zenith angle" at ground level
      Numeric za_g;

      if ( !refr )
      {
        za_g = RAD2DEG * asin( (r_geoid+z_tan) / (r_geoid+z_ground) );   

        los_geometric( z, psi, l_step, z_ground, za_g, atm_limit, r_geoid );

        psi0 = za + za_g - 180.0;
      }
      else
      {
        za_g = RAD2DEG * asin( c / ( (r_geoid+z_ground) * 
                  n_for_z( z_ground, p_abs, z_abs, refr_index, atm_limit ) ) );

        los_refraction( z, psi, l_step, z_ground, za_g, atm_limit, r_geoid, 
                                p_abs, z_abs, refr, refr_lfac, refr_index, c );

        // Determine the "zenith angle" of the LOS at the top of the atmosphere
        Numeric zmax = last( z );
        Numeric n = n_for_z( zmax, p_abs, z_abs, refr_index, atm_limit );
        Numeric theta = RAD2DEG * asin( c / ((r_geoid+zmax)*n)  );

        psi0 = theta + za - 180.0 + last(psi);
      }

      ground = 1;
      nz     = z.nelem();
    }

    // Add psi0 to all elements:
    if ( psi0 != 0 )
      psi += psi0;			// Matpack can add element-vise like this.
  
    start = stop = nz - 1;
  }


  //=== Inside the atmosphere looking upwards =================================
  else if ( za <= 90 )
  {
    if ( !refr )
      los_geometric( z, psi, l_step, z_plat, za, atm_limit, r_geoid );
    else
      los_refraction( z, psi, l_step, z_plat, za, atm_limit, r_geoid, 
                                p_abs, z_abs, refr, refr_lfac, refr_index, c );
    ground = 0;
    stop   = 0;
    start  = z.nelem() - 1;
  }

  //=== Inside the atmosphere looking downwards ===============================
  else
  {
    // Some temporary values are always double to avoid numerical problems
    // (this is especially a problem where r_geoid is squared).
    double   l1;     // Distance between platform and tangent point or ground
    double   a, b;   // Temporary values

    out3 << " (z_tan = " << z_tan/1e3 << " km)";

    // Only through the atmosphere
    if ( z_tan >= z_ground )
    {
      if ( !refr )
      {
	// Calculate the distance platform-tangent point
	a  = r_geoid + z_plat;
	b  = r_geoid + z_tan; 
	l1 = sqrt(a*a-b*b);   

        // A sufficient large distance between platform and tangent point
        if ( l1 > l_step_max/10 )
	{
  	  // Adjust l_step downwards to get an integer number of steps
	  stop  = Index( ceil( l1 / l_step_max + 1.0 ) - 1 );  
	  l_step = l1 / Numeric(stop);
	}

        // Ignore downwrd part if platform is very close to tangent point
        else
	{
          l_step = l_step_max;
          stop   = 0;
	}
	  
        los_geometric( z, psi, l_step, z_tan, 90.0, atm_limit, r_geoid );
      }
      else
      {
        // Calculate a first LOS from the tangent point and up to the sensor
        // using l_step/refr_lfac as step length
        Numeric   l = l_step / refr_lfac;
        los_refraction( z, psi, l, z_tan, 90.0, z_plat+l, r_geoid, 
                                        p_abs, z_abs, refr, 1, refr_index, c );

        // Determine the distance along the LOS between the tangent point and
        // the sensor by an interpolation
	{
	  Vector ls;
	  linspace( ls, 0, l*(z.nelem()-1) , l );
	  l1 = interp_lin( z, ls, z_plat );
	}

        // A sufficient large distance between platform and tangent point
        if ( l1 > l_step_max/10 )
	{
  	  // Adjust l_step downwards to get an integer number of steps
	  stop  = Index( ceil( l1 / l_step_max + 1.0 ) - 1 );  
	  l_step = l1 / Numeric(stop);
	}
  
        // Ignore downwrd part if platform is very close to tangent point
        else
	{
          l_step = l_step_max;
          stop   = 0;
	}

        los_refraction( z, psi, l_step, z_tan, 90.0, atm_limit, r_geoid, 
                                p_abs, z_abs, refr, refr_lfac, refr_index, c );
      }

      ground = 0;
      start  = z.nelem() - 1;

      // The angular distance between the sensor and the tangent point 
      // is psi[stop]
      psi += psi[stop];		// Matpack can add element-vise like this.
    }

    // Intersection with the ground
    else
    {
      Numeric za_g;       // The "zenith angle" at ground level

      if ( !refr )
      {
	// Calculate the distance platform-ground
	a  = r_geoid + z_plat;
	b  = r_geoid + z_tan;
	b  = b * b; 
	l1 = sqrt(a*a-b);   
	a  = r_geoid + z_ground;
	l1 = l1 - sqrt(a*a-b);   
  
        // A sufficient large distance between platform and ground
        if ( l1 > l_step_max/10 )
	{
	  // Adjust l_step downwards to get an integer number of steps
	  stop  = Index( ceil( l1 / l_step_max + 1.0 ) - 1 );  
	  l_step = l1 / Numeric(stop);
	}
  
        // Ignore downwrd part if platform is very close to ground
        else
	{
          l_step = l_step_max;
          stop   = 0;
	}

	za_g = RAD2DEG * asin( (r_geoid+z_tan) / (r_geoid+z_ground) );
  
	los_geometric( z, psi, l_step, z_ground, za_g, atm_limit, r_geoid );
      }

      else
      {
        za_g = RAD2DEG * asin( c / ( (r_geoid+z_ground) * 
                  n_for_z( z_ground, p_abs, z_abs, refr_index, atm_limit ) ) );

        Numeric   l = l_step / refr_lfac;
        los_refraction( z, psi, l, z_ground, za_g, z_plat+l, r_geoid, 
                                        p_abs, z_abs, refr, 1, refr_index, c );

        // Determine the distance along the LOS between the tangent point and
        // the sensor by an interpolation
	{
	  Vector ls;
	  linspace( ls, 0, l*(z.nelem()-1) , l );
	  l1 = interp_lin( z, ls, z_plat );
	}

        // A sufficient large distance between platform and ground
        if ( l1 > l_step_max/10 )
	{
	  // Adjust l_step downwards to get an integer number of steps
	  stop  = Index( ceil( l1 / l_step_max + 1.0 ) - 1 );  
	  l_step = l1 / Numeric(stop);
	}

        // Ignore downwrd part if platform is very close to ground
        else
	{
          l_step = l_step_max;
          stop   = 0;
	}

        los_refraction( z, psi, l_step, z_ground, za_g, atm_limit, r_geoid, 
                                p_abs, z_abs, refr, refr_lfac, refr_index, c );
      }

      ground = 1;
      start  = z.nelem() - 1;

      // The angular distance between the sensor and the ground
      // is psi[stop]
      psi += psi[stop];		// Matpack can add element-vise like this.
    }
  }
}



////////////////////////////////////////////////////////////////////////////
//   The sub-function to yCalc
////////////////////////////////////////////////////////////////////////////

/**
   \author Patrick Eriksson
   \date   2000-??-??
*/
void y_rte (
                    Vector&          y,
              const Los&             los,   
              const Vector&          f_mono,
              const Vector&          y_space,
              const ArrayOfMatrix&   source,
              const ArrayOfMatrix&   trans,
              const Vector&          e_ground,
              const Numeric&         t_ground )
{
  // Some variables
  const Index   n=los.start.nelem();  // Number of zenith angles 
  const Index   nf=f_mono.nelem();    // Number of frequencies 
        Vector   y_tmp(nf);           // Temporary storage for spectra
        Index   iy0=0;               // Reference index for output vector

  out2 << "  Integrating the radiative transfer eq. with emission.\n";

  // Resize y
  y.resize( nf*n );
        
  // Set up vector for ground blackbody radiation if any ground intersection
  // Check also if the ground emission vector has the correct length
  Vector   y_ground(f_mono.nelem()); 
  if ( any_ground(los.ground) )  
  {
    if ( t_ground <= 0 )
      throw runtime_error(
          "There are intersections with the ground, but the ground\n"
          "temperature is set to be <=0 (are dummy values used?).");
    if ( e_ground.nelem() != nf )
      throw runtime_error(
          "There are intersections with the ground, but the frequency and\n"
          "ground emission vectors have different lengths (are dummy values\n"
          "used?).");
    out2 << "  There are intersections with the ground.\n";
    planck( y_ground, f_mono, t_ground );
  }

  // Loop zenith angles
  out3 << "    Zenith angle nr:      ";
  for ( Index i=0; i<n; i++ )
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Iteration is done in seperate function    
    rte( y_tmp, los.start[i], los.stop[i], trans[i], 
                 source[i], y_space, los.ground[i], e_ground, y_ground);

    // Move values to output vector
    y[Range(iy0,nf)] = y_tmp;	// Copy to nf elements to y, starting at iy0.

    // Update iy0
    iy0 += nf;   
  }
  out3 << "\n";
}



/**
   \author Patrick Eriksson
   \date   2000-??-??
*/
void y_tau (
                    Vector&          y,
              const Los&             los,   
              const ArrayOfMatrix&   trans,
              const Vector&          e_ground )
{
  // Some variables
  const Index   n=los.start.nelem();    // Number of zenith angles 
  const Index   nf=trans[0].nrows();    // Number of frequencies 
        Index   iy, iy0=0;              // Index for output vector
        Vector  y_tmp;                  // Temporary storage for spectra

  out2 << "  Calculating optical thicknesses.\n";

  // Resize y and set to 1
  y.resize( nf*n );
  y = 1.0;			// Matpack can set all elements like this.

  // Check if the ground emission vector has the correct length
  if ( any_ground(los.ground) )  
  {
    if ( e_ground.nelem() != nf )
      throw runtime_error(
          "There are intersections with the ground, but the frequency and\n"
          "ground emission vectors have different lengths (are dummy values\n"
          "used?).");
    out2 << "  There are intersections with the ground.\n";
  }
        
  // Loop zenith angles
  out3 << "    Zenith angle nr:     ";
  for ( Index i=0; i<n; i++ )
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Iteration is done in seperate function    
    bl( y_tmp, los.start[i], los.stop[i], trans[i], los.ground[i], e_ground );

    // Convert to optical thicknesses and move values to output vector
    for ( iy=0; iy<nf; iy++ )
      y[iy0+iy] = -log( y_tmp[iy] );

    // Update iy0
    iy0 += nf;           
  }
  out3 << "\n";
}




////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void r_geoidStd( Numeric&    r_geoid )
{
  extern const Numeric EARTH_RADIUS;
  r_geoid = EARTH_RADIUS;
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void r_geoidWGS84( 
              Numeric&   r_geoid,
        const Numeric&   latitude,
        const Numeric&   obsdirection )
{
  check_if_in_range( -90, 90, latitude, "latitude" );
  check_if_in_range( -360, 360, obsdirection, "obsdirection" );

  const Numeric rq = 6378.138e3, rp = 6356.752e3;
        Numeric a, b, rns, rew;

  // Calculate NS and EW radius
  a   = cos(latitude*DEG2RAD);
  b   = sin(latitude*DEG2RAD);
  rns = rq*rq*rp*rp/pow(rq*rq*a*a+rp*rp*b*b,1.5);
  rew = rq*rq/sqrt(rq*rq*a*a+rp*rp*b*b);

  // Calculate the radius in the observation direction
  a       = cos(obsdirection*DEG2RAD);
  b       = sin(obsdirection*DEG2RAD);
  r_geoid = 1/(a*a/rns+b*b/rew);
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void groundSet( 
              Numeric&   z_ground,
              Numeric&   t_ground,
              Vector&    e_ground,
        const Vector&    p_abs,
        const Vector&    t_abs,
        const Vector&    z_abs,
        const Vector&    f_mono,
	const Numeric&   z,
	const Numeric&   e )
{
  check_if_in_range( 0, 1, e, "e" );
  z_ground = z;
  t_ground = interpz( p_abs, z_abs, t_abs, z );
  e_ground.resize( f_mono.nelem() );
  e_ground = e;		
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void groundAtBottom( 
              Numeric&   z_ground,
              Numeric&   t_ground,
              Vector&    e_ground,
        const Vector&    t_abs,
        const Vector&    z_abs,
        const Vector&    f_mono,
	const Numeric&   e )
{
  check_if_in_range( 0, 1, e, "e" );
  z_ground = z_abs[0];
  t_ground = t_abs[0];
  e_ground.resize( f_mono.nelem() );
  e_ground = e;			// Matpack can set all elements like this.
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-01-22
*/
void groundOff( 
              Numeric&   z_ground,
              Numeric&   t_ground,
              Vector&    e_ground,
        const Vector&    z_abs )
{
  z_ground = z_abs[0];
  t_ground = 0;
  e_ground.resize( 0 );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void emissionOn( Index&   emission )
{
  emission = 1;
}


/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-04-19
*/
void emissionOff( Index&   emission )
{
  emission = 0;
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-02-15
*/
void losCalc(       Los&        los,
                    Vector&     z_tan,
              const Numeric&    z_plat,
              const Vector&     za,
              const Numeric&    l_step,
              const Vector&     p_abs,
              const Vector&     z_abs,
              const Index&      refr,
              const Index&      refr_lfac,
              const Vector&     refr_index,
              const Numeric&    z_ground,
              const Numeric&    r_geoid )
{     
  Index   n = za.nelem();  // number of zenith angles

  // Some checks                
  check_if_bool( refr, "refr" );                                      
  if ( z_ground < z_abs[0] )
    throw runtime_error(
      "There is a gap between the ground and the lowest absorption altitude.");
  if ( z_plat < z_ground )
    throw runtime_error("Your platform altitude is below the ground.");
  if ( z_plat < z_abs[0] )  
    throw runtime_error(
      "The platform cannot be below the lowest absorption altitude.");
  if ( refr && ( p_abs.nelem() != refr_index.nelem() ) )
    throw runtime_error(
      "Refraction is turned on, but the length of refr_index does not match\n"
      "the length of p_abs. Are dummy vales used?");
  if ( refr && ( refr_lfac < 1 ) )
    throw runtime_error(
      "Refraction is turned on, but the refraction length factor is < 1. \n"
      "Are dummy vales used?");
    
  // Reallocate the los structure and z_tan
  los.p.resize(      n	);
  los.psi.resize(    n	);
  los.z.resize(      n	);
  los.l_step.resize( n	);
  los.ground.resize( n	);
  los.start.resize(  n	);
  los.stop.resize(   n	);
  z_tan.resize(      n	);

  // Print messages
  if ( refr == 0 )
    out2 << "  Calculating line of sights WITHOUT refraction.\n";
  else if ( refr == 1 )
    out2 << "  Calculating line of sights WITH refraction.\n";
  //
  out3 << "     z_plat: " << z_plat/1e3 << " km\n";

  // Loop the zenith angles
  for ( Index i=0; i<n; i++ )
  {
    out3 << "         za: " << za[i] << " degs.";

    los_1za( los.z[i], los.psi[i], los.l_step[i], los.ground[i], los.start[i],
             los.stop[i], z_tan[i], z_plat, za[i], l_step, z_ground, r_geoid,
                                   p_abs, z_abs, refr, refr_lfac, refr_index );
    out3 << "\n";

    // Convert altitudes to pressures
    los.p[i].resize( los.z[i].nelem() );
    z2p( los.p[i], z_abs, p_abs, los.z[i] );
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Carlos Jimenez
   \date   2000-03-27
*/
void zaFromZtan(
        // WS Goutput
              Vector&       za,
        const String&       za_name,
	 // WS input
	const Vector&       z_tan,
        const Numeric&      z_plat,
        const Vector&       p_abs,
        const Vector&       z_abs,
	const Index&        refr,
	const Vector&       refr_index,
	const Numeric&      r_geoid,
        const Numeric&      z_ground )
{
  check_lengths( p_abs, "p_abs", z_abs, "z_abs" );  
  check_lengths( p_abs, "p_abs", refr_index, "refr_index" );  

  const Numeric atm_limit = last(z_abs);
  const Index          nz = z_tan.nelem();

  za.resize(nz);

  for (Index i=0; i<nz; i++)
  {
    if (za[i]>z_plat)
      throw runtime_error(
        "One tangent altitude is larger than the platform altitude");      

    // No refraction
    if (!refr)    
       za[i] = 90 + RAD2DEG*acos ( (r_geoid + z_tan[i]) / (r_geoid + z_plat) );
 
    // Refraction
    else
    {
      Numeric nz_plat =  n_for_z(z_plat,p_abs,z_abs,refr_index,atm_limit);  
      if (z_tan[i]>=0)
      { 
	// Calculating constant
        Numeric nza =  n_for_z(z_tan[i],p_abs,z_abs,refr_index,atm_limit);
        Numeric c   = (r_geoid + z_tan[i]) * nza;
        za[i]       =  180 - RAD2DEG * asin( c / nz_plat / (r_geoid + z_plat));
      }
      else
      {
	// inside the Earth, looking for hitting point
	Numeric ze  = RAD2DEG * acos((r_geoid + z_tan[i]) / r_geoid);
        // from hitting point to platform
        Numeric nze =  n_for_z(z_ground,p_abs,z_abs,refr_index,atm_limit);
        Numeric c   =  r_geoid * sin(DEG2RAD * (90-ze)) * nze;
        za[i]       =  180 - RAD2DEG * asin( c / nz_plat / (r_geoid + z_plat));
      } 
    }  
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Carlos Jimenez
   \date   2000-04-09
*/
void zaFromDeltat(
        // WS Generic Output:
        Vector&             za,
        // WS Generic Output Names:
        const String&       za_name,
        // WS Input:
        const Numeric&      z_plat,
        const Vector&       p_abs,
        const Vector&       z_abs,
        const Numeric&      l_step,
	const Index&        refr,
	const Index&        refr_lfac,
	const Vector&       refr_index,
	const Numeric&      r_geoid,
        const Numeric&      z_ground,
        // Control Parameters:
        const Numeric&      delta_t,
        const Vector&       z_tan_lim )

{
  // Checking stuff
  check_lengths( p_abs, "p_abs", z_abs, "z_abs" );  
  check_lengths( p_abs, "p_abs", refr_index, "refr_index" );  
  if ( z_tan_lim[0] > z_tan_lim[1] )
    throw runtime_error(
       "The lower tangent latitude is larger than the higher (in z_tan_lim).");

  // No refraction
  if (!refr)     
  {
    // Geometric calculations
    Vector phi(2);
    Vector za_lim(2);
    String zastr = "za_lim";

    zaFromZtan(za_lim, zastr, z_tan_lim, z_plat, p_abs, z_abs, refr, 
                                                refr_index, r_geoid, z_ground);

    phi[0] = za_lim[0] - 90;
    phi[1] = za_lim[1] - 90;

    const Numeric ang_step  = RAD2DEG * delta_t * 
                             sqrt (EARTH_GRAV_CONST / pow(r_geoid + z_plat,3));

    if (((phi[0]-phi[1])/ang_step) <=0)
      throw runtime_error("The time given is too large to have any cross-link in the given altitude range");     

    const Index n=Index(floor((phi[0]-phi[1])/ang_step));

    za.resize(n);
    for ( Index j=0;j<n;j++ )
      za[j] = 90 + phi[0] - (j * ang_step);
  }

  // Refraction
  else
  {  
    const Index ztanstep = 100;  // 100 meters step
    const Index n=Index(floor((z_tan_lim[1]-z_tan_lim[0])/ztanstep));
    
    Vector z_tan_1(n);
    Vector z_tan_2(n);
    za.resize(n);

    // ztan altitudes for later doing the interpolation
    for ( Index j=0;j<n;j++ )
      z_tan_1[j]=z_tan_lim[0]+j*ztanstep;
    
    // corresponding zenith angles
    String za_str = "za";
    zaFromZtan(za, za_str, z_tan_1, z_plat, p_abs, z_abs, refr, refr_index, 
                                                            r_geoid, z_ground);
    
    // corresponding psi
    Los los;
    losCalc(los,z_tan_2,z_plat,za,l_step,p_abs,z_abs,refr,refr_lfac,
                                                  refr_index,z_ground,r_geoid);
    // psi corresponding to the ztan defined for interpolation
    Vector psizb(n);
    for ( Index j=0;j<n;j++ )
      psizb[j]=los.psi[j][0];
    // lower and higer psi
    const Numeric psitop = psizb[n-1];
    const Numeric psibot = psizb[0];       

    // vel * deltat 
    const Numeric ang_step = RAD2DEG * delta_t * 
                             sqrt (EARTH_GRAV_CONST / pow(r_geoid + z_plat,3));

    // number of cross links in the ztan range specified for the given deltat
    if (((psibot-psitop)/ang_step)<=0)
      throw runtime_error("The time given is too large to have any cross-link in the given altitude range");  
    const Index np=Index(floor((psibot-psitop)/ang_step));
    // corresponding psi (in that z_tan_lim[1]-z_tan_lim[0])
    Vector psit(n);
    for ( Index j=0;j<np;j++ )
      psit[j]=psibot-j*ang_step;


    // ztan for psi for deltat from  psi and ztan for interpolation        
    z_tan_1.resize(np);
    for ( Index j=0;j<np;j++ )
    {
      z_tan_1[j]=interp_lin(psizb,z_tan_2,psit[j]);
    }
 
    // corresponding zenith angles
    zaFromZtan(za, za_str, z_tan_1, z_plat, p_abs, z_abs, refr, refr_index, 
                                                            r_geoid, z_ground);
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void sourceCalc(
                    ArrayOfMatrix&   source,
	      const Index&           emission,
              const Los&             los,   
              const Vector&          p_abs,
              const Vector&          t_abs,
              const Vector&          f_mono )
{
  check_if_bool( emission, "emission" );
  check_lengths( p_abs, "p_abs", t_abs, "t_abs" );  

  if ( emission == 0 )
  {
    out2 << "  Setting the source array to be empty.\n";
    source.resize( 0 );
  }

  else
  {     
	  Vector   tlos;                   // temperatures along the LOS
    const Index    nza=los.start.nelem();  // the number of zenith angles  
    const Index    nf=f_mono.nelem();      // the number of frequencies
	  Index    nlos;                   // the number of pressure points
	  Matrix   b;                      // the Planck function for TLOS  
	  Index    iv, ilos;               // frequency and LOS point index
  
    out2 << "  Calculating the source function for LTE and no scattering.\n";
   
    // Resize the source array
    source.resize(nza);
  
    // Loop the zenith angles and:
    //  1. interpolate the temperature
    //  2. calculate the Planck function for the interpolated temperatures
    //  3. take the mean of neighbouring Planck values
    out3 << "    Zenith angle nr:      ";
    for (Index i=0; i<nza; i++ ) 
    {
      if ( (i%20)==0 )
	out3 << "\n      ";
      out3 << " " << i; cout.flush();
  
      if ( los.p[i].nelem() > 0 )
      {
	nlos = los.p[i].nelem();
	tlos.resize( nlos );
	interpp( tlos, p_abs, t_abs, los.p[i] );
	b.resize( nf, nlos );
	planck( b, f_mono, tlos );
	source[i].resize(nf,nlos-1);
	for ( ilos=0; ilos<(nlos-1); ilos++ )
	{
	  for ( iv=0; iv<nf; iv++ )
	    source[i](iv,ilos) = ( b(iv,ilos) + b(iv,ilos+1) ) / 2.0;
	}
      }
    }  
    out3 << "\n";
  }
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void transCalc(
                    ArrayOfMatrix&   trans,
              const Los&             los,   
              const Vector&          p_abs,
              const Matrix&          abs )
{    
  check_length_ncol( p_abs, "p_abs", abs, "abs" );

  // Some variables
  const Index     n = los.start.nelem(); // the number of zenith angles
  const Index     nf = abs.nrows();      // the number of frequencies
        Index     np;                    // the number of pressure points
        Index     row, col;              // counters
        Matrix    abs2 ;                 // matrix to store interpolated abs.
        Numeric   w;                     // = -l_step/2

  out2 << "  Calculating transmissions WITHOUT scattering.\n";
 
  // Resize the transmission array
  trans.resize(n);

  // Loop the zenith angles and:
  //  1. interpolate the absorption
  //  2. calculate the transmission using the mean absorption between points
  out3 << "    Zenith angle nr:     ";
  for (Index i=0; i<n; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    np = los.p[i].nelem();
    if ( np > 0 )
    {
      abs2.resize( nf, np );
      interp_lin_matrix( abs2, p_abs, abs, los.p[i] );
      trans[i].resize( nf, np-1 );
      w  =  -0.5*los.l_step[i];
      for ( row=0; row<nf; row++ )
      {
        for ( col=0; col<(np-1); col++ )
          trans[i](row,col) = exp( w * ( abs2(row,col)+abs2(row,col+1) ) );
      }
    }
  }    
  out3 << "\n";
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-?-?
*/
void y_spaceStd(
                    Vector&   y_space,
              const Vector&   f,
              const String&   choice )
{
  y_space.resize( f.nelem() );

  if ( choice == "zero" )
  {
    y_space = 0.0;		// Matpack can set all elements like this.
    out2 << "  Setting y_space to zero.\n";
  }
  else if ( choice == "cbgr" )
  {
    planck( y_space, f, COSMIC_BG_TEMP );
    out2 << "  Setting y_space to cosmic background radiation.\n";
  }
  else if ( choice == "sun" )
  {
    planck( y_space, f, SUN_TEMP );
    out2 << "  Setting y_space to blackbody radiation corresponding to "
         << "the Sun temperature\n";
  }
  else
    throw runtime_error(
      "Possible choices for y_space are \"zero\", \"cbgr\" and \"sun\".");

}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-03-30
*/
void yCalc (
                    Vector&          y,
	      const Index&           emission,
              const Los&             los,   
              const Vector&          f_mono,
              const Vector&          y_space,
              const ArrayOfMatrix&   source,
              const ArrayOfMatrix&   trans,
              const Vector&          e_ground,
              const Numeric&         t_ground )
{
  // Check input
  // (ground emission and temperature are checked in the sub-functions)
  //
  check_if_bool( emission, "emission" );                                      
  check_length_nrow( f_mono, "f_mono", trans[0], 
                                     "the transmission matrices (in trans)" );
  if ( los.p.nelem() != trans.nelem() )
    throw runtime_error(
      "The number of zenith angles of *los* and *trans* are different.");
  if ( emission )
  {
    check_lengths( f_mono, "f_mono", y_space, "y_space" );  
    check_length_nrow( f_mono, "f_mono", source[0], 
                                          "the source matrices (in source)" );
    if ( los.p.nelem() != source.nelem() )
      throw runtime_error(
        "The number of zenith angles of *los* and *source* are different.");
  }

  if ( emission == 0 )
    y_tau( y, los, trans, e_ground );
  else
    y_rte( y, los, f_mono, y_space, source, trans, e_ground, t_ground );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2000-08-31
*/
void yTB (
                    Vector&          y,
              const Vector&          f_mono,
              const Vector&          za_pencil )
{
  out2 << "  Converts the values of *y* to Planck temperatures.\n";
  invplanck( y, f_mono, za_pencil );
}



/**
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2001-08-28
*/
void yTRJ (
                    Vector&          y,
              const Vector&          f_mono,
              const Vector&          za_pencil )
{
  out2 << "  Converts the values of *y* to Rayleigh-Jean temperatures.\n";
  invrayjean( y, f_mono, za_pencil );
}



/**
   Convert a matrix containing radiances to Rayleigh-Jeans BTs. 

   \author Stefan, Viju, Sreerekha
   \date   2001-05-02
*/
void MatrixTRJ (// WS Generic Output:
                Matrix& kout,
                // WS Generic Output Names:
                const String& kout_name,
                // WS Input:
                const Vector& f_mono,
                const Vector& za_pencil,
                // WS Generic Input:
                const Matrix& kin,
                // WS Generic Input Names:
                const String& kin_name)
{
  out2 << "  Converts the values of *" << kin_name << "* to Rayleigh-Jean "
       << "temperatures,\n  and stores the result in *" << kout_name 
       << "*.\n";

  // Resize kout if necessary:
  if ( kout.nrows()!=kin.nrows() || kout.ncols()!=kin.ncols() )
    kout.resize(kin.nrows(),kin.ncols());
  
  Vector y(kin.nrows());
  for ( Index i=0; i<kin.ncols(); i++ )
  {
    y = kin(Range(joker),i);	    // Copy ith column of kin to y.
    invrayjean( y, f_mono, za_pencil );
    kout(Range(joker),i) = y;	    // Copy to ith column of kout.
  }
}



/**
   Convert a matrix containing radiances to Planck BTs. 

   \author Stefan, Viju, Sreerekha
   \date   2001-05-02
*/
void MatrixTB (// WS Generic Output:
                Matrix& kout,
                // WS Generic Output Names:
                const String& kout_name,
                // WS Input:
                const Vector& f_mono,
                const Vector& za_pencil,
                // WS Generic Input:
                const Matrix& kin,
                // WS Generic Input Names:
                const String& kin_name)
{
  out2 << "  Converts the values of *" << kin_name << "* to Planck "
       << "temperatures,\n  and stores the result in *" << kout_name 
       << "*.\n";

  // Resize kout if necessary:
  if (kout.nrows()!=kin.nrows() ||
      kout.ncols()!=kin.ncols())
    kout.resize(kin.nrows(),kin.ncols());
  
  Vector y(kin.nrows());
  for ( Index i=0; i<kin.ncols(); i++ )
  {
    y = kin(Range(joker),i);	// Copy ith column of kin to y.
    invplanck( y, f_mono, za_pencil );
    kout(Range(joker),i) = y;	// Copy to ith column of kout.
  }
}



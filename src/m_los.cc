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
   \file   m_los.cc

   This file contains functions associated with 1D LOS calculations.

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

#include "arts.h"
#include "atm_funcs.h"          
#include "vecmat.h"
#include "los.h"
#include "math_funcs.h"          
#include "messages.h"          
#include "wsv.h"          
extern const Numeric EARTH_RADIUS;
extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;
extern const Numeric COSMIC_BG_TEMP;
extern const Numeric SUN_TEMP;
extern const Numeric PLANCK_CONST;
extern const Numeric BOLTZMAN_CONST;
extern const Numeric SPEED_OF_LIGHT;



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
bool any_ground( const ARRAY<int>& ground )  
{
  for ( size_t i=0; i<ground.size(); i++ )
  {
    if ( ground[i] >= 0 )
      return 1;
  }  
  return 0;
}


//// geom2refr /////////////////////////////////////////////////////////////
/**
   Determines a refractive LOS using the prolongation factor.

   This function takes a geometrical LOS that is parallel with the refractive
   LOS at the lowest point and, by using the prolongation factor, calculates
   the distance along the LOS for the altitudes of the geometrical LOS. These
   distances are interpolated to get the refractive LOS.

   See further AUG.

   \retval   p           pressures of the refractive LOS
   \retval   l_step      step length along the final LOS
   \param    zs          vertical altitudes of the geometrical LOS
   \param    c           constant for the refractive LOS
   \param    z_tan       tangent altitude
   \param    r_geoid     local geoid curvature
   \param    atm_limit   upper atmospheric limit
   \param    p_abs       absorption pressure grid
   \param    z_abs       absorption altitude grid
   \param    refr_index  refrective index corresponding to p_abs
   \param    l_step_refr step length along the input geometrical LOS

   \author Patrick Eriksson
   \date   2000-10-02

   Adapted to MTL. We have to set the size of p before calling
   interp_lin at the end of the function.  
   \date   2001-01-05
   \author Stefan Buehler
*/
void geom2refr( 
              VECTOR&       p,
              Numeric&      l_step,
        const VECTOR&       zs,
        const Numeric&      c,
        const Numeric&      z_tan,
        const Numeric&      r_geoid,
        const Numeric&      atm_limit,
        const VECTOR&       p_abs,
        const VECTOR&       z_abs,
        const VECTOR&       refr_index,
        const Numeric&      l_step_refr )
{
  const size_t nz = zs.size();
  VECTOR n(nz), a(nz), b(nz), r(nz), l(nz), ps(nz);

  // Get refractive index at the points of the geometrical LOS
  interpz( n, p_abs, z_abs, refr_index, zs );

  // Calculate the prolongation factor

  setto( a, r_geoid );
  add( zs, a );			//   a = (r_geoid+zs);
  ele_mult( a, a, a );		//   a = emult( a, a );
  // b is just a dummy I have introduced.
  setto( b, -c*c/n[0]/n[0] );
  add( a, b );
  transf( b, sqrt, b );
  ele_mult( n, b, r );		// r = emult( n, sqrt(a-(c*c/n(1)/n(1))) );

  // now we can use b as a dummy for something else. Note also that we
  // modify a here, since it is not used below.
  ele_mult(n,a,a);
  ele_mult(n,a,a);
  setto(b,-c*c);
  add(b,a);
  transf(a,sqrt,a);
  ele_div( r, a, r );	     // r = ediv( r, sqrt( emult(a,emult(n,n))-c*c ) );

  // Handle tangent points
  if ( abs(zs[0]-z_tan) < 0.01 )
    r[0] = sqrt(n[0]/(n[0]+(r_geoid+zs[0])*(n[1]-n[0])/(zs[1]-zs[0])));

  // Safety check
  if ( isnan(r[1]) != 0 )
    throw logic_error(
        "NaN obtained for the prolongation factor (check consistency for c!)");

  // Calculate distance from the lowest altitude along the LOS using r
  l[0] = 0.0;
  for ( size_t i=1; i<nz; i++ )
    l[i] = l[i-1] + l_step_refr * (r[i-1]+r[i])/2;
  
  // Interpolate to get LOS with steps of L_STEP
  if ( l[nz-1]<l_step )
    l_step = l[nz-1];
  z2p( ps, z_abs, p_abs, zs );
  VECTOR dummy = linspace(0,l[nz-1],l_step);
  resize( p, dummy.size() );
  //interp_lin( p, l, ps, dummy );
  interp_lin_vector( p, l, ps, dummy );
}



//// refr_constant ///////////////////////////////////////////////////////////
/**
   Determines the constant for a refractive LOS.

   Calculates (Re+z)*n(z)*sin(theta) at the platform.

   \return               LOS constant
   \param    r_geoid     local geoid curvature
   \param    za          zenith angle
   \param    z_plat      platform altitude
   \param    p_abs       absorption pressure grid
   \param    z_abs       absorption altitude grid
   \param    refr_index  refrective index corresponding to p_abs

   \author Patrick Eriksson
   \date   2000-10-02
*/
Numeric refr_constant( 
        const Numeric&      r_geoid,
        const Numeric&      za,
        const Numeric&      z_plat,
        const VECTOR&       p_abs,
        const VECTOR&       z_abs,
        const VECTOR&       refr_index )
{
  Numeric n_plat;

  if ( z_plat > last(z_abs) )
    n_plat = 1;
  else
    n_plat = interpz( p_abs, z_abs, refr_index, z_plat );

  return (r_geoid + z_plat) * sin(DEG2RAD*za) * n_plat;
}



//// upward_geom ///////////////////////////////////////////////////////////
/**
   A help function to calculate LOS for upward observations without refraction.

   This function handles only a single zenith angle.

   \retval   z           vertical altitudes for the LOS
   \retval   l_step      step length along the LOS
   \param    r_geoid     local geoid curvature
   \param    z_plat      platform altitude
   \param    za          zenith angle
   \param    atm_limit   maximum altitude of the absorption grid
   \param    fit_limit   flag to match perfectly atmospheric limit

   \author Patrick Eriksson
   \date   2000-09-14
*/
void upward_geom(
              VECTOR&       z,
              Numeric&      l_step,
        const Numeric&      r_geoid,
        const Numeric&      z_plat,
        const Numeric&      za,
        const Numeric&      atm_limit,
        const size_t&       fit_limit )
{
  Numeric  a, b;          // temporary values
  Numeric  llim;          // distance to atmospheric limit
  VECTOR   l;             // distance from sensor

  if ( za > 90 )
    throw logic_error("Upward function used for zenith angle > 90 deg."); 

  a     = r_geoid + atm_limit;
  b     = (r_geoid + z_plat)*sin(DEG2RAD*za);
  llim  = sqrt(a*a-b*b) - (r_geoid+z_plat)*cos(DEG2RAD*za) ;

  if ( !fit_limit ) 
  {
    if ( llim < l_step )    // Handle the rare case that llim < l_step
      l_step = llim*0.9999; // *0.9999 to avoid problem in interpolations
  }
  else
  {
    double n = ceil( llim / l_step + 1.0 );  
    l_step = 0.9999*llim/(n-1); // *0.9999 to avoid problem in interpolations
  }
  linspace( l, 0, llim, l_step );

  b = r_geoid + z_plat;  
  //  z = sqrt(b*b+emult(l,l)+(2.0*b*cos(DEG2RAD*za))*l) - r_geoid;
  VECTOR d1(l.size()), d2(l.size());  // just two dummies 
  ele_mult(l,l,d1);

  copy( scaled( l, 2.0*b*cos(DEG2RAD*za) ), d2 );

  resize( z, l.size() );
  setto( z, b*b );
  add( d1, z );
  add( d2, z );
  add( VECTOR( l.size(), -r_geoid ), z );
}



//// upward_refr ///////////////////////////////////////////////////////////
/**
   A help function to calculate LOS for upward observations with refraction.

   This function handles only a single zenith angle.

   \retval   p           pressures the LOS
   \retval   l_step      step length along the LOS
   \param    c           LOS constant
   \param    z_plat      platform altitude
   \param    za          zenith angle
   \param    r_geoid     local geoid curvature
   \param    atm_limit   maximum altitude of the absorption grid
   \param    p_abs       absorption pressure grid
   \param    z_abs       absorption altitude grid
   \param    refr_index  refrective index corresponding to p_abs
   \param    l_step_refr step length along the input geometrical LOS

   \author Patrick Eriksson
   \date   2000-10-02
*/
void upward_refr(
              VECTOR&       p,
              Numeric&      l_step,
        const Numeric&      c,
        const Numeric&      z_plat,
        const Numeric&      za,
        const Numeric&      r_geoid,
        const Numeric&      atm_limit,
        const VECTOR&       p_abs,
        const VECTOR&       z_abs,
        const VECTOR&       refr_index,
              Numeric       l_step_refr )
{
  VECTOR  zs;

  if ( za > 90 )
    throw logic_error("Upward function used for zenith angle > 90 deg."); 

  upward_geom( zs, l_step_refr, r_geoid, z_plat, za, atm_limit, 1 );
  geom2refr( p, l_step, zs, c, -99e3, r_geoid, atm_limit, p_abs, z_abs, 
                                                   refr_index, l_step_refr );
}



//// space_geom ////////////////////////////////////////////////////////////
/**
   Calculates the LOS for observations from space without refraction.

   This function handles only a single zenith angle.

   \retval   z           vertical altitudes for the LOS
   \retval   l_step      step length along the LOS
   \param    z_tan       tangent altitude
   \param    r_geoid     local geoid curvature
   \param    atm_limit   maximum altitude of the absorption grid
   \param    fit_limit   flag to match perfectly atmospheric limit
   \param    z_ground    altitude of the ground

   \author Patrick Eriksson
   \date   2000-09-14
*/
void space_geom(
               VECTOR&     z,
               Numeric&    l_step,
         const Numeric&    z_tan,
         const Numeric&    r_geoid,
         const Numeric&    atm_limit,
         const size_t&     fit_limit,
         const Numeric&    z_ground )
{
  Numeric  a, b;          // temporary values
  Numeric  llim;          // distance to atmospheric limit
  VECTOR   l;             // length from the tangent point

  //  cout << "l_step = " << l_step << "\n";

  // If LOS outside the atmosphere, return empty vector
  if ( z_tan >= atm_limit )
    resize(z,0);

  // Only through the atmosphere
  else if ( z_tan >= z_ground )
  {
    a    = r_geoid + atm_limit;
    b    = r_geoid + z_tan;
    llim = sqrt( a*a - b*b );        

    if ( !fit_limit ) 
    {
      if ( llim < l_step )   // Handle the rare case that llim < l_step
        l_step = llim*0.9999; // *0.9999 to avoid problem in interpolations
    }
    else
    {
      double n = ceil( llim / l_step + 1.0 );  
      l_step = 0.9999*llim/(n-1); // *0.9999 to avoid problem in interpolations
    }
    linspace( l, 0, llim, l_step );

    //    cout << "l, " << l.size() << ": " << l << "\n";

    //    z = sqrt(b*b+emult(l,l)) - r_geoid;
    ele_mult(l,l,l);
    //    cout << "l^2, " << l.size() << ": " << l << "\n";

    // Add to l a vector of same size with contents b*b:
    add( VECTOR( l.size(), b*b ), l ); 
    //    cout << "l added, " << l.size() << ": " << l << "\n";

    transf( l, sqrt, l );
    //    cout << "l sqrt, " << l.size() << ": " << l << "\n";

    resize(z,l.size());
    setto( z, -r_geoid );
    add( l, z ); 
    //    cout << "z, " << z.size() << ": " << z << "\n";
  }   

  // Intersection with the ground
  else
  {
    // Determine the "zenith angle" at ground level and call upward function 
    Numeric za = RAD2DEG*asin((r_geoid+z_tan)/(r_geoid+z_ground));
    upward_geom( z, l_step, r_geoid, z_ground, za, atm_limit, 0 );    
  }
}



//// space_refr ////////////////////////////////////////////////////////////
/**
   Calculates the LOS for observations from space with refraction.

   This function handles only a single zenith angle.

   \retval   p           pressures along the LOS
   \retval   l_step      step length along the LOS
   \param    c           LOS constant
   \param    z_tan       tangent altitude
   \param    r_geoid     local geoid curvature
   \param    atm_limit   maximum altitude of the absorption grid
   \param    z_ground    altitude of the ground
   \param    p_abs       absorption pressure grid
   \param    z_abs       absorption altitude grid
   \param    refr_index  refrective index corresponding to p_abs
   \param    l_step_refr step length for intermediate LOS

   \author Patrick Eriksson
   \date   2000-10-02
*/
void space_refr(
               VECTOR&     p,
               Numeric&    l_step,
         const Numeric&    c,
         const Numeric&    z_tan,
         const Numeric&    r_geoid,
         const Numeric&    atm_limit,
         const Numeric&    z_ground,
         const VECTOR&     p_abs,
         const VECTOR&     z_abs,
         const VECTOR&     refr_index,
               Numeric     l_step_refr )
{
  // If LOS outside the atmosphere, return empty vector
  if ( z_tan >= atm_limit )
    resize(p,0);

  // Only through the atmosphere
  else if ( z_tan >= z_ground )
  {
    VECTOR zs;
    space_geom( zs, l_step_refr, z_tan, r_geoid, atm_limit, 1, z_ground );
    geom2refr( p, l_step, zs, c, z_tan, r_geoid, atm_limit, p_abs, z_abs, 
                                                    refr_index, l_step_refr );
  }   

  // Intersection with the ground
  else
  {
    // Determine the "zenith angle" at ground level and call upward function 
    Numeric n_ground = interpz( p_abs, z_abs, refr_index, z_ground );
    Numeric za_ground = RAD2DEG*asin(c/n_ground/(r_geoid+z_ground));
    upward_refr( p, l_step, c, z_ground, za_ground, r_geoid, atm_limit, p_abs,
                                              z_abs, refr_index, l_step_refr );
  }
}



////////////////////////////////////////////////////////////////////////////
//   Sub-functions to losBasic 
////////////////////////////////////////////////////////////////////////////

//// los_space /////////////////////////////////////////////////////////////
/**
   Sets up the LOS structure for observations from space.

   Both limb and nadir observations are handled.   

   \retval   los         los structure (for all zenith angles)
   \param    z_plat      platform altitude
   \param    za          zentith angles
   \param    l_step      maximum step length along the LOS
   \param    atm_limit   maximum altitude of the absorption grid
   \param    z_ground    altitude of the ground
   \param    r_geoid     the local geoid radius
   \param    refr        flag for refraction (0=no refraction)
   \param    z_abs       absorption altitude grid
   \param    p_abs       absorption pressure grid
   \param    refr_index  refrective index corresponding to p_abs
   \param    l_step_refr step length for intermediate LOS

   \author Patrick Eriksson
   \date   2000-09-14
*/
void los_space(
                    LOS&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const Numeric&    atm_limit,
              const Numeric&    z_ground,
              const Numeric&    r_geoid,
              const size_t&     refr,
              const VECTOR&     z_abs,
              const VECTOR&     p_abs,
              const VECTOR&     refr_index,
              const Numeric&    l_step_refr )
{
  Numeric z_tan;           // the tangent altitude
  Numeric c;               // LOS constant

  // Set all step lengths to the user defined value as a first guess.
  // Note that the step length can be changed in SPACE_GEOM/REFR.
  setto( los.l_step, l_step );

  // Loop the zenith angles
  for ( size_t i=0; i<za.size(); i++ )
  { 
    if ( refr )   // Refraction
    {
      c = refr_constant( r_geoid, za[i], z_plat, p_abs, z_abs, refr_index );
      z_tan = ztan_refr( c, za[i], z_plat, z_ground, p_abs, z_abs, refr_index);
      space_refr( los.p[i], los.l_step[i], c, z_tan, r_geoid, atm_limit, 
                             z_ground, p_abs, z_abs, refr_index, l_step_refr );
    }
    else          // Geometrical calculations
    {
      z_tan = ztan_geom( za[i], z_plat );
      space_geom( los.p[i], los.l_step[i], z_tan, r_geoid, atm_limit, 0, 
                                                                   z_ground );
      // The geometrical functions return altitudes, not pressures
      z2p( los.p[i], z_abs, p_abs, los.p[i] );
    }

    los.start[i]  = los.p[i].size() - 1;
    los.stop[i]   = los.start[i];
    if ( z_tan >= z_ground )
      los.ground[i] = -1;          // no ground intersection
    else
      los.ground[i] = 0;           // ground at index 0
  }
}



//// los_inside ////////////////////////////////////////////////////
/**
   Sets up the LOS for measurements from points inside the atmosphere.

   Both upward and downward observations are handled.

   \retval   los         los structure (for all zenith angles)
   \param    z_plat      platform altitude
   \param    za          zentith angles
   \param    l_step      maximum step length along the LOS
   \param    atm_limit   maximum altitude of the absorption grid
   \param    z_ground    altitude of the ground
   \param    r_geoid     the local geoid radius
   \param    refr        flag for refraction (0=no refraction)
   \param    z_abs       absorption altitude grid
   \param    p_abs       absorption pressure grid
   \param    refr_index  refrective index corresponding to p_abs
   \param    l_step_refr step length for intermediate LOS

   \author Patrick Eriksson
   \date   2000-09-14
*/
void los_inside(
                    LOS&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const Numeric&    atm_limit,
              const Numeric&    z_ground,
              const Numeric&    r_geoid,
              const bool&       refr,
              const VECTOR&     z_abs,
              const VECTOR&     p_abs,
              const VECTOR&     refr_index,
              const Numeric&    l_step_refr )
{
  Numeric z_tan, za_ground; 
  Numeric a, b, c, l1;         // see below (c can be the LOS constant)
  size_t  n = za.size();        // the number of zenith angles

  // Set all step lengths to the user defined value as a first guess.
  // Note that the step length can be changed in the sub-functions.
  setto( los.l_step, l_step);

  // Loop the zenith angles
  for ( size_t i=0; i<n; i++ )
  { 
    // Calculate the LOS constant
    if ( refr )
      c = refr_constant( r_geoid, za[i], z_plat, p_abs, z_abs, refr_index );

    // Upward
    if ( za[i] <= 90 )
    {
      if ( refr )
        upward_refr( los.p[i], los.l_step[i], c, z_plat, za[i], r_geoid, 
                           atm_limit, p_abs, z_abs, refr_index, l_step_refr );
      else
        upward_geom( los.p[i], los.l_step[i], r_geoid, z_plat, za[i], 
                                                               atm_limit, 0 );
      los.start[i]  = los.p[i].size() - 1;   // length of LOS
      los.stop[i]   = 0;                     // stop index here always 0
      los.ground[i] = -1;                    // no ground intersection
    }

    // Downward
    else
    {
      // Calculate the tangent altitude
      if ( refr )
        z_tan = ztan_refr( c, za[i], z_plat, z_ground, p_abs,z_abs,refr_index);
      else
	z_tan = ztan_geom( za[i], z_plat );

      // Only through the atmosphere
      if ( z_tan >= z_ground )
      {
        // Calculate the distance between the platform and the tangent point
        if ( refr )
        {
          // Calculate a first refractive LOS with short step length and 
          // interpolate to determine the distance
	  VECTOR  zs, ps, ls;  
	  Numeric l_temp=l_step_refr;
	  space_geom( zs, l_temp, z_tan, r_geoid, atm_limit, 1, z_ground );
	  geom2refr( ps, l_temp, zs, c, z_tan, r_geoid, atm_limit, p_abs, 
                                                   z_abs, refr_index, l_temp );
	  linspace( ls, 0, l_temp*(ps.size()-1), l_temp );
	  l1 = interp_lin( ps, ls, z_plat );  
	}
        else
	{
          // Use geometry
	  a      = r_geoid + z_plat;   // help variable
	  b      = r_geoid + z_tan;    // help variable
	  l1     = sqrt(a*a-b*b);           // distance platform-tangent point
        }

	// Adjust l_step downwards to get an integer number of steps
	los.stop[i]   = (size_t) ceil( l1 / l_step + 1.0 ) - 1;  
	los.l_step[i] = l1 / ( (Numeric)los.stop[i] );
        if ( refr )
	  space_refr( los.p[i], los.l_step[i], c, z_tan, r_geoid, atm_limit, 
                             z_ground, p_abs, z_abs, refr_index, l_step_refr );
	else
	  space_geom( los.p[i], los.l_step[i], z_tan, r_geoid, atm_limit, 0, 
                                                                    z_ground );
        los.start[i]  = los.p[i].size() - 1;
        los.ground[i] = -1;                // no gound intersection
      }   
    
      // Intersection with the ground
      else
      {
        if ( refr )
        {
          // Calculate a first refractive LOS with short step length and 
          // interpolate to determine the distance
	  VECTOR  zs, ps, ls;  
	  Numeric l_temp=l_step_refr;
          // Determine the "zenith angle" at ground and call upward function 
          Numeric n_ground = interpz( p_abs, z_abs, refr_index, z_ground );
          za_ground = RAD2DEG*asin(c/n_ground/(r_geoid+z_ground));
          upward_geom( zs, l_temp, r_geoid, z_ground, za_ground, atm_limit, 1);
	  geom2refr( ps, l_temp, zs, c, z_tan, r_geoid, atm_limit, p_abs, 
                                                  z_abs, refr_index, l_temp );
	  linspace( ls, 0, l_temp*(ps.size()-1), l_temp );
	  l1 = interp_lin( ps, ls, z_plat );  
	}
        else
	{
          a      = r_geoid + z_ground;
	  b      = r_geoid + z_plat;
	  c      = r_geoid + z_tan;
	  l1     = sqrt(b*b-c*c) - sqrt(a*a-c*c); // distance platform-ground
        }
        // Adjust l_step downwards to get an integer number of steps
	los.stop[i]   = (size_t) ceil( l1 / l_step );
	los.l_step[i] = l1 / ( (double)los.stop[i] );
        if ( refr )
	  upward_refr( los.p[i], los.l_step[i], c, z_ground, za_ground, 
                   r_geoid, atm_limit, p_abs, z_abs, refr_index, l_step_refr );
        else
          upward_geom( los.p[i], los.l_step[i], r_geoid, z_ground, 
                                             RAD2DEG*asin(c/a), atm_limit, 0 );
        los.start[i]  = los.p[i].size() - 1;
        los.ground[i] = 0;                // ground at index 0
      }
    }

    // The geometrical functions return altitudes, not pressures
    if ( !refr )
      z2p( los.p[i], z_abs, p_abs, los.p[i] );
  }
}



////////////////////////////////////////////////////////////////////////////
//   Workspace methods
////////////////////////////////////////////////////////////////////////////

void r_geoidStd( Numeric&    r_geoid )
{
   r_geoid = EARTH_RADIUS;
}



void r_geoidWGS84( 
              Numeric&   r_geoid,
        const Numeric&   latitude,
        const Numeric&   obsdirection )
{
  const Numeric rq = 6378.138e3, rp = 6356.752e3;
        Numeric a, b, rns, rew;

  // Calculate NS and EW radius
  a   = cos(latitude*DEG2RAD);
  b   = sin(latitude*DEG2RAD);
  rns = rq*rq*rp*rp/pow(rq*rq*a*a+rp*rp*b*b,1.5);
  rew = rq*rq/sqrt(rq*rq*a*a+rp*rp*b*b);

  // Calculate the radius in the observation direction
  a   = cos(obsdirection*DEG2RAD);
  b   = sin(obsdirection*DEG2RAD);
  r_geoid   = 1/(a*a/rns+b*b/rew);
}



void losCalc(       LOS&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const int&        refr,
              const Numeric&    l_step_refr,
              const VECTOR&     refr_index,
              const Numeric&    z_ground,
              const Numeric&    r_geoid )
{     
  size_t n = za.size();  // number of zenith angles

  //  cout << "losCalc/l_step = " << l_step << "\n";

  // Some checks                                                      
  if ( z_ground < z_abs[0] )
    throw runtime_error(
      "There is a gap between the ground and the lowest absorption altitude");
  if ( z_plat < z_ground )
    throw runtime_error("Your platform altitude is below the ground");
  if ( z_plat < z_abs[0] )  
    throw runtime_error(
      "The platform cannot be below the lowest absorption altitude");

  // Reallocate the l_step, ground, start and stop vectors
  resize( los.p,      n	);
  resize( los.l_step, n	);
  resize( los.ground, n	);
  resize( los.start,  n	);
  resize( los.stop,   n	);

  // Print messages
  if ( refr == 0 )
    out2 << "  Calculating line of sights WITHOUT refraction.\n";
  else if ( refr == 1 )
    out2 << "  Calculating line of sights WITH refraction.\n";
  else
    throw runtime_error("The refraction flag can only be 0 or 1.");
  out3 << "     z_plat: " << z_plat/1e3 << " km\n";
  if ( n == 1 )
    out3 << "         za: " << za[0] << " degs.\n";
  else
  {
    out3 << "     min za: " << min(za) << " degs.\n";
    out3 << "     max za: " << max(za) << " degs.\n";
  }

  // Determine the upper limit of the atmosphere
  const Numeric atm_limit = last(z_abs);

  // Two cases, from space or from within the atmosphere
  if ( z_plat >= atm_limit )
    los_space( los, z_plat, za, l_step, atm_limit, z_ground, r_geoid, refr, 
                                     z_abs, p_abs, refr_index, l_step_refr );
  else
    los_inside( los, z_plat, za, l_step, atm_limit, z_ground, r_geoid, refr, 
                                     z_abs, p_abs, refr_index, l_step_refr );
}



void losNoRefraction(
                    LOS&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const Numeric&    z_ground,
              const Numeric&    r_geoid )
{
  losCalc( los, z_plat, za, l_step, p_abs, z_abs, 0, 0, VECTOR(0) , z_ground,
                                                                    r_geoid ); 
}



void losUpward(
                    LOS&        los,
              const Numeric&    z_plat,
              const VECTOR&     za,
              const Numeric&    l_step,
              const VECTOR&     p_abs,
              const VECTOR&     z_abs,
              const Numeric&    r_geoid )
{
  if ( max(za) > 90 )
    throw runtime_error("At least one zenith angle > 90 degrees, that is, not upwards. Use losCalc or losNoRefraction.");

  losCalc( los, z_plat, za, l_step, p_abs, z_abs, 0, 0, VECTOR(0) , z_plat,
                                                                    r_geoid ); 
}



void sourceCalc(
                    ARRAYofMATRIX&   source,
              const LOS&             los,   
              const VECTOR&          p_abs,
              const VECTOR&          t_abs,
              const VECTOR&          f_mono )
{     
        VECTOR   tlos;                  // temperatures along the LOS
  const size_t   nza=los.start.size();  // the number of zenith angles  
  const size_t   nf=f_mono.size();      // the number of frequencies
        size_t   nlos;                  // the number of pressure points
        MATRIX   b;                     // the Planck function for TLOS  
        size_t   iv, ilos;              // frequency and LOS point index

  out2 << "  Calculating the source function for LTE and no scattering.\n";
 
  // Resize the source array
  resize(source,nza);

  // Loop the zenith angles and:
  //  1. interpolate the temperature
  //  2. calculate the Planck function for the interpolated temperatures
  //  3. take the mean of neighbouring Planck values
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=0; i<nza; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();

    if ( los.p[i].size() > 0 )
    {
      resize( tlos, los.p[i].size() );
      interpp( tlos, p_abs, t_abs, los.p[i] );
      nlos = tlos.size();
      resize( b, f_mono.size(), tlos.size() );
      planck( b, f_mono, tlos );
      resize(source[i],nf,nlos-1);
      for ( ilos=0; ilos<(nlos-1); ilos++ )
      {
        for ( iv=0; iv<nf; iv++ )
          source[i][iv][ilos] = ( b[iv][ilos] + b[iv][ilos+1] ) / 2.0;
      }
    }
  }  
  out3 << "\n";
}



void transCalc(
                    ARRAYofMATRIX&   trans,
              const LOS&             los,   
              const VECTOR&          p_abs,
              const MATRIX&          abs )
{    
  // Some variables
  const size_t   n = los.start.size(); // the number of zenith angles
  const size_t   nf = abs.nrows();     // the number of frequencies
        size_t   np;                   // the number of pressure points
        size_t   row, col;             // counters
        MATRIX   abs2 ;                // matrix to store interpolated absorp.
       Numeric   w;                    // = -l_step/2

  out2 << "  Calculating transmissions WITHOUT scattering.\n";
 
  // Resize the transmission array
  resize(trans,n);

  // Loop the zenith angles and:
  //  1. interpolate the absorption
  //  2. calculate the transmission using the mean absorption between points
  out3 << "    Zenith angle nr:\n      ";
  for (size_t i=0; i<n; i++ ) 
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    np = los.p[i].size();
    if ( np > 0 )
    {
      resize( abs2, abs.nrows(), los.p[i].size() );
      //interp_lin_row( abs2, p_abs, abs, los.p[i] );
      resize( abs2, nf, np );
      interp_lin_matrix( abs2, p_abs, abs, los.p[i] );
      resize(trans[i], nf, np-1 );
      w  =  -0.5*los.l_step[i];
      for ( row=0; row<nf; row++ )
      {
        for ( col=0; col<(np-1); col++ )
          trans[i][row][col] = exp( w * ( abs2[row][col]+abs2[row][col+1]) );
      }
    }
  }    
  out3 << "\n";
}



void y_spaceStd(
                    VECTOR&   y_space,
              const VECTOR&   f,
              const string&   choice )
{
  resize( y_space, f.size() );

  if ( choice == "zero" )
  {
    setto(y_space,0.0);
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
    out2 << "  Setting y_space to blackbody radiation corresponding to the Sun temperature\n";
  }
  else
    throw runtime_error(
      "Possible choices for Y_SPACE are \"zero\", \"cbgr\" and \"sun\".");

}



void yRte (
                    VECTOR&          y,
              const LOS&             los,   
              const VECTOR&          f_mono,
              const VECTOR&          y_space,
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground,
              const Numeric&         t_ground )
{
  // Some variables
  const size_t   n=los.start.size();  // Number of zenith angles 
  const size_t   nf=f_mono.size();    // Number of frequencies 
        VECTOR   y_tmp(nf);           // Temporary storage for spectra
        size_t   iy0=0;               // Reference index for output vector
        size_t   iy;                  // Frequency index

  out2 << "  Integrating the radiative transfer eq. WITHOUT scattering.\n";

  // Resize y
  resize( y, nf*n );
        
  // Set up vector for ground blackbody radiation if any ground intersection
  // Check also if the ground emission vector has the correct length
  VECTOR   y_ground(f_mono.size()); 
  if ( any_ground(los.ground) )  
  {
    out2 << "  There are intersections with the ground.\n";
    planck( y_ground, f_mono, t_ground );
    if ( e_ground.size() != nf )
      throw runtime_error("The frequency and ground emission vectors have different lengths.");
  }

  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for ( size_t i=0; i<n; i++ )
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Iteration is done in seperate function    
    rte( y_tmp, los.start[i], los.stop[i], trans[i], 
                 source[i], y_space, los.ground[i], e_ground, y_ground);

    // Move values to output vector
    for ( iy=0; iy<nf; iy++ )    
      y[iy0+iy] = y_tmp[iy];
    iy0 += nf;                    // update iy0
  }
  out3 << "\n";
}



void yRteNoGround (
                    VECTOR&          y,
              const LOS&             los,   
              const VECTOR&          f_mono,
              const VECTOR&          y_space,
              const ARRAYofMATRIX&   source,
              const ARRAYofMATRIX&   trans )
{
  if ( any_ground(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  yRte( y, los, f_mono, y_space, source, trans, VECTOR(0), 0.0 );
}



void yBl (
                    VECTOR&          y,
              const LOS&             los,   
              const ARRAYofMATRIX&   trans,
              const VECTOR&          e_ground )
{
  // Some variables
  const size_t   n=los.start.size();    // Number of zenith angles 
  const size_t   nf=trans[0].nrows();    // Number of frequencies 
        size_t   iy0=0;                 // Reference index for output vector
        size_t   iy;                    // Frequency index
        VECTOR   y_tmp;                 // Temporary storage for spectra

  out2 << "  Calculating total transmission WITHOUT scattering.\n";

  // Resize y and set to 1
  resize( y, nf*n );
  setto( y, 1.0 );

  // Check if the ground emission vector has the correct length
  if ( any_ground(los.ground) )  
  {
    out2 << "  There are intersections with the ground.\n";
    if ( e_ground.size() != nf )
      throw runtime_error("The frequency and ground emission vectors have different lengths.");
  }
        
  // Loop zenith angles
  out3 << "    Zenith angle nr:\n      ";
  for ( size_t i=0; i<n; i++ )
  {
    if ( (i%20)==0 )
      out3 << "\n      ";
    out3 << " " << i; cout.flush();
    
    // Iteration is done in seperate function    
    bl( y_tmp, los.start[i], los.stop[i], trans[i], los.ground[i], e_ground );

    // Move values to output vector
    for ( iy=0; iy<nf; iy++ )    
      y[iy0+iy] = y_tmp[iy];
    iy0 += nf;                    // update iy0
  }
  out3 << "\n";
}



void yBlNoGround (
                    VECTOR&          y,
              const LOS&             los,   
              const ARRAYofMATRIX&   trans )
{
  if ( any_ground(los.ground) )  
    throw runtime_error("There is at least one intersection with the ground and this function cannot be used.");

  yBl( y, los, trans, VECTOR(0) );
}



void yTB (
                    VECTOR&          y,
              const VECTOR&          f_sensor,
              const VECTOR&          za_sensor )
{
  const size_t   nf  = f_sensor.size();
  const size_t   nza = za_sensor.size();
  const size_t   ny  = y.size();
        size_t   i0;
  const Numeric  a = PLANCK_CONST/BOLTZMAN_CONST;
  const Numeric  b = 2*PLANCK_CONST/(SPEED_OF_LIGHT*SPEED_OF_LIGHT);
        Numeric  c,d;

  if ( max(y) > 1e-4 )  
    throw runtime_error("The spectrum is not in expected intensity unit (impossible value detected).");

  if ( nf*nza != ny )  
    throw runtime_error("The length of y does not match f_sensor and za_sensor.");

  out2 << "  Converts the spectrum to brightness (Planck) temperature.\n";

  for ( size_t i=0; i<nf; i++ )
  {
    c = a*f_sensor[i];
    d = b*f_sensor[i]*f_sensor[i]*f_sensor[i];
    for ( size_t j=0; j<nza; j++ )    
    {
      i0 = j*nf + i;
      y[i0] = c / ( log(d/y[i0]+1) );
    }
  }
}



void yTRJ (
                    VECTOR&          y,
              const VECTOR&          f_sensor,
              const VECTOR&          za_sensor )
{
  const size_t   nf  = f_sensor.size();
  const size_t   nza = za_sensor.size();
  const size_t   ny  = y.size();
        size_t   i0;
  const Numeric  a = SPEED_OF_LIGHT*SPEED_OF_LIGHT/(2*BOLTZMAN_CONST);
        Numeric  b;

  if ( max(y) > 1e-4 )  
    throw runtime_error("The spectrum is not in expected intensity unit (impossible value detected).");

  if ( nf*nza != ny )  
    throw runtime_error("The length of y does not match f_sensor and za_sensor.");

  out2 << "  Converts the spectrum to Rayleigh-Jean temperature.\n";

  for ( size_t i=0; i<nf; i++ )
  {
    b = a/(f_sensor[i]*f_sensor[i]);
    for ( size_t j=0; j<nza; j++ )    
    {
      i0 = j*nf + i;
      y[i0] = b * y[i0];
    }
  }
}



void yLoadCalibration (
                    VECTOR&          y,
              const VECTOR&          i_cal1,
              const VECTOR&          i_cal2,
              const VECTOR&          y_cal1,
              const VECTOR&          y_cal2,
              const VECTOR&          za_sensor )
{
  const size_t   nf  = i_cal1.size();
  const size_t   nza = za_sensor.size();
  const size_t   ny  = y.size();
        size_t   i0;
        Numeric  a;

  if ( max(y) > 1e-4 )  
    throw runtime_error("The spectrum is not in expected intensity unit (impossible value detected).");

  if ( (nf!=i_cal2.size()) || (nf!=y_cal1.size()) || (nf!=y_cal2.size()) )
    throw runtime_error("All the calibration vectors must have the same length.");
  if ( nf*nza != ny )  
    throw runtime_error("The length of y does not match za_sensor and the calibration data.");

  out2 << "  Performs a load switching calibration procedure.\n";

  for ( size_t i=0; i<nf; i++ )
  {
    a = (i_cal2[i]-i_cal1[i])/(y_cal2[i]-y_cal1[i]);
    for ( size_t j=0; j<nza; j++ )    
    {
      i0 = j*nf + i;
      y[i0] = i_cal1[i] + a * ( y[i0] - y_cal1[i] );
    }
  }
}



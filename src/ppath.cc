/* Copyright (C) 2002 Patrick Eriksson <patrick@rss.chalmers.se>

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
  \file   ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-02
  
  \brief  Functions releated to calculation of propagation paths.
  
  Functions to determine propagation paths for different atmospheric
  dimensionalities.

  The term propagation path is here shortened to ppath.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "mystring.h"
#include "logic.h"
#include "poly_roots.h"
#include "ppath.h"
#include "special_interp.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;





/*===========================================================================
  === Global variables
  ===========================================================================*/

// When calculating propagation path steps we want to check if the given start
// point is inside the grid range/cell described by the grid positions. 
// The end points of the path steps are normally at the end of a grid range,
// but numerical inaccuracy leads easily to that a point is found to be 
// slightly above or below the end of the grid range. If the point then is
// found to be outside the grid range, this will result in an error.
// The problem is especially tricky for the radius when the pressure surfaces
// are tilted (their radius changes with latitude and longitude).
// The problem is handled with the variable R_EPS, that gives the accepted 
// inaccuracy in radius when checking if a position is inside a grid cell. 
//
const Numeric R_EPS = 1e-3;





/*===========================================================================
  === Functions related to geometrical propagation paths
  ===========================================================================*/

//! geometrical_ppc
/*! 
   Calculates the propagation path constant for pure geometrical calculations.

   Both positive and negative zenith angles are handled.

   \return         Path constant.
   \param   r      Radius of the sensor position.
   \param   za     Zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Numeric geometrical_ppc( const Numeric& r, const Numeric& za )
{
  assert( r > 0 );
  assert( fabs(za) <= 180 );

  return r * sin( DEG2RAD * fabs(za) );
}



//! geompath_za_at_r
/*! 
   Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   For downlooking cases, the two points must be on the same side of 
   the tangent point.

   Both positive and negative zenith angles are handled.

   \return         Zenith angle at the point of interest.
   \param   ppc    Propagation path constant.
   \param   a_za   A zenith angle along the path on the same side of the 
                   tangent point as the point of interest.  
   \param   r      Radius of the point of interest.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Numeric geompath_za_at_r(
       const Numeric&   ppc,
       const Numeric&   a_za,
       const Numeric&   r )
{
  assert( ppc >= 0 );
  assert( fabs(a_za) <= 180 );
  assert( r >= ppc );

  Numeric za = RAD2DEG * asin( ppc / r );
  if( fabs(a_za) > 90 )
    { za = 180 - za; }
  if( a_za < 0 )
    { za = -za; }
  return za;
}



//! geompath_r_at_za
/*! 
   Calculates the zenith angle for a given radius along a geometrical 
   propagation path.

   Both positive and negative zenith angles are handled.

   \return         Radius at the point of interest.
   \param   ppc    Propagation path constant.
   \param   za     Zenith angle at the point of interest.

   \author Patrick Eriksson
   \date   2002-06-05
*/
Numeric geompath_r_at_za(
       const Numeric&   ppc,
       const Numeric&   za )
{
  assert( ppc >= 0 );
  assert( fabs(za) <= 180 );

  return ppc / sin( DEG2RAD * fabs(za) );
}



//! geompath_lat_at_za
/*!
   Calculates the latitude for a given zenith angle along a geometrical 
   propagation path.

   Positive and negative zenith angles are handled. A positive zenith angle
   means a movement towards higher latitudes.

   \return         The latitude of the second point.
   \param   za0    The zenith angle of the starting point.
   \param   lat0   The latitude of the starting point.
   \param   za     The zenith angle of the second point.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Numeric geompath_lat_at_za(
       const Numeric&   za0,
       const Numeric&   lat0,
       const Numeric&   za )
{
  assert( fabs(za0) <= 180 );
  assert( fabs(za) <= 180 );
  assert( ( za0 >= 0 && za >= 0 )  ||  ( za0 < 0 && za < 0 ) );

  return lat0 + za0 - za;
}



//! geompath_za_at_lat
/*!
   Calculates the zenith angle for a given latitude along a geometrical 
   propagation path.

   Positive and negative zenith angles are handled. A positive zenith angle
   means a movement towards higher latitudes.

   \return         The zenith angle for the second point.
   \param   za0    The zenith angle of the starting point.
   \param   lat0   The latitude of the starting point.
   \param   lat     The latitude of the second point.

   \author Patrick Eriksson
   \date   2002-07-03
*/
/* Not used:
Numeric geompath_za_at_lat(
       const Numeric&   za0,
       const Numeric&   lat0,
       const Numeric&   lat )
{
  assert( fabs(za0) <= 180 );
  assert( fabs(za) <= 180 );
  assert( ( za0 >= 0 && lat >= lat0 )  ||  ( lat < lat0 ) );

  return za0 - lat0 - lat;
}
*/


//! geompath_l_at_r
/*!
   Calculates the length from the tangent point for the given radius.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Length along the path from the tangent point.
   \param   ppc    Propagation path constant.
   \param   r      Radius of the point of concern.
g
   \author Patrick Eriksson
   \date   2002-05-20
*/
Numeric geompath_l_at_r(
       const Numeric&   ppc,
       const Numeric&   r )
{
  assert( ppc >= 0 );
  assert( r >= ppc );
  
  // Double is hard-coded here to improve accuracy
  double a=ppc*ppc, b=r*r;

  return sqrt( b - a );
}



//! geompath_r_at_l
/*!
   Calculates the radius for a distance from the tangent point.

   The tangent point is either rwal or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Radius. 
   \param   ppc    Propagation path constant.
   \param   l      Length from the tangent point.

   \author Patrick Eriksson
   \date   2002-05-20
*/
Numeric geompath_r_at_l(
       const Numeric&   ppc,
       const Numeric&   l )
{
  assert( ppc >= 0 );
  assert( l >= 0 );
  
  // Double is hard-coded here to improve accuracy
  double a=ppc*ppc, b=l*l;

  return sqrt( b + a );
}



//! geompath_r_at_lat
/*!
   Calculates the radius for a given latitude.

   \return         Radius at the point of interest.
   \param   ppc    Propagation path constant.
   \param   lat0   Latitude at some other point of the path.
   \param   za0    Zenith angle for the point with latitude lat0.
   \param   lat    Latitude of the point of interest.

   \author Patrick Eriksson
   \date   2002-06-05
*/
Numeric geompath_r_at_lat(
       const Numeric&   ppc,
       const Numeric&   lat0,
       const Numeric&   za0,
       const Numeric&   lat )
{
  assert( ppc >= 0 );
  assert( fabs(za0) <= 180 );
  assert( ( za0 >= 0 && lat >= lat0 )  ||  ( za0 <= 0 && lat <= lat0 ) );

  // Zenith angle at the new latitude
  const Numeric za = za0 + lat0 -lat;

  return geompath_r_at_za( ppc, za );
}



//! geompath_from_r1_to_r2
/*!
   Determines radii, latitudes and zenith angles between two points of a 
   propagation path.

   Both start and end point are included in the returned vectors.

   \param   r      Output: Radius of propagation path points.
   \param   lat    Output: Latitude of propagation path points.
   \param   za     Output: Zenith angle of propagation path points.
   \param   lstep  Output: Distance along the path between the points. 
   \param   ppc    Propagation path constant.
   \param   r1     Radius for first point.
   \param   lat1   Latitude for first point.
   \param   za1    Zenith angle for first point.
   \param   r2     Radius for second point.
   \param   lmax   Length criterion for distance between path points.
                   A negative value means no length criterion.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void geompath_from_r1_to_r2( 
	     Vector&    r,
	     Vector&    lat,
	     Vector&    za,
	     Numeric&   lstep,
       const Numeric&   ppc,
       const Numeric&   r1,
       const Numeric&   lat1,
       const Numeric&   za1,
       const Numeric&   r2,
       const Numeric&   lmax )
{
  // Calculate length along the path for point 1 and 2.
  const Numeric l1 =  geompath_l_at_r( ppc, r1 );
  const Numeric l2 =  geompath_l_at_r( ppc, r2 );
  
  // Calculate needed number of steps, considering a possible length criterion
  Index n;
  if( lmax > 0 )
    {
      // The absolute value of the length distance is needed here
      n = Index( ceil( fabs( l2 - l1 ) / lmax ) );
    }
  else
    { n = 1; }

  // Length of path steps (note that lstep here can be negative)
  lstep = ( l2 - l1 ) / n;

  // Allocate vectors and put in point 1
  r.resize(n+1);
  lat.resize(n+1);
  za.resize(n+1);
  r[0]   = r1;
  lat[0] = lat1;
  za[0]  = za1;

  // Loop steps (beside last) and calculate radius and zenith angle
  for( Index i=1; i<n; i++ )
    {
      r[i]   = geompath_r_at_l( ppc, l1 + lstep * i );
      za[i]  = geompath_za_at_r( ppc, za1, r[i] );
    }

  // For maximum accuracy, set last radius to be exactly r2.
  r[n]   = r2;
  za[n]  = geompath_za_at_r( ppc, za1, r[n] );

  // Ensure that zenith and nadir observations keep their zenith angle
  if( za1 == 0  ||  fabs(za1) == 180 )
    { za = za1; }

  // Calculate latitudes
  for( Index i=1; i<=n; i++ )
    { lat[i] = geompath_lat_at_za( za1, lat1, za[i] ); }

  // Take absolute value of lstep
  lstep = fabs( lstep );
}



//! za_geom2other_point
/*!
   Calculates the zenith angle for the geometrical propagation path between
   two specified points.

   The returned zenith angle is valid for point 1. That is, the propagation
   path goes from point 1 to point 2.

   \return         Zenith angle.
   \param   r1     Radius for point 1.
   \param   lat1   Latiytude for point 1.
   \param   r2     Radius for point 2.
   \param   lat2   Latitude for point 2.

   \author Patrick Eriksson
   \date   2002-07-03
*/
Numeric za_geom2other_point(
       const Numeric&   r1,
       const Numeric&   lat1,
       const Numeric&   r2,
       const Numeric&   lat2 )
{
  if( lat2 == lat1 )
    {
      if( r1 <= r2 )
	{ return 0; }
      else
	{ return 180; }
    }
  else
    {
      // Absolute value of the latitude difference
      const Numeric dlat = fabs( lat2 - lat1 );

      // The zenith angle is obtained by a combination of the lawes of sine
      // and cosine.
      Numeric za = dlat + RAD2DEG * asin( r1 * sin( DEG2RAD * dlat ) / 
		 sqrt( r1*r1 + r2*r2 - 2 * r1 * r2 * cos( DEG2RAD * dlat ) ) );

      // Consider viewing direction
      if( lat2 < lat1 )
	{ za = -za; }

      return za;
    }
}





/*===========================================================================
  === Functions related to propagation paths with refraction
  ===========================================================================*/

//! refraction_ppc
/*! 
   Calculates the propagation path constant for cases with refraction.

   Both positive and negative zenith angles are handled.

   \return         Path constant.
   \param   r      Radius of the sensor position.
   \param   za     Zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Numeric refraction_ppc( 
        const Numeric& r, 
	const Numeric& za, 
	const Numeric& refr_index )
{
  assert( r > 0 );
  assert( fabs(za) <= 180 );

  return r * refr_index * sin( DEG2RAD * fabs(za) );
}



//! refraction_gradient_2d
/*! 
   Calculates the radial and latitudinal derivative of the refractive
   index in a 2D grid cell.

   The refractive index is assumed to be a bi-linear function (as a
   function of radius and latitude).

   \param   dndr     Out: Radial derivative.
   \param   dndlat   Out: Latitudinal derivative.
   \param   r1       Radius of lower-left corner.
   \param   r4       Radius of upper-left corner (r4 > r1).
   \param   c2       Slope of lower pressure surface [m/deg].
   \param   c4       Slope of upper pressure surface [m/deg].
   \param   lat1     Latitude of left end face.
   \param   lat3     Latitude of right end face (lat3 > lat1).
   \param   n1       Refractive index at lower-left corner.
   \param   n2       Refractive index at lower-right corner.
   \param   n3       Refractive index at upper-right corner.
   \param   n4       Refractive index at upper-left corner.
   \param   r        Radius for point of concern.
   \param   lat      Latitude for point of concern.

   \author Patrick Eriksson
   \date   2002-11-18
*/
void refraction_gradient_2d( 
	      Numeric&   dndr,
	      Numeric&   dndlat,
        const Numeric&   r1, 
	const Numeric&   r4, 
	const Numeric&   c2,
	const Numeric&   c4,
	const Numeric&   lat1,
	const Numeric&   lat3,
        const Numeric&   n1, 
	const Numeric&   n2, 
	const Numeric&   n3,
	const Numeric&   n4,
	const Numeric&   r,
        const Numeric&   lat )
{
  // Help variables to avoid duplication of calculations
  const Numeric   xlat  = lat - lat1;
  const Numeric   rlow  = r1 + xlat * c2;
  const Numeric   rhigh = r4 + xlat * c4;
  const Numeric   dr    = rhigh - rlow;
  const Numeric   dlat  = lat3 - lat1;

  assert( r1 < r4 );
  assert( lat1 < lat3 );
  assert( lat >= lat1  &&  lat <= lat3 );
  assert( r >= rlow  &&  r <= rhigh );

  // Fractional distance for latitude
  Numeric   fd = xlat / dlat;

  // Derivative in the radius direction
  dndr = ( ( fd*n3 + (1-fd)*n4 ) - ( fd*n2 + (1-fd)*n1 ) ) / dr; 

  // Fractional distance for radius
  fd   = ( r - rlow ) / dr; 

  // Derivative in the latitude direction
  dndr = ( ( fd*n3 + (1-fd)*n2 ) - ( fd*n4 + (1-fd)*n1 ) ) / dlat; 
}





/*===========================================================================
  === Functions related to slope and tilt of the ground and pressure surfaces
  ===========================================================================*/

//! psurface_slope_2d
/*!
   Calculates the radial slope of the ground or a pressure surface for 2D.

   The radial slope is here the derivative of the radius with respect to the
   latitude. The unit is accordingly m/degree.

   Note that the radius is defined to change linearly between grid points.

   Note also that the slope is always calculated with respect to increasing
   latitudes, independently of the upwards argument. The upwards argument is
   only used to determine which grid range that is of interest when the
   position is exactly on top of a grid point.

   For a point exactly on a grid value it is not clear if it is the range 
   below or above that is of interest. The input argument upward is used to 
   resolve such cases, where upward == 1 means that it is the range above
   that is of interest.

   \return              The radial slope [m/degree]
   \param   lat_grid    The latitude grid.
   \param   r_geoid     Radius of the geoid for the latitude dimension.
   \param   z_surf      Geometrical altitude of the ground, or the pressure
                        surface of interest, for the latitide dimension
   \param   gp          Latitude grid position for the position of interest
   \param   upwards     See above.

   \author Patrick Eriksson
   \date   2002-06-03
*/
Numeric psurface_slope_2d(
	ConstVectorView   lat_grid,	      
	ConstVectorView   r_geoid,
	ConstVectorView   z_surf,
        const GridPos&    gp,
        const Index&      upwards )
{
  assert( is_bool( upwards ) );

  Index i1 = gridpos2gridrange( gp, upwards );
  const Numeric r1 = r_geoid[i1] + z_surf[i1];
  const Numeric r2 = r_geoid[i1+1] + z_surf[i1+1];
  return ( r2 - r1 ) / ( lat_grid[i1+1] - lat_grid[i1] );
}



//! psurface_tilt_2d
/*!
   Calculates the angular tilt of the ground or a pressure surface for 2D.

   Note that the tilt value is a local value. The tilt for a constant
   slope value, is different for different radii.

   \return        The angular tilt.
   \param    r    The radius for the surface at the point of interest.
   \param    c    The radial slope, as returned by psurface_slope_2d.

   \author Patrick Eriksson
   \date   2002-06-03
*/
Numeric psurface_tilt_2d(
        const Numeric&   r,
        const Numeric&   c )
{
  // The tilt (in radians) is c/r if c is converted to m/radian. So we get
  // conversion RAD2DEG twice
  return   RAD2DEG * RAD2DEG * c / r;
}



//! is_los_downwards_2d
/*!
   Determines if a line-of-sight is downwards compared to the angular tilt
   of the ground or a pressure surface.

   For example, this function can be used to determine if the line-of-sight
   goes into the ground for a starting point exactly on the ground radius.
  
   As the radius of the ground and pressure surfaces varies as a function of
   latitude, it is not clear if a zenith angle of 90 is above or below e.g.
   the ground.
 
   \return         Boolean that is true if LOS is downwards.
   \param   za     Zenith angle of line-of-sight.
   \param   tilt   Angular tilt of the ground or the pressure surface (as
                   returned by psurface_tilt_2d)

   \author Patrick Eriksson
   \date   2002-06-03
*/
bool is_los_downwards_2d( 
        const Numeric&   za,
        const Numeric&   tilt )
{
  assert( fabs(za) <= 180 );

  // Yes, it shall be -tilt in both cases, if you wonder.
  if( za > (90-tilt)  ||  za < (-90-tilt) )
    { return true; }
  else
    { return false; }
}



//! psurface_crossing_2d
/*!
   Calculates the angular distance to a crossing of a pressure surface
   or the ground.

   The function solves the problem mentioned above for a pressure
   surface, or the ground, where the radius changes quadratically as a
   function of latitude. No analytical solution to the original
   problem has been found. The problem involves sine and cosine of the
   latitude difference and these functions are replaced with their
   Taylor expansions where the two first terms are kept. This should
   be OK as in practical situations, the latitude difference inside a
   grid cell should not exceed 2 degrees, and the accuracy should be
   sufficient for values up to 3 degrees.

   The problem and its solution is further described in AUG. See the
   chapter on propagation paths.

   Both positive and negative zenith angles are handled. 

   The function only looks for crossings in the forward direction of
   the given zenith angle. This means that if r>r0 and the absolute
   value of the zenith angle is < 90, no crossing will be found (if
   not the slope of the pressure surface happen to be very strong).
   
   If the given path point is on the pressure surface (r=r0), the
   solution 0 is rejected.
 
   The latitude difference is set to 999 of no possible value is found, but
   there will normally be some value < 360 degrees. 

   The variable names below are the same as in AUG.

   \return         The angular distance to the crossing.
   \param   rp     Radius of a point of the path inside the grid cell
   \param   za     Zenith angle of the path at r.
   \param   r0     Radius of the pressure surface or the ground at the
                   latitude of r.
   \param   c1     Linear slope term, as returned by psurface_slope_2d.
   \param   c2     Quadratic slope term.

   \author Patrick Eriksson
   \date   2002-06-07
*/
Numeric psurface_crossing_2d(
        const Numeric&   rp,
        const Numeric&   za,
        const Numeric&   r0,
              Numeric    c1,
              Numeric    c2 )
{
  assert( fabs(za) <= 180 );

  const Numeric no_crossing = 999;

  // Handle the cases of za=0 and za=180. 
  if( za == 0 )
    {
      if( rp < r0 )
	{ return 0; }
      else
	{ return no_crossing; }
    }
  if( fabs(za) == 180 )
    {
      if( rp > r0 )
	{ return 0; }
      else
	{ return no_crossing; }
    }

  // Check if the given LOS goes in the direction towards the pressure surface.
  // If not, return 999.
  //
  const bool downwards = is_los_downwards_2d( za, psurface_tilt_2d( r0, c1 ) );
  //
  if( ( rp < r0  &&  downwards )  ||  ( rp >= r0  &&  !downwards ) )
    { return no_crossing; }


  // The case with c1=c2=0 can be handled analytically
  if( c1 == 0  &&  c2 == 0 )
    {
      return geompath_lat_at_za( za, 0, 
		       geompath_za_at_r( geometrical_ppc( rp, za ), za, r0 ) );
    }


  // Approximative solution:
  else
    {

      // The nadir angle in radians, and cosine and sine of that angle
      const Numeric   beta = DEG2RAD * ( 180 - fabs(za) );
      const Numeric   cv = cos( beta );
      const Numeric   sv = sin( beta );

      // Convert slope to m/radian and consider viewing direction
      c1 *= RAD2DEG;
      c2 *= RAD2DEG;
      if( za < 0 )
        { 
          c1 = -c1; 
          c2 = -c2; 
        }
      
      // The case when c2=0 must be treated seperately as the root solving 
      // function does not accept that the highest polynomial has 
      // coefficient = 0.
      Index n = 6;
      if( c2 == 0 )
        { n = 5; }
          
      // The vector of polynomial coefficients
      Vector p(n);
      //
      p[0] = ( r0 - rp ) * sv;
      p[1] = r0 * cv + c1 * sv;
      p[2] = -r0 * sv / 2 + c1 * cv + c2 * sv;
      p[3] = -r0 * cv / 6 - c1 * sv / 2 + c2 * cv;
      p[4] = -c1 * cv / 6 - c2 * sv / 2;
      //
      if( n == 6 )
        { p[5] = -c2 * cv / 6; }

      // Calculate roots of the polynomial
      Matrix roots(n-1,2);
      poly_root_solve( roots, p );

      //MatrixPrint(roots,"roots");

      // Find first root with imaginary part = 0, and real part > 0 or >= 0.
      // The solution dlat = 0 is not of interest if rp = r1.
      Numeric dlat = no_crossing;
      for( Index i=0; i<n-1; i++ )
        {
          if( roots(i,1) == 0  &&   roots(i,0) >= 0 )
    	{
    	  if( ( rp != r0 || roots(i,0) > 0 )  &&  RAD2DEG * roots(i,0) < dlat )
    	    { dlat = RAD2DEG * roots(i,0); }
    	}
        }  

      // Change sign if zenith angle is negative
      if( dlat < no_crossing  &&  za < 0 )
        { dlat = -dlat; }

      return dlat;
    }
}




/*===========================================================================
  === Functions operating on the Ppath structure
  ===========================================================================*/

//! ppath_init_structure
/*!
   Initiates a Ppath structure to hold the given number of points.

   All fields releated with the ground, symmetry and tangent point are set
   to 0 or empty. The background field is set to background case 0. The
   constant field is set to -1. The refraction field is set to 0.

   The length of the l_step field is set to np-1.

   \param   ppath            Output: A Ppath structure.
   \param   atmosphere_dim   The atmospheric dimensionality.
   \param   np               Number of points of the path.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_init_structure( 
	      Ppath&      ppath,
	const Index&      atmosphere_dim,
        const Index&      np )
{
  assert( atmosphere_dim >= 1 );
  assert( atmosphere_dim <= 3 );

  ppath.dim        = atmosphere_dim;
  ppath.np         = np;
  ppath.refraction = 0;
  ppath.method     = "-";
  ppath.constant   = -1;   
  if( atmosphere_dim < 3 )
    {
      ppath.pos.resize( np, 2 );
      ppath.los.resize( np, 1 );
    }
  else
    {
      ppath.pos.resize( np, atmosphere_dim );
      ppath.los.resize( np, 2 );
      ppath.gp_lon.resize( np );
    }
  ppath.gp_p.resize( np );
  if( atmosphere_dim >= 2 )
    { ppath.gp_lat.resize( np ); }
  ppath.z.resize( np );
  if( np > 0 )
    { ppath.l_step.resize( np-1 ); }
  else
    { ppath.l_step.resize( 0 ); }
  ppath_set_background( ppath, 0 );
  ppath.tan_pos.resize(0);
  ppath.geom_tan_pos.resize(0);
}



//! ppath_set_background 
/*!
   Sets the background field of a Ppath structure.

   The different background cases have a number coding to simplify a possible
   change of the strings and checking of the what case that is valid.

   The case numbers are:                    <br>
      0. Not yet set.                       <br>
      1. Space.                             <br>
      2. The ground.                        <br>
      3. The surface of the cloud box.      <br>
      4. The interior of the cloud box.     

   \param   ppath            Output: A Ppath structure.
   \param   case_nr          Case number (see above)

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_set_background( 
	      Ppath&      ppath,
        const Index&      case_nr )
{
  switch ( case_nr )
    {
    case 0:
      ppath.background = "";
      break;
    case 1:
      ppath.background = "space";
      break;
    case 2:
      ppath.background = "ground";
      break;
    case 3:
      ppath.background = "cloud box surface";
      break;
    case 4:
      ppath.background = "cloud box interior";
      break;
    default:
      ostringstream os;
      os << "Case number " << case_nr << " is not defined.";
      throw runtime_error(os.str());
    }
}



//! ppath_what_background
/*!
   Returns the case number for the radiative background.

   See further the function *ppath_set_background*.

   \return                   The case number.
   \param   ppath            A Ppath structure.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Index ppath_what_background( const Ppath&   ppath )
{
  if( ppath.background == "" )
    { return 0; }
  else if( ppath.background == "space" )
    { return 1; }
  else if( ppath.background == "ground" )
    { return 2; }
  else if( ppath.background == "cloud box surface" )
    { return 3; }
  else if( ppath.background == "cloud box interior" )
    { return 4; }
  else
    {
      ostringstream os;
      os << "The string " << ppath.background 
	 << " is not a valid background case.";
      throw runtime_error(os.str());
    }
}



//! ppath_copy
/*!
   Copy the content in ppath2 to ppath1.

   The ppath1 structure must be allocated before calling the function. The
   structure can be allocated to hold more points than found in ppath2.
   The data of ppath2 is placed in the first positions of ppath1.

   \param   ppath1    Output: PPath structure.
   \param   ppath2    The Ppath structure to be copied.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_copy(
	   Ppath&      ppath1,
     const Ppath&      ppath2 )
{
  assert( ppath1.np >= ppath2.np ); 

  // The field np shall not be copied !!!

  ppath1.dim        = ppath2.dim;
  ppath1.refraction = ppath2.refraction;
  ppath1.method     = ppath2.method;
  ppath1.constant   = ppath2.constant;
  ppath1.background = ppath2.background;

  for( Index i=0; i<ppath2.np; i++ )
    {
      ppath1.pos(i,0)      = ppath2.pos(i,0);
      ppath1.pos(i,1)      = ppath2.pos(i,1);
      ppath1.los(i,0)      = ppath2.los(i,0);
      ppath1.z[i]          = ppath2.z[i];
      gridpos_copy( ppath1.gp_p[i], ppath2.gp_p[i] );
      
      if( ppath1.dim >= 2 )
	{ gridpos_copy( ppath1.gp_lat[i], ppath2.gp_lat[i] ); }
      
      if( ppath1.dim == 3 )
	{
	  ppath1.pos(i,2)        = ppath2.pos(i,2);
	  ppath1.los(i,1)        = ppath2.los(i,1);
	  gridpos_copy( ppath1.gp_lon[i], ppath2.gp_lon[i] ); 
	}
      
      if( i > 0 )
	{ ppath1.l_step[i-1] = ppath2.l_step[i-1]; }
     }

  if( ppath2.tan_pos.nelem() )
    {
      ppath1.tan_pos.resize( ppath2.tan_pos.nelem() );
      ppath1.tan_pos = ppath2.tan_pos; 
    }
  if( ppath2.geom_tan_pos.nelem() )
    {
      ppath1.geom_tan_pos.resize( ppath2.geom_tan_pos.nelem() );
      ppath1.geom_tan_pos = ppath2.geom_tan_pos; 
    }
}



//! ppath_append
/*!
   Combines two Ppath structures   

   The function appends a Ppath structure to another structure. 
 
   All the data of ppath1 is kept.

   The first point in ppath2 is assumed to be the same as the last in ppath1.
   Only data in ppath from the fields pos, los, z, l_step, gp_xxx and 
   background are considered.

   \param   ppath1    Output: Ppath structure to be expanded.
   \param   ppath2    The Ppath structure to include in ppath.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_append(
	   Ppath&   ppath1,
     const Ppath&   ppath2 )
{
  const Index n1 = ppath1.np;
  const Index n2 = ppath2.np;

  Ppath   ppath;
  ppath_init_structure( ppath, ppath1.dim, n1 );
  ppath_copy( ppath, ppath1 );

  ppath_init_structure( ppath1, ppath1.dim, n1 + n2 - 1 );
  ppath_copy( ppath1, ppath );

  // Append data from ppath2
  Index i1;
  for( Index i=1; i<n2; i++ )
    {
      i1 = n1 + i - 1;

      ppath1.pos(i1,0)      = ppath2.pos(i,0);
      ppath1.pos(i1,1)      = ppath2.pos(i,1);
      ppath1.los(i1,0)      = ppath2.los(i,0);
      ppath1.z[i1]          = ppath2.z[i];
      gridpos_copy( ppath1.gp_p[i1], ppath2.gp_p[i] );

      if( ppath1.dim >= 2 )
	{ gridpos_copy( ppath1.gp_lat[i1], ppath2.gp_lat[i] ); }
      
      if( ppath1.dim == 3 )
	{
	  ppath1.pos(i1,2)        = ppath2.pos(i,2);
	  ppath1.los(i1,1)        = ppath2.los(i,1);
	  gridpos_copy( ppath1.gp_lon[i1], ppath2.gp_lon[i] ); 
	}
      
      ppath1.l_step[i1-1] = ppath2.l_step[i-1];
    }

  if( ppath_what_background( ppath2 ) )
    { ppath1.background = ppath2.background; }
}



//! ppath_fill_1d
/*!
   Fills a 1D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step and gp_p.
   The pressure grid positions (gp_p) are filtered through gridpos_check_fd.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath      Output: Ppath structure.
   \param   r          Vector with radius for the path points.
   \param   lat        Vector with latitude for the path points.
   \param   za         Vector with zenith angle for the path points.
   \param   lstep      Length along the path between the points.
   \param   r_geoid    Geoid radii.
   \param   z_field    Geometrical altitudes.
   \param   ip         Pressure grid range.

   \author Patrick Eriksson
   \date   2002-07-18
*/
void ppath_fill_1d(
	   Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   za,
     ConstVectorView   lstep,
     const Numeric&    r_geoid,
     ConstVectorView   z_field,
     const Index&      ip )
{
  // Help variables that are common for all points.
  const Numeric   r1 = r_geoid + z_field[ip];
  const Numeric   dr = z_field[ip+1] - z_field[ip];

  for( Index i=0; i<r.nelem(); i++ )
    {
      ppath.pos(i,0) = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.los(i,0) = za[i];
      
      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - r1 ) / dr;
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      ppath.z[i] = r[i] - r_geoid;

      if( i > 0 )
	{ ppath.l_step[i-1] = lstep[i-1]; }
    }
}



//! ppath_fill_2d
/*!
   Fills a 2D Ppath structure with position and LOS values.

   The function fills the fields: pos, los, z, l_step, gp_p and gp_lat.

   The structure fields must be allocated to correct size before calling the 
   function. The field size must be at least as large as the length of r,
   lat and za vectors.

   The length along the path shall be the same between all points.

   \param   ppath      Output: Ppath structure.
   \param   r          Vector with radius for the path points.
   \param   lat        Vector with latitude for the path points.
   \param   za         Vector with zenith angle for the path points.
   \param   lstep      Length along the path between the points.
   \param   r_geoid    Geoid radii.
   \param   z_field    Geometrical altitudes
   \param   lat_grid   Latitude grid.
   \param   ip         Pressure grid range.
   \param   ilat       Latitude grid range.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_fill_2d(
	   Ppath&      ppath,
     ConstVectorView   r,
     ConstVectorView   lat,
     ConstVectorView   za,
     const Numeric&    lstep,
     ConstVectorView   r_geoid,
     ConstMatrixView   z_field,
     ConstVectorView   lat_grid,
     const Index&      ip,
     const Index&      ilat )
{
  // Help variables that are common for all points.
  const Numeric   dlat  = lat_grid[ilat+1] - lat_grid[ilat];
  const Numeric   r1low = r_geoid[ilat] + z_field(ip,ilat);
  const Numeric   drlow = r_geoid[ilat+1] + z_field(ip,ilat+1) -r1low;
  const Numeric   r1upp = r_geoid[ilat] + z_field(ip+1,ilat);
  const Numeric   drupp = r_geoid[ilat+1] + z_field(ip+1,ilat+1) - r1upp;
  const Numeric   z1low =  z_field(ip,ilat);
  const Numeric   dzlow =  z_field(ip,ilat+1) -z1low;
  const Numeric   z1upp =  z_field(ip+1,ilat);
  const Numeric   dzupp =  z_field(ip+1,ilat+1) - z1upp;

  for( Index i=0; i<r.nelem(); i++ )
    {
      ppath.pos(i,0) = r[i];
      ppath.pos(i,1) = lat[i];
      ppath.los(i,0) = za[i];
      
      // Weigt in the latitude direction
      Numeric w = ( lat[i] - lat_grid[ilat] ) / dlat;

      // Radius of lower and upper face at present latitude
      const Numeric rlow = r1low + w * drlow;
      const Numeric rupp = r1upp + w * drupp;

      // Geometrical altitude of lower and upper face at present latitude
      const Numeric zlow = z1low + w * dzlow;
      const Numeric zupp = z1upp + w * dzupp;

      ppath.gp_p[i].idx   = ip;
      ppath.gp_p[i].fd[0] = ( r[i] - rlow ) / ( rupp - rlow );
      ppath.gp_p[i].fd[1] = 1 - ppath.gp_p[i].fd[0];
      gridpos_check_fd( ppath.gp_p[i] );

      ppath.z[i] = zlow + ppath.gp_p[i].fd[0] * ( zupp -zlow );

      ppath.gp_lat[i].idx   = ilat;
      ppath.gp_lat[i].fd[0] = ( lat[i] - lat_grid[ilat] ) / dlat;
      ppath.gp_lat[i].fd[1] = 1 - ppath.gp_lat[i].fd[0];
      gridpos_check_fd( ppath.gp_lat[i] );

      if( i > 0 )
	{ ppath.l_step[i-1] = lstep; }
    }
}





/*===========================================================================
  === A help functions to ppathCalc
  ===========================================================================*/

//! ppath_start_stepping
/*!
   Initiates a Ppath structure for calculation of a path with *ppath_step*.

   The function performs two main tasks. As mentioned above, it initiates
   a Ppath structure (a), but it also checks that the end point of the path is
   at an allowed location (b).

   (a): The Ppath structure is set to hold the position and LOS of the last
   point of the path inside the atmosphere. This point is either the
   sensor position, or the point where the path leaves the model atmosphere.
   If the path is totally outside the atmosphere, no point is put into the
   structure. If the (practical) end and start points are identical, such
   as when the sensor is inside the cloud box, the background field is set.

   (b): If it is found that the end point of the path is at an illegal position
   a detailed error message is given. Not allowed cases are: <br>  
      1. The sensor is placed below ground level. <br> 
      2. For 2D and 3D, the path leaves the model atmosphere at a latitude or
         longitude end face. <br> 
      3. For 2D and 3D, the path is totally outside the atmosphere and the 
         latitude and longitude of the tangent point is outside the range of
         the corresponding grids. 

   All input variables are identical with the WSV with the same name.
   The output variable is here called ppath for simplicity, but is in
   fact *ppath_step*.

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   r_geoid           The geoid radius.
   \param   z_ground          Ground altitude.
   \param   cloudbox_on       Flag to activate the cloud box.
   \param   cloudbox_limits   Index limits of the cloud box.
   \param   a_pos             The position of the sensor.
   \param   a_los             The line-of-sight of the sensor.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void ppath_start_stepping(
              Ppath&          ppath,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        ConstTensor3View      z_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        ConstVectorView       a_pos,
        ConstVectorView       a_los )
{
  // This function contains no checks or asserts as it is only a sub-function
  // to ppathCalc where the input is checked carefully.

  // Allocate the ppath structure
  ppath_init_structure(  ppath, atmosphere_dim, 1 );

  // Number of pressure levels
  const Index np = p_grid.nelem();

  // The different atmospheric dimensionalities are handled seperately

  //-- 1D ---------------------------------------------------------------------
  if( atmosphere_dim == 1 )
    {
      // Radius for the ground
      const Numeric r_ground = r_geoid(0,0) + z_ground(0,0);

      // Radius for the top of the atmosphere
      const Numeric r_top = r_geoid(0,0) + z_field(np-1,0,0);

      // The only forbidden case here is that the sensor is below the ground
      if( a_pos[0] < r_ground )
	{
          ostringstream os;
          os << "The sensor is placed " 
             << (r_ground - a_pos[0])/1e3 << " km below ground level.\n"
             << "The sensor must be above the ground.";
	  throw runtime_error(os.str());
	}

      out2 << "  sensor altitude        : " << (a_pos[0]-r_geoid(0,0))/1e3 
	   << " km\n";

      // If downwards, calculate geometrical tangent position
      Vector geom_tan_pos(0);
      if( a_los[0] >= 90 )
	{
	  geom_tan_pos.resize(2);
	  geom_tan_pos[0] = geometrical_ppc( a_pos[0], a_los[0] );
          geom_tan_pos[1] = geompath_lat_at_za( a_los[0], 0, 90 );
	  out2 << "  geom. tangent radius   : " << geom_tan_pos[0]/1e3 
	       <<" km\n";
	  out2 << "  geom. tangent latitude : " << geom_tan_pos[1] << "\n";
	  out2 << "  geom. tangent altitude : " 
	       << (geom_tan_pos[0]-r_geoid(0,0))/1e3 << " km\n";
	}

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = a_pos[0];
      ppath.los(0,0) = a_los[0];
      ppath.pos(0,1) = 0; 
      ppath.z[0]     = ppath.pos(0,0) - r_geoid(0,0);
      

      // The sensor is inside the model atmosphere, 1D ------------------------
      if( a_pos[0] < r_top )
	{
	  // Use below the values in ppath (instead of a_pos and a_los) as 
	  // they can be modified on the way.
     
	  // Is the sensor on the ground looking down?
	  // If yes and the sensor is inside the cloudbox, the background will
	  // be changed below.
	  if( ppath.pos(0,0) == r_ground  &&  ppath.los(0,0) > 90 )
	    { ppath_set_background( ppath, 2 ); }

	  // Check sensor position with respect to cloud box.
	  if( cloudbox_on )
	    {
	      // Is the sensor inside the cloud box?
	      if( ppath.z[0] > z_field(cloudbox_limits[0],0,0)  && 
		                 ppath.z[0] < z_field(cloudbox_limits[1],0,0) )
		{ ppath_set_background( ppath, 4 ); }

	      else if( ( ppath.z[0] == z_field(cloudbox_limits[0],0,0)  && 
                                                        ppath.los(0,0) <= 90 )
                     || 
		  ( ppath.z[0] == z_field(cloudbox_limits[1],0,0)  && 
                                                        ppath.los(0,0) > 90 ) )
		{
		  ppath_set_background( ppath, 3 );
		}
	    }
	}

      // The sensor is outside the model atmosphere, 1D -----------------------
      else
	{
	  // Upward observations are not allowed here
	  if( a_los[0] <= 90 )
	      throw runtime_error("When the sensor is placed outside the model"
                         " atmosphere, upward observations are not allowed." );

	  // We can here set the path constant, that equals the radius of the
	  // geometrical tangent point.
	  ppath.constant = geom_tan_pos[0];
 
	  // Path is above the atmosphere
	  if( ppath.constant >= r_top )
	    {
	      ppath_set_background( ppath, 1 );
	      out1 << "  --- WARNING ---, path is totally outside of the "
		   << "model atmosphere\n";
	    }

	  // Path enters the atmosphere
	  else
	    {
	      ppath.z[0]     = z_field(np-1,0,0);
              ppath.pos(0,0) = r_top;
	      ppath.los(0,0) = geompath_za_at_r( ppath.constant, a_los[0], 
                                                                       r_top );
	      ppath.pos(0,1) = geompath_lat_at_za( a_los[0], 0, 
                                                              ppath.los(0,0) );
	    }
	}

      // Get grid position for the end point, if it is inside the atmosphere.
      if( ppath.z[0] <= z_field(np-1,0,0) )
	{ gridpos( ppath.gp_p, z_field(Range(joker),0,0), ppath.z ); }

      // Set geometrical tangent point position
      if( geom_tan_pos.nelem() == 2 )
	{
	  ppath.geom_tan_pos.resize(2);
	  ppath.geom_tan_pos = geom_tan_pos;
	}

    }  // End 1D


  //-- 2D ---------------------------------------------------------------------
  else if( atmosphere_dim == 2 )
    {
      // Number of points in the pressure and latitude grids
      const Index   np = p_grid.nelem();
      const Index   nlat = lat_grid.nelem();

      // Is the sensor inside the latitude range of the model atmosphere,
      // and below the top of the atmosphere? If yes, is_inside = 1.
      // Store geoid and ground radii, grid position and interpolation weights
      // for later use.
      //
      Index   is_inside = 0;   
      Numeric rv_geoid=-1, rv_ground=-1;  // -1 to avoid compiler warnings
      ArrayOfGridPos gp_lat(1);
      Matrix itw(1,2);
      //
      if( a_pos[1] >= lat_grid[0]  &&  a_pos[1] <= lat_grid[nlat-1] )
	{
	  gridpos( gp_lat, lat_grid, Vector(1,a_pos[1]) );
	  interpweights( itw, gp_lat );

	  Vector v_rgeoid(1), v_zground(1), v_ztop(1);
	  interp( v_rgeoid, itw, r_geoid(Range(joker),0), gp_lat );
	  interp( v_zground, itw, z_ground(Range(joker),0), gp_lat );
	  interp( v_ztop, itw, z_field(np-1,Range(joker),0), gp_lat );

	  rv_geoid  = v_rgeoid[0];
	  rv_ground = rv_geoid + v_zground[0];

	  out2 << "  sensor altitude        : " << (a_pos[0]-rv_geoid)/1e3 
	       << " km\n";

	  if( a_pos[0] < ( v_rgeoid[0] + v_ztop[0] ) )
	    { is_inside = 1; }
	}

      // If downwards, calculate geometrical tangent position. If the tangent
      // point is inside the covered latitude range, calculate also the 
      // geometrical altitude of the tangent point and the top of atmosphere.
      //
      Vector  geom_tan_pos(0);
      Numeric geom_tan_z=-1, geom_tan_atmtop=-1;  // -1 to avoid warnings
      //
      if( fabs(a_los[0]) >= 90 )
	{
	  geom_tan_pos.resize(2);
	  geom_tan_pos[0] = geometrical_ppc( a_pos[0], a_los[0] );
	  if( a_los[0] > 0 )
	    { geom_tan_pos[1] = geompath_lat_at_za( a_los[0], a_pos[1], 90 ); }
	  else
	    { geom_tan_pos[1] = geompath_lat_at_za( a_los[0], a_pos[1], -90 );}
	  out2 << "  geom. tangent radius         : " << geom_tan_pos[0] / 1e3
	       <<" km\n";
	  out2 << "  geom. tangent latitude       : " << geom_tan_pos[1] 
	       << "\n";
	  //
	  if( geom_tan_pos[1] >= lat_grid[0]  &&  
	                                  geom_tan_pos[1] <= lat_grid[nlat-1] )
	    {
	      ArrayOfGridPos gp_tmp(1);
	      Matrix itw_tmp(1,2);
	      Vector v1_tmp(1), v2_tmp(1);
	      gridpos( gp_tmp, lat_grid, Vector(1,geom_tan_pos[1]) );
	      interpweights( itw_tmp, gp_tmp );
	      interp( v1_tmp, itw_tmp, r_geoid(Range(joker),0), gp_tmp );
	      geom_tan_z = geom_tan_pos[0] - v1_tmp[0];
	      interp( v2_tmp, itw_tmp, z_field(np-1,Range(joker),0), gp_tmp );
	      geom_tan_atmtop = v2_tmp[0];
	      out2 << "  geom. tangent altitude       : " << geom_tan_z/1e3 
		   << " km\n";
	    }
	}

      // Put sensor position and LOS in ppath as first guess
      ppath.pos(0,0) = a_pos[0];
      ppath.pos(0,1) = a_pos[1];
      ppath.los(0,0) = a_los[0];

      // The sensor is inside the model atmosphere, 2D ------------------------
      if( is_inside )
	{
	  // Check that the sensor is above the ground
	  if( a_pos[0] < rv_ground )
	    {
	      ostringstream os;
	      os << "The sensor is placed " 
		 << (rv_ground - a_pos[0])/1e3 << " km below ground level.\n"
		 << "The sensor must be above the ground.";
	      throw runtime_error(os.str());
	    }

	  // Check that not at latitude end point and looks out
	  if( ( a_pos[1] == lat_grid[0]  &&   a_los[0] < 0 ) )
	    throw runtime_error( "The sensor is at the lower latitude end "
                                            "point and the zenith angle < 0" );
	  if( a_pos[1] == lat_grid[nlat-1]  &&   a_los[0] > 0 ) 
	    throw runtime_error( "The sensor is at the upper latitude end "
                                            "point and the zenith angle > 0" );
	  
	  // Geometrical altitude
	  ppath.z[0] = ppath.pos(0,0) - rv_geoid;

	  // Use below the values in ppath (instead of a_pos and a_los) as 
	  // they can be modified on the way.
     
	  // Grid positions
	  ppath.gp_lat[0].idx   = gp_lat[0].idx;
	  ppath.gp_lat[0].fd[0] = gp_lat[0].fd[0];
	  ppath.gp_lat[0].fd[1] = gp_lat[0].fd[1];
	  // Create a vector with the geometrical altitude of the pressure 
	  // surfaces for the sensor latitude and use it to get ppath.gp_p.
	  Vector z_grid(np);
          z_at_lat_2d( z_grid, p_grid, lat_grid, z_field, gp_lat );
	  gridpos( ppath.gp_p, z_grid, ppath.z );

	  // Is the sensor on the ground looking down?
	  if( ppath.pos(0,0) == rv_ground )
	    {
	      // Calculate radial slope of the ground
	      const Numeric rslope = psurface_slope_2d( lat_grid, 
                        r_geoid(Range(joker),0), z_ground(Range(joker),0), 
                                              gp_lat[0], ppath.los(0,0) >= 0 );

	      // Calculate angular tilt of the ground
	      const Numeric atilt = psurface_tilt_2d( rv_ground, rslope);

	      // Are we looking down into the ground?
	      // If yes and the sensor is inside the cloudbox, the background 
	      // will be changed below.
	      if( is_los_downwards_2d( ppath.los(0,0), atilt ) )
		{ ppath_set_background( ppath, 2 ); }
	    }

	  // Check sensor position with respect to cloud box.
	  if( cloudbox_on )
	    {
	      // To check all possible cases here when the sensor is at the
	      // surface and can either look into or out from the box needs
	      // a lot of coding.
	      // So we are instead sloppy and set all cases when the sensor
	      // is inside or at the surface to be inside the box.
	      // The neglected cases should be very unlikely in for real
	      // situations.

	      if( ppath.pos(0,1) >= lat_grid[cloudbox_limits[2]]  &&
		               ppath.pos(0,1) <= lat_grid[cloudbox_limits[3]] )
		{
		  // Calculate the lower and upper altitude radius limit for
		  // the cloud box at the latitude of the sensor
	          Vector v_zlim(1);
		  interp( v_zlim, itw, 
                          z_field(cloudbox_limits[0],Range(joker),0), gp_lat );
		  Numeric rv_low = rv_geoid + v_zlim[0];
		  interp( v_zlim, itw, 
                          z_field(cloudbox_limits[1],Range(joker),0), gp_lat );
		  Numeric rv_upp = rv_geoid + v_zlim[0];

		  if( ppath.pos(0,0) >= rv_low  &&  ppath.pos(0,0) <= rv_upp )
		    { ppath_set_background( ppath, 4 ); }	
		}
	    }
	}

      // The sensor is outside the model atmosphere, 2D -----------------------
      else
	{
	  // Upward observations are not allowed here
	  if( fabs(a_los[0]) <= 90 )
	    {
	      ostringstream os;
	      os << "When the sensor is placed outside the model atmosphere,\n"
		 << "upward observations are not allowed.";
	      throw runtime_error( os.str() );
	    }
	  
	  // We can here set the path constant, that equals the radius of the
	  // geometrical tangent point.
	  ppath.constant = geom_tan_pos[0];

	  // Handle cases when the sensor appears to look the wrong way
	  if( ( a_pos[1] <= lat_grid[0]  &&  a_los[0] <= 0 )  || 
	                  ( a_pos[1] >= lat_grid[nlat-1]  &&  a_los[0] >= 0 ) )
	    {
	      ostringstream os;
	      os << "The sensor is outside (or at the limit) of the model "
		 << "atmosphere but\nlooks in the wrong direction (wrong sign "
		 << "for the zenith angle?).\nThis case includes nadir "
		 << "looking exactly at the latitude end points.";
	      throw runtime_error( os.str() );
	    }

	  // If the sensor is outside the latitude range, check that path is
	  // above the closest corner of the model atmosphere
	  if( a_pos[1] < lat_grid[0]  ||  a_pos[1] > lat_grid[nlat-1] )
	    {
	      Index   ic = 0;
	      String  sc = "lower";
	      if( a_pos[1] > lat_grid[0] )
		{ ic = nlat - 1;   sc = "upper"; }
	      const Numeric rv = geompath_r_at_lat( ppath.constant, a_pos[1], 
                                                      a_los[0], lat_grid[ic] );
	      if( rv < ( r_geoid(ic,0) + z_field(np-1,ic,0) ) )
		{
		  ostringstream os;
		  os << "The sensor is outside of the model atmosphere and "
		     << "looks in the\n" << sc << " latitude end face.\n"
		     << "The geometrical altitude of the corner point is "
		     << z_field(np-1,ic,0)/1e3 << " km.\n"
		     << "The geometrical altitude of the entrance point is "
		     << (rv-r_geoid(ic,0))/1e3 << " km.";
		  throw runtime_error( os.str() );
		}
	    }

	  // If the tangent point is inside covered latitude range, everything
	  // is OK. If not, the path must be below the corner of the model
	  // atmosphere. 
	  if( ( geom_tan_pos[1] < lat_grid[0]  ||  
		                         geom_tan_pos[1] > lat_grid[nlat-1] ) )
	    {
	      Index   ic = 0;
	      String  sc = "lower";
	      if( a_los[0] >= 0 )
		{ ic = nlat - 1;   sc = "upper"; }
	      const Numeric rv = geompath_r_at_lat( ppath.constant, a_pos[1], 
                                                      a_los[0], lat_grid[ic] );
	      if( rv >= ( r_geoid(ic,0) + z_field(np-1,ic,0) ) )
		{
		  ostringstream os;
		  os << "The combination of sensor position and line-of-sight "
		     << "gives a\npropagation path that goes above the model "
		     << "atmosphere, with\na tangent point outside the covered"
		     << " latitude range.\nThe latitude of the tangent point "
		     << "is " << geom_tan_pos[1] << " degrees.";
		  throw runtime_error( os.str() );
		}
	    }

	  // That should be all needed checks. We know now that the path is
	  // either totally outside the atmosphere, with a tangent point 
	  // inside lat_grid, or enters the atmosphere from the top 
	  // somewhere inside lat_grid. In the latter case we need to
	  // determine the latitude of the entrance point.
	  
	  // Path is above the atmosphere:
	  // Requieres that tangent point is inside lat_grid and above the
	  // top of the atmosphere.
	  if( geom_tan_pos[1] >= lat_grid[0]  &&  
                           geom_tan_pos[1] <= lat_grid[nlat-1]   &&  
                                                geom_tan_z >= geom_tan_atmtop )
	    {
	      ppath_set_background( ppath, 1 );
	      out1 << "  --- WARNING ---: path is totally outside of the "
		   << "model atmosphere\n";
	    }

	  // The path enters the atmosphere
	  else
	    {
	      // Find the latitude where the path passes top of the atmosphere.

	      // We are handling this in a rather dumb way. A test is performed
	      // for each latitude range using psurface_crossing_2d.
	      // A bit smarter algorithm was considered but that made the code 
	      // more messy.
	      // The case when the sensor is placed inside lat_grid must be
	      // hanled seperetaly.

	      // Determine first latitude range of interest, search direction
	      // and first test latitude.
	      Numeric lat0;
	      Index   ilat0, istep;
	      //
	      if( a_pos[1] <= lat_grid[0] )
		{
		  lat0  = lat_grid[0]; 
		  ilat0 = 0;
		  istep = 1;
		}
	      else if( a_pos[1] >= lat_grid[nlat-1] )
		{
		  lat0  = lat_grid[nlat-1]; 
		  ilat0 = nlat-1;
		  istep = -1;
		}
	      else
		{
		  lat0  = a_pos[1]; 
		  if( a_los[0] >= 0 )
		    {  
		      ilat0 = gridpos2gridrange( gp_lat[0], 1 );
		      istep = 1;
		    }
		  else
		    { 
		      ilat0 = gridpos2gridrange( gp_lat[0], 0 ) + 1;
		      istep = -1; 
		    }
		}

	      // Loop until entrance point is found
	      Index ready = 0;
	      while( !ready )
		{
		  // Calculate radius and zenith angle of path at lat0
		  Numeric r0  = geompath_r_at_lat( ppath.constant, a_pos[1], 
                                                              a_los[0], lat0 );
		  Numeric za0 = geompath_za_at_r( ppath.constant, a_los[0], 
                                                                          r0 );

		  // Calculate radius and slope to use in psurface_crossing_2d
		  Numeric rv1 = r_geoid(ilat0,0) + z_field(np-1,ilat0,0);
		  Numeric rv2 = r_geoid(ilat0+istep,0) + 
		                                   z_field(np-1,ilat0+istep,0);
		  Numeric latstep = lat_grid[ilat0+istep] - lat_grid[ilat0];
		  Numeric c = istep * ( rv2 - rv1 ) / latstep;

		  if( lat0 != lat_grid[ilat0] )
		    { rv1 = rv1 + c * ( lat0 - lat_grid[ilat0] ); } 

		  Numeric dlat = psurface_crossing_2d( r0, za0, rv1, c, 0 );

		  if( fabs(dlat) <= fabs(latstep) )
		    {
		      ready = 1;
		      ppath.pos(0,1) = lat0 + dlat;
		      ppath.pos(0,0) = rv1 + c * dlat;
		      ppath.los(0,0) = geompath_za_at_r( ppath.constant, 
                                                    a_los[0], ppath.pos(0,0) );
		      // Re-use some variables from above
		      rv1 = r_geoid(ilat0,0);
		      rv2 = r_geoid(ilat0+istep,0);
		      c   = ( rv2 - rv1 ) / latstep;
		      ppath.z[0] = ppath.pos(0,0) - ( rv1 + istep * c *
                                        ( ppath.pos(0,1) - lat_grid[ilat0] ) );
		      ppath.gp_p[0].idx = np - 2;
		      ppath.gp_p[0].fd[0] = 1;
		      ppath.gp_p[0].fd[1] = 0;
		      gridpos( ppath.gp_lat, lat_grid, 
                                                   ppath.pos(Range(joker),1) );
		    } 
		  else
		    {
		      ilat0 += istep;
		      lat0   = lat_grid[ilat0];
		    }
		} 
	    }
	}      

      // Set geometrical tangent point position
      if( geom_tan_pos.nelem() == 2 )
	{
	  ppath.geom_tan_pos.resize(2);
	  ppath.geom_tan_pos = geom_tan_pos;
	}

    }  // End 2D


  //-- 3D ---------------------------------------------------------------------
  else if( atmosphere_dim == 3 )
    {
      throw runtime_error("The function handles not yet 3D.");
    }  // End 3D
}





/*===========================================================================
  === Help functions for the *ppath_step* functions found below
  === These functions are mainly pieces of code that are common for at least
  === two functions (or two places in some function) and for this reason 
  === there is not much documentation. 
  ===========================================================================*/

// This function is copied from arts-1 as a temporary solution.

//// refr_index_BoudourisDryAir ///////////////////////////////////////////////
/**
   Calculates the refractive index for dry air at microwave frequncies 
   following Boudouris 1963.

   The expression is also found in Chapter 5 of the Janssen book.

   The atmosphere is assumed to have no water vapour.

   \retval   refr_index  refractive index
   \param    p_abs       absorption pressure grid
   \param    t_abs       temperatures at p_abs

   \author Patrick Eriksson
   \date   2001-02-16
*/
void refr_index_BoudourisDryAir (
             Vector&     refr_index,
       ConstVectorView   p_abs,
       ConstVectorView   t_abs )
{
  const Index   n = p_abs.nelem();
  refr_index.resize( n );

  assert ( n == t_abs.nelem() );

  // N = 77.593e-2 * p / t ppm
  for ( Index i=0; i<n; i++ )
    refr_index[i] = 1.0 + 77.593e-8 * p_abs[i] / t_abs[i];
}



//! get_refr_index_1d
/*! 
   A temporary function to get refractive index for 1D cases.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-13
*/
Numeric get_refr_index_1d(
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        const Numeric&    z )
{      
  ArrayOfGridPos gp(1);
  gridpos( gp, z_field, Vector(1,z) );

  Matrix itw(1,2);
  interpweights( itw, gp );

  Vector p_value(1), t_value(1);

  itw2p( p_value, p_grid, gp, itw );
  interp( t_value, itw, t_field, gp );

  Vector refr_index(1);

  refr_index_BoudourisDryAir( refr_index, p_value[0], t_value[0] );

  //return 1;
  return refr_index[0];
}



//! get_refr_index_2d
/*! 
   A temporary function to get refractive index for 2D cases.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-18
*/
Numeric get_refr_index_2d(
        ConstVectorView    p_grid,
        ConstVectorView    lat_grid,
        ConstVectorView    r_geoid,
        ConstMatrixView    z_field,
        ConstMatrixView    t_field,
        const Numeric&     z,
        const Numeric&     lat )
{     
  const Index      np = p_grid.nelem();
  ArrayOfGridPos   gp_p(1), gp_lat(1);
  Vector           z_grid(np), r_value(1), p_value(1), t_value(1);
  Matrix           itw(1,2);

  gridpos( gp_lat, lat_grid, Vector(1,lat) );
  interpweights( itw, gp_lat );
  interp( r_value, itw, r_geoid, gp_lat );

  z_at_lat_2d( z_grid, p_grid, lat_grid, z_field, gp_lat );

  gridpos( gp_p, z_grid, Vector(1,z-r_value[0]) );
  interpweights( itw, gp_p );
  itw2p( p_value, p_grid, gp_p, itw );

  itw.resize(1,4);
  interpweights( itw, gp_p, gp_lat );
  interp( t_value, itw, t_field, gp_p, gp_lat );

  Vector refr_index(1);

  refr_index_BoudourisDryAir( refr_index, p_value[0], t_value[0] );

  return refr_index[0];

}



//! ppath_start_1d
/*! 
   Internal help function for 1D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-13
*/
void ppath_start_1d(
	      Index&      imax,
	      Index&      npl,
	      Index&      ip,
              Numeric&    r_start,
              Numeric&    lat_start,
              Numeric&    za_start,
	const Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        const Numeric&    r_geoid,
        const Numeric&    z_ground )
{
  // Number of points in the incoming ppath
  imax = ppath.np - 1;

  // Number of pressure levels
  npl = p_grid.nelem();

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.pos(imax,0);
  lat_start = ppath.pos(imax,1);
  za_start  = ppath.los(imax,0);

  // Asserts
  assert( npl >= 2 );
  assert( is_decreasing( p_grid ) );
  assert( is_size( z_field, npl ) );
  assert( is_increasing( z_field ) );
  assert( r_geoid > 0 );
  //
  assert( ppath.dim == 1 );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= ( npl - 2 ) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );
  //
  assert( r_start >= r_geoid + z_ground );
  assert( za_start >= 0  &&  za_start <= 180 );
  assert( !( is_gridpos_at_index_i( ppath.gp_p[imax], 0 ) && za_start > 90 ) );
  assert( !( is_gridpos_at_index_i( ppath.gp_p[imax], npl-1 )  && 
                                                            za_start <= 90 ) );

  // Determine index of the pressure surface being the lower limit for the
  // grid range of interest.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], za_start<=90 );
}



//! ppath_start_2d
/*! 
   Internal help function for 1D path calculations.

   The function does the asserts and determined some variables that are common
   for geometrical and refraction calculations.

   See the code fo details.

   \author Patrick Eriksson
   \date   2002-11-18
*/
void ppath_start_2d(
	      Numeric&    r_start,
	      Numeric&    lat_start,
	      Numeric&    za_start,
	      Index&      ip,
	      Index&      ilat,
	      Numeric&    r1,
	      Numeric&    r2,
	      Numeric&    r3,
	      Numeric&    r4,
 	      Numeric&    lat1,
	      Numeric&    lat3,
	      Numeric&    c2,
	      Numeric&    c4,
	const Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground )
{
  // Number of points in the incoming ppath
  const Index imax = ppath.np - 1;

  // Number of pressure levels and latitudes
  const Index npl = p_grid.nelem();
  const Index nlat = lat_grid.nelem();

  // Extract starting radius, zenith angle and latitude
  r_start   = ppath.pos(imax,0);
  lat_start = ppath.pos(imax,1);
  za_start  = ppath.los(imax,0);

  // First asserts (more below)
  assert( npl >= 2 );
  assert( is_decreasing( p_grid ) );
  assert( nlat >= 2 );
  assert( is_increasing( lat_grid ) );
  assert( is_size( z_field, npl, nlat ) );
  assert( is_size( r_geoid, nlat ) );
  assert( is_size( z_ground, nlat ) );
  //
  assert( ppath.dim == 2 );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= ( npl - 2 ) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );
  assert( ppath.gp_lat[imax].idx >= 0 );
  assert( ppath.gp_lat[imax].idx <= ( nlat - 2 ) );
  assert( ppath.gp_lat[imax].fd[0] >= 0 );
  assert( ppath.gp_lat[imax].fd[0] <= 1 );
  //
  assert( za_start >= -180  &&  za_start <= 180 );
  assert( !( is_gridpos_at_index_i( ppath.gp_p[imax], 0 )  &&  
                                                       fabs(za_start) > 90 ) );
  assert( !( is_gridpos_at_index_i( ppath.gp_p[imax], npl-1 )  && 
                                                      fabs(za_start) <= 90 ) );
  assert( !( is_gridpos_at_index_i( ppath.gp_lat[imax], 0 )  &&  
                                                              za_start < 0 ) );
  assert( !( is_gridpos_at_index_i( ppath.gp_lat[imax], nlat-1 )  &&  
                                                              za_start > 0 ) );
  // more asserts below ...


  // The corners and the faces of the grid cell are numbered in anti-clockwise
  // direction. The lower left corner is number 1. The left face is number 1.
  // For the coding of end point, the ground is given number 5 and a tangent
  // point 6.


  // Determine interesting latitude grid range and latitude end points of 
  // the range.
  //
  ilat = gridpos2gridrange( ppath.gp_lat[imax], za_start >= 0 );
  //
  lat1 = lat_grid[ilat];
  lat3 = lat_grid[ilat+1];

  // Latitude distance between start point and left grid cell boundary
  const Numeric dlat_left  = lat_start - lat1;

  // Determine interesting pressure grid range. Do this first assuming that
  // the pressure surfaces are not tilted (that is, fabs(za_start<=90) always
  // mean upward observation). 
  // Set radius for the corners of the grid cell and the radial slope of 
  // pressure surface limits of the grid cell to match the found ip.
  //
  ip = gridpos2gridrange( ppath.gp_p[imax], fabs(za_start) <= 90);
  //
  r1 = r_geoid[ilat] + z_field(ip,ilat);        // lower-left
  r2 = r_geoid[ilat+1] + z_field(ip,ilat+1);    // lower-right
  r3 = r_geoid[ilat+1] + z_field(ip+1,ilat+1);  // upper-right
  r4 = r_geoid[ilat] + z_field(ip+1,ilat);      // upper-left
  c2 = psurface_slope_2d( lat_grid, r_geoid, 
         z_field(ip,Range(joker)), ppath.gp_lat[imax], za_start >= 0 );
  c4 = psurface_slope_2d( lat_grid, r_geoid, 
       z_field(ip+1,Range(joker)), ppath.gp_lat[imax], za_start >= 0 );

  // Check if the LOS zenith angle happen to be between 90 and the zenith angle
  // of the pressure surface (that is, 90 + tilt of pressure surface), and in
  // that case if ip must be changed. This check is only needed when the
  // start point is on a pressure surface. We can then take the oppertunity
  // to assert that the start radius is then consistent with gp_p.
  //
  if( is_gridpos_at_index_i( ppath.gp_p[imax], ip )  )
    {
      assert( fabs( r_start - ( r1 + c2 * dlat_left ) ) <= R_EPS );
      Numeric tilt = psurface_tilt_2d( r_start, c2 );
      if( is_los_downwards_2d( za_start, tilt ) )
	{
	  ip--;
	  r4 = r1;   r3 = r2;   c4 = c2;
	  r1 = r_geoid[ilat] + z_field(ip,ilat);
          r2 = r_geoid[ilat+1] + z_field(ip,ilat+1);
          c2 = psurface_slope_2d( lat_grid, r_geoid, 
                 z_field(ip,Range(joker)), ppath.gp_lat[imax], za_start >= 0 );
	}
    }
  else if( is_gridpos_at_index_i( ppath.gp_p[imax], ip+1 )  )
    {
      assert( fabs( r_start - ( r4 + c4 * dlat_left ) ) <= R_EPS );
      Numeric tilt = psurface_tilt_2d( r_start, c4 );
      if( !is_los_downwards_2d( za_start, tilt ) )
	{
	  ip++;
	  r1 = r4;   r2 = r3;   c2 = c4;
	  r3 = r_geoid[ilat+1] + z_field(ip+1,ilat+1);
	  r4 = r_geoid[ilat] + z_field(ip+1,ilat);    
	  c4 = psurface_slope_2d( lat_grid, r_geoid, 
               z_field(ip+1,Range(joker)), ppath.gp_lat[imax], za_start >= 0 );
	}
    }

  out3 << "  pressure grid range  : " << ip << "\n";
  out3 << "  latitude grid range  : " << ilat << "\n";

  // As a double check: 
  // assert that start position is inside the found grid cell.
  assert( lat_start >= lat1  &&  lat_start <= lat3 );
  assert( r1 < r4 );
  assert( r2 < r3 );
  assert( r_start >= r1 + c2 * dlat_left - R_EPS );
  assert( r_start <= r4 + c4 * dlat_left + R_EPS );
}



/*===========================================================================
  === Core functions for geometrical *ppath_step* functions
  ===========================================================================*/

//! ppath_step_geom_1d
/*! 
   Calculates 1D geometrical propagation path steps.

   This is the core function to determine 1D propagation path steps by pure
   geometrical calculations. Path points are included for crossings with the 
   grids, tangent points and points of ground intersections. In addition,
   points are included in the propgation path to ensure that the distance
   along the path between the points does not exceed the selected maximum 
   length (lmax). If lmax is <= 0, this means that no length criterion shall
   be applied.

   Note that the input variables are here compressed to only hold data for
   a 1D atmosphere. For example, z_field is z_field(:,0,0).

   For more information read the chapter on propagation paths in AUG.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   z_field            Geometrical altitudes corresponding to p_grid.
   \param   r_geoid           Geoid radius.
   \param   z_ground          Ground altitude.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-05-20
*/
void ppath_step_geom_1d(
	      Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        const Numeric&    r_geoid,
        const Numeric&    z_ground,
	const Numeric&    lmax )
{
  // Number of pressure levels and points in the incoming ppath
  Index npl, imax;

  // Starting radius, zenith angle and latitude
  Numeric r_start, lat_start, za_start;

  // Index of the pressure surface being the lower limit for the
  // grid range of interest.
  Index ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( imax, npl, ip, r_start, lat_start, za_start, 
                                   ppath, p_grid, z_field, r_geoid, z_ground );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  Numeric ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // Get end radius of the path step (r_end). If looking downwards, it must 
  // be checked if:
  //    a tangent point is passed
  //    there is an intersection with the ground
  //
  Numeric r_end;
  bool    tanpoint = false, ground = false;
  //
  if( za_start <= 90 )
    { r_end = r_geoid + z_field[ ip + 1 ]; }
  else
    {
      // Lowest possible radius for the path step
      Numeric r_lowest  = r_geoid + z_field[ ip ];

      // Ground radius
      Numeric r_ground = r_geoid + z_ground;

      // The tangent radius equals here ppc.

      if( ( r_lowest > r_ground )  &&  ( r_lowest > ppc ) )
	{
	  r_end    = r_lowest;
	}
      else if( ppc >= r_ground )
	{
	  r_end    = ppc;
	  tanpoint = true;
	}
      else
	{
	  r_end    = r_ground;
	  ground   = true;
	}
    }


  // Calculate basic variables from r_start to r_end.
  //
  Vector    r_v, lat_v, za_v;
  Numeric   lstep;
  //
  geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                                                       za_start, r_end, lmax );
  //
  const Index ilast = r_v.nelem() - 1;


  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure(  ppath, 1, r_v.nelem() );
  //
  if( lmax < 0 )
    { ppath.method     = "1D geometrical"; }
  else
    { ppath.method     = "1D geometrical with length criterion"; }
  ppath.refraction = 0;
  ppath.constant   = ppc;
  //
  ppath_fill_1d( ppath, r_v, lat_v, za_v, Vector(ilast,lstep), r_geoid, 
                                                                 z_field, ip );


  // Different options depending on position of end point of step:

  //--- End point is the ground
  if( ground )
    { ppath_set_background( ppath, 2 ); }

  //--- End point is a tangent point
  else if( tanpoint )
    {
      ppath.tan_pos.resize(2);
      ppath.tan_pos[0] = r_v[ilast];
      ppath.tan_pos[1] = lat_v[ilast];

      // Make part from tangent point and up to the starting pressure level.
      //
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      ppath_step_geom_1d( ppath2, p_grid, z_field, r_geoid, z_ground, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }

  //--- End point is on top of a pressure surface
  else
    {
      gridpos_force_end_fd( ppath.gp_p[ilast] );
    }
}



//! ppath_step_geom_2d
/*! 
   Calculates 2D geometrical propagation path steps.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   z_field           Geometrical altitudes
   \param   r_geoid           Geoid radii.
   \param   z_ground          Ground altitudes.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_step_geom_2d(
	      Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
	const Numeric&    lmax )
{
  // Radius, zenith angle and latitude of start point.
  Numeric   r_start, lat_start, za_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat;

  // Radius for corner points, and latitude and slope of faces of the grid cell
  //
  // The corners and the faces of the grid cell are numbered in anti-clockwise
  // direction. The lower left corner is number 1. The left face is number 1.
  // For the coding of end point, the ground is given number 5 and a tangent
  // point 6.
  Numeric   r1, r2, r3, r4, lat1, lat3, c2, c4;

  // Determine the variables defined above and make all possible asserts
  ppath_start_2d( r_start, lat_start, za_start, ip, ilat, 
                  r1, r2, r3, r4, lat1, lat3, c2, c4, 
                         ppath, p_grid, lat_grid, z_field, r_geoid, z_ground );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  Numeric ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // Define some useful variables:
  //
  // Number of points in the incoming ppath
  const Index imax = ppath.np - 1;
  //
  // Latitude distance between start point and left grid cell boundary
  const Numeric dlat_left  = lat_start - lat1;
  //
  // Latitude distance to latitude end face in the viewing direction
  Numeric dlat_endface;
  if( za_start >= 0 )
    { dlat_endface = lat3 - lat_start; }
  else
    { dlat_endface = -dlat_left; }


  // We shall now determine at what face of the grid cell the path step ends.
  // If the step passes a tangent point or the ground, we call this function
  // iteratively for the second part.

  // The end point is most easily determined by the latitude difference to 
  // the start point. For some cases there exist several crossings and we want 
  // the one closest in latitude to the start point. The latitude distance 
  // for the crossing shall not exceed dlat2end.
  // The variable endface is a number coding of at what cell face (including 
  // the ground and tangent points) the end point is found. This is needed to
  // set values for the end point in most accurate way.
  //
  Numeric dlat2end = 999;
  Index   endface  = 999;   

  // --- Lower face (pressure surface ip).
  //
  // This face is tricky as there can two crossings with the pressure surface
  // before the next latitude grid point is reached. This is the case as the
  // face is bended inwards. 
  // The zenith angles to the corner points of the cell cannot be used to 
  // determine if there is a crossing or not. Instead we have to call 
  // psurface_crossing_2d for all cases to test if there is a crossing.
  // We don't need to consider the face if we are standing on the pressure 
  // surface.
  //
  if( !is_gridpos_at_index_i( ppath.gp_p[imax], ip )  )
    {
      Numeric r_surface = r1+c2*dlat_left;

      // We must check that numerical inaccuracies not have caused that the
      // starting radius is below the pressure surface.
      // Similar checks are done below for the ground and the upper pressure
      // surface.
      //
      Numeric r_tmp = r_start;
      //
      if( r_tmp < r_surface )
	{ r_tmp = r_surface; }

      dlat2end = psurface_crossing_2d( r_tmp, za_start, r_surface, c2, 0 );
      endface = 2;  // This variable will be re-set if there was no crossing
    }

  // --- The ground.
  //
  // Ground radius at latitude end points, and ground slope
  const Numeric rground1 = r_geoid[ilat] + z_ground[ilat];
  const Numeric rground2 = r_geoid[ilat+1] + z_ground[ilat+1];
  const Numeric cground = psurface_slope_2d( lat_grid, r_geoid, 
		                 z_ground, ppath.gp_lat[imax], za_start >= 0 );
  //
  {
    Numeric r_ground = rground1 + cground * dlat_left;

    assert( r_start >= r_ground - R_EPS );

    // Check shall be done only if the ground is, at least partly, inside 
    //the grid cell.
    //
    if( rground1 >= r1  ||  rground2 >= r2 )
      {
	Numeric r_tmp = r_start;
	if( r_tmp < r_ground )
	  { r_tmp = r_ground; }

	Numeric dlat2ground = psurface_crossing_2d( r_tmp, za_start, r_ground,
						                  cground, 0 );
	if( fabs(dlat2ground) <= fabs(dlat2end) )
	  {
	    dlat2end = dlat2ground;
	    endface  = 5;
	  }
      }
  }

  // If dlat2end <= dlat_endface we are ready. Otherwise we have to check
  // remaining cell faces. The same applies after testing upper face.

  // --- Upper face  (pressure surface ip+1).
  //
  if( fabs(dlat2end) > fabs(dlat_endface) )
    {
      // We can here determine by zenith angles if there is a crossing with
      // the pressure surface. This should save some time compared to call
      // psurface_crossing_2d blindly.
      // Note that the path step can both start and end at the face.

      Numeric r_surface = r4 + c4 * dlat_left;

      Numeric r_tmp = r_start;
      if( is_gridpos_at_index_i( ppath.gp_p[imax], ip+1 ) )
	{ r_tmp = r_surface; }
      else if( r_tmp > r_surface )
	{ r_tmp = r_surface; }

      if( za_start >= za_geom2other_point( r_tmp, lat_start, r4, lat1 )  &&
                za_start <= za_geom2other_point( r_tmp, lat_start, r3, lat3 ) )
	{
	  dlat2end = psurface_crossing_2d( r_tmp, za_start, r_surface, c4, 0 );
	  assert( dlat2end < 999 );
	  endface  = 4;
	}
    }

  // Left or right end face
  if( fabs(dlat2end) > fabs(dlat_endface) )
    { 
      dlat2end = dlat_endface; 
      if( za_start >= 0 )
	{ endface  = 3; }
      else
	{ endface  = 1; }
    }

  // Check if a tangent point is passed before dlat2end is reached.
  if( fabs(za_start) > 90  &&  ( fabs(za_start) - fabs(dlat2end) ) < 90 ) 
    { endface  = 6; }


  // Calculate radius for end point.
  // To obtain best possible accuracy it is calculated folowing found end face,
  // and not based on dlat2end.
  //
  Numeric r_end = -1;
  //
  out3 << "  end face number code : " << endface << "\n";
  //
  if( endface == 1 )
    { r_end = geompath_r_at_lat( ppc, lat_start, za_start, lat1 ); }
  else if( endface == 2 )
    { r_end = r1 + c2 * ( dlat_left + dlat2end ); }
  else if( endface == 3 )
    { r_end = geompath_r_at_lat( ppc, lat_start, za_start, lat3 ); }
  else if( endface == 4 )
    { r_end = r4 + c4 * ( dlat_left + dlat2end ); }
  else if( endface == 5 )
    { r_end = rground1 + cground * ( dlat_left + dlat2end ); }
  else if( endface == 6 )
    { r_end = geompath_r_at_za( ppc, sign(za_start) * 90 ); }


  // Calculate basic variables from r_start to r_end.
  //
  Vector    r_v, lat_v, za_v;
  Numeric   lstep;
  //
  geompath_from_r1_to_r2( r_v, lat_v, za_v, lstep, ppc, r_start, lat_start, 
                                                       za_start, r_end, lmax );
  //
  const Index ilast = r_v.nelem() - 1;

  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure(  ppath, 2, r_v.nelem() );
  //
  if( lmax < 0 )
    { ppath.method     = "2D basic geometrical"; }
  else
    { ppath.method     = "2D geometrical with length criterion"; }
  ppath.refraction = 0;
  ppath.constant   = ppc;
  //
  ppath_fill_2d( ppath, r_v, lat_v, za_v, lstep, r_geoid, z_field, lat_grid, 
		                                                    ip, ilat );
  //
  if( endface == 5 )
    { ppath_set_background( ppath, 2 ); }
  else if( endface == 6 )
    {
      ppath.tan_pos.resize(2);
      ppath.tan_pos[0] = r_v[ilast];
      ppath.tan_pos[1] = lat_v[ilast];
    }

  // To avoid numerical inaccuracy as far as possible, we set values to match
  // the end of grid ranges when possible (radius of end point is handled
  // by geompath_from_r1_to_r2):
  //
  if( endface == 1 )
    {
      ppath.pos(ilast,1) = lat1;
      gridpos_force_end_fd( ppath.gp_lat[ilast] );
    }
  else if( endface == 2  ||  endface == 4 )
    { gridpos_force_end_fd( ppath.gp_p[ilast] ); }
  else if( endface == 3 )
    {
      ppath.pos(ilast,1) = lat3;
      gridpos_force_end_fd( ppath.gp_lat[ilast] );
    }
  else if( endface == 6 )
    { ppath.los(ilast,0) = sign(za_start)*90; }


  // Make part after a tangent point.
  //
  if( endface == 6 )
    {
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      // Call this function recursively
      ppath_step_geom_2d( ppath2, p_grid, lat_grid, z_field,
			                             r_geoid, z_ground, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }
}



/*===========================================================================
  === Core functions for refraction *ppath_step* functions
  ===========================================================================*/

//! ppath_step_refr_euler_1d
/*! 
   Calculates 1D propagation path steps, with refraction, using an Euler
   approach.

   This function works as the function *ppath_step_geom_1d* but considers
   also refraction. The upper length of the ray tracing steps is set by
   the argument *lraytrace*. This argument controls only the internal
   calculations. The maximum distance between the path points is still
   determined by *lmax*.

   A geometrical step with length of *lraytrace* is taken from each
   point.  The zenith angle for the end point of that step is
   calculated exactly by the expression c = r*n*sin(theta), and a new
   step is taken. The length of the last ray tracing step to reach the
   end radius is adopted to the distance to the end radius. The
   calculations are always performed from the lowest to the highest
   radius.

   For more information read the chapter on propagation paths in AUG.
   The algorithm used is described in that part of AUG.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   z_field           Geometrical altitudes corresponding to p_grid.
   \param   t_field           Temperatures corresponding to p_grid.
   \param   r_geoid           Geoid radius.
   \param   z_ground          Ground altitude.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-11-26
*/
void ppath_step_refr_euler_1d(
	      Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        const Numeric&    r_geoid,
        const Numeric&    z_ground,
	const Numeric&    lraytrace,
	const Numeric&    lmax )
{
  // Number of pressure levels and points in the incoming ppath
  Index npl, imax;

  // Starting radius, zenith angle and latitude
  Numeric r_start, lat_start, za_start;

  // Index of the pressure surface being the lower limit for the
  // grid range of interest.
  Index ip;

  // Determine the variables defined above, and make asserts of input
  ppath_start_1d( imax, npl, ip, r_start, lat_start, za_start, 
                                   ppath, p_grid, z_field, r_geoid, z_ground );

  // Assert not done for geometrical calculations
  assert( t_field.nelem() == npl );
  assert( lraytrace > 0 );
  assert( lmax < 0  ||  lmax >= lraytrace );

  // Refractive index at start point
  Numeric nvalue = get_refr_index_1d( p_grid, z_field, t_field, 
                                                           r_start - r_geoid );


  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  Numeric ppc;
  if( ppath.constant < 0 )
    { 
      // If the sensor is placed outside the atmosphere, the constant is
      // already set.
      ppc = refraction_ppc( r_start, za_start, nvalue ); 
    }
  else
    { ppc = ppath.constant; }


  // Get end radius of the path step (r_end). If looking downwards, it must 
  // be checked if:
  //    a tangent point is passed
  //    there is an intersection with the ground
  //
  Numeric r_end;
  bool    tanpoint = false, ground = false;
  //
  if( za_start <= 90 )
    { r_end = r_geoid + z_field[ ip + 1 ]; }

  else
    {
      // Lowest possible radius for the path step
      const Numeric r_lowest  = r_geoid + z_field[ ip ];

      // Refractive index at r_lowest
      const Numeric n_lowest = get_refr_index_1d( p_grid, z_field, t_field, 
                                                          r_lowest - r_geoid );
      // Ground radius
      const Numeric r_ground = r_geoid + z_ground;

      // Is the pressure surface below the end point?
      if( ( r_lowest > r_ground )  &&  ( r_lowest*n_lowest > ppc ) )
	{
	  r_end    = r_lowest;
	}

      // Is the tangent point is the end point?
      else if( ppc >= r_ground * get_refr_index_1d( p_grid, z_field, t_field, 
                                                         r_ground - r_geoid ) )
	{
	  	Numeric   r_tan = r_lowest;
	  	          r_end = r_tan + 999e3;   
	  	Numeric   n_end = n_lowest;
	        bool      first = true;
	  const Numeric   accuracy = 0.1;       // 0.1m 

	  // We loop until it is sure that accuracy is better then specified
	  while( fabs( r_end - r_tan ) > accuracy )
	    {
	      r_end = r_tan;

	      if( !first )
		{
		  n_end = get_refr_index_1d( p_grid, z_field, t_field, 
                                                             r_end - r_geoid );
		}
	      first = false;

	      const Numeric dn = ( nvalue - n_end ) / ( r_start - r_end );
	  
	      if( fabs(dn) < 1e-15 )        // An arbitrary threshold! 
		{                           // dn at ground level is about 3e-8
		  r_tan = ppc / n_end;
		}
	      else
		{
		  // See AUG for algorith used here
		  const Numeric x = ( n_end - dn * r_end ) / ( 2 * dn );
		  r_tan = -x + sign(dn) * sqrt( x*x + ppc/dn );
		}
	    }

	  tanpoint = true;
	}

      // Ground must be the end point!
      else
	{
	  r_end    = r_ground;
	  ground   = true;
	}
    }


  //--- Ray tracing
  //
  // Create some variables used for the ray tracing.
  // These variables are needed as the ray tracing is always done from the
  // lowest to the highest radius.
  //
  Numeric   r, rstop, za, lat=0;
  //
  if( r_start < r_end )
    {
      r     = r_start;
      rstop = r_end;
      za    = za_start;
    }
  else
    {
      r     = r_end;
      rstop = r_start;
      if( tanpoint )    // Here we don't know ZA for the start point.
	{ za = 90; }
      else
	{ 
	  const Numeric n_end = get_refr_index_1d( p_grid, z_field, t_field, 
                                                             r_end - r_geoid );
	  za = geompath_za_at_r( ppc/n_end, za_start, r_end ); 
	}
    }
  //
  // Perform the ray tracing. 
  // As we don't know how many points we will get, the result is first stored
  // in numeric arrays, and is later copied to vectors.
  //
  Array<Numeric>   r_rt1, lat_rt1, za_rt1, l_rt1;
  //
  r_rt1.push_back( r );
  lat_rt1.push_back( lat );
  za_rt1.push_back( za );
  //
  bool ready = false;
  //
  while( !ready )
    {
      const Numeric ppc_step = r * sin( DEG2RAD * za );
            Numeric l_step = lraytrace;

      Numeric r_new = geompath_r_at_l( ppc_step, 
                                     geompath_l_at_r( ppc_step, r ) + l_step );

      if( r_new >= rstop )
	{
	  r_new   = rstop;
	  l_step = geompath_l_at_r( ppc_step, r_new ) - 
                                                geompath_l_at_r( ppc_step, r );
	  ready  = true;
	}
      
      const Numeric n_new = get_refr_index_1d( p_grid, z_field, t_field, 
                                                             r_new - r_geoid );

      const Numeric   ppc_local = ppc / n_new; 
      
      if( ppc_local < r_new )
	{ za = geompath_za_at_r( ppc/n_new, za_start, r_new ); }
      else
	{ za = 90; }
      
      lat += RAD2DEG * acos( ( r_new*r_new + r*r - l_step*l_step ) / 
                                                             ( 2 * r_new*r ) );
      r = r_new;

      r_rt1.push_back( r );
      l_rt1.push_back( l_step );
      lat_rt1.push_back( lat );
      za_rt1.push_back( za );
    }
  //
  // Number of ray tracing points
  const Index   nrt=r_rt1.nelem();
  //
  // Copy to vectors
  // The copying differ depending on the true direction (upward/downward) of
  // the path step.
  Vector    r_rt(nrt), lat_rt(nrt), za_rt(nrt), l_rt(nrt-1);    
  //
  if( r_start < r_end )
    {
      for( Index i=0; i<nrt; i++ )
	{
	  r_rt[i]   = r_rt1[i];
	  lat_rt[i] = lat_start + lat_rt1[i]; 
	  za_rt[i]  = za_rt1[i];
	  if( i < (nrt-1) )
	    { l_rt[i]   = l_rt1[i]; }
	}
    }
  else
    {
      for( Index i=0; i<nrt; i++ )
	{
	  const Index i1 = nrt - 1 - i;
	  r_rt[i]   = r_rt1[i1];
	  lat_rt[i] = lat_start + lat_rt1[nrt-1] - lat_rt1[i1]; 
	  za_rt[i]  = za_rt1[i1];
	  if( i < (nrt-1) )
	    { l_rt[i]   = l_rt1[i1-1]; }
	}
    }

  //VectorPrint(r_rt,"r_rt");
  //VectorPrint(lat_rt,"lat_rt");
  //VectorPrint(za_rt,"za_rt");
  //VectorPrint(l_rt,"l_rt");


  // Interpolate the radii, zenith angles and latitudes to a set of points
  // linearly spaced along the path. If *lmax* <= 0, then only the end points
  // shall be kept.
  //
  Index     n = 2;
  Numeric   ltotsum = l_rt.sum();
  //
  if( lmax > 0 )
    { n = Index( ceil( ltotsum / lmax ) ) + 1; }
  //
  Vector    r_v(n), lat_v(n), za_v(n), l_v(n-1);
  //
  if( n == 2 )
    {
      r_v.resize(2);   r_v[0] = r_rt[0];     r_v[1] = r_rt[nrt-1];
      za_v.resize(2);  za_v[0] = za_rt[0];   za_v[1] = za_rt[nrt-1];
      lat_v.resize(2); lat_v[0] = lat_rt[0]; lat_v[1] = lat_rt[nrt-1];
      l_v.resize(1);   l_v[0] = ltotsum;
    }
  else
    {
      Vector           lsum(nrt), llin(n);
      ArrayOfGridPos   gp(n);
      Matrix           itw(n,2);
      
      lsum[0] = 0;
      for( Index i=1; i<nrt; i++ )
	{ lsum[i] = lsum[i-1] + l_rt[i-1]; }

      nlinspace( llin, 0, ltotsum, n );

      gridpos( gp, lsum, llin );
      gridpos_force_end_fd( gp[0] );
      gridpos_force_end_fd( gp[n-1] );

      interpweights( itw, gp );

      interp( r_v, itw, r_rt, gp );
      interp( za_v, itw, za_rt, gp );
      interp( lat_v, itw, lat_rt, gp );
      l_v = llin[1] - llin[0];
    }


  // Re-allocate ppath for return results and fill the structure
  //
  ppath_init_structure(  ppath, 1, n );
  //
  if( lmax < 0 )
    { ppath.method     = "1D with refraction (Euler)"; }
  else
    { ppath.method     = "1D with refraction and length criterion (Euler)"; }
  ppath.refraction = 1;
  ppath.constant   = ppc;
  //
  ppath_fill_1d( ppath, r_v[Range(0,n)], lat_v[Range(0,n)], 
                   za_v[Range(0,n)], l_v[Range(0,n-1)], r_geoid, z_field, ip );


  // Different options depending on position of end point of step:

  //--- End point is the ground
  if( ground )
    { ppath_set_background( ppath, 2 ); }

  //--- End point is a tangent point
  else if( tanpoint )
    {
      ppath.tan_pos.resize(2);
      ppath.tan_pos[0] = r_v[n-1];
      ppath.tan_pos[1] = lat_v[n-1];

      // Make sure that last zenith angle is exactly 90 degrees.
      // Otherwise "another" tangent point can be found in the recursive call.
      ppath.los(n-1,0) = 90;

      // Make part from tangent point and up to the starting pressure level.
      //
      Ppath ppath2;
      ppath_init_structure( ppath2, ppath.dim, ppath.np );
      ppath_copy( ppath2, ppath );

      out3 << "  --- Recursive step to include tangent point --------\n"; 

      ppath_step_refr_euler_1d( ppath2, p_grid, z_field, t_field, r_geoid, 
                                                   z_ground, lraytrace, lmax );

      out3 << "  ----------------------------------------------------\n"; 

      // Combine ppath and ppath2
      ppath_append( ppath, ppath2 );
    }

  //--- End point is on top of a pressure surface
  else
    {
      gridpos_force_end_fd( ppath.gp_p[n-1] );
    }
}



//! ppath_step_refr_euler_2d
/*! 
   Calculates 1D propagation path steps, with refraction, using a simple
   and fast ray tracing scheme.

   Works as the same function for 1D despite that some input arguments are
   of different type.

   \param   ppath             Output: A Ppath structure.
   \param   p_grid            Pressure grid.
   \param   lat_grid          Latitude grid.
   \param   z_field           Geometrical altitudes.
   \param   t_field           Atmospheric temperatures.
   \param   r_geoid           Geoid radii.
   \param   z_ground          Ground altitudes.
   \param   lraytrace         Maximum allowed length for ray tracing steps.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-07-03
*/
void ppath_step_refr_euler_2d(
	      Ppath&      ppath,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstVectorView   r_geoid,
        ConstVectorView   z_ground,
	const Numeric&    lraytrace,
	const Numeric&    lmax )
{
  // Radius, zenith angle and latitude of start point.
  Numeric   r_start, lat_start, za_start;

  // Lower grid index for the grid cell of interest.
  Index   ip, ilat;

  // Radius for corner points, and latitude and slope of faces of the grid cell
  //
  // The corners and the faces of the grid cell are numbered in anti-clockwise
  // direction. The lower left corner is number 1. The left face is number 1.
  // For the coding of end point, the ground is given number 5 and a tangent
  // point 6.
  Numeric   r1, r2, r3, r4, lat1, lat3, c2, c4;

  // Determine the variables defined above and make all possible asserts
  ppath_start_2d( r_start, lat_start, za_start, ip, ilat, 
                  r1, r2, r3, r4, lat1, lat3, c2, c4, 
                         ppath, p_grid, lat_grid, z_field, r_geoid, z_ground );

  // No constant for the path is valid here.


  // Refractive index for the corner points of the grid cell
  const Numeric   n1 = get_refr_index_2d( p_grid, lat_grid, r_geoid, 
                                                  z_field, t_field, r1, lat1 );
  const Numeric   n2 = get_refr_index_2d( p_grid, lat_grid, r_geoid,
                                                  z_field, t_field, r2, lat3 );
  const Numeric   n3 = get_refr_index_2d( p_grid, lat_grid, r_geoid,
                                                  z_field, t_field, r3, lat3 );
  const Numeric   n4 = get_refr_index_2d( p_grid, lat_grid, r_geoid,
                                                  z_field, t_field, r4, lat1 );


  //--- Ray tracing
  //
  Array<Numeric>   r_rt1, lat_rt1, za_rt1, l_rt1;
  Numeric          r = r_start, lat = lat_start, za = za_start; 
  //
  r_rt1.push_back( r );
  lat_rt1.push_back( lat );
  za_rt1.push_back( za );
  //
  bool   ready = false;
  //
  while( !ready )
    {
      // Radial and latitudinal gradient of refractive index
      Numeric   dndr, dndlat;
      refraction_gradient_2d( dndr, dndlat, r1,  r4, c2, c4, lat1, lat3, 
			                              n1, n2, n3, n4, r, lat );

      // Perform ray a tracing step of length *lraytrace*.
      //
      const Numeric ppc_step = r * sin( DEG2RAD * za );
            Numeric l_step = lraytrace;

      Numeric r_new = geompath_r_at_l( ppc_step, 
                                     geompath_l_at_r( ppc_step, r ) + l_step );

      Numeric lat_new = lat + RAD2DEG * acos( ( r_new*r_new + r*r - 
                                           l_step*l_step ) / ( 2 * r_new*r ) );

      // Check if the new point is outside of the grid cell. 
      // If yes, set *r_new* and *lat_new* to the crossong point.

      // Store found point
      r_rt1.push_back( r );
      lat_rt1.push_back( lat );
      za_rt1.push_back( za );
    }
}

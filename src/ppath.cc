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



/*****************************************************************************
 ***  File description 
 *****************************************************************************/

/*!
  \file   ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-02
  
  \brief  Functions releated to calculation of propagation paths.
  
  Functions to determine propagation paths for different atmospheric
  dimensionalities.

  The term propagation path is here shortened to ppath.
*/



/*****************************************************************************
 *** External declarations
 *****************************************************************************/

#include <math.h>
#include <stdexcept>
#include "array.h"
#include "check_input.h"
#include "math_funcs.h"
#include "mystring.h"
#include "interpolation.h"
#include "ppath.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



/*****************************************************************************
 *** Functions related to geometrical propagation paths
 *****************************************************************************/

//! geometrical_tangent_radius
/*! 
   Calculates the radius of the geometrical tangent point. 

   If the zenith angle is smaller than 90 degrees, the tangent point is 
   an imaginary point behind the sensor.

   Positive and negative zenith angles are handled.

   \return         The tangent point radius.
   \param   r      The radius of the sensor position.
   \param   za     The zenith angle of the sensor line-of-sight.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Numeric geometrical_tangent_radius( const Numeric& r, const Numeric& za )
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

   Positive and negative zenith angles are handled.

   \return         The zenith angle at the second point.
   \param   r0     The radius of the starting point.
   \param   za0    The zenith angle of the starting point.
   \param   r      The radius of the second point.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Numeric geompath_za_at_r(
       const Numeric&   r_tan,
       const Numeric&   a_za,
       const Numeric&   r )
{
  assert( r_tan > 0 );
  assert( fabs(a_za) <= 180 );
  assert( r >= r_tan );
  Numeric za;
  za = RAD2DEG * asin( r_tan / r );
  if( fabs(a_za) > 90 )
    za = 180 - za;
  if( a_za < 0 )
    za = -za;
  return za;
}



//! geompath_lat_at_za
/*!
   Calculates the latitude for a given zenith angle along a geometrical 
   propagation path.

   Positive and negative zenith angles are handled.

   \return         The latitude difference for the second point.
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
  assert( fabs(za0) >= fabs(za) );
  assert( (za0>=0&&za>=0) || (za0<0&&za<0) );
  const Numeric dlat = za0 - za;
  if( za0 > 0 )
    return lat0 + dlat;
  else
    return lat0 - dlat;
}



//! geompath_l_at_r
/*!
   Calculates the length from the tangent point for the given radius.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         The length along the path from the tangent point.
   \param   r_tan  Tangent radius
   \param   r      The radius of the point of concern.

   \author Patrick Eriksson
   \date   2002-05-20
*/
Numeric geompath_l_at_r(
       const Numeric&   r_tan,
       const Numeric&   r )
{
  assert( r_tan > 0 );
  assert( r >= r_tan );
  
  // Double is hard-coded here to improve accuracy
  double a=r_tan*r_tan, b=r*r;

  return sqrt( b - a );
}



//! geompath_r_at_l
/*!
   Calculates the radius for a distance from the tangent point.

   The tangent point is either rwal or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         The radius. 
   \param   r_tan  Tangent radius
   \param   l      The length from the tangent point.

   \author Patrick Eriksson
   \date   2002-05-20
*/
Numeric geompath_r_at_l(
       const Numeric&   r_tan,
       const Numeric&   l )
{
  assert( r_tan > 0 );
  assert( l >= 0 );
  
  // Double is hard-coded here to improve accuracy
  double a=r_tan*r_tan, b=l*l;

  return sqrt( b + a );
}





/*****************************************************************************
 *** Angle(s) after a ground reflection
 *****************************************************************************/

//! angle_after_ground_1D 
/*!
   Calculates the LOS after a ground reflection for 1D.

   The zenith angle before the reflection msut be > 90 degrees.

   \param   za     Output: Zenith angle.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void angle_after_ground_1D( Numeric& za )
{
  assert( za > 90 );
  za = 180 - za;
}





/*****************************************************************************
 *** Functions operating on the Ppath structure
 *****************************************************************************/

//! ppath_init_structure
/*!
   Initiates a Ppath structure to hold the given number of points.

   All fields releated with the ground, symmetry and tangent point are set
   to 0 or empty. The background field is set to background case 0.

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
      ppath.gp_lat.resize( np );
  ppath.z.resize( np );
  if( np > 1 )
    ppath.l_step.resize( np-1 );
  else
    ppath.l_step.resize( 0 );
  ppath_set_background( ppath, 0 );
  ppath.ground     = 0;
  ppath.i_ground   = 0;
  ppath.tan_pos.resize(0);
  ppath.symmetry   = 0;
  ppath.i_symmetry = 0;
}



//! ppath_set_background 
/*!
   Sets the background field of a Ppath structure.

   The different background cases have a number coding to simplify a possible
   change of the strings and checking of the what case that is valid.

   The case numbers are:
   \begin{verbatim}
      0. Not yet set.
      1. Space.
      2. Blackbody ground.
      3. The surface of the cloud box.
      4. The interior of the cloud box.
   \end{verbatim}

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
      ppath.background = "blackbody ground";
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
   \param   ppath            Output: A Ppath structure.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Index ppath_what_background( 
	      Ppath&      ppath )
{
  if( ppath.background == "" )
    return 0;
  else if( ppath.background == "space" )
    return 1;
  else if( ppath.background == "blackbody ground" )
    return 2;
  else if( ppath.background == "cloud box surface" )
    return 3;
  else if( ppath.background == "cloud box interior" )
    return 4;
  else
    {
      ostringstream os;
      os << "The string " << ppath.background 
	 << " is not a valid background case.";
      throw runtime_error(os.str());
    }
}





/*****************************************************************************
 *** Help functions to ppathCalc
 *****************************************************************************/

//! ppath_start_stepping
/*!
   Initiates a Ppath structure for calculation of a path with *ppath_step*.

   The function performs two main tasks. As mentiuoned above, it initiates
   a Ppath structure (a), but it also checks that the end point of the path is
   at an allowed location (b).

   (a): The Ppath structure is set to hold the position and LOS of the last
   point of the path inside the atmsophere. This is point is either the
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
   fact *ppath_partial*.

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   r_geoid           The geoid radius.
   \param   z_ground          Ground altitude.
   \param   blackbody_ground  Flag for treating the ground as a blackbody.
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
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Index&          blackbody_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         a_pos,
        const Vector&         a_los )
{
  // This function contains no checks or asserts as it is only a sub-function
  // to ppathCalc where the input is checked carefully.

  // Assume that there is a point to put into ppath.
  // If the path is totally outside the atmosphere, the function shall be
  // called again with 0 as last argument.
  ppath_init_structure(  ppath, atmosphere_dim, 1 );

  // Number of pressure levels
  const Index np = p_grid.nelem();


  // The different atmospheric dimensionalities are handled seperately

  //-- 1D ---------------------------------------------------------------------
  if( atmosphere_dim == 1 )
    {
      // The only forbidden case here is that the sensor is below the ground
      if( a_pos[0] < r_geoid(0,0) + z_ground(0,0) )
	{
          ostringstream os;
          os << "The sensor is placed " 
             << (r_geoid(0,0) + z_ground(0,0) - a_pos[0])/1e3 << " km below "
             << "ground level.\nThe sensor must be above the ground.";
	  throw runtime_error(os.str());
	}

      // The sensor is inside the model atmosphere ----------------------------
      if( a_pos[0] < r_geoid(0,0) + z_field(np-1,0,0) )
	{
	  // Put some values into ppath. Use these values below (instead of
	  // a_pos and a_los) as they can be modified on the way.
	  ppath.pos(0,1) = 0; 
	  ppath.z[0]     = a_pos[0] - r_geoid(0,0);
          ppath.los(0,0) = a_los[0];
     
	  // Is the sensor on the ground looking down?
	  if( ( ppath.pos(0,1) < r_geoid(0,0)+z_ground(0,0) ) 
                                                   && ( ppath.los(0,0) > 90 ) )
	    {
	      angle_after_ground_1D( ppath.los(0,0) );
	      ppath.ground = 1;
	      ppath.i_ground = 0;
	      if( blackbody_ground )
		ppath_set_background( ppath, 2 );
	    }

	  // Check sensor position with respect to cloud box, if activated
          // and not background set to blackbody ground
	  if( cloudbox_on && !ppath_what_background(ppath) )
	    {
	      // Radius of lower and upper limit of the cloud box
	      Numeric r_low = r_geoid(0,0) + z_field(cloudbox_limits[0],0,0);
	      Numeric r_upp = r_geoid(0,0) + z_field(cloudbox_limits[1],0,0);

	      // Is the sensor inside the cloud box?
	      if( ppath.z[0] > z_field(cloudbox_limits[0],0,0) && 
		                 ppath.z[0] < z_field(cloudbox_limits[1],0,0) )
		{
		  ppath_set_background( ppath, 4 );
		}

	      // Is the sensor on the surface of cloud box and looks into box?
	      if( ( ppath.z[0] == z_field(cloudbox_limits[0],0,0) && 
                                                        ppath.los(0,0) <= 90 )
                     || 
		  ( ppath.z[0] == z_field(cloudbox_limits[1],0,0) && 
                                                        ppath.los(0,0) > 90 ) )
		{
		  ppath_set_background( ppath, 3 );
		}
	    }
	}

      // The sensor is outside the model atmosphere ---------------------------
      else
	{
	  // Path is above the atmosphere
	  if( a_los[0] <= 90 || 
                        geometrical_tangent_radius(a_pos[0],a_los[0]) >=
  	                      (r_geoid(0,0) + z_field(np-1,0,0)) )
	    {
	      ppath_init_structure(  ppath, atmosphere_dim, 0 );
	      ppath_set_background( ppath, 1 );
	    }

	  // Path enters the atmosphere
	  else
	    {
	      ppath.z[0]     = z_field(np-1,0,0);
	      ppath.los(0,0) = geompath_za_at_r( 
                        geometrical_tangent_radius( a_pos[0], a_los[0] ),
                                         a_los[0], r_geoid(0,0) + ppath.z[0] );
	      ppath.pos(0,1) = geompath_lat_at_za( a_los[0], 0, 
                                                              ppath.los(0,0) );
	    }
	}
    }


  //-- 2D ---------------------------------------------------------------------
  else if( atmosphere_dim == 2 )
    {
      throw runtime_error("The function handles not yet 2D.");
    }


  //-- 3D ---------------------------------------------------------------------
  else if( atmosphere_dim == 3 )
    {
      throw runtime_error("The function handles not yet 3D.");
    }
  //-- End of 1-3D ------------------------------------------------------------


  // Get grid positions and pressure for the end point, if there is one.
  if( ppath.np == 1 )
    {
      // Pressure
      gridpos( ppath.gp_p, z_field(Range(joker),0,0), ppath.z );
      Matrix itw( 1, 2 );
      interpweights( itw, ppath.gp_p );
      interp( ppath.pos(Range(joker),0), itw, p_grid, ppath.gp_p );

      // Latitude
      if( atmosphere_dim >= 2 )
	gridpos( ppath.gp_lat, lat_grid, ppath.pos(Range(joker),1) );

      // Longitude
      if( atmosphere_dim == 3 )
	gridpos( ppath.gp_lon, lon_grid, ppath.pos(Range(joker),2) );
    }
}





/*****************************************************************************
 *** Core functions for *ppath_step* functions
 *****************************************************************************/

//! ppath_step_1d_geom
/*! 
 This is just a test text.

 Some more text.

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   z_grid            The geometrical altitude corresponding to p_grid.
   \param   r_geoid           The geoid radius.
   \param   z_ground          Ground altitude.
   \param   blackbody_ground  Flag for treating the ground as a blackbody.
   \param   lmax              Maximum allowed length between the path points.

   \author Patrick Eriksson
   \date   2002-05-20
*/
void ppath_step_1d_geom(
	      Ppath&      ppath,
        const Index&      atmosphere_dim,
        ConstVectorView   p_grid,
        ConstVectorView   z_grid,
        const Numeric&    r_geoid,
        const Numeric&    z_ground,
        const Index&      blackbody_ground,
	const Numeric&    lmax )
{
  // Number of points in the incoming ppath
  const Index imax = ppath.np - 1;

  // First asserts (more below)
  assert( ppath.dim == atmosphere_dim );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= (p_grid.nelem()-2) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );

  // Start and end point mean here the point where the calculations start and
  // end (which is the reversed order compared to definition of a path).

  // Extract starting radius and zenith angle
  const Numeric r_start  = ppath.pos(imax,0);
  const Numeric za_start = ppath.los(imax,0);

  // More asserts, checking not at any end point of grid and looks out etc.
  assert( !( ppath.gp_p[imax].idx == 0 && za_start > 90 ) );
  assert( !( ppath.gp_p[imax].idx == (p_grid.nelem()-2) && 
                                  ppath.gp_p[imax].fd[0] && za_start <= 90 ) );
  assert( r_start >= r_geoid + z_ground );

  // Determine index of the pressure surface being the lower limit for the
  // pressure grid step of interest.
  const Index ilow = gridpos2gridrange( ppath.gp_p[imax], za_start<=90 );

  // Tangent radius 
  Numeric r_tan = geometrical_tangent_radius( r_start, za_start );

  // Get end/lowest radius of the path step (r_end) and check:
  // if a tangent point is passed
  // if there is an intersection with the ground
  Numeric r_end;
  bool    tanpoint = false, ground = false;
  if( za_start <= 90 )
    r_end = r_geoid + z_grid[ ilow + 1 ];
  else
    {
      // Ground radius
      Numeric r_ground = r_geoid + z_ground;

      // Lowest possible radius for the path step
      Numeric r_lowest  = r_geoid + z_grid[ ilow ];
      bool    ground_is_low = false;
      if( r_ground >= r_lowest )
	{
	  r_lowest = r_ground;
	  ground_is_low = true;
	}

      // Looking down but not crossing r_lower
      if( r_tan >= r_lowest )
	{
	  r_end = r_tan;
	  tanpoint = true;
	}
      else
	{
	  r_end = r_lowest;
	  if( r_ground == r_lowest )
	    ground = true;
	}
    }

  // Create a vector with a minimum set of radi to describe the path step +
  // the start radius.
  //
  Vector rs;
  //
  // Upward, down to next pressure surface or intersection with blackbody 
  // ground.
  if( za_start<=90 || (!tanpoint && !ground) || (ground && blackbody_ground) ) 
    {
      rs.resize(2);
      rs[1] = r_end;
    }
  // Ground reflection, or passing a tangent point.
  else
    {
      // If the starting point is not at the upper pressure surface, that 
      // surface must be put in as last radius.
      if( r_start == r_geoid + z_grid[ ilow + 1 ] )
	rs.resize(3);
      else
	{
	  rs.resize(4);
          rs[3] = r_geoid + z_grid[ ilow + 1 ];
	}
      rs[1] = r_end;
      rs[2] = r_start;
    }
  rs[0] = r_start;

  // Create a vector with the lengths along the path between the points in rs
  // and the number of sub-steps needed to go from one radius in rs to next.
  //
  Vector         ls( rs.nelem() );
  ArrayOfIndex   ns( rs.nelem()-1 );
  Index          np = 0;                // The sum of ns;
  //
  const Numeric l0 = geompath_l_at_r( r_tan, rs[0] );
  // 
  for( Index i=0; i<rs.nelem(); i++ )
    {
      ls[i] = geompath_l_at_r( r_tan, rs[i] ) - l0;
      if( i > 0 )
	{
	  ns[i-1] = Index( ceil( (ls[i]-ls[i-1]) / lmax ) );
	  np   += ns[i];
	}
    }

  // Re-allocate ppath for return results
  ppath_init_structure(  ppath, atmosphere_dim, np );

  // Go from one radius to next and calculate radii, zenith angles etc.
  //
  np = 0;   // Now used as counter for points done
  //
  for( Index i=0; i<ns.nelem(); i++ )
    {
      Numeric dl   = ( ls[i+1] - ls[i] ) / ns[i]; 
      Numeric dz   = z_grid[ilow+1] - z_grid[ilow];
      Numeric a_za = za_start;   // A zenith angle valid for the path on the
                                 // present side of the tangent point
      for( Index j=1; j<=ns[i]; j++ )
	{
	  ppath.pos(np,0)      = geompath_r_at_l( r_tan, ls[i]+dl*j );
	  ppath.z[np]          = ppath.pos(np,0) - r_geoid;
          ppath.gp_p[np].idx   = ilow;
	  ppath.gp_p[np].fd[0] = ( ppath.z[np] - z_grid[ilow] ) / dz;
	  ppath.gp_p[np].fd[1] = 1 - ppath.gp_p[np].fd[0];
          ppath.l_step[np]     = dl;
          ppath.los(np,0)      = geompath_za_at_r( r_tan, a_za, 
                                                             ppath.pos(np,0) );
	  np ++;
	}
      
      // Handle tangent points
      if( i==0 )
	{
	  if( tanpoint )
	    {
	      ppath.symmetry   = 1;
	      ppath.i_symmetry = np - 1;
              ppath.tan_pos.resize(1);
	      ppath.tan_pos[0] = r_tan;
	      a_za = 80;   // Any value below 90 degrees is OK
	    }
	  else if( ground )
	    {
	      // Re-calculate last zenith angle
	      angle_after_ground_1D( ppath.los(np-1,0) );
	      a_za = ppath.los(np-1,0);

	      if( blackbody_ground )
		ppath_set_background( ppath, 2 );
	      else
		{
		  ppath.symmetry   = 1;
		  ppath.i_symmetry = np - 1;
		  ppath.ground     = 1;
		  ppath.i_ground   = np - 1;
		}
	    }
	}
    }

}

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////

/*!
  \file   ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   Sun May  5  2002
  
  \brief  Functions releated to calculation of propagation paths.
  
  Functions to determine propagation paths for different atmospheric
  dimensionalities.

  The term propagation path is here shortened to ppath.
*/



///////////////////////////////////////////////////////////////////////////////
//   External declarations
///////////////////////////////////////////////////////////////////////////////

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



///////////////////////////////////////////////////////////////////////////////
//   Functions related to geometrical propagation paths.
///////////////////////////////////////////////////////////////////////////////

//// geometrical_tangent_radius ///////////////////////////////////////////////
/**
   Calculates the radius of the geometrical tangent point.

   The absolute value of the zenith angle must be >= 90 degrees.

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
  assert( fabs(za) >= 90 );
  assert( fabs(za) <= 180 );
  return r * sin( DEG2RAD * fabs(za) );
}



//// geompath_za_at_r /////////////////////////////////////////////////////////
/**
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
       const Numeric&   r0,
       const Numeric&   za0,
       const Numeric&   r )
{
  Numeric za = fabs(za0);
  assert( r > 0 );
  assert( za <= 180 );
  assert( (za<90&&r>r0) || (za==90&&r>=r0) || (za>90&&r0>r) );
  if( za > 90 )
    assert( geometrical_tangent_radius( r0, za ) <= r );
  za = RAD2DEG * asin( (r0/r) * sin( DEG2RAD * za ) );
  if( fabs(za0) > 90 )
    za = 180 - za;
  if( za0 < 0 )
    za = -za;
  return za;
}



//// geompath_lat_at_za ///////////////////////////////////////////////////////
/**
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





///////////////////////////////////////////////////////////////////////////////
//   Angle(s) after a ground reflection.
///////////////////////////////////////////////////////////////////////////////

//// angle_after_ground_1D ////////////////////////////////////////////////////
/**
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





///////////////////////////////////////////////////////////////////////////////
//   Functions operating on the Ppath structure
///////////////////////////////////////////////////////////////////////////////

//// ppath_init_structure /////////////////////////////////////////////////////
/**
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



//// ppath_set_background /////////////////////////////////////////////////////
/**
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



//// ppath_what_background ////////////////////////////////////////////////////
/**
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





///////////////////////////////////////////////////////////////////////////////
//   Help functions to ppathCalc
///////////////////////////////////////////////////////////////////////////////

//// ppath_start_stepping /////////////////////////////////////////////////////
/**
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
   a detailed error message is given. Not allowed cases are:
   \begin{verbatim}
      1. The sensor is placed below ground level.
      2. For 2D and 3D, the path leaves the model atmosphere at a latitude or
         longitude end face.
      3. For 2D and 3D, the path is totally outside the atmosphere and the 
         latitude and longitude of the tangent point is outside the range of
         the corresponding grids.
   \end{verbatim}

   All input variables are identical with the WSV with the same name.
   The output variable is here called ppath for simplicity, but is in fact
   *ppath_partial*.

   \param   ppath             Output: A Ppath structure.
   \param   atmosphere_dim    The atmospheric dimensionality.
   \param   p_grid            The pressure grid.
   \param   lat_grid          The latitude grid.
   \param   lon_grid          The longitude grid.
   \param   z_field           The field of geometrical altitudes.
   \param   r_geoid           The geoid radius.
   \parm    z_ground          Ground altitude.
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
	      ppath.los(0,0) = geompath_za_at_r( a_pos[0], a_los[0], 
                                                   r_geoid(0,0) + ppath.z[0] );
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



///////////////////////////////////////////////////////////////////////////////
//   Core functions for *ppath_step* functions
///////////////////////////////////////////////////////////////////////////////

//// ppath_step_1d_geom ///////////////////////////////////////////////////////

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
  assert( ppath.gp_p[imax].idx <= (p_grid.nelem()-1) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] < 1 );

  // Start and end point mean here the point where the calculations start and
  // end (which is the reversed order compared to definition of a path).

  // Extract starting radius and zenith angle
  const Numeric r_start  = ppath.pos(imax,0);
  const Numeric za_start = ppath.los(imax,0);

  // More asserts (checking not at end point of grid and looks out)
  assert( !( ppath.gp_p[imax].idx == 0 && za_start > 90 ) );
  assert( !( ppath.gp_p[imax].idx == (p_grid.nelem()-1) && za_start <= 90 ) );

  // Tangent altitude, if downward obervation
  Numeric r_tan;
  if( za_start > 90 )
    r_tan = geometrical_tangent_radius( r_start, za_start );

  // Create a vector with a minimum set of radii 
  //
  Vector rs;
  //
  // Upward
  if( za_start >= 90 )
    {
      rs.resize(2);
      rs[0] = r_start;
      rs[1] = r_geoid + z_grid[ ppath.gp_p[imax].idx + 1 ];
    }
  // Looking down but not crossing the pressure surface below
  else if( r_tan >= ( r_geoid + z_grid[ppath.gp_p[imax].idx-1] ) )
    {
    }


  // Determine end radius 
  Numeric r_end;
  if( za_start >= 90 || r_tan >= ( r_geoid + z_grid[ppath.gp_p[imax].idx] ) )
    r_end = r_geoid + z_grid[ ppath.gp_p[imax].idx + 1 ];
  else if( r_tan >= ( r_geoid + z_grid[ppath.gp_p[imax].idx-1] ) )
    r_end = r_geoid + z_grid[ ppath.gp_p[imax].idx ];
  else
    r_end = r_geoid + z_grid[ ppath.gp_p[imax].idx - 1 ];

}

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

#include <math.h>
#include <stdexcept>
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "mystring.h"
#include "logic.h"
#include "interpolation.h"
#include "ppath.h"

extern const Numeric DEG2RAD;
extern const Numeric RAD2DEG;



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
  Numeric za;
  za = RAD2DEG * asin( ppc / r );
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
  assert( fabs(za) <= 180 );
  assert( (za0>=0&&za>=0) || (za0<0&&za<0) );
  return lat0 + za0 - za;
}



//! geompath_l_at_r
/*!
   Calculates the length from the tangent point for the given radius.

   The tangent point is either real or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

   \return         Length along the path from the tangent point.
   \param   ppc    Propagation path constant.
   \param   r      Radius of the point of concern.

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
   Calculates the radius for a distance from the tangent point.

   The tangent point is either rwal or imaginary depending on the zenith
   angle of the sensor. See geometrical_tangent_radius.

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
  assert( ( za0 >= 0 && lat >= lat0 )  ||  ( za0 < 0 && lat < lat0 ) );

  // Zenith angle at the new latitude
  const Numeric za = za0 + lat0 -lat;

  return geompath_r_at_za( ppc, za );
}






/*===========================================================================
  === Functions related to slope and tilt of the ground and pressure surfaces
  ===========================================================================*/

//! psurface_slope_2D
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
   \param   gp          Grid position for the position of interest
   \param   upwards     

   \author Patrick Eriksson
   \date   2002-06-03
*/
Numeric psurface_slope_2D(
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



//! psurface_tilt_2D
/*!
   Calculates the angular tilt of the ground or a pressure surface for 2D.

   Note that the tilt value is a local value. The tilt for a constant
   slope value, is different for different radiuses.

   \return        The angular tilt.
   \param    r    The radius for the surface at the point of interest.
   \param    c    The radial slope, as returned by psurface_slope_2D.

   \author Patrick Eriksson
   \date   2002-06-03
*/
Numeric psurface_tilt_2D(
        const Numeric&   r,
        const Numeric&   c )
{
  // The tilt (in radians) is c/r if c is converted to m/radian. So we get
  // conversion RAD2DEG twice
  return   RAD2DEG * RAD2DEG * c / r;
}



//! is_los_downwards_2D
/*!
   Determines if a line-of-sight is downwards compared to the angular tilt
   of the ground or a pressure surface.

   For example, this function can be used to determine if the line-of-sight
   goes into the ground for a starting point exactly on the ground radius.
  
   As the radius of the ground and pressure surfaces varies as a function of
   latitude, it is not clear if a zenith angle of 90 is above or below e.g.
   the ground.
 
   \return 
   \param   za     Zenith angle of line-of-sight.
   \param   tilt   Angular tilt of the ground or the pressure surface (as
                   returned by psurface_tilt_2D)

   \author Patrick Eriksson
   \date   2002-06-03
 */
bool is_los_downwards_2D( 
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





/*===========================================================================
  === Angle(s) after a ground reflection
  ===========================================================================*/

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



//! angle_after_ground_2D 
/*!
   Calculates the LOS after a ground reflection for 1D.

   The zenith angle before the reflection msut be > 90 degrees.

   \param   za     Output: Zenith angle.

   \author Patrick Eriksson
   \date   2002-05-17
*/
void angle_after_ground_2D( Numeric& za, const Numeric& tilt )
{
  assert( is_los_downwards_2D( za, tilt ) );
  
  const Numeric za_new = 180 - fabs(za) - 2 * tilt;

  if( za >= 0 )
    { za = za_new; }
  else
    { za = -za_new; }
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

   The case numbers are:                    <br>
      0. Not yet set.                       <br>
      1. Space.                             <br>
      2. Blackbody ground.                  <br>
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
   \param   ppath            A Ppath structure.

   \author Patrick Eriksson
   \date   2002-05-17
*/
Index ppath_what_background( 
	      Ppath&      ppath )
{
  if( ppath.background == "" )
    { return 0; }
  else if( ppath.background == "space" )
    { return 1; }
  else if( ppath.background == "blackbody ground" )
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





/*===========================================================================
  === Help functions to ppathCalc
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
      // Radius for the ground
      const Numeric r_ground = r_geoid(0,0) + z_ground(0,0);

      // Radius for the top of the atmosphere
      const Numeric r_top = r_geoid(0,0) + z_field(np-1,0,0);

      out2 << "  sensor altitude      : " << (a_pos[0]-r_geoid(0,0))/1e3 
	   << " km\n";

      // The only forbidden case here is that the sensor is below the ground
      if( a_pos[0] < r_ground )
	{
          ostringstream os;
          os << "The sensor is placed " 
             << (r_ground - a_pos[0])/1e3 << " km below ground level.\n"
             << "The sensor must be above the ground.";
	  throw runtime_error(os.str());
	}

      // If downwards, calculate geometrical tangent position
      Vector geom_tan_pos(0);
      if( a_los[0] >= 90 )
	{
	  geom_tan_pos.resize(2);
	  geom_tan_pos[0] = geometrical_ppc( a_pos[0], a_los[0] );
          geom_tan_pos[1] = geompath_lat_at_za( a_los[0], 0, 90 );
	  out2 << "  tangent radius       : " << geom_tan_pos[0]/1e3 <<" km\n";
	  out2 << "  tangent latitude     : " << geom_tan_pos[1] << "\n";
	  out2 << "  tangent altitude     : " 
	       << (geom_tan_pos[0]-r_geoid(0,0))/1e3 << " km\n";
	}


      // The sensor is inside the model atmosphere, 1D ------------------------
      if( a_pos[0] < r_top )
	{
	  // Put some values into ppath. Use these values below (instead of
	  // a_pos and a_los) as they can be modified on the way.
          ppath.pos(0,0) = a_pos[0];
          ppath.los(0,0) = a_los[0];
	  ppath.pos(0,1) = 0; 
	  ppath.z[0]     = ppath.pos(0,0) - r_geoid(0,0);
     
	  // Is the sensor on the ground looking down?
	  if( ppath.pos(0,0) == r_ground  &&  ppath.los(0,0) > 90 )
	    {
	      angle_after_ground_1D( ppath.los(0,0) );
	      ppath.ground = 1;
	      ppath.i_ground = 0;
	      if( blackbody_ground )
		{ ppath_set_background( ppath, 2 ); }
	    }

	  // Check sensor position with respect to cloud box, if activated
          // and not background set to blackbody ground
	  if( cloudbox_on )
	    {
	      // Is the sensor inside the cloud box?
	      if( ppath.z[0] > z_field(cloudbox_limits[0],0,0)  && 
		                 ppath.z[0] < z_field(cloudbox_limits[1],0,0) )
		{ ppath_set_background( ppath, 4 ); }

	      // Is the sensor on the surface of cloud box and looks into box?
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
	  if( fabs(a_los[0]) <= 90 )
	      throw runtime_error("When the sensor is placed outside the model"
                         " atmosphere, upward observations are not allowed." );

	  // We can here set the path constant, that equals the radius of the
	  // geometrical tangent point.
	  ppath.constant = geom_tan_pos[0];
 
	  // Path is above the atmosphere
	  if( a_los[0] <= 90  ||  ppath.constant  >= r_top )
	    {
	      ppath_init_structure(  ppath, atmosphere_dim, 0 );
	      ppath_set_background( ppath, 1 );
	      out2 << "  --- WARNING ---, path is totally outside of the "
		   << "model atmosphere";
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

      // Get grid position for the end point, if there is one.
      if( ppath.np == 1 )
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
      Numeric rv_geoid, rv_ground;
      ArrayOfGridPos gp(1);
      Matrix itw(1,2);
      //
      if( a_pos[1] >= lat_grid[0]  &&  a_pos[1] <= lat_grid[nlat-1] )
	{
	  gridpos( gp, lat_grid, Vector(1,a_pos[1]) );
	  interpweights( itw, gp );

	  Vector v_rgeoid(1), v_zground(1), v_ztop(1);
	  interp( v_rgeoid, itw, r_geoid(Range(joker),0), gp );
	  interp( v_zground, itw, z_ground(Range(joker),0), gp );
	  interp( v_ztop, itw, z_field(np-1,Range(joker),0), gp );

	  rv_geoid  = v_rgeoid[0];
	  rv_ground = rv_geoid + v_zground[0];

	  out2 << "  sensor altitude      : " << (a_pos[0]-rv_geoid)/1e3 
	       << " km\n";

	  if( a_pos[0] < ( v_rgeoid[0] + v_ztop[0] ) )
	    { is_inside = 1; }
	}

      // If downwards, calculate geometrical tangent position. If the tangent
      // point is inside the covered latitude range, calculate also the 
      // geometrical altitude of the tangent point and the top of atmosphere.
      Vector geom_tan_pos(0);
      Numeric geom_tan_z, geom_tan_atmtop;
      if( a_los[0] >= 90 )
	{
	  geom_tan_pos.resize(2);
	  geom_tan_pos[0] = geometrical_ppc( a_pos[0], a_los[0] );
	  if( a_los[0] > 0 )
	    { geom_tan_pos[1] = geompath_lat_at_za( a_los[0], a_pos[1], 90 ); }
	  else
	    { geom_tan_pos[1] = geompath_lat_at_za( a_los[0], a_pos[1], -90 );}
	  out2 << "  tangent radius       : " << geom_tan_pos[0]/1e3<<" km\n";
	  out2 << "  tangent latitude     : " << geom_tan_pos[1] << "\n";
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
	      out2 << "  tangent altitude     : " << geom_tan_z/1e3 << " km\n";
	    }
	}

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

	  
	  // Put some values into ppath. Use these values below (instead of
	  // a_pos and a_los) as they can be modified on the way.
          ppath.pos(0,0) = a_pos[0];
          ppath.pos(0,1) = a_pos[1];
          ppath.los(0,0) = a_los[0];
	  ppath.z[0]     = ppath.pos(0,0) - rv_geoid;
     
	  // Is the sensor on the ground looking down?
	  if( ppath.pos(0,0) == rv_ground )
	    {
	      // Calculate radial slope of the ground
	      const Numeric rslope = psurface_slope_2D( lat_grid, 
                        r_geoid(Range(joker),0), z_ground(Range(joker),0), 
                                                  gp[0], ppath.los(0,0) >= 0 );

	      // Calculate angular tilt of the ground
	      const Numeric atilt = psurface_tilt_2D( rv_ground, rslope);

	      // Are we looking down into the ground?
	      if( is_los_downwards_2D( ppath.los(0,0), atilt ) )
		{
		  angle_after_ground_2D( ppath.los(0,0), atilt );
		  ppath.ground = 1;
		  ppath.i_ground = 0;
		  if( blackbody_ground )
		    { ppath_set_background( ppath, 2 ); }
		}
	    }

	  // Check sensor position with respect to cloud box, if activated
          // and not background set to blackbody ground
	  if( cloudbox_on )
	    {
	      // Are we inside the interesting latitude range
	      if( ppath.pos(0,1) >= lat_grid[cloudbox_limits[2]]  &&
		               ppath.pos(0,1) <= lat_grid[cloudbox_limits[3]] )
		{
		  // Calculate the lower and upper altitude radius limit for
		  // the cloud box at the latitude of the sensor
	          Vector v_zlim(1);
		  interp( v_zlim, itw, 
                              z_field(cloudbox_limits[0],Range(joker),0), gp );
		  Numeric rv_low = rv_geoid + v_zlim[0];
		  interp( v_zlim, itw, 
                              z_field(cloudbox_limits[1],Range(joker),0), gp );
		  Numeric rv_upp = rv_geoid + v_zlim[0];

		  // Is the sensor inside the cloud box?
		  if( ppath.pos(0,0) > rv_low  &&  ppath.pos(0,0) < rv_upp )
		    { ppath_set_background( ppath, 4 ); }

		  // Is the sensor on the surface of cloud box and looks 
		  // into box?
		  else if( ppath.pos(0,0) == rv_low )
		    {
		      // Calculate radial slope etc. of the pressure surface at
		      // the lower cloud box limit
		      Numeric rslope = psurface_slope_2D( lat_grid, 
                              r_geoid(Range(joker),0), 
                                 z_field(cloudbox_limits[0],Range(joker),0), 
                                                  gp[0], ppath.los(0,0) >= 0 );
		      Numeric atilt = psurface_tilt_2D( rv_ground, rslope );
		      // Note the ! before the function command
		      if( !is_los_downwards_2D( ppath.los(0,0), atilt ) )
			{ ppath_set_background( ppath, 3 ); }
		    }
		  //
		  else if( ppath.pos(0,0) == rv_upp )
		    {
		      // Calculate radial slope etc. of the pressure surface at
		      // the upper cloud box limit
		      Numeric rslope = psurface_slope_2D( lat_grid, 
                              r_geoid(Range(joker),0), 
                                 z_field(cloudbox_limits[1],Range(joker),0), 
                                                  gp[0], ppath.los(0,0) >= 0 );
		      Numeric atilt = psurface_tilt_2D( rv_ground, rslope );
		      if( is_los_downwards_2D( ppath.los(0,0), atilt ) )
			{ ppath_set_background( ppath, 3 ); }
		    }
		}
	    }
	}

      // The sensor is outside the model atmosphere, 2D -----------------------
      else
	{
	  // Upward observations are not allowed here
	  if( fabs(a_los[0]) <= 90 )
	      throw runtime_error("When the sensor is placed outside the model"
                         " atmosphere, upward observations are not allowed." );
	  
	  // We can here set the path constant, that equals the radius of the
	  // geometrical tangent point.
	  ppath.constant = geom_tan_pos[0];

	  // Handle cases when the sensor appears to look the wrong way
	  if( ( a_pos[1] <= lat_grid[0]  &&  a_los[0] <= 0 )  || 
	                  ( a_pos[1] >= lat_grid[nlat-1]  &&  a_los[0] >= 0 ) )
	    {
	      ostringstream os;
	      os << "The sensor is outside (or at the limit) of the model "
		 << "atmosphere but looks in the wrong\ndirection (wrong sign "
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
		     << "looks in the " << sc << " latitude end face.\n"
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
		     << "gives a propagation path that goes above\nthe model "
		     << "atmosphere, with a tangent point outside the covered "
		     << "latitude range.\nThe latitude of the tangent point "
		     << "is " << geom_tan_pos[1] << " degress.";
		  throw runtime_error( os.str() );
		}
	    }

	  // That should be all needed checks. We know now that the path is
	  // either totally outside the atmosphere, with a tangent point 
	  // inside lat_grid, or enters the atmosphere from the top 
	  // somewhere inside lat_grid. In the latter case we need to
	  // determine the latitude of the entrence point.
	  
	  // Path is above the atmosphere
	  if( geom_tan_z >= geom_tan_atmtop )
	    {
	      ppath_init_structure(  ppath, atmosphere_dim, 0 );
	      ppath_set_background( ppath, 1 );
	      out2 << "  --- WARNING ---: path is totally outside of the "
		   << "model atmosphere\n";
	    }
	}      

      // Get grid position for the end point, if there is one.
      if( ppath.np == 1 )
	{ 
	  gridpos( ppath.gp_p, z_field(Range(joker),0,0), ppath.z ); 
	  gridpos( ppath.gp_lat, lat_grid, ppath.pos(Range(joker),1) ); 
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
  === Core functions for *ppath_step* functions
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
   a 1D atmosphere. For example, z_grid is z_field(:,0,0).

   For more information read the chapter on propagation paths in AUG.

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
void ppath_step_geom_1d(
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
  assert( p_grid.nelem() == z_grid.nelem() );
  assert( p_grid.nelem() >= 2 );
  assert( r_geoid > 0 );
  assert( is_bool( blackbody_ground ) );
  //
  assert( ppath.dim == atmosphere_dim );
  assert( ppath.np >= 1 );
  assert( ppath.gp_p[imax].idx >= 0 );
  assert( ppath.gp_p[imax].idx <= (p_grid.nelem()-2) );
  assert( ppath.gp_p[imax].fd[0] >= 0 );
  assert( ppath.gp_p[imax].fd[0] <= 1 );

  // Start and end point mean here the point where the calculations start and
  // end (which is the reversed order compared to definition of a path).

  // Extract starting radius, zenith angle and latitude
  const Numeric r_start   = ppath.pos(imax,0);
  const Numeric lat_start = ppath.pos(imax,1);
  const Numeric za_start  = ppath.los(imax,0);

  // More asserts, checking not at any end point of grid and looks out etc.
  assert( r_start >= r_geoid + z_ground );
  assert( r_start <= r_geoid + last( z_grid ) );
  assert( !( is_gridpos_at_index_i( ppath.gp_p[imax], 0 ) &&  za_start > 90 ));
  assert( !( is_gridpos_at_index_i( ppath.gp_p[imax], p_grid.nelem()-1 )  && 
                                                            za_start <= 90 ) );

  // If the field "constant" is negative, this is the first call of the
  // function and the path constant shall be calculated.
  Numeric ppc;
  if( ppath.constant < 0 )
    { ppc = geometrical_ppc( r_start, za_start ); }
  else
    { ppc = ppath.constant; }


  // Determine index of the pressure surface being the lower limit for the
  // pressure grid step of interest.
  const Index ilow = gridpos2gridrange( ppath.gp_p[imax], za_start<=90 );

  // Get end/lowest radius of the path step (r_end) and check:
  // if a tangent point is passed
  // if there is an intersection with the ground
  Numeric r_end;
  bool    tanpoint = false, ground = false;
  if( za_start <= 90 )
    { r_end = r_geoid + z_grid[ ilow + 1 ]; }
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

      // Looking down but not crossing r_lower. Note that ppc
      // here equals the tangent radius.
      if( ppc >= r_lowest )
	{
	  r_end = ppc;
	  tanpoint = true;
	}
      else
	{
	  r_end = r_lowest;
	  if( r_ground == r_lowest )
	    { ground = true; }
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
	{ rs.resize(3); }
      else
	{
	  rs.resize(4);
          rs[3] = r_geoid + z_grid[ ilow + 1 ];
	}
      rs[1] = r_end;
      rs[2] = r_start;
    }
  rs[0] = r_start;

  // Create a vector (ls) with the length between the points in rs and the 
  // tangent point, and another vector (ns) with the number of sub-steps needed
  // to go from one radius in rs to next. The last value in ns is one, which
  // corresponds to the last point.
  //
  const Index    nr = rs.nelem();
  Vector         ls( nr );
  ArrayOfIndex   ns( nr ); 
  Index          np = 0;           // The sum of ns. The number one corresponds
  //                                  to the starting point.
  for( Index i=0; i<nr; i++ )
    {
      ls[i] = geompath_l_at_r( ppc, rs[i] );
      if( i > 0 )
	{
	  if( lmax > 0 )
	    {
	      // The absolute value of the length distance is needed here
	      ns[i-1] = Index( ceil( fabs(ls[i]-ls[i-1]) / lmax ) );
	    }
	  else
	    { ns[i-1] = 1; }
	  np   += ns[i-1];
	}
    }
  // Last point
  ns[nr-1] = 1;
  np++;

  // If np is <= 1, something has gone wrong
  assert( np > 1 ); 

  // Re-allocate ppath for return results and put in some variables
  ppath_init_structure(  ppath, atmosphere_dim, np );
  //
  ppath.method     = "1D basic geometrical";
  ppath.refraction = 0;
  ppath.constant   = ppc;
  //
  np = 0;   // Now used as counter for points done
  //
  // Any angle valid for the path on the present side of the tangent point. 
  Numeric a_za = za_start; 
  //
  // A reference zenith angle and latitude for the calculation of latitudes.
  // This solution is needed as there is a step in the zenith angles at ground
  // reflections.
  Numeric za_lat0 = za_start;
  Numeric lat0    = lat_start;
  //
  // Go from one radius to next and calculate radii, zenith angles etc. and
  // fill the return structure.
  //
  for( Index i=0; i<nr; i++ )
    {
      // Note that dl is allowed to be negative
      Numeric dl;
      if( i < ( nr-1 ) )
	{ dl = ( ls[i+1] - ls[i] ) / ns[i]; }
      else
	{ dl = 0; }

      Numeric dz   = z_grid[ilow+1] - z_grid[ilow];

      for( Index j=0; j<ns[i]; j++ )
	{
	  if( j == 0 )
	    { ppath.pos(np,0)  = rs[i]; }
	  else
	    { ppath.pos(np,0)  = geompath_r_at_l( ppc, ls[i]+dl*j ); }
	  ppath.z[np]      = ppath.pos(np,0) - r_geoid;
	  if( i < ( nr-1 ) )
	    { ppath.l_step[np] = fabs( dl ); }
          ppath.los(np,0)  = geompath_za_at_r( ppc, a_za, ppath.pos(np,0) );
	  a_za             = ppath.los(np,0);
          ppath.pos(np,1)  = geompath_lat_at_za( za_lat0, lat0, a_za );

	  // Calculate grid position. Things to consider:
	  // Rounding errors can give fd-values below 0 or above 1, and the
	  // last point must have a fd-value of exactly 0 and 1, if not ending
	  // at a blackbody ground.
          ppath.gp_p[np].idx   = ilow;
	  ppath.gp_p[np].fd[0] = ( ppath.z[np] - z_grid[ilow] ) / dz;
	  ppath.gp_p[np].fd[1] = 1 - ppath.gp_p[np].fd[0];
          gridpos_check_fd( ppath.gp_p[np] );
	  if( i == (nr-1)  &&  j == (ns[i]-1)  && 
                                              !( ground && blackbody_ground ) )
	    { gridpos_force_end_fd( ppath.gp_p[np] ); }

	  np ++;
	}

      // Handle tangent points and ground reflections
      if( i == 0 )
	{
	  if( tanpoint )
	    {
	      ppath.symmetry   = 1;
	      ppath.i_symmetry = np - 1;
              ppath.tan_pos.resize(2);
	      ppath.tan_pos[0] = ppc;
	      ppath.tan_pos[1] = geompath_lat_at_za( za_start, lat_start, 90 );
	      a_za = 80;   // Any value below 90 degrees is OK
	    }
	  else if( ground )
	    {
	      // Only a_za must be set to the zenith angle after ground 
              // reflection. As the zenith angle shall be valid for the path
	      // when leaving the point, the value in ppath.los is correct.
	      a_za    = ppath.los(np-1,0);
	      angle_after_ground_1D( a_za );
	      za_lat0 = a_za;
	      lat0    = ppath.pos(np-1,1);

	      if( blackbody_ground )
		{ ppath_set_background( ppath, 2 ); }
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

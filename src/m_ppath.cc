/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                            
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
  ===  File description 
  ===========================================================================*/

/*!
  \file   m_ppath.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-08 

  \brief  Workspace functions releated to propagation paths variables.

  The file includes special functions to set the sensor position and LOS,
  and functions for calculation of propagation paths.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <iostream>
#include <cmath>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! a_losSet
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-30
*/
void a_losSet(
        // WS Output:
              Vector&    a_los,
        // WS Input:
        const Index&     atmosphere_dim,
        // Control Parameters:
        const Numeric&   za,
        const Numeric&   aa )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  if( atmosphere_dim == 1 )
    { a_los.resize(1); }
  else
    {
      a_los.resize(2);
      a_los[1] = aa;
    }
  a_los[0] = za;
}



//! a_posAddGeoidWGS84
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-20
*/
void a_posAddGeoidWGS84(
        // WS Output:
              Vector&    a_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        const Numeric&   latitude_1d,
        const Numeric&   meridian_angle_1d )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_vector_length( "a_pos", a_pos, atmosphere_dim );

  // Use *sensor_posAddGeoidWGS84* to perform the calculations.
  Matrix m(1,a_pos.nelem());
  m(0,Range(joker)) = a_pos;
  sensor_posAddGeoidWGS84( m, atmosphere_dim, latitude_1d, meridian_angle_1d);
  a_pos[0] = m(0,0);
}



//! a_posAddRgeoid
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-20
*/
void a_posAddRgeoid(
        // WS Output:
              Vector&    a_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Matrix&    r_geoid )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_vector_length( "a_pos", a_pos, atmosphere_dim );

  // Use *sensor_posAddRgeoid* to perform the calculations.
  Matrix m(1,a_pos.nelem());
  m(0,Range(joker)) = a_pos;
  sensor_posAddRgeoid( m, atmosphere_dim, lat_grid, lon_grid, r_geoid );
  a_pos[0] = m(0,0);
}



//! a_posSet
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-30
*/
void a_posSet(
        // WS Output:
              Vector&    a_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        // Control Parameters:
        const Numeric&   r_or_z,
        const Numeric&   lat,
        const Numeric&   lon )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  a_pos.resize(atmosphere_dim);
  a_pos[0] = r_or_z;
  if( atmosphere_dim >= 2 )
    { a_pos[1] = lat; }
  if( atmosphere_dim == 3 )
    { a_pos[2] = lon; }
}



//! ppathCalc
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void ppathCalc(
        // WS Output:
              Ppath&          ppath,
	      Ppath&          ppath_step,
        // WS Input:
	const Agenda&         ppath_step_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         a_pos,
        const Vector&         a_los )
{
  // This function is a WSM but it is normally only called from RteCalc. 
  // For that reason, this function does not repeat input checks that are
  // performed in RteCalc, it only performs checks regarding the sensor 
  // position and LOS.

  //--- Check input -----------------------------------------------------------

  // Sensor position and LOS
  //
  chk_vector_length( "a_pos", a_pos, atmosphere_dim );
  chk_if_over_0( "sensor radius", a_pos[0] );
  if( atmosphere_dim == 1 )
    {
      chk_vector_length( "a_los", a_los, 1 );
      chk_if_in_range( "sensor zenith angle", a_los[0], 0, 180 );
    }
  else if( atmosphere_dim == 2 )
    {
      chk_vector_length( "a_los", a_los, 1 );
      chk_if_in_range( "sensor zenith angle", a_los[0], -180, 180 );
    }
  else
    {
      chk_if_in_range( "sensor latitude", a_pos[1], -90, 90 );
      chk_if_in_range( "sensor longitude", a_pos[2], -360, 360 );
      chk_vector_length( "a_los", a_los, 2 );
      chk_if_in_range( "sensor zenith angle", a_los[0], 0, 180 );
      chk_if_in_range( "sensor azimuth angle", a_los[1], -180, 180 );
    }
  
  //--- End: Check input ------------------------------------------------------


  // Some messages
  out2 << "  -------------------------------------\n";
  out2 << "  sensor radius          : " << a_pos[0]/1e3 << " km\n";
  if( atmosphere_dim >= 2 )    
    out2 << "  sensor latitude        : " << a_pos[1] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor longitude       : " << a_pos[2] << "\n";
  out2 << "  sensor zenith angle    : " << a_los[0] << "\n";
  if( atmosphere_dim == 3 )
    out2 << "  sensor azimuth angle   : " << a_los[1] << "\n";


  // Initiate the partial Ppath structure. 
  // The function doing the work sets ppath_step to the point of the path
  // inside the atmosphere closest to the sensor, if the path is at all inside
  // the atmosphere.
  // If the background field is set by the function this flags that there is no
  // path to follow (for example when the sensor is inside the cloud box).
  // The function checks also that the sensor position and LOS give an
  // allowed path.
  //
  ppath_start_stepping( ppath_step, atmosphere_dim, p_grid, lat_grid, 
                        lon_grid, z_field, r_geoid, z_ground,
                        cloudbox_on, cloudbox_limits, a_pos, a_los );

  out2 << "  -------------------------------------\n";

  // Perform propagation path steps until the starting point is found, which
  // is flagged by ppath_step by setting the background field.
  //
  // The results of each step, returned by ppath_step_agenda as a new 
  // ppath_step, are stored as an array of Ppath structures.
  //
  Array<Ppath>   ppath_array;
  ppath_array.push_back( ppath_step );
  // 
  Index   np = ppath_step.np;   // Counter for number of points of the path
  Index   istep = 0;            // Counter for number of steps
  //
  const Index imax_p   = p_grid.nelem() - 1;
  const Index imax_lat = lat_grid.nelem() - 1;
  const Index imax_lon = lon_grid.nelem() - 1;
  //
  
  while( !ppath_what_background( ppath_step ) )
    {

      // Call ppath_step agenda
      //
      istep++;
      out3 << "  path step number     : " << istep << "\n";
      //
      
      // (CE:) Included istep here to execute the agenda silently.
      ppath_step_agenda.execute(istep);
      
      //PpathPrint(ppath_step,"ppath_step");

      // Number of points in returned path step
      const Index n = ppath_step.np;

      // Increase the total number
      np += n - 1;

      // Put new ppath_step in ppath_array
      ppath_array.push_back( ppath_step );

      // Check if the top of the atmosphere is reached
      if( is_gridpos_at_index_i( ppath_step.gp_p[n-1], imax_p ) )
	{ ppath_set_background( ppath_step, 1 ); }

      // Check that path does not exit at a latitude or longitude end face
      if( atmosphere_dim >= 2 )
	{
	  if( is_gridpos_at_index_i( ppath_step.gp_lat[n-1], 0 ) )
	    {
	      ostringstream os;
	      os << "The path exits the atmosphere through the lower latitude"
		 << " end face.\nThe exit point is at an altitude of " 
                 << ppath_step.z[n-1]/1e3 << " km.";
	      throw runtime_error( os.str() );
	    }
	  if( is_gridpos_at_index_i( ppath_step.gp_lat[n-1], imax_lat ) )
	    {
	      ostringstream os;
	      os << "The path exits the atmosphere through the upper latitude"
		 << " end face.\nThe exit point is at an altitude of " 
                 << ppath_step.z[n-1]/1e3 << " km.";
	      throw runtime_error( os.str() );
	    }

	  if( atmosphere_dim == 3 )
	    {
	      if( is_gridpos_at_index_i( ppath_step.gp_lon[n-1], 0 ) )
		{
		  ostringstream os;
		  os << "The path exits the atmosphere through the lower " 
		     << "longitude end face.\nThe exit point is at an altitude"
		     << "of " << ppath_step.z[n-1]/1e3 << " km.";
		  throw runtime_error( os.str() );
		}
	      if( is_gridpos_at_index_i( ppath_step.gp_lon[n-1], imax_lon ))
		{
		  ostringstream os;
		  os << "The path exits the atmosphere through the upper "
		     << "longitude end face.\nThe exit point is at an altitude"
		     << "of " << ppath_step.z[n-1]/1e3 << " km.";
		  throw runtime_error( os.str() );
		}
	    }
	}
      
    
      // Check if there is an intersection with an active cloud box
      if( cloudbox_on )
	{
	  Numeric ipos = Numeric( ppath_step.gp_p[n-1].idx ) + 
                                                    ppath_step.gp_p[n-1].fd[0];
	  if( ipos >= Numeric( cloudbox_limits[0] )  && 
	                                ipos <= Numeric( cloudbox_limits[1] ) )
	    {
	      if( atmosphere_dim == 1 )
		{ ppath_set_background( ppath_step, 3 ); }
	      else
		{
		  ipos = Numeric( ppath_step.gp_lat[n-1].idx ) + 
                                                  ppath_step.gp_lat[n-1].fd[0];
		  if( ipos >= Numeric( cloudbox_limits[2] )  && 
	                                ipos <= Numeric( cloudbox_limits[3] ) )
		    {
		      if( atmosphere_dim == 2 )
			{ ppath_set_background( ppath_step, 3 ); }
		      else
			{
			  ipos = Numeric( ppath_step.gp_lon[n-1].idx ) + 
                                                  ppath_step.gp_lon[n-1].fd[0];
			  if( ipos >= Numeric( cloudbox_limits[4] )  && 
	                                ipos <= Numeric( cloudbox_limits[5] ) )
			    { ppath_set_background( ppath_step, 3 ); } 
			}
		    }
		}
	    }
	}
    } // End path steps
  
 
  // Combine all structures in ppath_array to form the return Ppath structure.
  //
  ppath_init_structure( ppath, atmosphere_dim, np );
  //
  np = 0;   // Now used as counter for points moved to ppath
  //
  for( Index i=0; i<ppath_array.nelem(); i++ )
    {
      // For the first structure, the first point shall be included, but the
      // first structure can also be empty. 
      // For later structures, the first point shall not be included, but
      // there will always be at least two points.
      // Only the first structure can be empty.

      Index n = ppath_array[i].np;

      if( n )
	{
	  // First index to include
	  Index i1 = 1;
	  if( i == 0 )
	    { i1 = 0; }

	  // Vectors and matrices that can be handled by ranges.
	  ppath.z[ Range(np,n-i1) ] = ppath_array[i].z[ Range(i1,n-i1) ];
	  ppath.pos( Range(np,n-i1), Range(joker) ) = 
                            ppath_array[i].pos( Range(i1,n-i1), Range(joker) );
	  ppath.los( Range(np,n-i1), Range(joker) ) = 
	                    ppath_array[i].los( Range(i1,n-i1), Range(joker) );

	  // For i==1, there is no defined l_step. For higher i, all 
	  // values in l_step shall be copied.
	  if( i > 0 )
	    { ppath.l_step[ Range(np-1,n-1) ] = ppath_array[i].l_step; }

	  // Grid positions must be handled by a loop
	  for( Index j=i1; j<n; j++ )
	    { ppath.gp_p[np+j-i1] = ppath_array[i].gp_p[j]; }
	  if( atmosphere_dim >= 2 )
	    {
	      for( Index j=i1; j<n; j++ )
		{ ppath.gp_lat[np+j-i1] = ppath_array[i].gp_lat[j]; }
	    }
	  if( atmosphere_dim == 3 )
	    {
	      for( Index j=i1; j<n; j++ )
		{ ppath.gp_lon[np+j-i1] = ppath_array[i].gp_lon[j]; }
	    }

	  // Fields just set once
	  if( ppath_array[i].tan_pos.nelem() )
	    {
	      ppath.tan_pos.resize( ppath_array[i].tan_pos.nelem() );
	      ppath.tan_pos               = ppath_array[i].tan_pos; 
	    }
	  if( ppath_array[i].geom_tan_pos.nelem() )
	    {
	      ppath.geom_tan_pos.resize( ppath_array[i].tan_pos.nelem() );
	      ppath.geom_tan_pos          = ppath_array[i].geom_tan_pos; 
	    }

	  // Increase number of points done
	  np += n - i1;
	 
	}
    }  
  ppath.method     = ppath_step.method;
  ppath.refraction = ppath_step.refraction;
  ppath.constant   = ppath_step.constant;
  ppath.background = ppath_step.background;
}



//! ppath_stepGeometric
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void ppath_stepGeometric(
        // WS Output:
              Ppath&     ppath_step,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Tensor3&   z_field,
        const Matrix&    r_geoid,
        const Matrix&    z_ground )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here. This commented in the on-line information.

  // Note that lmax is here set to -1.

  if( atmosphere_dim == 1 )
    { ppath_step_geom_1d( ppath_step, p_grid, z_field(Range(joker),0,0), 
                                             r_geoid(0,0), z_ground(0,0), -1 );
    }
  else if( atmosphere_dim == 2 )
    { ppath_step_geom_2d( ppath_step, p_grid, lat_grid,
             z_field(Range(joker),Range(joker),0), r_geoid(Range(joker),0), 
                                                z_ground(Range(joker),0), -1 );
    }
  else
    {
      throw runtime_error( "3D propagation path steps are not yet handled." );
    }
}



//! ppath_stepGeometricWithLmax
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void ppath_stepGeometricWithLmax(
        // WS Output:
              Ppath&     ppath_step,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Tensor3&   z_field,
        const Matrix&    r_geoid,
        const Matrix&    z_ground,
        // Control Parameters:
        const Numeric&   lmax )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here. This commented in the on-line information for ppath_stepGeometric.

  if( atmosphere_dim == 1 )
    { ppath_step_geom_1d( ppath_step, p_grid, z_field(Range(joker),0,0), 
                                           r_geoid(0,0), z_ground(0,0), lmax );
    }
  else if( atmosphere_dim == 2 )
    { ppath_step_geom_2d( ppath_step, p_grid, lat_grid,
             z_field(Range(joker),Range(joker),0), r_geoid(Range(joker),0), 
                                              z_ground(Range(joker),0), lmax );
    }
  else
    {
      throw runtime_error( "3D propagation path steps are not yet handled." );
    }
}



//! ppath_stepRefractionEuler
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-11-14
*/
void ppath_stepRefractionEuler(
        // WS Output:
              Ppath&     ppath_step,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Tensor3&   z_field,
        const Tensor3&   t_field,
        const Matrix&    r_geoid,
        const Matrix&    z_ground,
        // Control Parameters:
	const Numeric&   lraytrace )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here. This commented in the on-line information.

  // Note that lmax is here set to -1.

  if( atmosphere_dim == 1 )
    { ppath_step_refr_1d( ppath_step, p_grid, z_field(Range(joker),0,0), 
                t_field(Range(joker),0,0), r_geoid(0,0), z_ground(0,0), 
                                               "linear_euler", lraytrace, -1 );
    }
  else if( atmosphere_dim == 2 )
    { ppath_step_refr_2d( ppath_step, p_grid, lat_grid,
                        z_field(Range(joker),Range(joker),0), 
                        t_field(Range(joker),Range(joker),0), 
                        r_geoid(Range(joker),0), z_ground(Range(joker),0), 
                                               "linear_euler", lraytrace, -1 );
    }
  else
    {
      throw runtime_error( "3D propagation path steps are not yet handled." );
    }
}



//! ppath_stepRefractionEulerWithLmax
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void ppath_stepRefractionEulerWithLmax(
        // WS Output:
              Ppath&     ppath_step,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    p_grid,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Tensor3&   z_field,
        const Tensor3&   t_field,
        const Matrix&    r_geoid,
        const Matrix&    z_ground,
        // Control Parameters:
	const Numeric&   lraytrace,
        const Numeric&   lmax )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here. This commented in the on-line information for ppath_stepGeometric.

  if( atmosphere_dim == 1 )
    { ppath_step_refr_1d( ppath_step, p_grid, z_field(Range(joker),0,0), 
                   t_field(Range(joker),0,0), r_geoid(0,0), z_ground(0,0), 
                                             "linear_euler", lraytrace, lmax );
    }
  else if( atmosphere_dim == 2 )
    { ppath_step_refr_2d( ppath_step, p_grid, lat_grid,
                       z_field(Range(joker),Range(joker),0), 
                       t_field(Range(joker),Range(joker),0), 
                       r_geoid(Range(joker),0), z_ground(Range(joker),0), 
                                             "linear_euler", lraytrace, lmax );
    }
  else
    {
      throw runtime_error( "3D propagation path steps are not yet handled." );
    }
}



//! PpathPrint
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void PpathPrint(
        // WS Generic Input:
        const Ppath&     ppath,
        // WS Generic Input Names:
        const String&    ppath_name )
{
  cout << "  The fields of *" << ppath_name <<"*:\n";
  IndexPrint( ppath.dim, "dim" );
  IndexPrint( ppath.np, "np" );
  IndexPrint( ppath.refraction, "refraction" );
  StringPrint( ppath.method, "method" );
  NumericPrint( ppath.constant, "constant" );
  MatrixPrint( ppath.pos, "pos" );
  VectorPrint( ppath.z, "z" );
  VectorPrint( ppath.l_step, "l_step" );
  ArrayOfGridPosPrint( ppath.gp_p, "gp_p" );
  if( ppath.dim >= 2 )
    ArrayOfGridPosPrint( ppath.gp_lat, "gp_lat" );
  if( ppath.dim == 3 )
    ArrayOfGridPosPrint( ppath.gp_lon, "gp_lon" );
  MatrixPrint( ppath.los, "los" );
  StringPrint( ppath.background, "background" );
  if( ppath.tan_pos.nelem() )
    VectorPrint( ppath.tan_pos, "tan_pos" );
  if( ppath.geom_tan_pos.nelem() )
    VectorPrint( ppath.geom_tan_pos, "geom_tan_pos" );
}



//! sensor_posAddGeoidWGS84
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-16
*/
void sensor_posAddGeoidWGS84(
        // WS Output:
              Matrix&    sensor_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        const Numeric&   latitude_1d,
        const Numeric&   meridian_angle_1d )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_matrix_ncols( "sensor_pos", sensor_pos, atmosphere_dim );

  // Number of positions
  const Index npos = sensor_pos.nrows();
  if( npos == 0 )
    throw runtime_error("The number of positions is 0, must be at least 1.");

  // The function *r_geoidWGS84 is used to get the geoid radius, but some
  // tricks are needed as this a WSF to set *r_geoid* for all crossings of the
  // latitude and longitude grids.
  // For 2D and 3D, the function is always called with the atmospheric 
  // dimensionality set to 2, and the latitudes of concern are put into 
  // the latitude grid. An extra dummy value is needed if there is only one 
  // position in *sensor_pos*.

  if( atmosphere_dim == 1 )
    {
      Vector lats(0);

      // The size of the r-matrix is set inside the function.
      Matrix r;
      r_geoidWGS84( r, 1, lats, Vector(0), latitude_1d, meridian_angle_1d );
      
      // Add the geoid radius to the geometric altitudes
      sensor_pos(Range(joker),0) += r(0,0);
    }

  else
    {
      Vector lats;
      if( npos == 1 )
	{
	  lats.resize(2);
	  lats[0] = sensor_pos(0,1);
	  lats[1] = lats[0] + 1;
	}
      else
	{
	  lats.resize(npos);
          lats = sensor_pos( Range(joker), 1 );
	}

      // The size of the r-matrix is set inside the function.
      Matrix r;
      r_geoidWGS84( r, 2, lats, Vector(0), -999, -999 );
      
      // Add the geoid radius to the geometric altitudes
      for( Index i=0; i<npos; i++ )
	sensor_pos(i,0) += r(i,0);
    }
}



//! sensor_posAddRgeoid
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-06-03
*/
void sensor_posAddRgeoid(
        // WS Output:
              Matrix&    sensor_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Matrix&    r_geoid )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_matrix_ncols( "sensor_pos", sensor_pos, atmosphere_dim );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );

  // Number of positions
  const Index npos = sensor_pos.nrows();
  //
  if( npos == 0 )
    throw runtime_error("The number of positions is 0, must be at least 1.");

  if( atmosphere_dim == 1 )
    { sensor_pos(Range(joker),0) += r_geoid(0,0); }

  else
    {
      // Check that positions in sensor_pos are inside the lat and lon grids
      if( min(sensor_pos(Range(joker),1)) < lat_grid[0]  || 
                             max(sensor_pos(Range(joker),1)) > last(lat_grid) )
	throw runtime_error(
             "You have given a position with a latitude outside *lat_grid*." );
      if( atmosphere_dim == 3 )
	{
	  if( min(sensor_pos(Range(joker),2)) < lon_grid[0]  || 
                            max(sensor_pos(Range(joker),2)) >= last(lon_grid) )
	    throw runtime_error(
            "You have given a position with a longitude outside *lat_grid*." );
	}

      if( atmosphere_dim == 2 )
	{
	  ArrayOfGridPos gp(npos);
	  Matrix itw(npos,2);
	  gridpos( gp, lat_grid, sensor_pos(Range(joker),1) );
	  interpweights( itw, gp );
	  Vector v_rgeoid(npos);
	  interp( v_rgeoid, itw, r_geoid(Range(joker),0), gp );
	  for( Index i=0; i<npos; i++ )
	    { sensor_pos(i,0) += v_rgeoid[i]; } 
	}
      else
	{
	  ArrayOfGridPos gp_lat(npos), gp_lon(npos);
	  Matrix itw(npos,4);
	  gridpos( gp_lat, lat_grid, sensor_pos(Range(joker),1) );
	  gridpos( gp_lon, lon_grid, sensor_pos(Range(joker),2) );
	  interpweights( itw, gp_lat, gp_lon );
	  Vector v_rgeoid(npos);
	  interp( v_rgeoid, itw, r_geoid, gp_lat, gp_lon );
	  for( Index i=0; i<npos; i++ )
	    { sensor_pos(i,0) += v_rgeoid[i]; } 
	}
    }
}




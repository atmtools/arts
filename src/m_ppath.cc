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

#include <cmath>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "math_funcs.h"
#include "messages.h"
#include "ppath.h"
#include "special_interp.h"
#include "xml_io.h"

extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;
extern const Numeric EARTH_GRAV_CONST;

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
  m(0,joker) = a_pos;
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
  m(0,joker) = a_pos;
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
        const Tensor3&        t_field,
        const Matrix&         r_geoid,
        const Matrix&         z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         a_pos,
        const Vector&         a_los )
{
  ppath_calc( ppath, ppath_step, ppath_step_agenda, atmosphere_dim, 
             p_grid, lat_grid, lon_grid, z_field, t_field, r_geoid, z_ground, 
                               cloudbox_on, cloudbox_limits, a_pos, a_los, 0 );
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
        const Matrix&    z_ground,
        // Control Parameters:
        const Numeric&   lmax )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here.

  if( atmosphere_dim == 1 )
    { ppath_step_geom_1d( ppath_step, p_grid, z_field(joker,0,0), 
                                         r_geoid(0,0), z_ground(0,0), lmax ); }

  else if( atmosphere_dim == 2 )
    { ppath_step_geom_2d( ppath_step, p_grid, lat_grid,
         z_field(joker,joker,0), r_geoid(joker,0), z_ground(joker,0), lmax ); }


  else if( atmosphere_dim == 3 )
    { ppath_step_geom_3d( ppath_step, p_grid, lat_grid, lon_grid,
                                          z_field, r_geoid, z_ground, lmax ); }

  else
    { throw runtime_error( "The atmospheric dimensionality must be 1-3." ); }
}



//! ppath_stepRefractionEuler
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-05-28
*/
void ppath_stepRefractionEuler(
        // WS Output:
              Ppath&      ppath_step,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
              Numeric&    refr_index,
        // WS Input:
        const Agenda&     refr_index_agenda,
        const Index&      atmosphere_dim,
        const Vector&     p_grid,
        const Vector&     lat_grid,
        const Vector&     lon_grid,
        const Tensor3&    z_field,
        const Tensor3&    t_field,
        const Tensor4&    vmr_field,
        const Matrix&     r_geoid,
        const Matrix&     z_ground,
        // Control Parameters:
        const Numeric&    lraytrace,
        const Numeric&    lmax )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we do only asserts. The keywords are checked here,
  // other input in the sub-functions to make them as simple as possible.

  assert( lraytrace > 0 );

  if( atmosphere_dim == 1 )
    { ppath_step_refr_1d( ppath_step, a_pressure, a_temperature, a_vmr_list,
                          refr_index, refr_index_agenda,
                          p_grid, z_field(joker,0,0), t_field(joker,0,0), 
                       vmr_field(joker,joker,0,0), r_geoid(0,0), z_ground(0,0),
                                           "linear_euler", lraytrace, lmax ); }

  else if( atmosphere_dim == 2 )
    { ppath_step_refr_2d( ppath_step, a_pressure, a_temperature, a_vmr_list,
                          refr_index, refr_index_agenda,
                          p_grid, lat_grid, z_field(joker,joker,0),
                       t_field(joker,joker,0), vmr_field(joker, joker,joker,0),
                          r_geoid(joker,0), z_ground(joker,0), 
                                           "linear_euler", lraytrace, lmax ); }

  else if( atmosphere_dim == 3 )
    { ppath_step_refr_3d( ppath_step, a_pressure, a_temperature, a_vmr_list,
                          refr_index, refr_index_agenda,
                          p_grid, lat_grid, lon_grid, z_field, 
                          t_field, vmr_field, r_geoid, z_ground, 
                                           "linear_euler", lraytrace, lmax ); }

  else
    { throw runtime_error( "The atmospheric dimensionality must be 1-3." ); }
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
  // tricks are needed as this WSM sets *r_geoid* for all crossings of the
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
      sensor_pos(joker,0) += r(0,0);
    }

  else
    {
      for( Index i=0; i<npos; i++ )
        {
          Vector lats(2);
          Index  pos_used;

          if( sensor_pos(i,1) < 90 )
            {
              lats[0]  = sensor_pos(i,1);
              lats[1]  = 90;
              pos_used = 0;
            }
          else
            {
              lats[0]  = -90;
              lats[1]  = sensor_pos(i,1);
              pos_used = 1;
            }

          // The size of the r-matrix is set inside the function.
          Matrix r;
          r_geoidWGS84( r, 2, lats, Vector(0), -999, -999 );
          
          sensor_pos(i,0) += r(i,pos_used);
        }
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
    { sensor_pos(joker,0) += r_geoid(0,0); }

  else
    {
      // Check that positions in sensor_pos are inside the lat and lon grids
      if( min(sensor_pos(joker,1)) < lat_grid[0]  || 
                                    max(sensor_pos(joker,1)) > last(lat_grid) )
        throw runtime_error(
             "You have given a position with a latitude outside *lat_grid*." );
      if( atmosphere_dim == 3 )
        {
          if( min(sensor_pos(joker,2)) < lon_grid[0]  || 
                                    max(sensor_pos(joker,2)) > last(lon_grid) )
            throw runtime_error(
            "You have given a position with a longitude outside *lon_grid*." );
        }

      if( atmosphere_dim == 2 )
        {
          ArrayOfGridPos gp(npos);
          Matrix itw(npos,2);
          gridpos( gp, lat_grid, sensor_pos(joker,1) );
          interpweights( itw, gp );
          Vector v_rgeoid(npos);
          interp( v_rgeoid, itw, r_geoid(joker,0), gp );
          for( Index i=0; i<npos; i++ )
            { sensor_pos(i,0) += v_rgeoid[i]; } 
        }
      else
        {
          ArrayOfGridPos gp_lat(npos), gp_lon(npos);
          Matrix itw(npos,4);
          gridpos( gp_lat, lat_grid, sensor_pos(joker,1) );
          gridpos( gp_lon, lon_grid, sensor_pos(joker,2) );
          interpweights( itw, gp_lat, gp_lon );
          Vector v_rgeoid(npos);
          interp( v_rgeoid, itw, r_geoid, gp_lat, gp_lon );
          for( Index i=0; i<npos; i++ )
            { sensor_pos(i,0) += v_rgeoid[i]; } 
        }
    }
}

//! VectorZtanToZa
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson and Mattias Ekström
   \date   2003-01-27
*/
void VectorZtanToZa(
        // WS Output:
              Vector&    za_vector,
	const String&    za_v_name,
        // WS Input:
        const Matrix&    sensor_pos,
        const Vector&    ztan_vector,
        const String&    ztan_v_name,
        const Numeric&   r_geoid )
{
  const Index   npos = sensor_pos.nrows();

  if( ztan_vector.nelem() != npos )
    {
      ostringstream os;
      os << "The number of altitudes in *" << ztan_v_name << "* must\n"
	 << "match the number of positions in *sensor_pos*.";
      throw runtime_error( os.str() );
    }

  out2 << "   Filling *" << za_v_name << "* with zenith angles, based on\n"
       << "   tangent altitudes from *" << ztan_v_name << ".\n";

  za_vector.resize(npos);

  for( Index i=0; i<npos; i++ )
    { za_vector[i] = geompath_za_at_r( r_geoid + ztan_vector[i], 100,
                                                           sensor_pos(i,0) ); }

}

//! ZaSatOccultation
/*!
   See the the online help (arts -d FUNCTION_NAME)

   \author Mattias Ekström
   \date   2003-02-10
*/
void ZaSatOccultation(
        // WS Output:
        Ppath&                  ppath_step,
        // WS Generic Output:
        Vector&                 za_out,
        // WS Generic Output Names:
        const String&           za_out_name,
        // WS Input:
        const Agenda&           ppath_step_agenda,
        const Index&            atmosphere_dim,
        const Vector&           p_grid,
        const Vector&           lat_grid,
        const Vector&           lon_grid,
        const Tensor3&          z_field,
        const Tensor3&          t_field,
        const Matrix&           r_geoid,
        const Matrix&           z_ground,
        // WS Control Parameters:
        const Numeric&          z_recieve,
        const Numeric&          z_send,
        const Numeric&          t_sample,
        const Numeric&          z_scan_low,
        const Numeric&          z_scan_high )
{
  //Check that indata is valid
  assert( atmosphere_dim ==1 );
  if ( z_scan_low >= z_scan_high ) {
    ostringstream os;
    os << "The lowest scan altitude *z_scan_low*number is higher or equal\n"
    << "to the highest scan altitude *z_scan_high*.";
    throw runtime_error( os.str() );
  }
  chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
  chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid,
                                                                    lon_grid );
  chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid,
                                                                    lon_grid );
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_ground", z_ground, atmosphere_dim, lat_grid, lon_grid );

  //Convert heights to radius
  const Numeric     r_recieve = r_geoid(0,0) + z_recieve;
  const Numeric     r_send = r_geoid(0,0) + z_send;

  //Extend upper and lower scan limits to include the limits
  //and convert to radius
  const Numeric     d = 3e3;
  const Numeric     r_scan_max = z_scan_high + d + r_geoid(0,0);
  const Numeric     r_scan_min = z_scan_low - d + r_geoid(0,0);

  //Calculate max and min zenith angles to cover z_scan_min to z_scan_max
  const Numeric     za_ref_max = geompath_za_at_r( r_scan_min, 100, r_recieve );
  const Numeric     za_ref_min = geompath_za_at_r( r_scan_max, 100, r_recieve );

  //Create vector with equally spaced zenith angles
  const Numeric     za_step = 600.0;
  Vector za_ref;
  linspace( za_ref, za_ref_min, za_ref_max, 1/za_step );

  //Calculate ppaths for all angles in za_ref
  Vector lat_ref( za_ref.nelem() );
  Vector z_ref( za_ref.nelem() );
  Index i;
  for ( i=0; i<za_ref.nelem(); i++ ) {
    Ppath ppath;
    ppath_calc( ppath, ppath_step, ppath_step_agenda, atmosphere_dim,
                p_grid, lat_grid, lon_grid, z_field, t_field, r_geoid,
                z_ground, 0, ArrayOfIndex(0), Vector(1,r_recieve),
                Vector(1,za_ref[i]), 0 );
    if ( ppath.tan_pos.nelem() != 0 ) {
      //Calculate the latitude from the recieving sensor, concisting of
      //the latitude of the atmospheric exist point and the geometric
      //path from that point to the transmitting satellite
      lat_ref[i] = geompath_lat_at_za( ppath.los(ppath.np-1,0),
        ppath.pos(ppath.np-1,1), geompath_za_at_r(
            ppath.pos(ppath.np-1,0)*sin(DEG2RAD*ppath.los(ppath.np-1,0)),
            10, r_send ));
      z_ref[i] = ppath.tan_pos[0] - r_geoid(0,0);
    } else {
      //If ppath_calc hits the ground, use linear extrapolation to find the
      //last point
      lat_ref[i] = 2*lat_ref[i-1] - lat_ref[i-2];
      z_ref[i] = 2*z_ref[i-1] - z_ref[i-2];
      break;
    }
  }

  //Create vectors for ppaths that passed above ground
  Vector za_ref_above = za_ref[Range(0,i)];
  Vector z_ref_above = z_ref[Range(0,i)];
  Vector lat_ref_above = lat_ref[Range(0,i)];

  //Retrieve the interpolated latitudes that corresponds to z_scan_high
  //and z_scan_low
  GridPos gp_high, gp_low;
  gridpos( gp_high, z_ref_above, z_scan_high );
  gridpos( gp_low, z_ref_above, z_scan_low );

  Vector itw_high(2), itw_low(2);
  interpweights( itw_high, gp_high );
  interpweights( itw_low, gp_low );

  Numeric lat_min = interp( itw_high, lat_ref_above, gp_high );
  Numeric lat_max = interp( itw_low, lat_ref_above, gp_low );

  //Create latitude vector with steps corresponding to t_sample
  //First calculate the lat_sample corresponding to t_sample
  //NB: satellite velocities defined as opposite
  Numeric w_send, w_recieve, lat_dot, lat_sample;
  w_send = sqrt(EARTH_GRAV_CONST/r_send) / r_send;
  w_recieve = sqrt(EARTH_GRAV_CONST/r_recieve) / r_recieve;
  lat_dot = (w_send + w_recieve) * RAD2DEG;
  lat_sample = t_sample * lat_dot;

  Vector lat_out;
  linspace( lat_out, lat_min, lat_max, lat_sample );

  //Interpolate values in za_out corresponding to latitudes in lat_out
  ArrayOfGridPos gp(lat_out.nelem());
  gridpos( gp, lat_ref_above, lat_out );
  Matrix itw(lat_out.nelem(), 2);
  interpweights( itw, gp );
  za_out.resize( lat_out.nelem() );
  interp( za_out, itw, za_ref_above, gp );

}


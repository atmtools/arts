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



//! DoGridcell2D
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-12-21
*/
void DoGridcell2D(
        const Numeric&   r_start,
        const Numeric&   lat_start,
        const Numeric&   za_start,
        const Numeric&   lmax,
        const Numeric&   lat1,
        const Numeric&   lat3,
        const Numeric&   r1a,
        const Numeric&   r3a,
        const Numeric&   r3b,
        const Numeric&   r1b,
        const Numeric&   rground1,
        const Numeric&   rground3 )
{
  const Numeric   ppc = geometrical_ppc( r_start, za_start );

  Vector    r, lat, za;
  Numeric   lstep;
  Index     endface, tanpoint;

  do_gridcell_2d( r, lat, za, lstep, endface, tanpoint, r_start, lat_start, 
                  za_start, ppc, lmax, lat1, lat3, r1a, r3a, r3b, r1b, 
                                                          rground1, rground3 );

  String filename;

  filename = "";
  filename_xml( filename, "r" );
  xml_write_to_file ( filename, r );
  filename = "";
  filename_xml( filename, "lat" );
  xml_write_to_file ( filename, lat );
  filename = "";
  filename_xml( filename, "za" );
  xml_write_to_file ( filename, za );
  filename = "";
  filename_xml( filename, "lstep" );
  xml_write_to_file ( filename, lstep );
  filename = "";
  filename_xml( filename, "endface" );
  xml_write_to_file ( filename, endface );
}



//! DoGridcell3D
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-12-25
*/
void DoGridcell3D(
        const Numeric&   r_start,
        const Numeric&   lat_start,
        const Numeric&   lon_start,
        const Numeric&   za_start,
        const Numeric&   aa_start,
        const Numeric&   lmax,
        const Numeric&   lat1,
        const Numeric&   lat3,
        const Numeric&   lon5,
        const Numeric&   lon6,
        const Numeric&   r15a,
        const Numeric&   r35a,
        const Numeric&   r36a,
        const Numeric&   r16a,
        const Numeric&   r15b,
        const Numeric&   r35b,
        const Numeric&   r36b,
        const Numeric&   r16b,
        const Numeric&   rground15,
        const Numeric&   rground35,
        const Numeric&   rground36,
        const Numeric&   rground16 )
{
  const Numeric   ppc = geometrical_ppc( r_start, za_start );

  Vector    r, lat, lon, za, aa;
  Numeric   lstep;
  Index     endface, tanpoint;

  do_gridcell_3d( r, lat, lon, za, aa, lstep, endface, tanpoint, 
                  r_start, lat_start, lon_start, za_start, aa_start, ppc, lmax,
                  lat1, lat3,  lon5, lon6, 
                  r15a, r35a, r36a, r16a, r15b, r35b, r36b, r16b,
                                  rground15, rground35, rground36, rground16 );

  String filename;

  filename = "";
  filename_xml( filename, "r" );
  xml_write_to_file ( filename, r );
  filename = "";
  filename_xml( filename, "lat" );
  xml_write_to_file ( filename, lat );
  filename = "";
  filename_xml( filename, "lon" );
  xml_write_to_file ( filename, lon );
  filename = "";
  filename_xml( filename, "za" );
  xml_write_to_file ( filename, za );
  filename = "";
  filename_xml( filename, "aa" );
  xml_write_to_file ( filename, aa );
  filename = "";
  filename_xml( filename, "lstep" );
  xml_write_to_file ( filename, lstep );
  filename = "";
  filename_xml( filename, "endface" );
  xml_write_to_file ( filename, endface );
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



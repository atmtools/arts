/* Copyright (C) 2002-2008 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                            
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
#include "refraction.h"
#include "m_general.h"

extern const Numeric RAD2DEG;
extern const Numeric DEG2RAD;
extern const Numeric EARTH_GRAV_CONST;

/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losSet(
        // WS Output:
              Vector&    rte_los,
        // WS Input:
        const Index&     atmosphere_dim,
        // Control Parameters:
        const Numeric&   za,
        const Numeric&   aa )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  if( atmosphere_dim == 1 )
    { rte_los.resize(1); }
  else
    {
      rte_los.resize(2);
      rte_los[1] = aa;
    }
  rte_los[0] = za;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_posAddGeoidWGS84(
        // WS Output:
              Vector&    rte_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        const Numeric&   latitude_1d,
        const Numeric&   meridian_angle_1d )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_vector_length( "rte_pos", rte_pos, atmosphere_dim );

  // Use *sensor_posAddGeoidWGS84* to perform the calculations.
  Matrix m(1,rte_pos.nelem());
  m(0,joker) = rte_pos;
  sensor_posAddGeoidWGS84( m, atmosphere_dim, latitude_1d, meridian_angle_1d);
  rte_pos[0] = m(0,0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_posAddRgeoid(
        // WS Output:
              Vector&    rte_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        const Vector&    lat_grid,
        const Vector&    lon_grid,
        const Matrix&    r_geoid )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  chk_vector_length( "rte_pos", rte_pos, atmosphere_dim );

  // Use *sensor_posAddRgeoid* to perform the calculations.
  Matrix m(1,rte_pos.nelem());
  m(0,joker) = rte_pos;
  sensor_posAddRgeoid( m, atmosphere_dim, lat_grid, lon_grid, r_geoid );
  rte_pos[0] = m(0,0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_posSet(
        // WS Output:
              Vector&    rte_pos,
        // WS Input:
        const Index&     atmosphere_dim,
        // Control Parameters:
        const Numeric&   r_or_z,
        const Numeric&   lat,
        const Numeric&   lon )
{
  // Check input
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  rte_pos.resize(atmosphere_dim);
  rte_pos[0] = r_or_z;
  if( atmosphere_dim >= 2 )
    { rte_pos[1] = lat; }
  if( atmosphere_dim == 3 )
    { rte_pos[2] = lon; }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_pos_and_losFromTangentPressure(
                    // WS Output:
                    Vector&          rte_pos,
                    Vector&          rte_los,
                    Ppath&           ppath,
                    Ppath&           ppath_step,
                    // WS Input:
                    const Index&     atmosphere_dim,
                    const Vector&    p_grid,
                    const Tensor3&   z_field,
                    const Vector &     lat_grid,      
                    const Vector &     lon_grid,
                    const Agenda &     ppath_step_agenda,
                    const Matrix &     r_geoid,
                    const Matrix &     z_surface,
                    // Control Parameters:
                    const Numeric&   tan_p
                    )
{
  //This function is only ready for 1D calculations
  if (atmosphere_dim > 1)
    throw runtime_error(
      "Sorry, this function currently only works for atmosphere_dim==1");
  

  Index np=p_grid.nelem();
  Vector log_p_grid(np);
  Numeric log_p=log10(tan_p);
  GridPos gp;
  Vector itw(2);
  Numeric z;
  rte_pos.resize(atmosphere_dim);
  rte_los.resize(atmosphere_dim);
  
  //find z for given tangent pressure
  for (Index i=0;i<np;i++){log_p_grid[i]=log10(p_grid[i]);}
  
  gridpos(gp,log_p_grid,log_p);
  out1 << "gp.index = " << gp.idx << "\n";
  interpweights(itw,gp);
  z = interp(itw,z_field(Range(joker),0,0),gp);
  out1 <<  " Tangent pressure corresponds to an altitude of " << z << "m.\n";
  
  //Use ppath to find the desired point on the edge of the atmosphere
  rte_pos[0]=z+r_geoid(0,0);
  rte_los[0]=90;

  Index cloudbox_on=0;
  ArrayOfIndex cloudbox_limits(2);
  bool outside_cloudbox=true;
  ppath_calc(ppath,ppath_step,ppath_step_agenda,atmosphere_dim,p_grid,lat_grid,
         lon_grid,z_field, r_geoid, z_surface,cloudbox_on, cloudbox_limits,
         rte_pos, rte_los, outside_cloudbox);
  rte_pos = ppath.pos(ppath.np-1,Range(0,atmosphere_dim));
  rte_los = ppath.los(ppath.np-1,joker);
  rte_los[0]=180.0-rte_los[0];
}


/* Workspace method: Doxygen documentation will be auto-generated */
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
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         rte_pos,
        const Vector&         rte_los )
{
  ppath_calc( ppath, ppath_step, ppath_step_agenda, atmosphere_dim, 
              p_grid, lat_grid, lon_grid, z_field, r_geoid, z_surface, 
              cloudbox_on, cloudbox_limits, rte_pos, rte_los, 1 );
}


/* Workspace method: Doxygen documentation will be auto-generated */
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
        const Matrix&    z_surface,
        // Control Parameters:
        const Numeric&   lmax )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we perform asserts in the sub-functions, but no checks 
  // here.

  if( atmosphere_dim == 1 )
    { ppath_step_geom_1d( ppath_step, p_grid, z_field(joker,0,0), 
                                         r_geoid(0,0), z_surface(0,0), lmax ); }

  else if( atmosphere_dim == 2 )
    { ppath_step_geom_2d( ppath_step, p_grid, lat_grid,
         z_field(joker,joker,0), r_geoid(joker,0), z_surface(joker,0), lmax ); }


  else if( atmosphere_dim == 3 )
    { ppath_step_geom_3d( ppath_step, p_grid, lat_grid, lon_grid,
                                          z_field, r_geoid, z_surface, lmax ); }

  else
    { throw runtime_error( "The atmospheric dimensionality must be 1-3." ); }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ppath_stepRefractionEuler(
        // WS Output:
              Ppath&      ppath_step,
              Numeric&    rte_pressure,
              Numeric&    rte_temperature,
              Vector&     rte_vmr_list,
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
        const Matrix&     z_surface,
        // Control Parameters:
        const Numeric&    lraytrace,
        const Numeric&    lmax )
{
  // Input checks here would be rather costly as this function is called
  // many times. So we do only asserts. The keywords are checked here,
  // other input in the sub-functions to make them as simple as possible.

  assert( lraytrace > 0 );

  if( atmosphere_dim == 1 )
    { ppath_step_refr_1d( ppath_step, rte_pressure, rte_temperature, 
                          rte_vmr_list, refr_index, refr_index_agenda,
                          p_grid, z_field(joker,0,0), t_field(joker,0,0), 
                       vmr_field(joker,joker,0,0), r_geoid(0,0), z_surface(0,0),
                                           "linear_euler", lraytrace, lmax ); }

  else if( atmosphere_dim == 2 )
    { ppath_step_refr_2d( ppath_step, rte_pressure, rte_temperature, 
                          rte_vmr_list, refr_index, refr_index_agenda,
                          p_grid, lat_grid, z_field(joker,joker,0),
                       t_field(joker,joker,0), vmr_field(joker, joker,joker,0),
                          r_geoid(joker,0), z_surface(joker,0), 
                                           "linear_euler", lraytrace, lmax ); }

  else if( atmosphere_dim == 3 )
    { ppath_step_refr_3d( ppath_step, rte_pressure, rte_temperature, 
                          rte_vmr_list, refr_index, refr_index_agenda,
                          p_grid, lat_grid, lon_grid, z_field, 
                          t_field, vmr_field, r_geoid, z_surface, 
                                           "linear_euler", lraytrace, lmax ); }

  else
    { throw runtime_error( "The atmospheric dimensionality must be 1-3." ); }
}


/* Workspace method: Doxygen documentation will be auto-generated */
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


/* Workspace method: Doxygen documentation will be auto-generated */
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


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZa1D(
                        // WS Generic Output:
                        Vector&             za_vector,
                        // WS Input:
                        const Matrix&       sensor_pos,
                        const Matrix&       r_geoid,
                        const Index&        atmosphere_dim,
                        // WS Generic Input:
                        const Vector&       ztan_vector)
{
  if( atmosphere_dim != 1 ) {
    throw runtime_error( "The function can only be used for 1D atmospheres." );
  }

  const Index   npos = sensor_pos.nrows();

  if( ztan_vector.nelem() != npos )
    {
      ostringstream os;
      os << "The number of altitudes in the geometric tangent altitude vector must\n"
         << "match the number of positions in *sensor_pos*.";
      throw runtime_error( os.str() );
    }

  za_vector.resize(npos);

  for( Index i=0; i<npos; i++ )
    { za_vector[i] = geompath_za_at_r( r_geoid(0,0) + ztan_vector[i], 100,
                                                           sensor_pos(i,0) ); }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void VectorZtanToZaRefr1D(
                        // WS Output:
                        Numeric&            refr_index,
                        Numeric&            rte_pressure,
                        Numeric&            rte_temperature,
                        Vector&             rte_vmr_list,
                        // WS Generic Output:
                        Vector&             za_vector,
                        // WS Input:
                        const Agenda&       refr_index_agenda,
                        const Matrix&       sensor_pos,
                        const Vector&       p_grid,
                        const Tensor3&      t_field,
                        const Tensor3&      z_field,
                        const Tensor4&      vmr_field,
                        const Matrix&       r_geoid,
                        const Index&        atmosphere_dim,
                        // WS Generic Input:
                        const Vector&       ztan_vector)
{
  if( atmosphere_dim != 1 ) {
    throw runtime_error( "The function can only be used for 1D atmospheres." );
  }

  if( ztan_vector.nelem() != sensor_pos.nrows() ) {
    ostringstream os;
    os << "The number of altitudes in true tangent altitude vector must\n"
       << "match the number of positions in *sensor_pos*.";
    throw runtime_error( os.str() );
  }

  //No output from get_refr_index_1d
  Index agenda_verb = 1;

  //Set za_vector's size equal to ztan_vector
  za_vector.resize( ztan_vector.nelem() );

  //Calculate refractive index for the tangential altitudes
  for( Index i=0; i<ztan_vector.nelem(); i++ ) {
    get_refr_index_1d( refr_index, rte_pressure, rte_temperature, rte_vmr_list, 
      refr_index_agenda, agenda_verb, p_grid, r_geoid(0,0), 
      z_field(joker,0,0), t_field(joker,0,0), vmr_field(joker,joker,0,0),
      ztan_vector[i]+r_geoid(0,0) );

    //Calculate zenith angle
    za_vector[i] = 180-asin(refr_index*(ztan_vector[i]+r_geoid(0,0)) /
                            sensor_pos(i,0))*RAD2DEG;
  }

}


/* Workspace method: Doxygen documentation will be auto-generated */
void ZaSatOccultation(
        // WS Output:
        Ppath&                  ppath_step,
        // WS Generic Output:
        Vector&                 za_out,
        // WS Input:
        const Agenda&           ppath_step_agenda,
        const Index&            atmosphere_dim,
        const Vector&           p_grid,
        const Vector&           lat_grid,
        const Vector&           lon_grid,
        const Tensor3&          z_field,
        const Matrix&           r_geoid,
        const Matrix&           z_surface,
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
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_surface", z_surface, atmosphere_dim, lat_grid, lon_grid );

  //Convert heights to radius
  const Numeric     r_recieve = r_geoid(0,0) + z_recieve;
  const Numeric     r_send = r_geoid(0,0) + z_send;

  //Extend upper and lower scan limits to include the limits
  //and convert to radius
  const Numeric    d = 3e3;
  const Numeric    r_scan_max = z_scan_high + d + r_geoid(0,0);
  const Numeric    r_scan_min = z_scan_low - d + r_geoid(0,0);

  //Calculate max and min zenith angles to cover z_scan_min to z_scan_max
  const Numeric    za_ref_max = geompath_za_at_r( r_scan_min, 100, r_recieve );
  const Numeric    za_ref_min = geompath_za_at_r( r_scan_max, 100, r_recieve );

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
                p_grid, lat_grid, lon_grid, z_field, r_geoid,
                z_surface, 0, ArrayOfIndex(0), Vector(1,r_recieve),
                Vector(1,za_ref[i]), 1 );
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
      //If ppath_calc hits the surface, use linear extrapolation to find the
      //last point
      lat_ref[i] = 2*lat_ref[i-1] - lat_ref[i-2];
      z_ref[i] = 2*z_ref[i-1] - z_ref[i-2];
      break;
    }
  }

  //Create vectors for ppaths that passed above surface
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

/* Copyright (C) 2002 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
                      Stefan Buehler   <sbuehler@uni-bremen.de>
                            
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
  \file   rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-29

  \brief  Functions to solve radiative transfer tasks.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "auto_md.h"
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"
#include "lin_alg.h"

extern const Numeric DEG2RAD;



/*===========================================================================
  === The functions in alphabetical order
  ===========================================================================*/

//! get_radiative_background
/*!
    Sets *i_rte* to the radiative background for a propagation path.

    The function uses *ppath* to determine the radiative background for a 
    propagation path. The possible backgrounds are listed in the header of
    the function ppath_set_background (in ppath.cc).

    The main purpose of the function is set *i_rte*. It is NOT needed to set 
    *i_rte* to the correct size before calling the function. The size is set
    to be [f_grid.nelem(),stokes_dim].

    \param   i_rte              Output: As the WSV with the same name.
    \param   ppath_step         Output: As the WSV with the same name.
    \param   rte_pos            Output: As the WSV with the same name.
    \param   rte_los            Output: As the WSV with the same name.
    \param   rte_gp_p           Output: As the WSV with the same name.
    \param   rte_gp_lat         Output: As the WSV with the same name.
    \param   rte_gp_lon         Output: As the WSV with the same name.
    \param   i_space            Output: As the WSV with the same name.
    \param   surface_emission    Output: As the WSV with the same name.
    \param   surface_los         Output: As the WSV with the same name.
    \param   surface_refl_coeffs Output: As the WSV with the same name.
    \param   ppath              Input: As the WSV with the same name.
    \param   ppath_step_agenda  Input: As the WSV with the same name.
    \param   rte_agenda         Input: As the WSV with the same name.
    \param   i_space_agenda     Input: As the WSV with the same name.
    \param   surface_agenda Input: As the WSV with the same name.
    \param   atmosphere_dim     Input: As the WSV with the same name.
    \param   p_grid             Input: As the WSV with the same name.
    \param   lat_grid           Input: As the WSV with the same name.
    \param   lon_grid           Input: As the WSV with the same name.
    \param   z_field            Input: As the WSV with the same name.
    \param   t_field            Input: As the WSV with the same name.
    \param   r_geoid            Input: As the WSV with the same name.
    \param   z_surface           Input: As the WSV with the same name.
    \param   cloudbox_on        Input: As the WSV with the same name.
    \param   cloudbox_limits    Input: As the WSV with the same name.
    \param   scat_i_p           Input: As the WSV with the same name.
    \param   scat_i_lat         Input: As the WSV with the same name.
    \param   scat_i_lon         Input: As the WSV with the same name.
    \param   scat_za_grid       Input: As the WSV with the same name.
    \param   scat_aa_grid       Input: As the WSV with the same name.
    \param   f_grid             Input: As the WSV with the same name.
    \param   stokes_dim         Input: As the WSV with the same name.
    \param   antenna_dim        Input: As the WSV with the same name.
    \param   agenda_verb        Input: Argument handed to agendas to control
                                verbosity.

    \author Patrick Eriksson 
    \date   2002-09-17
*/
void get_radiative_background(
              Matrix&         i_rte,
              Ppath&          ppath_step,
              Vector&         rte_pos,
              Vector&         rte_los,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Matrix&         i_space,
              Matrix&         surface_emission,
              Matrix&         surface_los, 
              Tensor4&        surface_refl_coeffs,
              Ppath&          ppath,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         i_space_agenda,
        const Agenda&         surface_agenda,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Tensor7&        scat_i_p,
        const Tensor7&        scat_i_lat,
        const Tensor7&        scat_i_lon,
        ConstVectorView       scat_za_grid,
        ConstVectorView       scat_aa_grid,
        ConstVectorView       f_grid,
        const Index&          stokes_dim,
        const Index&          agenda_verb,
        const Index&          scat_za_interp)
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Resize i_rte to have the correct the size
  i_rte.resize( nf, stokes_dim );

  // Set rte_pos, rte_gp_XXX and rte_los to match the last point in ppath.
  // The agendas below use different combinations of these varaibles.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  rte_pos.resize( atmosphere_dim );
  rte_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  rte_los.resize( ppath.los.ncols() );
  rte_los = ppath.los(np-1,joker);
  gridpos_copy( rte_gp_p, ppath.gp_p[np-1] );
  if( atmosphere_dim > 1 )
    { gridpos_copy( rte_gp_lat, ppath.gp_lat[np-1] ); }
  if( atmosphere_dim > 2 )
    { gridpos_copy( rte_gp_lon, ppath.gp_lon[np-1] ); }

  out3 << "Radiative background: " << ppath.background << "\n";


    // Initialize i_rte to the radiative background
  switch ( ppath_what_background( ppath ) )
    {

    case 1:   //--- Space ---------------------------------------------------- 
      {
        chk_not_empty( "i_space_agenda", i_space_agenda );
        i_space_agenda.execute( agenda_verb );

        if( i_space.nrows() != nf  ||  i_space.ncols() != stokes_dim )
          throw runtime_error( "The size of *i_space* returned from "
                                           "*i_space_agenda* is not correct.");

        i_rte = i_space;
      }
      break;


    case 2:   //--- The surface -----------------------------------------------
      {
        chk_not_empty( "surface_agenda", surface_agenda );
        surface_agenda.execute( agenda_verb );

        // Check returned variables
        if( surface_emission.nrows() != nf  ||  
                                        surface_emission.ncols() != stokes_dim )
          throw runtime_error(
                  "The size of the created *surface_emission* is not correct.");

        // If *surface_los* is empty, the upwelling radiation is just
        // the surface emission. Else, the downwelling radiation shall
        // be calculated and be added with the surface emission.
        //
        if( surface_los.nrows() == 0 )
          {
            i_rte = surface_emission;
          }

        else
          {
            // Make a copy of *ppath* as the WSV will be changed below, but the
            // original data is needed when going back to *rte_calc*.
            Ppath   pp_copy;
            ppath_init_structure( pp_copy, atmosphere_dim, ppath.np );
            ppath_copy( pp_copy, ppath );

            // Copy surface data to local variables the surface variables can be
            // changed below by call of *rte_calc*. There is no need to copy
            // *surface_los*.
            Index     nlos = surface_los.nrows();
            Matrix    surface_emission_local(nf,stokes_dim);
            Tensor4   surface_refl_coeffs_local(nlos,nf,stokes_dim,stokes_dim);
            //
            surface_emission_local    = surface_emission;
            surface_refl_coeffs_local = surface_refl_coeffs;

            // Sum of reflected radiation
            Matrix   i_sum(nf,stokes_dim,0);

            // Each surface los is here treated as a 
            // measurement block, with no za and aa grids.
            Matrix   sensor_pos( 1, rte_pos.nelem() );
                     sensor_pos( 0, joker ) = rte_pos;
            Matrix   sensor_los( 1, rte_los.nelem() );

            // Use some local variables to avoid unwanted side effects
            Vector   y_local;
            
            // Use a dummy variables for the sensor
            Vector   sensor_rot;
            Sparse   sensor_response;
            Index    antenna_dim = 1;
            Vector   mblock_za_grid(1,0);
            Vector   mblock_aa_grid(0);


            for( Index ilos=0; ilos<nlos; ilos++ )
              {
                sensor_los( 0, joker ) = surface_los( ilos, joker);

                rte_calc( y_local, ppath, ppath_step, i_rte,
                   rte_pos, rte_los, rte_gp_p, rte_gp_lat, rte_gp_lon,
                   i_space, surface_emission, surface_los, surface_refl_coeffs, 
                   ppath_step_agenda, rte_agenda, i_space_agenda, 
                   surface_agenda, atmosphere_dim, p_grid, lat_grid, 
                   lon_grid, z_field, t_field, r_geoid, z_surface, 
                   cloudbox_on, cloudbox_limits, scat_i_p, scat_i_lat,
                   scat_i_lon, scat_za_grid, scat_aa_grid, 
                   sensor_response, sensor_pos, sensor_los,
                          sensor_rot, f_grid, stokes_dim, antenna_dim,
                          mblock_za_grid, mblock_aa_grid, false, false, 
                          agenda_verb, 0 );

                // Add the calculated spectrum (*i_rte*) to *i_sum*,
                // considering the reflection coeff. matrix.
                //
                for( Index iv=0; iv<nf; iv++ )
                  {
                    if( stokes_dim == 1 )
                      {
                        i_sum(iv,0) += 
                                 surface_refl_coeffs(ilos,iv,0,0) * i_rte(iv,0);
                      }
                    else
                      {
                        for( Index is0=0; is0<stokes_dim; is0++ )
                          { 
                            for( Index is=0; is<stokes_dim; is++ )
                              {
                                i_sum(iv,is0) += 
                                      surface_refl_coeffs(ilos,iv,is0,is) * 
                                                                  i_rte(iv,is);
                              }
                          }
                      }
                  }
              }

            // Copy from *i_sum* to *i_rte*, and add the surface emission
            for( Index iv=0; iv<nf; iv++ )
              {
                for( Index is=0; is<stokes_dim; is++ )
                  { 
                    i_rte(iv,is) = i_sum(iv,is) + surface_emission_local(iv,is);
                  }
              }

            // Copy data back to *ppath*.
            ppath_init_structure( ppath, atmosphere_dim, pp_copy.np );
            ppath_copy( ppath, pp_copy );
          }
      }
      break;


    case 3:   //--- Cloudbox surface -----------------------------------------

      {
        if(scat_za_interp == 0)
          {
            CloudboxGetOutgoing( i_rte, "i_rte", scat_i_p, scat_i_lat, 
                                 scat_i_lon, 
                                 rte_gp_p, rte_gp_lat, rte_gp_lon,
                                 rte_los, 
                                 cloudbox_on, cloudbox_limits,
                                 atmosphere_dim,
                                 stokes_dim, scat_za_grid, scat_aa_grid,
                                 f_grid );
          }
        else
          {
            CloudboxGetOutgoingCubic( i_rte, "i_rte", scat_i_p, scat_i_lat, 
                                      scat_i_lon, 
                                      rte_gp_p, rte_gp_lat, rte_gp_lon,
                                      rte_los, 
                                      cloudbox_on, cloudbox_limits,
                                      atmosphere_dim,
                                      stokes_dim, scat_za_grid, scat_aa_grid,
                                      f_grid );
          }
      }
      break;


    case 4:   //--- Inside of cloudbox ---------------------------------------
      {
        throw runtime_error( "RTE calculations starting inside the cloud box "
                                                       "are not yet handled" );
      }
      break;


    default:  //--- ????? ----------------------------------------------------
      // Are we here, the coding is wrong somewhere
      assert( false );
    }
}



//! surface_specular_los
/*!
    Calculates the LOS for a specular surface reflection.

    The function calculates the LOS for the downwelling radiation to consider
    when the surface is treated to give no scattering. The tilt of the surface
    (that is, a change in radius as a function of latitude etc.) is considered.

    \param   los                Output: As the WSV with the same name.
    \param   atmosphere_dim     Input: As the WSV with the same name.
    \param   r_geoid            Input: As the WSV with the same name.
    \param   z_surface           Input: As the WSV with the same name.
    \param   lat_grid           Input: As the WSV with the same name.
    \param   lon_grid           Input: As the WSV with the same name.
    \param   rte_gp_lat         Input: As the WSV with the same name.
    \param   rte_gp_lon         Input: As the WSV with the same name.
    \param   rte_los            Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-09-22
*/
void surface_specular_los(
              VectorView   los,
        const Index&       atmosphere_dim,
        ConstMatrixView    r_geoid,
        ConstMatrixView    z_surface,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        const GridPos&     rte_gp_lat,
        const GridPos&     rte_gp_lon,
        ConstVectorView    rte_los )
{
  assert( atmosphere_dim >= 1  &&  atmosphere_dim <= 3 );

  if( atmosphere_dim == 1 )
    {
      assert( r_geoid.nrows() == 1 );
      assert( r_geoid.ncols() == 1 );
      assert( z_surface.nrows() == 1 );
      assert( z_surface.ncols() == 1 );
      assert( rte_los.nelem() == 1 );
      assert( rte_los[0] > 90 );      // Otherwise surface refl. not possible
      assert( rte_los[0] <= 180 ); 
      assert( los.nelem() == 1 );

      los[0] = 180 - rte_los[0];
    }

  else if( atmosphere_dim == 2 )
    {
      assert( r_geoid.ncols() == 1 );
      assert( z_surface.ncols() == 1 );
      assert( rte_los.nelem() == 1 );
      assert( los.nelem() == 1 );
      assert( abs(rte_los[0]) <= 180 ); 

      
      // Calculate LOS neglecting any tilt of the surface
      los[0] = sign( rte_los[0] ) * 180 - rte_los[0];

      // Interpolate to get radius for the surface at reflection point.
      const Numeric r_surface =
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                              r_geoid, "r_geoid", rte_gp_lat, rte_gp_lon ) +
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                                z_surface, "z_surface", rte_gp_lat, rte_gp_lon );

      // Calculate surface slope (unit is m/deg).
      Numeric slope = plevel_slope_2d( lat_grid, r_geoid(joker,0), 
                                   z_surface(joker,0), rte_gp_lat, rte_los[0] );

      // Calculate surface (angular) tilt (unit is deg).
      Numeric tilt = plevel_angletilt( r_surface, slope );

      // Check that rte_los contains a downward LOS
      assert( is_los_downwards( rte_los[0], tilt ) );
      
      // Include surface tilt
      los[0] -= 2 * tilt;
    }

  else if( atmosphere_dim == 3 )
    {
      assert( rte_los.nelem() == 2 );
      assert( los.nelem() == 2 );
      assert( rte_los[0] >= 0 ); 
      assert( rte_los[0] <= 180 ); 
      assert( abs( rte_los[1] ) <= 180 ); 

      assert( rte_gp_lat.idx >= 0 );
      assert( rte_gp_lat.idx <= ( lat_grid.nelem() - 2 ) );
      assert( rte_gp_lon.idx >= 0 );
      assert( rte_gp_lon.idx <= ( lon_grid.nelem() - 2 ) );

      // Calculate LOS neglecting any tilt of the surface
      los[0] = 180 - rte_los[0];
      los[1] = rte_los[1];

      // Below you find a first version to include the surface tilt.
      // Is the solution correct? 
      // At least, it does not work for za = 180.
      if( 0 )
        {
      // Interpolate to get radius for the surface at reflection point.
      const Numeric r_surface =
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                              r_geoid, "r_geoid", rte_gp_lat, rte_gp_lon ) +
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                                z_surface, "z_surface", rte_gp_lat, rte_gp_lon );

      // Restore latitude and longitude values
      Vector   itw(2);
      Numeric  lat, lon;
      interpweights( itw, rte_gp_lat );
      lat = interp( itw, lat_grid, rte_gp_lat );
      interpweights( itw, rte_gp_lon );
      lon = interp( itw, lon_grid, rte_gp_lon );

      // Calculate surface slope along the viewing direction (unit is m/deg).
      //
      Index   ilat = gridpos2gridrange( rte_gp_lat, abs( rte_los[1] ) <= 90 ); 
      Index   ilon = gridpos2gridrange( rte_gp_lon, rte_los[1] >= 0 );
      //
      Numeric   lat1 = lat_grid[ilat];
      Numeric   lat3 = lat_grid[ilat+1];
      Numeric   lon5 = lon_grid[ilon];
      Numeric   lon6 = lon_grid[ilon+1];
      Numeric   r15  = r_geoid(ilat,ilon) + z_surface(ilat,ilon);
      Numeric   r35  = r_geoid(ilat+1,ilon) + z_surface(ilat+1,ilon);
      Numeric   r16  = r_geoid(ilat,ilon+1) + z_surface(ilat,ilon+1);
      Numeric   r36  = r_geoid(ilat+1,ilon+1) + z_surface(ilat+1,ilon+1);
      //
      Numeric slope = plevel_slope_3d( lat1, lat3, lon5, lon6, 
                                    r15, r35, r36, r16, lat, lon, rte_los[1] );

      // Calculate surface (angular) tilt (unit is deg).
      Numeric tilt = plevel_angletilt( r_surface, slope );

      // Check that rte_los contains a downward LOS
      assert( is_los_downwards( rte_los[0], tilt ) );
      
      // Include surface tilt in zenith angle
      los[0] -= 2 * tilt;

      // Calculate surface slope along the viewing direction (unit is m/deg).
      Numeric   aa = rte_los[1] + 90;
      if( aa > 180 )
        { aa -= 360; }
      //
      slope = plevel_slope_3d( lat1, lat3, lon5, lon6, 
                                            r15, r35, r36, r16, lat, lon, aa );

      // Calculate surface (angular) tilt (unit is deg).
      tilt = plevel_angletilt( r_surface, slope );

      // Include surface tilt in azimuth angle
      los[1] -= tilt;
        }
    }
}



//! rte_calc
/*!
    Implementation of RteCalc.

    This function performs the actual calculations for RteCalc. This
    function was introduced to make it possible to turn off input
    checks as they are rather costly, if repeated many times. This
    feature is controlled by *check_input*. This variable shall be set
    to br true when called from RteCalc. When called from
    CloudboxGetIncoming this variable can be set to false, except for
    the first of the many repeated calls of this function.

    The LOS angles are always checked. Those checks are performed in
    ppath_calc.

    If no sensor is applied (*apply_sensor* = false), data are not copied to 
    *y*. If the function is called to obtain monochromatic pencil beam values
    this has to be done for one line-of-sight at the time, and the result
    is returned by *i_rte*. This is done by setting *sensor_pos* and
    *sensor_los* to contain one measurement block matching the wanted position
    and line-of-sight. The antenna dimension (*antenna_dim*) can then be set 
    to 1, *mblock_za_grid* to have length 1 with 0 as only value and 
    *mblock_aa_grid* to be empty.

    All variables, beside the onses listed below, corresponde to the
    WSV with the same name.

    \param   check_input    Input: Boolean to perform check of function input.
                            True means that checka are performed.
    \param   apply_sensor   Input: Boolean to apply sensor responses. False
                            means no sensor.
    \param   agenda_verb    Input: Boolean to control the verbosity of called
                            agendas.

    \author Patrick Eriksson 
    \date   2003-01-18
*/
void rte_calc(
              Vector&         y,
              Ppath&          ppath,
              Ppath&          ppath_step,
              Matrix&         i_rte,
              Vector&         rte_pos,
              Vector&         rte_los,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Matrix&         i_space,
              Matrix&         surface_emission, 
              Matrix&         surface_los, 
              Tensor4&        surface_refl_coeffs,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         i_space_agenda,
        const Agenda&         surface_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Tensor7&        scat_i_p,
        const Tensor7&        scat_i_lat,
        const Tensor7&        scat_i_lon,
        const Vector&         scat_za_grid,
        const Vector&         scat_aa_grid,
        const Sparse&         sensor_response,
        const Matrix&         sensor_pos,
        const Matrix&         sensor_los,
        const Vector&         sensor_rot,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Index&          antenna_dim,
        const Vector&         mblock_za_grid,
        const Vector&         mblock_aa_grid,
        const bool&           check_input,
        const bool&           apply_sensor,
        const bool&           agenda_verb,
        const Index&          scat_za_interp)
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index nmblock = sensor_pos.nrows();
  const Index nza     = mblock_za_grid.nelem();

  // Number of azimuthal direction for pencil beam calculations
  Index naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }

  // Resize *y* to have the correct length.
//   y.resize( nmblock*nf*nza*naa*stokes_dim );
  // FIXME: this is an ugly solution, but works if appropriate checks are
  // preformed
  y.resize( nmblock*sensor_response.nrows() );

  // Create vector for MPB radiances for 1 measurement block.
  Vector ib( nf*nza*naa*stokes_dim );

  // Number of elements of *y* for one mblock
  Index    nblock = sensor_response.nrows();

  //--- Check input -----------------------------------------------------------

  if( check_input )
    {
      // Agendas (agendas not always used are checked when actually used)
      //
      chk_not_empty( "ppath_step_agenda", ppath_step_agenda );
      chk_not_empty( "rte_agenda", rte_agenda );

      // Stokes
      //
      chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

      // Basic checks of atmospheric, geoid and surface variables
      //  
      chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
      chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
      chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
      chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
      chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, 
                                                                    lon_grid );
      chk_atm_surface( "z_surface", z_surface, atmosphere_dim, lat_grid, 
                                                                    lon_grid );

      // Check that z_field has strictly increasing pages.
      //
      for( Index row=0; row<z_field.nrows(); row++ )
        {
          for( Index col=0; col<z_field.ncols(); col++ )
            {
              ostringstream os;
              os << "z_field (for latitude nr " << row << " and longitude nr " 
                 << col << ")";
              chk_if_increasing( os.str(), z_field(joker,row,col) ); 
            }
        }

      // Check that there is no gap between the surface and lowest pressure 
      // surface
      //
      for( Index row=0; row<z_surface.nrows(); row++ )
        {
          for( Index col=0; col<z_surface.ncols(); col++ )
            {
              if( z_surface(row,col)<z_field(0,row,col) ||
                       z_surface(row,col)>=z_field(z_field.npages()-1,row,col) )
                {
                  ostringstream os;
                  os << "The surface altitude (*z_surface*) cannot be outside "
                     << "of the altitudes in *z_field*.";
                  if( atmosphere_dim > 1 )
                    os << "\nThis was found to be the case for:\n"
                       << "latitude " << lat_grid[row];
                  if( atmosphere_dim > 2 )
                    os << "\nlongitude " << lon_grid[col];
                  throw runtime_error( os.str() );
                }
            }
        }

      // Cloud box
      //  
      chk_cloudbox( atmosphere_dim, p_grid, lat_grid, lon_grid,
                                                cloudbox_on, cloudbox_limits );

      // Frequency grid
      //
      if( nf == 0 )
        throw runtime_error( "The frequency grid is empty." );
      chk_if_increasing( "f_grid", f_grid );

      // Antenna
      //
      chk_if_in_range( "antenna_dim", antenna_dim, 1, 2 );
      if( nza == 0 )
        throw runtime_error( 
                         "The measurement block zenith angle grid is empty." );
      chk_if_increasing( "mblock_za_grid", mblock_za_grid );
      if( antenna_dim == 1 )
        {
          if( mblock_aa_grid.nelem() != 0 )
            throw runtime_error( 
              "For antenna_dim = 1, the azimuthal angle grid must be empty." );
        }
      else
        {
          if( atmosphere_dim < 3 )
            throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
          if( mblock_aa_grid.nelem() == 0 )
            throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
          chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
        }

      // Sensor
      //
      if( apply_sensor ) {
        if( sensor_response.ncols() != nf * nza * naa * stokes_dim ) {
          ostringstream os;
          os << "The *sensor_response* matrix does not have the right size, \n"
             << "either the method *sensor_responseInit* has not been run \n"
             << "prior to the call to *RteCalc* or some of the other sensor\n"
             << "response methods has not been correctly configured.";
          throw runtime_error( os.str() );
        }
      if( sensor_rot.nelem() != nmblock )
          throw runtime_error( 
                       "The length of *sensor_rot* must match *sensor_pos*." );
      }

      // Sensor position and LOS. 
      //
      // That the angles are inside OK ranges are checked inside ppathCalc.
      //
      if( sensor_pos.ncols() != atmosphere_dim )
        throw runtime_error( "The number of columns of sensor_pos must be "
                                  "equal to the atmospheric dimensionality." );
      if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
        throw runtime_error( 
                          "For 1D and 2D, sensor_los shall have one column." );
      if( atmosphere_dim == 3  &&  sensor_los.ncols() != 2 )
        throw runtime_error( "For 3D, sensor_los shall have two columns." );
      if( sensor_los.nrows() != nmblock )
        {
          ostringstream os;
          os << "The number of rows of sensor_pos and sensor_los must be "
             << "identical, but sensor_pos has " << nmblock << " rows,\n"
             << "while sensor_los has " << sensor_los.nrows() << " rows.";
          throw runtime_error( os.str() );
        }
    }
  //--- End: Check input ------------------------------------------------------


  // Loop:  measurement block / zenith angle / azimuthal angle
  //
  Index    nydone = 0;                 // Number of positions in y done
  Index    nbdone;                     // Number of positions in ib done
  Vector   los( sensor_los.ncols() );  // LOS of interest
  Index    iaa, iza, ip;
  //
  for( Index mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      nbdone = 0;

      for( iza=0; iza<nza; iza++ )
        {
          for( iaa=0; iaa<naa; iaa++ )
            {
              // Argument for verbosity of agendas
              bool  ag_verb = !( agenda_verb == 0  && 
                                             (iaa + iza + mblock_index) == 0 );

              // LOS of interest
              los     = sensor_los( mblock_index, joker );
              los[0] += mblock_za_grid[iza];
              if( antenna_dim == 2 )
                { los[1] += mblock_aa_grid[iaa]; }

              // Determine propagation path
              ppath_calc( ppath, ppath_step, ppath_step_agenda, 
                          atmosphere_dim, p_grid, lat_grid, lon_grid, 
                          z_field, r_geoid, z_surface,
                          cloudbox_on, cloudbox_limits,
                          sensor_pos(mblock_index,joker), los, 1, ag_verb );

              // Determine the radiative background
              get_radiative_background( i_rte, ppath_step, rte_pos, rte_los, 
                      rte_gp_p, rte_gp_lat, rte_gp_lon, i_space, 
                      surface_emission, surface_los, surface_refl_coeffs, ppath, 
                      ppath_step_agenda, rte_agenda, 
                      i_space_agenda, surface_agenda, atmosphere_dim, 
                      p_grid, lat_grid, lon_grid, z_field, t_field, 
                      r_geoid, z_surface, 
                      cloudbox_on, cloudbox_limits, scat_i_p, scat_i_lat, 
                      scat_i_lon, scat_za_grid, scat_aa_grid, f_grid, 
                                        stokes_dim, ag_verb, scat_za_interp );

	      
              // Execute the *rte_agenda*
              rte_agenda.execute( ag_verb );

              // If the sensor should be applied *i_rte* is reshaped
              // to a column vector to match *sensor_response*
              if( apply_sensor ) {
                for( ip=0; ip<stokes_dim; ip++ )
                  ib[Range(nbdone+ip,nf,stokes_dim)] = i_rte(joker,ip);

              // Increase nbdone
              nbdone += nf*stokes_dim;
              }
            }
        }

      // Apply sensor response matrix
      //
      if( apply_sensor )
        mult( y[Range(nydone,nblock)], sensor_response, ib );
      /* FIXME: Should *y* be created even if no sensor is applied?
         NOTE: sensorOff doesn't turn off apply_sensor.
      else
          // Copy ib to *y*
          y[Range(nydone,nblock)] = ib;
      */

      // Increase nydone
      nydone += nblock;
    }
}



//! rte_step
/*!
    Solves monochromatic VRTE for an atmospheric slab with constant 
    conditions.

    The function can be used for clearsky and cloudbox calculations.

    The function is best explained by considering a homogenous layer. That is,
    the physical conditions inside the layer are constant. In reality they
    are not constant, so in practical all coefficients have to be averaged 
    before calling this function. 
    Total extinction and absorption inside the layer are described by
    *ext_mat_av* and *abs_vec_av* respectively,
    the blackbdody radiation of the layer is given by *rte_planck_value*
    and the propagation path length through the layer is *l_step*.

    There is an additional scattering source term in the 
    VRTE, the scattering integral term. For this function a constant
    scattering term is assumed. The radiative transfer sterp is only a part 
    the iterative solution of the scattering problem, for more 
    information consider AUG. In the clearsky case this variable has to be
    set to 0.

    When calling the function, the vector *stokes_vec* shall contain the
    Stokes vector for the incoming radiation. The function returns this
    vector, then containing the outgoing radiation on the other side of the 
    layer.

    The function performs the calculations differently depending on the
    conditions to improve the speed. There are three cases: <br>
       1. Scalar absorption (stokes_dim = 1). <br>
       2. The matrix ext_mat_gas is diagonal (unpolarised absorption). <br>
       3. The total general case.

    \param   stokes_vec         Input/Output: A Stokes vector.
    \param   ext_mat_av         Input: Averaged extinction matrix.
    \param   abs_vec_av         Input: Averaged absorption vector.
    \param   sca_vec_av         Input: averaged scattering vector.
    \param   l_step             Input: The length of the RTE step.
    \param   rte_planck_value     Input: Blackbody radiation.

    \author Patrick Eriksson, Claudia Emde 
    \date   2002-11-22
*/
void rte_step(//Output and Input:
              VectorView stokes_vec,
              //Input
              ConstMatrixView ext_mat_av,
              ConstVectorView abs_vec_av,
              ConstVectorView sca_vec_av, 
              const Numeric& l_step,
              const Numeric& rte_planck_value )
{

  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim)); 
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( rte_planck_value >= 0 );
  assert( l_step >= 0 );

  // Check, if only the first component of abs_vec is non-zero:
  bool unpol_abs_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (abs_vec_av[i] != 0)
      unpol_abs_vec = false;
  
  bool unpol_sca_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (sca_vec_av[i] != 0)
      unpol_sca_vec = false;


  const Numeric   trans = exp(-ext_mat_av(0,0) * l_step);

  // Scalar case: 
  if( stokes_dim == 1 )
    { 
      stokes_vec[0] = stokes_vec[0] * trans + 
        (abs_vec_av[0] * rte_planck_value + sca_vec_av[0]) / ext_mat_av(0,0) 
        * (1 - trans);
    }
  // Vector case: 
    
      // We have here two cases, diagonal or non-diagonal ext_mat_gas
      // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
      // non-zero value in position 1.
    
      // Unpolarised
  else if( is_diagonal(ext_mat_av) && unpol_abs_vec && unpol_sca_vec )
    {
      // Stokes dim 1
      //    assert( ext_mat_av(0,0) == abs_vec_av[0] );
      //   Numeric transm = exp( -l_step * abs_vec_av[0] );
      stokes_vec[0] = stokes_vec[0] * trans + 
        (abs_vec_av[0] * rte_planck_value + sca_vec_av[0]) / ext_mat_av(0,0) 
        * (1 - trans);
      
      // Stokes dims > 1
      for( Index i=1; i<stokes_dim; i++ )
        {
          //      assert( abs_vec_av[i] == 0.);
          stokes_vec[i] = stokes_vec[i] * trans + 
            sca_vec_av[i] / ext_mat_av(i,i)  * (1 - trans);
        }
    }
  
  
  //General case
  else
    {
      stokes_vecGeneral(stokes_vec, ext_mat_av, abs_vec_av, sca_vec_av,
                        l_step, rte_planck_value);
    }
}




//! Calculate vector radiative transfer with fixed scattering integral
/*! 
  This function computes the radiative transfer for a thin layer.
  It is a general function which works for both, the vector and the scalar 
  RTE. 
  All coefficients and the scattered field vector are assumed to be constant
  inside the grid cell/layer.
  Then an analytic solution can be found (see AUG for details). 
  
  \param stokes_vec Output and Input: Stokes Vector after traversing a grid
                 cell/layer. 
  \param ext_mat_av Input: Extinction coefficient matrix.
  \param abs_vec_av Input: Absorption coefficient vector.
  \param sca_vec_av Input: Scattered field vector.
  \param l_step  Input: Pathlength through a grid cell/ layer.
  \param rte_planck_value  Input: Planck function.

  \author Claudia Emde
  \date 2002-06-08
*/
void stokes_vecGeneral(
                  //WS Output and Input:
                  VectorView stokes_vec,
                  //Input
                  ConstMatrixView ext_mat_av,
                  ConstVectorView abs_vec_av,
                  ConstVectorView sca_vec_av, 
                  const Numeric& l_step,
                  const Numeric& rte_planck_value )
{ 

  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim)); 
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( rte_planck_value >= 0 );
  // assert( l_step > 0 );

  //Initialize internal variables:

  // Matrix LU used for LU decompostion and as dummy variable:
  Matrix LU(stokes_dim, stokes_dim); 
  ArrayOfIndex indx(stokes_dim); // index for pivoting information 
  Vector b(stokes_dim); // dummy variable 
  Vector x(stokes_dim); // solution vector for K^(-1)*b
  Matrix I(stokes_dim, stokes_dim);

  Vector B_abs_vec(stokes_dim);
  B_abs_vec = abs_vec_av;
  B_abs_vec *= rte_planck_value;
  
  for (Index i=0; i<stokes_dim; i++) 
    b[i] = B_abs_vec[i] + sca_vec_av[i];  // b = abs_vec * B + sca_vec

  // solve K^(-1)*b = x
  ludcmp(LU, indx, ext_mat_av);
  lubacksub(x, LU, b, indx);

  Matrix ext_mat_ds(stokes_dim, stokes_dim);
  ext_mat_ds = ext_mat_av;
  ext_mat_ds *= -l_step; // ext_mat_ds = -ext_mat*ds
  
  Index q = 10;  // index for the precision of the matrix exponential function
  Matrix exp_ext_mat(stokes_dim, stokes_dim);
  matrix_exp(exp_ext_mat, ext_mat_ds, q);

  Vector term1(stokes_dim);
  Vector term2(stokes_dim);
  
  id_mat(I);
  for(Index i=0; i<stokes_dim; i++)
    {
      for(Index j=0; j<stokes_dim; j++)
        LU(i,j) = I(i,j) - exp_ext_mat(i,j); // take LU as dummy variable
    }
  mult(term2, LU, x); // term2: second term of the solution of the RTE with
                      //fixed scattered field

  
  mult(term1, exp_ext_mat, stokes_vec);// term1: first term of the solution of
                                    // the RTE with fixed scattered field
  
  for (Index i=0; i<stokes_dim; i++) 
    stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector

}

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
#include "rte.h"
#include "special_interp.h"
#include "lin_alg.h"

//! get_radiative_background
/*!
    Sets *i_rte* to the radiative background for a propagation path.

    The function uses *ppath* to determine the radiative background for a 
    propagation path. The possible backgrounds are listed in the header of
    the function ppath_set_background (in ppath.cc).

    The input variables beside *ppath* are the variables needed to calculate
    the magnitude of the radiative background. 

    The main purpose of the function is set *i_rte*. It is NOT needed to set 
    *i_rte* to the correct size before calling the function. The size is set
    to be [f_grid.nelem(),stokes_dim].

    \param   i_rte              Output: As the WSV with the same name.
    \param   ppath_step         Output: As the WSV with the same name.
    \param   a_pos              Output: As the WSV with the same name.
    \param   a_los              Output: As the WSV with the same name.
    \param   a_gp_p             Output: As the WSV with the same name.
    \param   a_gp_lat           Output: As the WSV with the same name.
    \param   a_gp_lon           Output: As the WSV with the same name.
    \param   i_space            Output: As the WSV with the same name.
    \param   ground_emission    Output: As the WSV with the same name.
    \param   ground_los         Output: As the WSV with the same name.
    \param   ground_refl_coeffs Output: As the WSV with the same name.
    \param   ppath              Input: As the WSV with the same name.
    \param   mblock_index       Input: As the WSV with the same name.
    \param   ppath_step_agenda  Input: As the WSV with the same name.
    \param   rte_agenda         Input: As the WSV with the same name.
    \param   i_space_agenda     Input: As the WSV with the same name.
    \param   ground_refl_agenda Input: As the WSV with the same name.
    \param   atmosphere_dim     Input: As the WSV with the same name.
    \param   p_grid             Input: As the WSV with the same name.
    \param   lat_grid           Input: As the WSV with the same name.
    \param   lon_grid           Input: As the WSV with the same name.
    \param   z_field            Input: As the WSV with the same name.
    \param   r_geoid            Input: As the WSV with the same name.
    \param   z_ground           Input: As the WSV with the same name.
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
              Vector&         a_pos,
              Vector&         a_los,
              GridPos&        a_gp_p,
              GridPos&        a_gp_lat,
              GridPos&        a_gp_lon,
              Matrix&         i_space,
              Matrix&         ground_emission, 
              Matrix&         ground_los, 
              Tensor4&        ground_refl_coeffs,
        const Ppath&          ppath,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         i_space_agenda,
        const Agenda&         ground_refl_agenda,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Tensor7&        scat_i_p,
        const Tensor7&        scat_i_lat,
        const Tensor7&        scat_i_lon,
        ConstVectorView       scat_za_grid,
        ConstVectorView       scat_aa_grid,
        ConstVectorView       f_grid,
        const Index&          stokes_dim,
        const Index&          antenna_dim,
        const Index&          agenda_verb )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index np      = ppath.np;

  // Resize i_rte to have the correct the size
  i_rte.resize( nf, stokes_dim );

  // Set a_pos, a_gp_XXX and a_los to match the last point in ppath.
  // The agendas below use different combinations of these varaibles.
  //
  // Note that the Ppath positions (ppath.pos) for 1D have one column more
  // than expected by most functions. Only the first atmosphere_dim values
  // shall be copied.
  //
  a_pos.resize( atmosphere_dim );
  a_pos = ppath.pos(np-1,Range(0,atmosphere_dim));
  a_los.resize( ppath.los.ncols() );
  a_los = ppath.los(np-1,joker);
  gridpos_copy( a_gp_p, ppath.gp_p[np-1] );
  if( atmosphere_dim > 1 )
    { gridpos_copy( a_gp_lat, ppath.gp_lat[np-1] ); }
  if( atmosphere_dim > 2 )
    { gridpos_copy( a_gp_lon, ppath.gp_lon[np-1] ); }


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


    case 2:   //--- The ground -----------------------------------------------
      {
        chk_not_empty( "ground_refl_agenda", ground_refl_agenda );
        ground_refl_agenda.execute( agenda_verb );

        // Check returned variables
        if( ground_emission.nrows() != nf  ||  
                                        ground_emission.ncols() != stokes_dim )
          throw runtime_error(
                  "The size of the created *ground_emission* is not correct.");

        // Copy to local variables as the ground variables can be changed by
        // by call of RteCalc. There is no need to copy *ground_los*.
        Index     nlos = ground_los.nrows();
        Matrix    ground_emission_local(nf,stokes_dim);
        Tensor4   ground_refl_coeffs_local(nlos,nf,stokes_dim,stokes_dim);
        //
        ground_emission_local    = ground_emission;
        ground_refl_coeffs_local = ground_refl_coeffs;


        // Calculate the spectra hitting the ground, if any reflection
        // directions have been set (*ground_los*). If *ground_los* is empty,
        // we are ready.
        //
        if( nlos == 0 )
          {
            i_rte = 0;    // This operation is also performed inside
          }               // the else-block
        else
          {
            // Set zenith and azimuthal angles (note that these variables 
            // are local variables). Each ground los is here treated as a 
            // measurement block, with no za and aa grids.

            Matrix   sensor_pos( nlos, a_pos.nelem() );
            Matrix   sensor_los( nlos, a_los.nelem() );
            for( Index ilos=0; ilos < nlos; ilos++ )
              {
                sensor_pos( ilos, joker ) = a_pos;
                sensor_los( ilos, joker ) = ground_los( ilos, joker);
              }
            Vector   mblock_za_grid(1,0);
            Vector   mblock_aa_grid(0);
            Index    mblock_index_local;
            if( antenna_dim > 1 )
              {
                mblock_aa_grid.resize(1);
                mblock_aa_grid = 0;
              }

            // Use some local variables to avoid unwanted  side effects
            Ppath    ppath_local;
            Matrix   y_rte_local;

            RteCalc( y_rte_local, ppath_local, ppath_step, i_rte, 
                 mblock_index_local, a_pos, a_los, a_gp_p, a_gp_lat, a_gp_lon,
                 i_space, ground_emission, ground_los, ground_refl_coeffs, 
                 ppath_step_agenda, rte_agenda, i_space_agenda, 
                 ground_refl_agenda, atmosphere_dim, p_grid, lat_grid, 
                 lon_grid, z_field, t_field, r_geoid, z_ground, 
                 cloudbox_on, cloudbox_limits, scat_i_p, scat_i_lat,
                 scat_i_lon, scat_za_grid, scat_aa_grid, sensor_pos, 
                 sensor_los, f_grid, stokes_dim, antenna_dim, 
                 mblock_za_grid, mblock_aa_grid );


            // Add the the calculated spectra (y_rte_local) to *i_rte*, 
            // considering the reflection coeff. matrix.
            //
            i_rte = 0;
            //
            for( Index ilos=0; ilos < nlos; ilos++ )
              {
                Index   i0 = ilos*nf;
            
                for( Index iv=0; iv<nf; iv++ )
                  {
                    if( stokes_dim == 1 )
                      {
                        i_rte(iv,0) += ground_refl_coeffs(ilos,iv,0,0) * 
                                                          y_rte_local(i0+iv,0);
                      }
                    else
                      {
                        Vector stokes_vec(stokes_dim);
                        mult( stokes_vec, 
                              ground_refl_coeffs(ilos,iv,joker,joker), 
                              y_rte_local(i0+iv,joker) );
                        for( Index is=0; is < stokes_dim; is++ )
                          { i_rte(iv,is) += stokes_vec[is]; }
                      }
                  }
              }
          }
        
        // Add the ground emission (the local copy) to *i_rte*
        for( Index iv=0; iv<nf; iv++ )
          {
            for( Index is=0; is < stokes_dim; is++ )
              { i_rte(iv,is) += ground_emission_local(iv,is); }
          }
      }
      break;


    case 3:   //--- Cloudbox surface -----------------------------------------

      {
        CloudboxGetOutgoing( i_rte, "i_rte", scat_i_p, scat_i_lat, scat_i_lon, 
                           a_gp_p, a_gp_lat, a_gp_lon, a_los, 
                           cloudbox_on, cloudbox_limits, atmosphere_dim,
                              stokes_dim, scat_za_grid, scat_aa_grid, f_grid );
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



//! ground_specular_los
/*!
    Calculates the LOS for a specular ground reflection.

    The function calculates the LOS for the downwelling radiation to consider
    when the ground is treated to give no scattering. The tilt of the ground
    (that is, a change in radius as a function of latitude etc.) is considered.

    \param   los                Output: As the WSV with the same name.
    \param   atmosphere_dim     Input: As the WSV with the same name.
    \param   r_geoid            Input: As the WSV with the same name.
    \param   z_ground           Input: As the WSV with the same name.
    \param   lat_grid           Input: As the WSV with the same name.
    \param   lon_grid           Input: As the WSV with the same name.
    \param   a_gp_lat           Input: As the WSV with the same name.
    \param   a_gp_lon           Input: As the WSV with the same name.
    \param   a_los              Input: As the WSV with the same name.

    \author Patrick Eriksson 
    \date   2002-09-22
*/
void ground_specular_los(
              VectorView   los,
        const Index&       atmosphere_dim,
        ConstMatrixView    r_geoid,
        ConstMatrixView    z_ground,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        const GridPos&     a_gp_lat,
        const GridPos&     a_gp_lon,
        ConstVectorView    a_los )
{
  assert( atmosphere_dim >= 1  &&  atmosphere_dim <= 3 );

  if( atmosphere_dim == 1 )
    {
      assert( r_geoid.nrows() == 1 );
      assert( r_geoid.ncols() == 1 );
      assert( z_ground.nrows() == 1 );
      assert( z_ground.ncols() == 1 );
      assert( a_los.nelem() == 1 );
      assert( a_los[0] > 90 );      // Otherwise ground refl. not possible
      assert( a_los[0] <= 180 ); 
      assert( los.nelem() == 1 );

      los[0] = 180 - a_los[0];
    }

  else if( atmosphere_dim == 2 )
    {
      assert( r_geoid.ncols() == 1 );
      assert( z_ground.ncols() == 1 );
      assert( a_los.nelem() == 1 );
      assert( los.nelem() == 1 );
      assert( abs(a_los[0]) <= 180 ); 

      
      // Calculate LOS neglecting any tilt of the ground
      los[0] = sign( a_los[0] ) * 180 - a_los[0];

      // Interpolate to get radius for the ground at reflection point.
      const Numeric r_ground =
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                              r_geoid, "r_geoid", a_gp_lat, a_gp_lon ) +
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                                    z_ground, "z_ground", a_gp_lat, a_gp_lon );

      // Calculate ground slope (unit is m/deg).
      Numeric slope = psurface_slope_2d( lat_grid, r_geoid(joker,0), 
                                       z_ground(joker,0), a_gp_lat, a_los[0] );

      // Calculate ground (angular) tilt (unit is deg).
      Numeric tilt = psurface_angletilt( r_ground, slope );

      // Check that a_los contains a downward LOS
      assert( is_los_downwards( a_los[0], tilt ) );
      
      // Include ground tilt
      los[0] -= 2 * tilt;
    }

  else if( atmosphere_dim == 3 )
    {
      assert( a_los.nelem() == 2 );
      assert( los.nelem() == 2 );
      assert( a_los[0] >= 0 ); 
      assert( a_los[0] <= 180 ); 
      assert( abs( a_los[1] ) <= 180 ); 

      assert( a_gp_lat.idx >= 0 );
      assert( a_gp_lat.idx <= ( lat_grid.nelem() - 2 ) );
      assert( a_gp_lon.idx >= 0 );
      assert( a_gp_lon.idx <= ( lon_grid.nelem() - 2 ) );

      // Calculate LOS neglecting any tilt of the ground
      los[0] = 180 - a_los[0];
      los[1] = a_los[1];

      // Below you find a first version to include the ground tilt.
      // Is the solution correct? 
      // At least, it does not work for za = 180.
      if( 0 )
        {
      // Interpolate to get radius for the ground at reflection point.
      const Numeric r_ground =
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                              r_geoid, "r_geoid", a_gp_lat, a_gp_lon ) +
          interp_atmsurface_by_gp( atmosphere_dim, lat_grid, lon_grid,
                                    z_ground, "z_ground", a_gp_lat, a_gp_lon );

      // Restore latitude and longitude values
      Vector   itw(2);
      Numeric  lat, lon;
      interpweights( itw, a_gp_lat );
      lat = interp( itw, lat_grid, a_gp_lat );
      interpweights( itw, a_gp_lon );
      lon = interp( itw, lon_grid, a_gp_lon );

      // Calculate ground slope along the viewing direction (unit is m/deg).
      //
      Index   ilat = gridpos2gridrange( a_gp_lat, abs( a_los[1] ) <= 90 ); 
      Index   ilon = gridpos2gridrange( a_gp_lon, a_los[1] >= 0 );
      //
      Numeric   lat1 = lat_grid[ilat];
      Numeric   lat3 = lat_grid[ilat+1];
      Numeric   lon5 = lon_grid[ilon];
      Numeric   lon6 = lon_grid[ilon+1];
      Numeric   r15  = r_geoid(ilat,ilon) + z_ground(ilat,ilon);
      Numeric   r35  = r_geoid(ilat+1,ilon) + z_ground(ilat+1,ilon);
      Numeric   r16  = r_geoid(ilat,ilon+1) + z_ground(ilat,ilon+1);
      Numeric   r36  = r_geoid(ilat+1,ilon+1) + z_ground(ilat+1,ilon+1);
      //
      Numeric slope = psurface_slope_3d( lat1, lat3, lon5, lon6, 
                                      r15, r35, r36, r16, lat, lon, a_los[1] );

      // Calculate ground (angular) tilt (unit is deg).
      Numeric tilt = psurface_angletilt( r_ground, slope );

      // Check that a_los contains a downward LOS
      assert( is_los_downwards( a_los[0], tilt ) );
      
      // Include ground tilt in zenith angle
      los[0] -= 2 * tilt;

      // Calculate ground slope along the viewing direction (unit is m/deg).
      Numeric   aa = a_los[1] + 90;
      if( aa > 180 )
        { aa -= 360; }
      //
      slope = psurface_slope_3d( lat1, lat3, lon5, lon6, 
                                            r15, r35, r36, r16, lat, lon, aa );

      // Calculate ground (angular) tilt (unit is deg).
      tilt = psurface_angletilt( r_ground, slope );

      // Include ground tilt in azimuth angle
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
    ppathCalc.

    All variables, beside *check_input*, corresponde to the WSV with the
    same name.

    \author Patrick Eriksson 
    \date   2003-01-18
*/
void rte_calc(
              Matrix&         y_rte,
              Ppath&          ppath,
              Ppath&          ppath_step,
              Matrix&         i_rte,
              Index&          mblock_index,
              Vector&         a_pos,
              Vector&         a_los,
              GridPos&        a_gp_p,
              GridPos&        a_gp_lat,
              GridPos&        a_gp_lon,
              Matrix&         i_space,
              Matrix&         ground_emission, 
              Matrix&         ground_los, 
              Tensor4&        ground_refl_coeffs,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         i_space_agenda,
        const Agenda&         ground_refl_agenda,
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
        const Tensor7&        scat_i_p,
        const Tensor7&        scat_i_lat,
        const Tensor7&        scat_i_lon,
        const Vector&         scat_za_grid,
        const Vector&         scat_aa_grid,
        const Matrix&         sensor_pos,
        const Matrix&         sensor_los,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Index&          antenna_dim,
        const Vector&         mblock_za_grid,
        const Vector&         mblock_aa_grid,
        const bool&           check_input )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index nmblock = sensor_pos.nrows();
  const Index nza     = mblock_za_grid.nelem();


  //--- Check input -----------------------------------------------------------

  if( check_input )
    {

      // Agendas (agendas not always used are checked when actually used)
      chk_not_empty( "ppath_step_agenda", ppath_step_agenda );
      chk_not_empty( "rte_agenda", rte_agenda );

      // Basic checks of atmospheric, geoid and ground variables
      //  
      chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
      chk_atm_grids( atmosphere_dim, p_grid, lat_grid, lon_grid );
      chk_atm_field( "z_field", z_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
      chk_atm_field( "t_field", t_field, atmosphere_dim, p_grid, lat_grid, 
                                                                    lon_grid );
      chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, 
                                                                    lon_grid );
      chk_atm_surface( "z_ground", z_ground, atmosphere_dim, lat_grid, 
                                                                    lon_grid );

      // Check that z_field has strictly increasing pages.
      //
      for( Index row=0; row<z_field.nrows(); row++ )
        {
          for( Index col=0; col<z_field.ncols(); col++ )
            {
              // I could not get the compliler to accept a solution without 
              // dummy!!
              Vector dummy(z_field.npages());
              dummy = z_field(joker,row,col);
              ostringstream os;
              os << "z_field (for latitude nr " << row << " and longitude nr " 
                 << col << ")";
              chk_if_increasing( os.str(), dummy ); 
            }
        }

      // Check that there is no gap between the ground and lowest pressure 
      // surface
      //
      for( Index row=0; row<z_ground.nrows(); row++ )
        {
          for( Index col=0; col<z_ground.ncols(); col++ )
            {
              if( z_ground(row,col)<z_field(0,row,col) || 
                       z_ground(row,col)>=z_field(z_field.npages()-1,row,col) )
                {
                  ostringstream os;
                  os << "The ground altitude (*z_ground*) cannot be outside "
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

      // Stokes
      //
      chk_if_in_range( "stokes_dim", stokes_dim, 1, 4 );

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


  // Number of azimuthal direction for pencil beam calculations
  Index naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }
  
  // Resize *y_rte* to have the correct size.
  // This size must be changed later when Hb is introduced
  y_rte.resize( nmblock * nf * nza * naa, stokes_dim );

  // Create a matrix to hold the spectra for one measurement block
  Matrix ib( nf * nza * naa, stokes_dim );

  // Loop:  measurement block / zenith angle / azimuthal angle
  //
  Index    nydone = 0;                 // Number of positions in y_rte done
  Index    nbdone;                     // Number of positions in ib done
  Index    nblock = nf*nza*naa;        // Number of spectral values of 1 block
  Vector   los( sensor_los.ncols() );  // LOS of interest
  Index    iaa, iza;
  //
  for( mblock_index=0; mblock_index<nmblock; mblock_index++ )
    {
      nbdone = 0;

      for( iza=0; iza<nza; iza++ )
        {
          for( iaa=0; iaa<naa; iaa++ )
            {
              // Argument for verbosity of agendas
              const Index  agenda_verb = iaa + iza + mblock_index;

              // LOS of interest
              los     = sensor_los( mblock_index, joker );
              los[0] += mblock_za_grid[iza];
              if( antenna_dim == 2 )
                { los[1] += mblock_aa_grid[iaa]; }

              // Determine propagation path
              ppath_calc( ppath, ppath_step, ppath_step_agenda, 
                          atmosphere_dim, p_grid, lat_grid, lon_grid, 
                          z_field, r_geoid, z_ground,
                          cloudbox_on, cloudbox_limits, 
                          sensor_pos(mblock_index,joker), los, agenda_verb );

              // Determine the radiative background
              get_radiative_background( i_rte, ppath_step, 
                      a_pos, a_los, a_gp_p, a_gp_lat, a_gp_lon, i_space, 
                      ground_emission, ground_los, ground_refl_coeffs, ppath, 
                      ppath_step_agenda, rte_agenda, 
                      i_space_agenda, ground_refl_agenda, atmosphere_dim, 
                      p_grid, lat_grid, lon_grid, z_field, t_field, 
                      r_geoid, z_ground, 
                      cloudbox_on, cloudbox_limits, scat_i_p, scat_i_lat, 
                      scat_i_lon, scat_za_grid, scat_aa_grid, f_grid, 
                      stokes_dim, antenna_dim, agenda_verb );

              // Execute the *rte_agenda*
              rte_agenda.execute( agenda_verb );
              
              // Copy i_rte to ib
              ib(Range(nbdone,nf),joker) = i_rte;

              // Increase nbdone
              nbdone += nf;
            }
        }

      // Copy ib to y_rte
      // The matrix Hb shall be applied here
      y_rte(Range(nydone,nblock),joker) = ib;

      // Increase nbdone
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
    the blackbdody radiation of the layer is given by *a_planck_value*
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
    \param   a_planck_value     Input: Blackbody radiation.

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
              const Numeric& a_planck_value )
{

  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim)); 
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( a_planck_value >= 0 );
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
        (abs_vec_av[0] * a_planck_value + sca_vec_av[0]) / ext_mat_av(0,0) 
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
        (abs_vec_av[0] * a_planck_value + sca_vec_av[0]) / ext_mat_av(0,0) 
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
                        l_step, a_planck_value);
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
  \param a_planck_value  Input: Planck function.

  \author Claudia Emde
  \date 2002-06-08
*/

void
stokes_vecGeneral(//WS Output and Input:
                  VectorView stokes_vec,
                  //Input
                  ConstMatrixView ext_mat_av,
                  ConstVectorView abs_vec_av,
                  ConstVectorView sca_vec_av, 
                  const Numeric& l_step,
                  const Numeric& a_planck_value )
{ 

  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim)); 
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( a_planck_value >= 0 );
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
  B_abs_vec *= a_planck_value; 
  
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

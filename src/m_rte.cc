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
  \file   m_rte.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-11 

  \brief  Workspace functions for solving clear sky radiative transfer.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <cmath>
#include <stdexcept>
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "physics_funcs.h"
#include "ppath.h"
#include "rte.h"
#include "special_interp.h"



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

//! RteCalc
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-09-17
*/
void RteCalc(
        // WS Output:
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
        // WS Input:
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
        const Vector&         mblock_aa_grid )
{
  // Some sizes
  const Index nf      = f_grid.nelem();
  const Index nmblock = sensor_pos.nrows();
  const Index nza     = mblock_za_grid.nelem();


  //--- Check input -----------------------------------------------------------

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
  chk_atm_surface( "r_geoid", r_geoid, atmosphere_dim, lat_grid, lon_grid );
  chk_atm_surface( "z_ground", z_ground, atmosphere_dim, lat_grid, lon_grid );

  // Check that z_field has strictly increasing pages.
  //
  for( Index row=0; row<z_field.nrows(); row++ )
    {
      for( Index col=0; col<z_field.ncols(); col++ )
        {
          // I could not get the compliler to accept a solution without dummy!!
          Vector dummy(z_field.npages());
          dummy = z_field(Range(joker),row,col);
          ostringstream os;
          os << "z_field (for latitude nr " << row << " and longitude nr " 
             << col << ")";
          chk_if_increasing( os.str(), dummy ); 
        }
    }

  // Check that there is no gap between the ground and lowest pressure surface
  //
  for( Index row=0; row<z_ground.nrows(); row++ )
    {
      for( Index col=0; col<z_ground.ncols(); col++ )
        {
          if( z_ground(row,col)<z_field(0,row,col) || 
                       z_ground(row,col)>=z_field(z_field.npages()-1,row,col) )
            {
              ostringstream os;
              os << "The ground altitude (*z_ground*) cannot be outside of the"
                 << " altitudes in *z_field*.";
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
    throw runtime_error( "The measurement block zenith angle grid is empty." );
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
        throw runtime_error( 
         "2D antennas (antenna_dim=2) can only be used with 3D atmospheres." );
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
    throw runtime_error( "The number of columns of sensor_pos must be equal "
                                        "to the atmospheric dimensionality." );
  if( atmosphere_dim <= 2  &&  sensor_los.ncols() != 1 )
    throw runtime_error( "For 1D and 2D, sensor_los shall have one column." );
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
              // LOS of interest
              los     = sensor_los( mblock_index, Range(joker) );
              los[0] += mblock_za_grid[iza];
              if( antenna_dim == 2 )
                { los[1] += mblock_aa_grid[iaa]; }

              // Determine propagation path
              ppathCalc( ppath, ppath_step, ppath_step_agenda, atmosphere_dim, 
                        p_grid, lat_grid, lon_grid, z_field, t_field, 
                              r_geoid, z_ground,
                              cloudbox_on, cloudbox_limits, 
                                  sensor_pos(mblock_index,Range(joker)), los );

              // Determine the radiative background
              get_radiative_background( i_rte, ppath_step, 
                      a_pos, a_los, a_gp_p, a_gp_lat, a_gp_lon, i_space, 
                      ground_emission, ground_los, ground_refl_coeffs, ppath, 
                      mblock_index, ppath_step_agenda, rte_agenda, 
                      i_space_agenda, ground_refl_agenda, atmosphere_dim, 
                      p_grid, lat_grid, lon_grid, z_field, t_field, 
                      r_geoid, z_ground, 
                      cloudbox_on, cloudbox_limits, scat_i_p, scat_i_lat, 
                      scat_i_lon, scat_za_grid, scat_aa_grid, f_grid, 
                      stokes_dim, antenna_dim );

              // Execute the *rte_agenda*
              rte_agenda.execute();
              
              // Copy i_rte to ib
              ib(Range(nbdone,nf),Range(joker)) = i_rte;
            
              // Increase nbdone
              nbdone += nf;
            }
        }

      // Copy ib to y_rte
      // The matrix Hb shall be applied here
      y_rte(Range(nydone,nblock),Range(joker)) = ib;

      // Increase nbdone
      nydone += nblock;
    } 
}



//! RteEmissionStd
/*! 
 
   This function does a clearsky radiative transfer calculation for a given 
   propagation path. 
   Gaseous emission and absorption is calculated for each propagation
   path point using the agenda *gas_absorption_agenda*. 
   Absorption vector and extinction matrix are created using 
   *opt_prop_part_agenda*.
   The coefficients for the radiative transfer are averaged between two
   successive propagation path points.
 

   \author Claudia Emde (first version by Patrick Eriksson)
   \date   2003-01-07
*/
void RteEmissionStd(
      // WS Output:
             Matrix&    i_rte,
             Matrix&    abs_vec,
             Tensor3&   ext_mat,
             Numeric&   a_pressure,
             Numeric&   a_temperature,
             Vector&    a_vmr_list,
             Index&     f_index,
             Matrix&    abs_scalar_gas,
       // WS Input:
       const Ppath&     ppath,
       const Vector&    f_grid,
       const Index&     stokes_dim,
       const Index&     atmosphere_dim,
       const Vector&    p_grid,
       const Vector&    lat_grid,
       const Vector&    lon_grid,
       const Tensor3&   t_field,
       const Tensor4&   vmr_field,
       const Agenda&    scalar_gas_absorption_agenda,
       const Agenda&    opt_prop_gas_agenda )
{
  // Relevant checks are assumed to be done in RteCalc

  // Some sizes
  const Index   nf = f_grid.nelem();
  const Index   np = ppath.np;
  // Number of species:
  const Index   ns = vmr_field.nbooks();


  // If the number of propagation path points is 0 or 1, we are already ready,
  // the observed spectrum equals then the radiative background.
  if( np > 1 )
    {
      // Determine the pressure at each propagation path point
      Vector   p_ppath(np);
      Matrix   itw_p(np,2);
      //
      interpweights( itw_p, ppath.gp_p );      
      itw2p( p_ppath, p_grid, ppath.gp_p, itw_p );

      // Determine the atmospheric temperature at each propagation path point
      Vector t_ppath(np);
      //
      interp_atmfield_by_itw( t_ppath,  atmosphere_dim, p_grid, lat_grid, 
                               lon_grid, t_field, "t_field", 
                               ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_p );

      // Determine the vmrs at each propagation path point
      Matrix   vmr_ppath(ns,np), itw_field;
      //
      interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid, lat_grid, 
                            lon_grid, ppath.gp_p, ppath.gp_lat, ppath.gp_lon );

      // We have to loop over all gaseous species to obtain the 
      // full vmr_field.
      for( Index is = 0; is < ns; is ++)
        {
          interp_atmfield_by_itw( vmr_ppath(is, Range(joker)), atmosphere_dim,
                    p_grid, lat_grid, lon_grid, 
                    vmr_field(is, Range(joker), Range(joker),  Range(joker)), 
              "vmr_field", ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
        }
      
      // Calculate absorption vector and extinction matrix for all points 
      // along the propagation path.

      // Define a variable which holds the scalar gas absorption for 
      // all pressure values and all frequencies.
      // Dimensions: [ # frequencies, # pressure levels, # species ]
      Tensor3 scalar_gas_array;

      // Variables for extinction matrix and absorption vector at each 
      // propagation path point.
      ArrayOfTensor3   ext_mat_ppath(np);
      ArrayOfMatrix    abs_vec_ppath(np);
    
      for ( Index ip = 0; ip < np; ip ++)
        { 
          a_pressure    = p_ppath[ip];
          a_temperature = t_ppath[ip];
          a_vmr_list    = vmr_ppath(joker,ip);
          
          // If f_index < 0, scalar gas absorption is calculated for 
          // all frequencies in f_grid.
          f_index = -1;
          scalar_gas_absorption_agenda.execute();
          
          // We don't know how many gas species we have before the first
          // call of scalar_gas_absorption_agenda. So we have to resize
          // scalar_gas_array after the first call:
          //
          if ( ip == 0 )
            { scalar_gas_array.resize( nf, np, abs_scalar_gas.ncols() ); }
          //
          scalar_gas_array(joker,ip,joker) = abs_scalar_gas(joker,joker);
            

          // Calculate extinction matrix and absorption vector
          // for all propagation path points.
          
          abs_vec_ppath[ip].resize(nf,stokes_dim);
          ext_mat_ppath[ip].resize(nf,stokes_dim,stokes_dim);
          
          opt_prop_gas_agenda.execute();
              
          abs_vec_ppath[ip] = abs_vec;
          ext_mat_ppath[ip] = ext_mat;
          
        }

      
      // Loop the propagation path steps
      //
      // The number of path steps is np-1.
      // The path points are stored in such way that index 0 corresponds to
      // the point closest to the sensor.
      //
          
      // Define variables which hold averaged coefficients:
      Matrix   ext_mat_av(stokes_dim, stokes_dim,0.);
      Vector   abs_vec_av(stokes_dim,0.);
              
      for( Index ip=np-1; ip>0; ip-- )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              // Calculate averaged values for extinction matrix and 
              // absorption vector.
              for (Index i = 0; i < stokes_dim; i++)
                {
                  // Extinction matrix requires a second loop over 
                  // stokes_dim:
                  for (Index j = 0; j < stokes_dim; j++)
                    {
                      ext_mat_av(i, j) = 0.5*( ext_mat_ppath[ip](iv, i, j) +
                                                ext_mat_ppath[ip-1](iv, i, j));
                    }
                 
                  // Absorption vector
                  abs_vec_av[i] = 0.5*( abs_vec_ppath[ip](iv,i) +
                                                   abs_vec_ppath[ip-1](iv,i) );
                }
                  
              // Calculate an effective blackbody radiation for the step
              // The mean of the temperature at the end points is used.
              Numeric planck_value = 
                           planck( f_grid[iv], (t_ppath[ip]+t_ppath[ip-1])/2 );
                  
              // Dummy vector for scattering integral. It has to be 
              // set to 0 for clear sky calculations.
              Vector sca_vec_dummy(stokes_dim, 0.);

              assert (!is_singular( ext_mat_av ));   

              // Perform the RTE step.
              rte_step( i_rte(iv, Range(joker)), ext_mat_av, abs_vec_av,
                             sca_vec_dummy, ppath.l_step[ip-1], planck_value );
            }
        }
    }
}




//! yNoPolarisation
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Patrick Eriksson
   \date   2002-08-09
*/
void yNoPolarisation(
        // WS Output:
              Vector&         y,
        // WS Input:
        const Matrix&         y_rte )
{
  out2 << "   Creates *y* from *y_rte* assuming equal sensitivity to all "
       << "polarisations.\n";

  y.resize( y_rte.nrows() );

  y = y_rte( Range(joker), 0 );
}

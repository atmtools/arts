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
              Vector&         y,
              Ppath&          ppath,
              Ppath&          ppath_step,
              Matrix&         iy,
              Vector&         rte_pos,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Vector&         rte_los,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         iy_surface_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Sparse&         sensor_response,
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

  // Number of azimuthal direction for pencil beam calculations
  Index naa = mblock_aa_grid.nelem();
  if( antenna_dim == 1 )
    { naa = 1; }


  //--- Check input -----------------------------------------------------------
  //---------------------------------------------------------------------------

  // Agendas (agendas not always used are checked elsewhere when used)
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
        throw runtime_error( "2D antennas (antenna_dim=2) can only be "
                                                 "used with 3D atmospheres." );
      if( mblock_aa_grid.nelem() == 0 )
        throw runtime_error(
                      "The measurement block azimuthal angle grid is empty." );
      chk_if_increasing( "mblock_aa_grid", mblock_aa_grid );
    }

  // Sensor
  //
  if( sensor_response.ncols() != nf * nza * naa * stokes_dim ) 
    {
      ostringstream os;
      os << "The *sensor_response* matrix does not have the right size, \n"
         << "either the method *sensor_responseInit* has not been run \n"
         << "prior to the call to *RteCalc* or some of the other sensor\n"
         << "response methods has not been correctly configured.";
      throw runtime_error( os.str() );
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
  //---------------------------------------------------------------------------
  //--- End: Check input ------------------------------------------------------


  // Resize *y* to have correct length.
  y.resize( nmblock*sensor_response.nrows() );

  // Create vector for MPB radiances for 1 measurement block.
  Vector ib( nf*nza*naa*stokes_dim );

  // Number of elements of *y* for one mblock
  Index    nblock = sensor_response.nrows();


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
              bool  ag_verb = ( (iaa + iza + mblock_index) != 0 );

              // LOS of interest
              los     = sensor_los( mblock_index, joker );
              los[0] += mblock_za_grid[iza];
              if( antenna_dim == 2 )
                { los[1] += mblock_aa_grid[iaa]; }

              // Calculate *iy*
              iy_calc( iy, ppath, ppath_step, rte_pos, rte_gp_p, rte_gp_lat, 
                 rte_gp_lon, rte_los, 
                 ppath_step_agenda, rte_agenda, iy_space_agenda, 
                 iy_surface_agenda, iy_cloudbox_agenda, atmosphere_dim, 
                 p_grid, lat_grid, lon_grid, z_field, r_geoid, z_surface, 
                 cloudbox_on,  cloudbox_limits, sensor_pos(mblock_index,joker),
                 los, f_grid, stokes_dim, ag_verb );

              // Copy *iy* to *ib*
              for( ip=0; ip<stokes_dim; ip++ )
                { ib[Range(nbdone+ip,nf,stokes_dim)] = iy(joker,ip); }

              // Increase nbdone
              nbdone += nf*stokes_dim;
            }
        }

      // Apply sensor response matrix
      mult( y[Range(nydone,nblock)], sensor_response, ib );

      // Increase nydone
      nydone += nblock;
    }
}



//! RteEmissionStd
/*! 
   See the the online help (arts -d FUNCTION_NAME)

   \author Claudia Emde and Patrick Eriksson
   \date   2003-01-07
*/
void RteEmissionStd(
      // WS Output:
             Matrix&    iy,
             Matrix&    abs_vec,
             Tensor3&   ext_mat,
             Numeric&   rte_pressure,
             Numeric&   rte_temperature,
             Vector&    rte_vmr_list,
             Index&     f_index,
             Index&     ppath_index,
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

      // Determine the atmospheric temperature and species VMR at 
      // each propagation path point
      Vector   t_ppath(np);
      Matrix   vmr_ppath(ns,np), itw_field;
      rte_vmr_list.resize(ns);
      //
      interp_atmfield_gp2itw( itw_field, atmosphere_dim, p_grid, lat_grid, 
                            lon_grid, ppath.gp_p, ppath.gp_lat, ppath.gp_lon );
      //
      interp_atmfield_by_itw( t_ppath,  atmosphere_dim, p_grid, lat_grid, 
                               lon_grid, t_field, "t_field", 
                           ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
      // 
      for( Index is=0; is<ns; is++ )
        {
          interp_atmfield_by_itw( vmr_ppath(is, joker), atmosphere_dim,
                    p_grid, lat_grid, lon_grid, 
                    vmr_field(is, joker, joker,  joker), 
              "vmr_field", ppath.gp_p, ppath.gp_lat, ppath.gp_lon, itw_field );
        }

      // Calculate absorption vector and extinction matrix for all points 
      // along the propagation path.

      // Variables for extinction matrix and absorption vector at each 
      // propagation path point.
      ArrayOfTensor3   ext_mat_ppath(np);
      ArrayOfMatrix    abs_vec_ppath(np);

      // If f_index < 0, scalar gas absorption is calculated for 
      // all frequencies in f_grid.
      f_index = -1;
      
      // Loop the propagation path steps
      //
      // The number of path steps is np-1.
      // The path points are stored in such way that index 0 corresponds to
      // the point closest to the sensor.
      //
          
      // Define variables which hold averaged coefficients:
      Matrix   ext_mat_av(stokes_dim, stokes_dim,0.);
      Vector   abs_vec_av(stokes_dim,0.);

      // Dummy vector for scattering integral. It has to be 
      // set to 0 for clear sky calculations.
      Vector sca_vec_dummy(stokes_dim, 0.);
              
      for( Index ip=np-1; ip>0; ip-- )
        {
          rte_pressure    = 0.5*(p_ppath[ip] + p_ppath[ip-1]);
          rte_temperature = 0.5*(t_ppath[ip] + t_ppath[ip-1]);

          for( Index is = 0; is < ns; is ++)
            {
              rte_vmr_list[is] = 0.5*(vmr_ppath(is,ip) + vmr_ppath(is, ip-1));
            }
          
          // The absO2ZeemanModel needs the position in the propagation
          // path. 
          ppath_index = ip;
          
          scalar_gas_absorption_agenda.execute( ip );

          opt_prop_gas_agenda.execute( ip ); 

          for( Index iv=0; iv<nf; iv++ )
            {
              // Calculate averaged values for extinction matrix and 
              // absorption vector.
              ext_mat_av = ext_mat(iv, joker, joker);
              abs_vec_av = abs_vec(iv, joker);
              
              // Calculate an effective blackbody radiation for the step
              // The mean of the temperature at the end points is used.
              Numeric planck_value = 
                           planck( f_grid[iv], rte_temperature);
                  
              assert (!is_singular( ext_mat_av ));   

              // Perform the RTE step.
              rte_step_std( iy(iv, joker), ext_mat_av, abs_vec_av,
                             sca_vec_dummy, ppath.l_step[ip-1], planck_value );
	      
            }
        }
    }
}



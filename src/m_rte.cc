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
        const Sparse&         sensor_response,
        const Matrix&         sensor_pos,
        const Matrix&         sensor_los,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Index&          antenna_dim,
        const Vector&         mblock_za_grid,
        const Vector&         mblock_aa_grid )
{
  const bool   check_input = true;
  const bool   apply_sensor = true;

  rte_calc( y_rte, ppath, ppath_step, i_rte, mblock_index, a_pos, a_los,
            a_gp_p, a_gp_lat, a_gp_lon, i_space, 
            ground_emission, ground_los, ground_refl_coeffs, 
            ppath_step_agenda, rte_agenda, i_space_agenda, ground_refl_agenda,
            atmosphere_dim, p_grid, lat_grid, lon_grid, z_field, t_field, 
            r_geoid, z_ground, cloudbox_on,  cloudbox_limits, 
            scat_i_p, scat_i_lat, scat_i_lon, scat_za_grid, scat_aa_grid, 
            sensor_response, sensor_pos, sensor_los, f_grid, stokes_dim,
            antenna_dim, mblock_za_grid, mblock_aa_grid, 
            check_input, apply_sensor );
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
    
      for ( Index ip = 0; ip < np; ip ++)
        { 
          a_pressure    = p_ppath[ip];
          a_temperature = t_ppath[ip];
          a_vmr_list    = vmr_ppath(joker,ip);
          
          scalar_gas_absorption_agenda.execute( ip );
          
          // Calculate extinction matrix and absorption vector
          // for all propagation path points.
          
          abs_vec_ppath[ip].resize(nf,stokes_dim);
          ext_mat_ppath[ip].resize(nf,stokes_dim,stokes_dim);
          
          opt_prop_gas_agenda.execute( ip );
              
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

      // Dummy vector for scattering integral. It has to be 
      // set to 0 for clear sky calculations.
      Vector sca_vec_dummy(stokes_dim, 0.);
              
      for( Index ip=np-1; ip>0; ip-- )
        {
          for( Index iv=0; iv<nf; iv++ )
            {
              // Calculate averaged values for extinction matrix and 
              // absorption vector.
              for (Index i = 0; i < stokes_dim; i++)
                {
                  // Extinction matrix requires a second loop over stokes_dim:
                  for (Index j = 0; j < stokes_dim; j++)
                    { ext_mat_av(i, j) = 0.5*( ext_mat_ppath[ip](iv, i, j) +
                                              ext_mat_ppath[ip-1](iv, i, j)); }
                 
                  // Absorption vector
                  abs_vec_av[i] = 0.5*( abs_vec_ppath[ip](iv,i) +
                                                   abs_vec_ppath[ip-1](iv,i) );
                }
                  
              // Calculate an effective blackbody radiation for the step
              // The mean of the temperature at the end points is used.
              Numeric planck_value = 
                           planck( f_grid[iv], (t_ppath[ip]+t_ppath[ip-1])/2 );
                  
              assert (!is_singular( ext_mat_av ));   

              // Perform the RTE step.
              rte_step( i_rte(iv, joker), ext_mat_av, abs_vec_av,
                             sca_vec_dummy, ppath.l_step[ip-1], planck_value );
            }
        }
    }
}



//! RteEmissionTest
/*! 
   \author Patrick Eriksson
   \date   2003-04-06
*/
void RteEmissionTest(
      // WS Output:
             Matrix&    i_rte,
             Matrix&    abs_vec,
             Tensor3&   ext_mat,
             Numeric&   a_pressure,
             Numeric&   a_temperature,
             Vector&    a_vmr_list,
             Index&     f_index,
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

  assert( stokes_dim == 1 );

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
    
      for ( Index ip = 0; ip < np; ip ++)
        { 
          a_pressure    = p_ppath[ip];
          a_temperature = t_ppath[ip];
          a_vmr_list    = vmr_ppath(joker,ip);
          
          scalar_gas_absorption_agenda.execute( ip );
          
          // Calculate extinction matrix and absorption vector
          // for all propagation path points.
          
          abs_vec_ppath[ip].resize(nf,stokes_dim);
          ext_mat_ppath[ip].resize(nf,stokes_dim,stokes_dim);
          
          opt_prop_gas_agenda.execute( ip );
              
          abs_vec_ppath[ip] = abs_vec;
          ext_mat_ppath[ip] = ext_mat;
          
        }

      
      // Loop the propagation path steps
      //
      // The number of path steps is np-1.
      // The path points are stored in such way that index 0 corresponds to
      // the point closest to the sensor.
          
      // Define variables which hold averaged coefficients:
      Matrix   ext_mat_av(stokes_dim, stokes_dim,0.);
      Vector   abs_vec_av(stokes_dim,0.);

      // Dummy vector for scattering integral. It has to be 
      // set to 0 for clear sky calculations.
      Vector sca_vec_dummy(stokes_dim, 0.);
              

      for( Index iv=0; iv<nf; iv++ )
        {
          Numeric   b1 = planck( f_grid[iv], t_ppath[np-1] );

          for( Index ip=np-1; ip>0; ip-- )
            {
              Numeric   tau = ppath.l_step[ip-1] * 0.5 *
                       ( abs_vec_ppath[ip](iv,0) + abs_vec_ppath[ip-1](iv,0) );
                  
              Numeric   b2 = planck( f_grid[iv], t_ppath[ip-1] );

              if( tau < 0.005 )
                {
                  Numeric   xhalf   = tau / 2;
                  Numeric   xsquare = tau * tau;

                  i_rte(iv,0) = i_rte(iv,0) * exp( -tau ) + 
                    b1 * ( xhalf - xsquare/3 ) + b2 * ( xhalf - xsquare/6 );
                }
              else
                {
                  Numeric   trans = exp( -tau );

                  i_rte(iv,0) = i_rte(iv,0) * trans + 
                    b2 + (b1-b2) * (1-trans) / tau - b1 * trans;
                }

              b1 = b2;
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

  y = y_rte( joker, 0 );
}

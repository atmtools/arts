/* Copyright (C) 2002,2003 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
   USA. 
*/
  
/*!
  \file   scatrte.cc
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Wed Jun 04 11:03:57 2003
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox.
  
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "array.h"
#include "auto_md.h"
#include "matpackVII.h"
#include "ppath.h"
#include "agenda_class.h"
#include "physics_funcs.h"
#include "lin_alg.h"
#include "math_funcs.h"
#include "messages.h"
#include "xml_io.h"
#include "rte.h"
#include "special_interp.h"
#include "scatrte.h"
#include "logic.h"
#include "check_input.h"
#include "sorting.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;

//! cloud_fieldsCalc
/*! 
  Calculate ext_mat, abs_vec for all points inside the cloudbox for one 
  propagation direction.
  sca_vec can be obtained from the workspace variable doit_scat_field.
  As we need the average for each layer, it makes sense to calculte
  the coefficients once and store them in an array instead of 
  calculating at each point the coefficient of the point above and 
  the point below. 

  // Output
  \param ext_mat_field extinction matrix field
  \param abs_vec_field absorption vector field
  // Input
  \param spt_calc_agenda Agenda for calculation of single scattering properties
  \param opt_prop_part_agenda Agenda for summing over all hydrometeor species
  \param scat_za_index Indices for
  \param scat_aa_index    propagation direction
  \param cloudbox_limits Cloudbox limits.
  \param t_field Temperature field
  \param pnd_field Particle number density field. 

  \author Claudia Emde
  \date 2002-06-03
*/
void cloud_fieldsCalc(// Output and Input:
                        Tensor5View ext_mat_field,
                        Tensor4View abs_vec_field,
                        // Input:
                        const Agenda& spt_calc_agenda,
                        const Agenda& opt_prop_part_agenda,
                        const Index& scat_za_index, 
                        const Index& scat_aa_index,
                        const ArrayOfIndex& cloudbox_limits,
                        const Tensor3View t_field, 
                        const Tensor4View pnd_field
                        )
{
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from 
  // where this function is called.
  
  out3 << "Calculate scattering properties in cloudbox \n";
  
  const Index atmosphere_dim = cloudbox_limits.nelem()/2;
  const Index N_pt = pnd_field.nbooks();
  const Index stokes_dim = ext_mat_field.ncols(); 
  
  assert( atmosphere_dim == 1 || atmosphere_dim ==3 );
  assert( ext_mat_field.ncols() == ext_mat_field.nrows() &&
          ext_mat_field.ncols() == abs_vec_field.ncols());
  
  const Index Np_cloud = cloudbox_limits[1]-cloudbox_limits[0]+1;
  
  // If atmosohere_dim == 1
  Index Nlat_cloud = 1;
  Index Nlon_cloud = 1;
  
  if (atmosphere_dim == 3)
    {
      Nlat_cloud = cloudbox_limits[3]-cloudbox_limits[2]+1;
      Nlon_cloud = cloudbox_limits[5]-cloudbox_limits[4]+1;
    }
  
  // Initialize ext_mat(_spt), abs_vec(_spt)
  // Resize and initialize variables for storing optical properties
  // of cloud particles
  Matrix abs_vec_spt_local(N_pt, stokes_dim, 0.);
  Tensor3 ext_mat_spt_local(N_pt, stokes_dim, stokes_dim, 0.);
  Matrix abs_vec_local;
  Tensor3 ext_mat_local;
  Numeric rte_temperature_local;
  
  // Calculate ext_mat, abs_vec for all points inside the cloudbox.
  // sca_vec can be obtained from the workspace variable doit_scat_field.
  // As we need the average for each layer, it makes sense to calculte
  // the coefficients once and store them in an array instead of 
  // calculating at each point the coefficient of the point above and 
  // the point below. 
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:               
  
  // Loop over all positions inside the cloudbox defined by the 
  // cloudbox_limits.
  for(Index scat_p_index_local = 0; scat_p_index_local < Np_cloud; 
      scat_p_index_local ++)
    {
      for(Index scat_lat_index_local = 0; scat_lat_index_local < Nlat_cloud; 
          scat_lat_index_local ++)
        {
          for(Index scat_lon_index_local = 0; 
              scat_lon_index_local < Nlon_cloud; 
              scat_lon_index_local ++)
            {
              if (atmosphere_dim == 1)
                rte_temperature_local = 
                  t_field(scat_p_index_local + cloudbox_limits[0], 0, 0);
              else
                rte_temperature_local = 
                  t_field(scat_p_index_local + cloudbox_limits[0],
                          scat_lat_index_local + cloudbox_limits[2],
                          scat_lon_index_local + cloudbox_limits[4]);
              
              //Calculate optical properties for single particle types:
              //( Execute agendas silently. )
              spt_calc_agendaExecute(ext_mat_spt_local, 
                                     abs_vec_spt_local,
                                     scat_p_index_local,
                                     scat_lat_index_local,
                                     scat_lon_index_local, 
                                     rte_temperature_local,
                                     scat_za_index,
                                     scat_aa_index,
                                     spt_calc_agenda,
                                     true);

              opt_prop_part_agendaExecute(ext_mat_local, abs_vec_local, 
                                          ext_mat_spt_local, 
                                          abs_vec_spt_local,
                                          scat_p_index_local,
                                          scat_lat_index_local,
                                          scat_lon_index_local,
                                          opt_prop_part_agenda,
                                          true);
           
              // Store coefficients in arrays for the whole cloudbox.
              abs_vec_field(scat_p_index_local, scat_lat_index_local,
                            scat_lon_index_local,
                            joker) = abs_vec_local(0, joker);
              
              ext_mat_field(scat_p_index_local, scat_lat_index_local,
                            scat_lon_index_local,
                            joker, joker) = ext_mat_local(0, joker, joker);
            } 
        }
    }
}
  




//! cloud_ppath_update1D
/*! 
  This function calculates the radiation field along a propagation path 
  step for specified zenith direction and pressure level.
  This function is used in the sequential update and called inside a loop over
  the pressure grid. 
  In the function the intersection point of the propagation path with the 
  next layer is calculated and all atmospheric properties are 
  interpolated an the intersection point. Then a radiative transfer step is 
  performed.

  WS Output:
  \param doit_i_field Updated radiation field inside the cloudbox. 
  WS Input:
  \param p_index // Pressure index
  \param scat_za_index // Index for proagation direction
  \param scat_za_grid
  \param cloudbox_limits 
  \param doit_scat_field Scattered field.
  Calculate scalar gas absorption:
  \param abs_scalar_gas_agenda
  \param vmr_field
  Scalar gas absorption:
  \param opt_prop_gas_agenda
  Propagation path calculation:
  \param ppath_step_agenda
  \param p_grid
  \param z_field
  \param r_geoid
  Calculate thermal emission:
  \param t_field
  \param f_grid
  \param f_index
  Optical properties of particles
  \param ext_mat_field
  \param abs_vec_field
  \param surface_prop_agenda
  \param scat_za_interp

  \author Claudia Emde
  \date 2003-06-04
*/
void cloud_ppath_update1D(
                          // Input and output
                          Tensor6View doit_i_field,
                          // ppath_step_agenda:
                          const Index& p_index,
                          const Index& scat_za_index,
                          ConstVectorView scat_za_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstTensor6View doit_scat_field,
                          // Calculate scalar gas absorption:
                          const Agenda& abs_scalar_gas_agenda,
                          ConstTensor4View vmr_field,
                          // Gas absorption: 
                          const Agenda& opt_prop_gas_agenda,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          ConstVectorView  p_grid,
                          ConstTensor3View z_field,
                          ConstMatrixView r_geoid,
                          ConstMatrixView z_surface,
                          // Calculate thermal emission:
                          ConstTensor3View t_field,
                          ConstVectorView f_grid,
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          const Agenda& surface_prop_agenda,
                          //const Agenda& surface_prop_agenda, 
                          const Index& scat_za_interp)
{
  Matrix iy;
  Matrix surface_emission;
  Matrix surface_los;
  Tensor4 surface_rmatrix;
  Ppath ppath_step;
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from 
  // where this function is called.
  
  //Initialize ppath for 1D.
  ppath_init_structure(ppath_step, 1, 1);
  // See documentation of ppath_init_structure for understanding
  // the parameters.
  
  // Assign value to ppath.pos:
  ppath_step.z[0]     = z_field(p_index,0,0);
  ppath_step.pos(0,0) = r_geoid(0,0) + ppath_step.z[0];
  
  // Define the direction:
  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
  
  
  // Define the grid positions:
  ppath_step.gp_p[0].idx   = p_index;
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;
  
  // Call ppath_step_agenda: 
  Vector unused_lat_grid(0);
  Vector unused_lon_grid(0);
  ppath_step_agendaExecute(ppath_step, 1, p_grid,
                           unused_lat_grid, unused_lon_grid,
                           z_field, r_geoid, z_surface,
                           ppath_step_agenda, true);
  
  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  
  if((cloudbox_limits[0] <= ppath_step.gp_p[1].idx) &&
     cloudbox_limits[1] > ppath_step.gp_p[1].idx ||
     (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
      abs(ppath_step.gp_p[1].fd[0]) < 1e-6))
    {
      // Stokes dimension
      const Index stokes_dim = doit_i_field.ncols();
      // Number of species
      const Index N_species = vmr_field.nbooks();
      
      // Ppath_step normally has 2 points, the starting
      // point and the intersection point.
      // But there can be points in between, when a maximum 
      // l_step is given. We have to interpolate on all the 
      // points in the ppath_step.

      // Initialize variables for interpolated values
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np, 0.);
      Matrix abs_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix doit_i_field_int(stokes_dim, ppath_step.np, 0.);
      Vector t_int(ppath_step.np, 0.);
      Matrix vmr_list_int(N_species, ppath_step.np, 0.);
      Vector p_int(ppath_step.np, 0.);
      
      interp_cloud_coeff1D(ext_mat_int, abs_vec_int, sca_vec_int,
                           doit_i_field_int, t_int, vmr_list_int, p_int, 
                           ext_mat_field, abs_vec_field, doit_scat_field, 
                           doit_i_field, t_field, vmr_field, p_grid, 
                           ppath_step, cloudbox_limits, scat_za_grid, 
                           scat_za_interp);
     
          
      // ppath_what_background(ppath_step) tells the 
      // radiative background.  More information in the 
      // function get_radiative_background.
      // if there is no background we proceed the RT
      Index bkgr = ppath_what_background(ppath_step);

      // Radiative transfer from one layer to the next, starting
      // at the intersection with the next layer and propagating
      // to the considered point.
      cloud_RT_no_background(doit_i_field, 
                             abs_scalar_gas_agenda,
                                 opt_prop_gas_agenda, ppath_step, 
                             t_int, vmr_list_int,
                             ext_mat_int, abs_vec_int, sca_vec_int,
                                 doit_i_field_int,
                             p_int, cloudbox_limits, 
                                 f_grid, f_index, p_index, 0, 0,
                             scat_za_index, 0);
      
      // bkgr=2 indicates that the background is surface
      if (bkgr == 2)
        {
          // cout << "hit surface "<< ppath_step.gp_p << endl;
          cloud_RT_surface(
                           doit_i_field, surface_prop_agenda, 
                           f_index, stokes_dim, ppath_step, cloudbox_limits, 
                           scat_za_grid, scat_za_index); 
          
        }
      
    }//end if inside cloudbox
}

//! cloud_ppath_update1D_noseq
/*
  Basically the same as cloud_ppath_update1D, the only difference is that
  i_field_old is always used as incoming Stokes vector.

  \author Claudia Emde
  \date 2005-05-04
*/
void cloud_ppath_update1D_noseq(
                          // Output
                          Tensor6View doit_i_field,
                          // ppath_step_agenda:
                          const Index& p_index,
                          const Index& scat_za_index,
                          ConstVectorView scat_za_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstTensor6View doit_i_field_old,
                          ConstTensor6View doit_scat_field,
                          // Calculate scalar gas absorption:
                          const Agenda& abs_scalar_gas_agenda,
                          ConstTensor4View vmr_field,
                          // Gas absorption: 
                          const Agenda& opt_prop_gas_agenda,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          ConstVectorView  p_grid,
                          ConstTensor3View z_field,
                          ConstMatrixView r_geoid,
                          ConstMatrixView z_surface,
                          // Calculate thermal emission:
                          ConstTensor3View t_field,
                          ConstVectorView f_grid,
                          // used for surface ?
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          //const Agenda& surface_prop_agenda,
                          const Index& scat_za_interp
                         )
{
  Matrix iy;
  Matrix surface_emission;
  Matrix surface_los;
  Tensor4 surface_rmatrix;
  Ppath ppath_step;
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from 
  // where this function is called.
  
  //Initialize ppath for 1D.
  ppath_init_structure(ppath_step, 1, 1);
  // See documentation of ppath_init_structure for understanding
  // the parameters.
  
  // Assign value to ppath.pos:
  ppath_step.z[0]     = z_field(p_index,0,0);
  ppath_step.pos(0,0) = r_geoid(0,0) + ppath_step.z[0];
  
  // Define the direction:
  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
  
  
  // Define the grid positions:
  ppath_step.gp_p[0].idx   = p_index;
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;
  
  // Call ppath_step_agenda: 
  Vector unused_lat_grid(0);
  Vector unused_lon_grid(0);
  ppath_step_agendaExecute(ppath_step, 1, p_grid,
                           unused_lat_grid, unused_lon_grid,
                           z_field, r_geoid, z_surface,
                           ppath_step_agenda, true);
  
  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  
  if((cloudbox_limits[0] <= ppath_step.gp_p[1].idx) &&
     cloudbox_limits[1] > ppath_step.gp_p[1].idx ||
     (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
      abs(ppath_step.gp_p[1].fd[0]) < 1e-6))
    {
      // Stokes dimension
      const Index stokes_dim = doit_i_field.ncols();
      // Number of species
      const Index N_species = vmr_field.nbooks();
      
      // Ppath_step normally has 2 points, the starting
      // point and the intersection point.
      // But there can be points in between, when a maximum 
      // l_step is given. We have to interpolate on all the 
      // points in the ppath_step.

      // Initialize variables for interpolated values
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np, 0.);
      Matrix abs_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix doit_i_field_int(stokes_dim, ppath_step.np, 0.);
      Vector t_int(ppath_step.np, 0.);
      Matrix vmr_list_int(N_species, ppath_step.np, 0.);
      Vector p_int(ppath_step.np, 0.);
      
      interp_cloud_coeff1D(ext_mat_int, abs_vec_int, sca_vec_int,
                           doit_i_field_int, t_int, vmr_list_int, p_int, 
                           ext_mat_field, abs_vec_field, doit_scat_field, 
                           doit_i_field_old, t_field, vmr_field, p_grid, 
                           ppath_step, cloudbox_limits, scat_za_grid, 
                           scat_za_interp);
      
      // ppath_what_background(ppath_step) tells the 
      // radiative background.  More information in the 
      // function get_radiative_background.
      // if there is no background we proceed the RT
      Index bkgr = ppath_what_background(ppath_step);
      
      // if 0, there is no background
      if (bkgr == 0)
        {
          // Radiative transfer from one layer to the next, starting
          // at the intersection with the next layer and propagating
          // to the considered point.
          cloud_RT_no_background(doit_i_field,
                                 abs_scalar_gas_agenda,
                                 opt_prop_gas_agenda, ppath_step, 
                                 t_int, vmr_list_int,
                                 ext_mat_int, abs_vec_int, sca_vec_int,
                                 doit_i_field_int,
                                 p_int, cloudbox_limits, 
                                 f_grid, f_index, p_index, 0, 0, 
                                 scat_za_index, 0);
        }// if loop end - for non_ground background
      
      // bkgr=2 indicates that the background is surface
      else if (bkgr == 2)
        {
          throw runtime_error("Surface reflections only implemented for sequential update.\n" );
        }//end else loop over surface
    }//end if inside cloudbox
}




//! Radiative transfer calculation along a path inside the cloudbox (3D).
/*! 
  This function calculates the radiation field along a propagation path 
  step for a specified zenith direction. This function is used for the 
  sequential update if the radiation field and called inside a loop over
  the pressure grid. 
  In the function the intersection point of the propagation path with the 
  next layer is calculated and all atmospheric properties are 
  interpolated an the intersection point. Then a radiative transfer step is 
  performed using the stokes vector as output and input. The inermediate
  Stokes vectors are stored in the WSV doit_i_field.

 WS Output:
  \param doit_i_field Updated radiation field inside the cloudbox. 
  WS Input:
  \param p_index // Pressure index
  \param lat_index
  \param lon_index
  \param scat_za_index // Index for proagation direction
  \param scat_aa_index
  \param scat_za_grid
  \param scat_aa_grid
  \param cloudbox_limits 
  \param doit_scat_field Scattered field.
  Calculate scalar gas absorption:
  \param abs_scalar_gas_agenda
  \param vmr_field
  Scalar gas absorption:
  \param opt_prop_gas_agenda
  Propagation path calculation:
  \param ppath_step_agenda
  \param p_grid
  \param lat_grid
  \param lon_grid
  \param z_field
  \param r_geoid
  Calculate thermal emission:
  \param t_field
  \param f_grid
  \param f_index
  \param ext_mat_field
  \param abs_vec_field

  \author Claudia Emde
  \date 2003-06-04
*/
void cloud_ppath_update3D(
                          Tensor6View doit_i_field,
                          // ppath_step_agenda:
                          const Index& p_index,
                          const Index& lat_index,
                          const Index& lon_index,
                          const Index& scat_za_index,
                          const Index& scat_aa_index,
                          ConstVectorView scat_za_grid,
                          ConstVectorView scat_aa_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstTensor6View doit_scat_field,
                          // Calculate scalar gas absorption:
                          const Agenda& abs_scalar_gas_agenda,
                          ConstTensor4View vmr_field,
                          // Gas absorption:
                          const Agenda& opt_prop_gas_agenda,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          ConstVectorView p_grid,
                          ConstVectorView lat_grid,
                          ConstVectorView lon_grid,
                          ConstTensor3View z_field,
                          ConstMatrixView r_geoid,
                          ConstMatrixView z_surface,
                          // Calculate thermal emission:
                          ConstTensor3View t_field,
                          ConstVectorView f_grid,
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          const Index& //scat_za_interp
                          )
{
  Ppath ppath_step;
  const Index stokes_dim = doit_i_field.ncols();
  
  Vector sca_vec_av(stokes_dim,0);
  Vector aa_grid(scat_aa_grid.nelem());

  for(Index i = 0; i<scat_aa_grid.nelem(); i++)
    aa_grid[i] = scat_aa_grid[i]-180.;

   //Initialize ppath for 3D.
  ppath_init_structure(ppath_step, 3, 1);
  // See documentation of ppath_init_structure for
  // understanding the parameters.
              
  // Assign value to ppath.pos:
  ppath_step.z[0] = z_field(p_index,lat_index,
                            lon_index);

  // The first dimension of pos are the points in 
  // the propagation path. 
  // Here we initialize the first point.
  // The second is: radius, latitude, longitude

  ppath_step.pos(0,0) =
    r_geoid(lat_index, lon_index) + ppath_step.z[0];
  ppath_step.pos(0,1) = lat_grid[lat_index];
  ppath_step.pos(0,2) = lon_grid[lon_index];
              
  // Define the direction:
  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
  ppath_step.los(0,1) = aa_grid[scat_aa_index];
              
  // Define the grid positions:
  ppath_step.gp_p[0].idx   = p_index;
  ppath_step.gp_p[0].fd[0] = 0.;
  ppath_step.gp_p[0].fd[1] = 1.;

  ppath_step.gp_lat[0].idx   = lat_index;
  ppath_step.gp_lat[0].fd[0] = 0.;
  ppath_step.gp_lat[0].fd[1] = 1.;
                    
  ppath_step.gp_lon[0].idx   = lon_index;
  ppath_step.gp_lon[0].fd[0] = 0.;
  ppath_step.gp_lon[0].fd[1] = 1.;

  // Call ppath_step_agenda: 
  ppath_step_agendaExecute(ppath_step, 3, p_grid,
                           lat_grid, lon_grid, z_field, r_geoid, z_surface,
                           ppath_step_agenda, true);

    // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  if (is_inside_cloudbox(ppath_step, cloudbox_limits, true))
    {      
      // Gridpositions inside the cloudbox.
      // The optical properties are stored only inside the
      // cloudbox. For interpolation we use grids
      // inside the cloudbox.
      
      ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;
      ArrayOfGridPos cloud_gp_lat = ppath_step.gp_lat;
      ArrayOfGridPos cloud_gp_lon = ppath_step.gp_lon;
       
      for(Index i = 0; i<ppath_step.np; i++ )
        {
          cloud_gp_p[i].idx -= cloudbox_limits[0];  
          cloud_gp_lat[i].idx -= cloudbox_limits[2];
          cloud_gp_lon[i].idx -= cloudbox_limits[4];
          
          // This is necessary because the ppath_step_agenda sometimes 
          // returns nan values. (FIXME: Report this problem to Patrick)
        //   const Numeric TOL = 1e-6;
//           if (abs(ppath_step.los(0,0)) < TOL)
//             ppath_step.los(i,0) = 0.;
//           if (abs(ppath_step.los(0,0)-180.) < TOL)
//             ppath_step.los(i,0) = 180.;
          
        }
      
      fix_gridpos_at_boundary(cloud_gp_p, cloudbox_limits[1] - cloudbox_limits[0] +1);
      fix_gridpos_at_boundary(cloud_gp_lat, cloudbox_limits[3] - cloudbox_limits[2] +1); 
      fix_gridpos_at_boundary(cloud_gp_lon, cloudbox_limits[5] - cloudbox_limits[4] +1);
      
      Matrix itw(ppath_step.np, 8);
      interpweights(itw, cloud_gp_p, cloud_gp_lat, cloud_gp_lon);

      Matrix itw_p(ppath_step.np, 2);
      interpweights(itw_p, cloud_gp_p);
 
      // The zenith angles and azimuth of the propagation path are
      // needed as we have to 
      // interpolate the intensity field and the scattered field on the 
      // right angles.
      VectorView los_grid_za = ppath_step.los(joker,0);
      VectorView los_grid_aa = ppath_step.los(joker,1);

      for(Index i = 0; i<los_grid_aa.nelem(); i++)
        los_grid_aa[i] = los_grid_aa[i] + 180.;
  
      ArrayOfGridPos gp_za(los_grid_za.nelem()); 
      gridpos(gp_za, scat_za_grid, los_grid_za);

      ArrayOfGridPos gp_aa(los_grid_aa.nelem()); 
      gridpos(gp_aa, scat_aa_grid, los_grid_aa);

      Matrix itw_p_za(ppath_step.np, 32);
      interpweights(itw_p_za, cloud_gp_p, cloud_gp_lat, cloud_gp_lon, 
                    gp_za, gp_aa);
            
      // Ppath_step normally has 2 points, the starting
      // point and the intersection point.
      // But there can be points in between, when a maximum 
      // l_step is given. We have to interpolate on all the 
      // points in the ppath_step.
      
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np);
      Matrix abs_vec_int(stokes_dim, ppath_step.np);
      Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix doit_i_field_int(stokes_dim, ppath_step.np, 0.);
      Vector t_int(ppath_step.np);
      Vector vmr_int(ppath_step.np);
      Vector p_int(ppath_step.np);
      Vector stokes_vec(stokes_dim);
      //Tensor3 ext_mat_gas(stokes_dim, stokes_dim, ppath_step.np);
      //Matrix abs_vec_gas(stokes_dim, ppath_step.np);
      
            // Calculate the average of the coefficients for the layers
      // to be considered in the 
      // radiative transfer calculation.
      
      for (Index i = 0; i < stokes_dim; i++)
        {
          // Extinction matrix requires a second loop 
          // over stokes_dim
          out3 << "Interpolate ext_mat:\n";
          for (Index j = 0; j < stokes_dim; j++)
            {
              //
              // Interpolation of ext_mat
              //
              interp( ext_mat_int(i, j, joker), itw, 
                      ext_mat_field(joker, joker, joker, i, j), cloud_gp_p,
                      cloud_gp_lat, cloud_gp_lon); 
            }
          // Absorption vector:
          //
          // Interpolation of abs_vec
          //
           interp( abs_vec_int(i,joker), itw, 
                  abs_vec_field(joker, joker, joker, i),
                   cloud_gp_p, cloud_gp_lat, cloud_gp_lon); 
          //
          // Scattered field:
          //
          // Interpolation of sca_vec:
          //
          out3 << "Interpolate doit_scat_field:\n";
          interp( sca_vec_int(i, joker), itw_p_za, 
                  doit_scat_field(joker, joker, joker, joker, joker, i),
                  cloud_gp_p,
                  cloud_gp_lat, cloud_gp_lon, gp_za, gp_aa);
          out3 << "Interpolate doit_i_field:\n";
          interp( doit_i_field_int(i, joker), itw_p_za, 
                  doit_i_field(joker, joker, joker, joker, joker, i), 
                  cloud_gp_p,
                  cloud_gp_lat, cloud_gp_lon, gp_za, gp_aa);
        }
      //
      // Planck function
      // 
      // Interpolate temperature field
      //
      out3 << "Interpolate temperature field\n";
      interp( t_int, itw, 
              t_field(joker, joker, joker), ppath_step.gp_p, 
              ppath_step.gp_lat, ppath_step.gp_lon);
      
      // 
      // The vmr_field is needed for the gaseous absorption 
      // calculation.
      //
      const Index N_species = vmr_field.nbooks();
      //
      // Interpolated vmr_list, holds a vmr_list for each point in 
      // ppath_step.
      //
      Matrix vmr_list_int(N_species, ppath_step.np);
      
      for (Index i = 0; i < N_species; i++)
        {
          out3 << "Interpolate vmr field\n";
          interp( vmr_int, itw, 
                  vmr_field(i, joker, joker, joker), ppath_step.gp_p,
                  ppath_step.gp_lat, ppath_step.gp_lon );
          
          vmr_list_int(i, joker) = vmr_int;
        }
      
      // Presssure (needed for the calculation of gas absorption)
      itw2p( p_int, p_grid, ppath_step.gp_p, itw_p);
      
      out3 << "Calculate radiative transfer inside cloudbox.\n";
      cloud_RT_no_background(doit_i_field, 
                             abs_scalar_gas_agenda,
                             opt_prop_gas_agenda, ppath_step, 
                             t_int, vmr_list_int,
                             ext_mat_int, abs_vec_int, sca_vec_int,
                             doit_i_field_int,
                             p_int, cloudbox_limits, 
                             f_grid, f_index, p_index, lat_index, lon_index, 
                             scat_za_index, scat_aa_index);
    }//end if inside cloudbox
}

//! cloud_RT_no_background
/*
  This function calculates RT in the cloudbox (1D) if the intersection 
  point with the next layer is in the atmosphere (not on the surface). 
  It is used inside the functions cloud_ppath_update1DXXX.
  
  Output: 
  \param doit_i_field Radiation field in cloudbox. This variable is filled
  with the updated values for a given zenith angle (scat_za_index) and 
  pressure (p_index).
  Input:
  \param abs_scalar_gas_agenda Calculate scalar gas absorption. 
  \param opt_prop_gas_agenda Calculate absorption vector and extiction matrix
  due to gas. 
  \param ppath_step Propagation path step from one pressure level to the next 
  (this can include several points)
  \param t_int Temperature values interpolated on propagation path points.
  \param vmr_list_int Interpolated volume mixing ratios. 
  \param ext_mat_int Interpolated particle extinction matrix.
  \param abs_vec_int Interpolated particle absorption vector. 
  \param sca_vec_int Interpolated particle scattering vector. 
  \param doit_i_field_int Interpolated radiances. 
  \param p_int Interpolated pressure values. 
  \param cloudbox_limits Cloudbox_limits. 
  \param f_grid Frequency grid.
  \param f_index Frequency index of (monochromatic) scattering calculation. 
  \param p_index Pressure index in *doit_i_field*.
  \param scat_za_index Zenith angle index in *doit_i_field*.
  
  \author Claudia Emde
  \date 2005-05-13
*/
void cloud_RT_no_background(//Output
                            Tensor6View doit_i_field,
                            // Input
                            const Agenda& abs_scalar_gas_agenda,
                            const Agenda& opt_prop_gas_agenda,
                            const Ppath& ppath_step, 
                            ConstVectorView t_int,
                            ConstMatrixView vmr_list_int,
                            ConstTensor3View ext_mat_int,
                            ConstMatrixView abs_vec_int,
                            ConstMatrixView sca_vec_int,
                            ConstMatrixView doit_i_field_int,
                            ConstVectorView p_int,
                            const ArrayOfIndex& cloudbox_limits,
                            ConstVectorView f_grid,
                            const Index& f_index,
                            const Index& p_index,
                            const Index& lat_index,
                            const Index& lon_index, 
                            const Index& scat_za_index,
                            const Index& scat_aa_index)
{
  const Index N_species = vmr_list_int.nrows();
  const Index stokes_dim = doit_i_field.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem()/2;

  Vector sca_vec_av(stokes_dim,0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rte_vmr_list_local(N_species,0.); 
  
  Matrix abs_scalar_gas_local(1, N_species, 0.);
  Tensor3 ext_mat_local;
  Matrix abs_vec_local;  

  // Incoming stokes vector
  stokes_vec = doit_i_field_int(joker, ppath_step.np-1);

  for( Index k= ppath_step.np-1; k > 0; k--)
    {
      // Length of the path between the two layers.
      Numeric l_step = ppath_step.l_step[k-1];
      // Average temperature
      Numeric rte_temperature_local = 0.5 * (t_int[k] + t_int[k-1]);
      //
      // Average pressure
      Numeric rte_pressure_local = 0.5 * (p_int[k] + p_int[k-1]);
      //
      // Average vmrs
      for (Index i = 0; i < N_species; i++)
        rte_vmr_list_local[i] = 0.5 * (vmr_list_int(i,k) +
                                       vmr_list_int(i,k-1));
      //
      // Calculate scalar gas absorption and add it to abs_vec 
      // and ext_mat.
      //
              
      abs_scalar_gas_agendaExecute(abs_scalar_gas_local, 
                                          f_index, 
                                          rte_pressure_local, 
                                          rte_temperature_local, 
                                          rte_vmr_list_local,
                                          abs_scalar_gas_agenda,
                                          true);
              
      opt_prop_gas_agendaExecute(ext_mat_local, abs_vec_local, 
                                  abs_scalar_gas_local,
                                  opt_prop_gas_agenda,
                                  true);
              
      //
      // Add average particle extinction to ext_mat. 
      //
      for (Index i = 0; i < stokes_dim; i++)
        {
          for (Index j = 0; j < stokes_dim; j++)
            {
              ext_mat_local(0,i,j) += 0.5 *
                (ext_mat_int(i,j,k) + ext_mat_int(i,j,k-1));
            }
          //
          // Add average particle absorption to abs_vec.
          //
          abs_vec_local(0,i) += 0.5 * 
            (abs_vec_int(i,k) + abs_vec_int(i,k-1));
          
          //
          // Averaging of sca_vec:
          //
          sca_vec_av[i] =  0.5 *
            (sca_vec_int(i, k) + sca_vec_int(i, k-1));
                  
        }
      // Frequency
      Numeric f = f_grid[f_index];
      //
      // Calculate Planck function
      //
      Numeric rte_planck_value = planck(f, rte_temperature_local);
              
      // Some messages:
      out3 << "-----------------------------------------\n";
      out3 << "Input for radiative transfer step \n"
           << "calculation inside"
           << " the cloudbox:" << "\n";
      out3 << "Stokes vector at intersection point: \n" 
           << stokes_vec 
           << "\n"; 
      out3 << "l_step: ..." << l_step << "\n";
      out3 << "------------------------------------------\n";
      out3 << "Averaged coefficients: \n";
      out3 << "Planck function: " << rte_planck_value << "\n";
      out3 << "Scattering vector: " << sca_vec_av << "\n"; 
      out3 << "Absorption vector: " << abs_vec_local(0,joker) << "\n"; 
      out3 << "Extinction matrix: " << ext_mat_local(0,joker,joker) << "\n"; 
              
              
      assert (!is_singular( ext_mat_local(0,joker,joker)));
              
      // Radiative transfer step calculation. The Stokes vector
      // is updated until the considered point is reached.
      rte_step_std(stokes_vec, Matrix(stokes_dim,stokes_dim), 
                   ext_mat_local(0,joker,joker), abs_vec_local(0,joker), 
                   sca_vec_av, l_step, rte_planck_value);
              
    }// End of loop over ppath_step. 
  // Assign calculated Stokes Vector to doit_i_field. 
  if (atmosphere_dim == 1)
    doit_i_field(p_index - cloudbox_limits[0], 0, 0, scat_za_index, 0, joker)
      = stokes_vec;
  else if (atmosphere_dim == 3)
    doit_i_field(p_index - cloudbox_limits[0],
                 lat_index - cloudbox_limits[2],
                 lon_index - cloudbox_limits[4],
                 scat_za_index, scat_aa_index,
                 joker) = stokes_vec;

}

//! cloud_RT_surface
/*
  This function calculates RT in the cloudbox if the intersection 
  point with the next layer is the surface. 

  CE (2006-05-29) Included surface_prop_agenda here.  

  \author Claudia Emde
  \date 2005-05-13
*/
void cloud_RT_surface(
                      //Output
                      Tensor6View doit_i_field,
                      //Input
                      const Agenda& surface_prop_agenda,
                      const Index& f_index,
                      const Index& stokes_dim,
                      const Ppath& ppath_step,
                      const ArrayOfIndex& cloudbox_limits, 
                      ConstVectorView scat_za_grid, 
                      const Index& scat_za_index
                     )
{
  chk_not_empty( "surface_prop_agenda", surface_prop_agenda );

  Matrix iy; 
  
  // Local output of surface_prop_agenda.
  Matrix surface_emission;
  Matrix surface_los; 
  Tensor4 surface_rmatrix;


  //Set rte_pos, rte_gp_p and rte_los to match the last point
  //in ppath.
  
  Index np = ppath_step.np;
  
  Vector rte_los; 
  rte_los.resize( ppath_step.los.ncols() );
  rte_los = ppath_step.los(np-1,joker);
  
  GridPos dummy_ppath_step_gp_lat;
  GridPos dummy_ppath_step_gp_lon;
  
  //Execute the surface_prop_agenda which gives the surface 
  //parameters.
  
  surface_prop_agendaExecute(surface_emission, surface_los, 
                             surface_rmatrix, ppath_step.gp_p[np-1],
                             dummy_ppath_step_gp_lat, 
                             dummy_ppath_step_gp_lon, rte_los,
                             surface_prop_agenda, true);
  
  iy = surface_emission;

  Vector rtmp(stokes_dim); // Reflected Stokes vector for 1 frequency
  Index nlos = surface_los.nrows();
  
  if( nlos > 0 )
    {
      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          
          mult( rtmp, 
                surface_rmatrix(ilos,f_index,joker,joker), 
                doit_i_field(cloudbox_limits[0], 0, 0,
                             (scat_za_grid.nelem() -1 - scat_za_index), 0,
                             joker) );
          iy(f_index,joker) += rtmp;
          
          doit_i_field(cloudbox_limits[0], 0, 0,
                       scat_za_index, 0,
                       joker) = iy(0, joker);
        }
    }  
}


/*! Calculates for a given point and a given direction one
  propagation path step.

  This function initializes the ppath structure and 
  executes ppath_step_agenda. Output of the fuinction is 
  a propagation path with two points. The starting point 
  and the next point.

  The function is needed in the sequential update
  (doit_i_fieldUpdateSeq3D).
  
  \param ppath_step Propagation path step
  \param ppath_step_agenda Agenda for calculating propagation paths
  \param p Pressure index 
  \param lat Latitude index
  \param lon Longitude index
  \param z_field Altitude field
  \param r_geoid Geoid
  \param scat_za_grid Zenith angle grid
  \param aa_grid Azimuth angle grid
  \param scat_za_index Zenith angle index
  \param scat_aa_index Azimuth angle index
  \param lat_grid Latitude grid
  \param lon_grid Longitude grid

  \author Claudia Emde
  \date 2003-06-06
*/
void ppath_step_in_cloudbox(//Output:
                            Ppath& ppath_step,
                            //Input:
                            const Agenda& ppath_step_agenda,
                            const Index& p,
                            const Index& lat, 
                            const Index& lon,
                            ConstTensor3View z_field,
                            ConstMatrixView r_geoid,
                            ConstMatrixView z_surface,
                            ConstVectorView scat_za_grid,
                            ConstVectorView aa_grid,
                            const Index& scat_za_index,
                            const Index& scat_aa_index,
                            ConstVectorView p_grid,
                            ConstVectorView lat_grid,
                            ConstVectorView lon_grid)
{
  //Initialize ppath for 3D.
  ppath_init_structure(ppath_step, 3, 1);
  // See documentation of ppath_init_structure for
  // understanding the parameters.
    
  // Assign value to ppath.pos:
  //
  ppath_step.z[0] = z_field(p, lat, lon);
                                  
  // The first dimension of pos are the points in 
  // the propagation path. 
  // Here we initialize the first point.
  // The second is: radius, latitude, longitude

  ppath_step.pos(0,0) = r_geoid(lat, lon) + ppath_step.z[0];
  ppath_step.pos(0,1) = lat_grid[lat];
  ppath_step.pos(0,2) = lon_grid[lon];
              
  // Define the direction:
  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
  ppath_step.los(0,1) = aa_grid[scat_aa_index];
              
  // Define the grid positions:
  ppath_step.gp_p[0].idx   = p;
  ppath_step.gp_p[0].fd[0] = 0.;
  ppath_step.gp_p[0].fd[1] = 1.;

  ppath_step.gp_lat[0].idx   = lat;
  ppath_step.gp_lat[0].fd[0] = 0.;
  ppath_step.gp_lat[0].fd[1] = 1.;
                    
  ppath_step.gp_lon[0].idx   = lon;
  ppath_step.gp_lon[0].fd[0] = 0.;
  ppath_step.gp_lon[0].fd[1] = 1.;
              
  // Call ppath_step_agenda: 
  ppath_step_agendaExecute(ppath_step, 3, p_grid,
                           lat_grid, lon_grid, z_field, r_geoid, z_surface,
                           ppath_step_agenda, true);
}

//! interp_cloud_coeff1D 
/*!  
  Interpolate all inputs of the VRTE on a propagation path step. 
  Used in the WSM cloud_ppath_update1D. 
  
  \author Claudia Emde
  \date 2003-06-06
*/
void interp_cloud_coeff1D(//Output
                          Tensor3View ext_mat_int,
                          MatrixView abs_vec_int,
                          MatrixView sca_vec_int,
                          MatrixView doit_i_field_int,
                          VectorView t_int, 
                          MatrixView vmr_list_int,
                          VectorView p_int,
                          //Input
                          ConstTensor5View ext_mat_field, 
                          ConstTensor4View abs_vec_field,
                          ConstTensor6View doit_scat_field,
                          ConstTensor6View doit_i_field,
                          ConstTensor3View t_field, 
                          ConstTensor4View vmr_field, 
                          ConstVectorView p_grid,
                          const Ppath& ppath_step,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstVectorView scat_za_grid,
                          const Index& scat_za_interp)
{
  // Stokes dimension
  const Index stokes_dim = doit_i_field.ncols();
  
  // Gridpositions inside the cloudbox.
  // The optical properties are stored only inside the
  // cloudbox. For interpolation we use grids
  // inside the cloudbox.
  ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;
  
  for(Index i = 0; i < ppath_step.np; i++ )
    cloud_gp_p[i].idx -= cloudbox_limits[0];
  
  Matrix itw(cloud_gp_p.nelem(),2);
  interpweights( itw, cloud_gp_p );
  
  // The zenith angles of the propagation path are needed as we have to 
  // interpolate the intensity field and the scattered field on the 
  // right angles.
  Vector los_grid = ppath_step.los(joker,0);
  
  ArrayOfGridPos gp_za(los_grid.nelem()); 
  gridpos(gp_za, scat_za_grid, los_grid);
  
  Matrix itw_p_za(cloud_gp_p.nelem(), 4);
  interpweights(itw_p_za, cloud_gp_p, gp_za );
  
  // Calculate the average of the coefficients for the layers
  // to be considered in the 
  // radiative transfer calculation.
      
  for (Index i = 0; i < stokes_dim; i++)
    {
      // Extinction matrix requires a second loop 
      // over stokes_dim
      out3 << "Interpolate ext_mat:\n";
      for (Index j = 0; j < stokes_dim; j++)
        {
          //
          // Interpolation of ext_mat
          //
          interp( ext_mat_int(i, j, joker), itw, 
                  ext_mat_field(joker, 0, 0, i, j), cloud_gp_p ); 
        }
      // Particle absorption vector:
      //
      // Interpolation of abs_vec
      //  //
      out3 << "Interpolate abs_vec:\n";
      interp( abs_vec_int(i,joker), itw, 
              abs_vec_field(joker, 0, 0, i), cloud_gp_p ); 
      //
      // Scattered field:
      //
      //
      
      out3 << "Interpolate doit_scat_field and doit_i_field:\n";
      if(scat_za_interp == 0) // linear interpolation
        {
          interp( sca_vec_int(i, joker), itw_p_za, 
                  doit_scat_field(joker, 0, 0, joker, 0, i), 
                  cloud_gp_p, gp_za);
          interp( doit_i_field_int(i, joker), itw_p_za, 
                  doit_i_field(joker, 0, 0, joker, 0, i), cloud_gp_p, gp_za);
        }
      else if (scat_za_interp == 1) //polynomial interpolation
        {
          // These intermediate variables are needed because polynomial 
          // interpolation is not implemented as multidimensional 
          // interpolation.
          Tensor3 sca_vec_int_za(stokes_dim, ppath_step.np,
                                 scat_za_grid.nelem(), 0.);
          Tensor3 doit_i_field_int_za(stokes_dim, ppath_step.np, 
                                      scat_za_grid.nelem(), 0.);
          for (Index za = 0; za < scat_za_grid.nelem(); za++)
            {
              interp( sca_vec_int_za(i, joker, za), itw, 
                      doit_scat_field(joker, 0, 0, za, 0, i), cloud_gp_p);
              out3 << "Interpolate doit_i_field:\n";
              interp( doit_i_field_int_za(i, joker, za), itw, 
                      doit_i_field(joker, 0, 0, za, 0, i), cloud_gp_p);
            }
          for (Index ip = 0; ip < ppath_step.np; ip ++)
            {
              sca_vec_int(i, ip) = 
                interp_poly(scat_za_grid, sca_vec_int_za(i, ip, joker),
                            los_grid[ip], gp_za[ip]);
              doit_i_field_int(i, ip) = 
                interp_poly(scat_za_grid, doit_i_field_int_za(i, ip, joker),
                            los_grid[ip], gp_za[ip]);
            }
        }
    }
  //
  // Planck function
  // 
  // Interpolate temperature field
  //
  out3 << "Interpolate temperature field\n";
  interp( t_int, itw, t_field(joker, 0, 0), ppath_step.gp_p );
  // 
  // The vmr_field is needed for the gaseous absorption 
  // calculation.
  //
  const Index N_species = vmr_field.nbooks();
  //
  // Interpolated vmr_list, holds a vmr_list for each point in 
  // ppath_step.
  //
  Vector vmr_int(ppath_step.np);
      
  for (Index i_sp = 0; i_sp < N_species; i_sp ++)
    {
      out3 << "Interpolate vmr field\n";
      interp( vmr_int, itw, vmr_field(i_sp, joker, 0, 0), ppath_step.gp_p );
      vmr_list_int(i_sp, joker) = vmr_int;
    }
  // 
  // Interpolate pressure
  //
  itw2p( p_int, p_grid, ppath_step.gp_p, itw);
}

/*! Checks, whether a gridpoint is inside the cloudbox.

    \return true is returned if the point is inside the 
          cloudbox.
          
  \param gp_p  pressure GridPos
  \param gp_lat latitude GridPos
  \param gp_lon longitude GridPos
  \param cloudbox_limits The limits of the cloudbox.
  \param include_boundaries boolean: determines whther or not points on the 
  boundary are considered to be inside the cloudbox.

  \author Claudia Emde (rewritten by Cory Davis 2005-07-03)
  \date 2003-06-06

*/
bool is_gp_inside_cloudbox(const GridPos& gp_p,
                           const GridPos& gp_lat,
                           const GridPos& gp_lon,
                           const ArrayOfIndex& cloudbox_limits,
                           const bool include_boundaries)
                        
{
  bool result=true;
  // Pressure dimension
  double ipos = fractional_gp( gp_p );
  if (include_boundaries){
    if( ipos < double( cloudbox_limits[0] )  ||
        ipos > double( cloudbox_limits[1] ) )
      { result=false; }
    
    else {
      // Latitude dimension
      ipos = fractional_gp( gp_lat );
      if( ipos < double( cloudbox_limits[2] )  || 
          ipos > double( cloudbox_limits[3] ) )
        { result=false; }
      
      else
        {
          // Longitude dimension
          ipos = fractional_gp( gp_lon );
          if( ipos < double( cloudbox_limits[4] )  || 
              ipos > double( cloudbox_limits[5] ) )
            { result=false; } 
        }
    }
  }
  else
    {
      if( ipos <= double( cloudbox_limits[0] )  ||
          ipos >= double( cloudbox_limits[1] ) )
        { result=false; }
      
      else {
        // Latitude dimension
        ipos = fractional_gp( gp_lat );
        if( ipos <= double( cloudbox_limits[2] )  || 
            ipos >= double( cloudbox_limits[3] ) )
          { result=false; }
        
        else
          {
            // Longitude dimension
            ipos = fractional_gp( gp_lon );
            if( ipos <= double( cloudbox_limits[4] )  || 
                ipos >= double( cloudbox_limits[5] ) )
              { result=false; }           
          }
      }
    }
  return result;
  
}


/*! Checks, whether the last point of a propagation path 
  is inside the cloudbox.

    \return true is returned if the point is inside the 
          cloudbox.
          
  \param ppath_step Propagation path step.
  \param cloudbox_limits The limits of the cloudbox.
  \param include_boundaries boolean: determines whther or not points on the 
  boundary are considered to be inside the cloudbox.

  \author Claudia Emde (rewritten by Cory Davis 2005-07-03)
  \date 2003-06-06

*/
bool is_inside_cloudbox(const Ppath& ppath_step,
                        const ArrayOfIndex& cloudbox_limits,
                        const bool include_boundaries)
                        
{
  const Index np=ppath_step.np;
  
  return is_gp_inside_cloudbox(ppath_step.gp_p[np-1],ppath_step.gp_lat[np-1],
                               ppath_step.gp_lon[np-1],cloudbox_limits,include_boundaries);
  
}


//! Radiative transfer calculation inside cloudbox for planeparallel case.
/*! 
  This function calculates the radiation field along a line of sight. This 
  function is used for the sequential update of the radiation field and 
  called inside a loop over the pressure grid. 
  
  The function gets all the atmospheric points on the pressure grid. 
  Then a radiative transfer step is 
  performed using the stokes vector as output and input. The inermediate
  Stokes vectors are stored in the WSV doit_i_field.

 WS Output:
  \param doit_i_field Updated radiation field inside the cloudbox. 
  Variables used in opt_prop_xxx_agenda:
  \param ext_mat
  \param abs_vec  
  Variables used in ppath_step_agenda:
  \param ppath_step
  WS Input:
  \param p_index // Pressure index
  \param scat_za_index // Index for proagation direction
  \param scat_za_grid
  \param cloudbox_limits 
  \param doit_scat_field Scattered field.
  Calculate scalar gas absorption:
  \param abs_scalar_gas_agenda
  \param vmr_field
  Scalar gas absorption:
  \param opt_prop_gas_agenda
  Propagation path calculation:
  \param ppath_step_agenda
  \param p_grid
  \param z_field
  \param r_geoid
  Calculate thermal emission:
  \param t_field
  \param f_grid
  \param f_index
  \param ext_mat_field
  \param abs_vec_field

  \author Sreerekha Ravi
  \date 2003-11-17
*/
void cloud_ppath_update1D_planeparallel(
                                        Tensor6View doit_i_field,
                                        const Index& p_index,
                                        const Index& scat_za_index,
                                        ConstVectorView scat_za_grid,
                                        const ArrayOfIndex& cloudbox_limits,
                                        ConstTensor6View doit_scat_field,
                                        // Calculate scalar gas absorption:
                                        const Agenda& abs_scalar_gas_agenda,
                                        ConstTensor4View vmr_field,
                                        // Gas absorption:
                                        const Agenda& opt_prop_gas_agenda,
                                        // Propagation path calculation:
                                        const Agenda& /* ppath_step_agenda */,
                                        ConstVectorView p_grid,
                                        ConstTensor3View z_field,
                                        ConstMatrixView r_geoid,
                                        // Calculate thermal emission:
                                        ConstTensor3View t_field,
                                        ConstVectorView f_grid,
                                        const Index& f_index,
                                        //particle opticla properties
                                        ConstTensor5View ext_mat_field,
                                        ConstTensor4View abs_vec_field
                                        )
{
  const Index N_species = vmr_field.nbooks();
  const Index stokes_dim = doit_i_field.ncols();
  const Index atmosphere_dim = 1;   
  Matrix abs_scalar_gas(1, N_species, 0.);
  Tensor3 ext_mat;
  Matrix abs_vec;
  Vector rte_vmr_list(N_species,0.); 
  Vector sca_vec_av(stokes_dim,0);
 
  // Radiative transfer from one layer to the next, starting
  // at the intersection with the next layer and propagating
  // to the considered point.
  Vector stokes_vec(stokes_dim);
  Index bkgr;
  if((p_index == 0) && (scat_za_grid[scat_za_index] > 90))
    {
      bkgr = 2;
    }
  else
    {
      bkgr = 0;
    }
      
      // if 0, there is no background
      if (bkgr == 0)
        {
          if(scat_za_grid[scat_za_index] <= 90.0)
            {
              stokes_vec = doit_i_field(p_index-cloudbox_limits[0] +1, 0, 0, scat_za_index, 0, joker);         
              Numeric z_field_above = z_field(p_index +1, 0, 0);
              Numeric z_field_0 = z_field(p_index, 0, 0);
             
              Numeric cos_theta, l_step;
              if (scat_za_grid[scat_za_index] == 90.0)
                {
                  //cos_theta = cos((scat_za_grid[scat_za_index] -1)*PI/180.);
                  //                  cos_theta = cos(89.999999999999*PI/180.);
                  cos_theta = 1e-20;
                 
                }
              else
                {
                  cos_theta = cos(scat_za_grid[scat_za_index]* PI/180.0);
                }
              Numeric dz = (z_field_above -  z_field_0);
              
              l_step =  abs(dz/cos_theta) ;
              
              // Average temperature
              Numeric rte_temperature =   0.5 * (t_field(p_index,0,0)
                                                 + t_field(p_index + 1,0,0));
              //
              // Average pressure
              Numeric rte_pressure = 0.5 * (p_grid[p_index]
                                            + p_grid[p_index + 1]);
              
              // Average vmrs
              for (Index j = 0; j < N_species; j++)
                rte_vmr_list[j] = 0.5 * (vmr_field(j,p_index,0,0) + 
                                         vmr_field(j,p_index + 1,0,0));
              //
              // Calculate scalar gas absorption and add it to abs_vec 
              // and ext_mat.
              //
              
              abs_scalar_gas_agendaExecute(abs_scalar_gas, 
                                                  f_index, 
                                                  rte_pressure, 
                                                  rte_temperature, 
                                                  rte_vmr_list,
                                                  abs_scalar_gas_agenda,
                                                  (p_index != 0));
              
              opt_prop_gas_agendaExecute(ext_mat, abs_vec, abs_scalar_gas,
                                         opt_prop_gas_agenda, (p_index != 0));
              
              //
              // Add average particle extinction to ext_mat. 
        
              for (Index k = 0; k < stokes_dim; k++)
                {
                  for (Index j = 0; j < stokes_dim; j++)
                    {
                      
                      ext_mat(0,k,j) += 0.5 *
                        (ext_mat_field(p_index - cloudbox_limits[0],0,0,k,j) + 
                         ext_mat_field(p_index - cloudbox_limits[0]+ 1,0,0,k,j));
                    }
                  //
                  // Add average particle absorption to abs_vec.
                  //
                  abs_vec(0,k) += 0.5 * 
                    (abs_vec_field(p_index- cloudbox_limits[0],0,0,k) + 
                     abs_vec_field(p_index - cloudbox_limits[0]+ 1,0,0,k));
                  
                  //
                  // Averaging of sca_vec:
                  //
                  sca_vec_av[k] += 0.5 * 
                    (doit_scat_field(p_index- cloudbox_limits[0],0,0,scat_za_index,0, k) + 
                     doit_scat_field(p_index - cloudbox_limits[0]+ 1,0,0,scat_za_index,0,k));
                   
                }
              // Frequency
              Numeric f = f_grid[f_index];
              //
              // Calculate Planck function
              //
              Numeric rte_planck_value = planck(f, rte_temperature);
              
              // Some messages:
              out3 << "-----------------------------------------\n";
              out3 << "Input for radiative transfer step \n"
                   << "calculation inside"
                   << " the cloudbox:" << "\n";
              out3 << "Stokes vector at intersection point: \n" 
                   << stokes_vec 
                   << "\n"; 
              out3 << "l_step: ..." << l_step << "\n";
              out3 << "------------------------------------------\n";
              out3 << "Averaged coefficients: \n";
              out3 << "Planck function: " << rte_planck_value << "\n";
              out3 << "Scattering vector: " << sca_vec_av << "\n"; 
              out3 << "Absorption vector: " << abs_vec(0,joker) << "\n"; 
              out3 << "Extinction matrix: " << ext_mat(0,joker,joker) << "\n"; 
              
              
              assert (!is_singular( ext_mat(0,joker,joker)));
              
              // Radiative transfer step calculation. The Stokes vector
              // is updated until the considered point is reached.
              rte_step_std(stokes_vec, Matrix(stokes_dim,stokes_dim), 
                           ext_mat(0,joker,joker), abs_vec(0,joker), 
                           sca_vec_av, l_step, rte_planck_value);
              
              // Assign calculated Stokes Vector to doit_i_field. 
              doit_i_field(p_index - cloudbox_limits[0],
                      0, 0,
                      scat_za_index, 0,
                      joker) = stokes_vec;
            }
          if(scat_za_grid[scat_za_index] > 90)
            {
              stokes_vec = doit_i_field(p_index-cloudbox_limits[0] -1, 0, 0, scat_za_index, 0, joker);
              Numeric z_field_below = z_field(p_index -1, 0, 0);
              Numeric z_field_0 = z_field(p_index, 0, 0);
              
              Numeric cos_theta, l_step;
              if (scat_za_grid[scat_za_index] == 90.0)
                {
                  cos_theta = cos((scat_za_grid[scat_za_index] -1)*PI/180.);
                  //cos_theta = cos(90.00000001*PI/180.);
                  //cout<<"cos_theta"<<cos_theta<<endl;
                }
              else
                {
                  cos_theta = cos(scat_za_grid[scat_za_index]* PI/180.0);
                }
              Numeric dz = ( z_field_0 - z_field_below);
              l_step =  abs(dz/cos_theta) ;
              
              // Average temperature
              Numeric rte_temperature =   0.5 * (t_field(p_index,0,0)
                                                 + t_field(p_index - 1,0,0));
              //
              // Average pressure
              Numeric rte_pressure = 0.5 * (p_grid[p_index]
                                            + p_grid[p_index - 1]);
              
              //
              // Average vmrs
              for (Index k = 0; k < N_species; k++)
                rte_vmr_list[k] = 0.5 * (vmr_field(k,p_index,0,0) + 
                                         vmr_field(k,p_index - 1,0,0));
              //
              // Calculate scalar gas absorption and add it to abs_vec 
              // and ext_mat.
              //

              abs_scalar_gas_agendaExecute(abs_scalar_gas, 
                                                  f_index, 
                                                  rte_pressure, 
                                                  rte_temperature, 
                                                  rte_vmr_list,
                                                  abs_scalar_gas_agenda,
                                                  (p_index != 0));

              opt_prop_gas_agendaExecute(ext_mat, abs_vec, abs_scalar_gas,
                                         opt_prop_gas_agenda, (p_index != 0));

              //
              // Add average particle extinction to ext_mat. 
              //
             
              // cout<<"cloudbox_limits[1]"<<cloudbox_limits[1]<<endl;
//               cout<<"p_index - cloudbox_limits[0]"<<p_index - cloudbox_limits[0]<<endl;
              for (Index k = 0; k < stokes_dim; k++)
                {
                  for (Index j = 0; j < stokes_dim; j++)
                    {
                      
                      ext_mat(0,k,j) += 0.5 *
                        (ext_mat_field(p_index - cloudbox_limits[0],0,0,k,j) + 
                         ext_mat_field(p_index  - cloudbox_limits[0]- 1,0,0,k,j));
                    }
                  //
                  // Add average particle absorption to abs_vec.
                  //
                  abs_vec(0,k) += 0.5 * 
                    (abs_vec_field(p_index - cloudbox_limits[0],0,0,k) + 
                     abs_vec_field(p_index  - cloudbox_limits[0]- 1,0,0,k));
                  
                  //
                  // Averaging of sca_vec:
                  //
                  sca_vec_av[k] += 0.5 * 
                    (doit_scat_field(p_index- cloudbox_limits[0],0,0,scat_za_index,0, k) + 
                     doit_scat_field(p_index - cloudbox_limits[0]- 1,0,0,scat_za_index,0,k));
                  
                }
              // Frequency
              Numeric f = f_grid[f_index];
              //
              // Calculate Planck function
              //
              Numeric rte_planck_value = planck(f, rte_temperature);
              
              // Some messages:
              out3 << "-----------------------------------------\n";
              out3 << "Input for radiative transfer step \n"
                   << "calculation inside"
                   << " the cloudbox:" << "\n";
              out3 << "Stokes vector at intersection point: \n" 
                   << stokes_vec 
                   << "\n"; 
              out3 << "l_step: ..." << l_step << "\n";
              out3 << "------------------------------------------\n";
              out3 << "Averaged coefficients: \n";
              out3 << "Planck function: " << rte_planck_value << "\n";
              out3 << "Scattering vector: " << sca_vec_av << "\n"; 
              out3 << "Absorption vector: " << abs_vec(0,joker) << "\n"; 
              out3 << "Extinction matrix: " << ext_mat(0,joker,joker) << "\n"; 
              
              
              assert (!is_singular( ext_mat(0,joker,joker)));
              
              // Radiative transfer step calculation. The Stokes vector
              // is updated until the considered point is reached.
              rte_step_std(stokes_vec, Matrix(stokes_dim,stokes_dim), 
                           ext_mat(0,joker,joker), abs_vec(0,joker), 
                           sca_vec_av, l_step, rte_planck_value);
              
              // Assign calculated Stokes Vector to doit_i_field. 
              doit_i_field(p_index - cloudbox_limits[0],
                      0, 0,
                      scat_za_index, 0,
                      joker) = stokes_vec;
            }
          
        
        }// if loop end - for non_ground background
      
      // bkgr=2 indicates that the background is surface
      else if (bkgr == 2)
        {
          //Set rte_pos, rte_gp_p and rte_los to match the last point
          //in ppath.
          //pos
          Vector rte_pos( atmosphere_dim );
          Numeric z_field_0 = z_field(0, 0, 0);
          rte_pos = z_field_0 + r_geoid(0,0);//ppath_step.pos(np-1,Range(0,atmosphere_dim));
          //los
          Vector rte_los(1);
          rte_los = scat_za_grid[scat_za_index];//ppath_step.los(np-1,joker);
          //gp_p
          GridPos rte_gp_p;
          rte_gp_p.idx   = p_index;
          rte_gp_p.fd[0] = 0;
          rte_gp_p.fd[1] = 1;
          //gridpos_copy( rte_gp_p, ppath_step.gp_p[np-1] ); 
          // Executes the surface agenda
          // FIXME: Convert to new agenda scheme before using
          // surface_prop_agenda.execute();

      throw runtime_error( 
                     "Surface reflections inside cloud box not yet handled." );
      /*
        See comment in function above
          // Check returned variables
          if( surface_emission.nrows() != f_grid.nelem()  ||  
              surface_emission.ncols() != stokes_dim )
            throw runtime_error(
                  "The size of the created *surface_emission* is not correct.");

          Index nlos = surface_los.nrows();

          // Define a local vector doit_i_field_sum which adds the 
          // products of groudnd_refl_coeffs with the downwelling 
          // radiation for each elements of surface_los
          Vector doit_i_field_sum(stokes_dim,0);
          // Loop over the surface_los elements
          for( Index ilos=0; ilos < nlos; ilos++ )
            {
              if( stokes_dim == 1 )
                {
                  doit_i_field_sum[0] += surface_refl_coeffs(ilos,f_index,0,0) *
                    doit_i_field(cloudbox_limits[0],
                            0, 0,
                            (scat_za_grid.nelem() -1 - scat_za_index), 0,
                            0);
                }
              else 
                {
                  Vector stokes_vec2(stokes_dim);
                  mult( stokes_vec2, 
                        surface_refl_coeffs(ilos,0,joker,joker), 
                        doit_i_field(cloudbox_limits[0],
                                0, 0,
                                (scat_za_grid.nelem() -1 - scat_za_index), 0,
                                joker));
                  for( Index is=0; is < stokes_dim; is++ )
                    { 
                      doit_i_field_sum[is] += stokes_vec2[is];
                    }
                  
                }
            }
          // Copy from *doit_i_field_sum* to *doit_i_field*, and add the surface emission
          for( Index is=0; is < stokes_dim; is++ )
            {
              doit_i_field (cloudbox_limits[0],
                       0, 0,
                       scat_za_index, 0,
                       is) = doit_i_field_sum[is] + surface_emission(f_index,is);
            }
         
          //cout<<"scat_za_grid"<<scat_za_grid[scat_za_index]<<endl;
          //cout<<"p_index"<<p_index<<endl;
          //cout<<"cloudboxlimit[0]"<<cloudbox_limits[0]<<endl;
          // now the RT is done to the next point in the path.
          // 
          Vector stokes_vec_local;
          stokes_vec_local = doit_i_field (0,
                                      0, 0,
                                      scat_za_index, 0,
                                      joker);
          
          Numeric z_field_above = z_field(p_index +1, 0, 0);
          //Numeric z_field_0 = z_field(p_index, 0, 0);
          Numeric cos_theta;
          if (scat_za_grid[scat_za_index] == 90.0)
            {
              //cos_theta = cos((scat_za_grid[scat_za_index] -1)*PI/180.);
              cos_theta = cos(90.00000001*PI/180.);
            }
          else
            {
              cos_theta = cos(scat_za_grid[scat_za_index]* PI/180.0);
            }
          Numeric dz = (z_field_above -  z_field_0);
          
          Numeric l_step =  abs(dz/cos_theta) ;
          
          // Average temperature
          Numeric rte_temperature =   0.5 * (t_field(p_index,0,0)
                                             + t_field(p_index + 1,0,0));
          
          //
          // Average pressure
          Numeric rte_pressure = 0.5 * (p_grid[p_index]
                                        + p_grid[p_index + 1]);
          
          //
          const Index N_species = vmr_field.nbooks();
          // Average vmrs
          for (Index k = 0; k < N_species; k++)
            {
              rte_vmr_list[k] = 0.5 * (vmr_field(k,p_index,0,0) + 
                                       vmr_field(k,p_index + 1,0,0));
            }
          //
          // Calculate scalar gas absorption and add it to abs_vec 
          // and ext_mat.
          //
          
          // FIXME: Convert to new agenda scheme before using
          //abs_scalar_gas_agenda.execute(p_index);
          
          opt_prop_gas_agendaExecute(ext_mat, abs_vec, abs_scalar_gas,
                                     opt_prop_gas_agenda, (p_index != 0));
          
          //
          // Add average particle extinction to ext_mat. 
          //
          //cout<<"Reached Here ????????????????????????????????????????????????";
          for (Index k = 0; k < stokes_dim; k++)
            {
              for (Index j = 0; j < stokes_dim; j++)
                {
                  ext_mat(0,k,j) += 0.5 *
                    (ext_mat_field(p_index - cloudbox_limits[0],0,0,k,j) + 
                     ext_mat_field(p_index - cloudbox_limits[0]+ 1,0,0,k,j));
                }
              
              
              //
              //
              // Add average particle absorption to abs_vec.
              //
              abs_vec(0,k) += 0.5 * 
                (abs_vec_field(p_index - cloudbox_limits[0],0,0,k) + 
                 abs_vec_field(p_index - cloudbox_limits[0]+1,0,0,k));
              //
              // Averaging of sca_vec:
              //
              sca_vec_av[k] = 0.5*( doit_scat_field(p_index- cloudbox_limits[0], 
                                               0, 0, scat_za_index, 0, k)
                                    + doit_scat_field(p_index- cloudbox_limits[0]+1,
                                                 0, 0, scat_za_index, 0, k)) ;
              
            }
          // Frequency
          Numeric f = f_grid[f_index];
          //
          // Calculate Planck function
          //
          Numeric rte_planck_value = planck(f, rte_temperature);
          
          assert (!is_singular( ext_mat(0,joker,joker)));
          
          // Radiative transfer step calculation. The Stokes vector
          // is updated until the considered point is reached.
          rte_step_std(stokes_vec_local, ext_mat(0,joker,joker), 
                   abs_vec(0,joker), 
                   sca_vec_av, l_step, rte_planck_value);
          // Assign calculated Stokes Vector to doit_i_field. 
          doit_i_field(p_index - cloudbox_limits[0],
                  0, 0,
                  scat_za_index, 0,
                  joker) = stokes_vec_local;
      */  
        }//end else loop over surface
}




/*! Optimize the zenith angle grid, 

  This method optimizes the zenith angle grid. For optimization it uses the 
  interpolation method given by *scat_za_interp* (0 - linear interpolation,
  1 - polynomial interpolation). 
  As input it needs the intensity field calculated on a very fine zenith angle 
  grid (*za_grid_fine*). The function picks out as many grid points as required
  to achieve the required accuracy (*acc* [%]). This methods optimizes only 
  the intensity (first Stokes component) for 1D cases (first latitude and 
  longitude of the intensity field.  
  (Could be modified to optimize all Stokes components at the same time, if we
  don't want to use the clearsky field for grid optimuzation.)

  \param za_grid_opt Optimized zenith angle grid.
  \param doit_i_field_opt Optimized intensity field. 
  \param za_grid_fine Fine zenith angle grid.
  \param doit_i_field Radiation field calculated on a very fine za grid. 
  \param acc Accuracy of optimization [%].
  \param scat_za_interp Interpolation method. 

  \author Claudia Emde
  \date 2004-04-05
*/
void za_gridOpt(//Output:
                Vector& za_grid_opt,
                Matrix& doit_i_field_opt,
                // Input
                ConstVectorView za_grid_fine,
                ConstTensor6View doit_i_field,
                const Numeric& acc,
                const Index& scat_za_interp)
{
  Index N_za = za_grid_fine.nelem();

  assert(doit_i_field.npages() == N_za);
  
  Index N_p = doit_i_field.nvitrines();
  
  Vector i_approx_interp(N_za);
  Vector za_reduced(2);

  ArrayOfIndex idx;
  idx.push_back(0);
  idx.push_back(N_za-1);
  ArrayOfIndex idx_unsorted;

  Numeric max_diff = 100;

  ArrayOfGridPos gp_za(N_za);
  Matrix itw(za_grid_fine.nelem(), 2);

  ArrayOfIndex i_sort;
  Vector diff_vec(N_za);
  Vector max_diff_za(N_p);
  ArrayOfIndex ind_za(N_p);
  Numeric max_diff_p;
  Index ind_p=0;
  
  while( max_diff > acc )
    {
      za_reduced.resize(idx.nelem());
      doit_i_field_opt.resize(N_p, idx.nelem());
      max_diff_za = 0.;
      max_diff_p = 0.;

      // Interpolate reduced internsity field on fine za_grid for 
      // all pressure levels
      for( Index i_p = 0; i_p < N_p; i_p++ )
        {
          for( Index i_za_red = 0; i_za_red < idx.nelem(); i_za_red ++)
            {
              za_reduced[i_za_red] = za_grid_fine[idx[i_za_red]];
              doit_i_field_opt(i_p, i_za_red) = doit_i_field(i_p, 0, 0, idx[i_za_red], 
                                                   0, 0);
            }
          // Calculate grid positions
          gridpos(gp_za, za_reduced, za_grid_fine); 
          //linear interpolation 
          if(scat_za_interp == 0 || idx.nelem() < 3)
            {
              interpweights(itw, gp_za);
              interp(i_approx_interp, itw, doit_i_field_opt(i_p, joker), gp_za);
            }
          else if(scat_za_interp == 1)
            {
              for(Index i_za = 0; i_za < N_za; i_za ++)
                {
                  i_approx_interp[i_za] = 
                    interp_poly(za_reduced, doit_i_field_opt(i_p, joker),
                                 za_grid_fine[i_za],
                                 gp_za[i_za]);
                }
            }
          else
            // Interpolation method not defined
            assert(false);
          
          // Calculate differences between approximated i-vector and 
          // exact i_vector for the i_p pressure level
          for (Index i_za = 0; i_za < N_za; i_za ++)
            {
              diff_vec[i_za]  =  abs( doit_i_field(i_p, 0, 0, i_za, 0 ,0)
                                      -  i_approx_interp[i_za]);
              if( diff_vec[i_za] > max_diff_za[i_p] )
                {
                  max_diff_za[i_p] = diff_vec[i_za];
                  ind_za[i_p] = i_za;
                }
            }
          // Take maximum value of max_diff_za
          if( max_diff_za[i_p] > max_diff_p )
            {
              max_diff_p = max_diff_za[i_p];
              ind_p = i_p;
            }
        }
      
      
      //Transform in %
      max_diff = max_diff_p/doit_i_field(ind_p, 0, 0, ind_za[ind_p], 0, 0)*100.;
      
      idx.push_back(ind_za[ind_p]);
      idx_unsorted = idx;

      i_sort.resize(idx_unsorted.nelem());
      get_sorted_indexes(i_sort, idx_unsorted);
      
      for (Index i = 0; i<idx_unsorted.nelem(); i++)
        idx[i] = idx_unsorted[i_sort[i]];
  
      za_reduced.resize(idx.nelem());
    }
   
  za_grid_opt.resize(idx.nelem());
  doit_i_field_opt.resize(N_p, idx.nelem());
  for(Index i = 0; i<idx.nelem(); i++)
    {
      za_grid_opt[i] = za_grid_fine[idx[i]];
      doit_i_field_opt(joker, i) = doit_i_field(joker, 0, 0, idx[i], 0, 0);
    }
}
          


            
                    

  
  

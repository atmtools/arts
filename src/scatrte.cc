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
  \file   m_scatrte.cc
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




//! Calculation of scattering properties in the cloudbox.
/*! 
  Calculate ext_mat, abs_vec for all points inside the cloudbox for one 
  propagation direction.
  sca_vec can be obtained from the workspace variable scat_field.
  As we need the average for each layer, it makes sense to calculte
  the coefficients once and store them in an array instead of 
  calculating at each point the coefficient of the point above and 
  the point below. 

  \param ext_mat_field extinction matrix field
  \param abs_vec_field absorption vector field
  \param scat_p_index pressure index in cloudbox
  \param scat_lat_index latitude index in cloudbox
  \param scat_lon_index longitude index in cloudbox
  \param scat_za_index Index for propagation direction
  \param spt_calc_agenda Agenda for calculation of single scattering properties
  \param opt_prop_part_agenda Agenda for summing over all hydrometeor species
  \param cloudbox_limits Cloudbox limits.

  \author Claudia Emde
  \date 2002-06-03
*/
void cloud_fieldsCalc(// Output:
                        Tensor5View ext_mat_field,
                        Tensor4View abs_vec_field,
                        Index& scat_p_index,
                        Index& scat_lat_index,
                        Index& scat_lon_index, 
                        Tensor3& ext_mat,
                        Matrix& abs_vec,  
                        // Input:
                        const Index& scat_za_index,
                        const Index& scat_aa_index,
                        const Agenda& spt_calc_agenda,
                        const Agenda& opt_prop_part_agenda,
                        const ArrayOfIndex& cloudbox_limits
                        )
{
  
  //Calculate optical properties for single particle types:
  spt_calc_agenda.execute(scat_za_index || scat_aa_index);
  
  // Calculate ext_mat, abs_vec for all points inside the cloudbox.
  // sca_vec can be obtained from the workspace variable scat_field.
  // As we need the average for each layer, it makes sense to calculte
  // the coefficients once and store them in an array instead of 
  // calculating at each point the coefficient of the point above and 
  // the point below. 
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:               

  // Loop over all positions inside the cloudbox defined by the 
  // cloudbox_limits.
  for(scat_p_index = cloudbox_limits[0]; scat_p_index
        <= cloudbox_limits[1]; scat_p_index ++)
    {
      for(scat_lat_index = cloudbox_limits[2]; scat_lat_index
            <= cloudbox_limits[3]; scat_lat_index ++)
        {
          for(scat_lon_index = cloudbox_limits[4]; scat_lon_index
                <= cloudbox_limits[5]; scat_lon_index ++)
            {
              // Execute agendas silently, only the first call is
              // output on the screen (no other reason for argument 
              // in agenda.execute).
              opt_prop_part_agenda.execute(scat_za_index ||
                                           scat_aa_index ||
                                           scat_p_index-cloudbox_limits[0] ||
                                           scat_lat_index-cloudbox_limits[2] ||
                                           scat_lon_index-cloudbox_limits[4]);
                      
              // Store coefficients in arrays for the whole cloudbox.
              abs_vec_field(scat_p_index - cloudbox_limits[0], 
                            scat_lat_index - cloudbox_limits[2],
                            scat_lon_index- cloudbox_limits[4],
                            joker) = 
                abs_vec(0, joker);
          
              ext_mat_field(scat_p_index - cloudbox_limits[0], 
                            scat_lat_index - cloudbox_limits[2],
                            scat_lon_index- cloudbox_limits[4],
                            joker, joker) = 
                ext_mat(0, joker, joker);
            } 
        }
    }
}
  




//! Radiative transfer calculation along a path inside the cloudbox.
/*! 
  This function calculates the radiation field along a propagation path 
  step for a specified zenith direction. This function is used for the 
  sequential update if the radiation field and called inside a loop over
  the pressure grid. 
  In the function the intersection point of the propagation path with the 
  next layer is calculated and all atmospheric properties are 
  interpolated an the intersection point. Then a radiative transfer step is 
  performed using the stokes vector as output and input. The inermediate
  Stokes vectors are stored in the WSV i_field.

 WS Output:
  \param i_field Updated radiation field inside the cloudbox. 
  Variables used in scalar_gas_abs_agenda:
  \param abs_scalar_gas
  \param a_pressure
  \param a_temperature
  \param a_vmr_list
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
  \param scat_field Scattered field.
  Calculate scalar gas absorption:
  \param scalar_gas_absorption_agenda
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

  \author Claudia Emde
  \date 2002-06-04
*/
void cloud_ppath_update1D(
                  Tensor6& i_field,
                  Vector& stokes_vec,
                   // scalar_gas_abs_agenda:
                  Numeric& a_pressure,
                  Numeric& a_temperature,
                  Vector& a_vmr_list,
                  // opt_prop_xxx_agenda:
                  Tensor3& ext_mat,
                  Matrix& abs_vec,  
                  // ppath_step_agenda:
                  Ppath& ppath_step, 
                  const Index& p_index,
                  const Index& scat_za_index,
                  const Vector& scat_za_grid,
                  const ArrayOfIndex& cloudbox_limits,
                  const Tensor6& scat_field,
                  // Calculate scalar gas absorption:
                  const Agenda& scalar_gas_absorption_agenda,
                  const Tensor4& vmr_field,
                  // Gas absorption:
                  const Agenda& opt_prop_gas_agenda,
                  // Propagation path calculation:
                  const Agenda& ppath_step_agenda,
                  const Vector& p_grid,
                  const Tensor3& z_field,
                  const Matrix& r_geoid,
                  // Calculate thermal emission:
                  const Tensor3& t_field,
                  const Vector& f_grid,
                  const Index& f_index,
                  //particle opticla properties
                  ConstTensor5View ext_mat_field,
                  ConstTensor4View abs_vec_field
                  )
{

  const Index atmosphere_dim = 1;
  const Index stokes_dim = stokes_vec.nelem();

  // Grid ranges inside cloudbox:
  const Range p_range = Range(cloudbox_limits[0],
                              (cloudbox_limits[1] - cloudbox_limits[0]+1) );
  
  Vector sca_vec_av(stokes_dim,0); 

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
  ppath_step_agenda.execute((scat_za_index + 
                             (p_index - cloudbox_limits[0])));
               
  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  if ((cloudbox_limits[0] <= ppath_step.gp_p[1].idx) &&
      cloudbox_limits[1] > ppath_step.gp_p[1].idx ||
      (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
       fabs(ppath_step.gp_p[1].fd[0]) < 1e-6)) 
    {
                  
      // If the intersection points lies exactly on a 
      // upper boundary the gridposition index is 
      // reduced by one and the first interpolation weight 
      // is set to 1.
                        
      for( Index i = 0; i<2; i++)
        { 
          if (cloudbox_limits[1] == ppath_step.gp_p[i].idx &&
              fabs(ppath_step.gp_p[i].fd[0]) < 1e-6)
            {
              ppath_step.gp_p[i].idx -= 1;
              ppath_step.gp_p[i].fd[0] = 1;
              ppath_step.gp_p[i].fd[1] = 0;
            }
        }
                  
                  
                  
      // Check if the agenda has returned ppath.step with 
      // reasonable values. 
      // PpathPrint( ppath_step, "ppath");
                  
      // Gridpositions inside the cloudbox.
      // The optical properties are stored only inside the
      // cloudbox. For interpolation we use grids
      // inside the cloudbox.
      ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;
      ArrayOfGridPos dummy_gp;
      Vector dummy_grid(0);
                  
                                    
      for(Index i = 0; i < ppath_step.np; i++ )
        cloud_gp_p[i].idx -= cloudbox_limits[0];
                  
      Matrix itw_field;
                  
      interp_atmfield_gp2itw(
                             itw_field, atmosphere_dim,
                             p_grid[ Range(p_range)], dummy_grid,
                             dummy_grid,
                             cloud_gp_p, dummy_gp, dummy_gp);
                  
      // Ppath_step has 2 points, the starting
      // point and the intersection point.
      // But there can be points in between, when a maximum 
      // l_step is given. We have to interpolate on all the 
      // points in the ppath_step.
                  
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np);
      Matrix abs_vec_int(stokes_dim, ppath_step.np);
      Matrix sca_vec_int(stokes_dim, ppath_step.np);
      Vector t_int(ppath_step.np);
      Vector vmr_int(ppath_step.np);
      Vector p_int(ppath_step.np);
                  
                  
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
              interp_atmfield_by_itw
                (ext_mat_int(i, j, joker),
                 atmosphere_dim,
                 p_grid[p_range], dummy_grid, dummy_grid,
                 ext_mat_field(joker, joker, joker, i, j),
                 "ext_mat_array",
                 cloud_gp_p, dummy_gp, dummy_gp,
                 itw_field);
            }
          // Particle absorption vector:
          //
          // Interpolation of abs_vec
          //  //
          out3 << "Interpolate abs_vec:\n";
          interp_atmfield_by_itw
            (abs_vec_int(i,joker),
             atmosphere_dim,
             p_grid[p_range], dummy_grid, dummy_grid, 
             abs_vec_field(joker, joker, joker, i),
             "abs_vec_array",
             cloud_gp_p, dummy_gp, dummy_gp,
             itw_field);
          //
          // Scattered field:
          //
          // Interpolation of sca_vec:
          //
          out3 << "Interpolate scat_field:\n";
          interp_atmfield_by_itw
            (sca_vec_int(i, joker),
             atmosphere_dim,
             p_grid[p_range], dummy_grid, dummy_grid,
             scat_field(joker, joker, joker, scat_za_index,
                        0, i),
             "scat_field",
             cloud_gp_p,
             dummy_gp, dummy_gp,
             itw_field);
        }
                    
      //
      // Planck function
      // 
      // Interpolate temperature field
      //
      out3 << "Interpolate temperature field\n";
      interp_atmfield_by_itw
        (t_int,
         atmosphere_dim,
         p_grid, dummy_grid, dummy_grid,
         t_field(joker, joker, joker),
         "t_field",
         ppath_step.gp_p,
         dummy_gp, dummy_gp,
         itw_field);

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
          interp_atmfield_by_itw
            (vmr_int,
             atmosphere_dim,
             p_grid, dummy_grid, dummy_grid,
             vmr_field(i, joker, joker, joker),
             "vmr_field",
             ppath_step.gp_p,
             dummy_gp, dummy_gp,
             itw_field);
                  
          vmr_list_int(i, joker) = vmr_int;
        }
              
      // 
      // Interpolate pressure
      //
      itw2p( p_int, p_grid, ppath_step.gp_p, itw_field); 
              
      // Radiative transfer from one layer to the next, starting
      // at the intersection with the next layer and propagating
      // to the considered point.
              
      for( Index k= ppath_step.np-1; k > 0; k--)
        {
                  
          // Length of the path between the two layers.
          Numeric l_step = ppath_step.l_step[k-1];
          // Average temperature
          a_temperature =   0.5 * (t_int[k] + t_int[k-1]);
          //
          // Average pressure
          a_pressure = 0.5 * (p_int[k] + p_int[k-1]);
          //
          // Average vmrs
          for (Index i = 0; i < N_species; i++)
            a_vmr_list[i] = 0.5 * (vmr_list_int(i,k) + 
                                   vmr_list_int(i,k-1));
          //
          // Calculate scalar gas absorption and add it to abs_vec 
          // and ext_mat.
          //
		    
          scalar_gas_absorption_agenda.execute(p_index);
		      
          opt_prop_gas_agenda.execute(p_index);
          //
          // Add average particle extinction to ext_mat. 
          //
          for (Index i = 0; i < stokes_dim; i++)
            {
              for (Index j = 0; j < stokes_dim; j++)
                {
                  ext_mat(0,i,j) += 0.5 *
                    (ext_mat_int(i,j,k) + ext_mat_int(i,j,k-1));
			  
                }
              //
              // Add average particle absorption to abs_vec.
              //
              abs_vec(0,i) += 0.5 * 
                (abs_vec_int(i,k) + abs_vec_int(i,k-1));
			  
              //
              // Averaging of sca_vec:
              //
              sca_vec_av[i] =  0.5*
                (sca_vec_int(i,k) + sca_vec_int(i,k-1));
              // 
            }
                  
                  
          // Frequency
          Numeric f = f_grid[f_index];
          //
          // Calculate Planck function
          //
          Numeric a_planck_value = planck(f, a_temperature);
                  
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
          out3 << "Planck function: " << a_planck_value << "\n";
          out3 << "Scattering vector: " << sca_vec_av << "\n"; 
          out3 << "Absorption vector: " << abs_vec(0,joker) << "\n"; 
          out3 << "Extinction matrix: " << ext_mat(0,joker,joker) << "\n"; 
                      
                  
          assert (!is_singular( ext_mat(0,joker,joker)));
                  
          // Radiative transfer step calculation. The Stokes vector
          // is updated until the considered point is reached.
          rte_step(stokes_vec, ext_mat(0,joker,joker), 
                   abs_vec(0,joker), 
                   sca_vec_av, l_step, a_planck_value);
                  
        }// End of loop over ppath_step. 
      // Assign calculated Stokes Vector to i_field. 
      i_field(p_index - cloudbox_limits[0],
              0, 0,
              scat_za_index, 0,
              joker) = stokes_vec;
                   

    } //end if inside cloudbox
  // 
  // If the intersection point is outside the cloudbox
  // no radiative transfer step is performed.
  // The value on the cloudbox boundary remains unchanged.
  //
              
}
          
 



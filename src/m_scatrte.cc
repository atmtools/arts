/* Copyright (C) 2002,2003 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                           Sreerekha T.R. <rekha@uni-bremen.de>
                           
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
  \author Sreerekha T.R. <rekha@uni-bremen.de>
  \date   Wed Jun 19 11:03:57 2002
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox.
  
  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <stdexcept>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include "arts.h"
#include "array.h"
#include "auto_md.h"
#include "check_input.h"
#include "matpackVII.h"
#include "logic.h"
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
#include "m_general.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
  
/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


//! Convergence test (maximum absolute difference). 
/*! 
  The function calculates the absolute differences for two successive 
  iteration fields. It picks out the maximum values for each Stokes 
  component separately. The convergence test is fullfilled under the following
  conditions:
  |I(m+1) - I(m)| < epsilon_1                Intensity.
  |Q(m+1) - Q(m)| < epsilon_2                The other Stokes components.
  |U(m+1) - U(m)| < epsilon_3 
  |V(m+1) - V(m)| < epsilon_4
  The limits for convergence have to be set in the controlfile by setting the 
  vector *epsilon* to appropriate values.

  The conditions have to be valid for all positions in the cloudbox and
  for all directions.
  Then convergence_flag is set to 1. 

  WS Output:
  \param convergence_flag Flag for convergence test.
  WS Input:
  \param i_field Radiation field.
  \param i_field_old Old radiation field.
  Keyword : 
  \param epsilon   Limiting values for the convergence test

  \author Claudia Emde
  \date 2002-06-17

*/
void convergence_flagAbs(//WS Output:
                      Index& convergence_flag,
                      // WS Input:
                      const Tensor6& i_field,
                      const Tensor6& i_field_old,
                      // Keyword:
                      const Vector& epsilon)
{
  //Check the input:
  assert( convergence_flag == 0 );

  const Index N_p = i_field.nvitrines();
  const Index N_lat = i_field.nshelves();
  const Index N_lon = i_field.nbooks();
  const Index N_za = i_field.npages();
  const Index N_aa = i_field.nrows();
  const Index stokes_dim = i_field.ncols();
   
  // Check if i_field and i_field_old have the same dimensions:
  assert(is_size( i_field_old, 
                  N_p, N_lat, N_lon, N_za, N_aa, stokes_dim));
  
  // Check keyword "epsilon":
  if ( epsilon.nelem() != stokes_dim )
    throw runtime_error(
                        "You have to specify limiting values for the "
                        "convergence test for each Stokes component "
                        "separately. That means that *epsilon* must "
                        "have *stokes_dim* elements!"
                        );
  
   for (Index p_index = 0; p_index < N_p; p_index++)
    { 
      for (Index lat_index = 0; lat_index < N_lat; lat_index++)
        {
          for (Index lon_index = 0; lon_index <N_lon; lon_index++)
            {
              for (Index scat_za_index = 0; scat_za_index < N_za;
                   scat_za_index++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < N_aa;
                       scat_aa_index++)
                    {
                      for (Index stokes_index = 0; stokes_index <
                             stokes_dim; stokes_index ++) 
                        {
                          Numeric diff =
                            fabs( i_field(p_index, lat_index, lon_index, 
                                          scat_za_index, scat_aa_index, 
                                          stokes_index) -
                                  i_field_old(p_index, lat_index, 
                                              lon_index, scat_za_index,
                                              scat_aa_index, 
                                              stokes_index ));
                          
                          // If the absolute difference of the components
                          // is larger than the pre-defined values, return
                          // to *i_fieldIterarte* and do next iteration
                          
                          if( diff > epsilon[stokes_index])
                            return;
                          
                          
                        }// End loop stokes_dom.
                    }// End loop scat_aa_grid. 
                }// End loop scat_za_grid.
            }// End loop lon_grid. 
        }// End loop lat_grid.
    } // End p_grid.
  
  // Convergence test has been successful, convergence_flag can be set to 1.
  convergence_flag = 1;
}

//! Convergence test in BT(maximum absolute difference in BT(RJ)). 
/*! 
  The function calculates the absolute differences for two successive 
  iteration fields in BT units. This function works exactly similar to
  convergence_flagAbs.  Additionally, we need frequency as input.

  Note that we use Rayleigh Jeans Brightness temperature for epsilon.
  This is because epsilon is a difference of intensity and Planck BT
  is not linear for small radiance values.  For higher stokes components
  also Planck BT cannot be used because of the same reason.
 
  WS Output:
  \param convergence_flag Fag for convergence test.
  WS Input:
  \param i_field Radiation field.
  \param i_field_old Old radiation field.
  \param f_grid frequency grid
  \param f_index frequency index
  Keyword : 
  \param epsilon   Limiting values for the convergence test in 
                   Rayleigh-Jeans brightness temperature unit.

  \author Sreerekha T.R.
  \date 2003-04-01

*/
void convergence_flagAbs_BT(//WS Output:
                            Index& convergence_flag,
                            // WS Input:
                            const Tensor6& i_field,
                            const Tensor6& i_field_old,
                            const Vector& f_grid,
                            const Index& f_index, 
                            // Keyword:
                            const Vector& epsilon)
{
  //Check the input:
  assert( convergence_flag == 0 );

  const Index N_p = i_field.nvitrines();
  const Index N_lat = i_field.nshelves();
  const Index N_lon = i_field.nbooks();
  const Index N_za = i_field.npages();
  const Index N_aa = i_field.nrows();
  const Index stokes_dim = i_field.ncols();
   
  // Check if i_field and i_field_old have the same dimensions:
  assert(is_size( i_field_old, 
                  N_p, N_lat, N_lon, N_za, N_aa, stokes_dim));
  
  // Check keyword "epsilon":
  if ( epsilon.nelem() != stokes_dim )
    throw runtime_error(
                        "You have to specify limiting values for the "
                        "convergence test for each Stokes component "
                        "separately. That means that *epsilon* must "
                        "have *stokes_dim* elements!"
                        );
  
  //3D atmosphere:
  for (Index p_index = 0; p_index < N_p; p_index++)
    { 
      for (Index lat_index = 0; lat_index < N_lat; lat_index++)
        {
          for (Index lon_index = 0; lon_index <N_lon; lon_index++)
            {
              for (Index scat_za_index = 0; scat_za_index < N_za;
                   scat_za_index++)
                {
                  for (Index scat_aa_index = 0; scat_aa_index < N_aa;
                       scat_aa_index++)
                    {
                      for (Index stokes_index = 0; stokes_index <
                             stokes_dim; stokes_index ++) 
                        {
                          Numeric diff =
                            fabs( i_field(p_index, lat_index, lon_index, 
                                          scat_za_index, scat_aa_index, 
                                          stokes_index) -
                                  i_field_old(p_index, lat_index, 
                                              lon_index, scat_za_index,
                                              scat_aa_index, 
                                              stokes_index ));
                          
                          // If the absolute difference of the components
                          // is larger than the pre-defined values, return
                          // to *i_fieldIterarte* and do next iteration
                          Numeric diff_bt = invrayjean(diff, f_grid[f_index]);
                          if( diff_bt > epsilon[stokes_index])
                            return;
                          
                          
                        }// End loop stokes_dom.
                    }// End loop scat_aa_grid. 
                }// End loop scat_za_grid.
            }// End loop lon_grid. 
        }// End loop lat_grid.
    } // End p_grid.
  
  // Convergence test has been successful, convergence_flag can be set to 1.
  convergence_flag = 1;
}



//! Ckecking grid_stepsize.
/*!
  This method calculates the value of the constant gridstepsize of
  the workspacevariables *scat_za_grid* and *scat_aa_grid* and stores it
  in the workspace-Vector *grid_stepsize*

  If one of the grid stepsizes isn't constant, -1 will be stored

  grid_stepsize[0] <-> scat_za_grid
  grid_stepsize[1] <-> scat_aa_grid

  \param grid_stepsize
  \param scat_za_grid
  \param scat_aa_grid

  \author Claas Teichmann
  \date 2003-05-26
*/

void grid_stepsizeCheck(Vector& grid_stepsize,
                      const Vector& scat_za_grid,
                      const Vector& scat_aa_grid)
{
  grid_stepsize.resize(2);
  grid_stepsize[0] = 0;
  grid_stepsize[1] = 0;

  Index n_za = scat_za_grid.nelem();
  Index n_aa = scat_aa_grid.nelem();

  assert(n_za > 1);
  assert(n_aa > 1);

  Numeric diff_za = scat_za_grid[1] - scat_za_grid[0];
  Numeric diff_aa = scat_aa_grid[1] - scat_aa_grid[0];

  for (Index i = 0; i < n_za - 1; i++)
    {
      if (fabs(diff_za - fabs(scat_za_grid[i] - scat_za_grid[i+1])) > 1e-6) grid_stepsize[0] = -1;
    }
  if (grid_stepsize[0] == 0) grid_stepsize[0]=diff_za;

  for (Index i = 0; i < n_aa - 1; i++)
    {
      if (fabs(diff_aa - fabs(scat_aa_grid[i] - scat_aa_grid[i+1])) > 1e-6) grid_stepsize[1] = -1;
    }
  if (grid_stepsize[1] == 0) grid_stepsize[1]=diff_aa;
  return;
}


//! Set grid_stepsize for scattering integral.
/*!
  In this method the grid stepsizes (zenith angle grid and azimuth angle grid)
  for the scattering integral calculation
  can be specified. For the calculation of the scattering integral equidistant
  grids are appropriate as the peak of the phase function can be anywhere, 
  depending on incoming and scattered directions. 

  grid_stepsize[0] <-> scat_za_grid
  grid_stepsize[1] <-> scat_aa_grid

  WS output:
  \param za_grid_size Number of points in zenith angle grid for 
                      scattering integral calculation 
  \param scat_aa_grid Azimuth angle grid.
  Keywords:
  \param N_za_grid  Number of grid points for scattering    
  \param N_aa_grid  integral calculation 

  \author Claudia Emde
  \date 2003-11-17

  */
void grid_sizeSet(// WS Output:
                  Index& za_grid_size,
                  Vector& scat_aa_grid,
                  // Keywords:
                  const Index& N_za_grid,
                  const Index& N_aa_grid)
{

  // Number of zenith angle grid points (only for scattering integral): 
  za_grid_size = N_za_grid; 
    
  // Azimuth angle grid (the same is used for the scattering integral and
  // for the radiative transfer.
  nlinspace(scat_aa_grid, 0, 360, N_aa_grid);
}



//! Iterative solution of the RTE.
/*! 
  A solution for the RTE with scattering is found using an iterative scheme:

  1. Calculate scattering integral using *scat_field_agenda*.
  2. Calculate RT with fixed scattered field using *scat_rte_agenda*.
  3. Convergence test using *convergence_test_agenda*.

  Note: The atmospheric dimensionality *atmosphere_dim* can be either 1 or 3. 
        To these dimensions the method adapts automatically. 
        If *atmosphere_dim* equals 2, it returns an error message, as 2D
        scattering calculations can not be performed.

  WS Input and Output:
  \param i_field Intensity field.
  \param i_field_old Intensity field from previous iteration. 
  WS Output:
  \param convergence_flag Flag for convergence test
  WS Input:
  \param scat_field_agenda Agenda for calculation of scattering integral
  \param scat_rte_agenda Agenda for RT in cloudbox 
  \param convergence_test_agenda Agenda for convergence test.
  
  \author Claudia Emde
  \date 2002-05-29
*/
void
i_fieldIterate(
               // WS Input and Output:
               Tensor6& i_field,
               Tensor6& i_field_old,
               // WS Output
               Index& convergence_flag,
               // WS Input:  
               const Agenda& scat_field_agenda,
               const Agenda& scat_rte_agenda,
               const Agenda& convergence_test_agenda
               )
{
  convergence_flag = 0;
  while(convergence_flag == 0) {
    
    // 1. Copy i_field to i_field_old.
    i_field_old = i_field;
    
    // 2.Calculate scattered field vector for all points in the cloudbox.
    
    // Calculate the scattered field.
    out2 << "Execute scat_field_agenda. \n";
    scat_field_agenda.execute(true);
    
    // Update i_field.
    out2 << "Execute scat_rte_agenda. \n";
    scat_rte_agenda.execute(true);

    //Convergence test.
    out2 << "Execute convergence_test_agenda. \n";
    convergence_test_agenda.execute(true);

  }//end of while loop, convergence is reached.
}


//! 1D RT calculation inside the cloud box.
/*! 
  This function loops over all grid points and all directions and performs 
  the RT calculation with a fixed scattering integral for one frequency 
  of the frequency grid specified by *f_index*. 
  
  The loop over directions is the outermost loop. Here the optical properties
  for single particle types are calculated as they are not depending on the 
  position of the particles. 
  The inner loop is the loop over the positions. Inside this loop the total 
  optical properties including the partile types as well as the gaseous
  species are calculated. Then the radiative transfer equation can be computed.

  It is not recommended to use this method as it can be very slow. The 
  method i_fieldUpdateSeq1D is much more efficient (see AUG).

  WS Output:
  \param i_field Updated radiation field inside the cloudbox. 
  Variables used in scalar_gas_abs_agenda:
  \param rte_pressure
  \param rte_temperature
  \param rte_vmr_list
  Variables used in spt_calc_agenda:
  \param scat_za_index
  Variables used in opt_prop_xxx_agenda:
  \param ext_mat
  \param abs_vec  
  \param scat_p_index
  Variables used in ppath_step_agenda:
  \param ppath_step
  WS Input:
  \param i_field_old Old radiation field.
  \param scat_field Scattered field.
  \param cloudbox_limits 
  Calculate scalar gas absorption:
  \param scalar_gas_absorption_agenda
  \param vmr_field
  Optical properties for single particle type:
  \param spt_calc_agenda
  \param scat_za_grid
  Optical properties for gases and particles:
  \param opt_prop_part_agenda
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
  \date 2002-05-30
*/
void
i_fieldUpdate1D(// WS Output:
                Tensor6& i_field,
                // scalar_gas_abs_agenda:
                Numeric& rte_pressure,
                Numeric& rte_temperature,
                Vector& rte_vmr_list,
                // spt_calc_agenda:
                Index& scat_za_index ,
                // opt_prop_xxx_agenda:
                Tensor3& ext_mat,
                Matrix& abs_vec,  
                Index& scat_p_index,
                // ppath_step_agenda:
                Ppath& ppath_step, 
                // WS Input:
                const Tensor6& i_field_old,
                const Tensor6& scat_field,
                const ArrayOfIndex& cloudbox_limits,
                // Calculate scalar gas absorption:
                const Agenda& scalar_gas_absorption_agenda,
                const Tensor4& vmr_field,
                // Optical properties for single particle type:
                const Agenda& spt_calc_agenda,
                const Vector& scat_za_grid,
                // Optical properties for gases and particles:
                const Agenda& opt_prop_part_agenda,
                const Agenda& opt_prop_gas_agenda,
                // Propagation path calculation:
                const Agenda& ppath_step_agenda,
                const Vector& p_grid,
                const Tensor3& z_field,
                const Matrix& r_geoid,
                 // Calculate thermal emission:
                const Tensor3& t_field,
                const Vector& f_grid,
                const Index& f_index
                )
{

  out2 << "i_fieldUpdate1D: Radiative transfer calculatiuon in cloudbox. \n";
  out2 << "------------------------------------------------------------- \n";
  
  const Index stokes_dim = scat_field.ncols();
  const Index atmosphere_dim = 1;

  //Check the input
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));

  assert( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));  
  
  // Is the frequency index valid?
  assert( f_index <= f_grid.nelem() );

  // End of checks



  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();

  // Grid ranges inside cloudbox:
  const Range p_range = Range(cloudbox_limits[0],
                              (cloudbox_limits[1] - cloudbox_limits[0]+1) );
  
  //=======================================================================
  // Calculate scattering coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are 
  // calulated for interpolated VMRs, temperature and pressure.
       
  //Loop over all directions, defined by scat_za_grid 
  for(scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      //Calculate optical properties for single particle types:
      spt_calc_agenda.execute(scat_za_index);
    

      // Calculate ext_mat, abs_vec for all points inside the cloudbox.
      // sca_vec can be obtained from the workspace variable scat_field.
      // As we need the average for each layer, it makes sense to calculte
      // the coefficients once and store them in an array instead of 
      // calculating at each point the coefficient of the point above and 
      // the point below. 
      // To use special interpolation functions for atmospheric fields we 
      // use ext_mat_field and abs_vec_field:
      Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                            stokes_dim, stokes_dim, 0.);
      Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                            stokes_dim, 0.);

      // Loop over all positions inside the cloudbox defined by the 
      // cloudbox_limits.
      for(Index p_index = cloudbox_limits[0]; p_index
        <= cloudbox_limits[1]; p_index ++)
        {
          out3 << "-------------------------------------------------\n";
          out3 << " Cloudbox limits: "<< cloudbox_limits[0] << " " 
               << cloudbox_limits[1] << "\n";
          out3 << " Zenith angle index: "<< scat_za_index << "\n";
          out3 << " Pressure index: "<< p_index << "\n";
          
          // The required workspace variable is scat_p_index.
          scat_p_index = p_index; 
          
          // Execute agendas silently, only the first call is output on
          // the screen (no other reason for argument in agenda.execute).
          //  opt_prop_gas_agenda.execute((scat_za_index + 
          //                          (p_index - cloudbox_limits[0])) );
          opt_prop_part_agenda.execute((scat_za_index + 
                                       (p_index - cloudbox_limits[0])) );
          
          // Store coefficients in arrays for the whole cloudbox.
          abs_vec_field(p_index-cloudbox_limits[0], 0, 0, joker) = 
            abs_vec(0, joker);

          ext_mat_field(p_index-cloudbox_limits[0], 0, 0, joker, joker) = 
            ext_mat(0, joker, joker);
                    
        }//End of p_grid loop over the cloudbox
  
      

      //======================================================================
      // Radiative transfer inside the cloudbox
      //=====================================================================
      
      // Define variables which hold averaged coefficients:
      
      Vector stokes_vec(stokes_dim,0.);
      Vector sca_vec_av(stokes_dim,0.);

      // Loop over all positions inside the cloud box defined by the 
      // cloudbox_limits exculding the upper boundary.
      for(Index p_index = cloudbox_limits[0]; p_index
            <= cloudbox_limits[1]; p_index ++)
        {
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

              
              
              // Check if the agenda has returned ppath.step with reasonable 
              // values. 
              //Print( ppath_step, "ppath", 0);
              
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
              Matrix sto_vec_int(stokes_dim, ppath_step.np);
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
                  //
                  // Stokes vector:
                  //
                  // Interpolation of i_field_old.
                  out3 << "Interpolate i_field_old:\n";
                  interp_atmfield_by_itw
                    (sto_vec_int(i, joker),
                     atmosphere_dim,
                     p_grid[p_range], dummy_grid, dummy_grid,
                     i_field_old(joker, joker, joker, scat_za_index,
                             0, i),
                     "i_field_old",
                     cloud_gp_p,
                     dummy_gp, dummy_gp,
                     itw_field);
                  // 
                  // For the radiative transfer equation we 
                  // need the Stokes vector at the intersection
                  // point with the next layer.
                  //
                  stokes_vec[i] = sto_vec_int(i,ppath_step.np-1);
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
                  rte_temperature =   0.5 * (t_int[k] + t_int[k-1]);
                  //
                  // Average pressure
                  rte_pressure = 0.5 * (p_int[k] + p_int[k-1]);
                  //
                  // Average vmrs
                  for (Index i = 0; i < N_species; i++)
                    rte_vmr_list[i] = 0.5 * (vmr_list_int(i,k) + 
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
                  rte_step(stokes_vec, ext_mat(0,joker,joker), 
                           abs_vec(0,joker), 
                           sca_vec_av, l_step, rte_planck_value);
                  
                } // End of loop over ppath_step. 
              // Assign calculated Stokes Vector to i_field. 
              i_field(p_index - cloudbox_limits[0],
                      0, 0,
                      scat_za_index, 0,
                      joker) = stokes_vec;
            } //end if
              // 
              // If the intersection point is outside the cloudbox
              // no radiative transfer step is performed.
              // The value on the cloudbox boundary remains unchanged.
              //
              
        }// Close loop over p_grid (inside cloudbox).
    }// Closes loop over scat_za_grid.
  
} // End of the function.


//! 1D RT calculation inside the cloud box.
/*! 
  This function loops over all grid points and all directions and performs 
  the RT calculation with a fixed scattering integral for one frequency 
  of the frequency grid specified by *f_index*. 
  
  The loop over directions is the outermost loop. Here the optical properties
  for single particle types are calculated as they are not depending on the 
  position of the particles. 
  The inner loop is the loop over the positions. Inside this loop the total 
  optical properties including the partile types as well as the gaseous
  species are calculated. Then the radiative transfer equation can be computed.
  The looping is arranged in such a way, that the radiation field is updated 
  sequentially.
  
  It is recommended to use the sequential update as it is much more efficient
  (see AUG) than the "normal" update. 

  WS Output:
  \param i_field Updated radiation field inside the cloudbox. 
  Variables used in scalar_gas_abs_agenda:
  \param rte_pressure
  \param rte_temperature
  \param rte_vmr_list
  Variables used in spt_calc_agenda:
  \param scat_za_index
  Variables used in opt_prop_xxx_agenda:
  \param ext_mat
  \param abs_vec  
  \param scat_p_index
  Variables used in ppath_step_agenda:
  \param ppath_step
  Ground related variables:
  \param ground_los
  \param ground_emission
  \param ground_refl_coeffs
  \param rte_los
  \param rte_pos
  \param rte_gp_p
  WS Input:
  \param scat_field Scattered field.
  \param cloudbox_limits 
  Calculate scalar gas absorption:
  \param scalar_gas_absorption_agenda
  \param vmr_field
  Optical properties for single particle type:
  \param spt_calc_agenda
  \param scat_za_grid
  Optical properties for gases and particles:
  \param opt_prop_part_agenda
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
  Ground reflection
  \param ground_refl_agenda

  \author Claudia Emde
  \date 2002-05-30
*/
void
i_fieldUpdateSeq1D(// WS Output:
                   Tensor6& i_field,
                   // scalar_gas_abs_agenda:
                   Numeric& rte_pressure,
                   Numeric& rte_temperature,
                   Vector& rte_vmr_list,
                   // spt_calc_agenda:
                   Index& scat_za_index ,
                   // opt_prop_xxx_agenda:
                   Tensor3& ext_mat,
                   Matrix& abs_vec,  
                   Index& scat_p_index,
                   // ppath_step_agenda:
                   Ppath& ppath_step, 
                   // ground related variables STR
                   Matrix& ground_los,
                   Matrix& ground_emission,
                   Tensor4& ground_refl_coeffs,
                   Vector& rte_los,
                   Vector& rte_pos,
                   GridPos& rte_gp_p,
                   // WS Input:
                   const Tensor6& scat_field,
                   const ArrayOfIndex& cloudbox_limits,
                   // Calculate scalar gas absorption:
                   const Agenda& scalar_gas_absorption_agenda,
                   const Tensor4& vmr_field,
                   // Optical properties for single particle type:
                   const Agenda& spt_calc_agenda,
                   const Vector& scat_za_grid,
                   // Optical properties for gases and particles:
                   const Agenda& opt_prop_part_agenda,
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
                   const Agenda& ground_refl_agenda //STR
                   )
{
  
  out2 << "i_fieldUpdateSeq1D: Radiative transfer calculatiuon in cloudbox.\n";
  out2 << "------------------------------------------------------------- \n";
  
  const Index stokes_dim = scat_field.ncols();
 
  //Check the input
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));

  assert( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));  
  
  // Is the frequency index valid?
  assert( f_index <= f_grid.nelem() );

  // End of checks



  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();

 
  
  //=======================================================================
  // Calculate scattering coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are 
  // calulated for interpolated VMRs, temperature and pressure.
      
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, stokes_dim, 0.);
  Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, 0.);
 
  //Loop over all directions, defined by scat_za_grid 
  for(scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      
      //Only dummy variables:
      Index scat_lat_index = 0;
      Index scat_lon_index = 0;

      // This function has to be called inside the angular loop, as
      // it spt_calc_agenda takes *scat_za_index* and *scat_aa_index* 
      // from the workspace.
      cloud_fieldsCalc(ext_mat_field, abs_vec_field, scat_p_index,
                       scat_lat_index, scat_lon_index,
                       ext_mat, abs_vec, 
                       spt_calc_agenda, 
                       opt_prop_part_agenda, cloudbox_limits);
      
      //======================================================================
      // Radiative transfer inside the cloudbox
      //=====================================================================
      
      // Calculate the tangent point for g 

      // Lower boundary of cloudbox
      Index p_low = cloudbox_limits[0];
      //
      // Upper boundary of cloudbox
      Index p_up = cloudbox_limits[1];
      //
      Vector stokes_vec(stokes_dim,0.);
      
      Numeric theta_lim = 180. - asin((r_geoid(0,0)+z_field(p_low,0,0))/
                                     (r_geoid(0,0)+z_field(p_up,0,0)))*RAD2DEG;

      // Sequential update for uplooking angles
      if ( scat_za_grid[scat_za_index] <= 90.) 
        {
          // Loop over all positions inside the cloud box defined by the 
          // cloudbox_limits exculding the upper boundary. For uplooking
          // directions, we start from cloudbox_limits[1]-1 and go down
          // to cloudbox_limits[0] to do a sequential update of the
          // aradiation field
          for(Index p_index = cloudbox_limits[1]-1; p_index
                >= cloudbox_limits[0]; p_index --)
            {
              cloud_ppath_update1D(i_field, 
                                   rte_pressure, rte_temperature, rte_vmr_list,
                                   ext_mat, abs_vec, ground_los,
                                   ground_emission, ground_refl_coeffs,
                                   rte_los, rte_pos, rte_gp_p, ppath_step, 
                                   p_index, scat_za_index, scat_za_grid,
                                   cloudbox_limits, scat_field,
                                   scalar_gas_absorption_agenda, vmr_field,
                                   opt_prop_gas_agenda, ppath_step_agenda,
                                   p_grid,  z_field, r_geoid, t_field, 
                                   f_grid, f_index, ext_mat_field,
                                   abs_vec_field,ground_refl_agenda); 
            }
        }
      else if ( scat_za_grid[scat_za_index] > theta_lim) 
        {
          //
          // Sequential updating for downlooking angles
          //
          for(Index p_index = cloudbox_limits[0]+1; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
              cloud_ppath_update1D(i_field,  
                                   rte_pressure, rte_temperature, rte_vmr_list,
                                   ext_mat, abs_vec, ground_los,
                                    ground_emission, ground_refl_coeffs,
                                    rte_los, rte_pos, rte_gp_p, ppath_step, 
                                   p_index, scat_za_index, scat_za_grid,
                                   cloudbox_limits, scat_field,
                                   scalar_gas_absorption_agenda, vmr_field,
                                   opt_prop_gas_agenda, ppath_step_agenda,
                                   p_grid,  z_field, r_geoid, t_field, 
                                   f_grid, f_index, ext_mat_field, 
                                   abs_vec_field,ground_refl_agenda); 
            }// Close loop over p_grid (inside cloudbox).
        } // end if downlooking.
      
      //
      // Limb looking:
      // We have to include a special case here, as we may miss the endpoints
      // when the intersection point is at the same level as the aactual point.
      // To be save we loop over the full cloudbox. Inside the function 
      // cloud_ppath_update1D it is checked whether the intersection point is 
      // inside the cloudbox or not.
      else if (  scat_za_grid[scat_za_index] > 90 &&
                 scat_za_grid[scat_za_index] < theta_lim ) 
        {
          for(Index p_index = cloudbox_limits[0]; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
              // For this case the cloudbox goes down to the ground an we
              // look downwards. These cases are outside the cloudbox and 
              // not needed. Switch is included here, as ppath_step_agenda 
              // gives an error for such cases.
              if (!(p_index == 0 && scat_za_grid[scat_za_index] > 90.))
                {
                  cloud_ppath_update1D(i_field,  
                                       rte_pressure, rte_temperature, 
                                       rte_vmr_list,
                                       ext_mat, abs_vec, ground_los,
                                       ground_emission, ground_refl_coeffs,
                                       rte_los, rte_pos, rte_gp_p, ppath_step, 
                                       p_index, scat_za_index, scat_za_grid,
                                       cloudbox_limits, scat_field,
                                       scalar_gas_absorption_agenda, vmr_field,
                                       opt_prop_gas_agenda, ppath_step_agenda,
                                       p_grid,  z_field, r_geoid, t_field, 
                                       f_grid, f_index, ext_mat_field, 
                                       abs_vec_field,ground_refl_agenda);
                }
            }
        } 
    }// Closes loop over scat_za_grid.
} // End of the function.




//! 3D RT calculation inside the cloud box.
/*! 
  This function loops over all grid points and all directions and performs 
  the RT calculation with a fixed scattering integral for one frequency 
  of the frequency grid specified by *f_index*. 
  
  The loop over directions is the outermost loop. Here the optical properties
  for single particle types are calculated as they are not depending on the 
  position of the particles. 
  The inner loop is the loop over the positions. Inside this loop the total 
  optical properties including the partile types as well as the gaseous
  species are calculated. Then the radiative transfer equation can be computed.

  It is recommended to use i_fieldUpdateSeq3D as it is much more 
  efficient (see AUG).

  Note: Ground reflection is not yet implemented in 3D scattering 
  calculations.

  WS Output:
  \param i_field Updated radiation field inside the cloudbox. 
  Variables used in scalar_gas_abs_agenda:
  \param rte_pressure
  \param rte_temperature
  \param rte_vmr_list
  Variables used in spt_calc_agenda:
  \param scat_za_index
  \param scat_aa_index
  Variables used in opt_prop_xxx_agenda:
  \param ext_mat
  \param abs_vec  
  \param scat_p_index
  \param scat_lat_index
  \param scat_lon_index
  Variables used in ppath_step_agenda:
  \param ppath_step
  WS Input:
  \param i_field_old Old radiation field.
  \param scat_field Scattered field.
  \param cloudbox_limits 
  Calculate scalar gas absorption:
  \param scalar_gas_absorption_agenda
  \param vmr_field
  Optical properties for single particle type:
  \param spt_calc_agenda
  \param scat_za_grid
  \param scat_aa_grid
  Optical properties for gases and particles:
  \param opt_prop_part_agenda
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

  \author Claudia Emde
  \date 2003-01-07
*/
void i_fieldUpdate3D(// WS Output:
                     Tensor6& i_field,
                     // scalar_gas_abs_agenda:
                     Numeric& rte_pressure,
                     Numeric& rte_temperature,
                     Vector& rte_vmr_list,
                     // spt_calc_agenda:
                     Index& scat_za_index,
                     Index& scat_aa_index,
                     // opt_prop_xxx_agenda:
                     Tensor3& ext_mat,
                     Matrix& abs_vec,  
                     Index& scat_p_index,
                     Index& scat_lat_index,
                     Index& scat_lon_index,
                     // ppath_step_agenda:
                     Ppath& ppath_step, 
                     // WS Input:
                     const Tensor6& i_field_old,
                     const Tensor6& scat_field,
                     const ArrayOfIndex& cloudbox_limits,
                     // Calculate scalar gas absorption:
                     const Agenda& scalar_gas_absorption_agenda,
                     const Tensor4& vmr_field,
                     // Optical properties for single particle type:
                     const Agenda& spt_calc_agenda,
                     const Vector& scat_za_grid,
                     const Vector& scat_aa_grid,
                     // Optical properties for gases and particles:
                     const Agenda& opt_prop_part_agenda,
                     const Agenda& opt_prop_gas_agenda,
                     // Propagation path calculation:
                     const Agenda& ppath_step_agenda,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& z_field,
                     const Matrix& r_geoid,
                     // Calculate thermal emission:
                     const Tensor3& t_field,
                     const Vector& f_grid,
                     const Index& f_index
                     )
{
  
  out2 << "i_fieldUpdate3D: Radiative transfer calculatiuon in cloudbox. \n";
  out2 << "------------------------------------------------------------- \n";

  const Index stokes_dim = scat_field.ncols();
  const Index atmosphere_dim = 3;

  //Check the input
  
  assert( atmosphere_dim == 3);
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) + 1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
                       scat_za_grid.nelem(), 
                       scat_aa_grid.nelem(),
                       stokes_dim));

  assert( is_size( scat_field, 
                     (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                     (cloudbox_limits[3] - cloudbox_limits[2]) + 1, 
                     (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
   
  // Is the frequency index valid?
  assert( f_index <= f_grid.nelem() );

  // End of checks



  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  // Number of azimuth angles
  const Index N_scat_aa = scat_aa_grid.nelem();

  // Grid ranges inside cloudbox:
  const Range p_range = Range(cloudbox_limits[0],
                              (cloudbox_limits[1] - cloudbox_limits[0]+1) );
  const Range lat_range = Range(cloudbox_limits[2],
                                (cloudbox_limits[3] - cloudbox_limits[2]+1) );
  const Range lon_range = Range(cloudbox_limits[4],
                                (cloudbox_limits[5] - cloudbox_limits[4]+1) );
  
  // The definition of the azimth angle grids is different for clearsky and
  // cloudbox. (SHOULD BE FIXED!!!!)
  Vector aa_grid(scat_aa_grid.nelem());
  for(Index i = 0; i<scat_aa_grid.nelem(); i++)
    aa_grid[i] = scat_aa_grid[i] - 180;



  //=======================================================================
  // Calculate coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are 
  // calulated for interpolated VMRs, temperature and pressure.
  
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, stokes_dim, 0.);
  Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, 0.);

  //Loop over all directions, defined by scat_za_grid 
  for(scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
     //Loop over azimuth directions (scat_aa_grid)
      for(scat_aa_index = 0; scat_aa_index < N_scat_aa; scat_aa_index ++)
        {
          cloud_fieldsCalc(ext_mat_field, abs_vec_field, scat_p_index,
                           scat_lat_index, scat_lon_index,
                           ext_mat, abs_vec, 
                           spt_calc_agenda, 
                           opt_prop_part_agenda, cloudbox_limits);
        
          //==================================================================
          // Radiative transfer inside the cloudbox
          //==================================================================

          // Define variables which hold averaged coefficients:

          Vector stokes_vec(stokes_dim,0.);
          Vector sca_vec_av(stokes_dim,0.);
          
          // Loop over all positions inside the cloud box defined by the 
          // cloudbox_limits.
          for(Index p_index = cloudbox_limits[0]; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
              for(Index lat_index = cloudbox_limits[2]; lat_index
                    <= cloudbox_limits[3]; lat_index ++)
                {
                  for(Index lon_index = cloudbox_limits[4]; lon_index
                        <= cloudbox_limits[5]; lon_index ++)
                  {
                    
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
                    ppath_step_agenda.execute(scat_za_index &&
                                              scat_aa_index &&
                                              p_index - cloudbox_limits[0] &&
                                              lat_index - cloudbox_limits[2]&&
                                              lon_index - cloudbox_limits[4]);
              
                    
                    // Check if the agenda has returned ppath.step with 
                    // reasonable values. 
                    // Print( ppath_step, "ppath", 0);
              
                    // Check whether the next point is inside or outside the
                    // cloudbox. Only if the next point lies inside the
                    // cloudbox a radiative transfer step caclulation has to
                    // be performed.
                    
                    // Tolerance value for checking if a point is exactly on
                    // a grid point.
                    const Numeric TOL = 1e-2;
                    
                    if (
                        // inside pressure boundaries
                        (cloudbox_limits[0] <= ppath_step.gp_p[1].idx) &&
                        (cloudbox_limits[1] > ppath_step.gp_p[1].idx ||
                         (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
                          fabs(ppath_step.gp_p[1].fd[0]) < TOL)) &&
                        // inside latitude boundaries 
                        (cloudbox_limits[2] <= ppath_step.gp_lat[1].idx) &&
                        (cloudbox_limits[3] > ppath_step.gp_lat[1].idx ||
                         (cloudbox_limits[3] == ppath_step.gp_lat[1].idx &&
                          fabs(ppath_step.gp_lat[1].fd[0]) < TOL)) &&
                        // inside longitude boundaries 
                        (cloudbox_limits[4] <= ppath_step.gp_lon[1].idx) &&
                         (cloudbox_limits[5] > ppath_step.gp_lon[1].idx ||
                         (cloudbox_limits[5] == ppath_step.gp_lon[1].idx &&
                          fabs(ppath_step.gp_lon[1].fd[0]) < TOL )) 
                        )
                      {

                        // If the intersection points lies exactly on a 
                        // upper boundary the gridposition index is 
                        // reduced by one and the first interpolation weight 
                        // is set to 1.
                        
                        for( Index i = 0; i<2; i++)
                          { 
                            if (cloudbox_limits[1] == ppath_step.gp_p[i].idx &&
                                fabs(ppath_step.gp_p[i].fd[0]) < TOL)
                              {
                                ppath_step.gp_p[i].idx -= 1;
                                ppath_step.gp_p[i].fd[0] = 1;
                                ppath_step.gp_p[i].fd[1] = 0;
                              }
                            
                            if (cloudbox_limits[3]==ppath_step.gp_lat[i].idx &&
                                fabs(ppath_step.gp_lat[i].fd[0]) < TOL)
                              {
                                ppath_step.gp_lat[i].idx -= 1;
                                ppath_step.gp_lat[i].fd[0] = 1;
                                ppath_step.gp_lat[i].fd[1] = 0;
                              }
                            
                            if (cloudbox_limits[5]==ppath_step.gp_lon[i].idx &&
                                fabs(ppath_step.gp_lon[i].fd[0]) < TOL)
                              {
                                ppath_step.gp_lon[i].idx -= 1;
                                ppath_step.gp_lon[i].fd[0] = 1;
                                ppath_step.gp_lon[i].fd[1] = 0;
                              }
                          }
                        
                        // Gridpositions inside the cloudbox.
                        // The optical properties are stored only inside the
                        // cloudbox. For interpolation we use grids
                        // inside the cloudbox.
                        
                        ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;
                        ArrayOfGridPos cloud_gp_lat = ppath_step.gp_lat;
                        ArrayOfGridPos cloud_gp_lon = ppath_step.gp_lon;
                        
                        for(Index i = 0; i<2; i++ )
                          {
                           cloud_gp_p[i].idx -= cloudbox_limits[0];  
                           cloud_gp_lat[i].idx -= cloudbox_limits[2];
                           cloud_gp_lon[i].idx -= cloudbox_limits[4];
                          }

                        Matrix itw_field;
                        
                        interp_atmfield_gp2itw
                          (itw_field, atmosphere_dim,
                           p_grid[ Range( cloudbox_limits[0], 
                                          (cloudbox_limits[1]-
                                          cloudbox_limits[0]+1))],
                           lat_grid[ Range( cloudbox_limits[2], 
                                            (cloudbox_limits[3]-
                                             cloudbox_limits[2]+1))],
                           lon_grid[ Range( cloudbox_limits[4], 
                                            (cloudbox_limits[5]-
                                             cloudbox_limits[4]+1))],
                           cloud_gp_p, cloud_gp_lat, cloud_gp_lon );

                        // Ppath_step always has 2 points, the starting
                        // point and the intersection point.
                        Tensor3 ext_mat_int(stokes_dim, stokes_dim, 
                                            ppath_step.np);
                        Matrix abs_vec_int(stokes_dim, ppath_step.np);
                        Matrix sca_vec_int(stokes_dim, ppath_step.np);
                        Matrix sto_vec_int(stokes_dim, ppath_step.np);
                        Vector t_int(ppath_step.np);
                        Vector vmr_int(ppath_step.np);
                        Vector p_int(ppath_step.np);

                        // Interpolate ext_mat, abs_vec and sca_vec on the
                        // intersection point.
                        
                        // Calculate the average of the coefficients for the 
                        // layers to be considered in the 
                        // radiative transfer calculation.
                        
                        for (Index i = 0; i < stokes_dim; i++)
                          {
                            
                            // Extinction matrix requires a second loop 
                            // over stokes_dim
                            out3 << "Interpolate ext_mat:\n";
                            for (Index j = 0; j < stokes_dim; j++)
                              {
                                // Interpolation of ext_mat
                                //
                                interp_atmfield_by_itw
                                  (ext_mat_int(i, j, joker),
                                   atmosphere_dim,
                                   p_grid[p_range],
                                   lat_grid[lat_range],
                                   lon_grid[lon_range],
                                   ext_mat_field( joker, joker, joker, i , j),
                                   "ext_mat_field",
                                   cloud_gp_p, cloud_gp_lat, cloud_gp_lon,
                                   itw_field);
                                                            }
                            // Absorption vector:
                            //
                            // Interpolation of abs_vec
                            //
                            out3 << "Interpolate abs_vec:\n";
                            interp_atmfield_by_itw
                                  (abs_vec_int(i,joker),
                                   atmosphere_dim,
                                   p_grid[p_range],
                                   lat_grid[lat_range],
                                   lon_grid[lon_range],
                                   abs_vec_field( joker, joker, joker, i),
                                   "abs_vec_field",
                                   cloud_gp_p, cloud_gp_lat, cloud_gp_lon,
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
                               p_grid[p_range], lat_grid[lat_range],
                               lon_grid[lon_range],
                               scat_field(joker, joker, joker, scat_za_index,
                                         scat_aa_index, i),
                               "scat_field",
                               cloud_gp_p,
                               cloud_gp_lat,
                               cloud_gp_lon,
                               itw_field);
                            //
                            // Stokes vector:
                            //
                            // Interpolation of i_field_old.
                            out3 << "Interpolate i_field_old:\n";
                            interp_atmfield_by_itw
                              (sto_vec_int(i, joker),
                               atmosphere_dim,
                               p_grid[p_range],
                               lat_grid[lat_range],
                               lon_grid[lon_range],
                               i_field_old(joker, joker, joker, scat_za_index,
                                       scat_aa_index, i),
                               "i_field_old",
                               cloud_gp_p,
                               cloud_gp_lat,
                               cloud_gp_lon,
                               itw_field);
                            // 
                            // For the radiative transfer equation we 
                            // need the Stokes vector at the intersection
                            // point.
                            //
                            stokes_vec[i] = sto_vec_int(i,ppath_step.np);
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
                           p_grid, lat_grid,
                           lon_grid,
                           t_field(joker, joker, joker),
                           "t_field",
                           ppath_step.gp_p,
                           ppath_step.gp_lat,
                           ppath_step.gp_lon,
                           itw_field);

                         // 
                        // The vmr_field is needed for the gaseous absorption 
                        // calculation.
                        //
                        const Index N_species = vmr_field.nbooks();
                        //
                        // Interpolated vmr_list, holds a vmr_list for
                        //each point in 
                        // ppath_step.
                        //
                        Matrix vmr_list_int(N_species, ppath_step.np);

                        for (Index i = 0; i < N_species; i++)
                          {
                            out3 << "Interpolate vmr field\n";
                            interp_atmfield_by_itw
                              (vmr_int,
                               atmosphere_dim,
                               p_grid, lat_grid, lon_grid,
                               vmr_field(i, joker, joker, joker),
                               "vmr_field",
                               ppath_step.gp_p,
                               ppath_step.gp_lat,
                               ppath_step.gp_lon,
                               itw_field);
                  
                            vmr_list_int(i, joker) = vmr_int;
                          }
                        // 
                        // Interpolate pressure, latitude, longitude
                        //
                        itw2p( p_int, p_grid, ppath_step.gp_p, itw_field); 
                       
                        // Radiative transfer from one layer to the next,
                        // starting at the intersection with the next layer 
                        // and propagating to the considered point.
                        for( Index k= ppath_step.np-1; k > 0; k--)
                          {
                            // Length of the path between the two layers.
                            Numeric l_step = ppath_step.l_step[k-1];
                    
                            
                            // Average temperature
                            rte_temperature =   0.5 * (t_int[k] + t_int[k-1]);
                            //
                            // Average pressure
                            rte_pressure = 0.5 * (p_int[k] + p_int[k-1]);
                            //
                            // Average vmrs
                            for (Index i = 0; i < N_species; i++)
                              rte_vmr_list[i] = 0.5 * (vmr_list_int(i,k) + 
                                                     vmr_list_int(i,k-1));
                            //
                            // Calculate scalar gas absorption and add it to 
                            // abs_vec and ext_mat.
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
                                      (ext_mat_int(i,j,k) +
                                       ext_mat_int(i,j,k-1));
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
                            Numeric rte_planck_value = planck(f, rte_temperature);
                            
                            // Some messages:
                            out3 << "-------------------------------------\n";
                            out3 << "Input for radiative transfer step \n"
                                 << "calculation inside"
                                 << " the cloudbox:" << "\n";
                            out3 << "Stokes vector at intersection point: \n" 
                                 << stokes_vec 
                                 << "\n"; 
                            out3 << "l_step: ..." << l_step << "\n";
                            out3 << "--------------------------------------\n";
                            out3 << "Averaged coefficients: \n";
                            out3 << "Planck function: " 
                                 << rte_planck_value << "\n";
                            out3 << "Scattering vector: " 
                                 << sca_vec_av << "\n"; 
                            out3 << "Absorption vector: " 
                                 << abs_vec(0,joker) << "\n"; 
                            out3 << "Extinction matrix: " 
                                 << ext_mat(0,joker,joker) << "\n"; 
                      
                  
                            assert (!is_singular( ext_mat(0,joker,joker)));
                            
                            // Radiative transfer step calculation. 
                            // The Stokes vector is
                            // updated until the considered point is reached.
                            rte_step(stokes_vec, ext_mat(0,joker,joker), 
                                     abs_vec(0,joker), 
                                     sca_vec_av, l_step, rte_planck_value);
                          }
                        
                        
                        // Assign calculated Stokes Vector to i_field. 
                        i_field(p_index - cloudbox_limits[0],
                                lat_index - cloudbox_limits[2],
                                lon_index - cloudbox_limits[4],
                                scat_za_index, scat_aa_index,
                                joker) = stokes_vec;
                        //
                      } //end if
                    // 
                    // If the intersection point is outside the cloudbox
                    // no radiative transfer step is performed.
                    // The value on the cloudbox boundary remains unchanged.
                    //
                    
                    //
                  } // end of loop over lon_grid
                }  // end of loop over lat_grid
            } // end loop over p_grid
        }// end loop over scat_aa_grid
    }// end loop over scat_za_grid
  
  out2 << "Finished i_fieldUpdate3D.\n";
}// end of function.

                        
//! 3D RT calculation inside the cloud box.
/*! 
  This function loops over all grid points and all directions and performs 
  the RT calculation with a fixed scattering integral for one frequency 
  of the frequency grid specified by *f_index*. 
  
  The loop over directions is the outermost loop. Here the optical properties
  for single particle types are calculated as they are not depending on the 
  position of the particles. 
  The inner loop is the loop over the positions. Inside this loop the total 
  optical properties including the partile types as well as the gaseous
  species are calculated. Then the radiative transfer equation can be computed.

  Note: Ground reflection is not yet implemented in 3D scattering 
  calculations.

  WS Output:
  \param i_field Updated radiation field inside the cloudbox. 
  Variables used in scalar_gas_abs_agenda:
  \param rte_pressure
  \param rte_temperature
  \param rte_vmr_list
  Variables used in spt_calc_agenda:
  \param scat_za_index
  \param scat_aa_index
  Variables used in opt_prop_xxx_agenda:
  \param ext_mat
  \param abs_vec  
  \param scat_p_index
  \param scat_lat_index
  \param scat_lon_index
  Variables used in ppath_step_agenda:
  \param ppath_step
  WS Input:
  \param i_field_old Old radiation field.
  \param scat_field Scattered field.
  \param cloudbox_limits 
  Calculate scalar gas absorption:
  \param scalar_gas_absorption_agenda
  \param vmr_field
  Optical properties for single particle type:
  \param spt_calc_agenda
  \param scat_za_grid
  \param scat_aa_grid
  Optical properties for gases and particles:
  \param opt_prop_part_agenda
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

  \author Claudia Emde
  \date 2003-01-07
*/
void i_fieldUpdateSeq3D(// WS Output:
                     Tensor6& i_field,
                     // scalar_gas_abs_agenda:
                     Numeric& rte_pressure,
                     Numeric& rte_temperature,
                     Vector& rte_vmr_list,
                     // spt_calc_agenda:
                     Index& scat_za_index,
                     Index& scat_aa_index,
                     // opt_prop_xxx_agenda:
                     Tensor3& ext_mat,
                     Matrix& abs_vec,  
                     Index& scat_p_index,
                     Index& scat_lat_index,
                     Index& scat_lon_index,
                     // ppath_step_agenda:
                     Ppath& ppath_step, 
                     // WS Input:
                     const Tensor6& scat_field,
                     const ArrayOfIndex& cloudbox_limits,
                     // Calculate scalar gas absorption:
                     const Agenda& scalar_gas_absorption_agenda,
                     const Tensor4& vmr_field,
                     // Optical properties for single particle type:
                     const Agenda& spt_calc_agenda,
                     const Vector& scat_za_grid,
                     const Vector& scat_aa_grid,
                     // Optical properties for gases and particles:
                     const Agenda& opt_prop_part_agenda,
                     const Agenda& opt_prop_gas_agenda,
                     // Propagation path calculation:
                     const Agenda& ppath_step_agenda,
                     const Vector& p_grid,
                     const Vector& lat_grid,
                     const Vector& lon_grid,
                     const Tensor3& z_field,
                     const Matrix& r_geoid,
                     // Calculate thermal emission:
                     const Tensor3& t_field,
                     const Vector& f_grid,
                     const Index& f_index
                     )
{
  
  out2 << "i_fieldUpdate3D: Radiative transfer calculatiuon in cloudbox. \n";
  out2 << "------------------------------------------------------------- \n";

  const Index stokes_dim = scat_field.ncols();
    
  //Check the input
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) + 1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
                       scat_za_grid.nelem(), 
                       scat_aa_grid.nelem(),
                       stokes_dim));

  assert( is_size( scat_field, 
                     (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                     (cloudbox_limits[3] - cloudbox_limits[2]) + 1, 
                     (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
   
  // Is the frequency index valid?
  assert( f_index <= f_grid.nelem() );

  // End of checks



  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  // Number of azimuth angles
  const Index N_scat_aa = scat_aa_grid.nelem();

 
  //=======================================================================
  // Calculate coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are 
  // calulated for interpolated VMRs, temperature and pressure.
  
   // Define shorter names for cloudbox_limits.

  const Index p_low = cloudbox_limits[0];
  const Index p_up = cloudbox_limits[1];
  const Index lat_low = cloudbox_limits[2];
  const Index lat_up = cloudbox_limits[3];
  const Index lon_low = cloudbox_limits[4];
  const Index lon_up = cloudbox_limits[5];

  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(p_up-p_low+1, lat_up-lat_low+1, lon_up-lon_low+1,
                        stokes_dim, stokes_dim, 0.);
  Tensor4 abs_vec_field(p_up-p_low+1, lat_up-lat_low+1, lon_up-lon_low+1,
                        stokes_dim, 0.);
 
  
  //Loop over all directions, defined by scat_za_grid 
  for(scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      //Loop over azimuth directions (scat_aa_grid). First and last point in 
      // azimuth angle grid are euqal. Start with second element.
      for(scat_aa_index = 1; scat_aa_index < N_scat_aa; scat_aa_index ++)
        {
         //==================================================================
          // Radiative transfer inside the cloudbox
          //==================================================================

          // This function has to be called inside the angular loop, as
          // it spt_calc_agenda takes *scat_za_index* and *scat_aa_index* 
          // from the workspace.
          cloud_fieldsCalc(ext_mat_field, abs_vec_field, scat_p_index,
                           scat_lat_index, scat_lon_index,
                           ext_mat, abs_vec, 
                           spt_calc_agenda, 
                           opt_prop_part_agenda, cloudbox_limits);
          

          Vector stokes_vec(stokes_dim,0.);
          
          Numeric theta_lim = 180. - asin((r_geoid(0,0)+z_field(p_low,0,0))
                                         /(r_geoid(0,0)+z_field(p_up,0,0)))
            *RAD2DEG;

          // Sequential update for uplooking angles
          if ( scat_za_grid[scat_za_index] <= 90.) 
            {
              // Loop over all positions inside the cloud box defined by the 
              // cloudbox_limits exculding the upper boundary. For uplooking
              // directions, we start from cloudbox_limits[1]-1 and go down
              // to cloudbox_limits[0] to do a sequential update of the
              // aradiation field
              for(Index p_index = p_up-1; p_index >= p_low; p_index --)
                {
                  for(Index lat_index = lat_low; lat_index <= lat_up; 
                      lat_index ++)
                    {
                      for(Index lon_index = lon_low; lon_index <= lon_up; 
                          lon_index ++)
                        {
                          cloud_ppath_update3D(i_field, 
                                               rte_pressure, rte_temperature, 
                                               rte_vmr_list, ext_mat, abs_vec,
                                               ppath_step, p_index, lat_index, 
                                               lon_index, scat_za_index, 
                                               scat_aa_index, scat_za_grid, 
                                               scat_aa_grid, cloudbox_limits, 
                                               scat_field, 
                                               scalar_gas_absorption_agenda,
                                               vmr_field, 
                                               opt_prop_gas_agenda,
                                               ppath_step_agenda, p_grid, 
                                               lat_grid, lon_grid, z_field, 
                                               r_geoid, t_field, f_grid,
                                               f_index,
                                               ext_mat_field, abs_vec_field);
                        }
                    }
                }
            }// close up-looking case
          else if ( scat_za_grid[scat_za_index] > theta_lim) 
            {
              //
              // Sequential updating for downlooking angles
              //
              for(Index p_index = p_low+1; p_index <= p_up; p_index ++)
                {
                  for(Index lat_index = lat_low; lat_index <= lat_up; 
                      lat_index ++)
                    {
                      for(Index lon_index = lon_low; lon_index <= lon_up; 
                          lon_index ++)
                        {
                          cloud_ppath_update3D(i_field, 
                                               rte_pressure, rte_temperature, 
                                               rte_vmr_list, ext_mat, abs_vec,
                                               ppath_step, p_index, lat_index, 
                                               lon_index, scat_za_index, 
                                               scat_aa_index, scat_za_grid, 
                                               scat_aa_grid, cloudbox_limits, 
                                               scat_field, 
                                               scalar_gas_absorption_agenda,
                                               vmr_field, 
                                               opt_prop_gas_agenda,
                                               ppath_step_agenda, p_grid, 
                                               lat_grid, lon_grid, z_field, 
                                               r_geoid, t_field, f_grid,
                                               f_index,
                                               ext_mat_field, abs_vec_field);
                        }
                    }
                }
            } // end if downlooking.
      
          //
          // Limb looking:
          // We have to include a special case here, as we may miss the endpoints
          // when the intersection point is at the same level as the actual point.
          // To be save we loop over the full cloudbox. Inside the function 
          // cloud_ppath_update3D it is checked whether the intersection point is 
          // inside the cloudbox or not.
          else if (  scat_za_grid[scat_za_index] > 90. &&
                     scat_za_grid[scat_za_index] < theta_lim ) 
            {
              for(Index p_index = p_low; p_index <= p_up; p_index ++)
                {
                  // For this case the cloudbox goes down to the ground an we
                  // look downwards. These cases are outside the cloudbox and 
                  // not needed. Switch is included here, as ppath_step_agenda 
                  // gives an error for such cases.
                  if (!(p_index == 0 && scat_za_grid[scat_za_index] > 90.))
                    {
                      for(Index lat_index = lat_low; lat_index <= lat_up; 
                          lat_index ++)
                        {
                          for(Index lon_index = lon_low; lon_index <= lon_up; 
                              lon_index ++)
                            {
                              cloud_ppath_update3D(i_field, 
                                                   rte_pressure,
                                                   rte_temperature, 
                                                   rte_vmr_list,
                                                   ext_mat, abs_vec,
                                                   ppath_step, p_index, 
                                                   lat_index, 
                                                   lon_index, scat_za_index, 
                                                   scat_aa_index,
                                                   scat_za_grid, 
                                                   scat_aa_grid,
                                                   cloudbox_limits, 
                                                   scat_field, 
                                                   scalar_gas_absorption_agenda,
                                                   vmr_field, 
                                                   opt_prop_gas_agenda,
                                                   ppath_step_agenda, p_grid, 
                                                   lat_grid, lon_grid,
                                                   z_field, 
                                                   r_geoid, t_field, f_grid,
                                                   f_index,
                                                   ext_mat_field,
                                                   abs_vec_field); 
                            }
                        }
                    }
                }
            }
        } //  Closes loop over aa_grid.
    }// Closes loop over scat_za_grid.

  i_field(joker, joker, joker, joker, 0, joker) = 
    i_field(joker, joker, joker, joker, N_scat_aa-1, joker);

  
  
}// end of function.          
      



//! 1D RT calculation inside the cloud box in a plane parallel geometry.
/*! 
  This function is baseically same as i_fieldUpdateSeq1D in that it updates the
  i_field.  The difference is that it assumes that inside the cloudbox the 
  atmosphere is planeparallel.  This is included with the intention that it 
  can be faster compared to the spherical.  Also it will be good for comparisons.
  

  WS Output:
  \param i_field       Updated intensity field. 
  Variables used in scalar_gas_abs_agenda:
  \param rte_pressure    FIXME: Add documentation.
  \param rte_temperature FIXME: Add documentation.
  \param rte_vmr_list    FIXME: Add documentation.
  Variables used in spt_calc_agenda:
  \param scat_za_index
   Variables used in opt_prop_xxx_agenda:
  \param ext_mat
  \param abs_vec  
  \param scat_p_index
  Variables used in ppath_step_agenda:
  \param ppath_step
  \param ppath_step
  WS Input:
  \param i_field_old Old radiation field.
  \param scat_field Scattered field.
  \param cloudbox_limits 
  Calculate scalar gas absorption:
  \param scalar_gas_absorption_agenda
  \param vmr_field
  Optical properties for single particle type:
  \param spt_calc_agenda
  \param scat_za_grid
  Optical properties for gases and particles:
  \param opt_prop_part_agenda
  \param pnd_field
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

  \author Sreerekha Ravi
  \date 2002-11-28
*/
void
i_fieldUpdateSeq1D_PlaneParallel(// WS Output:
                Tensor6& i_field,
                // scalar_gas_abs_agenda:
                Numeric& rte_pressure,
                Numeric& rte_temperature,
                Vector& rte_vmr_list,
                // spt_calc_agenda:
                Index& scat_za_index ,
                // opt_prop_xxx_agenda:
                Tensor3& ext_mat,
                Matrix& abs_vec,  
                Index& scat_p_index,
                // ppath_step_agenda:
                Ppath& ppath_step, 
                // ground related variables STR
                Matrix& ground_los,
                Matrix& ground_emission,
                Tensor4& ground_refl_coeffs,
                Vector& rte_los,
                Vector& rte_pos,
                GridPos& rte_gp_p,
                // WS Input:
                const Tensor6& scat_field,
                const ArrayOfIndex& cloudbox_limits,
                // Calculate scalar gas absorption:
                const Agenda& scalar_gas_absorption_agenda,
                const Tensor4& vmr_field,
                // Optical properties for single particle type:
                const Agenda& spt_calc_agenda,
                const Vector& scat_za_grid,
                // Optical properties for gases and particles:
                const Agenda& opt_prop_part_agenda,
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
                const Agenda& ground_refl_agenda //STR
                )
{

  out2 << "i_fieldUpdatePlaneparallel: Radiative transfer calculatiuon in cloudbox.\n";
  out2 << "------------------------------------------------------------- \n";
  
  const Index stokes_dim = scat_field.ncols();
  //  const Index atmosphere_dim = 1;

  //Check the input
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));

  assert( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));  
  
  // Is the frequency index valid?
  assert( f_index <= f_grid.nelem() );

  // End of checks



  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();

  
  
  //=======================================================================
  // Calculate scattering coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are 
  // calulated for interpolated VMRs, temperature and pressure.
  
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:
     
      Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                            stokes_dim, stokes_dim, 0.);
      Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                            stokes_dim, 0.);

      //Loop over all directions, defined by scat_za_grid 
  for(scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      
      //Only dummy variables:
      Index scat_lat_index = 0;
      Index scat_lon_index = 0;
      
      cloud_fieldsCalc(ext_mat_field, abs_vec_field, scat_p_index,
                       scat_lat_index, scat_lon_index,
                       ext_mat, abs_vec, 
                       spt_calc_agenda, 
                       opt_prop_part_agenda, cloudbox_limits);

      //======================================================================
      // Radiative transfer inside the cloudbox
      //=====================================================================
      
      Vector stokes_vec(stokes_dim,0.);
       // Sequential update for uplooking angles
      if ( scat_za_grid[scat_za_index] <= 90) 
        {
          // Loop over all positions inside the cloud box defined by the 
          // cloudbox_limits exculding the upper boundary. For uplooking
          // directions, we start from cloudbox_limits[1]-1 and go down
          // to cloudbox_limits[0] to do a sequential update of the
          // aradiation field
     
          // Loop over all positions inside the cloudbox defined by the 
          // cloudbox_limits.
          for(Index p_index = cloudbox_limits[1] -1; p_index
                >= cloudbox_limits[0]; p_index --)
            {
              cloud_ppath_update1D_planeparallel(i_field, 
                                                 rte_pressure, rte_temperature, rte_vmr_list,
                                                 ext_mat, abs_vec, ground_los,
                                                 ground_emission, ground_refl_coeffs,
                                                 rte_los, rte_pos, rte_gp_p, ppath_step, 
                                                 p_index, scat_za_index, scat_za_grid,
                                                 cloudbox_limits, scat_field,
                                                 scalar_gas_absorption_agenda, vmr_field,
                                                 opt_prop_gas_agenda, ppath_step_agenda,
                                                 p_grid,  z_field, r_geoid, t_field, 
                                                 f_grid, f_index, ext_mat_field,
                                                 abs_vec_field,ground_refl_agenda); 
            }   
        }
      else if ( scat_za_grid[scat_za_index] > 90) 
        {
          //
          // Sequential updating for downlooking angles
          //
          for(Index p_index = cloudbox_limits[0]+1; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
              cloud_ppath_update1D_planeparallel(i_field,  
                                                 rte_pressure, rte_temperature, rte_vmr_list,
                                                 ext_mat, abs_vec, ground_los,
                                                 ground_emission, ground_refl_coeffs,
                                                 rte_los, rte_pos, rte_gp_p, ppath_step, 
                                                 p_index, scat_za_index, scat_za_grid,
                                                 cloudbox_limits, scat_field,
                                                 scalar_gas_absorption_agenda, vmr_field,
                                                 opt_prop_gas_agenda, ppath_step_agenda,
                                                 p_grid,  z_field, r_geoid, t_field, 
                                                 f_grid, f_index, ext_mat_field, 
                                                 abs_vec_field,ground_refl_agenda); 
            }// Close loop over p_grid (inside cloudbox).
        } // end if downlooking.
      
    }// Closes loop over scat_za_grid.
}// End of the function.


//! This method computes the scattering integral
/*! 
    By scattering integral, we mean the field generated by integrating
    the product of intensity field and phase matrix over all incident 
    angles.

    This method uses the same angular grids as the radiative transfer methods
    *i_fieldUpdateXXX*. 
    
    WS Output:
    \param scat_field Scattering integral field.
    WS Output and Input:
    \param pha_mat Phase matrix.
    \param pha_mat_spt Phase matrix for single particle type
    \param scat_za_index Indices for propagation directions,
    \param scat_aa_index needed for communication with pha_mat_spt_agenda.
    WS Input:
    \param pha_mat_spt_agenda Agenda for calculation of pha_mat_spt. 
    \param i_field Intensity field 
    \param pnd_field Particle number density field.
    \param scat_za_grid Zenith angle grid
    \param scat_aa_grid Azimuth angle grid
    \param atmosphere_dim Atmospheric dimension
    \param cloudbox_limits Limits of the cloudbox.
    \param za_gridsize Number of grid points in zenith angle grid.

\author Sreerekha Ravi, Claudia Emde
\date 2002-06-20

\date 2003-10-24 Modifications to gain efficiency (CE).  

*/

void
scat_fieldCalc(//WS Output:
               Tensor6& scat_field,
               Tensor4& pha_mat,
               Tensor5& pha_mat_spt,
               Index& scat_za_index, 
               Index& scat_aa_index,
               //WS Input:
               const Agenda& pha_mat_spt_agenda,
               const Tensor6& i_field,
               const Tensor4& pnd_field,
               const Vector& scat_za_grid,
               const Vector& scat_aa_grid,
               const Index& atmosphere_dim,
               const ArrayOfIndex& cloudbox_limits,
               const Index& za_grid_size
               )
  
{
  Index stokes_dim = scat_field.ncols();

  // Some useful indices :
  Index Nza = scat_za_grid.nelem();
  Index Naa = scat_aa_grid.nelem();

  if (za_grid_size != Nza)
    throw runtime_error("Error in *scat_fieldCalc*: \n"
                        "Zenith angle grids for radiative transfer \n"
                        "(*scat_za_grid*) and computation of scattering  \n"
                        "integral (specified in *grid_stepsizeSet*) must \n"
                        "be equal.\n"
                        "If you need different grids (for Limb calculations)\n"
                        "you have to use *scat_fieldCalcLimb*"
                        );

  Vector grid_stepsize(2);
  grid_stepsize[0] = 180./(za_grid_size - 1);
  grid_stepsize[1] = 360./(Naa - 1);     
  
  Tensor3 product_field(Nza, Naa, stokes_dim, 0);
 
  out2 << "Calculate the scattered field\n";
  
  if  ( atmosphere_dim == 1 )
    {
      assert ( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));
    assert ( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));
   
    scat_aa_index = 0;
    
    // Get pha_mat at the grid positions
    // Since atmosphere_dim = 1, there is no loop over lat and lon grids
    for (Index p_index = 0; p_index <= cloudbox_limits[1]- cloudbox_limits[0] ;
         p_index++)
      {
        //There is only loop over zenith angle grid ; no azimuth angle grid.
        for (scat_za_index = 0; scat_za_index < Nza; scat_za_index ++)
          {
            // Calculate the phase matric of a single particle type
            pha_mat_spt_agenda.execute(scat_za_index || p_index );

            // Sum over all particle types
            pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
                         atmosphere_dim, p_index + cloudbox_limits[0], 0, 
                         0);

             product_field = 0;

            // za_in and aa_in are for incoming zenith and azimutha 
            //angle direction for which pha_mat is calculated
            for (Index za_in = 0; za_in < Nza; ++ za_in)
              { 
                for (Index aa_in = 0; aa_in < Naa; ++ aa_in)
                  {
                    // Multiplication of phase matrix with incoming 
                    // intensity field.

                    for ( Index i = 0; i < stokes_dim; i++)
                       {
                       for (Index j = 0; j< stokes_dim; j++)
                        {
                            product_field(za_in, aa_in, i) +=
                              pha_mat(za_in, aa_in, i, j) * 
                              i_field(p_index, 0, 0, za_in, 0, j);
                        }
                       }
                    
                  }//end aa_in loop
              }//end za_in loop
            //integration of the product of ifield_in and pha
            //  over zenith angle and azimuth angle grid. It calls
            for (Index i = 0; i < stokes_dim; i++)
              {
                
                MatrixView product_field_mat =
                  product_field( Range(joker),
                                 Range(joker),
                                 i);
                // scat_field is also defined for all points inside the cloud
                //box for each propagion angle
                scat_field( p_index,
                            0,
                            0,
                            scat_za_index, 
                            0,
                            i)  =   AngIntegrate_trapezoid_opti(product_field_mat,
                                                                scat_za_grid,
                                                                scat_aa_grid,
                                                                grid_stepsize);
                
              }//end i loop
          }//end za_prop loop
      }//end p_index loop
    
  }//end atmosphere_dim = 1

  
  //When atmospheric dimension , atmosphere_dim = 3
  if( atmosphere_dim == 3 ){

     assert ( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
    assert ( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
    
    /*there is a loop over pressure, latitude and longitudeindex
      when we calculate the pha_mat from pha_mat_spt and pnd_field
      using the method pha_matCalc.  */
    
    for (Index p_index = 0; p_index <= cloudbox_limits[1] - cloudbox_limits[0];
         p_index++)
      {
        for (Index lat_index = 0; lat_index <= 
               cloudbox_limits[3] - cloudbox_limits[2]; lat_index++)
          {
            for (Index lon_index = 0; lon_index <= 
                   cloudbox_limits[5]-cloudbox_limits[4]; lon_index++)
              {
                for (scat_aa_index = 1; scat_aa_index < Naa; 
                     scat_aa_index++)
                  {
                    for (scat_za_index = 0; scat_za_index < Nza; 
                         scat_za_index ++)
                      {
                        pha_mat_spt_agenda.execute( true );

                        pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
                                    atmosphere_dim, 
                                    p_index + cloudbox_limits[0], 
                                    lat_index + cloudbox_limits[2], 
                                    lon_index + cloudbox_limits[4]);
                   
                        product_field = 0;

                        //za_in and aa_in are the incoming directions
                        //for which pha_mat_spt is calculated
                        for (Index za_in = 0; za_in < Nza; ++ za_in)
                          {
                            for (Index aa_in = 0; aa_in < Naa; ++ aa_in)
                              { 
                                // Multiplication of phase matrix
                                // with incloming intensity field.
                                for ( Index i = 0; i < stokes_dim; i++)
                                  {
                                    for (Index j = 0; j< stokes_dim; j++)
                                      {
                                        product_field(za_in, aa_in, i) +=
                                          pha_mat(za_in, aa_in, i, j) * 
                                          i_field(p_index, lat_index, 
                                                  lon_index, scat_za_index,
                                                  scat_aa_index, j);
                                      }
                                  }
                              }//end aa_in loop
                          }//end za_in loop
                        //integration of the product of ifield_in and pha
                        //over zenith angle and azimuth angle grid. It 
                        //calls here the integration routine 
                        //AngIntegrate_trapezoid_opti
                        for (Index i = 0; i < stokes_dim; i++)
                          {
                            scat_field( p_index,
                                        lat_index,
                                        lon_index,
                                        scat_za_index, 
                                        scat_aa_index,
                                        i)  =  
                              AngIntegrate_trapezoid_opti(product_field
                                                          ( joker,
                                                            joker, i),
                                                          scat_za_grid,
                                                          scat_aa_grid,
                                                          grid_stepsize);
                          }//end i loop
                      }//end aa_prop loop
                  }//end za_prop loop
              }//end lon loop
          }// end lat loop
      }// end p loop
    scat_field(joker, joker, joker, joker, 0, joker) =
      scat_field(joker, joker, joker, joker, Naa-1, joker);
  }// end atmosphere_dim = 3
}



//! This method computes the scattering integral (Limb)
/*! 
    By scattering integral, we mean the field generated by integrating
    the product of intensity field and phase matrix over all incident 
    angles. 
    
    This methods used the grids specified in *grid_stepsizeSet*. The zenith
    angle grid used in the radiative transfer part has to be very fine around
    90°, i.e., different grids are used in the two parts of the scattering
    calculation. 
    
    In the case of different grids the intensity field and the scattered field 
    have to be interpolated.
 
    WS Output:
    \param scat_field Scattering integral field.
    WS Output and Input:
    \param pha_mat Phase matrix.
    \param pha_mat_spt Phase matrix for single particle type
    \param scat_za_index Indices for propagation directions,
    \param scat_aa_index needed for communication with pha_mat_spt_agenda.
    WS Input:
    \param pha_mat_spt_agenda Agenda for calculation of pha_mat_spt. 
    \param i_field Intensity field 
    \param pnd_field Particle number density field.
    \param scat_za_grid Zenith angle grid
    \param scat_aa_grid Azimuth angle grid
    \param atmosphere_dim Atmospheric dimension
    \param cloudbox_limits Limits of the cloudbox.
    \param za_gridsize Number of grid points in zenith angle grid.

\author Claudia Emde
\date 2003-11-28

*/
void
scat_fieldCalcLimb(//WS Output:
               Tensor6& scat_field,
               Tensor4& pha_mat,
               Tensor5& pha_mat_spt,
               Index& scat_za_index, 
               Index& scat_aa_index,
               //WS Input:
               const Agenda& pha_mat_spt_agenda,
               const Tensor6& i_field,
               const Tensor4& pnd_field,
               const Vector& scat_za_grid,
               const Vector& scat_aa_grid,
               const Index& atmosphere_dim,
               const ArrayOfIndex& cloudbox_limits,
               const Index& za_grid_size
               )
  
{
  Index stokes_dim = scat_field.ncols();

  // Some useful indices :
  Index Nza = scat_za_grid.nelem();
  Index Naa = scat_aa_grid.nelem();

   // Create the grids for the calculation of the scattering integral.
  Vector za_grid;
  nlinspace(za_grid, 0, 180, za_grid_size);
 
  // Two interpolations are required. First we have to interpolate the 
  // intensity field on the equidistant grid: 
  ArrayOfGridPos gp_za_i(za_grid_size);
  gridpos(gp_za_i, scat_za_grid, za_grid);
  
  Matrix itw_za_i(za_grid_size, 2);
  interpweights(itw_za_i, gp_za_i);

  // Intensity field interpolated on equidistant grid.
  Matrix i_field_int(za_grid_size, stokes_dim, 0);

  // Second, we have to interpolate the scattering integral on the RT
  // zenith angle grid.
  ArrayOfGridPos gp_za(Nza);
  gridpos(gp_za, za_grid, scat_za_grid);

  Matrix itw_za(Nza, 2);
  interpweights(itw_za, gp_za);

  // Original scattered field, on equidistant zenith angle grid.
  Matrix scat_field_org(za_grid_size, stokes_dim, 0);
  
  //  Grid stepsize of zenith and azimuth angle grid, these are needed for the 
  // integration function. 
  Vector grid_stepsize(2);
  grid_stepsize[0] = 180./(za_grid_size - 1);
  grid_stepsize[1] = 360./(Naa - 1);
    
  Tensor3 product_field(za_grid_size, Naa, stokes_dim, 0);
  
  if  ( atmosphere_dim == 1 )
    {
    
      assert ( is_size( i_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        1, 
                        1,
                        scat_za_grid.nelem(), 
                        1,
                        stokes_dim));
      assert ( is_size( scat_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        1, 
                        1,
                        scat_za_grid.nelem(), 
                        1,
                        stokes_dim));
    
      scat_aa_index = 0;
  
      // Get pha_mat at the grid positions
      // Since atmosphere_dim = 1, there is no loop over lat and lon grids
      for (Index p_index = 0; p_index <= cloudbox_limits[1]-cloudbox_limits[0];
           p_index++)
        {
          // Interpolate intensity field:
          for (Index i = 0; i < stokes_dim; i++)
            {
              interp(i_field_int(joker, i), itw_za_i, 
                     i_field(p_index, 0, 0, joker, 0, i), gp_za_i);
            }       
          
          //There is only loop over zenith angle grid; no azimuth angle grid.
          for( scat_za_index = 0; scat_za_index < za_grid_size;
               scat_za_index ++)
            {
            // Calculate the phase matrix of a single particle type
            out3 << "Calculate the phase matrix \n"; 
            pha_mat_spt_agenda.execute(true);

            // Sum over all particle types
            pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
                        atmosphere_dim, p_index + cloudbox_limits[0], 0, 
                        0);

            out3 << "Multiplication of phase matrix with incoming" << 
              " intensities \n";
            
            product_field = 0;

            // za_in and aa_in are for incoming zenith and azimuth 
            // angle direction for which pha_mat is calculated
            for( Index za_in = 0; za_in < za_grid_size; za_in ++)
              {
                for (Index aa_in = 0; aa_in < Naa; ++ aa_in)
                  {
                    // Multiplication of phase matrix with incoming 
                    // intensity field.
                    
                    for ( Index i = 0; i < stokes_dim; i++)
                      {
                        for (Index j = 0; j< stokes_dim; j++)
                          {
                            product_field(za_in, aa_in, i) +=
                              pha_mat(za_in, aa_in, i, j) * 
                              i_field_int(za_in, j);
                          }
                      }
                    
                  }//end aa_in loop
              }//end za_in loop
            
            out3 << "Compute integral. \n"; 
            for (Index i = 0; i < stokes_dim; i++)
              {
                
                MatrixView product_field_mat =
                  product_field( Range(joker),
                                 Range(joker),
                                 i);
                // scat_field is also defined for all points inside the cloud
                //box for each propagion angle
                scat_field_org(scat_za_index, i)=
                  AngIntegrate_trapezoid_opti(product_field_mat,
                                              za_grid,
                                              scat_aa_grid,
                                              grid_stepsize);
                
              }//end i loop
          }//end za_prop loop
        
        // Interpolation on scat_za_grid, which is used in radiative transfer 
        // part.
        for (Index i = 0; i < stokes_dim; i++)
          { 
            interp(scat_field(p_index,
                              0,
                              0,
                              joker,
                              0,
                              i),
                   itw_za,
                   scat_field_org(joker, i),
                   gp_za);
          }
      }//end p_index loop
    
  }//end atmosphere_dim = 1

  
  if( atmosphere_dim == 3 ){
    

    assert ( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
    assert ( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
      
    // Loop over all positions
    for (Index p_index = 0; p_index <= cloudbox_limits[1] - cloudbox_limits[0];
         p_index ++)
      {
        for (Index lat_index = 0; lat_index <= 
               cloudbox_limits[3] - cloudbox_limits[2]; lat_index++)
          {
            for (Index lon_index = 0; lon_index <= 
                   cloudbox_limits[5] - cloudbox_limits[4]; lon_index++)
              {
               
                // Loop over scattered directions
                for (scat_aa_index = 1; scat_aa_index < Naa; 
                     scat_aa_index ++)
                  {
                   // Interpolate intensity field:
                    for (Index i = 0; i < stokes_dim; i++)
                      {
                        interp(i_field_int(joker, i), itw_za_i, 
                               i_field(p_index, lat_index, lon_index,
                                       joker, scat_aa_index, i), gp_za_i);
                      }       
                    
                    for (scat_za_index = 0; scat_za_index < 
                           za_grid_size; scat_za_index ++)
                      {
                        out3 << "Calculate phase matrix \n";
                        pha_mat_spt_agenda.execute(true);
                        
                        pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
                                    atmosphere_dim, 
                                    p_index + cloudbox_limits[0], 
                                    lat_index + cloudbox_limits[2], 
                                    lon_index + cloudbox_limits[4]);
                        
                        product_field = 0;
                        
                        
                        //za_in and aa_in are the incoming directions
                        //for which pha_mat_spt is calculated
                        out3 << "Multiplication of phase matrix with" << 
                          "incoming intensity \n";
                        
                        for( Index za_in = 0; za_in < za_grid_size; za_in ++)
                          {
                            for (Index aa_in = 0; aa_in < Naa; ++ aa_in)
                              { 
                                // Multiplication of phase matrix
                                // with incloming intensity field.
                                for ( Index i = 0; i < stokes_dim; i++)
                                  {
                                    for (Index j = 0; j< stokes_dim; j++)
                                      {
                                        product_field(za_in, aa_in, i) +=
                                          pha_mat(za_in, aa_in, i, j) * 
                                          i_field_int(za_in, j);
                                      }
                                  }
                              }//end aa_in loop
                          }//end za_in loop
                        
                        out3 << "Compute the integral \n";

                        for (Index i = 0; i < stokes_dim; i++)
                          {
                            scat_field_org(scat_za_index, i)  =  
                              AngIntegrate_trapezoid_opti(product_field
                                                     ( joker,
                                                       joker, i),
                                                        za_grid,
                                                        scat_aa_grid,
                                                        grid_stepsize
                                                     );
                          }//end stokes_dim loop

                      }//end za_prop loop
                    //Interpolate on original za_grid. 
                    for (Index i = 0; i < stokes_dim; i++)
                      {
                        interp(scat_field(p_index,
                                          lat_index,
                                          lon_index,
                                          joker,
                                          scat_aa_index,
                                          i),
                               itw_za,
                               scat_field_org(joker, i),
                               gp_za);
                      }
                  } // end aa_prop loop
              }//end lon loop
          }//end lat loop
      }// end p loop
    scat_field(joker, joker, joker, joker, 0, joker) =
      scat_field(joker, joker, joker, joker, Naa-1, joker);
  }// end atm_dim=3

  out2 << "Finished scattered field.\n"; 
 
}

//! This method computes the scattering integral

/*! 
By scattering integral, we mean the field generated by integrating
the product of intensity field and phase matrix over all incident 
angles.  
  
\param scat_field Output : scattering integraal
\param pha_mat Output : phase matrix.
\param amp_mat Input : Amplitude matrix.
\param i_field Input : the intensity field 
\param pha_mat_spt Input : the phase matrix for single particle type
\param pnd_field Input : the particle number density field.
\param scat_za_grid zenith angle grid
\param scat_aa_grid azimuth angle grid
\param p_grid pressure grid
\param lat_grid latitude grid
\param lon_grid longitude grid
\param atmosphere_dim atmospheric dimension
\param cloudbox_limits Limits of the cloudbox.
\param grid_stepsize FIXME: Add documentation.

\author Sreerekha Ravi
\date 2002-06-20
*/

void
scat_fieldCalcFromAmpMat(//WS Output:
                         Tensor6& scat_field,
                         Tensor4& pha_mat,
                         Tensor5& pha_mat_spt,
                         //WS Input: 
                         const Tensor6& amp_mat,
                         const Tensor6& i_field,
                         const Tensor4& pnd_field,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         // FIXME parameters only used in assertions
#ifndef NDEBUG
                         const Vector& p_grid,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
#else
                         const Vector&,
                         const Vector&,
                         const Vector&,
#endif
                         const Index& atmosphere_dim,
                         const ArrayOfIndex& cloudbox_limits,
                         const Vector& grid_stepsize
                         )

  
{
  // Some useful indices :
// FIXME variable only used in assertion
#ifndef NDEBUG
  Index N_pt = pha_mat_spt.nshelves();
#endif
  Index Nza = scat_za_grid.nelem();
  Index Nza_prop = i_field.npages();
  Index Naa = scat_aa_grid.nelem();
  Index Naa_prop = i_field.ncols();
  
  const Index stokes_dim = pha_mat_spt.ncols();

  cout<<"Naa in the scattering integral routine"<<" "<<Naa<<"\n";
  //Tensor4 product_field(Np,Nza, Naa, stokes_dim);//earlier tensor3; added pressure index STR
  // now tensor5 after za_index_in index
  Tensor3 product_field(Nza, Naa, stokes_dim,0);
 
  
  /*
  cout<<"N_pt"<<" "<<N_pt;
  cout<<"N_p_grid"<<" "<<p_grid.nelem();
  cout<<"N_lat_grid"<<" "<<lat_grid.nelem();
  cout<<"N_lon_grid"<<" "<<lon_grid.nelem();*/

  //Note that the size of i_field and scat_field are the same 

  if (atmosphere_dim ==3){
    
    assert ( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
    assert ( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                      (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                      scat_za_grid.nelem(), 
                      scat_aa_grid.nelem(),
                      stokes_dim));
    assert ( is_size( pnd_field, 
                      N_pt,
                      p_grid.nelem(),
                      lat_grid.nelem(), 
                      lon_grid.nelem()));
    
  }
  
  else if (atmosphere_dim == 1 ){
    assert ( is_size( i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));
    assert ( is_size( scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));
    assert ( is_size( pnd_field, 
                      N_pt,
                      p_grid.nelem(),
                      1, 
                      1));
  }
  
   //When atmospheric dimension , atmosphere_dim = 1
  if( atmosphere_dim == 1 ){
  
    // Get pha_mat at the grid positions
    // Since atmosphere_dim = 1, there is no loop over lat and lon grids
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
         p_index++)
      {
       
        //There is only loop over zenith angle grid ; no azimuth angle grid.
        for (Index za_prop = 0; za_prop < Nza_prop; ++ za_prop)
          {
            pha_mat_sptCalc(pha_mat_spt,
                            amp_mat,
                            za_prop,
                            0);
            
            pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
                        atmosphere_dim, p_index, 0, 
                        0);
            // za_in and aa_in are for incoming zenith and azimutha 
            //angle direction for which pha_mat is calculated
            
            for (Index za_in = 0; za_in < Nza; ++ za_in)
              { 
                for (Index aa_in = 0; aa_in < Naa; ++ aa_in)
                  {
                    //This is a matrix with last two dimensions 
                    //being that of stokes_dim
                    ConstMatrixView pha = pha_mat(za_in,
                                                  aa_in,
                                                  Range(joker),
                                                  Range(joker));

                    //the incoming field expressed as a vector for all
                    //pressure points, and looking angles
                    ConstVectorView i_field_in = 
                      i_field((p_index - cloudbox_limits[0]),
                              0,
                              0,
                              za_prop,
                              0,
                              Range(joker));

                    // multiplication of intensity field vector and 
                    //pha_mat matrix
                    mult(product_field( za_in,
                                        aa_in,
                                        Range(joker)),
                         pha, 
                         i_field_in);
                            
                  }//end aa_in loop
              }//end za_in loop
            /*integration of the product of ifield_in and pha
              over zenith angle and azimuth angle grid. It calls
              here the integration routine AngIntegrate_trapezoid_opti*/
           
            for (Index i = 0; i < stokes_dim; i++)
              {
                
                MatrixView product_field_mat =
                  product_field( Range(joker),
                                 Range(joker),
                                 i);
                // scat_field is also defined for all points inside the cloud
                //box for each propagion angle
                scat_field( (p_index - cloudbox_limits[0]),
                            0,
                            0,
                            za_prop, 
                            0,
                            i)  =   AngIntegrate_trapezoid_opti(product_field_mat,
                                                                scat_za_grid,
                                                                scat_aa_grid,
                                                                grid_stepsize);
                
              }//end i loop
          }//end za_prop loop
      }//end p_index loop
    
  }//end atmosphere_dim = 1
  //  xml_write_to_file("scat_field.xml", scat_field);    
  //cout<<"scat_field"<<scat_field<<"\n";       
  //arts_exit ();
  
  
  //When atmospheric dimension , atmosphere_dim = 3
  if( atmosphere_dim == 3 ){
    
    /*there is a loop over pressure, latitude and longitudeindex
      when we calculate the pha_mat from pha_mat_spt and pnd_field
      using the method pha_matCalc.  */
    
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
         p_index++)
      {
        for (Index lat_index = cloudbox_limits[2]; lat_index <= 
               cloudbox_limits[3]; lat_index++)
          {
            for (Index lon_index = cloudbox_limits[4]; lon_index <= 
                   cloudbox_limits[5]; lon_index++)
              {
                pha_matCalc(pha_mat, pha_mat_spt, pnd_field, 
                            atmosphere_dim, p_index, lat_index, 
                            lon_index);
              }
          }
        
      }
    
    // Get pha_mat at the grid positions
    for (Index p_index = cloudbox_limits[0]; p_index <= cloudbox_limits[1];
         p_index++)
      {
        for (Index lat_index = cloudbox_limits[2]; lat_index <= 
               cloudbox_limits[3]; lat_index++)
          {
            for (Index lon_index = cloudbox_limits[4]; lon_index <= 
                   cloudbox_limits[5]; lon_index++)
              {
                //za_prop and aa_prop are the propagation directions
                for (Index za_prop = 0; za_prop < Nza_prop; ++ za_prop)
                  {
                    for (Index aa_prop = 0; aa_prop < Naa_prop; ++ aa_prop)
                      {
                        //za_in and aa_in are the incoming directions
                        //for which pha_mat_spt is calculated
                        for (Index za_in = 0; za_in < Nza; ++ za_in)
                          {
                            for (Index aa_in = 0; aa_in < Naa; ++ aa_in)
                              { 
                                
                                ConstMatrixView pha = pha_mat(za_in,
                                                  aa_in,
                                                  Range(joker),
                                                  Range(joker));

                                ConstVectorView i_field_in = 
                                  i_field((p_index - cloudbox_limits[0]),
                                          (lat_index - cloudbox_limits[1]),
                                          (lon_index - cloudbox_limits[2]),
                                          za_prop,
                                          aa_prop,
                                          Range(joker));
                                
                                mult(product_field( za_in,
                                                    aa_in,
                                                    Range(joker)),
                                     pha, 
                                     i_field_in);
                                
                              }//end aa_in loop
                          }//end za_in loop
                        //integration of the product of ifield_in and pha
                        //over zenith angle and azimuth angle grid. It 
                        //calls here the integration routine 
                        //AngIntegrate_trapezoid_opti
                        for (Index i = 0; i < stokes_dim; i++)
                          {
                            MatrixView product_field_mat =
                              product_field( Range(joker),
                                             Range(joker),
                                             i);
                            scat_field( (p_index - cloudbox_limits[0]),
                                        (lat_index - cloudbox_limits[1]),
                                        (lon_index - cloudbox_limits[2]),
                                        za_prop, 
                                        aa_prop,
                                        i)  =  
                              AngIntegrate_trapezoid_opti(product_field_mat,
                                                          scat_za_grid,
                                                          scat_aa_grid,
                                                          grid_stepsize);
                          }//end i loop
                      }//end aa_prop loop
                  }//end za_prop loop
              }//end lon loop
          }// end lat loop
      }// end p loop
  }// end atmosphere_dim = 3
}



//! Initialize variables for a scattering calculation
/*! 
  Variables needed in the scattering calculations are initialzed here. This 
  method has to be executed before using *ScatteringMain*.

  \param scat_p_index: Index for a position inside the cloudbox
  \param scat_lat_index: Index for a position inside the cloudbox
  \param scat_lon_index: Index for a position inside the cloudbox
  \param scat_za_index: Index for a zenith direction.
  \param scat_aa_index: Index for an azimuth direction.
  \param iteration_counter: Index for counting iterations.
  \param pha_mat : Phase matrix.
  \param pha_mat_spt: Phase matrix for a single particle type.
  \param ext_mat_spt: Extinction matrix for a single particle type.
  \param abs_vec_spt: Absorption vector for a single particle type.
  \param i_field : Radiation field
  \param scat_field: Scattered field.
  \param stokes_dim: Stokes dimension.
  \param atmosphere_dim: Atmospheric dimension.
  \param scat_za_grid: Zenith angle grid for scattering calculation.
  \param scat_aa_grid: Azimuth angle grid for scattering calculation.
  \param cloudbox_limits: Limits of cloudbox.
  \param scat_data_raw FIXME: Add documentation.
*/
void ScatteringInit(
                    Index& scat_p_index,
                    Index& scat_lat_index,
                    Index& scat_lon_index,
                    Index& scat_za_index,
                    Index& scat_aa_index,
                    Tensor4& pha_mat,
                    Tensor5& pha_mat_spt,
                    Tensor3& ext_mat_spt,
                    Matrix& abs_vec_spt,
                    Tensor6& scat_field,
                    Tensor6& i_field,
                    const Index& stokes_dim,
                    const Index& atmosphere_dim,
                    const Vector& scat_za_grid,
                    const Vector& scat_aa_grid,
                    const Index& za_grid_size,
                    const ArrayOfIndex& cloudbox_limits,
                    const ArrayOfSingleScatteringData& scat_data_raw
                    )
{
  // Check the input:
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*"); 
  scat_p_index = 0;
  scat_lat_index = 0;
  scat_lon_index = 0;
  scat_za_index = 0;
  scat_aa_index = 0;
  
  // Number of particle types:
  const Index N_pt = scat_data_raw.nelem();
  
  pha_mat.resize(za_grid_size, scat_aa_grid.nelem(), stokes_dim,
                 stokes_dim);
  
  pha_mat_spt.resize(N_pt, za_grid_size, 
                     scat_aa_grid.nelem(), stokes_dim, stokes_dim);
  
  abs_vec_spt.resize(N_pt, stokes_dim);
  abs_vec_spt = 0.;

  ext_mat_spt.resize(N_pt, stokes_dim, stokes_dim);
  ext_mat_spt = 0.;  

             
  if (atmosphere_dim == 1)
    {
      i_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                     1, 
                     1,
                     scat_za_grid.nelem(), 
                     1,
                     stokes_dim);
      
      scat_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        1, 
                        1,
                        scat_za_grid.nelem(), 
                        1,
                        stokes_dim);  
    }
  else if (atmosphere_dim == 3)
    {
      i_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                     (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                     (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                     scat_za_grid.nelem(), 
                     scat_aa_grid.nelem(),
                     stokes_dim);
      
      scat_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                        (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                        scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(),
                        stokes_dim);
    }
  else 
    {
      throw runtime_error(
                        "Scattering calculations are not possible for a 2D"
                        "atmosphere. If you want to do scattering calculations"
                        "*atmosphere_dim* has to be either 1 or 3"
                          );
    }
  
  i_field = 0.;
  scat_field = 0.;

}

//! Initialize variables for a scattering calculation
/*! 
  Variables needed in the scattering calculations are initialzed here. This 
  method has to be executed before using *ScatteringMain*.

  \param scat_p_index: Index for a position inside the cloudbox
  \param scat_lat_index: Index for a position inside the cloudbox
  \param scat_lon_index: Index for a position inside the cloudbox
  \param scat_za_index: Index for a zenith direction.
  \param scat_aa_index: Index for an azimuth direction.
  \param iteration_counter: Index for counting iterations.
  \param pha_mat : Phase matrix.
  \param pha_mat_spt: Phase matrix for a single particle type.
  \param ext_mat_spt: Extinction matrix for a single particle type.
  \param abs_vec_spt: Absorption vector for a single particle type.
  \param i_field : Radiation field
  \param scat_field: Scattered field.
  \param stokes_dim: Stokes dimension.
  \param atmosphere_dim: Atmospheric dimension.
  \param scat_za_grid: Zenith angle grid for scattering calculation.
  \param scat_aa_grid: Azimuth angle grid for scattering calculation.
  \param cloudbox_limits: Limits of cloudbox.
  \param amp_mat_raw FIXME: Add documentation.
*/
void ScatteringInitAmpMat(
                    Index& scat_p_index,
                    Index& scat_lat_index,
                    Index& scat_lon_index,
                    Index& scat_za_index,
                    Index& scat_aa_index,
                    Index& iteration_counter,
                    Tensor4& pha_mat,
                    Tensor5& pha_mat_spt,
                    Tensor3& ext_mat_spt,
                    Matrix& abs_vec_spt,
                    Tensor6& scat_field,
                    Tensor6& i_field,
                    const Index& stokes_dim,
                    const Index& atmosphere_dim,
                    const Vector& scat_za_grid,
                    const Vector& scat_aa_grid,
                    const ArrayOfIndex& cloudbox_limits,
                    const ArrayOfArrayOfTensor6& amp_mat_raw
                    )
{
  // Check the input:
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*"); 
  scat_p_index = 0;
  scat_lat_index = 0;
  scat_lon_index = 0;
  scat_za_index = 0;
  scat_aa_index = 0;
  iteration_counter = 0;

  // Number of particle types:
  const Index N_pt = amp_mat_raw.nelem();
  
  pha_mat.resize(scat_za_grid.nelem(), scat_aa_grid.nelem(), stokes_dim,
                 stokes_dim);
  pha_mat = 0.;
  
  pha_mat_spt.resize(N_pt, scat_za_grid.nelem(), 
                     scat_aa_grid.nelem(), stokes_dim, stokes_dim);
  pha_mat_spt = 0.;

  abs_vec_spt.resize(N_pt, stokes_dim);
  abs_vec_spt = 0.;

  ext_mat_spt.resize(N_pt, stokes_dim, stokes_dim);
  ext_mat_spt = 0.;  

             
  if (atmosphere_dim == 1)
    {
      i_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                     1, 
                     1,
                     scat_za_grid.nelem(), 
                     1,
                     stokes_dim);
      
      scat_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        1, 
                        1,
                        scat_za_grid.nelem(), 
                        1,
                        stokes_dim);  
    }
  else if (atmosphere_dim == 3)
    {
      i_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                     (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                     (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                     scat_za_grid.nelem(), 
                     scat_aa_grid.nelem(),
                     stokes_dim);
      
      scat_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                        (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                        scat_za_grid.nelem(), 
                        scat_aa_grid.nelem(),
                        stokes_dim);
    }
  else 
    {
      throw runtime_error(
                        "Scattering calculations are not possible for a 2D"
                        "atmosphere. If you want to do scattering calculations"
                        "*atmosphere_dim* has to be either 1 or 3"
                          );
    }
  
  i_field = 0.;
  scat_field = 0.;

}

//! Main function for the radiative transfer in cloudbox.  
/*!
  This function executes *scat_mono_agenda* for each frequency in *f_grid*.  
  
  \param f_index Frequency index
  \param f_grid Frequency grid
  \param scat_mono_agenda Agenda for monochromatic scattering calculation. 
  
  \author Claudia Emde
*/
void ScatteringMain(
                    Index& f_index,
                    const Vector& f_grid,
                    const Agenda& scat_mono_agenda
                    )
                  
{

  Index Nf = f_grid.nelem();

  // Is the frequency index valid?
  assert( f_index <= f_grid.nelem() );

  for ( f_index = 0; f_index < Nf; ++ f_index)
    {
      out1 << "----------------------------------------------\n";
      out1 << "Frequency: " << f_grid[f_index]/1e9 <<" GHz \n" ;
      out1 << "----------------------------------------------\n";
      scat_mono_agenda.execute(); 
    }
}



//! Increase iteration counter.  
/*!
  This function has to be used if you want to write the separate iteration 
  fields into differtent files using *Tensor6WriteIteration*.

  \param iteration_counter Counter for the number of iterations.
*/
void iteration_counterIncrease(Index& iteration_counter)
{
  iteration_counter += 1;
}   

 
 
//! Write iterated fields.
/*!
  This function writed intermediate resultes, the iterations of fields to xml 
  files. It can be used to check the solution method for the RTE with 
  scattering integral, which is an iterative numerical method. It is useful
  to look how the radiation field *i_field* and the scattered field
  *scat_field* behave.

  The user can give an array containing the iterations which shall be written 
  to files as a keyword to the method. E.g. if "iterations = [3, 6, 9]",
  the 3rd, 6th and 9th iterations are stored in the files
  "iteration_field_3.xml",  "iteration_field_6.xml" ...
  
  If you want to save all the iterations the array has to contain 
  just one element set to 0: "iterations = [0]". 
  
  Note: The workspace variable iteration_counter has to be set as 0 in the
  control file before using this method.
      
  WS Input/Output: 
  \param iteration_counter Counter for the iterations. 
  \param field Iterated field
  \param field_name Name of intensity field.
  \param iterations Array containing the iteration numbers to be written

  \author Claudia Emde
  \date 2002-08-26
     
*/ 
void Tensor6WriteIteration(//WS input 
                           const Index& iteration_counter,
                           //Global  Input :
                           const Tensor6& field,
                           const String& field_name,
                           //Keyword:
                           const ArrayOfIndex& iterations)
{
  // if(iteration_counter>1000) return;
   
  ostringstream os;
  os << iteration_counter;
  
  // All iterations are written to files
  if( iterations[0] == 0 )
    {
      xml_write_to_file(field_name + os.str() + ".xml", field);  
   }
  // Only the iterations given by the keyword are written to a file
  else
    {
      for (Index i = 0; i<iterations.nelem(); i++)
        {
          if (iteration_counter == iterations[i])
            xml_write_to_file(field_name + os.str() + ".xml", field);  
        }
    }
}





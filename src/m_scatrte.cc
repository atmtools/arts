/* Copyright (C) 2002-2008
   Claudia Emde <claudia.emde@dlr.de>
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
  \author Claudia Emde <claudia.emde@dlr.de>
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
#include "wsv_aux.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;
  
/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


/* Workspace method: Doxygen documentation will be auto-generated */
void DoitAngularGridsSet(// WS Output:
                  Index& doit_za_grid_size,
                  Vector& scat_aa_grid,
                  Vector& scat_za_grid,
                  // Keywords:
                  const Index& N_za_grid,
                  const Index& N_aa_grid,
                  const String& za_grid_opt_file)
{
  // Check input
  //
  // The recommended values were found by testing the accuracy and the speed of 
  // 1D DOIT calculations for different grid sizes. For 3D calculations it can 
  // be necessary to use more grid points. 
  if (N_za_grid < 16)
    throw runtime_error("N_za_grid must be greater than 15 for accurate results");
  else if (N_za_grid > 100)
    out1 << "Warning: N_za_grid is very large which means that the \n"
         << "calculation will be very slow. The recommended value is 19.\n";
  
  if (N_aa_grid < 6)
    throw runtime_error("N_aa_grid must be greater than 5 for accurate results");
  else if (N_aa_grid > 100)
    out1 << "Warning: N_aa_grid is very large which means that the \n"
         << "calculation will be very slow. The recommended value is 10.\n";
  
  // Azimuth angle grid (the same is used for the scattering integral and
  // for the radiative transfer.
  nlinspace(scat_aa_grid, 0, 360, N_aa_grid);
  
  // Zenith angle grid: 
  // Number of zenith angle grid points (only for scattering integral): 
  doit_za_grid_size = N_za_grid; 

  if( za_grid_opt_file == "" ) 
    nlinspace(scat_za_grid, 0, 180, N_za_grid);
  else
    xml_read_from_file(za_grid_opt_file, scat_za_grid);

}


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_conv_flagAbs(//WS Input and Output:
                      Index& doit_conv_flag,
                      Index& doit_iteration_counter,
                      // WS Input:
                      const Tensor6& doit_i_field,
                      const Tensor6& doit_i_field_old,
                      // Keyword:
                      const Vector& epsilon)
{
  //------------Check the input---------------------------------------
  if( doit_conv_flag != 0 )
    throw runtime_error("Convergence flag is non-zero, which means that this\n"
                        "WSM is not used correctly. *doit_conv_flagAbs* should\n"
                        "be used only in *doit_conv_test_agenda*\n");
  
  
  if (doit_iteration_counter > 100)
    {
      out1 <<"Warning in DOIT calculation:" 
           <<"Method does not converge (number of iterations \n"
           <<"is > 100). Either the cloud particle number density \n"
           <<"is too large or the numerical setup for the DOIT \n"
           <<"calculation is not correct. In case of limb \n"
           <<"simulations please make sure that you use an \n"
           <<"optimized zenith angle grid. \n"
           <<"*doit_i_field* might be wrong.\n";
      doit_conv_flag = 1;
    }
  
  const Index N_p = doit_i_field.nvitrines();
  const Index N_lat = doit_i_field.nshelves();
  const Index N_lon = doit_i_field.nbooks();
  const Index N_za = doit_i_field.npages();
  const Index N_aa = doit_i_field.nrows();
  const Index stokes_dim = doit_i_field.ncols();
  
  // Check keyword "epsilon":
  if ( epsilon.nelem() != stokes_dim )
    throw runtime_error(
                        "You have to specify limiting values for the "
                        "convergence test for each Stokes component "
                        "separately. That means that *epsilon* must "
                        "have *stokes_dim* elements!"
                        );

  // Check if doit_i_field and doit_i_field_old have the same dimensions:
  if(!is_size( doit_i_field_old, 
                  N_p, N_lat, N_lon, N_za, N_aa, stokes_dim))
    throw runtime_error("The fields (Tensor6) *doit_i_field* and \n"
                        "*doit_i_field_old* which are compared in the \n"
                        "convergence test do not have the same size.\n");
  
  //-----------End of checks-------------------------------------------------
                        

  out2 << "  Number of DOIT iteration: " << doit_iteration_counter << "\n";
  doit_iteration_counter +=1;


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
                             (doit_i_field(p_index, lat_index, lon_index, 
                                               scat_za_index, scat_aa_index, 
                                               stokes_index) -
                              doit_i_field_old(p_index, lat_index, 
                                               lon_index, scat_za_index,
                                               scat_aa_index, 
                                               stokes_index ));
                            
                          // If the absolute difference of the components
                          // is larger than the pre-defined values, return
                          // to *doit_i_fieldIterarte* and do next iteration
                          
                          if( abs(diff) > epsilon[stokes_index])
                            {
                              out1 << "difference: " << diff <<"\n";
                              return;
                            }
                          
                          
                        }// End loop stokes_dom.
                    }// End loop scat_aa_grid. 
                }// End loop scat_za_grid.
            }// End loop lon_grid. 
        }// End loop lat_grid.
    } // End p_grid.
  
  // Convergence test has been successful, doit_conv_flag can be set to 1.
  doit_conv_flag = 1;
  out2 << "  Number of DOIT-iterations: " << doit_iteration_counter << "\n";
}
      

/* Workspace method: Doxygen documentation will be auto-generated */
void doit_conv_flagAbsBT(//WS Input and Output:
                         Index& doit_conv_flag,
                         Index& doit_iteration_counter,
                         // WS Input:
                         const Tensor6& doit_i_field,
                         const Tensor6& doit_i_field_old,
                         const Vector& f_grid,
                         const Index& f_index, 
                         // Keyword:
                         const Vector& epsilon)
{
   //------------Check the input---------------------------------------
  
  if( doit_conv_flag != 0 )
    throw runtime_error("Convergence flag is non-zero, which means that this \n"
                        "WSM is not used correctly. *doit_conv_flagAbs* should\n"
                        "be used only in *doit_conv_test_agenda*\n");
  
  if (doit_iteration_counter > 100)
    {
      out1 <<"Warning in DOIT calculation at frequency" << f_grid[f_index] 
           << "GHz: \n"
           <<"Method does not converge (number of iterations \n"
           <<"is > 100). Either the cloud particle number density \n"
           <<"is too large or the numerical setup for the DOIT \n"
           <<"calculation is not correct. In case of limb \n"
           <<"simulations please make sure that you use an \n"
           <<"optimized zenith angle grid. \n"
           <<"*doit_i_field* might be wrong.\n";
      doit_conv_flag = 1;
    }
  
  const Index N_p = doit_i_field.nvitrines();
  const Index N_lat = doit_i_field.nshelves();
  const Index N_lon = doit_i_field.nbooks();
  const Index N_za = doit_i_field.npages();
  const Index N_aa = doit_i_field.nrows();
  const Index stokes_dim = doit_i_field.ncols();
  
  // Check keyword "epsilon":
  if ( epsilon.nelem() != stokes_dim )
    throw runtime_error(
                        "You have to specify limiting values for the "
                        "convergence test for each Stokes component "
                        "separately. That means that *epsilon* must "
                        "have *stokes_dim* elements!"
                        );
  
  // Check if doit_i_field and doit_i_field_old have the same dimensions:
  if(!is_size( doit_i_field_old, 
               N_p, N_lat, N_lon, N_za, N_aa, stokes_dim))
    throw runtime_error("The fields (Tensor6) *doit_i_field* and \n"
                        "*doit_i_field_old* which are compared in the \n"
                        "convergence test do not have the same size.\n");

  // Frequency grid
  //
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );

  // Is the frequency index valid?
  if ( f_index >= f_grid.nelem() )
    throw runtime_error("*f_index* is greater than number of elements in the\n"
                        "frequency grid.\n");
  
  //-----------End of checks--------------------------------

  out2 << "  Number of DOIT iteration: " << doit_iteration_counter << "\n";
  doit_iteration_counter +=1;

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
                            abs( doit_i_field(p_index, lat_index, lon_index, 
                                          scat_za_index, scat_aa_index, 
                                          stokes_index) -
                                  doit_i_field_old(p_index, lat_index, 
                                              lon_index, scat_za_index,
                                              scat_aa_index, 
                                              stokes_index ));
                          
                          // If the absolute difference of the components
                          // is larger than the pre-defined values, return
                          // to *doit_i_fieldIterarte* and do next iteration
                          Numeric diff_bt = invrayjean(diff, f_grid[f_index]);
                          if( diff_bt > epsilon[stokes_index])
                            {
                              out1 << "BT difference: " << diff_bt <<"\n";
                              return;
                            }
                        }// End loop stokes_dom.
                    }// End loop scat_aa_grid. 
                }// End loop scat_za_grid.
            }// End loop lon_grid. 
        }// End loop lat_grid.
    } // End p_grid.
  
  // Convergence test has been successful, doit_conv_flag can be set to 1.
  doit_conv_flag = 1;
  out1 << "Number of DOIT-iterations:" << doit_iteration_counter << "\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_conv_flagLsq(//WS Output:
                      Index& doit_conv_flag,
                      Index& doit_iteration_counter,
                      // WS Input:
                      const Tensor6& doit_i_field,
                      const Tensor6& doit_i_field_old,
                      const Vector& f_grid,
                      const Index& f_index,
                      // Keyword:
                      const Vector& epsilon)
{
  //------------Check the input---------------------------------------
  
  if( doit_conv_flag != 0 )
    throw runtime_error("Convergence flag is non-zero, which means that this \n"
                        "WSM is not used correctly. *doit_conv_flagAbs* should\n"
                        "be used only in *doit_conv_test_agenda*\n");
  
  
  if (doit_iteration_counter > 100)
    throw runtime_error("Error in DOIT calculation: \n"
                        "Method does not converge (number of iterations \n"
                        "is > 100). Either the cloud particle number density \n"
                        "is too large or the numerical setup for the DOIT \n"
                        "calculation is not correct. In case of limb \n"
                        "simulations please make sure that you use an \n"
                        "optimized zenith angle grid. \n");
  
  const Index N_p = doit_i_field.nvitrines();
  const Index N_lat = doit_i_field.nshelves();
  const Index N_lon = doit_i_field.nbooks();
  const Index N_za = doit_i_field.npages();
  const Index N_aa = doit_i_field.nrows();
  const Index stokes_dim = doit_i_field.ncols();
  
  // Check keyword "epsilon":
  if ( epsilon.nelem() != stokes_dim )
    throw runtime_error(
                        "You have to specify limiting values for the "
                        "convergence test for each Stokes component "
                        "separately. That means that *epsilon* must "
                        "have *stokes_dim* elements!"
                        );

  // Check if doit_i_field and doit_i_field_old have the same dimensions:
  if(!is_size( doit_i_field_old, 
               N_p, N_lat, N_lon, N_za, N_aa, stokes_dim))
    throw runtime_error("The fields (Tensor6) *doit_i_field* and \n"
                        "*doit_i_field_old* which are compared in the \n"
                        "convergence test do not have the same size.\n");
  
  // Frequency grid
  //
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );

  // Is the frequency index valid?
  if ( f_index >= f_grid.nelem() )
    throw runtime_error("*f_index* is greater than number of elements in the\n"
                        "frequency grid.\n");
  
  //-----------End of checks--------------------------------

 
  out2 << "  Number of DOIT iteration: " << doit_iteration_counter << "\n";
  doit_iteration_counter +=1;                   
  
  Vector lqs(4, 0.);
  
  // Will be set to zero if convergence not fullfilled
  doit_conv_flag = 1;
  for (Index i = 0; i < epsilon.nelem(); i ++)
    {
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
                          lqs[i] 
                            += pow(
                                   doit_i_field(p_index, lat_index, 
                                                lon_index, 
                                           scat_za_index, scat_aa_index, i) -
                                   doit_i_field_old(p_index, lat_index, 
                                               lon_index, scat_za_index,
                                               scat_aa_index, i) 
                                   , 2);
                        }// End loop scat_aa_grid. 
                    }// End loop scat_za_grid.
                }// End loop lon_grid. 
            }// End loop lat_grid.
        } // End p_grid.
      
      lqs[i] = sqrt(lqs[i]);
      lqs[i] /= (N_p*N_lat*N_lon*N_za*N_aa);

      // Convert difference to Rayleigh Jeans BT
      lqs[i] = invrayjean(lqs[i], f_grid[f_index]);
      
      if (lqs[i] >= epsilon[i] )
        doit_conv_flag = 0;
    }
  // end loop stokes_index
  out1 << "lqs [I]: " << lqs[0] << "\n";
  
  if (doit_conv_flag == 1)
    {
      // Convergence test has been successful,
      out1 << "Number of DOIT-iterations: " << doit_iteration_counter 
           << "\n";
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_i_fieldIterate(
               // WS Input and Output:
               Tensor6& doit_i_field,
               // WS Input:  
               const Agenda& doit_scat_field_agenda,
               const Agenda& doit_rte_agenda,
               const Agenda& doit_conv_test_agenda
               )
{
  //---------------Check input---------------------------------
  chk_not_empty( "doit_scat_field_agenda", doit_scat_field_agenda);
  chk_not_empty( "doit_rte_agenda", doit_rte_agenda);
  chk_not_empty( "doit_conv_test_agenda", doit_conv_test_agenda);
  
  //doit_i_field can not be checked here, because there is no way
  //to find out the size without including a lot more interface 
  //variables
  //-----------End of checks-------------------------------------- 

  Tensor6 doit_i_field_old_local;
  Index doit_conv_flag_local;
  Index doit_iteration_counter_local;

  // Resize and initialize doit_scat_field,
  // which  has the same dimensions as doit_i_field
  Tensor6 doit_scat_field_local
    (doit_i_field.nvitrines(), doit_i_field.nshelves(),
     doit_i_field.nbooks(), doit_i_field.npages(),
     doit_i_field.nrows(), doit_i_field.ncols(), 0.);

  doit_conv_flag_local = 0;
  doit_iteration_counter_local = 0;

  while(doit_conv_flag_local == 0) {
    
    // 1. Copy doit_i_field to doit_i_field_old.
    doit_i_field_old_local = doit_i_field;
    
    // 2.Calculate scattered field vector for all points in the cloudbox.
    
    // Calculate the scattered field.
    out2 << "  Execute doit_scat_field_agenda. \n";
    doit_scat_field_agendaExecute(doit_scat_field_local,
                                  doit_i_field,
                                  doit_scat_field_agenda,
                                  true);
    
    // Update doit_i_field.
    out2 << "  Execute doit_rte_agenda. \n";
    doit_rte_agendaExecute(doit_i_field, doit_scat_field_local, 
                                doit_rte_agenda, true);

    //Convergence test.
    doit_conv_test_agendaExecute(doit_conv_flag_local,
                                 doit_iteration_counter_local,
                                 doit_i_field,
                                 doit_i_field_old_local,
                                 doit_conv_test_agenda,
                                 true);

  }//end of while loop, convergence is reached.
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
doit_i_fieldUpdate1D(// WS Input and Output:
                   Tensor6& doit_i_field,
                   // WS Input:
                   const Tensor6& doit_i_field_old,
                   const Tensor6& doit_scat_field,
                   const ArrayOfIndex& cloudbox_limits,
                   // Calculate scalar gas absorption:
                   const Agenda& abs_scalar_gas_agenda,
                   const Tensor4& vmr_field,
                   // Optical properties for single particle type:
                   const Agenda& spt_calc_agenda,
                   const Vector& scat_za_grid,
                   const Tensor4& pnd_field,
                   // Optical properties for gases and particles:
                   const Agenda& opt_prop_part_agenda,
                   const Agenda& opt_prop_gas_agenda,
                   // Propagation path calculation:
                   const Agenda& ppath_step_agenda,
                   const Vector& p_grid,
                   const Tensor3& z_field,
                   const Matrix& r_geoid,
                   const Matrix& z_surface,
                   // Calculate thermal emission:
                   const Tensor3& t_field,
                   const Vector& f_grid,
                   const Index& f_index,
                   const Agenda&, //surface_prop_agenda,
                   const Index& doit_za_interp
                   )
{
  
  out2 << "  doit_i_fieldUpdateSeq1D: Radiative transfer calculation in cloudbox\n";
  out2 << "  ------------------------------------------------------------- \n";
  
  // ---------- Check the input ----------------------------------------
  
  // Agendas
  chk_not_empty( "spt_calc_agenda", spt_calc_agenda);
  chk_not_empty( "opt_prop_part_agenda", opt_prop_part_agenda);
  chk_not_empty( "opt_prop_gas_agenda", opt_prop_gas_agenda);
  chk_not_empty( "ppath_step_agenda", ppath_step_agenda);
  
  if (cloudbox_limits.nelem() != 2)
    throw runtime_error(
                        "The cloudbox dimension is not 1D! \n"
                        "Do you really want to do a 1D calculation? \n"
                        "If not, use *doit_i_fieldUpdateSeq3D*.\n"
                        );
  
  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  
  if (scat_za_grid[0] != 0. || scat_za_grid[N_scat_za-1] != 180.)
    throw runtime_error("The range of *scat_za_grid* must [0 180].");
  
  if( p_grid.nelem() < 2 )
    throw runtime_error( "The length of *p_grid* must be >= 2." );
  chk_if_decreasing( "p_grid", p_grid );

  chk_size("z_field", z_field, p_grid.nelem(), 1, 1);
  chk_size("t_field", t_field, p_grid.nelem(), 1, 1);
  
  const Vector dummy(1,0.);
  chk_atm_surface( "r_geoid", r_geoid, 1, dummy, 
                   dummy);
  
  // Frequency grid
  //
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  
  // Is the frequency index valid?
  if ( f_index >= f_grid.nelem() )
    throw runtime_error("*f_index* is greater than number of elements in the\n"
                        "frequency grid.\n");
  
  if( !(doit_za_interp == 0  ||  doit_za_interp == 1 ) )
    throw runtime_error( "Interpolation method is not defined. Use \n"
                         "*doit_za_interpSet*.\n");
  
  const Index stokes_dim = doit_scat_field.ncols();
  assert(stokes_dim > 0 || stokes_dim < 5);


  // These variables are calculated internally, so assertions should be o.k.
  assert( is_size( doit_i_field, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1, 1, 1, 
                   N_scat_za, 1, stokes_dim));
  
  assert( is_size( doit_i_field_old, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1, 1, 1, 
                   scat_za_grid.nelem(), 1, stokes_dim));
  
  assert( is_size( doit_scat_field, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1, 1, 1, 
                   N_scat_za, 1, stokes_dim));
  
  // FIXME: Check *vmr_field* 
  
  // -------------- End of checks --------------------------------------
 
  
  //=======================================================================
  // Calculate scattering coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are calulated for interpolated VMRs,
  // temperature and pressure.
      
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, stokes_dim, 0.);
  Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, 0.);
 
  //Only dummy variable:
  Index scat_aa_index_local = 0; 

  //Loop over all directions, defined by scat_za_grid 
  for( Index scat_za_index_local = 0; scat_za_index_local < N_scat_za; 
       scat_za_index_local ++)
    {
      // This function has to be called inside the angular loop, as
      // spt_calc_agenda takes *scat_za_index_local* and *scat_aa_index* 
      // from the workspace.
      // *scat_p_index* is needed for communication with agenda 
      // *opt_prop_part_agenda*.
      cloud_fieldsCalc(ext_mat_field, abs_vec_field,
                       spt_calc_agenda, 
                       opt_prop_part_agenda, scat_za_index_local, 
                       scat_aa_index_local,
                       cloudbox_limits, t_field, pnd_field);
      
      //======================================================================
      // Radiative transfer inside the cloudbox
      //=====================================================================
      
      for(Index p_index = cloudbox_limits[0]; p_index
            <= cloudbox_limits[1]; p_index ++)
        {
          cloud_ppath_update1D_noseq(doit_i_field, 
                                     p_index, scat_za_index_local, 
                                     scat_za_grid,
                                     cloudbox_limits, doit_i_field_old, 
                                     doit_scat_field,
                                     abs_scalar_gas_agenda, vmr_field,
                                     opt_prop_gas_agenda, ppath_step_agenda,
                                     p_grid,  z_field, r_geoid, z_surface,
                                     t_field, f_grid, f_index, ext_mat_field, 
                                     abs_vec_field,
                                     doit_za_interp);
        }
    }// Closes loop over scat_za_grid.
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
doit_i_fieldUpdateSeq1D(// WS Input and Output:
                   Tensor6& doit_i_field,
                   // WS Input:
                   const Tensor6& doit_scat_field,
                   const ArrayOfIndex& cloudbox_limits,
                   // Calculate scalar gas absorption:
                   const Agenda& abs_scalar_gas_agenda,
                   const Tensor4& vmr_field,
                   // Optical properties for single particle type:
                   const Agenda& spt_calc_agenda,
                   const Vector& scat_za_grid,
                   const Tensor4& pnd_field, 
                   // Optical properties for gases and particles:
                   const Agenda& opt_prop_part_agenda,
                   const Agenda& opt_prop_gas_agenda,
                   // Propagation path calculation:
                   const Agenda& ppath_step_agenda,
                   const Vector& p_grid,
                   const Tensor3& z_field,
                   const Matrix& r_geoid,
                   const Matrix& z_surface,
                   // Calculate thermal emission:
                   const Tensor3& t_field,
                   const Vector& f_grid,
                   const Index& f_index,
                   const Agenda& surface_prop_agenda, //STR
                   const Index& doit_za_interp
                   )
{
  
  out2<<"  doit_i_fieldUpdateSeq1D: Radiative transfer calculation in cloudbox\n";
  out2 << "  ------------------------------------------------------------- \n";

 // ---------- Check the input ----------------------------------------
  
  // Agendas
  chk_not_empty( "abs_scalar_gas_agenda", abs_scalar_gas_agenda);
  chk_not_empty( "spt_calc_agenda", spt_calc_agenda);
  chk_not_empty( "opt_prop_part_agenda", opt_prop_part_agenda);
  chk_not_empty( "opt_prop_gas_agenda", opt_prop_gas_agenda);
  chk_not_empty( "ppath_step_agenda", ppath_step_agenda);
  
  if (cloudbox_limits.nelem() != 2)
    throw runtime_error(
                        "The cloudbox dimension is not 1D! \n"
                        "Do you really want to do a 1D calculation? \n"
                        "For 3D use *doit_i_fieldUpdateSeq3D*.\n"
                        );
   
  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  
  if (scat_za_grid[0] != 0. || scat_za_grid[N_scat_za-1] != 180.)
    throw runtime_error("The range of *scat_za_grid* must [0 180].");
  
  if( p_grid.nelem() < 2 )
    throw runtime_error( "The length of *p_grid* must be >= 2." );
  chk_if_decreasing( "p_grid", p_grid );

  chk_size("z_field", z_field, p_grid.nelem(), 1, 1);
  chk_size("t_field", t_field, p_grid.nelem(), 1, 1);
  
  const Vector dummy(1,0.);
  chk_atm_surface( "r_geoid", r_geoid, 1, dummy, 
                  dummy);
  
  // Frequency grid
  //
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  
  // Is the frequency index valid?
  if ( f_index >= f_grid.nelem() )
    throw runtime_error("*f_index* is greater than number of elements in the\n"
                        "frequency grid.\n");
  
  if( !(doit_za_interp == 0  ||  doit_za_interp == 1 ) )
    throw runtime_error( "Interpolation method is not defined. Use \n"
                         "*doit_za_interpSet*.\n");
  
  const Index stokes_dim = doit_scat_field.ncols();
  assert(stokes_dim > 0 || stokes_dim < 5); 
  
  
  // These variables are calculated internally, so assertions should be o.k.
  assert( is_size( doit_i_field, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1, 1, 1, 
                   N_scat_za, 1, stokes_dim));
  
  assert( is_size( doit_scat_field, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1, 1, 1, 
                   N_scat_za, 1, stokes_dim));
  
  // FIXME: Check *vmr_field* 
  
  // -------------- End of checks --------------------------------------
  
     
  //=======================================================================
  // Calculate scattering coefficients for all positions in the cloudbox 
  //=======================================================================
  out3 << "Calculate single particle properties \n";

  // At this place only the particle properties are calculated. Gaseous
  // absorption is calculated inside the radiative transfer part. Inter-
  // polating absorption coefficients for gaseous species gives very bad
  // results, so they are calulated for interpolated VMRs,
  // temperature and pressure.
      
  // To use special interpolation functions for atmospheric fields we 
  // use ext_mat_field and abs_vec_field:
  Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, stokes_dim, 0.);
  Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                        stokes_dim, 0.);
  
     
  // If theta is between 90° and the limiting value, the intersection point
  // is exactly at the same level as the starting point (cp. AUG)
  Numeric theta_lim = 180. - asin((r_geoid(0,0)+
                                   z_field(cloudbox_limits[0],0,0))/
                                  (r_geoid(0,0)+
                                   z_field(cloudbox_limits[1],0,0)))*RAD2DEG;
  //Only dummy variables:
  Index scat_aa_index_local = 0;
  
  //Loop over all directions, defined by scat_za_grid 
  for(Index scat_za_index_local = 0; scat_za_index_local < N_scat_za; 
      scat_za_index_local ++)
    {
      // This function has to be called inside the angular loop, as
      // spt_calc_agenda takes *scat_za_index* and *scat_aa_index* 
      // from the workspace.
      // *scat_p_index* is needed for communication with agenda 
      // *opt_prop_part_agenda*.
      cloud_fieldsCalc(ext_mat_field, abs_vec_field, 
                       spt_calc_agenda, opt_prop_part_agenda, 
                       scat_za_index_local, scat_aa_index_local, 
                       cloudbox_limits, t_field, pnd_field);
  
      
     
      //  xml_write_to_file("ext_mat_field.xml", ext_mat_field );
      // xml_write_to_file("abs_vec_field.xml", ext_mat_field );
      
      //======================================================================
      // Radiative transfer inside the cloudbox
      //=====================================================================
   
      
      // Sequential update for uplooking angles
      if ( scat_za_grid[scat_za_index_local] <= 90.) 
        {
          // Loop over all positions inside the cloud box defined by the 
          // cloudbox_limits excluding the upper boundary. For uplooking
          // directions, we start from cloudbox_limits[1]-1 and go down
          // to cloudbox_limits[0] to do a sequential update of the
          // radiation field
          for(Index p_index = cloudbox_limits[1]-1; p_index
                >= cloudbox_limits[0]; p_index --)
            {
              cloud_ppath_update1D(doit_i_field, 
                                   p_index, scat_za_index_local, scat_za_grid,
                                   cloudbox_limits, doit_scat_field,
                                   abs_scalar_gas_agenda, vmr_field,
                                   opt_prop_gas_agenda, ppath_step_agenda,
                                   p_grid,  z_field, r_geoid, z_surface,
                                   t_field, f_grid, f_index, ext_mat_field,
                                   abs_vec_field, 
                                   surface_prop_agenda, doit_za_interp); 
            }
        }
      else if ( scat_za_grid[scat_za_index_local] >= theta_lim) 
        {
          //
          // Sequential updating for downlooking angles
          //
          for(Index p_index = cloudbox_limits[0]+1; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
              cloud_ppath_update1D(doit_i_field,  
                                   p_index, scat_za_index_local, scat_za_grid,
                                   cloudbox_limits, doit_scat_field,
                                   abs_scalar_gas_agenda, vmr_field,
                                   opt_prop_gas_agenda, ppath_step_agenda,
                                   p_grid,  z_field, r_geoid, z_surface,
                                   t_field, f_grid, f_index, ext_mat_field, 
                                   abs_vec_field, 
                                   surface_prop_agenda, doit_za_interp); 
            }// Close loop over p_grid (inside cloudbox).
        } // end if downlooking.
      
      //
      // Limb looking:
      // We have to include a special case here, as we may miss the endpoints
      // when the intersection point is at the same level as the aactual point.
      // To be save we loop over the full cloudbox. Inside the function 
      // cloud_ppath_update1D it is checked whether the intersection point is 
      // inside the cloudbox or not.
      else if (  scat_za_grid[scat_za_index_local] > 90. &&
                 scat_za_grid[scat_za_index_local] < theta_lim ) 
        {
          for(Index p_index = cloudbox_limits[0]; p_index
                <= cloudbox_limits[1]; p_index ++)
            {
              // For this case the cloudbox goes down to the surface and we
              // look downwards. These cases are outside the cloudbox and 
              // not needed. Switch is included here, as ppath_step_agenda 
              // gives an error for such cases.
              if (!(p_index == 0 && scat_za_grid[scat_za_index_local] > 90.))
                {
                  cloud_ppath_update1D(doit_i_field,  
                                       p_index, scat_za_index_local,
                                       scat_za_grid,
                                       cloudbox_limits, doit_scat_field,
                                       abs_scalar_gas_agenda, vmr_field,
                                       opt_prop_gas_agenda, ppath_step_agenda,
                                       p_grid,  z_field, r_geoid, z_surface,
                                       t_field, f_grid, f_index, ext_mat_field, 
                                       abs_vec_field,
                                       surface_prop_agenda, doit_za_interp);
                }
            }
        } 
    }// Closes loop over scat_za_grid.
} // End of the function.

                         
/* Workspace method: Doxygen documentation will be auto-generated */
void
doit_i_fieldUpdateSeq3D(// WS Output and Input:
                        Tensor6& doit_i_field,
                        // WS Input:
                        const Tensor6& doit_scat_field,
                        const ArrayOfIndex& cloudbox_limits,
                        // Calculate scalar gas absorption:
                        const Agenda& abs_scalar_gas_agenda,
                        const Tensor4& vmr_field,
                        // Optical properties for single particle type:
                        const Agenda& spt_calc_agenda,
                        const Vector& scat_za_grid,
                        const Vector& scat_aa_grid,
                        const Tensor4& pnd_field,
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
                        const Matrix& z_surface,
                        // Calculate thermal emission:
                        const Tensor3& t_field,
                        const Vector& f_grid,
                        const Index& f_index,
                        const Index& doit_za_interp
                     )
{
  out2<<"  doit_i_fieldUpdateSeq3D: Radiative transfer calculatiuon in cloudbox.\n";
  out2 << "  ------------------------------------------------------------- \n";
  
  // ---------- Check the input ----------------------------------------

   // Agendas
  chk_not_empty( "abs_scalar_gas_agenda", abs_scalar_gas_agenda);
  chk_not_empty( "spt_calc_agenda", spt_calc_agenda);
  chk_not_empty( "opt_prop_part_agenda", opt_prop_part_agenda);
  chk_not_empty( "opt_prop_gas_agenda", opt_prop_gas_agenda);
  chk_not_empty( "ppath_step_agenda", ppath_step_agenda);
  
  if (cloudbox_limits.nelem() != 6)
    throw runtime_error(
                        "The cloudbox dimension is not 3D! \n"
                        "Do you really want to do a 3D calculation? \n"
                        "For 1D use *doit_i_fieldUpdateSeq1D*.\n"
                        );

  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  
  if (scat_za_grid[0] != 0. || scat_za_grid[N_scat_za-1] != 180.)
    throw runtime_error("The range of *scat_za_grid* must [0 180].");
  
  // Number of azimuth angles.
  const Index N_scat_aa = scat_aa_grid.nelem();
  
  if (scat_aa_grid[0] != 0. || scat_aa_grid[N_scat_aa-1] != 360.)
    throw runtime_error("The range of *scat_za_grid* must [0 360]."); 

  // Check atmospheric grids
  chk_atm_grids(3, p_grid, lat_grid, lon_grid);

  // Check atmospheric fields
  chk_size("z_field", z_field, p_grid.nelem(), lat_grid.nelem(), 
           lon_grid.nelem());
  chk_size("t_field", t_field, p_grid.nelem(), lat_grid.nelem(), 
           lon_grid.nelem());

  chk_atm_surface( "r_geoid", r_geoid, 3, lat_grid, 
                   lon_grid);
  
  // Frequency grid
  //
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );
  
  // Is the frequency index valid?
  if ( f_index >= f_grid.nelem() )
    throw runtime_error("*f_index* is greater than number of elements in the\n"
                        "frequency grid.\n");
  
  if( !(doit_za_interp == 0  ||  doit_za_interp == 1 ) )
    throw runtime_error( "Interpolation method is not defined. Use \n"
                         "*doit_za_interpSet*.\n");
  
  const Index stokes_dim = doit_scat_field.ncols();
  assert(stokes_dim > 0 || stokes_dim < 5); 
  
  // These variables are calculated internally, so assertions should be o.k.
  assert( is_size( doit_i_field, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                   (cloudbox_limits[3] - cloudbox_limits[2]) + 1, 
                   (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
                   N_scat_za,
                   N_scat_aa,
                   stokes_dim));

  assert( is_size( doit_scat_field, 
                   (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                   (cloudbox_limits[3] - cloudbox_limits[2]) + 1, 
                   (cloudbox_limits[5] - cloudbox_limits[4]) + 1,
                   N_scat_za,
                   N_scat_aa,
                   stokes_dim));
  
  // FIXME: Check *vmr_field* 
  
  // ---------- End of checks ------------------------------------------

  
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
  for(Index scat_za_index = 0; scat_za_index < N_scat_za; scat_za_index ++)
    {
      //Loop over azimuth directions (scat_aa_grid). First and last point in 
      // azimuth angle grid are euqal. Start with second element.
      for(Index scat_aa_index = 1; scat_aa_index < N_scat_aa; scat_aa_index ++)
        {
         //==================================================================
          // Radiative transfer inside the cloudbox
          //==================================================================

          // This function has to be called inside the angular loop, as
          // it spt_calc_agenda takes *scat_za_index* and *scat_aa_index* 
          // from the workspace.
          cloud_fieldsCalc(ext_mat_field, abs_vec_field, 
                           spt_calc_agenda, 
                           opt_prop_part_agenda, scat_za_index, 
                           scat_aa_index, cloudbox_limits, t_field, 
                           pnd_field);
          

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
                          cloud_ppath_update3D(doit_i_field, 
                                               p_index, lat_index, 
                                               lon_index, scat_za_index, 
                                               scat_aa_index, scat_za_grid, 
                                               scat_aa_grid, cloudbox_limits, 
                                               doit_scat_field, 
                                               abs_scalar_gas_agenda,
                                               vmr_field, 
                                               opt_prop_gas_agenda,
                                               ppath_step_agenda, p_grid, 
                                               lat_grid, lon_grid, z_field, 
                                               r_geoid, z_surface, t_field,
                                               f_grid, f_index,
                                               ext_mat_field, abs_vec_field,
                                               doit_za_interp);
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
                          cloud_ppath_update3D(doit_i_field, 
                                               p_index, lat_index, 
                                               lon_index, scat_za_index, 
                                               scat_aa_index, scat_za_grid, 
                                               scat_aa_grid, cloudbox_limits, 
                                               doit_scat_field, 
                                               abs_scalar_gas_agenda,
                                               vmr_field, 
                                               opt_prop_gas_agenda,
                                               ppath_step_agenda, p_grid, 
                                               lat_grid, lon_grid, z_field, 
                                               r_geoid, z_surface, t_field,
                                               f_grid, f_index,
                                               ext_mat_field, abs_vec_field,
                                               doit_za_interp);
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
                  // For this case the cloudbox goes down to the surface an we
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
                              cloud_ppath_update3D(doit_i_field, 
                                                   p_index, 
                                                   lat_index, 
                                                   lon_index, scat_za_index, 
                                                   scat_aa_index,
                                                   scat_za_grid, 
                                                   scat_aa_grid,
                                                   cloudbox_limits, 
                                                   doit_scat_field, 
                                                   abs_scalar_gas_agenda,
                                                   vmr_field, 
                                                   opt_prop_gas_agenda,
                                                   ppath_step_agenda, p_grid, 
                                                   lat_grid, lon_grid,
                                                   z_field, 
                                                   r_geoid, z_surface,
                                                   t_field, f_grid,
                                                   f_index,
                                                   ext_mat_field,
                                                   abs_vec_field,
                                                   doit_za_interp); 
                            }
                        }
                    }
                }
            }
        } //  Closes loop over aa_grid.
    }// Closes loop over scat_za_grid.

  doit_i_field(joker, joker, joker, joker, 0, joker) = 
    doit_i_field(joker, joker, joker, joker, N_scat_aa-1, joker);
  
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
doit_i_fieldUpdateSeq1DPP(// WS Output:
                Tensor6& doit_i_field,
                // spt_calc_agenda:
                Index& scat_za_index ,
                // WS Input:
                const Tensor6& doit_scat_field,
                const ArrayOfIndex& cloudbox_limits,
                // Calculate scalar gas absorption:
                const Agenda& abs_scalar_gas_agenda,
                const Tensor4& vmr_field,
                // Optical properties for single particle type:
                const Agenda& spt_calc_agenda,
                const Vector& scat_za_grid,
                const Tensor4& pnd_field,
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

  out2 << "  doit_i_fieldUpdateSeq1DPP: Radiative transfer calculation in cloudbox.\n";
  out2 << "  --------------------------------------------------------------------- \n";
  
  const Index stokes_dim = doit_scat_field.ncols();
  //  const Index atmosphere_dim = 1;

  //Check the input
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");
  
  assert( is_size( doit_i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) + 1,
                      1, 
                      1,
                      scat_za_grid.nelem(), 
                      1,
                      stokes_dim));

  assert( is_size( doit_scat_field, 
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
      Index scat_aa_index = 0;
      
      cloud_fieldsCalc(ext_mat_field, abs_vec_field, 
                       spt_calc_agenda, 
                       opt_prop_part_agenda, scat_za_index, scat_aa_index, 
                       cloudbox_limits, t_field, 
                       pnd_field);

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
              cloud_ppath_update1D_planeparallel(doit_i_field, 
                                                 p_index, scat_za_index,
                                                 scat_za_grid,
                                                 cloudbox_limits,
                                                 doit_scat_field,
                                                 abs_scalar_gas_agenda,
                                                 vmr_field,
                                                 opt_prop_gas_agenda,
                                                 ppath_step_agenda,
                                                 p_grid, z_field, r_geoid,
                                                 t_field, 
                                                 f_grid, f_index,
                                                 ext_mat_field,
                                                 abs_vec_field); 
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
              cloud_ppath_update1D_planeparallel(doit_i_field,  
                                                 p_index, scat_za_index,
                                                 scat_za_grid,
                                                 cloudbox_limits,
                                                 doit_scat_field,
                                                 abs_scalar_gas_agenda,
                                                 vmr_field,
                                                 opt_prop_gas_agenda,
                                                 ppath_step_agenda,
                                                 p_grid, z_field, r_geoid,
                                                 t_field, 
                                                 f_grid, f_index,
                                                 ext_mat_field, 
                                                 abs_vec_field);  
            }// Close loop over p_grid (inside cloudbox).
        } // end if downlooking.
      
    }// Closes loop over scat_za_grid.
}


/* Workspace method: Doxygen documentation will be auto-generated */
void DoitInit(
              //WS Output
              Index& scat_p_index,
              Index& scat_lat_index,
              Index& scat_lon_index,
              Index& scat_za_index,
              Index& scat_aa_index,
              Tensor6& doit_scat_field,
              Tensor6& doit_i_field,
              Index& doit_za_interp,
              Index& doit_is_initialized,
              // WS Input
              const Index& stokes_dim,
              const Index& atmosphere_dim,
              const Vector& scat_za_grid,
              const Vector& scat_aa_grid,
              const Index& doit_za_grid_size,
              const ArrayOfIndex& cloudbox_limits,
              const ArrayOfSingleScatteringData& scat_data_raw
              )
{
  // -------------- Check the input ------------------------------
  
  if (stokes_dim < 0 || stokes_dim > 4)
    throw runtime_error(
                        "The dimension of stokes vector must be"
                        "1,2,3, or 4");

  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );
  
  // Number of zenith angles.
  const Index N_scat_za = scat_za_grid.nelem();
  
  if (scat_za_grid[0] != 0. || scat_za_grid[N_scat_za-1] != 180.)
    throw runtime_error("The range of *scat_za_grid* must [0 180].");
  
  // Number of azimuth angles.
  const Index N_scat_aa = scat_aa_grid.nelem();
  
  if (scat_aa_grid[0] != 0. || scat_aa_grid[N_scat_aa-1] != 360.)
    throw runtime_error("The range of *scat_za_grid* must [0 360]."); 
  
  if (doit_za_grid_size < 16)
    throw runtime_error(
     "*doit_za_grid_size* must be greater than 15 for accurate results");
  else if (doit_za_grid_size > 100)
    out1 << "Warning: doit_za_grid_size is very large which means that the \n"
         << "calculation will be very slow. The recommended value is 19.\n";
  
  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");

  if (scat_data_raw.nelem() == 0)
    throw runtime_error(
                         "No scattering data files have been added.\n"
                         "Please use the WSM *ParticleTypeAdd* or \n"
                         "*ParticleTypeAddAll* to define the cloud \n"
                         "properties for the scattering calculation.\n"
                         );

  //------------- end of checks ---------------------------------------
  
  
  // Initialize indices

  scat_p_index = 0;
  scat_lat_index = 0;
  scat_lon_index = 0;
  scat_za_index = 0;
  scat_aa_index = 0;
  
  
  // Resize and initialize radiation field in the cloudbox
  if (atmosphere_dim == 1)
    {
      doit_i_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                     1, 
                     1,
                     scat_za_grid.nelem(), 
                     1,
                     stokes_dim);
      
      doit_scat_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        1, 
                        1,
                        scat_za_grid.nelem(), 
                        1,
                        stokes_dim);
    }
  else if (atmosphere_dim == 3)
    {
      doit_i_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
                     (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                     (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                     scat_za_grid.nelem(), 
                     scat_aa_grid.nelem(),
                     stokes_dim);
      
      doit_scat_field.resize((cloudbox_limits[1] - cloudbox_limits[0]) +1,
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
  
  doit_i_field = 0.;
  doit_scat_field = 0.;

                                    
                                    
  // Default interpolation method is "linear"
  if (doit_za_interp != 1)
    doit_za_interp = 0;
  
  doit_is_initialized = 1;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void DoitWriteIterationFields(//WS input 
                           const Index& doit_iteration_counter,
                           const Tensor6& doit_i_field,
                           //Keyword:
                           const ArrayOfIndex& iterations)
{
  // Checks of doit_i_field have been done elsewhere, e.g. in
  // scat_fieldCalc(Limb).

  // If the number of iterations is less than a number specified in the 
  // keyword *iterations*, this number will be ignored.

  ostringstream os;
  os << doit_iteration_counter;
  
  // All iterations are written to files
  if( iterations[0] == 0 )
    {
      xml_write_to_file("doit_iteration_" + os.str() + ".xml", 
                        doit_i_field);  
    }
  
  // Only the iterations given by the keyword are written to a file
  else
    {
      for (Index i = 0; i<iterations.nelem(); i++)
        {
          if (doit_iteration_counter == iterations[i])
            xml_write_to_file("doit_iteration_" + os.str() + ".xml", 
                              doit_i_field);
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
doit_scat_fieldCalc(// WS Output and Input
                    Tensor6& doit_scat_field,
                    //WS Input:
                    const Agenda& pha_mat_spt_agenda,
                    const Tensor6& doit_i_field,
                    const Tensor4& pnd_field,
                    const Tensor3& t_field,
                    const Index& atmosphere_dim,
                    const ArrayOfIndex& cloudbox_limits,
                    const Vector& scat_za_grid,
                    const Vector& scat_aa_grid,
                    const Index& doit_za_grid_size
                    )
  
{
  // ------------ Check the input -------------------------------

  // Agenda for calculation of phase matrix
  chk_not_empty( "pha_mat_spt_agenda", pha_mat_spt_agenda);

  // Number of zenith angles.
  const Index Nza = scat_za_grid.nelem();
  
  if (scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180.)
    throw runtime_error("The range of *scat_za_grid* must [0 180].");
  
  // Number of azimuth angles.
  const Index Naa = scat_aa_grid.nelem();
  
  if (scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360.)
    throw runtime_error("The range of *scat_za_grid* must [0 360]."); 

  // Get stokes dimension from *doit_scat_field*:
  const Index stokes_dim = doit_scat_field.ncols();
  assert(stokes_dim > 0 || stokes_dim < 5); 

  // Size of particle number density field can not be checked here, 
  // because the function does not use the atmospheric grids.
  // Check will be included after fixing the size of pnd field, so 
  // that it is defined only inside the cloudbox. 
  
  // Check atmospheric dimension and dimensions of 
  // radiation field (*doit_i_field*) and scattering integral field
  // (*doit_scat_field*)
  if (atmosphere_dim == 1)
    {
      assert( is_size(doit_i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 1, Nza, 1, stokes_dim));
      assert( is_size(doit_scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 1, scat_za_grid.nelem(), 1, stokes_dim));
    }
  else if (atmosphere_dim == 3)
    {
      assert ( is_size( doit_i_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                        (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                        Nza, Naa, stokes_dim));
      assert ( is_size( doit_scat_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                        (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                        Nza, Naa, stokes_dim));
    }
  else
    {
      ostringstream os;
      os << "The atmospheric dimension must be 1D or 3D \n"
         << "for scattering calculations using the DOIT\n"
         << "module, but it is not. The value of *atmosphere_dim*\n"
         << "is " << atmosphere_dim << ".";
      throw runtime_error( os.str() );
    }

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");
  
  // This function should only be used for down-looking cases where no 
  // optimized zenith angle grid is required. 
  if (doit_za_grid_size != Nza)
    throw runtime_error(
                        "The zenith angle grids for the computation of\n"
                        "the scattering integral and the RT part must \n"
                        "be equal. Check definitions in \n"
                        "*DoitAngularGridsSet*. The keyword \n"
                        "'za_grid_opt_file' should be empty. \n"
                        );

  // ------ end of checks -----------------------------------------------

  // Initialize variables *pha_mat* and *pha_mat_spt*
  Tensor4 pha_mat_local(doit_za_grid_size, scat_aa_grid.nelem(), 
                        stokes_dim, stokes_dim, 0.);
  
  Tensor5 pha_mat_spt_local(pnd_field.nbooks(), doit_za_grid_size,
                            scat_aa_grid.nelem(), stokes_dim, stokes_dim, 0.);
  
  // Equidistant step size for integration
  Vector grid_stepsize(2);
  grid_stepsize[0] = 180./(doit_za_grid_size - 1);
  grid_stepsize[1] = 360./(Naa - 1);     
  
  Tensor3 product_field(Nza, Naa, stokes_dim, 0);
 
  out2 << "  Calculate the scattered field\n";
  
  if  ( atmosphere_dim == 1 )
    {
      Index scat_aa_index_local = 0;
      
      // Get pha_mat at the grid positions
      // Since atmosphere_dim = 1, there is no loop over lat and lon grids
      for (Index p_index = 0; p_index<=cloudbox_limits[1]-cloudbox_limits[0] ;
           p_index++)
        {
          Numeric rte_temperature_local =
            t_field(p_index + cloudbox_limits[0], 0, 0);
          //There is only loop over zenith angle grid ; no azimuth angle grid.
          for (Index scat_za_index_local = 0;
               scat_za_index_local < Nza; scat_za_index_local ++)
            {
              // Dummy index
              Index index_zero = 0;
              
              // Calculate the phase matric of a single particle type
              out3 << "Calculate the phase matrix \n"; 
              pha_mat_spt_agendaExecute(pha_mat_spt_local,
                                        scat_za_index_local,
                                        index_zero,
                                        index_zero,
                                        p_index,
                                        scat_aa_index_local,
                                        rte_temperature_local,
                                        pha_mat_spt_agenda,
                                        true);
              
              // Sum over all particle types
              pha_matCalc(pha_mat_local, pha_mat_spt_local, pnd_field, 
                          atmosphere_dim, p_index, 0, 
                          0);

              out3 << "Multiplication of phase matrix with incoming" << 
                " intensities \n";
              
              product_field = 0;
              
              // za_in and aa_in are for incoming zenith and azimuth 
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
                                pha_mat_local(za_in, aa_in, i, j) * 
                                doit_i_field(p_index, 0, 0, za_in, 0, j);
                          }
                      }
                      
                    }//end aa_in loop
                }//end za_in loop
              //integration of the product of ifield_in and pha
              //  over zenith angle and azimuth angle grid. It calls
              for (Index i = 0; i < stokes_dim; i++)
                {
                  doit_scat_field( p_index, 0, 0, scat_za_index_local, 0, i)
                    = AngIntegrate_trapezoid_opti
                    (product_field(joker, joker, i),
                     scat_za_grid,
                     scat_aa_grid,
                     grid_stepsize);
                  
                }//end i loop
            }//end za_prop loop
        }//end p_index loop
    }//end atmosphere_dim = 1
  
  
  //atmosphere_dim = 3
  else if( atmosphere_dim == 3 )
    {
      /*there is a loop over pressure, latitude and longitudeindex
        when we calculate the pha_mat from pha_mat_spt and pnd_field
        using the method pha_matCalc.  */
      
      for (Index p_index = 0; p_index <=
             cloudbox_limits[1] - cloudbox_limits[0];
           p_index++)
        {
          for (Index lat_index = 0; lat_index <= 
                 cloudbox_limits[3] - cloudbox_limits[2]; lat_index++)
            {
              for (Index lon_index = 0; lon_index <= 
                     cloudbox_limits[5]-cloudbox_limits[4]; lon_index++)
                {
                  Numeric rte_temperature_local = 
                    t_field(p_index + cloudbox_limits[0],
                            lat_index + cloudbox_limits[2],
                            lon_index + cloudbox_limits[4]);
                
                  for (Index scat_aa_index_local = 1; 
                       scat_aa_index_local < Naa; 
                       scat_aa_index_local++)
                    {
                      for (Index scat_za_index_local = 0; 
                           scat_za_index_local < Nza; 
                           scat_za_index_local ++)
                        {
                          out3 << "Calculate phase matrix \n";
                          pha_mat_spt_agendaExecute(pha_mat_spt_local,
                                                    scat_za_index_local,
                                                    lat_index,
                                                    lon_index,
                                                    p_index, 
                                                    scat_aa_index_local,
                                                    rte_temperature_local,
                                                    pha_mat_spt_agenda,
                                                    true);
                          
                          pha_matCalc(pha_mat_local, pha_mat_spt_local,
                                      pnd_field, 
                                      atmosphere_dim, 
                                      p_index, 
                                      lat_index, 
                                      lon_index);
                          
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
                                            pha_mat_local
                                            (za_in, aa_in, i, j) * 
                                            doit_i_field(p_index, lat_index, 
                                                         lon_index, 
                                                         scat_za_index_local,
                                                         scat_aa_index_local,
                                                         j);
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
                              doit_scat_field( p_index,
                                               lat_index,
                                               lon_index,
                                               scat_za_index_local, 
                                               scat_aa_index_local,
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
      // aa = 0 is the same as aa = 180:
      doit_scat_field(joker, joker, joker, joker, 0, joker) =
        doit_scat_field(joker, joker, joker, joker, Naa-1, joker);
    }// end atmosphere_dim = 3
}


/* Workspace method: Doxygen documentation will be auto-generated */
void
doit_scat_fieldCalcLimb(// WS Output and Input
                        Tensor6& doit_scat_field,
                        //WS Input:
                        const Agenda& pha_mat_spt_agenda,
                        const Tensor6& doit_i_field,
                        const Tensor4& pnd_field,
                        const Tensor3& t_field,
                        const Index& atmosphere_dim,
                        const ArrayOfIndex& cloudbox_limits,
                        const Vector& scat_za_grid,
                        const Vector& scat_aa_grid,
                        const Index& doit_za_grid_size,
                        const Index& doit_za_interp
                        )
{
  // ------------ Check the input -------------------------------
   
  // Agenda for calculation of phase matrix
  chk_not_empty( "pha_mat_spt_agenda", pha_mat_spt_agenda);
   
  // Number of zenith angles.
  const Index Nza = scat_za_grid.nelem();
   
  if (scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180.)
    throw runtime_error("The range of *scat_za_grid* must [0 180].");
   
  // Number of azimuth angles.
  const Index Naa = scat_aa_grid.nelem();
   
  if (scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360.)
    throw runtime_error("The range of *scat_za_grid* must [0 360]."); 

  // Get stokes dimension from *doit_scat_field*:
  const Index stokes_dim = doit_scat_field.ncols();
  assert(stokes_dim > 0 || stokes_dim < 5); 

  // Size of particle number density field can not be checked here, 
  // because the function does not use the atmospheric grids.
  // Check will be included after fixing the size of pnd field, so 
  // that it is defined only inside the cloudbox. 
  
  // Check atmospheric dimension and dimensions of 
  // radiation field (*doit_i_field*) and scattering integral field
  // (*doit_scat_field*)
  if (atmosphere_dim == 1)
    {
      assert( is_size(doit_i_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 1, Nza, 1, stokes_dim));
      assert( is_size(doit_scat_field, 
                      (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                      1, 1, scat_za_grid.nelem(), 1, stokes_dim));
    }
  else if (atmosphere_dim == 3)
    {
      assert ( is_size( doit_i_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                        (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                        Nza, Naa, stokes_dim));
      assert ( is_size( doit_scat_field, 
                        (cloudbox_limits[1] - cloudbox_limits[0]) +1,
                        (cloudbox_limits[3] - cloudbox_limits[2]) +1, 
                        (cloudbox_limits[5] - cloudbox_limits[4]) +1,
                        Nza, Naa, stokes_dim));
    }
  else
    {
      ostringstream os;
      os << "The atmospheric dimension must be 1D or 3D \n"
         << "for scattering calculations using the DOIT\n"
         << "module, but it is not. The value of *atmosphere_dim*\n"
         << "is " << atmosphere_dim << ".";
      throw runtime_error( os.str() );
    }
  
  if( !(doit_za_interp == 0  ||  doit_za_interp == 1 ) )
    throw runtime_error( "Interpolation method is not defined. Use \n"
                         "*doit_za_interpSet*.\n");

  if ( cloudbox_limits.nelem()!= 2*atmosphere_dim)
    throw runtime_error(
                        "*cloudbox_limits* is a vector which contains the"
                        "upper and lower limit of the cloud for all "
                        "atmospheric dimensions. So its dimension must"
                        "be 2 x *atmosphere_dim*");
  
  if (doit_za_grid_size < 16)
    throw runtime_error(
                        "*doit_za_grid_size* must be greater than 15 for"
                        "accurate results");
  else if (doit_za_grid_size > 100)
    out1 << "Warning: doit_za_grid_size is very large which means that the \n"
         << "calculation will be very slow. The recommended value is 19.\n";
  
  // ------ end of checks -----------------------------------------------
  
  // Initialize variables *pha_mat* and *pha_mat_spt*
  Tensor4 pha_mat_local(doit_za_grid_size, scat_aa_grid.nelem(), 
                        stokes_dim, stokes_dim, 0.);

  Tensor5 pha_mat_spt_local(pnd_field.nbooks(), doit_za_grid_size,
                            scat_aa_grid.nelem(), stokes_dim, stokes_dim, 0.);


  // Create the grids for the calculation of the scattering integral.
  Vector za_grid;
  nlinspace(za_grid, 0, 180, doit_za_grid_size);
 
  // Two interpolations are required. First we have to interpolate the 
  // intensity field on the equidistant grid: 
  ArrayOfGridPos gp_za_i(doit_za_grid_size);
  gridpos(gp_za_i, scat_za_grid, za_grid);
  
  Matrix itw_za_i(doit_za_grid_size, 2);
  interpweights(itw_za_i, gp_za_i);

  // Intensity field interpolated on equidistant grid.
  Matrix doit_i_field_int(doit_za_grid_size, stokes_dim, 0);

  // Second, we have to interpolate the scattering integral on the RT
  // zenith angle grid.
  ArrayOfGridPos gp_za(Nza);
  gridpos(gp_za, za_grid, scat_za_grid);

  Matrix itw_za(Nza, 2);
  interpweights(itw_za, gp_za);
  
  // Original scattered field, on equidistant zenith angle grid.
  Matrix doit_scat_field_org(doit_za_grid_size, stokes_dim, 0);
  
  //  Grid stepsize of zenith and azimuth angle grid, these are needed for the 
  // integration function. 
  Vector grid_stepsize(2);
  grid_stepsize[0] = 180./(doit_za_grid_size - 1);
  grid_stepsize[1] = 360./(Naa - 1);
    
  Tensor3 product_field(doit_za_grid_size, Naa, stokes_dim, 0);

  if  ( atmosphere_dim == 1 )
    {
      Index scat_aa_index_local = 0;
      
      // Get pha_mat at the grid positions
      // Since atmosphere_dim = 1, there is no loop over lat and lon grids
      for (Index p_index = 0;
           p_index <= cloudbox_limits[1]-cloudbox_limits[0];
           p_index++)
        {
          Numeric rte_temperature_local = 
            t_field(p_index + cloudbox_limits[0], 0, 0);
          // Interpolate intensity field:
          for (Index i = 0; i < stokes_dim; i++)
            {
              if (doit_za_interp == 0)
                {
                  interp(doit_i_field_int(joker, i), itw_za_i, 
                         doit_i_field(p_index, 0, 0, joker, 0, i), gp_za_i);
                } 
              else if (doit_za_interp == 1)
                {
                  // Polynomial
                  for(Index za = 0; za < za_grid.nelem(); za++)
                    {
                      doit_i_field_int(za, i) = 
                        interp_poly(scat_za_grid, 
                                     doit_i_field(p_index, 0, 0, joker, 0, i),
                                     za_grid[za],
                                     gp_za_i[za]);
                    }
                }
              // doit_za_interp must be 0 or 1 (linear or polynomial)!!!
              else assert(false);
            }       
          
          //There is only loop over zenith angle grid; no azimuth angle grid.
          for( Index scat_za_index_local = 0;
               scat_za_index_local < doit_za_grid_size;
               scat_za_index_local++)
            {
              // Dummy index
              Index index_zero = 0;
              
              // Calculate the phase matrix of a single particle type
              out3 << "Calculate the phase matrix \n"; 
              pha_mat_spt_agendaExecute(pha_mat_spt_local,
                                        scat_za_index_local,
                                        index_zero,
                                        index_zero,
                                        p_index,
                                        scat_aa_index_local,
                                        rte_temperature_local,
                                        pha_mat_spt_agenda,
                                        true);
              
              // Sum over all particle types
              pha_matCalc(pha_mat_local, pha_mat_spt_local, pnd_field, 
                          atmosphere_dim, p_index, 0, 
                          0);

              out3 << "Multiplication of phase matrix with incoming" << 
                " intensities \n";
            
              product_field = 0;

              // za_in and aa_in are for incoming zenith and azimuth 
              // angle direction for which pha_mat is calculated
              for( Index za_in = 0; za_in < doit_za_grid_size; za_in ++)
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
                                pha_mat_local(za_in, aa_in, i, j) * 
                                doit_i_field_int(za_in, j);
                            }
                        }
                      
                    }//end aa_in loop
                }//end za_in loop
            
              out3 << "Compute integral. \n"; 
              for (Index i = 0; i < stokes_dim; i++)
                {
                  doit_scat_field_org(scat_za_index_local, i)=
                    AngIntegrate_trapezoid_opti(product_field(joker, joker, i),
                                                za_grid,
                                                scat_aa_grid,
                                                grid_stepsize);
                }//end i loop
            }//end za_prop loop
          
          // Interpolation on scat_za_grid, which is used in 
          //radiative transferpart.
          for (Index i = 0; i < stokes_dim; i++)
            {
            if(doit_za_interp == 0) // linear interpolation
              {
                interp(doit_scat_field(p_index,
                                  0,
                                  0,
                                  joker,
                                  0,
                                  i),
                       itw_za,
                       doit_scat_field_org(joker, i),
                       gp_za);
              }
            else // polynomial interpolation
              {
                for(Index za = 0; za < scat_za_grid.nelem(); za++)
                  {
                    doit_scat_field(p_index, 0, 0, za, 0, i) = 
                      interp_poly(za_grid, 
                                   doit_scat_field_org(joker, i),
                                   scat_za_grid[za],
                                   gp_za[za]);
                  }
              }
          }
        }//end p_index loop
      
    }//end atmosphere_dim = 1
  
  
  else if( atmosphere_dim == 3 ){
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
                
                Numeric rte_temperature_local =
                  t_field(p_index + cloudbox_limits[0],
                          lat_index + cloudbox_limits[2],
                          lon_index + cloudbox_limits[4]);
                
                // Loop over scattered directions
                for (Index scat_aa_index_local = 1;
                     scat_aa_index_local < Naa; 
                     scat_aa_index_local++)
                  {
                   // Interpolate intensity field:
                    for (Index i = 0; i < stokes_dim; i++)
                      {
                        interp(doit_i_field_int(joker, i), itw_za_i, 
                               doit_i_field(p_index, lat_index, lon_index,
                                       joker, scat_aa_index_local, i), gp_za_i);
                      }       
                    
                    for (Index scat_za_index_local = 0;
                         scat_za_index_local < doit_za_grid_size;
                         scat_za_index_local++)
                      {
                        
                        out3 << "Calculate phase matrix \n";
                        pha_mat_spt_agendaExecute(pha_mat_spt_local,
                                                  scat_za_index_local,
                                                  lat_index,
                                                  lon_index,
                                                  p_index, 
                                                  scat_aa_index_local,
                                                  rte_temperature_local,
                                                  pha_mat_spt_agenda,
                                                  true);
  
                        pha_matCalc(pha_mat_local, pha_mat_spt_local,
                                    pnd_field, 
                                    atmosphere_dim, 
                                    p_index, 
                                    lat_index, 
                                    lon_index);
                        
                        product_field = 0;
                        
                        
                        //za_in and aa_in are the incoming directions
                        //for which pha_mat_spt is calculated
                        out3 << "Multiplication of phase matrix with" << 
                          "incoming intensity \n";
                        
                        for( Index za_in = 0; za_in < doit_za_grid_size; za_in ++)
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
                                          pha_mat_local(za_in, aa_in, i, j) * 
                                          doit_i_field_int(za_in, j);
                                      }
                                  }
                              }//end aa_in loop
                          }//end za_in loop
                        
                        out3 << "Compute the integral \n";

                        for (Index i = 0; i < stokes_dim; i++)
                          {
                            doit_scat_field_org(scat_za_index_local, i)  =  
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
                        interp(doit_scat_field(p_index,
                                               lat_index,
                                               lon_index,
                                               joker,
                                               scat_aa_index_local,
                                               i),
                               itw_za,
                               doit_scat_field_org(joker, i),
                               gp_za);
                      }
                  } // end aa_prop loop
              }//end lon loop
          }//end lat loop
      }// end p loop
    doit_scat_field(joker, joker, joker, joker, 0, joker) =
      doit_scat_field(joker, joker, joker, joker, Naa-1, joker);
  }// end atm_dim=3
  out2 << "  Finished scattered field.\n"; 
}


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_za_grid_optCalc(//WS Output
                          Vector& doit_za_grid_opt,
                          // WS Input:
                          const Tensor6& doit_i_field,
                          const Vector& scat_za_grid,
                          const Index& doit_za_interp,
                          //Keywords:
                          const Numeric& acc
                          )
{
  //-------- Check the input ---------------------------------
  
  // Here it is checked whether doit_i_field is 1D and whether it is 
  // consistent with scat_za_grid. The number of pressure levels and the 
  // number of stokes components does not matter. 
  chk_size("doit_i_field", doit_i_field, 
           doit_i_field.nvitrines() , 1, 1, scat_za_grid.nelem(), 1,
           doit_i_field.ncols());
  
  if(doit_i_field.ncols()<1 || doit_i_field.ncols()>4)
    throw runtime_error("The last dimension of *doit_i_field* corresponds\n"
                        "to the Stokes dimension, therefore the number of\n"
                        "columns in *doit_i_field* must be a number between\n"
                        "1 and 4, but it is not!");
  
  if(scat_za_grid.nelem() < 500)
    throw runtime_error("The fine grid (*scat_za_grid*) has less than \n"
                        "500 grid points which is not sufficient for \n"
                        "grid_optimization");

  if( !(doit_za_interp == 0  ||  doit_za_interp == 1 ) )
    throw runtime_error( "Interpolation method is not defined. Use \n"
                         "*doit_za_interpSet*.\n");

  // ------------- end of checks ------------------------------------- 
  
  // Here only used as dummy variable. 
  Matrix doit_i_field_opt_mat;
  doit_i_field_opt_mat = 0.;
  
  // Optimize zenith angle grid. 
  za_gridOpt(doit_za_grid_opt, doit_i_field_opt_mat,
             scat_za_grid, doit_i_field, acc,
             doit_za_interp);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void doit_za_interpSet(
                       Index& doit_za_interp,
                       const Index& atmosphere_dim,
                       //Keyword
                       const String& method
                       
                       )
{
  chk_if_in_range( "atmosphere_dim", atmosphere_dim, 1, 3 );

  if (atmosphere_dim != 1 && method == "polynomial")
    throw runtime_error(
                        "Polynomial interpolation is only implemented for\n"
                        "1D DOIT calculations as \n"   
                        "in 3D there can be numerical problems.\n"
                        "Please use 'linear' interpolation method."
                        );

  if(method == "linear")
    doit_za_interp = 0;
  else if (method == "polynomial")
    doit_za_interp = 1;
  else
    throw runtime_error("Possible interpolation methods are 'linear' "
                        "and 'polynomial'.\n");
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ScatteringDoit(
                    Tensor6& doit_i_field,
                    Tensor7& scat_i_p, 
                    Tensor7& scat_i_lat, 
                    Tensor7& scat_i_lon,
                    Tensor4& doit_i_field1D_spectrum,
                    const Vector& f_grid,
                    const Agenda& doit_mono_agenda,
                    const Index& doit_is_initialized
                    )
                  
{
  //-------- Check input -------------------------------------------
  
  chk_not_empty( "doit_mono_agenda", doit_mono_agenda );

  // Frequency grid
  //
  if( f_grid.nelem() == 0 )
    throw runtime_error( "The frequency grid is empty." );
  chk_if_increasing( "f_grid", f_grid );

  // Check whether DoitInit was executed
  if (!doit_is_initialized)
    throw runtime_error(
                        "Initialization method *DoitInit* has to be "
                        "put before\n"
                        "start of *ScatteringDoit*");

  //-------- end of checks ----------------------------------------

  for (Index f_index = 0; f_index < f_grid.nelem(); f_index ++)
    {
      out1 << "Frequency: " << f_grid[f_index]/1e9 <<" GHz \n" ;
      doit_mono_agendaExecute(doit_i_field, scat_i_p, scat_i_lat,
                              scat_i_lon, doit_i_field1D_spectrum,
                              f_index, doit_mono_agenda,
                              false); 
    }
}


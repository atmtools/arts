/* Copyright (C) 2002-2012 Claudia Emde <claudia.emde@dlr.de>
                      
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
  \file   doit.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Wed Jun 04 11:03:57 2003
  
  \brief  This file contains functions to calculate the radiative transfer
  inside the cloudbox using the DOIT method.
  
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
#include "doit.h"
#include "logic.h"
#include "check_input.h"
#include "sorting.h"
#include "cloudbox.h"
#include "geodetic.h"

extern const Numeric PI;
extern const Numeric RAD2DEG;



//! rte_step_doit
/*!
    Solves monochromatic VRTE for an atmospheric slab with constant 
    conditions.

    The function can be used for cloudbox calculations.

    The function is best explained by considering a homogenous layer. That is,
    the physical conditions inside the layer are constant. In reality they
    are not constant, so in practical all coefficients have to be averaged 
    before calling this function. 
    Total extinction and absorption inside the layer are described by
    *ext_mat_av* and *abs_vec_av* respectively,
    the blackbdody radiation of the layer is given by *rte_planck_value*
    and the propagation path length through the layer is *lstep*.

    There is an additional scattering source term in the 
    VRTE, the scattering integral term. For this function a constant
    scattering term is assumed. The radiative transfer step is only a part 
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
    \param   trans_mat          Output: Transmission matrix of slab.
    \param   ext_mat_av         Input: Averaged extinction matrix.
    \param   abs_vec_av         Input: Averaged absorption vector.
    \param   sca_vec_av         Input: averaged scattering vector.
    \param   lstep              Input: The length of the RTE step.
    \param   rtp_planck_value   Input: Blackbody radiation.

    \author Claudia Emde and Patrick Eriksson, 
    \date   2002-11-22
*/
void rte_step_doit(//Output and Input:
              VectorView stokes_vec,
              MatrixView trans_mat,
              //Input
              ConstMatrixView ext_mat_av,
              ConstVectorView abs_vec_av,
              ConstVectorView sca_vec_av,
              const Numeric& lstep,
              const Numeric& rtp_planck_value,
              const bool& trans_is_precalc )
{
  //Stokes dimension:
  Index stokes_dim = stokes_vec.nelem();

  //Check inputs:
  assert(is_size(trans_mat, stokes_dim, stokes_dim));
  assert(is_size(ext_mat_av, stokes_dim, stokes_dim));
  assert(is_size(abs_vec_av, stokes_dim));
  assert(is_size(sca_vec_av, stokes_dim));
  assert( rtp_planck_value >= 0 );
  assert( lstep >= 0 );
  assert (!is_singular( ext_mat_av ));

  // Check, if only the first component of abs_vec is non-zero:
  bool unpol_abs_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (abs_vec_av[i] != 0)
      unpol_abs_vec = false;

  bool unpol_sca_vec = true;

  for (Index i = 1; i < stokes_dim; i++)
    if (sca_vec_av[i] != 0)
      unpol_sca_vec = false;

  // Calculate transmission by general function, if not precalculated
  Index extmat_case = 0;
  if( !trans_is_precalc )
    { ext2trans( trans_mat, extmat_case, ext_mat_av, lstep ); }

  //--- Scalar case: ---------------------------------------------------------
  if( stokes_dim == 1 )
    {
      stokes_vec[0]  = stokes_vec[0] * trans_mat(0,0) +
        ( abs_vec_av[0] * rtp_planck_value + sca_vec_av[0] ) / 
        ext_mat_av(0,0) * (1 - trans_mat(0,0) );
    }


  //--- Vector case: ---------------------------------------------------------

  // We have here two cases, diagonal or non-diagonal ext_mat_gas
  // For diagonal ext_mat_gas, we expect abs_vec_gas to only have a
  // non-zero value in position 1.

  //- Unpolarised
  else if( extmat_case==1 && unpol_abs_vec && unpol_sca_vec )
    {
      // Stokes dim 1
      stokes_vec[0] = stokes_vec[0] * trans_mat(0,0) +
                      ( abs_vec_av[0] * rtp_planck_value + sca_vec_av[0] ) /
                      ext_mat_av(0,0) * (1 - trans_mat(0,0) );

      // Stokes dims > 1
      for( Index i=1; i<stokes_dim; i++ )
        {
          stokes_vec[i]  = stokes_vec[i] * trans_mat(i,i) + sca_vec_av[i] / 
                           ext_mat_av(i,i)  * (1 - trans_mat(i,i));
        }
    }


  //- General case
  else
    {
      //Initialize internal variables:

      // Matrix LU used for LU decompostion and as dummy variable:
      Matrix LU(stokes_dim, stokes_dim);
      ArrayOfIndex indx(stokes_dim); // index for pivoting information
      Vector b(stokes_dim); // dummy variable
      Vector x(stokes_dim); // solution vector for K^(-1)*b
      Matrix I(stokes_dim, stokes_dim);

      Vector B_abs_vec(stokes_dim);
      B_abs_vec = abs_vec_av;
      B_abs_vec *= rtp_planck_value;
      
      for (Index i=0; i<stokes_dim; i++)
        b[i] = B_abs_vec[i] + sca_vec_av[i];  // b = abs_vec * B + sca_vec

      // solve K^(-1)*b = x
      ludcmp(LU, indx, ext_mat_av);
      lubacksub(x, LU, b, indx);

      Matrix ext_mat_ds(stokes_dim, stokes_dim);
      ext_mat_ds = ext_mat_av;
      ext_mat_ds *= -lstep; // ext_mat_ds = -ext_mat*ds

      Vector term1(stokes_dim);
      Vector term2(stokes_dim);

      id_mat(I);
      for(Index i=0; i<stokes_dim; i++)
        {
          for(Index j=0; j<stokes_dim; j++)
            LU(i,j) = I(i,j) - trans_mat(i,j); // take LU as dummy variable
        }
      mult(term2, LU, x); // term2: second term of the solution of the RTE with
                          //fixed scattered field

      // term1: first term of solution of the RTE with fixed scattered field
      mult( term1, trans_mat, stokes_vec );

      for (Index i=0; i<stokes_dim; i++)
        stokes_vec[i] = term1[i] + term2[i];  // Compute the new Stokes Vector
    }
}




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
  \param ws Current Workspace
  \param ext_mat_field extinction matrix field
  \param abs_vec_field absorption vector field
  // Input
  \param spt_calc_agenda Agenda for calculation of single scattering properties
  \param opt_prop_part_agenda Agenda for summing over all scattering elements
  \param scat_za_index Indices for
  \param scat_aa_index    propagation direction
  \param cloudbox_limits Cloudbox limits.
  \param t_field Temperature field
  \param pnd_field Particle number density field. 

  \author Claudia Emde
  \date 2002-06-03
*/
void cloud_fieldsCalc(Workspace& ws,
                      // Output and Input:
                      Tensor5View ext_mat_field,
                      Tensor4View abs_vec_field,
                      // Input:
                      const Agenda& spt_calc_agenda,
                      const Agenda& opt_prop_part_agenda,
                      const Index& scat_za_index, 
                      const Index& scat_aa_index,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstTensor3View t_field, 
                      ConstTensor4View pnd_field,
                      const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Input variables are checked in the WSMs i_fieldUpdateSeqXXX, from 
  // where this function is called.
  
  out3 << "Calculate scattering properties in cloudbox \n";
  
  const Index atmosphere_dim = cloudbox_limits.nelem()/2;
  const Index N_se = pnd_field.nbooks();
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
  // of all scattering elements.
  Matrix abs_vec_spt_local(N_se, stokes_dim, 0.);
  Tensor3 ext_mat_spt_local(N_se, stokes_dim, stokes_dim, 0.);
  Matrix abs_vec_local;
  Tensor3 ext_mat_local;
  Numeric rtp_temperature_local;
  
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
                rtp_temperature_local = 
                  t_field(scat_p_index_local + cloudbox_limits[0], 0, 0);
              else
                rtp_temperature_local = 
                  t_field(scat_p_index_local + cloudbox_limits[0],
                          scat_lat_index_local + cloudbox_limits[2],
                          scat_lon_index_local + cloudbox_limits[4]);
              
              //Calculate optical properties for individual scattering elements:
              //( Execute agendas silently. )
              spt_calc_agendaExecute(ws, ext_mat_spt_local, 
                                     abs_vec_spt_local,
                                     scat_p_index_local,
                                     scat_lat_index_local,
                                     scat_lon_index_local, 
                                     rtp_temperature_local,
                                     scat_za_index,
                                     scat_aa_index,
                                     spt_calc_agenda);

              opt_prop_part_agendaExecute(ws, ext_mat_local, abs_vec_local, 
                                          ext_mat_spt_local, 
                                          abs_vec_spt_local,
                                          scat_p_index_local,
                                          scat_lat_index_local,
                                          scat_lon_index_local,
                                          opt_prop_part_agenda);
           
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

  \param[in,out] ws Current Workspace
  WS Output:
  \param[out] doit_i_field_mono Updated radiation field inside the cloudbox.
  WS Input:
  \param p_index // Pressure index
  \param scat_za_index // Index for proagation direction
  \param scat_za_grid
  \param cloudbox_limits 
  \param doit_scat_field Scattered field.
  Calculate scalar gas absorption:
  \param propmat_clearsky_agenda
  \param vmr_field
  Propagation path calculation:
  \param ppath_step_agenda
  \param p_grid
  \param z_field
  \param refellipsoid
  \param z_surface
  Calculate thermal emission:
  \param t_field
  \param f_grid
  \param f_index
  Optical properties of particles
  \param ext_mat_field
  \param abs_vec_field
  \param surface_rtprop_agenda
  \param scat_za_interp

  \author Claudia Emde
  \date 2003-06-04
*/
void cloud_ppath_update1D(Workspace& ws,
                          // Input and output
                          Tensor6View doit_i_field_mono,
                          // ppath_step_agenda:
                          const Index& p_index,
                          const Index& scat_za_index,
                          ConstVectorView scat_za_grid,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstTensor6View doit_scat_field,
                          // Calculate scalar gas absorption:
                          const Agenda& propmat_clearsky_agenda,
                          ConstTensor4View vmr_field,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          const Numeric&   ppath_lraytrace,
                          ConstVectorView  p_grid,
                          ConstTensor3View z_field,
                          ConstVectorView refellipsoid,
                          // Calculate thermal emission:
                          ConstTensor3View t_field,
                          ConstVectorView f_grid,
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          const Agenda& surface_rtprop_agenda,
                          //const Agenda& surface_rtprop_agenda, 
                          const Index& scat_za_interp,
                          const Verbosity& verbosity)
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
  ppath_step.pos(0,0) = z_field(p_index,0,0);
  ppath_step.r[0]     = refellipsoid[0] + z_field(p_index,0,0);
  
  // Define the direction:
  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
  
  
  // Define the grid positions:
  ppath_step.gp_p[0].idx   = p_index;
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;
  
  // Call ppath_step_agenda: 
  ppath_step_agendaExecute( ws, ppath_step, ppath_lraytrace, t_field, z_field, 
                            vmr_field, 
                            Vector(1,f_grid[f_index]), ppath_step_agenda );
  
  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  
  if((cloudbox_limits[0] <= ppath_step.gp_p[1].idx &&
     cloudbox_limits[1] > ppath_step.gp_p[1].idx) ||
     (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
      abs(ppath_step.gp_p[1].fd[0]) < 1e-6))
    {
      // Stokes dimension
      const Index stokes_dim = doit_i_field_mono.ncols();
      // Number of species
      const Index N_species = vmr_field.nbooks();
      
      // Ppath_step normally has 2 points, the starting
      // point and the intersection point.
      // But there can be points in between, when a maximum 
      // lstep is given. We have to interpolate on all the 
      // points in the ppath_step.

      // Initialize variables for interpolated values
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np, 0.);
      Matrix abs_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix doit_i_field_mono_int(stokes_dim, ppath_step.np, 0.);
      Vector t_int(ppath_step.np, 0.);
      Matrix vmr_list_int(N_species, ppath_step.np, 0.);
      Vector p_int(ppath_step.np, 0.);
      
      interp_cloud_coeff1D(ext_mat_int, abs_vec_int, sca_vec_int,
                           doit_i_field_mono_int, t_int, vmr_list_int, p_int, 
                           ext_mat_field, abs_vec_field, doit_scat_field, 
                           doit_i_field_mono, t_field, vmr_field, p_grid,
                           ppath_step, cloudbox_limits, scat_za_grid, 
                           scat_za_interp, verbosity);
     
          
      // ppath_what_background(ppath_step) tells the 
      // radiative background.  More information in the 
      // function get_iy_of_background.
      // if there is no background we proceed the RT
      Index bkgr = ppath_what_background(ppath_step);

      // Radiative transfer from one layer to the next, starting
      // at the intersection with the next layer and propagating
      // to the considered point.
      cloud_RT_no_background(ws, doit_i_field_mono, 
                             propmat_clearsky_agenda,
                                 ppath_step,
                             t_int, vmr_list_int,
                             ext_mat_int, abs_vec_int, sca_vec_int,
                                 doit_i_field_mono_int,
                             p_int, cloudbox_limits, 
                                 f_grid, f_index, p_index, 0, 0,
                             scat_za_index, 0, verbosity);
      
      // bkgr=2 indicates that the background is the surface
      if (bkgr == 2)
        {
          // cout << "hit surface "<< ppath_step.gp_p << endl;
          cloud_RT_surface(ws,
                           doit_i_field_mono, surface_rtprop_agenda, f_grid,
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
void cloud_ppath_update1D_noseq(Workspace& ws,
                                // Output
                                Tensor6View doit_i_field_mono,
                                // ppath_step_agenda:
                                const Index& p_index,
                                const Index& scat_za_index,
                                ConstVectorView scat_za_grid,
                                const ArrayOfIndex& cloudbox_limits,
                                ConstTensor6View doit_i_field_mono_old,
                                ConstTensor6View doit_scat_field,
                                // Calculate scalar gas absorption:
                                const Agenda& propmat_clearsky_agenda,
                                ConstTensor4View vmr_field,
                                // Propagation path calculation:
                                const Agenda& ppath_step_agenda,
                                const Numeric&   ppath_lraytrace,
                                ConstVectorView  p_grid,
                                ConstTensor3View z_field,
                                ConstVectorView refellipsoid,
                                // Calculate thermal emission:
                                ConstTensor3View t_field,
                                ConstVectorView f_grid,
                                // used for surface ?
                                const Index& f_index,
                                //particle optical properties
                                ConstTensor5View ext_mat_field,
                                ConstTensor4View abs_vec_field,
                                const Agenda& surface_rtprop_agenda,
                                const Index& scat_za_interp,
                                const Verbosity& verbosity)
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
  ppath_step.pos(0,0) = z_field(p_index,0,0);
  ppath_step.r[0]     = refellipsoid[0] + z_field(p_index,0,0);
  
  // Define the direction:
  ppath_step.los(0,0) = scat_za_grid[scat_za_index];
  
  
  // Define the grid positions:
  ppath_step.gp_p[0].idx   = p_index;
  ppath_step.gp_p[0].fd[0] = 0;
  ppath_step.gp_p[0].fd[1] = 1;
  
  // Call ppath_step_agenda: 
  ppath_step_agendaExecute( ws, ppath_step, ppath_lraytrace, t_field, z_field, 
                            vmr_field, 
                            Vector(1,f_grid[f_index]), ppath_step_agenda );
  
  // Check whether the next point is inside or outside the
  // cloudbox. Only if the next point lies inside the
  // cloudbox a radiative transfer step caclulation has to
  // be performed.
  
  if((cloudbox_limits[0] <= ppath_step.gp_p[1].idx &&
     cloudbox_limits[1] > ppath_step.gp_p[1].idx) ||
     (cloudbox_limits[1] == ppath_step.gp_p[1].idx &&
      abs(ppath_step.gp_p[1].fd[0]) < 1e-6))
    {
      // Stokes dimension
      const Index stokes_dim = doit_i_field_mono.ncols();
      // Number of species
      const Index N_species = vmr_field.nbooks();
      
      // Ppath_step normally has 2 points, the starting
      // point and the intersection point.
      // But there can be points in between, when a maximum 
      // lstep is given. We have to interpolate on all the 
      // points in the ppath_step.

      // Initialize variables for interpolated values
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np, 0.);
      Matrix abs_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix doit_i_field_mono_int(stokes_dim, ppath_step.np, 0.);
      Vector t_int(ppath_step.np, 0.);
      Matrix vmr_list_int(N_species, ppath_step.np, 0.);
      Vector p_int(ppath_step.np, 0.);
      
      interp_cloud_coeff1D(ext_mat_int, abs_vec_int, sca_vec_int,
                           doit_i_field_mono_int, t_int, vmr_list_int, p_int, 
                           ext_mat_field, abs_vec_field, doit_scat_field, 
                           doit_i_field_mono_old, t_field, vmr_field, p_grid, 
                           ppath_step, cloudbox_limits, scat_za_grid, 
                           scat_za_interp, verbosity);
      
      // ppath_what_background(ppath_step) tells the 
      // radiative background.  More information in the 
      // function get_iy_of_background.
      // if there is no background we proceed the RT
      Index bkgr = ppath_what_background(ppath_step);
      
      // if 0, there is no background
      // do this in any case. cause we need downwelling doit_i_field_mono
      // at the surface for calculating surface scattering

      // Radiative transfer from one layer to the next, starting
      // at the intersection with the next layer and propagating
      // to the considered point.
      cloud_RT_no_background(ws, doit_i_field_mono,
                             propmat_clearsky_agenda,
                             ppath_step,
                             t_int, vmr_list_int,
                             ext_mat_int, abs_vec_int, sca_vec_int,
                             doit_i_field_mono_int,
                             p_int, cloudbox_limits, 
                             f_grid, f_index, p_index, 0, 0, 
                             scat_za_index, 0, verbosity);

      if (bkgr == 2)
        {
          cloud_RT_surface( ws,
                            doit_i_field_mono, surface_rtprop_agenda, f_grid,
                            f_index, stokes_dim, ppath_step, cloudbox_limits, 
                            scat_za_grid, scat_za_index); 
        
        }
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
  Stokes vectors are stored in the WSV doit_i_field_mono.

  \param[in,out] ws Current workspace
 WS Output:
  \param doit_i_field_mono Updated radiation field inside the cloudbox.
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
  \param propmat_clearsky_agenda
  \param vmr_field
  Propagation path calculation:
  \param ppath_step_agenda
  \param p_grid
  \param lat_grid
  \param lon_grid
  \param z_field
  \param refellipsoid
  \param z_surface
  Calculate thermal emission:
  \param t_field
  \param f_grid
  \param f_index
  \param ext_mat_field
  \param abs_vec_field

  \author Claudia Emde
  \date 2003-06-04
*/
void cloud_ppath_update3D(Workspace& ws,
                          Tensor6View doit_i_field_mono,
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
                          const Agenda& propmat_clearsky_agenda,
                          ConstTensor4View vmr_field,
                          // Propagation path calculation:
                          const Agenda& ppath_step_agenda,
                          const Numeric&  ppath_lraytrace,
                          ConstVectorView p_grid,
                          ConstVectorView lat_grid,
                          ConstVectorView lon_grid,
                          ConstTensor3View z_field,
                          ConstVectorView refellipsoid,
                          // Calculate thermal emission:
                          ConstTensor3View t_field,
                          ConstVectorView f_grid,
                          const Index& f_index,
                          //particle optical properties
                          ConstTensor5View ext_mat_field,
                          ConstTensor4View abs_vec_field,
                          const Index&, //scat_za_interp
                          const Verbosity& verbosity
                          )
{
  CREATE_OUT3;
  
  Ppath ppath_step;
  const Index stokes_dim = doit_i_field_mono.ncols();
  
  Vector sca_vec_av(stokes_dim,0);
  Vector aa_grid(scat_aa_grid.nelem());

  for(Index i = 0; i<scat_aa_grid.nelem(); i++)
    aa_grid[i] = scat_aa_grid[i]-180.;

   //Initialize ppath for 3D.
  ppath_init_structure(ppath_step, 3, 1);
  // See documentation of ppath_init_structure for
  // understanding the parameters.
              
  // The first dimension of pos are the points in 
  // the propagation path. 
  // Here we initialize the first point.
  // The second is: radius, latitude, longitude

  
  ppath_step.pos(0,2) = lon_grid[lon_index];
  ppath_step.pos(0,1) = lat_grid[lat_index];
  ppath_step.pos(0,0) = z_field( p_index, lat_index, lon_index );
  // As always on top of the lat. grid positions, OK to call refell2r:
  ppath_step.r[0] = refell2r( refellipsoid, ppath_step.pos(0,1) ) + 
                    ppath_step.pos(0,0);

              
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
  ppath_step_agendaExecute( ws, ppath_step, ppath_lraytrace, t_field, z_field, 
                            vmr_field, 
                            Vector(1,f_grid[f_index]), ppath_step_agenda);

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
        }
      const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
      const Index n2 = cloudbox_limits[3] - cloudbox_limits[2];
      const Index n3 = cloudbox_limits[5] - cloudbox_limits[4];
      gridpos_upperend_check( cloud_gp_p[0], n1 );
      gridpos_upperend_check( cloud_gp_p[ppath_step.np-1],   n1 );
      gridpos_upperend_check( cloud_gp_lat[0], n2 );
      gridpos_upperend_check( cloud_gp_lat[ppath_step.np-1], n2);
      gridpos_upperend_check( cloud_gp_lon[0], n3 );
      gridpos_upperend_check( cloud_gp_lon[ppath_step.np-1], n3 );
      
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
      // lstep is given. We have to interpolate on all the 
      // points in the ppath_step.
      
      Tensor3 ext_mat_int(stokes_dim, stokes_dim, ppath_step.np);
      Matrix abs_vec_int(stokes_dim, ppath_step.np);
      Matrix sca_vec_int(stokes_dim, ppath_step.np, 0.);
      Matrix doit_i_field_mono_int(stokes_dim, ppath_step.np, 0.);
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
          out3 << "Interpolate doit_i_field_mono:\n";
          interp( doit_i_field_mono_int(i, joker), itw_p_za, 
                  doit_i_field_mono(joker, joker, joker, joker, joker, i), 
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
      cloud_RT_no_background(ws, doit_i_field_mono, 
                             propmat_clearsky_agenda,
                             ppath_step,
                             t_int, vmr_list_int,
                             ext_mat_int, abs_vec_int, sca_vec_int,
                             doit_i_field_mono_int,
                             p_int, cloudbox_limits, 
                             f_grid, f_index, p_index, lat_index, lon_index, 
                             scat_za_index, scat_aa_index, verbosity);
    }//end if inside cloudbox
}

//! cloud_RT_no_background
/*
  This function calculates RT in the cloudbox (1D) if the next intersected 
  level is an atmospheric level (in contrast to the surface). 
  It is used inside the functions cloud_ppath_update1DXXX.
  
  Output: 
  \param doit_i_field_mono Radiation field in cloudbox. This variable is filled
  with the updated values for a given zenith angle (scat_za_index) and 
  pressure (p_index).
  Input:
  \param propmat_clearsky_agenda Calculate gas absorption.
  \param ppath_step Propagation path step from one pressure level to the next 
  (this can include several points)
  \param t_int Temperature values interpolated on propagation path points.
  \param vmr_list_int Interpolated volume mixing ratios. 
  \param ext_mat_int Interpolated total particle extinction matrix.
  \param abs_vec_int Interpolated total particle absorption vector. 
  \param sca_vec_int Interpolated total particle scattering vector. 
  \param doit_i_field_int Interpolated radiances. 
  \param p_int Interpolated pressure values. 
  \param cloudbox_limits Cloudbox_limits. 
  \param f_grid Frequency grid.
  \param f_index Frequency index of (monochromatic) scattering calculation. 
  \param p_index Pressure index in *doit_i_field_mono*.
  \param scat_za_index Zenith angle index in *doit_i_field_mono*.
  
  \author Claudia Emde
  \date 2005-05-13
*/
void cloud_RT_no_background(Workspace& ws,
                            //Output
                            Tensor6View doit_i_field_mono,
                            // Input
                            const Agenda& propmat_clearsky_agenda,
                            const Ppath& ppath_step, 
                            ConstVectorView t_int,
                            ConstMatrixView vmr_list_int,
                            ConstTensor3View ext_mat_int,
                            ConstMatrixView abs_vec_int,
                            ConstMatrixView sca_vec_int,
                            ConstMatrixView doit_i_field_mono_int,
                            ConstVectorView p_int,
                            const ArrayOfIndex& cloudbox_limits,
                            ConstVectorView f_grid,
                            const Index& f_index,
                            const Index& p_index,
                            const Index& lat_index,
                            const Index& lon_index, 
                            const Index& scat_za_index,
                            const Index& scat_aa_index,
                            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  const Index N_species = vmr_list_int.nrows();
  const Index stokes_dim = doit_i_field_mono.ncols();
  const Index atmosphere_dim = cloudbox_limits.nelem()/2;

  Vector sca_vec_av(stokes_dim,0);
  Vector stokes_vec(stokes_dim, 0.);
  Vector rtp_temperature_nlte_dummy(0);
  Vector rtp_vmr_local(N_species,0.); 

  // Two propmat_clearsky to average between
  Tensor4 cur_propmat_clearsky;
  Tensor4 prev_propmat_clearsky;

  Tensor3 ext_mat_local;
  Matrix abs_vec_local;  

  // Incoming stokes vector
  stokes_vec = doit_i_field_mono_int(joker, ppath_step.np-1);

  for( Index k = ppath_step.np-1; k >= 0; k--)
    {
      // Save propmat_clearsky from previous level by
      // swapping it with current level
      swap(cur_propmat_clearsky, prev_propmat_clearsky);

      //
      // Calculate scalar gas absorption
      //
      const Vector rtp_mag_dummy(3,0);
      const Vector ppath_los_dummy;
      
      Tensor3 nlte_dummy; //FIXME: do this right?
      ArrayOfTensor3 partial_dummy; // This is right since there should be only clearsky partials
      ArrayOfMatrix partial_source_dummy,partial_nlte_dummy; // This is right since there should be only clearsky partials
      propmat_clearsky_agendaExecute( ws, cur_propmat_clearsky,
                                      nlte_dummy,partial_dummy,partial_source_dummy,partial_nlte_dummy,
                                      ArrayOfRetrievalQuantity(0),
                                    f_grid[Range(f_index, 1)],
                                    rtp_mag_dummy, ppath_los_dummy,
                                    p_int[k], 
                                    t_int[k], 
                                    rtp_temperature_nlte_dummy,
                                    vmr_list_int(joker,k),
                                    propmat_clearsky_agenda );

      // Skip any further calculations for the first point.
      // We need values at two ppath points before we can average.
      if (k == ppath_step.np-1)
          continue;
    
      // Average prev_propmat_clearsky with cur_propmat_clearsky
      prev_propmat_clearsky += cur_propmat_clearsky;
      prev_propmat_clearsky *= 0.5;
        
      opt_prop_sum_propmat_clearsky(ext_mat_local, abs_vec_local,
                                    prev_propmat_clearsky);
        
      //
      // Add average particle extinction to ext_mat. 
      //
      for (Index i = 0; i < stokes_dim; i++)
        {
          for (Index j = 0; j < stokes_dim; j++)
            {
              ext_mat_local(0,i,j) += 0.5 *
                (ext_mat_int(i,j,k) + ext_mat_int(i,j,k+1));
            }
          //
          // Add average particle absorption to abs_vec.
          //
          abs_vec_local(0,i) += 0.5 * 
            (abs_vec_int(i,k) + abs_vec_int(i,k+1));
          
          //
          // Averaging of sca_vec:
          //
          sca_vec_av[i] =  0.5 *
            (sca_vec_int(i, k) + sca_vec_int(i, k+1));
                  
        }
      // Frequency
      Numeric f = f_grid[f_index];
      //
      // Calculate Planck function
      //
      Numeric rte_planck_value = planck(f, 0.5 * (t_int[k] + t_int[k+1]));
              
      // Length of the path between the two layers.
      Numeric lstep = ppath_step.lstep[k];
        
      // Some messages:
      out3 << "-----------------------------------------\n";
      out3 << "Input for radiative transfer step \n"
           << "calculation inside"
           << " the cloudbox:" << "\n";
      out3 << "Stokes vector at intersection point: \n" 
           << stokes_vec 
           << "\n"; 
      out3 << "lstep: ..." << lstep << "\n";
      out3 << "------------------------------------------\n";
      out3 << "Averaged coefficients: \n";
      out3 << "Planck function: " << rte_planck_value << "\n";
      out3 << "Scattering vector: " << sca_vec_av << "\n"; 
      out3 << "Absorption vector: " << abs_vec_local(0,joker) << "\n"; 
      out3 << "Extinction matrix: " << ext_mat_local(0,joker,joker) << "\n"; 
              
              
      assert (!is_singular( ext_mat_local(0,joker,joker)));
              
      // Radiative transfer step calculation. The Stokes vector
      // is updated until the considered point is reached.
      rte_step_doit(stokes_vec, Matrix(stokes_dim,stokes_dim), 
                   ext_mat_local(0,joker,joker), abs_vec_local(0,joker), 
                   sca_vec_av, lstep, rte_planck_value);
              
    }// End of loop over ppath_step. 
  // Assign calculated Stokes Vector to doit_i_field_mono.
  if (atmosphere_dim == 1)
    doit_i_field_mono(p_index - cloudbox_limits[0], 0, 0, scat_za_index, 0, joker)
      = stokes_vec;
  else if (atmosphere_dim == 3)
    doit_i_field_mono(p_index - cloudbox_limits[0],
                 lat_index - cloudbox_limits[2],
                 lon_index - cloudbox_limits[4],
                 scat_za_index, scat_aa_index,
                 joker) = stokes_vec;

}



//! cloud_RT_surface
/*
  This function calculates RT in the cloudbox if the next intersected 
  level is the surface. 

  CE (2006-05-29) Included surface_rtprop_agenda here.  

  \author Claudia Emde
  \date 2005-05-13
*/
void cloud_RT_surface(Workspace& ws,
                      //Output
                      Tensor6View doit_i_field_mono,
                      //Input
                      const Agenda& surface_rtprop_agenda,
                      ConstVectorView f_grid,
                      const Index& f_index,
                      const Index& stokes_dim,
                      const Ppath& ppath_step,
                      const ArrayOfIndex& cloudbox_limits, 
                      ConstVectorView scat_za_grid, 
                      const Index& scat_za_index
                     )
{
  chk_not_empty( "surface_rtprop_agenda", surface_rtprop_agenda );

  Matrix iy; 
  
  // Local output of surface_rtprop_agenda.
  Matrix surface_emission;
  Matrix surface_los; 
  Tensor4 surface_rmatrix;


  //Set rte_pos and rte_los to match the last point in ppath.
  
  Index np = ppath_step.np;
  
  Vector rte_pos;    // ppath_step.pos contains two columns for 1D
  rte_pos.resize( ppath_step.dim );
  rte_pos = ppath_step.pos(np-1,Range(0,ppath_step.dim));

  Vector rte_los; 
  rte_los.resize( ppath_step.los.ncols() );
  rte_los = ppath_step.los(np-1,joker);
  
  //Execute the surface_rtprop_agenda which gives the surface 
  //parameters.
  
  surface_rtprop_agendaExecute( ws, surface_emission, surface_los, 
                                surface_rmatrix, Vector(1,f_grid[f_index]),
                                rte_pos, rte_los, surface_rtprop_agenda );
  
  iy = surface_emission;

  Index nlos = surface_los.nrows();
  
  if( nlos > 0 )
    {
      Vector rtmp(stokes_dim); // Reflected Stokes vector for 1 frequency

      for( Index ilos=0; ilos<nlos; ilos++ )
        {
          // Several things needs to be fixed here. As far as I understand it,
          // this works only for specular cases and if the lower cloudbox limit
          // is exactly at the surface (PE, 120828)
  
          mult( rtmp, surface_rmatrix(ilos,0,joker,joker), 
                doit_i_field_mono( cloudbox_limits[0], 0, 0,
                      (scat_za_grid.nelem() -1 - scat_za_index), 0, joker) );
          iy(0,joker) += rtmp;
          
          doit_i_field_mono( cloudbox_limits[0], 0, 0, scat_za_index, 0, joker ) = 
                                                                iy( 0, joker );
        }
    }  
}

void doit_i_field_ngAcceleration(Tensor6& doit_i_field_mono,
                                 const ArrayOfTensor6& acceleration_input,
                                 const Index& accelerated,
                                 const Verbosity& )
{
    const Index N_p = doit_i_field_mono.nvitrines();
    const Index N_za = doit_i_field_mono.npages();
    
    // Loop over 4 components of Stokes Vector
    for (Index i = 0; i < accelerated; ++i)
    {
        ConstMatrixView S1 = acceleration_input[0](joker,0,0,joker,0,i);
        ConstMatrixView S2 = acceleration_input[1](joker,0,0,joker,0,i);
        ConstMatrixView S3 = acceleration_input[2](joker,0,0,joker,0,i);
        ConstMatrixView S4 = acceleration_input[3](joker,0,0,joker,0,i);
       
        ConstMatrixView J = S4;
        Matrix Q1;
        Matrix Q2;
        Matrix Q3;
        Numeric A1 = 0;
        Numeric A2B1 = 0;
        Numeric B2 = 0;
        Numeric C1 = 0;
        Numeric C2 = 0;
        Numeric NGA = 0;
        Numeric NGB = 0;
        
        // Q1 = -2*S3 + S4 + S2
        
        Q1 = S3;
        Q1 *= -2;
        Q1 += S4;
        Q1 += S2;
        
        // Q2 = S4 - S3 - S2 + S1
        Q2 = S4;
        Q2 -= S3;
        Q2 -= S2;
        Q2 += S1;
        
        //Q3 = S4 - S3
        Q3 = S4;
        Q3 -= S3;
        
        for ( Index p_index = 0; p_index < N_p; ++p_index)
        {
            for ( Index za_index = 0; za_index < N_za; ++za_index)
            {
                A1 += Q1(p_index,za_index) * Q1(p_index,za_index) * J(p_index,za_index);
                A2B1 += Q2(p_index,za_index) * Q1(p_index,za_index) * J(p_index,za_index);
                B2 += Q2(p_index,za_index) * Q2(p_index,za_index) * J(p_index,za_index);
                C1 += Q1(p_index,za_index) * Q3(p_index,za_index) * J(p_index,za_index);
                C2 += Q2(p_index,za_index) * Q3(p_index,za_index) * J(p_index,za_index);
            }
        }
        
        NGA = (C1*B2 - C2*A2B1) / (A1*B2 - A2B1*A2B1);
        NGB = (C2*A1 - C1*A2B1) / (A1*B2 - A2B1*A2B1);
        

        // Calculating the accelerated field
        for ( Index p_index = 0; p_index < N_p; ++p_index)
        {
            for ( Index za_index = 0; za_index < N_za; ++za_index)
            {
                Q1(p_index,za_index) = (1-NGA-NGB)*S4(p_index,za_index)
                                      + NGA*S3(p_index,za_index) + NGB*S2(p_index,za_index);
            }
        }
        doit_i_field_mono(joker,0,0,joker,0,i) = Q1;
    }
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
                          MatrixView doit_i_field_mono_int,
                          VectorView t_int, 
                          MatrixView vmr_list_int,
                          VectorView p_int,
                          //Input
                          ConstTensor5View ext_mat_field, 
                          ConstTensor4View abs_vec_field,
                          ConstTensor6View doit_scat_field,
                          ConstTensor6View doit_i_field_mono,
                          ConstTensor3View t_field, 
                          ConstTensor4View vmr_field, 
                          ConstVectorView p_grid,
                          const Ppath& ppath_step,
                          const ArrayOfIndex& cloudbox_limits,
                          ConstVectorView scat_za_grid,
                          const Index& scat_za_interp,
                          const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  // Stokes dimension
  const Index stokes_dim = doit_i_field_mono.ncols();
  
  // Gridpositions inside the cloudbox.
  // The optical properties are stored only inside the
  // cloudbox. For interpolation we use grids
  // inside the cloudbox.
  ArrayOfGridPos cloud_gp_p = ppath_step.gp_p;
  
  for(Index i = 0; i < ppath_step.np; i++ )
    cloud_gp_p[i].idx -= cloudbox_limits[0];
  
  // Grid index for points at upper limit of cloudbox must be shifted
  const Index n1 = cloudbox_limits[1] - cloudbox_limits[0];
  gridpos_upperend_check( cloud_gp_p[0], n1 );
  gridpos_upperend_check( cloud_gp_p[ppath_step.np-1], n1 );

  
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
      
      out3 << "Interpolate doit_scat_field and doit_i_field_mono:\n";
      if(scat_za_interp == 0) // linear interpolation
        {
          interp( sca_vec_int(i, joker), itw_p_za, 
                  doit_scat_field(joker, 0, 0, joker, 0, i), 
                  cloud_gp_p, gp_za);
          interp( doit_i_field_mono_int(i, joker), itw_p_za, 
                  doit_i_field_mono(joker, 0, 0, joker, 0, i), cloud_gp_p, gp_za);
        }
      else if (scat_za_interp == 1) //polynomial interpolation
        {
          // These intermediate variables are needed because polynomial 
          // interpolation is not implemented as multidimensional 
          // interpolation.
          Tensor3 sca_vec_int_za(stokes_dim, ppath_step.np,
                                 scat_za_grid.nelem(), 0.);
          Tensor3 doit_i_field_mono_int_za(stokes_dim, ppath_step.np, 
                                      scat_za_grid.nelem(), 0.);
          for (Index za = 0; za < scat_za_grid.nelem(); za++)
            {
              interp( sca_vec_int_za(i, joker, za), itw, 
                      doit_scat_field(joker, 0, 0, za, 0, i), cloud_gp_p);
              out3 << "Interpolate doit_i_field_mono:\n";
              interp( doit_i_field_mono_int_za(i, joker, za), itw, 
                      doit_i_field_mono(joker, 0, 0, za, 0, i), cloud_gp_p);
            }
          for (Index ip = 0; ip < ppath_step.np; ip ++)
            {
              sca_vec_int(i, ip) = 
                interp_poly(scat_za_grid, sca_vec_int_za(i, ip, joker),
                            los_grid[ip], gp_za[ip]);
              doit_i_field_mono_int(i, ip) = 
                interp_poly(scat_za_grid, doit_i_field_mono_int_za(i, ip, joker),
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



//! Radiative transfer calculation inside cloudbox for planeparallel case.
/*! 
  This function calculates the radiation field along a line of sight. This 
  function is used for the sequential update of the radiation field and 
  called inside a loop over the pressure grid. 
  
  The function gets all the atmospheric points on the pressure grid. 
  Then a radiative transfer step is 
  performed using the stokes vector as output and input. The inermediate
  Stokes vectors are stored in the WSV doit_i_field_mono.

  \param[in,out] ws Current Workspace
 WS Output:
  \param[out] doit_i_field_mono Updated radiation field inside the cloudbox.
  Variables used in opt_prop_xxx_agenda:
  WS Input:
  \param p_index // Pressure index
  \param scat_za_index // Index for proagation direction
  \param scat_za_grid
  \param cloudbox_limits 
  \param doit_scat_field Scattered field.
  Calculate scalar gas absorption:
  \param propmat_clearsky_agenda
  \param vmr_field
  Propagation path calculation:
  \param ppath_step_agenda
  \param p_grid
  \param z_field
  \param refellipsoid
  Calculate thermal emission:
  \param t_field
  \param f_grid
  \param f_index
  \param ext_mat_field
  \param abs_vec_field

  \author Sreerekha Ravi
  \date 2003-11-17
*/
void cloud_ppath_update1D_planeparallel(Workspace& ws,
                                        Tensor6View doit_i_field_mono,
                                        const Index& p_index,
                                        const Index& scat_za_index,
                                        ConstVectorView scat_za_grid,
                                        const ArrayOfIndex& cloudbox_limits,
                                        ConstTensor6View doit_scat_field,
                                        // Calculate scalar gas absorption:
                                        const Agenda& propmat_clearsky_agenda,
                                        ConstTensor4View vmr_field,
                                        // Propagation path calculation:
                                        ConstVectorView p_grid,
                                        ConstTensor3View z_field,
                                        // Calculate thermal emission:
                                        ConstTensor3View t_field,
                                        ConstVectorView f_grid,
                                        const Index& f_index,
                                        //particle opticla properties
                                        ConstTensor5View ext_mat_field,
                                        ConstTensor4View abs_vec_field,
                                        const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  const Index N_species = vmr_field.nbooks();
  const Index stokes_dim = doit_i_field_mono.ncols();
  const Index atmosphere_dim = 1;   
  Tensor4 propmat_clearsky;
  Tensor3 ext_mat;
  Matrix abs_vec;
  Vector rtp_vmr(N_species,0.), rtp_temperature_nlte_dummy(0); 
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
              stokes_vec = doit_i_field_mono(p_index-cloudbox_limits[0] +1, 0, 0, scat_za_index, 0, joker);         
              Numeric z_field_above = z_field(p_index +1, 0, 0);
              Numeric z_field_0 = z_field(p_index, 0, 0);
             
              Numeric cos_theta, lstep;
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
              
              lstep =  abs(dz/cos_theta) ;
              
              // Average temperature
              Numeric rtp_temperature =   0.5 * (t_field(p_index,0,0)
                                                 + t_field(p_index + 1,0,0));
              //
              // Average pressure
              Numeric rtp_pressure = 0.5 * (p_grid[p_index]
                                            + p_grid[p_index + 1]);
              
              // Average vmrs
              for (Index j = 0; j < N_species; j++)
                rtp_vmr[j] = 0.5 * (vmr_field(j,p_index,0,0) + 
                                         vmr_field(j,p_index + 1,0,0));
              //
              // Calculate scalar gas absorption and add it to abs_vec 
              // and ext_mat.
              //

              const Vector rtp_mag_dummy(3,0);
              const Vector ppath_los_dummy;

              Tensor3 nlte_dummy; //FIXME: do this right?
              ArrayOfTensor3 partial_dummy; // This is right since there should be only clearsky partials
              ArrayOfMatrix partial_source_dummy,partial_nlte_dummy; // This is right since there should be only clearsky partials
              propmat_clearsky_agendaExecute(ws, propmat_clearsky, 
                                             nlte_dummy,partial_dummy,
                                             partial_source_dummy,
                                             partial_nlte_dummy,
                                             ArrayOfRetrievalQuantity(0),
                                             f_grid[Range(f_index, 1)],
                                             rtp_mag_dummy,ppath_los_dummy,
                                             rtp_pressure,
                                             rtp_temperature,
                                             rtp_temperature_nlte_dummy,
                                             rtp_vmr,
                                             propmat_clearsky_agenda);
              
              opt_prop_sum_propmat_clearsky(ext_mat, abs_vec, propmat_clearsky);
              
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
              Numeric rte_planck_value = planck(f, rtp_temperature);
              
              // Some messages:
              out3 << "-----------------------------------------\n";
              out3 << "Input for radiative transfer step \n"
                   << "calculation inside"
                   << " the cloudbox:" << "\n";
              out3 << "Stokes vector at intersection point: \n" 
                   << stokes_vec 
                   << "\n"; 
              out3 << "lstep: ..." << lstep << "\n";
              out3 << "------------------------------------------\n";
              out3 << "Averaged coefficients: \n";
              out3 << "Planck function: " << rte_planck_value << "\n";
              out3 << "Scattering vector: " << sca_vec_av << "\n"; 
              out3 << "Absorption vector: " << abs_vec(0,joker) << "\n"; 
              out3 << "Extinction matrix: " << ext_mat(0,joker,joker) << "\n"; 
              
              
              assert (!is_singular( ext_mat(0,joker,joker)));
              
              // Radiative transfer step calculation. The Stokes vector
              // is updated until the considered point is reached.
              rte_step_doit(stokes_vec, Matrix(stokes_dim,stokes_dim), 
                           ext_mat(0,joker,joker), abs_vec(0,joker), 
                           sca_vec_av, lstep, rte_planck_value);
              
              // Assign calculated Stokes Vector to doit_i_field_mono.
              doit_i_field_mono(p_index - cloudbox_limits[0],
                      0, 0,
                      scat_za_index, 0,
                      joker) = stokes_vec;
            }
          if(scat_za_grid[scat_za_index] > 90)
            {
              stokes_vec = doit_i_field_mono(p_index-cloudbox_limits[0] -1, 0, 0, scat_za_index, 0, joker);
              Numeric z_field_below = z_field(p_index -1, 0, 0);
              Numeric z_field_0 = z_field(p_index, 0, 0);
              
              Numeric cos_theta, lstep;
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
              lstep =  abs(dz/cos_theta) ;
              
              // Average temperature
              Numeric rtp_temperature =   0.5 * (t_field(p_index,0,0)
                                                 + t_field(p_index - 1,0,0));
              //
              // Average pressure
              Numeric rtp_pressure = 0.5 * (p_grid[p_index]
                                            + p_grid[p_index - 1]);
              
              //
              // Average vmrs
              for (Index k = 0; k < N_species; k++)
                rtp_vmr[k] = 0.5 * (vmr_field(k,p_index,0,0) + 
                                         vmr_field(k,p_index - 1,0,0));
              //
              // Calculate scalar gas absorption and add it to abs_vec 
              // and ext_mat.
              //
                
              const Vector rtp_mag_dummy(3,0);
              const Vector ppath_los_dummy;
              Tensor3 nlte_dummy; //FIXME: do this right?
              ArrayOfTensor3 partial_dummy; // This is right since there should be only clearsky partials
              ArrayOfMatrix partial_source_dummy,partial_nlte_dummy; // This is right since there should be only clearsky partials
              propmat_clearsky_agendaExecute( ws, propmat_clearsky,
                                              nlte_dummy,partial_dummy,
                                              partial_source_dummy,
                                              partial_nlte_dummy,
                                              ArrayOfRetrievalQuantity(0),
                                              f_grid[Range(f_index, 1)],
                                              rtp_mag_dummy,ppath_los_dummy,
                                              rtp_pressure, 
                                              rtp_temperature, 
                                              rtp_temperature_nlte_dummy,
                                              rtp_vmr,
                                              propmat_clearsky_agenda );

              opt_prop_sum_propmat_clearsky(ext_mat, abs_vec, propmat_clearsky);

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
              Numeric rte_planck_value = planck(f, rtp_temperature);
              
              // Some messages:
              out3 << "-----------------------------------------\n";
              out3 << "Input for radiative transfer step \n"
                   << "calculation inside"
                   << " the cloudbox:" << "\n";
              out3 << "Stokes vector at intersection point: \n" 
                   << stokes_vec 
                   << "\n"; 
              out3 << "lstep: ..." << lstep << "\n";
              out3 << "------------------------------------------\n";
              out3 << "Averaged coefficients: \n";
              out3 << "Planck function: " << rte_planck_value << "\n";
              out3 << "Scattering vector: " << sca_vec_av << "\n"; 
              out3 << "Absorption vector: " << abs_vec(0,joker) << "\n"; 
              out3 << "Extinction matrix: " << ext_mat(0,joker,joker) << "\n"; 
              
              
              assert (!is_singular( ext_mat(0,joker,joker)));
              
              // Radiative transfer step calculation. The Stokes vector
              // is updated until the considered point is reached.
              rte_step_doit(stokes_vec, Matrix(stokes_dim,stokes_dim), 
                           ext_mat(0,joker,joker), abs_vec(0,joker), 
                           sca_vec_av, lstep, rte_planck_value);
              
              // Assign calculated Stokes Vector to doit_i_field_mono.
              doit_i_field_mono(p_index - cloudbox_limits[0],
                      0, 0,
                      scat_za_index, 0,
                      joker) = stokes_vec;
            }
          
        
        }// if loop end - for non_ground background
      
      // bkgr=2 indicates that the background is the surface
      else if (bkgr == 2)
        {
          //Set rte_pos, rte_gp_p and rte_los to match the last point
          //in ppath.
          //pos
          Vector rte_pos( atmosphere_dim );
          Numeric z_field_0 = z_field(0, 0, 0);
          rte_pos = z_field_0;  //ppath_step.pos(np-1,Range(0,atmosphere_dim));
          //los
          Vector rte_los(1);
          rte_los = scat_za_grid[scat_za_index];//ppath_step.los(np-1,joker);
          //gp_p
          //GridPos rte_gp_p;
          //rte_gp_p.idx   = p_index;
          //rte_gp_p.fd[0] = 0;
          //rte_gp_p.fd[1] = 1;
          //gridpos_copy( rte_gp_p, ppath_step.gp_p[np-1] ); 
          // Executes the surface agenda
          // FIXME: Convert to new agenda scheme before using
          // surface_rtprop_agenda.execute();

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

          // Define a local vector doit_i_field_mono_sum which adds the
          // products of groudnd_refl_coeffs with the downwelling 
          // radiation for each elements of surface_los
          Vector doit_i_field_mono_sum(stokes_dim,0);
          // Loop over the surface_los elements
          for( Index ilos=0; ilos < nlos; ilos++ )
            {
              if( stokes_dim == 1 )
                {
                  doit_i_field_mono_sum[0] += surface_refl_coeffs(ilos,f_index,0,0) *
                    doit_i_field_mono(cloudbox_limits[0],
                            0, 0,
                            (scat_za_grid.nelem() -1 - scat_za_index), 0,
                            0);
                }
              else 
                {
                  Vector stokes_vec2(stokes_dim);
                  mult( stokes_vec2, 
                        surface_refl_coeffs(ilos,0,joker,joker), 
                        doit_i_field_mono(cloudbox_limits[0],
                                0, 0,
                                (scat_za_grid.nelem() -1 - scat_za_index), 0,
                                joker));
                  for( Index is=0; is < stokes_dim; is++ )
                    { 
                      doit_i_field_mono_sum[is] += stokes_vec2[is];
                    }
                  
                }
            }
          // Copy from *doit_i_field_mono_sum* to *doit_i_field_mono*, and add the surface emission
          for( Index is=0; is < stokes_dim; is++ )
            {
              doit_i_field_mono (cloudbox_limits[0],
                       0, 0,
                       scat_za_index, 0,
                       is) = doit_i_field_mono_sum[is] + surface_emission(f_index,is);
            }
         
          //cout<<"scat_za_grid"<<scat_za_grid[scat_za_index]<<endl;
          //cout<<"p_index"<<p_index<<endl;
          //cout<<"cloudboxlimit[0]"<<cloudbox_limits[0]<<endl;
          // now the RT is done to the next point in the path.
          // 
          Vector stokes_vec_local;
          stokes_vec_local = doit_i_field_mono (0,
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
          
          Numeric lstep =  abs(dz/cos_theta) ;
          
          // Average temperature
          Numeric rtp_temperature =   0.5 * (t_field(p_index,0,0)
                                             + t_field(p_index + 1,0,0));
          
          //
          // Average pressure
          Numeric rtp_pressure = 0.5 * (p_grid[p_index]
                                        + p_grid[p_index + 1]);
          
          //
          const Index N_species = vmr_field.nbooks();
          // Average vmrs
          for (Index k = 0; k < N_species; k++)
            {
              rtp_vmr[k] = 0.5 * (vmr_field(k,p_index,0,0) + 
                                       vmr_field(k,p_index + 1,0,0));
            }
          //
          // Calculate scalar gas absorption and add it to abs_vec 
          // and ext_mat.
          //
          
          // FIXME: Convert to new agenda scheme before using
          //abs_scalar_gas_agenda.execute(p_index);
          
          opt_prop_gas_agendaExecute(ext_mat, abs_vec, abs_scalar_gas,
                                     opt_prop_gas_agenda);
          
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
          Numeric rte_planck_value = planck(f, rtp_temperature);
          
          assert (!is_singular( ext_mat(0,joker,joker)));
          
          // Radiative transfer step calculation. The Stokes vector
          // is updated until the considered point is reached.
          rte_step_doit(stokes_vec_local, ext_mat(0,joker,joker), 
                   abs_vec(0,joker), 
                   sca_vec_av, lstep, rte_planck_value);
          // Assign calculated Stokes Vector to doit_i_field_mono.
          doit_i_field_mono(p_index - cloudbox_limits[0],
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
  \param doit_i_field_mono Radiation field calculated on a very fine za grid.
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
                ConstTensor6View doit_i_field_mono,
                const Numeric& acc,
                const Index& scat_za_interp)
{
  Index N_za = za_grid_fine.nelem();

  assert(doit_i_field_mono.npages() == N_za);
  
  Index N_p = doit_i_field_mono.nvitrines();
  
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

      // Interpolate reduced intensity field on fine za_grid for 
      // all pressure levels
      for( Index i_p = 0; i_p < N_p; i_p++ )
        {
          for( Index i_za_red = 0; i_za_red < idx.nelem(); i_za_red ++)
            {
              za_reduced[i_za_red] = za_grid_fine[idx[i_za_red]];
              doit_i_field_opt(i_p, i_za_red) = doit_i_field_mono(i_p, 0, 0, idx[i_za_red], 
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
              diff_vec[i_za]  =  abs( doit_i_field_mono(i_p, 0, 0, i_za, 0 ,0)
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
      max_diff = max_diff_p/doit_i_field_mono(ind_p, 0, 0, ind_za[ind_p], 0, 0)*100.;
      
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
      doit_i_field_opt(joker, i) = doit_i_field_mono(joker, 0, 0, idx[i], 0, 0);
    }
}
          

//! Interpolation of cloud box intensity field
/* 
   See WSM *iyInterpCloudboxField*.

   \param iy                Out: As the WSV with same name.
   \param doit_i_field      In: As the WSV with same name.
   \param rte_gp_p          In: As the WSV with same name.
   \param rte_gp_lat        In: As the WSV with same name.
   \param rte_gp_lon        In: As the WSV with same name.
   \param rte_los           In: As the WSV with same name.
   \param cloudbox_on       In: As the WSV with same name.
   \param cloudbox_limits   In: As the WSV with same name.
   \param atmosphere_dim    In: As the WSV with same name.
   \param stokes_dim        In: As the WSV with same name.
   \param scat_za_grid      In: As the WSV with same name. 
   \param scat_aa_grid      In: As the WSV with same name. 
   \param f_grid            In: As the WSV with same name.
   \param p_grid            In: As the WSV with same name.
   \param interpmeth        Interpolation method. Can be "linear" or 
   "polynomial".
   \param rigorous          Fail if incoming field is not safely interpolable.
   \param maxratio          Maximum allowed ratio of two radiances regarded as
   interpolable.

   \author Claudia Emde and Patrick Eriksson
   \date 2004-09-29
*/
void iy_interp_cloudbox_field(Matrix&               iy,
                              const Tensor7&        doit_i_field,
                              const GridPos&        rte_gp_p,
                              const GridPos&        rte_gp_lat,
                              const GridPos&        rte_gp_lon,
                              const Vector&         rte_los,
                              const Index&          cloudbox_on,
                              const ArrayOfIndex&   cloudbox_limits,
                              const Index&          atmosphere_dim,
                              const Index&          stokes_dim,
                              const Vector&         scat_za_grid,
                              const Vector&         scat_aa_grid,
                              const Vector&         f_grid,
                              const String&         interpmeth,
                              const Index&          rigorous,
                              const Numeric&        maxratio,
                              const Verbosity&      verbosity)
{
  CREATE_OUT3;
  
  //--- Check input -----------------------------------------------------------
  if( !(atmosphere_dim == 1  ||  atmosphere_dim == 3) )
    throw runtime_error( "The atmospheric dimensionality must be 1 or 3.");
  if( !cloudbox_on )
    throw runtime_error( "The cloud box is not activated and no outgoing "
                         "field can be returned." );
  if ( cloudbox_limits.nelem() != 2*atmosphere_dim )
    throw runtime_error(
       "*cloudbox_limits* is a vector which contains the upper and lower\n"
       "limit of the cloud for all atmospheric dimensions.\n"
       "So its length must be 2 x *atmosphere_dim*" ); 
  if( scat_za_grid.nelem() == 0 )
    throw runtime_error( "The variable *scat_za_grid* is empty. Are dummy "
                         "values from *cloudboxOff used?" );
  if( !( interpmeth == "linear"  ||  interpmeth == "polynomial" ) )
    throw runtime_error( "Unknown interpolation method. Possible choices are "
                         "\"linear\" and \"polynomial\"." );
  if( interpmeth == "polynomial"  &&  atmosphere_dim != 1  )
    throw runtime_error( "Polynomial interpolation method is only available "
                         "for *atmosphere_dim* = 1." );
  if( doit_i_field.nlibraries() != f_grid.nelem() )
    throw runtime_error( "Inconsistency in size between f_grid and doit_i_field! "
         "(This method does not yet handle dispersion type calculations.)" );
  //---------------------------------------------------------------------------


  //--- Determine if at border or inside of cloudbox (or outside!)
  //
  // Let us introduce a number coding for cloud box borders.
  // Borders have the same number as position in *cloudbox_limits*.
  // Inside cloud box is coded as 99, and outside as > 100.
  Index  border  = 999;
  //
  //- Check if at any border
  if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[0], false ) )
    { border = 0; }
  else if( is_gridpos_at_index_i( rte_gp_p, cloudbox_limits[1], false ) )
    { border = 1; }
  if( atmosphere_dim == 3  &&  border > 100 )
    {
      if( is_gridpos_at_index_i( rte_gp_lat, cloudbox_limits[2], false ) )
        { border = 2; }
      else if( is_gridpos_at_index_i( rte_gp_lat, cloudbox_limits[3], false ) )
        { border = 3; }
      else if( is_gridpos_at_index_i( rte_gp_lon, cloudbox_limits[4], false ) )
        { border = 4; }
      else if( is_gridpos_at_index_i( rte_gp_lon, cloudbox_limits[5], false ) )
        { border = 5; }
    }

  //
  //- Check if inside (till here border<100 means we are at some border)
  if( border > 100 )
    {
      // Assume inside as it is easiest to detect if outside (can be detected
      // check in one dimension at the time)
      bool inside = true;
      Numeric fgp;

      // Check in pressure dimension
      fgp = fractional_gp( rte_gp_p );
      if( fgp < Numeric(cloudbox_limits[0])  || 
          fgp > Numeric(cloudbox_limits[1]) )
        { inside = false; }

      // Check in lat and lon dimensions
     if( atmosphere_dim == 3  &&  inside )
       {
         fgp = fractional_gp( rte_gp_lat );
         if( fgp < Numeric(cloudbox_limits[2])  || 
             fgp > Numeric(cloudbox_limits[3]) )
           { inside = false; }
         fgp = fractional_gp( rte_gp_lon );
         if( fgp < Numeric(cloudbox_limits[4])  || 
             fgp > Numeric(cloudbox_limits[5]) )
           { inside = false; }
       }

     if( inside )
       { border = 99; }
    }

  // If outside, something is wrong
  if( border > 100 )
    {
      throw runtime_error( 
                 "Given position has been found to be outside the cloudbox." );
    }

  const Index border_index = border ? doit_i_field.nvitrines()-1 : 0;

  //- Sizes
  const Index   nf  = f_grid.nelem();
  DEBUG_ONLY (const Index   np  = cloudbox_limits[1] - cloudbox_limits[0] + 1);
  const Index   nza  = scat_za_grid.nelem();

  //- Resize *iy*
  iy.resize( nf, stokes_dim );


  // Sensor points inside the cloudbox
  if( border == 99 )
    {
      if (atmosphere_dim == 3)
        {
          throw runtime_error(
                              "3D DOIT calculations are not implemented\n"
                              "for observations from inside the cloudbox.\n"
                              );
        }
      else
        {
          assert(atmosphere_dim == 1);
          
          // *doit_i_field* is normally calculated internally:
          assert( is_size(doit_i_field, nf, np, 1, 1, nza, 1, stokes_dim) );
          
          out3 << "    Interpolating outgoing field:\n";
          out3 << "       zenith_angle: " << rte_los[0] << "\n";
          out3 << " Sensor inside cloudbox at position:  " << 
            rte_gp_p << "\n";
          
          // Grid position in *scat_za_grid*
          GridPos gp_za;
          gridpos( gp_za, scat_za_grid, rte_los[0] );
          
          // Corresponding interpolation weights
          Vector itw_za(2);
          interpweights( itw_za, gp_za );
          
          // Grid position in *p_grid* (only cloudbox part because 
          // doit_i_field1D_spectra is only defined inside the cloudbox
          GridPos gp_p;
          gp_p = rte_gp_p;
          gp_p.idx = rte_gp_p.idx - cloudbox_limits[0]; 
          gridpos_upperend_check( gp_p, cloudbox_limits[1] - 
                                        cloudbox_limits[0] );
          
          //          cout << gp_p << endl;

          Vector itw_p(2);
          interpweights( itw_p, gp_p );

          Vector iy_p(nza);

          if( interpmeth == "linear" )
            {
              for(Index is = 0; is < stokes_dim; is++ )
                {
                  for(Index iv = 0; iv < nf; iv++ )
                    {
                      for (Index i_za = 0; i_za < nza; i_za++)
                        {
                          iy_p[i_za] = interp
                          (itw_p, doit_i_field(iv, joker, 0, 0, i_za, 0, is),
                             gp_p);
                        }
                      iy(iv,is) = interp( itw_za, iy_p, gp_za);
                    }
                }
            }
          else
            {   
              for(Index is = 0; is < stokes_dim; is++ )
                {
                  for(Index iv = 0; iv < nf; iv++ )
                    {
                      for (Index i_za = 0; i_za < nza; i_za++)
                        {
                          iy_p[i_za] = interp
                          (itw_p, doit_i_field(iv, joker, 0, 0, i_za, 0, is),
                             gp_p);
                        }
                      iy(iv,is) =  interp_poly( scat_za_grid, iy_p, rte_los[0],
                                                gp_za );
                    }
                }
            }
          
        }
      
    }
  
  // Sensor outside the cloudbox

  // --- 1D ------------------------------------------------------------------
  else if( atmosphere_dim == 1 )
    {
      out3 << "    Interpolating outgoing field:\n";
      out3 << "       zenith_angle: " << rte_los[0] << "\n";
      if( border )
        out3 << "       top side\n";
      else
        out3 << "       bottom side\n";

      // Grid position in *scat_za_grid*
      GridPos gp;
      gridpos( gp, scat_za_grid, rte_los[0] );

      // Corresponding interpolation weights
      Vector itw(2);
      interpweights( itw, gp );

      if( interpmeth == "linear" )
        {
          for(Index is = 0; is < stokes_dim; is++ )
            {
              if( is==0 && rigorous )
                {
                  for(Index iv = 0; iv < nf; iv++ )
                    {
                      if( doit_i_field(iv,border_index,0,0,gp.idx,0,is)/
                          doit_i_field(iv,border_index,0,0,gp.idx+1,0,is) > 1/maxratio &&
                          doit_i_field(iv,border_index,0,0,gp.idx,0,is)/
                          doit_i_field(iv,border_index,0,0,gp.idx+1,0,is) < maxratio )
                        {
                          iy(iv,is) = interp( itw, 
                                              doit_i_field( iv, border_index, 0, 0, joker, 0, is ),
                                              gp );
                        }
                      else
                        {
                          ostringstream os;
                          os << "ERROR: Radiance difference between interpolation\n"
                             << "points is too large (factor " << maxratio << ") to\n"
                             << "safely interpolate. This might be due to za_grid\n"
                             << "being too coarse or the radiance field being a\n"
                             << "step-like function.\n";
                          os << "Happens at boundary " << border << " between zenith\n"
                             << "angels " << scat_za_grid[gp.idx] << " and "
                             << scat_za_grid[gp.idx+1] << "deg for frequency"
                             << "#" << iv << ", where radiances are "
                             << doit_i_field(iv,border_index,0,0,gp.idx,0,0)
                             << " and " << doit_i_field(iv,border_index,0,0,gp.idx+1,0,0)
                             << " W/(sr m2 Hz).";
                          throw runtime_error(os.str());
                        }
                    }
                }
              else
                {
                  for(Index iv = 0; iv < nf; iv++ )
                    {
                      iy(iv,is) = interp( itw, 
                                          doit_i_field( iv, border_index, 0, 0, joker, 0, is ),
                                          gp );
                    }                  
                }
            }
        }
      else
        {
          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp_poly( scat_za_grid, 
                       doit_i_field( iv, border_index, 0, 0, joker, 0, is ) , rte_los[0],
                                            gp );
                }
            }
        }
    }
  
  // --- 3D ------------------------------------------------------------------
  else
    {
      // Use asserts to check *scat_i_XXX* as these variables should to 99% be
      // calculated internally, and thus make it possible to avoid this check.
      assert ( is_size( doit_i_field, nf, doit_i_field.nvitrines(), doit_i_field.nshelves(),
                        doit_i_field.nbooks(), scat_za_grid.nelem(),
                        scat_aa_grid.nelem(), stokes_dim ));

      out3 << "    Interpolating outgoing field:\n";
      out3 << "       zenith angle : " << rte_los[0] << "\n";
      out3 << "       azimuth angle: " << rte_los[1]+180. << "\n";

      
      // Scattering angle grid positions
      GridPos gp_za, gp_aa;
      gridpos( gp_za, scat_za_grid, rte_los[0] );
      gridpos( gp_aa, scat_aa_grid, rte_los[1]+180. );

      // Interpolation weights (for 4D "red" interpolation)
      Vector   itw(16);

      // Outgoing from pressure level
      if( border <= 1 )
        {
          // Lat and lon grid positions with respect to cloud box 
          GridPos cb_gp_lat, cb_gp_lon;
          cb_gp_lat      = rte_gp_lat;
          cb_gp_lon      = rte_gp_lon;
          cb_gp_lat.idx -= cloudbox_limits[2];
          cb_gp_lon.idx -= cloudbox_limits[4]; 
          //          
          gridpos_upperend_check( cb_gp_lat, cloudbox_limits[3] - 
                                             cloudbox_limits[2] );
          gridpos_upperend_check( cb_gp_lon, cloudbox_limits[5] - 
                                             cloudbox_limits[4] );

          interpweights( itw, cb_gp_lat, cb_gp_lon, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                        doit_i_field( iv, border_index, joker, joker, joker, joker, is ),
                                      cb_gp_lat, cb_gp_lon, gp_za, gp_aa );
                }
            }
        }

      // Outgoing from latitude level
      else if( border <= 3 )
        {
          // Pressure and lon grid positions with respect to cloud box 
          GridPos cb_gp_p, cb_gp_lon;
          cb_gp_p        = rte_gp_p;
          cb_gp_lon      = rte_gp_lon;
          cb_gp_p.idx   -= cloudbox_limits[0];
          cb_gp_lon.idx -= cloudbox_limits[4]; 
          //          
          gridpos_upperend_check( cb_gp_p,   cloudbox_limits[1] - 
                                             cloudbox_limits[0] );
          gridpos_upperend_check( cb_gp_lon, cloudbox_limits[5] - 
                                             cloudbox_limits[4] );
          
          interpweights( itw, cb_gp_p, cb_gp_lon, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                    doit_i_field( iv, joker, border_index-2, joker, joker, joker, is ),
                                      cb_gp_p, cb_gp_lon, gp_za, gp_aa );
                }
            }
        }

      // Outgoing from longitude level
      else
        {
          // Pressure and lat grid positions with respect to cloud box 
          GridPos cb_gp_p, cb_gp_lat;
          cb_gp_p        = rte_gp_p;
          cb_gp_lat      = rte_gp_lat;
          cb_gp_p.idx   -= cloudbox_limits[0]; 
          cb_gp_lat.idx -= cloudbox_limits[2];
          //          
          gridpos_upperend_check( cb_gp_p,   cloudbox_limits[1] - 
                                             cloudbox_limits[0] );
          gridpos_upperend_check( cb_gp_lat, cloudbox_limits[3] - 
                                             cloudbox_limits[2] );
          
          interpweights( itw, cb_gp_p, cb_gp_lat, gp_za, gp_aa );

          for(Index is = 0; is < stokes_dim; is++ )
            {
              for(Index iv = 0; iv < nf; iv++ )
                {
                  iy(iv,is) = interp( itw, 
                     doit_i_field( iv, joker, joker, border_index-4, joker, joker, is ),
                                      cb_gp_p, cb_gp_lat, gp_za, gp_aa );
                }
            }
        }
    }
}


//! Normalization of scattered field
/*!
 Calculate the scattered extinction field and apply the
 derived correction factor to doit_scat_field.
 
 Only 1D is supported.
 
 \param[in,out] ws Current workspace
 \param[in,out] doit_scat_field Scattered field
 \param[in]     cloudbox_limits WS Input
 \param[in]     spt_calc_agenda WS Input
 \param[in]     atmopshere_dim WS Input
 \param[in]     scat_za_grid WS Input
 \param[in]     scat_aa_grid WS Input
 \param[in]     pnd_field WS Input
 \param[in]     opt_prop_part_agenda WS Input
 \param[in]     t_field WS Input
 \param[in]     norm_error_threshold  Normalization error threshold
 \param[in]     verbosity Verbosity
 
 \author Oliver Lemke
 \date 2013-01-17
 */
void
doit_scat_fieldNormalize(Workspace& ws,
                         Tensor6& doit_scat_field,
                         const Tensor6& doit_i_field_mono,
                         const ArrayOfIndex& cloudbox_limits,
                         const Agenda& spt_calc_agenda,
                         const Index& atmosphere_dim,
                         const Vector& scat_za_grid,
                         const Vector& scat_aa_grid,
                         const Tensor4& pnd_field,
                         const Agenda& opt_prop_part_agenda,
                         const Tensor3& t_field,
                         const Numeric& norm_error_threshold,
                         const Index& norm_debug,
                         const Verbosity& verbosity)
{
    if (atmosphere_dim != 1)
        throw runtime_error("Only 1D is supported here for now");

    CREATE_OUT0;
    CREATE_OUT2;

    // Number of zenith angles.
    const Index Nza = scat_za_grid.nelem();

    if (scat_za_grid[0] != 0. || scat_za_grid[Nza-1] != 180.)
        throw runtime_error("The range of *scat_za_grid* must [0 180].");

    // Number of azimuth angles.
    const Index Naa = scat_aa_grid.nelem();

    if (scat_aa_grid[0] != 0. || scat_aa_grid[Naa-1] != 360.)
        throw runtime_error("The range of *scat_aa_grid* must [0 360].");

    // Get stokes dimension from *doit_scat_field*:
    const Index stokes_dim = doit_scat_field.ncols();
    assert(stokes_dim > 0 || stokes_dim < 5);

    // To use special interpolation functions for atmospheric fields we
    // use ext_mat_field and abs_vec_field:
    Tensor5 ext_mat_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                          stokes_dim, stokes_dim, 0.);
    Tensor4 abs_vec_field(cloudbox_limits[1] - cloudbox_limits[0] + 1, 1, 1,
                          stokes_dim, 0.);

    const Index Np = doit_scat_field.nvitrines();

    Tensor5 doit_scat_ext_field(doit_scat_field.nvitrines(),
                                doit_scat_field.nshelves(),
                                doit_scat_field.nbooks(),
                                doit_scat_field.npages(),
                                doit_scat_field.nrows(),
                                0.);

    Index scat_aa_index_local = 0;

    // Calculate scattering extinction field
    for(Index scat_za_index_local = 0; scat_za_index_local < Nza;
        scat_za_index_local ++)
    {
        // This function has to be called inside the angular loop, as
        // spt_calc_agenda takes *scat_za_index* and *scat_aa_index*
        // from the workspace.
        // *scat_p_index* is needed for communication with agenda
        // *opt_prop_part_agenda*.
        cloud_fieldsCalc(ws, ext_mat_field, abs_vec_field,
                         spt_calc_agenda, opt_prop_part_agenda,
                         scat_za_index_local, scat_aa_index_local,
                         cloudbox_limits, t_field, pnd_field, verbosity);

        for(Index p_index = 0;
            p_index <= (cloudbox_limits[1] - cloudbox_limits[0]);
            p_index ++)
        {
            // For all in p_grid (in cloudbox):
            // I_ext = (ext_mat_field - abs_vec_field) * doit_i_field_mono
            // equivalent to:
            // I_ext = I * (K11-a1) + Q * (K12 - a2) + U * (K13 - a3) + V * (K14 - a4)
            for (Index i = 0; i < stokes_dim; i++)
            {
                doit_scat_ext_field(p_index, 0, 0, scat_za_index_local, 0)
                += doit_i_field_mono(p_index, 0, 0, scat_za_index_local, 0, i)
                * (ext_mat_field(p_index, 0, 0, 0, i) - abs_vec_field(p_index, 0, 0, i));
            }
        }
    }

    Numeric corr_max = .0;
    Index corr_max_p_index = -1;

    for (Index p_index = 0; p_index < Np; p_index++)
    {
        // Calculate scattering integrals
        const Numeric scat_int
        = AngIntegrate_trapezoid(doit_scat_field(p_index, 0, 0, joker, 0, 0),
                                 scat_za_grid);

        const Numeric scat_ext_int
        = AngIntegrate_trapezoid(doit_scat_ext_field(p_index, 0, 0, joker, 0),
                                 scat_za_grid);

        // Calculate factor between scattered extinction field integral
        // and scattered field integral
        const Numeric corr_factor = scat_ext_int / scat_int;

        // If no scattering is present, the correction factor can become
        // inf or nan. We just don't apply it for those cases.
        if (!isnan(corr_factor) && !isinf(corr_factor))
        {
            if (abs(corr_factor) > abs(corr_max))
            {
                corr_max = corr_factor;
                corr_max_p_index = p_index;
            }
            if (abs(1.-corr_factor) > norm_error_threshold)
            {
                ostringstream os;
                os <<   "ERROR: DOIT correction factor exceeds threshold: "
                << setprecision(4) <<  1.-corr_factor << " at p_index " << p_index << "\n";
                throw runtime_error(os.str());
            }
            else if (abs(1.-corr_factor) > norm_error_threshold/2.)
            {
                out0 << "  WARNING: DOIT correction factor above threshold/2: "
                << 1.-corr_factor << " at p_index " << p_index << "\n";
            }

            // Scale scattered field with correction factor
            doit_scat_field(p_index, 0, 0, joker, 0, joker) *= corr_factor;
        }
    }


    ArtsOut& norm_out = out2;
    if (norm_debug) norm_out = out0;

    if (corr_max_p_index != -1)
    {
        norm_out << "  Max. DOIT correction factor in this iteration: " << 1.-corr_max
        << " at p_index " << corr_max_p_index << "\n";
    }
    else
    {
        norm_out << "  No DOIT correction performed in this iteration.\n";
    }
}



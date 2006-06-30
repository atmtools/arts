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
  \file   rte.h
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2002-05-29

  \brief  Declaration of functions in rte.cc.
*/



#ifndef rte_h
#define rte_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "agenda_class.h"
#include "arts.h"
#include "complex.h"          
#include "ppath.h"
#include "matpackI.h"
#include "matpackIII.h"



/*===========================================================================
  === Functions in rte.cc
  ===========================================================================*/

void include_trans_in_diy_dq( 
            ArrayOfTensor4&   diy_dq,  
      const Index&            iv,
            bool              pol_trans,
      ConstMatrixView         trans,
      const ArrayOfPpath&     ppath_array, 
      const Index&            ppath_array_index );

void iy_calc(
              Matrix&                  iy,
              Ppath&                   ppath,
              Ppath&                   ppath_step,
              Index&                   ppath_array_index,
              ArrayOfPpath&            ppath_array,
              ArrayOfTensor4&          diy_dvmr,
              ArrayOfTensor4&          diy_dt,
        const Agenda&                  ppath_step_agenda,
        const Agenda&                  rte_agenda,
        const Agenda&                  iy_space_agenda,
        const Agenda&                  surface_prop_agenda,
        const Agenda&                  iy_cloudbox_agenda,
        const Index&                   atmosphere_dim,
        const Vector&                  p_grid,
        const Vector&                  lat_grid,
        const Vector&                  lon_grid,
        const Tensor3&                 z_field,
        const Tensor3&                 t_field,
        const Tensor4&                 vmr_field,
        const Matrix&                  r_geoid,
        const Matrix&                  z_surface,
        const Index&                   cloudbox_on, 
        const ArrayOfIndex&            cloudbox_limits,
        const Vector&                  pos,
        const Vector&                  los,
        const Vector&                  f_grid,
        const Index&                   stokes_dim,
        const Index&                   ppath_array_do,
        const ArrayOfIndex&            rte_do_vmr_jacs,
        const Index&                   rte_do_t_jacs,
        const bool&                    agenda_verb );

void iy_calc_no_jacobian(
              Matrix&         iy,
              Ppath&          ppath,
              Ppath&          ppath_step,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         iy_space_agenda,
        const Agenda&         surface_prop_agenda,
        const Agenda&         iy_cloudbox_agenda,
        const Index&          atmosphere_dim,
        const Vector&         p_grid,
        const Vector&         lat_grid,
        const Vector&         lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        const Tensor4&        vmr_field,
        const Matrix&         r_geoid,
        const Matrix&         z_surface,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Vector&         pos,
        const Vector&         los,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const bool&           agenda_verb );

void rte_step_std(
         //Output and Input:
         VectorView stokes_vec,
         MatrixView trans_mat,
         //Input
         ConstMatrixView ext_mat_av,
         ConstVectorView abs_vec_av,
         ConstVectorView sca_vec_av, 
         const Numeric& l_step,
         const Numeric& rte_planck_value );

void rte_std(
             Matrix&                  iy,
             Vector&                  emission,
             Matrix&                  abs_scalar_gas,
             Tensor4&                 ppath_transmissions,
             ArrayOfTensor4&          diy_dvmr,
             ArrayOfTensor4&          diy_dt,
       const Ppath&                   ppath,
       const ArrayOfPpath&            ppath_array, 
       const Index&                   ppath_array_index,
       const Vector&                  f_grid,
       const Index&                   stokes_dim,
       const Agenda&                  emission_agenda,
       const Agenda&                  abs_scalar_gas_agenda,
       const ArrayOfIndex&            rte_do_gas_jacs,
       const Index&                   rte_do_t_jacs,
       const bool&                    do_transmissions );

void surface_calc(
              Matrix&         iy,
        const Tensor3&        I,
        const Matrix&         surface_los,
        const Tensor4&        surface_rmatrix,
        const Matrix&         surface_emission );

void surface_specular_los(
              VectorView   los,
        const Index&       atmosphere_dim );

void surface_specular_R_and_b(
              MatrixView   surface_rmatrix,
              VectorView   surface_emission,
        const Complex&     Rv,
        const Complex&     Rh,
        const Numeric&     f,
        const Index&       stokes_dim,
        const Numeric&     surface_skin_t );


#endif  // rte_h

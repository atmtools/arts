/* Copyright (C) 2002-2008
   Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
   Stefan Buehler   <sbuehler@ltu.se>
                            
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
#include "jacobian.h"
#include "ppath.h"
#include "matpackI.h"
#include "matpackII.h"
#include "matpackIII.h"


/*===========================================================================
  === Functions in rte.cc
  ===========================================================================*/

void apply_y_unit( 
       MatrixView     iy, 
    const String&     y_unit, 
    ConstVectorView   f_grid );

void apply_y_unit2( 
    Tensor3View       J,
    ConstMatrixView   iy, 
    const String&     y_unit, 
    ConstVectorView   f_grid );

void get_ptvmr_for_ppath( 
        Vector&      ppath_p, 
        Vector&      ppath_t, 
        Matrix&      ppath_vmr, 
  const Ppath&       ppath,
  const Index&       atmosphere_dim,
  ConstVectorView    p_grid,
  ConstTensor3View   t_field,
  ConstTensor4View   vmr_field );

Range get_rowindex_for_mblock( 
  const Sparse&   sensor_response, 
  const Index&    imblock );

void get_step_vars_for_standardRT( 
        Workspace&   ws,
        Tensor3&     ppath_abs_scalar, 
        Matrix&      ppath_emission, 
        Matrix&      ppath_tau,
        Vector&      total_tau,
  const Agenda&      abs_scalar_agenda,
  const Agenda&      emission_agenda,
  const Ppath&       ppath,
  ConstVectorView    ppath_p, 
  ConstVectorView    ppath_t, 
  ConstMatrixView    ppath_vmr, 
  const Index&       nf,
  const Index&       emission_do );

void iy_transmission_for_scalar_tau( 
       Tensor3&     iy_transmission,
  const Index&      stokes_dim,
  ConstVectorView   tau );

void iy_transmission_mult( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstTensor3View   trans );

void iy_transmission_mult_scalar_tau( 
       Tensor3&      iy_trans_new,
  ConstTensor3View   iy_transmission,
  ConstVectorView    tau );

void iyb_calc(
        Workspace&                  ws,
        Vector&                     iyb,
        Vector&                     iyb_error,
        Index&                      iy_error_type,
        Matrix&                     iyb_aux,
        Index&                      n_aux,
        ArrayOfMatrix&              diyb_dx,
  const Index&                      imblock,
  const Index&                      atmosphere_dim,
  ConstVectorView                   p_grid,
  ConstVectorView                   lat_grid,
  ConstVectorView                   lon_grid,
  ConstTensor3View                  t_field,
  ConstTensor4View                  vmr_field,
  const Index&                      cloudbox_on,
  const Index&                      stokes_dim,
  ConstVectorView                   f_grid,
  ConstMatrixView                   sensor_pos,
  ConstMatrixView                   sensor_los,
  ConstVectorView                   mblock_za_grid,
  ConstVectorView                   mblock_aa_grid,
  const Index&                      antenna_dim,
  const Agenda&                     iy_clearsky_agenda,
  const Index&                      iy_aux_do,
  const String&                     y_unit,
  const Index&                      j_analytical_do,
  const ArrayOfRetrievalQuantity&   jacobian_quantities,
  const ArrayOfArrayOfIndex&        jacobian_indices );

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

void surface_calc(
              Matrix&         iy,
        const Tensor3&        I,
        const Matrix&         surface_los,
        const Tensor4&        surface_rmatrix,
        const Matrix&         surface_emission );



#endif  // rte_h

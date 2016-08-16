/* Copyright (C) 2006-2012 Claudia Emde <claudia.emde@dlr.de>

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

/**
 \file   disort.h
 \author Claudia Emde <claudia.emde@dlr.de>
 \date   Tue Feb  7 11:48:17 2006
  
 \brief  Functions for disort interface
 * 
 * 
 */

#ifndef disort_h
#define disort_h

#include "agenda_class.h"
#include "matpackIV.h"
#include "mystring.h"
#include "optproperties.h"


void dtauc_ssalbCalc( Workspace &ws,
                      VectorView dtauc,
                      VectorView ssalb,
                      const Agenda& propmat_clearsky_agenda,
                      const Agenda& spt_calc_agenda,
                      const Agenda& opt_prop_part_agenda,
                      ConstTensor4View pnd_field,
                      ConstTensor3View t_field,
                      ConstTensor3View z_field, 
                      ConstTensor4View vmr_field,
                      ConstVectorView p_grid,
                      const ArrayOfIndex& cloudbox_limits,
                      ConstVectorView f_mono );

void phase_functionCalc2( Workspace& ws,
                          //Output
                          MatrixView phase_function,
                          //Input
                          const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                          const Agenda& spt_calc_agenda,
                          const Agenda& opt_prop_part_agenda,
                          ConstTensor4View pnd_field,
                          ConstTensor3View t_field,
                          const ArrayOfIndex& cloudbox_limits,
                          const Index& pfct_za_grid_size,
                          const Verbosity& verbosity );

void phase_functionCalc( //Output
                         MatrixView phase_function,
                         //Input
                         const ArrayOfArrayOfSingleScatteringData& scat_data_mono,
                         ConstTensor4View pnd_field,
                         const ArrayOfIndex& cloudbox_limits,
                         const String pfct_method );

void pmomCalc2( //Output
                MatrixView pmom,
                //Input
                ConstMatrixView phase_function, 
                ConstVectorView scat_angle_grid,
                const Index n_legendre,
                const Verbosity& verbosity );

void pmomCalc( //Output
               MatrixView pmom,
               //Input
               ConstMatrixView phase_function, 
               ConstVectorView scat_angle_grid,
               const Index n_legendre,
               const Verbosity& verbosity );

#ifdef ENABLE_DISORT
void get_cb_inc_field( Workspace&      ws,
                       Matrix&         cb_inc_field,
                       const Agenda&   iy_main_agenda,
                       const Tensor3&  z_field,
                       const Tensor3&  t_field,
                       const Tensor4&  vmr_field,
                       const ArrayOfIndex&   cloudbox_limits,
                       const Vector&   f_grid,
                       const Vector&   scat_za_grid,
                       const Index&    nstreams );
#else /* ENABLE_DISORT */
void get_cb_inc_field( Workspace&,
                       Matrix&,
                       const Agenda&,
                       const Tensor3&,
                       const Tensor3&,
                       const Tensor4&,
                       const Index&,
                       const ArrayOfIndex&,
                       const Vector&,
                       const Vector&,
                       const Index& );

#endif /* ENABLE_DISORT */

#endif /* disort_h */


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

void iy_calc(
              Matrix&         iy,
              Ppath&          ppath,
              Ppath&          ppath_step,
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
        const Vector&         pos,
        const Vector&         los,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const bool&           agenda_verb );

void rte_step_std(
         //Output and Input:
         VectorView stokes_vec,
         //Input
         ConstMatrixView ext_mat_av,
         ConstVectorView abs_vec_av,
         ConstVectorView sca_vec_av, 
         const Numeric& l_step,
         const Numeric& rte_planck_value );

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

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
#include "ppath.h"
#include "matpackI.h"
#include "matpackIII.h"



/*===========================================================================
  === Functions in rte.cc
  ===========================================================================*/

void get_radiative_background(
              Matrix&         i_rte,
              Ppath&          ppath_step,
              Vector&         rte_pos,
              Vector&         rte_los,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Matrix&         i_space,
              Matrix&         ground_emission, 
              Matrix&         ground_los, 
              Tensor4&        ground_refl_coeffs,
              Ppath&          ppath,
        const Agenda&         ppath_step_agenda,
        const Agenda&         rte_agenda,
        const Agenda&         i_space_agenda,
        const Agenda&         ground_refl_agenda,
        const Index&          atmosphere_dim,
        ConstVectorView       p_grid,
        ConstVectorView       lat_grid,
        ConstVectorView       lon_grid,
        const Tensor3&        z_field,
        const Tensor3&        t_field,
        ConstMatrixView       r_geoid,
        ConstMatrixView       z_ground,
        const Index&          cloudbox_on, 
        const ArrayOfIndex&   cloudbox_limits,
        const Tensor7&        scat_i_p,
        const Tensor7&        scat_i_lat,
        const Tensor7&        scat_i_lon,
        ConstVectorView       scat_za_grid,
        ConstVectorView       scat_aa_grid,
        ConstVectorView       f_grid,
        const Index&          stokes_dim,
        const Index&          agenda_verb );

void ground_specular_los(
              VectorView   los,
        const Index&       atmosphere_dim,
        ConstMatrixView    r_geoid,
        ConstMatrixView    z_ground,
        ConstVectorView    lat_grid,
        ConstVectorView    lon_grid,
        const GridPos&     rte_gp_lat,
        const GridPos&     rte_gp_lon,
        ConstVectorView    rte_los );

void rte_step_clearsky_with_emission(
              VectorView    stokes_vec,                
        const Index&        stokes_dim,
        ConstMatrixView     ext_mat_gas,
        ConstVectorView     abs_vec_gas,
        const Numeric&      l_step,
        const Numeric&      rte_planck_value );

void rte_calc(
              Vector&         y,
              Ppath&          ppath,
              Ppath&          ppath_step,
              Matrix&         i_rte,
              Vector&         rte_pos,
              Vector&         rte_los,
              GridPos&        rte_gp_p,
              GridPos&        rte_gp_lat,
              GridPos&        rte_gp_lon,
              Matrix&         i_space,
              Matrix&         ground_emission, 
              Matrix&         ground_los, 
              Tensor4&        ground_refl_coeffs,
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
        const Matrix&         sensor_pol,
        const Vector&         sensor_rot,
        const Vector&         f_grid,
        const Index&          stokes_dim,
        const Index&          antenna_dim,
        const Vector&         mblock_za_grid,
        const Vector&         mblock_aa_grid,
        const bool&           check_input,
        const bool&           apply_sensor,
        const bool&           agenda_verb );

void
rte_step(//Output and Input:
         VectorView stokes_vec,
         //Input
         ConstMatrixView ext_mat_av,
         ConstVectorView abs_vec_av,
         ConstVectorView sca_vec_av, 
         const Numeric& l_step,
         const Numeric& rte_planck_value );


void
stokes_vecGeneral(//WS Output and Input:
                  VectorView stokes_vec,
                  //Input
                  ConstMatrixView ext_mat_av,
                  ConstVectorView abs_vec_av,
                  ConstVectorView sca_vec_av, 
                  const Numeric& l_step,
                  const Numeric& rte_planck_value );




#endif  // rte_h

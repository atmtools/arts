/* Copyright (C) 2023 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
    @file    m_poslos.cc
    @author  Patrick Eriksson <patrick.eriksson@chalmers.se>
    @date    2023-01-14
 
    @brief   Workspace methods for setting and extracting positions (pos)
             and line-of-sights (los).  
*/


/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "auto_md.h"
#include "check_input.h"
#include "geodetic.h"
#include "lin_alg.h"
#include "ppath.h"
#include "variousZZZ.h"


/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losGeometricToPosition(Vector& rte_los,
                                const Vector& refellipsoid,
                                const Vector& rte_pos,
                                const Vector& target_pos,
                                const Verbosity&)
{
    chk_rte_pos("rte_pos", rte_pos);
    chk_rte_pos("target_pos", target_pos);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);
    rte_los.resize(2);
    
    geodetic2ecef(ecef, rte_pos, refellipsoid);
    geodetic2ecef(ecef_target, target_pos, refellipsoid);
    ecef_vector_distance(decef, ecef, ecef_target);
    decef /= norm2(decef);
    ecef2geodetic_los(dummy, rte_los, ecef, decef, refellipsoid);    
}


/* Workspace method: Doxygen documentation will be auto-generated */
void rte_losRefractedToPosition(Workspace& ws,
                                Vector& rte_los,
                                Ppath& ppath,
                                const Agenda& refr_index_air_ZZZ_agenda,
                                const Numeric& ppath_lstep,
                                const Numeric& ppath_lraytrace,
                                const Vector& refellipsoid,
                                const GriddedField2& surface_elevation,
                                const Numeric& surface_search_accuracy,
                                const Vector& rte_pos,
                                const Vector& target_pos,
                                const Numeric& target_dl,
                                const String& algorithm,
                                const Index& max_iterations,
                                const Index& robust,
                                const Numeric& z_toa,
                                const Index& do_horizontal_gradients,
                                const Index& do_twosided_perturb,
                                const Verbosity&)
{
    chk_rte_pos("rte_pos", rte_pos);
    chk_rte_pos("target_pos", target_pos);

    if (algorithm == "basic") {
      refracted_link_basic(ws,
                           ppath,
                           refr_index_air_ZZZ_agenda,
                           ppath_lstep,
                           ppath_lraytrace,
                           refellipsoid,
                           surface_elevation,
                           surface_search_accuracy,
                           z_toa,
                           do_horizontal_gradients,
                           do_twosided_perturb,
                           rte_pos,
                           target_pos,
                           target_dl,
                           max_iterations,
                           robust);
    
    } else {
      ARTS_USER_ERROR("Allowed options for *algorithm* are: \"basic\n");
    }
  
    rte_los = ppath.start_los;
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricToPosition(Matrix& sensor_los,
                                   const Vector& refellipsoid,
                                   const Matrix& sensor_pos,
                                   const Vector& target_pos,
                                   const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_rte_pos("target_pos", target_pos);
    
    const Index n = sensor_pos.nrows();
    sensor_los.resize(n, 2);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);
    geodetic2ecef(ecef_target, target_pos, refellipsoid);

    for (Index i=0; i<n; ++i) {
      geodetic2ecef(ecef, sensor_pos(i, joker), refellipsoid);
      ecef_vector_distance(decef, ecef, ecef_target);
      decef /= norm2(decef);
      ecef2geodetic_los(dummy, sensor_los(i, joker), ecef, decef, refellipsoid);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losGeometricToPositions(Matrix& sensor_los,
                                    const Vector& refellipsoid,
                                    const Matrix& sensor_pos,
                                    const Matrix& target_pos,
                                    const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_sensor_pos("target_pos", target_pos);
    
    const Index n = sensor_pos.nrows();
    ARTS_USER_ERROR_IF(target_pos.nrows() != n,
        "*sensor_pos* and *target_pos* must have the same number of rows.");
    sensor_los.resize(n, 2);

    Vector ecef(3), ecef_target(3), decef(3), dummy(3);

    for (Index i=0; i<n; ++i) {
      geodetic2ecef(ecef, sensor_pos(i, joker), refellipsoid);
      geodetic2ecef(ecef_target, target_pos(i, joker), refellipsoid);
      ecef_vector_distance(decef, ecef, ecef_target);
      decef /= norm2(decef);
      ecef2geodetic_los(dummy, sensor_los(i, joker), ecef, decef, refellipsoid);
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losRefractedToPosition(Workspace& ws,
                                   Matrix& sensor_los,
                                   const Agenda& refr_index_air_ZZZ_agenda,
                                   const Numeric& ppath_lstep,
                                   const Numeric& ppath_lraytrace,
                                   const Vector& refellipsoid,
                                   const GriddedField2& surface_elevation,
                                   const Numeric& surface_search_accuracy,
                                   const Matrix& sensor_pos,
                                   const Vector& target_pos,
                                   const Numeric& target_dl,
                                   const String& algorithm,
                                   const Index& max_iterations,
                                   const Index& robust,
                                   const Numeric& z_toa,
                                   const Index& do_horizontal_gradients,
                                   const Index& do_twosided_perturb,
                                   const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_rte_pos("target_pos", target_pos);
    
    const Index n = sensor_pos.nrows();
    sensor_los.resize(n, 2);

    if (algorithm == "basic") {
      for (Index i=0; i<n; ++i) {
        Ppath ppath;
        refracted_link_basic(ws,
                             ppath,
                             refr_index_air_ZZZ_agenda,
                             ppath_lstep,
                             ppath_lraytrace,
                             refellipsoid,
                             surface_elevation,
                             surface_search_accuracy,
                             z_toa,
                             do_horizontal_gradients,
                             do_twosided_perturb,
                             sensor_pos(i, joker),
                             target_pos,
                             target_dl,
                             max_iterations,
                             robust);
        sensor_los(i, joker) = ppath.start_los;
      }
      
  } else {
    ARTS_USER_ERROR("Allowed options for *algorithm* are: \"basic\n");
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void sensor_losRefractedToPositions(Workspace& ws,
                                    Matrix& sensor_los,
                                    const Agenda& refr_index_air_ZZZ_agenda,
                                    const Numeric& ppath_lstep,
                                    const Numeric& ppath_lraytrace,
                                    const Vector& refellipsoid,
                                    const GriddedField2& surface_elevation,
                                    const Numeric& surface_search_accuracy,
                                    const Matrix& sensor_pos,
                                    const Matrix& target_pos,
                                    const Numeric& target_dl,
                                    const String& algorithm,
                                    const Index& max_iterations,
                                    const Index& robust,
                                    const Numeric& z_toa,
                                    const Index& do_horizontal_gradients,
                                    const Index& do_twosided_perturb,
                                    const Verbosity&)
{
    chk_sensor_pos("sensor_pos", sensor_pos);
    chk_sensor_pos("target_pos", target_pos);
    
    const Index n = sensor_pos.nrows();
    ARTS_USER_ERROR_IF(target_pos.nrows() != n,
        "*sensor_pos* and *target_pos* must have the same number of rows.");
    sensor_los.resize(n, 2);

    if (algorithm == "basic") {
      for (Index i=0; i<n; ++i) {
        Ppath ppath;
        refracted_link_basic(ws,
                             ppath,
                             refr_index_air_ZZZ_agenda,
                             ppath_lstep,
                             ppath_lraytrace,
                             refellipsoid,
                             surface_elevation,
                             surface_search_accuracy,
                             z_toa,
                             do_horizontal_gradients,
                             do_twosided_perturb,
                             sensor_pos(i, joker),
                             target_pos(i, joker),
                             target_dl,
                             max_iterations,
                             robust);
        sensor_los(i, joker) = ppath.start_los;
      }
      
  } else {
    ARTS_USER_ERROR("Allowed options for *algorithm* are: \"basic\n");
  }
}

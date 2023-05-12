/* Copyright (C) 2021 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
    @file   ppath_OLD.h
    @author Patrick Eriksson <patrick.eriksson@chalmers.se>
    @date   2023-01-01

    @brief  Functions releated to calculation of propagation paths.
*/

#ifndef ppath_h
#define ppath_h

#include "agenda_class.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "surf.h"


/*===========================================================================
  === Functions
  ===========================================================================*/

/** Calculates the geometrical length to the surface

    A negative length is returned if the geomtrical path has no intersection
    with the surface.

    The function also check that the observation position is actually above the
    surface. 

    @param[in]  rte_pos                  As the WSV with the same name.
    @param[in]  rte_los                  As the WSV with the same name.
    @param[in]  ecef                     rte_pos in ECEF.
    @param[in]  decef                    rte_los in ECEF.
    @param[in]  refellipsoid             As the WSV with same name.
    @param[in]  surface_elevation        As the WSV with same name.
    @param[in]  surface_search_accuracy  See WSM IntersectionGeometricalWithSurface.
                surface_search_safe      See WSM IntersectionGeometricalWithSurface.

    @return  Length to the surface.

    @author  Patrick Eriksson
    @date    2021-08-06
*/
Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const SurfaceField& surface_field,
                                     const Numeric& surface_search_accuracy,
                                     const Index& surface_search_safe);


/** Adding grid crossings to a ppath

    When called the variable ppath shall contain a description of the
    propagation path fine enough that it can be treated as piecewise
    linear. In addition, the spacing of the input ppath must be fine
    enough that altitudes, latitudes and longitides vary monotonically
    over each ppaths step. The function will not catch deviations from
    these criteria, neither crash, but the result will not be ideal.
   
    The output ppath contains all found grid crossings and other points
    added to meet l_step_max. Set a grid to be empty, to not add points
    for that dimension.

    @param[in,out]  ppath         As the WSV with the same name.
    @param[in]      refellipsoid  As the WSV with same name.
    @param[in]      z_grid        A set of sorted altitudes.
    @param[in]      lat_grid      A set of sorted latitudes.
    @param[in]      lon_grid      A set of sorted longitudes.
    @param[in]      ppath_lstep   As the WSV with same name.

    @author Patrick Eriksson
    @date   2022-09-29
*/
void ppath_add_grid_crossings(Ppath& ppath,
                              const Vector2 refellipsoid,
                              const Vector& z_grid,
                              const Vector& lat_grid,
                              const Vector& lon_grid,
                              const Numeric& l_step_max);


/** Extends a ppath with another one

    The propagation path described by *ppath2* is appended to *ppath*.

    The first point in *ppath2* shall be the same as the last one in
    *ppath*. This common point is repeated in the new *ppath*. The
    corresponding lstep value is (of course) zero.

    @param[in,out]  ppath   A first ppath, to be extended with *ppath2*
    @param[in]      ppath2  A second ppath

    @author Patrick Eriksson
    @date   2023-01-06
*/
void ppath_extend(Ppath& ppath,
                  const Ppath& ppath2);


/** Locates rte_pos with respect to the top of the atmosphere

    @param[out]  l2toa         Length to TOA. Set to zero if sensor is
                               inside  of atmosphere, and to -1 if TOA 
                               not reached from above.
    @param[in]   rte_pos       As the WSV with the same name.
    @param[in]   rte_los       As the WSV with the same name.
    @param[in]   ecef          rte_pos in ECEF.
    @param[in]   decef         rte_los in ECEF.
    @param[in]   refellipsoid  As the WSV with same name.
    @param[in]   z_toa         Top-of-the-atmosphere altitude

    @return  True if rte_pos is above z_toa, false otherwise

    @author Patrick Eriksson
    @date   2023-01-01
*/
bool ppath_l2toa_from_above(Numeric& l2toa,
                            ConstVectorView rte_pos,
                            ConstVectorView rte_los,
                            ConstVectorView ecef,
                            ConstVectorView decef,
                            const Vector2 refellipsoid,
                            const Numeric& z_toa);


/** Basic algorithm for finding refracted path between two points

    This is the implementation of the basic option of
    *ppathRefractedToPosition* and see built-in doc of that WSM for a
    short description of the algorithm and its GIN parameters.

    @param[in,out]  ws
    @param[out]     ppath
    @param[in]      refr_index_air_ZZZ_agenda  As the WSV with same name
    @param[in]      ppath_lstep                As the WSV with same name
    @param[in]      ppath_lraytrace            As the WSV with same name
    @param[in]      surface_elevation          As the WSV with same name
    @param[in]      surface_search_accuracy    As the WSV with same name
    @param[in]      z_toa                      As the GIN with same name
    @param[in]      do_horizontal_gradients    As the GIN with same name
    @param[in]      do_twosided_perturb        As the GIN with same name
    @param[in]      start_pos                  Matches WSV *rte_pos*
    @param[in]      target_pos                 As the GIN with same name
    @param[in]      target_dl                  As the GIN with same name
    @param[in]      max_iterations             As the GIN with same name
    @param[in]      robust                     As the GIN with same name

    @author Patrick Eriksson
    @date   2023-01-08
*/
void refracted_link_basic(Workspace& ws,
                          Ppath& ppath,
                          const Agenda& refr_index_air_ZZZ_agenda,
                          const Numeric& ppath_lstep,
                          const Numeric& ppath_lraytrace,
                          const SurfaceField& surface_field,
                          const Numeric& surface_search_accuracy,
                          const Numeric& z_toa,
                          const Index& do_horizontal_gradients,
                          const Index& do_twosided_perturb,
                          const Vector& start_pos,
                          const Vector& target_pos,
                          const Numeric& target_dl,
                          const Index& max_iterations,
                          const Index& robust);


/** Direction of specular reflection

    The surface is assumed to have no roughness, but topography is
    considered. That is, the local tilt of the surface due to
    variation in elevation is allowed to affect the specular direction.

    @param[out]  los_new            Specular direction as LOS vector
    @param[in]   refellipsoid       As the WSV with same name.
    @param[in]   surface_elevation  As the WSV with same name.
    @param[in]   pos2D              Vector with latitude and longitude
    @param[in]   los                Line-of-sight vector
    @param[in]   ignore_topography  Set to true to ignore possible
                                    surface tilt

    @author Patrick Eriksson
    @date   2023-01-07
*/
void specular_los_calc(VectorView los_new,
                       const SurfaceField& surface_field,
                       ConstVectorView pos2D,
                       ConstVectorView los,
                       const bool& ignore_topography = false);

#endif  // ppath_h

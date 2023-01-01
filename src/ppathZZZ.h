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
 * @file   ppath.h
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2021-07-28
 * 
 * @brief  Propagation path basic functions.
 * 
 * As these functions are called many times, very few checks and asserts are
 * included. Correctness of input shall be done by calling function. 
 */

#ifndef ppathZZZ_h
#define ppathZZZ_h

#include "agenda_class.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "ppath.h"

/*===========================================================================
  === Functions
  ===========================================================================*/

/** Calculates the geometrical length to the cloudbox

   A negative length is returned if the geomtrical path has no intersection
   with the cloudbox.

   @param[in]   rte_pos         As the WSV with the same name.
   @param[in]   rte_los         As the WSV with the same name.
   @param[in]   ecef            *rte_pos* in ECEF.
   @param[in]   decef           *rte_los* in ECEF.
   @param[in]   atmosphere_dim  As the WSV with the same name.
   @param[in]   refellipsoid    As the WSV with same name.
   @param[in]   z_grid          As the WSV with same name.
   @param[in]   lon_grid        As the WSV with same name.
   @param[in]   lat_grid        As the WSV with same name.
   @param[in]   cloudbox_on     As the WSV with same name.
   @param[in]   cloudbox_limits As the WSV with same name.
   @param[in]   is_outside      Shall be set to true/false if *rte_pos* is
                                outside/inside of cloudbox.

   @return   Length to the cloudbox.

   @author Patrick Eriksson
   @date   2021-08-12
 */
Numeric find_crossing_cloudbox(const Vector rte_pos,
                               const Vector rte_los,
                               const Vector ecef,
                               const Vector decef,
                               const Index& atmosphere_dim,
                               const Vector& refellipsoid,
                               const Vector& z_grid,
                               const Vector& lat_grid,
                               const Vector& lon_grid,
                               const Index& cloudbox_on,
                               const ArrayOfIndex& cloudbox_limits,
                               const bool& is_outside);

/** Checks if a position is outside of the atmosphere.

   For the return value, this coding is used:
     0: Inside the atmosphere
     1: Above the atmosphere
     2: Outside in the latitude dimension
     3: Outside in the longitide dimension

   The value 1 just ensures that the position is above the highest altitude in
   z_grid, the position could still be outside in latitude and longitude. 

   Values 2 and 3 should not be OK for sensor positions.

   There is no check against the surface altitude, i.e. 0 could still be a
   point below the surface.

   @param[in]   pos               Position vector.
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   z_toa             Top-of-the-atmosphere altitude.
   @param[in]   lat_grid          As the WSV with the same name.
   @param[in]   lon_grid          As the WSV with the same name.

   @return   See above.

   @author Patrick Eriksson
   @date   2021-07-28
 */
Index is_pos_outside_atmosphere(const Vector pos,
                                const Index& atmosphere_dim,
                                const Numeric& z_toa,
                                const Vector& lat_grid,
                                const Vector& lon_grid);

/** Geometric ppath, with same distance between all points

   Generates a full ppath structure holding a geomtrical path. The points of the
   ppath have equidistant spacing, according to l_step_max. The total length of
   the ppath does not exceed l_total_max (if > 0). 
   
   Intersections with TOA, surface and cloudbox are considered. The ppath is 
   checked to fully be inside the atmosphere.

   @param[out]  ppath             As the WSV with the same name.
   @param[in]   rte_pos           As the WSV with the same name.
   @param[in]   rte_los           As the WSV with the same name.
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with same name.
   @param[in]   z_grid            As the WSV with same name.
   @param[in]   lat_grid          As the WSV with same name.
   @param[in]   lon_grid          As the WSV with same name.
   @param[in]   cloudbox_on       As the WSV with same name.
   @param[in]   cloudbox_limits   As the WSV with same name.
   @param[in]   surface_elevation As the WSV with same name.
   @param[in]   l_step_max        Max distance between points of ppath.
   @param[in]   l_total_max       Max total length of ppath.
   @param[in]   l_accuracy        See WSM IntersectionGeometricalWithSurface.
   @param[in]   safe_search       See WSM IntersectionGeometricalWithSurface.

   @return  Radiative background. Set to be undefined if not found to be space
   or cloudbox.

   @author Patrick Eriksson
   @date   2021-08-11
 */
void ppath_geom_const_lstep(Ppath& ppath,
                            const Vector& rte_pos,
                            const Vector& rte_los,
                            const Index& atmosphere_dim,
                            const Vector& refellipsoid,
                            const Vector& z_grid,
                            const Vector& lat_grid,
                            const Vector& lon_grid,
                            const Index& cloudbox_on,
                            const ArrayOfIndex& cloudbox_limits,
                            const GriddedField2& surface_elevation,
                            const Numeric& l_step_max,
                            const Numeric& l_total_max,
                            const Numeric& l_accuracy,
                            const Index& safe_surface_search,
                            const Index& do_not_calc_gps);

/** Geometric ppath including grid crossings

  When called the variable ppath shall contain a description of the propagation 
  path fine enough that it can be treated as piecewise linear. In addition, the
  spacing of the input ppath must be fine enough that altitudes, latitudes and 
  longitides vary monotonically over each ppaths step. The function will not catch
  deviations from these criteria, neither crash, but the result will not be ideal. 
   
  The output ppath contains all found grid crossings and other points added to meet
  l_step_max.

   @param[in,out]  ppath          As the WSV with the same name.
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with same name.
   @param[in]   z_grid            As the WSV with same name.
   @param[in]   lat_grid          As the WSV with same name.
   @param[in]   lon_grid          As the WSV with same name.
   @param[in]   l_step_max        Max distance between points of ppath.

   @author Patrick Eriksson
   @date   2022-09-29
 */
void ppath_grid_crossings(Ppath& ppath,
                          const Index& atmosphere_dim,
                          const Vector& refellipsoid,
                          const Vector& z_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid,
                          const Numeric& l_step_max,
                          const Numeric& l_accuracy,
                          const Index& do_not_calc_gps);

/** Returns surface elevation at lat and lon of a position

   Length one grids and infinite extrapolation applied. That is, surface grids
   work as the retrieval grids.

   Can also be used for other GriddedField2 that are defined in the same way.

   ZZZ Move to file for surface stuff ZZZ

   @param[in]   pos               Position vector.
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   surface_elevation As the WSV with the same name.

   @return   Elevation.

   @author Patrick Eriksson
   @date   2021-07-28
 */
Numeric surface_z_at_pos(const Vector pos,
                         const Index& atmosphere_dim,
                         const GriddedField2& surface_elevation);

#endif  // ppathZZZ_h

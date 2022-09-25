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
 * @brief  Propagation path structure and basic functions.
 * 
 * This file contains the definition of the Ppath structure and the
 * functions in ppath_basic.cc that are of interest elsewhere.
 *
 * As these functions are called many times, very few checks and asserts are
 * included. Correctness of input shall be done by calling function. 
 */

#ifndef ppathZZZ_h
#define ppathZZZ_h

#include "agenda_class.h"
#include "interpolation.h"
#include "gridded_fields.h"
#include "ppath.h"



/*===========================================================================
  === Functions
  ===========================================================================*/

/** Calculates the geometrical length to the cloudbox

   A negative length is returned if the geomtrical path has no intersection
   with the cloudbox.

   @param[in]   rte_pos         As the WSV with the same name.
   @param[in]   rte_los         As the WSV with the same name.
   @param[in]   ecef            *rte_pos* in ECEF
   @param[in]   decef           *rte_los* in ECEF
   @param[in]   atmosphere_dim  As the WSV with the same name.
   @param[in]   refellipsoid    As the WSV with same name.
   @param[in]   z_grid          As the WSV with same name.
   @param[in]   lon_grid        As the WSV with same name.
   @param[in]   lat_grid        As the WSV with same name.
   @param[in]   is_outside      Shall be set to true/false if *rte_pos* is
                                outside/inside of cloudbox

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

/** Calculates the geometrical length to the surface

   A negative length is returned if the geomtrical path has no intersection
   with the surface.

   The function also check that the observation position is actually above the
   surface. 

   @param[in]   rte_pos           As the WSV with the same name.
   @param[in]   rte_los           As the WSV with the same name.
   @param[in]   ecef              rte_pos in ECEF
   @param[in]   decef             rte_los in ECEF
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with same name.
   @param[in]   surface_elevation As the WSV with same name.
   @param[in]   l_accuracy        See WSM IntersectionGeometricalWithSurface
   @param[in]   safe_search       See WSM IntersectionGeometricalWithSurface

   @return   Length to the surface.

   @author Patrick Eriksson
   @date   2021-08-06
 */
Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const Index& atmosphere_dim,
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& l_accuracy,
                                     const Index& safe_search);

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

/** Adjusts longitudes and calculates grid positions of a ppath
   
   A help function doing two things:

   If 3D, adjusts the longitudes in pos, start_pos and end_pos so they match
   lon_grid

   Calculates grid positions

   @param[in,out]  ppath      ppath-structure with at least fields dim, np,
                              pos, end_pos and start_pos set
   @param[in]   z_grid        As the WSV with the same name.
   @param[in]   lat_grid      As the WSV with the same name.
   @param[in]   lon_grid      As the WSV with the same name.


   @author Patrick Eriksson
   @date   2021-08-16
 */
void ppath_fix_lon_and_gp(Ppath& ppath,
                          const Vector& z_grid,
                          const Vector& lat_grid,
                          const Vector& lon_grid);

/** Some initial steps to determine a propagation path

   The function performs initial checks that are independent if the geometrical
   or refracted path shall be calculated. 

   If sensor is outside of the atmosphere, the function checks that the
   propagation path enters the atmosphere from the top (or is totally in
   space).

   The function also checks if rte_pos is inside the cloudbox. Inside here
   includes to be at the boundary.

   No check with respect to the surface is made.

   @param[out]  l2toa             Length to TOA. Set to be negative if sensor
                                  is inside  of atmosphere
   @param[out]  pos_toa           If l2toa >= 0, holds the position where the
                                  path enters the atmosphere.
   @param[out]  los_toa           If l2toa >= 0, holds the LOS at pos_toa
   @param[in]   rte_pos           As the WSV with the same name.
   @param[in]   rte_los           As the WSV with the same name.
   @param[in]   ecef              rte_pos in ECEF
   @param[in]   decef             rte_los in ECEF
   @param[in]   atmosphere_dim    As the WSV with the same name.
   @param[in]   refellipsoid      As the WSV with same name.
   @param[in]   surface_elevation As the WSV with same name.
   @param[in]   l_accuracy        See WSM IntersectionGeometricalWithSurface
   @param[in]   safe_search       See WSM IntersectionGeometricalWithSurface

   @return  Radiative background. Set to be undefined if not found to be space
   or cloudbox.

   @author Patrick Eriksson
   @date   2021-08-11
 */
enum PpathBackground ppath_init_calc(Numeric& l2toa,
                                     VectorView pos_toa,
                                     VectorView los_toa,
                                     const Vector rte_pos,
                                     const Vector rte_los,
                                     ConstVectorView ecef,
                                     ConstVectorView decef,
                                     const Index& atmosphere_dim,
                                     const Vector& refellipsoid,
                                     const Vector& z_grid,
                                     const Vector& lat_grid,
                                     const Vector& lon_grid,
                                     const Index& cloudbox_on,
                                     const ArrayOfIndex& cloudbox_limits);

/** Returns surface elevation at lat and lon of a position

   Length one grids and infinite extrapolation applied. That is, surface grids
   work as the retrieval grids.

   Can also be used for other GriddedField2 that are defined in the same way.

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

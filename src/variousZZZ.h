/* Copyright (C) 2022 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

#ifndef variousZZZ_h
#define variousZZZ_h

#include "gridded_fields.h"


/** Checks if a vector is a valid line-of-sight (direction vector)

    The function gives an error message if this is not the case.

    @param[in]  name  Name of the variable, if there is an error
    @param[in]  los   The vector

    @author  Patrick Eriksson 
    @date    2022-12-31
 */
void chk_rte_los(const String& name,
                 ConstVectorView los);

/** Checks if a vector is a valid position

    The function gives an error message if this is not the case.

    @param[in]  name  Name of the variable, if there is an error
    @param[in]  pos   The vector

    @author  Patrick Eriksson 
    @date    2022-12-31
 */
void chk_rte_pos(const String& name,
                 ConstVectorView pos);

/** Performs all needed checks of refellipsoid

    @param[in]  refellipsoid  As the WSV with the same name.

    @author  Patrick Eriksson 
    @date    2021-07-30
 */
void chk_refellipsoidZZZ(ConstVectorView refellipsoid);

/** Checks if a matrix is a valid sensor_pos

    @param[in]  name        Name of the variable, if there is an error
    @param[in]  sensor_pos  The matrix

    @author  Patrick Eriksson 
    @date    2021-07-30
 */
void chk_sensor_pos(const String& name,
                    ConstMatrixView sensor_pos);


/** Checks if a matrix is a valid sensor_los

    @param[in]  name        Name of the variable, if there is an error
    @param[in]  sensor_pos  The matrix

    @author  Patrick Eriksson 
    @date    2021-07-30
 */
void chk_sensor_los(const String& name,
                    ConstMatrixView sensor_los);

/** Checks if two matrices are valid combination of sensor_posLlos

    @param[in]  name1       Name of the pos variable, if there is an error
    @param[in]  sensor_pos  A first matrix
    @param[in]  name2       Name of the los variable, if there is an error
    @param[in]  sensor_los  A second matrix

    @author  Patrick Eriksson 
    @date    2021-07-30
 */
void chk_sensor_poslos(const String& name1,
                       ConstMatrixView sensor_pos,
                       const String& name2,
                       ConstMatrixView sensor_los);

/** Checks WSV *surface_elevation*

    It is NOT checked that grids are strictly increasing, as
    relatively costly.

    @param[in]  surface_elevation  As the WSV with the same name.

    @author  Patrick Eriksson 
    @date    2021-08-08
 */
void chk_surface_elevation(const GriddedField2& surface_elevation);

/** Calculates the geometrical length to the surface

   A negative length is returned if the geomtrical path has no intersection
   with the surface.

   The function also check that the observation position is actually above the
   surface. 

   @param[in]  rte_pos       As the WSV with the same name.
   @param[in]  rte_los       As the WSV with the same name.
   @param[in]  ecef          rte_pos in ECEF.
   @param[in]  decef         rte_los in ECEF.
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
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& surface_search_accuracy,
                                     const Index& surface_search_safe);

/** Interpolates a GriddedField2 to a position

    Interpolating a gridded 2D field allowing length 1 grids and
    extrapolation of nearest type.

    The values and units in *pos* shall match the grids in *G*.

    Longitudes need special attention (due to [-180,180] vs [0,360]).
    This is handled automatically, but longitude can only be the
    last (column) dimension.

   @param[in]  G    Gridded field to interpolate
   @param[in]  pos  Position vector, length 2

   @return  Result of interpolation

   @author Patrick Eriksson
   @date   2023-01-03
 */
Numeric interp_gfield2(const GriddedField2& G,
                       const Vector& pos);

/** Interpolates a GriddedField3 to a position

    Interpolating a gridded 3D field allowing length 1 grids and
    extrapolation of nearest type.

    The values and units in *pos* shall match the grids in *G*.

    Longitudes need special attention (due to [-180,180] vs [0,360]).
    This is handled automatically, but longitude can only be the
    last (column) dimension.

   @param[in]  G    Gridded field to interpolate
   @param[in]  pos  Position vector, length 3

   @return  Result of interpolation

   @author Patrick Eriksson
   @date   2023-01-03
 */
Numeric interp_gfield3(const GriddedField3& G,
                       const Vector& pos);

/** Calls refr_index_air_agenda to get n and its gradients

    For more accurate calculations, but slower, consider the two bool
    parameters.

    Default is to only determine the altitude gradients, as in general
    this is the dominating term. The latitude and longitude gradients
    are returned as zero. To calculate all three gradients, set 
    *do_vertical_gradients* to true.

    The gradients are obtained by shifting *pos* with small positive
    values. With *do_twosided_perturb* set to true, there is also a
    perturbation in the negative direction.

   @param[out]  refr_index_air           As the WSV with same name
   @param[out]  refr_index_air_group     As the WSV with same name
   @param[out]  dndz                     Altitude gradient of n
   @param[out]  dndlat                   Latitude gradient of n
   @param[out]  dndlon                   Longitude gradient of n
   @param[in]   ws                       The workspace
   @param[in]   refr_index_air_agenda    As the WSV with same name
   @param[in]   pos                      Position vector
   @param[in]   do_horizontal_gradients  See above
   @param[in]   do_twosided_perturb      See above

   @author Patrick Eriksson
   @date   2023-01-05
 */
void refr_index_and_its_gradients(Numeric& refr_index_air,
                                  Numeric& refr_index_air_group,
                                  Numeric& dndz,
                                  Numeric& dndlat,
                                  Numeric& dndlon,
                                  Workspace& ws,
                                  const Agenda& refr_index_air_agenda,
                                  ConstVectorView pos,
                                  const bool& do_horizontal_gradients,
                                  const bool& do_twosided_perturb);

#endif  // variousZZZ_h

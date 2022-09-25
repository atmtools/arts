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
 * @file   ppath_funcs.h
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2021-07-28
 * 
 * @brief  Functions to calculate propagation paths.
 * 
 *
 * As these functions are called many times, very few checks and asserts are
 * included. Correctness of input shall be done by calling function. 
 */

#ifndef ppath_funcs_h
#define ppath_funcs_h

#include "gridded_fields.h"
#include "ppath_struct.h"

Numeric find_crossing_with_surface_z(const Vector rte_pos,
                                     const Vector rte_los,
                                     const Vector ecef,
                                     const Vector decef,
                                     const Index& atmosphere_dim,
                                     const Vector& refellipsoid,
                                     const GriddedField2& surface_elevation,
                                     const Numeric& l_accuracy,
                                     const Index& safe_search);

enum PpathBackground ppath_init_calc(Numeric& l2toa,
                                     VectorView pos_toa,
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
#endif  // ppath_funcs_h

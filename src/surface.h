/* Copyright (C) 2012
   Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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

/**
   @file   surface.h
   @author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
   @date   2012-02-06 

   This file contains definitions of internal functions associated with 
   the surface.
 */

#ifndef surface_h
#define surface_h

#include "complex.h"
#include "matpackIV.h"
#include "mystring.h"
#include "ppath.h"

/**
    Calculates the incidence angle for a flat surface, based on rte_los and
    specular_los.

    @return                     Incidence angle.
    @param[in]   rte_los        As the WSV with the same name.
    @param[in]   specular_los   As the WSV with the same name.

    @author Patrick Eriksson 
    @date   2012-11-15
 */
Numeric calc_incang(ConstVectorView rte_los, ConstVectorView specular_los);

/**
    Lccates the surface with respect to pressure levels

    The function returns the index of the pressure level at or just below the
    surface.

    @return               The found index.
    @param[in] z_surface  A surface value (probably a value from WSV z_surface)
    @param[in] z_profile  A profile of geometrical altitudes (probably a
                          page of z_field) 

    @author Patrick Eriksson 
    @date   2019-10-26
 */
Index index_of_zsurface(const Numeric& z_surface,
                        ConstVectorView z_profile);

/**
    Weights together downwelling radiation and surface emission.

    *iy* must have correct size when function is called.

    @param[out]  iy                 Radiation matrix, amtching 
                                    the WSV with the same name.
    @param[in]   I                  Downwelling radiation, with dimensions
                                    (surface_los, f_grid, stokes_dim)
    @param[in]   surface_los        As the WSV with the same name.
    @param[in]   surface_rmatrix    As the WSV with the same name.
    @param[in]   surface_emission   As the WSV with the same name.

    @author Patrick Eriksson 
    @date   2005-04-07
 */
void surface_calc(Matrix& iy,
                  ConstTensor3View I,
                  ConstMatrixView surface_los,
                  ConstTensor4View surface_rmatrix,
                  ConstMatrixView surface_emission);

/**
    Sets up the surface reflection matrix and emission vector for
    the case of specular reflection.

    The function handles only one frequency at the time.

    See further the surface chapter in the user guide.

    @param[out]  surface_rmatrix   As the WSV with the same name, but slice
                                   for one direction and one frequency.
    @param[out]  surface_emission  As the WSV with the same name, but slice
                                   for one direction and one frequency.
    @param[in]   Rv                Complex amplitude relection coefficient
                                   vertical polarisation.
    @param[in]   Rh                Complex amplitude relection coefficient
                                   horisontal polarisation.
    @param[in]   f                 Frequency (a scalar).
    @param[in]   stokes_dim        As the WSV with the same name.
    @param[in]   surface_skin_t    As the WSV with the same name.

    @author Patrick Eriksson 
    @date   2004-09-24
 */
void surface_specular_R_and_b(MatrixView surface_rmatrix,
                              VectorView surface_emission,
                              const Complex& Rv,
                              const Complex& Rh,
                              const Numeric& f,
                              const Index& stokes_dim,
                              const Numeric& surface_skin_t);

/**
    Peforms basic checks of *surface_props_data* and *surface_props_names*

    No calculaions, just checks

    @param[in] atmosphere_dim       As the WVS with the same name.
    @param[in] lat_grid             As the WVS with the same name.
    @param[in] lon_grid             As the WVS with the same name.
    @param[in] surface_props_data   As the WVS with the same name.
    @param[in] surface_props_names  As the WVS with the same name.

    @author Patrick Eriksson 
    @date   2018-09-01
 */
void surface_props_check(const Index& atmosphere_dim,
                         const Vector& lat_grid,
                         const Vector& lon_grid,
                         const Tensor3& surface_props_data,
                         const ArrayOfString& surface_props_names);

/**
    Peforms an interpolation of *surface_props_data* 

    The function assumes that the intrpolatin gives a single value. The vector
    *v* must have lenght 1. 

    @param[in] v                    Interpolated value
    @param[in] vname                Name of surface variable to interpolate       
    @param[in] atmosphere_dim       As the WVS with the same name.
    @param[in] gp_lat               Pre-calculatd latitude grid positions
    @param[in] gp_lon               Pre-calculatd longitude grid positions
    @param[in] itw                  Pre-calculatd interpolation weight
    @param[in] surface_props_data   As the WVS with the same name.
    @param[in] surface_props_names  As the WVS with the same name.

    @author Patrick Eriksson 
    @date   2018-09-01
 */
void surface_props_interp(Vector& v,
                          const String& vname,
                          const Index& atmosphere_dim,
                          const ArrayOfGridPos& gp_lat,
                          const ArrayOfGridPos& gp_lon,
                          const Matrix& itw,
                          const Tensor3& surface_props_data,
                          const ArrayOfString& surface_props_names);

/**
    Peforms basic checks of the dsurface variables

    No calculaions, just checks

    @param[in] surface_props_data   As the WVS with the same name.
    @param[in] surface_props_names  As the WVS with the same name.
    @param[in] dsurface_names       As the WVS with the same name.
    @param[in] dsurface_ratrix_dx   As the WVS with the same name.
    @param[in] dsurface_emission_dx As the WVS with the same name.

    @author Patrick Eriksson 
    @date   2018-09-01
 */
void dsurface_check(const ArrayOfString& surface_props_names,
                    const ArrayOfString& dsurface_names,
                    const ArrayOfTensor4 dsurface_rmatrix_dx,
                    const ArrayOfMatrix& dsurface_emission_dx);
#endif  // surface_h

/* Copyright (C) 2003-2012 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
  ===  File description 
  ===========================================================================*/

/*!
  \file   refraction.h
  \author Patrick Eriksson <Patrick.Eriksson@chalmers.se>
  \date   2003-01-17
  
  \brief  Refraction functions.
  
   This file contains the definition of the functions in refraction.cc.
*/

#ifndef refraction_h
#define refraction_h

#include "agenda_class.h"
#include "arts.h"
#include "matpack_data.h"

void complex_n_water_liebe93(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t);

void complex_n_ice_matzler06(Matrix& complex_n,
                             const Vector& f_grid,
                             const Numeric& t);

void get_refr_index_1d(Workspace& ws,
                       Numeric& refr_index,
                       Numeric& refr_index_group,
                       const Agenda& refr_index_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r);

void get_refr_index_2d(Workspace& ws,
                       Numeric& refr_index,
                       Numeric& refr_index_group,
                       const Agenda& refr_index_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat);

void get_refr_index_3d(Workspace& ws,
                       Numeric& refr_index,
                       Numeric& refr_index_group,
                       const Agenda& refr_index_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView lon_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat,
                       const Numeric& lon);

void refr_gradients_1d(Workspace& ws,
                       Numeric& refr_index_air,
                       Numeric& refr_index_air_group,
                       Numeric& dndr,
                       const Agenda& refr_index_air_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r);

void refr_gradients_2d(Workspace& ws,
                       Numeric& refr_index,
                       Numeric& refr_index_group,
                       Numeric& dndr,
                       Numeric& dndlat,
                       const Agenda& refr_index_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat);

void refr_gradients_3d(Workspace& ws,
                       Numeric& refr_index,
                       Numeric& refr_index_group,
                       Numeric& dndr,
                       Numeric& dndlat,
                       Numeric& dndlon,
                       const Agenda& refr_index_agenda,
                       ConstVectorView p_grid,
                       ConstVectorView lat_grid,
                       ConstVectorView lon_grid,
                       ConstVectorView refellipsoid,
                       ConstTensor3View z_field,
                       ConstTensor3View t_field,
                       ConstTensor4View vmr_field,
                       ConstVectorView f_grid,
                       const Numeric& r,
                       const Numeric& lat,
                       const Numeric& lon);

/** Refractive index of water and steam for the optical and near infrared \n
 *
 * From:
 *  Revised formulation for the Refractive
 *  Index of Water and Steam as a Function of
 *  Wavelength, Temperature and Density
 *  \n
 *  Journal of Physical and Chemical Reference Data 27, 761 (1998);
 *  https://doi.org/10.1063/1.556029 27, 761
 *  \n
 *  see also
 *  http://www.iapws.org/release.html
 *
 * @param[out]n refractive index
 * @param[in] only_valid_range flag if true refractive index is calculated only
 *                             within range of validity. If false no check is made,
 *                             so you at your own risk.
 * @param[in] frequency frequency
 * @param[in] temperature temperature
 * @param[in] density density of water or steam
  */
void refractive_index_water_and_steam_VisNIR(Numeric& n,
                                             const Index& only_valid_range,
                                             const Numeric& frequency,
                                             const Numeric& temperature,
                                             const Numeric& density);

#endif  // refraction_h

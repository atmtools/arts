/* Copyright (C) 2003 Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>

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
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2003-01-17
  
  \brief  Refraction functions.
  
   This file contains the definition of the functions in refraction.cc.
*/


#ifndef refraction_h
#define refraction_h


#include "agenda_class.h"
#include "arts.h"
#include "matpackI.h"




void get_refr_index_1d(
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        const Index&      agenda_verb,
        ConstVectorView   p_grid,
        const Numeric&    r_geoid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        const Numeric&    r );

void get_refr_index_2d(
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        const Index&      agenda_verb,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   r_geoid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstTensor3View  vmr_field,
        const Numeric&    r,
        const Numeric&    lat );

void get_refr_index_3d(
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        const Index&      agenda_verb,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstMatrixView   r_geoid,
        ConstTensor3View  z_field,
        ConstTensor3View  t_field,
        ConstTensor4View  vmr_field,
        const Numeric&    r,
        const Numeric&    lat,
        const Numeric&    lon );

void refr_gradients_1d(
              Numeric&    refr_index,
              Numeric&    dndr,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        const Index&      agenda_verb,
        ConstVectorView   p_grid,
        const Numeric&    r_geoid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        const Numeric&    r );

void refr_gradients_2d(
              Numeric&    refr_index,
              Numeric&    dndr,
              Numeric&    dndlat,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        const Index&      agenda_verb,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   r_geoid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstTensor3View  vmr_field,
        const Numeric&    r,
        const Numeric&    lat );

void refr_gradients_3d(
              Numeric&    refr_index,
              Numeric&    dndr,
              Numeric&    dndlat,
              Numeric&    dndlon,
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
        const Agenda&     refr_index_agenda,
        const Index&      agenda_verb,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstVectorView   lon_grid,
        ConstMatrixView   r_geoid,
        ConstTensor3View  z_field,
        ConstTensor3View  t_field,
        ConstTensor4View  vmr_field,
        const Numeric&    r,
        const Numeric&    lat,
        const Numeric&    lon );

void refr_index_thayer_1974(
              Numeric&   refr_index,
        const Numeric&   p,
        const Numeric&   t,
        const Numeric&   h2o_vmr );

void refr_index_ir(
              Numeric&   refr_index,
        const Numeric&   p,
		const Numeric&   t );

#endif  // refraction_h

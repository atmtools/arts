/* Copyright (C) 2002 Patrick Eriksson <patrick@rss.chalmers.se>

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
  \file   refraction.cc
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   2003-01-17
  
  \brief  Functions releated to calculation of refractive index.
*/



/*===========================================================================
  === External declarations
  ===========================================================================*/

#include "interpolation.h"
#include "refraction.h"
#include "special_interp.h"



/*===========================================================================
  === The functions (in alphabetical order)
  ===========================================================================*/


//! get_refr_index_1d
/*! 
   Extracts the refractive index for 1D cases.

   The function interpolates the atmospheric pressure and fields, and
   calls *refr_index_agenda* to determine the refractive index for the
   given point.

   The atmosphere is given by its 1D view. That is, the latitude and
   longitude dimensions are removed from the atmospheric fields. For
   example, the temperature is given as a vector (the vertical profile).

   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   p_grid              As the WSV with the same name.
   \param   z_field             The geometric altitude of each pressure surface
   \param   t_field             The temperature profile.
   \param   vmr_field           The VMR profile for each species.
   \param   z                   The geometric altitude of the position of
                                interest.

   \author Patrick Eriksson
   \date   2003-01-16
*/
void get_refr_index_1d(
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   z_field,
        ConstVectorView   t_field,
        ConstMatrixView   vmr_field,
        const Numeric&    z )
{ 
  // Altitude (equal to pressure) grid position
  ArrayOfGridPos   gp(1);
  gridpos( gp, z_field, Vector(1,z) );

  // Altitude interpolation weights
  Matrix   itw(1,2);
  interpweights( itw, gp );

  // Pressure
  Vector   dummy(1);
  itw2p( dummy, p_grid, gp, itw );
  a_pressure = dummy[0];

  // Temperature
  interp( dummy, itw, t_field, gp );
  a_temperature = dummy[0];

  // VMR
  const Index   ns = vmr_field.nrows();
  //
  a_vmr_list.resize(ns);
  //
  for( Index is=0; is<ns; is++ )
    {
      interp( dummy, itw, vmr_field(is,joker), gp );
      a_vmr_list[is] = dummy[0];
    }

  refr_index_agenda.execute();
}



//! get_refr_index_2d
/*! 
   Extracts the refractive index for 2D cases.

   The function interpolates the atmospheric pressure and fields, and
   calls *refr_index_agenda* to determine the refractive index for the
   given point.

   The atmosphere is given by its 2D view. That is, the longitude
   dimension is removed from the atmospheric fields. For example,
   the temperature is given as a matrix.

   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   p_grid              As the WSV with the same name.
   \param   lat_grid            As the WSV with the same name.
   \param   z_field             The geometric altitude of each pressure surface
                                at each latitude.
   \param   t_field             The temperature 2D field.
   \param   vmr_field           The VMR 2D field for each species.
   \param   z                   The geometric altitude of the position of
                                interest.
   \param   lat                 The latitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-17
*/
void get_refr_index_2d(
              Numeric&    a_pressure,
              Numeric&    a_temperature,
              Vector&     a_vmr_list,
              Numeric&    refr_index,
        const Agenda&     refr_index_agenda,
        ConstVectorView   p_grid,
        ConstVectorView   lat_grid,
        ConstMatrixView   z_field,
        ConstMatrixView   t_field,
        ConstTensor3View  vmr_field,
        const Numeric&    z,
        const Numeric&    lat )
{ 
  // Determine the geometric altitudes at *lat*
  const Index      np = p_grid.nelem();
  Vector           z_grid(np);
  ArrayOfGridPos   gp_lat(1);
  //
  gridpos( gp_lat, lat_grid, lat );
  z_at_lat_2d( z_grid, p_grid, lat_grid, z_field, gp_lat[0] );
  
  // Altitude (equal to pressure) grid position
  ArrayOfGridPos   gp_p(1);
  gridpos( gp_p, z_grid, Vector(1,z) );

  // Altitude interpolation weights
  Matrix   itw(1,2);
  interpweights( itw, gp_p );

  // Pressure
  Vector   dummy(1);
  itw2p( dummy, p_grid, gp_p, itw );
  a_pressure = dummy[0];

  // Temperature
  itw.resize(1,4);
  interpweights( itw, gp_p, gp_lat );
  interp( dummy, itw, t_field, gp_p, gp_lat );
  a_temperature = dummy[0];

  // VMR
  const Index   ns = vmr_field.nrows();
  //
  a_vmr_list.resize(ns);
  //
  for( Index is=0; is<ns; is++ )
    {
      interp( dummy, itw, vmr_field(is,joker,joker), gp_p, gp_lat );
      a_vmr_list[is] = dummy[0];
    }

  refr_index_agenda.execute();
}



//! refr_index_thayer_1974
/*! 
   Calculation of microwave refractive index following Thayer 1974.

   Calculates the microwave refractive index due to gases in the
   Earth's atmosphere. The refractivity of dry air and water vapour is
   summed. All other gases has a negligible contribution.

   The parameterisation of Thayer (Radio Science, 9, 803-807, 1974)
   is used. See also Eq. 3 and 5 of Solheim et al. (JGR, 104, pp. 9664).

   \author Patrick Eriksson
   \date   2002-11-13
*/
void refr_index_thayer_1974(
              Numeric&   refr_index,
        const Numeric&   p,
        const Numeric&   t,
        const Numeric&   h2o_vmr )
{
  const Numeric   e = p * h2o_vmr;

  refr_index = 1 + ( 77.6e-8 * ( p - e ) + ( 72e-8 + 3.754e-3 / t ) * e ) / t;
}

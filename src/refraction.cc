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

   \param   refr_index          Output: As the WSV with the same name.
   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   agenda_verb         This argument is given as input to the agenda
                                above to control the verbosity.
   \param   p_grid              As the WSV with the same name.   
   \param   r_geoid             As the WSV with the same name.
   \param   z_field             The geometric altitude of each pressure surface
   \param   t_field             The temperature profile.
   \param   vmr_field           The VMR profile for each species.
   \param   r                   The radius of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-16
*/
void get_refr_index_1d(
              Numeric&    refr_index,
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
        const Numeric&    r )
{ 
  // Altitude (equal to pressure) grid position
  ArrayOfGridPos   gp(1);
  gridpos( gp, z_field, Vector( 1, r - r_geoid ) );

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

  refr_index_agenda.execute( agenda_verb );
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

   \param   refr_index          Output: As the WSV with the same name.
   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   agenda_verb         This argument is given as input to the agenda
                                above to control the verbosity.
   \param   p_grid              As the WSV with the same name.
   \param   lat_grid            As the WSV with the same name.
   \param   r_geoid             As the WSV with the same name.
   \param   z_field             The geometric altitude of each pressure surface
                                at each latitude.
   \param   t_field             The temperature 2D field.
   \param   vmr_field           The VMR 2D field for each species.
   \param   z                   The radius of the position of interest.
   \param   lat                 The latitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-14
*/
void get_refr_index_2d(
              Numeric&    refr_index,
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
        const Numeric&    lat )
{ 
  // Determine the geometric altitudes at *lat*
  const Index      np = p_grid.nelem();
  Vector           z_grid(np);
  ArrayOfGridPos   gp_lat(1);
  //
  gridpos( gp_lat, lat_grid, lat );
  z_at_lat_2d( z_grid, p_grid, lat_grid, z_field, gp_lat[0] );

  // Determine the geoid radius at *lat*
  Matrix   itw(1,2);
  Vector   dummy(1);
  interpweights( itw, gp_lat );
  interp( dummy, itw, r_geoid, gp_lat );
  const Numeric   rgeoid = dummy[0];

  // Altitude (equal to pressure) grid position
  ArrayOfGridPos   gp_p(1);
  gridpos( gp_p, z_grid, Vector( 1, r - rgeoid ) );

  // Altitude interpolation weights
  interpweights( itw, gp_p );

  // Pressure
  itw2p( dummy, p_grid, gp_p, itw );
  a_pressure = dummy[0];

  // Temperature
  itw.resize(1,4);
  interpweights( itw, gp_p, gp_lat );
  interp( dummy, itw, t_field, gp_p, gp_lat );
  a_temperature = dummy[0];

  // VMR
  const Index   ns = vmr_field.npages();
  //
  a_vmr_list.resize(ns);
  //
  for( Index is=0; is<ns; is++ )
    {
      interp( dummy, itw, vmr_field(is,joker,joker), gp_p, gp_lat );
      a_vmr_list[is] = dummy[0];
    }

  refr_index_agenda.execute( agenda_verb );
}



/*! get_refr_index_3d

   Extracts the refractive index for 3D cases.

   The function interpolates the atmospheric pressure and fields, and
   calls *refr_index_agenda* to determine the refractive index for the
   given point.

   \param   refr_index          Output: As the WSV with the same name.
   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   agenda_verb         This argument is given as input to the agenda
                                above to control the verbosity.
   \param   p_grid              As the WSV with the same name.
   \param   lat_grid            As the WSV with the same name.
   \param   lon_grid            As the WSV with the same name.
   \param   r_geoid             As the WSV with the same name.
   \param   z_field             As the WSV with the same name.
   \param   t_field             As the WSV with the same name.
   \param   vmr_field           As the WSV with the same name.
   \param   r                   The radius of the position of interest.
   \param   lat                 The latitude of the position of interest.
   \param   lon                 The longitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-17
*/
void get_refr_index_3d(
              Numeric&    refr_index,
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
        const Numeric&    lon )
{ 
  // Determine the geometric altitudes at *lat* and *lon*
  const Index      np = p_grid.nelem();
  Vector           z_grid(np);
  ArrayOfGridPos   gp_lat(1), gp_lon(1);
  //
  gridpos( gp_lat, lat_grid, lat );
  gridpos( gp_lon, lon_grid, lon );
  z_at_latlon( z_grid, p_grid, lat_grid, lon_grid, z_field, 
                                                        gp_lat[0], gp_lon[0] );
  
  // Determine the geoid radius at *lat*
  Matrix   itw(1,4);
  Vector   dummy(1);
  interpweights( itw, gp_lat, gp_lon );
  interp( dummy, itw, r_geoid, gp_lat, gp_lon );
  const Numeric   rgeoid = dummy[0];

  // Altitude (equal to pressure) grid position
  ArrayOfGridPos   gp_p(1);
  gridpos( gp_p, z_grid, Vector( 1, r - rgeoid ) );

  // Altitude interpolation weights
  itw.resize(1,2);
  interpweights( itw, gp_p );

  // Pressure
  itw2p( dummy, p_grid, gp_p, itw );
  a_pressure = dummy[0];

  // Temperature
  itw.resize(1,8);
  interpweights( itw, gp_p, gp_lat, gp_lon );
  interp( dummy, itw, t_field, gp_p, gp_lat, gp_lon );
  a_temperature = dummy[0];

  // VMR
  const Index   ns = vmr_field.nbooks();
  //
  a_vmr_list.resize(ns);
  //
  for( Index is=0; is<ns; is++ )
    {
      interp( dummy, itw, vmr_field(is,joker,joker,joker), gp_p, gp_lat, 
                                                                      gp_lon );
      a_vmr_list[is] = dummy[0];
    }

  refr_index_agenda.execute( agenda_verb );
}



//! refr_gradients_2d
/*! 
   Determines the refractive index, and its gradients, for the given position.

   The gradients are calculated in pure numerical way. That is, the
   refractive index is calculated for slightly shifted radius or
   latitude and the difference to the refractive index at the given
   point determines the gradient.

   The atmosphere is given by its 2D view. That is, the longitude
   dimension is removed from the atmospheric fields. For example,
   the temperature is given as a matrix.

   \param   refr_index          Output: As the WSV with the same name.
   \param   dndr                Output: Radial gradient of refractive index.
   \param   dndlat              Output: Latitude gradient of refractive index.
   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   agenda_verb         This argument is given as input to the agenda
                                above to control the verbosity.
   \param   p_grid              As the WSV with the same name.
   \param   lat_grid            As the WSV with the same name.
   \param   r_geoid             As the WSV with the same name.
   \param   z_field             The geometric altitude of each pressure surface
                                at each latitude.
   \param   t_field             The temperature 2D field.
   \param   vmr_field           The VMR 2D field for each species.
   \param   r                   The radius of the position of interest.
   \param   lat                 The latitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-14
*/
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
        const Numeric&    lat )
{ 
   get_refr_index_2d( refr_index, a_pressure,  a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, 
                      lat_grid, r_geoid, z_field, t_field, vmr_field, r, lat );

   const Numeric   n0 = refr_index;

   get_refr_index_2d( refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, 
                    lat_grid, r_geoid, z_field, t_field, vmr_field, r+1, lat );

   dndr = refr_index - n0;

   const Numeric   dlat = 1e-4;

   get_refr_index_2d( refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, 
                 lat_grid, r_geoid, z_field, t_field, vmr_field, r, lat+dlat );

   dndlat = ( refr_index - n0 ) / dlat; 

   refr_index = n0;
}



//! refr_gradients_3d
/*! 
   Determines the refractive index, and its gradients, for the given position.

   The gradients are calculated in pure numerical way. That is, the
   refractive index is calculated for slightly shifted radius,
   latitude or longitude and the difference to the refractive index at
   the given point determines the gradient.

   \param   refr_index          Output: As the WSV with the same name.
   \param   dndr                Output: Radial gradient of refractive index.
   \param   dndlat              Output: Latitude gradient of refractive index.
   \param   dndlat              Output: Longitude gradient of refractive index.
   \param   a_pressure          Output: As the WSV with the same name.
   \param   a_temperature       Output: As the WSV with the same name.
   \param   a_vmr_list          Output: As the WSV with the same name.
   \param   refr_index_agenda   As the WSV with the same name.
   \param   agenda_verb         This argument is given as input to the agenda
                                above to control the verbosity.
   \param   p_grid              As the WSV with the same name.
   \param   lat_grid            As the WSV with the same name.
   \param   lon_grid            As the WSV with the same name.
   \param   r_geoid             As the WSV with the same name.
   \param   z_field             As the WSV with the same name.
   \param   t_field             As the WSV with the same name.
   \param   vmr_field           As the WSV with the same name.
   \param   r                   The radius of the position of interest.
   \param   lat                 The latitude of the position of interest.
   \param   lon                 The longitude of the position of interest.

   \author Patrick Eriksson
   \date   2003-01-17
*/
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
        const Numeric&    lon )
{ 
   get_refr_index_3d( refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, lat_grid, 
                 lon_grid, r_geoid, z_field, t_field, vmr_field, r, lat, lon );

   const Numeric   n0 = refr_index;

   get_refr_index_3d( refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, lat_grid, 
               lon_grid, r_geoid, z_field, t_field, vmr_field, r+1, lat, lon );

   dndr = refr_index - n0;

   const Numeric   dlat = 1e-4;

   get_refr_index_3d( refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, lat_grid, 
             lon_grid, r_geoid, z_field, t_field, vmr_field, r, lat+dlat, lon);

   dndlat = ( refr_index - n0 ) / dlat; 

   const Numeric   dlon = 1e-4;

   get_refr_index_3d( refr_index, a_pressure, a_temperature, a_vmr_list, 
                      refr_index_agenda, agenda_verb, p_grid, lat_grid, 
             lon_grid, r_geoid, z_field, t_field, vmr_field, r, lat, lon+dlon);

   dndlon = ( refr_index - n0 ) / dlon; 

   refr_index = n0;
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

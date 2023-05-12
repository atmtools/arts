/* Copyright (C) 2021
   Jon Petersen <jon.petersen@studium.uni-hamburg.de>
   Manfred Brath  <manfred.brath@uni-hamburg.de>

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

#include "matpack_concepts.h"
#include "physics_funcs.h"
#include "arts.h"
#include "auto_md.h"
#include "geodetic.h"
#include "sun.h"
#include "surf.h"
#include <iostream>


/*!
  \file   m_sun.cc
  \author Jon Petersen  <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@.uni-hamburg.de>
  \date   2021-02-08

  \brief  Workspace functions related to simulation of radiation fluxes.

  These functions are listed in the doxygen documentation as entries of the
  file auto_md.h.
*/


using Constant::pi;

/*===========================================================================
  === The functions
  ===========================================================================*/

/* Workspace method: Doxygen documentation will be auto-generated */
void sunsAddSingleBlackbody(ArrayOfSun &suns,
                         Index &suns_do,
                         // Inputs:
                         const Vector &f_grid,
                         const Index &stokes_dim,
                         const Numeric &radius,
                         const Numeric &distance,
                         const Numeric &temperature,
                         const Numeric &latitude,
                         const Numeric &longitude) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the sun (",distance," m) \n"
                     " is smaller than the radius of the sun (", radius," m )")

  Sun& new_sun = suns.emplace_back();

  // spectrum
  new_sun.spectrum=Matrix(f_grid.nelem(), stokes_dim,0. );

  planck(new_sun.spectrum(joker,0), f_grid, temperature);
  new_sun.spectrum *= pi ; // outgoing flux at the surface of the sun.


  new_sun.description = "Blackbody sun" ;
  new_sun.radius = radius;
  new_sun.distance = distance;
  new_sun.latitude = latitude;
  new_sun.longitude = longitude;

  // set flag
  suns_do = 1;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sunsAddSingleFromGrid(ArrayOfSun &suns,
                         Index &suns_do,
                         // Inputs:
                         const Vector &f_grid,
                         const Index &stokes_dim,
                         const GriddedField2& sun_spectrum_raw,
                         const Numeric &radius,
                         const Numeric &distance,
                         const Numeric &temperature,
                         const Numeric &latitude,
                         const Numeric &longitude,
                         const String &description) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the sun (",distance," m) \n"
                     " is smaller than the radius of the sun (", radius," m )")

  // init sun
  Sun& new_sun = suns.emplace_back();
  new_sun.spectrum = regrid_sun_spectrum(sun_spectrum_raw, f_grid, stokes_dim, temperature); // set spectrum
  new_sun.description = description;
  new_sun.radius = radius;
  new_sun.distance = distance;
  new_sun.latitude = latitude;
  new_sun.longitude = longitude;

  // set flag
  suns_do = 1;

}

/* Workspace method: Doxygen documentation will be auto-generated */
void sunsAddSingleFromGridAtLocation(
                         ArrayOfSun &suns,
                         Index &suns_do,
                         // Inputs:
                         const Vector &f_grid,
                         const Index &stokes_dim,
                         const SurfaceField &surface_field,
                         const GriddedField2 &sun_spectrum_raw,
                         const Numeric &radius,
                         const Numeric &distance,
                         const Numeric &temperature,
                         const Numeric &zenith,
                         const Numeric &azimuth,
                         const String &description,
                         const Numeric &location_latitude,
                         const Numeric &location_longitude,
                         const Numeric &location_altitude) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the sun (",distance," m) \n"
                      "is smaller than the radius of the sun (", radius," m )")
  ARTS_USER_ERROR_IF (location_altitude<0.,
                      "The altitude of the solar spectrum should be positiv,\n"
                      "but is ",location_altitude," m) ")

  // from local position to global position
  Numeric toa_altitude = location_altitude;// + refell2r(refellipsoid, location_latitude);
  ARTS_USER_ERROR("ERROR")

  Numeric sun_altitude, sun_latitude, sun_longitude;
  if (zenith){// < ANGTOL){
    ARTS_USER_ERROR("ERROR")
    sun_altitude = distance + toa_altitude;
    sun_latitude = location_latitude;
    sun_longitude = location_longitude;
  } else if (zenith){// > 180 - ANGTOL) {
    ARTS_USER_ERROR("ERROR")
    sun_altitude = distance - toa_altitude;
    sun_latitude = -location_latitude;
    sun_longitude = location_longitude + 180 - 360.0 * Numeric(round((location_longitude - 0.0) / 360.0));
  } else {
    Numeric x, y, z, dx, dy, dz;
    ARTS_USER_ERROR("ERROR")
   // poslos2cart(x,
     //           y,
       //         z,
         //       dx,
           //     dy,
             //   dz,
               // toa_altitude,
   //             location_latitude,
     //           location_longitude,
       //         zenith,
         //       azimuth);

 //   cart2sph(sun_altitude, 
   //         sun_latitude,
     //       sun_longitude,
       //     x+distance*dx,
         //   y+distance*dy,
           // z+distance*dz,
     //       location_latitude,
       //     location_longitude,
         //   zenith, azimuth);
  }


  // Geometric scaling factor, scales the sun spectral irradiance at the given
  // location to the spectral irradiance of the suns surface.
  Numeric scale_factor = (radius*radius + distance*distance)/
                         (radius*radius);

  // init sun
  Sun& new_sun = suns.emplace_back();

  new_sun.spectrum = regrid_sun_spectrum(sun_spectrum_raw, f_grid, stokes_dim, temperature);
  new_sun.spectrum *= scale_factor; // scale to sun surface

  new_sun.description = description;
  new_sun.radius = radius;
  new_sun.distance = sun_altitude;
  new_sun.latitude = sun_latitude;
  new_sun.longitude = sun_longitude;

  // set flag
  suns_do = 1;

}

void sunsOff(Index &suns_do,
             ArrayOfSun &suns){

  // set flag to False (default)
  suns_do = 0;

  // create empty Array of Matrix for the sun_spectrum
  suns.resize(0);

}

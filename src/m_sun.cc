/*===========================================================================
  ===  File description
  ===========================================================================*/

#include "matpack_concepts.h"
#include "messages.h"
#include "physics_funcs.h"
#include "arts.h"
#include "auto_md.h"
#include "geodetic.h"
#include "sun.h"
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
                         const Numeric &longitude,
                         const Verbosity &) {

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
                         const String &description,
                         const Verbosity &verbosity) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the sun (",distance," m) \n"
                     " is smaller than the radius of the sun (", radius," m )")

  // init sun
  Sun& new_sun = suns.emplace_back();
  new_sun.spectrum = regrid_sun_spectrum(sun_spectrum_raw, f_grid, stokes_dim, temperature, verbosity); // set spectrum
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
                         const Vector &refellipsoid,
                         const GriddedField2 &sun_spectrum_raw,
                         const Numeric &radius,
                         const Numeric &distance,
                         const Numeric &temperature,
                         const Numeric &zenith,
                         const Numeric &azimuth,
                         const String &description,
                         const Numeric &location_latitude,
                         const Numeric &location_longitude,
                         const Numeric &location_altitude,
                         const Verbosity &verbosity) {

  // some sanity checks
  ARTS_USER_ERROR_IF (distance<radius,
                      "The distance to the center of the sun (",distance," m) \n"
                      "is smaller than the radius of the sun (", radius," m )")
  ARTS_USER_ERROR_IF (location_altitude<0.,
                      "The altitude of the solar spectrum should be positiv,\n"
                      "but is ",location_altitude," m) ")

  // from local position to global position
  Numeric toa_altitude = location_altitude + refell2r(refellipsoid, location_latitude);
  
  Numeric sun_altitude, sun_latitude, sun_longitude;
  if (zenith < ANGTOL){
    sun_altitude = distance + toa_altitude;
    sun_latitude = location_latitude;
    sun_longitude = location_longitude;
  } else if (zenith > 180 - ANGTOL) {
    sun_altitude = distance - toa_altitude;
    sun_latitude = -location_latitude;
    sun_longitude = location_longitude + 180 - 360.0 * Numeric(round((location_longitude - 0.0) / 360.0));
  } else {
    Numeric x, y, z, dx, dy, dz;
    poslos2cart(x,
                y,
                z,
                dx,
                dy,
                dz,
                toa_altitude,
                location_latitude,
                location_longitude,
                zenith,
                azimuth);

    cart2sph(sun_altitude, 
            sun_latitude,
            sun_longitude,
            x+distance*dx,
            y+distance*dy,
            z+distance*dz,
            location_latitude,
            location_longitude,
            zenith, azimuth);
  }


  // Geometric scaling factor, scales the sun spectral irradiance at the given
  // location to the spectral irradiance of the suns surface.
  Numeric scale_factor = (radius*radius + distance*distance)/
                         (radius*radius);

  // init sun
  Sun& new_sun = suns.emplace_back();

  new_sun.spectrum = regrid_sun_spectrum(sun_spectrum_raw, f_grid, stokes_dim, temperature, verbosity);
  new_sun.spectrum *= scale_factor; // scale to sun surface

  new_sun.description = description;
  new_sun.radius = radius;
  new_sun.distance = sun_altitude;
  new_sun.latitude = sun_latitude;
  new_sun.longitude = sun_longitude;

  // set flag
  suns_do = 1;

}

void sunsChangeGeometry(ArrayOfSun &suns,
                        // Inputs:
                        const Numeric &radius,
                        const Numeric &distance,
                        const Numeric &latitude,
                        const Numeric &longitude,
                        const Index &sun_index,
                        const Verbosity &)
{
  if (sun_index == -999) return;

  ARTS_USER_ERROR_IF(sun_index+1 > suns.nelem(),
                     "The sun index ", sun_index, " is out of range. \n"
                     "The sun array has only ", suns.nelem(), " elements.")

  //some sanity checks
  ARTS_USER_ERROR_IF(distance<radius,
                     "The distance to the center of the sun (", distance, " m) \n"
                     " is smaller than the radius of the sun (", radius, " m )")

  if (radius > 0) suns[sun_index].radius = radius;
  if (distance > 0) suns[sun_index].distance = distance;

  if (latitude != -999)
  {
    ARTS_USER_ERROR_IF(latitude < -90 || latitude > 90,
                       "The latitude of the sun should be between -90 and 90 degrees,\n"
                       "but is ", latitude, " degrees")

    suns[sun_index].latitude = latitude;
  }

  if (longitude != -999)
  {
    ARTS_USER_ERROR_IF(longitude < -180 || longitude > 360,
                       "The longitude of the sun should be between -180 and 360 degrees,\n"
                       "but is ", longitude, " degrees")

    suns[sun_index].longitude = longitude;
  }
}

void sunsOff(Index &suns_do,
             ArrayOfSun &suns,
             const Verbosity &){

  // set flag to False (default)
  suns_do = 0;

  // create empty Array of Matrix for the sun_spectrum
  suns.resize(0);

}

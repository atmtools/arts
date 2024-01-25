/*===========================================================================
  ===  File description
  ===========================================================================*/

#include <workspace.h>

#include "physics_funcs.h"
#include "sun.h"

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
                            // Inputs:
                            const AscendingGrid &f_grid,
                            const Numeric &radius,
                            const Numeric &distance,
                            const Numeric &temperature,
                            const Numeric &latitude,
                            const Numeric &longitude) {
  // some sanity checks
  ARTS_USER_ERROR_IF(distance < radius,
                     "The distance to the center of the sun (",
                     distance,
                     " m) \n"
                     " is smaller than the radius of the sun (",
                     radius,
                     " m )")

  Sun &new_sun = suns.emplace_back();

  // spectrum
  new_sun.spectrum = Matrix(f_grid.nelem(), 4, 0.);

  planck(new_sun.spectrum(joker, 0), f_grid, temperature);
  new_sun.spectrum *= pi;  // outgoing flux at the surface of the sun.

  new_sun.description = "Blackbody sun";
  new_sun.radius = radius;
  new_sun.distance = distance;
  new_sun.latitude = latitude;
  new_sun.longitude = longitude;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void sunsAddSingleFromGrid(ArrayOfSun &suns,
                           // Inputs:
                           const AscendingGrid &f_grid,
                           const GriddedField2 &sun_spectrum_raw,
                           const Numeric &radius,
                           const Numeric &distance,
                           const Numeric &temperature,
                           const Numeric &latitude,
                           const Numeric &longitude,
                           const String &description) {
  // some sanity checks
  ARTS_USER_ERROR_IF(distance < radius,
                     "The distance to the center of the sun (",
                     distance,
                     " m) \n"
                     " is smaller than the radius of the sun (",
                     radius,
                     " m )")

  // init sun
  Sun &new_sun = suns.emplace_back();
  new_sun.spectrum = regrid_sun_spectrum(
      sun_spectrum_raw, f_grid, temperature);  // set spectrum
  new_sun.description = description;
  new_sun.radius = radius;
  new_sun.distance = distance;
  new_sun.latitude = latitude;
  new_sun.longitude = longitude;
}

void sunsOff(ArrayOfSun &suns) {
  // create empty Array of Matrix for the sun_spectrum
  suns.resize(0);
}

/*!
  \file   sun.h
  \author Jon Petersen <jon.petersen@studium.uni-hamburg.de>
          Manfred Brath  <manfred.brath@uni-hamburg.de>
  \date   2021-02-22

  \brief  Declaration of functions in star.cc.
*/

#ifndef star_h
#define star_h

/*===========================================================================
  === External declarations
  ===========================================================================*/

#include <matpack.h>
#include <rtepack.h>

#include "atm.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "optproperties.h"
#include "ppath_struct.h"
#include "surf.h"

/*===========================================================================
  === structs/classes  in sun.h
  ===========================================================================*/

/** The structure to describe a propagation path and releated quantities.
 *
 *  The fields of the structure are described more in detail inside the ARTS
 *  user guide (AUG).
 */
struct Sun {
  /** Sun description */
  String description;
  /** Sun spectrum, monochrmatic radiance spectrum at the surface of the sun*/
  Matrix spectrum;
  /** Sun radius */
  Numeric radius;
  /** distance from center of planet to center of sun*/
  Numeric distance;
  /** latitude of the sun in the sky of the planet */
  Numeric latitude;
  /** longitude of the sun in the sky of the planet */
  Numeric longitude;

  friend std::ostream& operator<<(std::ostream& os, const Sun& sun);
};



/** An array of sun. */
using ArrayOfSun = Array<Sun>;

std::ostream& operator<<(std::ostream& os, const ArrayOfSun& a);
 
 /** regrid_sun_spectrum
 *
 * Regrids a given spectrum from a griddedfield2 to the f_grid.
 * if the f_grid covers a larger range as the given one, one
 * can choose between two padding options:
 * zeros: Intensities outside the given spectrum are set to zero 
 * planck: Intensities outside the given spectrum are initilizied 
 *        with the black body value at that frequency.
 *
 * @param[in]  sun_spectrum_raw  gf2 of the given spectrum.
 * @param[in]  f_grid  f_grid for the calculation.
 * @param[in]  temperature  Temperature for the planck padding.
 *
 * @return     interpolated spectrum
 *
 * @author Jon Petersen
 * @date   2022-01-19
 */
Matrix regrid_sun_spectrum(const GriddedField2& sun_spectrum_raw,
                          const Vector &f_grid,
                          const Numeric &temperature);

#endif /* star_h */

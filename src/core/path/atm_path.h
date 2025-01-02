#pragma once

#include <atm.h>

#include "path_point.h"

/** Helper function that resizes the input ArrayOfAtmPoint.
 *
 * @param[out] ppvar_atm As WSV
 * @param[in] rad_path As WSV
 * @return ArrayOfAtmPoint& As ppvar_atm WSV
 */
ArrayOfAtmPoint &atm_path_resize(ArrayOfAtmPoint &ppvar_atm,
                                 const ArrayOfPropagationPathPoint &rad_path);

/** Fills the propagation path atmospheric point variable, size is known
 *
 * @param[out] ppvar_atm As WSV
 * @param[in] rad_path As WSV
 * @param[in] atm_field As WSV
 */
void forward_atm_path(ArrayOfAtmPoint &ppvar_atm,
                      const ArrayOfPropagationPathPoint &rad_path,
                      const AtmField &atm_field);

/** Outputs the propagation path atmospheric point variable
 *
 * @param[in] rad_path As WSV
 * @param[in] atm_field As WSV
 * @return ArrayOfAtmPoint As ppvar_atm WSV
 */
ArrayOfAtmPoint forward_atm_path(const ArrayOfPropagationPathPoint &rad_path,
                                 const AtmField &atm_field);

/** Frequency shift at a single ray path ppoint
 * 
 * @param path_freq Shifted frequency grid
 * @param main_freq Original frequency grid
 * @param rad_path path point of the radiation
 * @param atm_path path point of the atmosphere
 */
void forward_path_freq(AscendingGrid &path_freq,
                       const AscendingGrid &main_freq,
                       const PropagationPathPoint &rad_path,
                       const AtmPoint &atm_path);

/** Set frequency grid along the atmospheric path
 * 
 * @param[out] ppvar_f As WSV
 * @param[in] f_grid As WSV
 * @param[in] rad_path As WSV
 * @param[in] ppvar_atm As WSV
 * @return ArrayOfVector& As ppvar_f WSV
 */
void forward_path_freq(ArrayOfAscendingGrid &ppvar_f,
                       const AscendingGrid &f_grid,
                       const ArrayOfPropagationPathPoint &rad_path,
                       const ArrayOfAtmPoint &ppvar_atm);

/** Extracts a 1D atmospheric "path" from a 3D atmospheric field
 *
 * @param[in] atm_field As WSV
 * @param[in] z_grid An altitude grid (or 1-long vector for single altitude)
 * @param[in] lat_grid A latitude grid (or 1-long vector for single latitude)
 * @param[in] lon_grid A longitude grid (or 1-long vector for single longitude)
 * @return ArrayOfAtmPoint A "path" through the 3D atmospheric field
 */
ArrayOfAtmPoint extract1D(const AtmField &atm_field,
                          const Vector &z_grid,
                          const Vector &lat_grid = {0},
                          const Vector &lon_grid = {0});

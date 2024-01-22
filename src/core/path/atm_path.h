#pragma once

#include "atm.h"
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

/** Set the size of the output
 * 
 * @param[out] ppvar_f As WSV
 * @param[in] f_grid As WSV
 * @param[in] ppvar_atm As WSV
 * @return ArrayOfVector& As ppvar_f WSV
 */
ArrayOfVector &path_freq_resize(ArrayOfVector &ppvar_f,
                                const Vector &f_grid,
                                const ArrayOfAtmPoint &ppvar_atm);

/** Set frequency grid along the atmospheric path
 * 
 * @param[out] ppvar_f As WSV
 * @param[in] f_grid As WSV
 * @param[in] rad_path As WSV
 * @param[in] ppvar_atm As WSV
 * @param[in] rte_alonglos_v As WSV
 * @return ArrayOfVector& As ppvar_f WSV
 */
void forward_path_freq(ArrayOfVector &ppvar_f,
                       const Vector &f_grid,
                       const ArrayOfPropagationPathPoint &rad_path,
                       const ArrayOfAtmPoint &ppvar_atm,
                       const Numeric rte_alonglos_v);

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

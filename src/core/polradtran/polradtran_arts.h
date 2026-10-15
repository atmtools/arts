#pragma once

#include <atm.h>
#include <matpack.h>
#include <path_point.h>
#include <rtepack.h>

namespace polradtran {
//! The layers of an ARTS propagation path, as layers_from_path makes them
struct path_layers {
  //! freq_grid[freq_index] [Hz]
  Numeric frequency;
  //! [nlev] the altitudes of ray_path [m], top first, strictly decreasing
  Vector height;
  //! [nlev] the temperatures of atm_path [K]
  Vector temperature;
  //! [nlev - 1] per layer, the mean of the A elements of the gas
  //! propagation matrices at its two levels [m-1]
  Vector gas_extinction;
};

/** The layers of a propagation path as RT3 and RT4 take them, for
 * rt3::problem_from_path and rt4::problem_from_path.
 *
 * The conventions are those of the DISORT workspace methods
 * (disort_settingsOpticalThicknessFromPath and
 * disort_settingsLayerThermalEmissionLinearInTau): ray_path, atm_path and
 * spectral_propmat_path have one entry per level, top first, as for a
 * down-looking path, with at least 2 levels.  The altitudes of ray_path must
 * decrease strictly; they are the heights [m], so layer l lies between
 * levels l and l + 1.  Only the altitudes are used, not the lines of sight
 * or the horizontal positions.  spectral_propmat_path is the gas
 * propagation matrix only (no particles), per metre, with freq_grid.size()
 * entries per level; polarized gas propagation matrices are rejected, as
 * the gas extinction of RT3 and RT4 is scalar.
 */
path_layers layers_from_path(const ArrayOfPropagationPathPoint& ray_path,
                             const ArrayOfAtmPoint&             atm_path,
                             const ArrayOfPropmatVector&        spectral_propmat_path,
                             const AscendingGrid&               freq_grid,
                             Index                              freq_index);
}  // namespace polradtran

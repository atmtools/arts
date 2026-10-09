#include "polradtran_arts.h"

#include <debug.h>

#include <algorithm>
#include <numeric>

namespace polradtran {
path_layers layers_from_path(const ArrayOfPropagationPathPoint& ray_path,
                             const ArrayOfAtmPoint&             atm_path,
                             const ArrayOfPropmatVector&        spectral_propmat_path,
                             const AscendingGrid&               freq_grid,
                             Index                              freq_index) {
  const Index nlev = static_cast<Index>(ray_path.size());
  const Index nf   = static_cast<Index>(freq_grid.size());

  ARTS_USER_ERROR_IF(nlev < 2, "ray_path needs at least 2 points (1 layer), got {}", nlev);
  ARTS_USER_ERROR_IF(
      static_cast<Index>(atm_path.size()) != nlev or static_cast<Index>(spectral_propmat_path.size()) != nlev,
      "ray_path, atm_path and spectral_propmat_path must have one entry per level; they have {}, {} "
      "and {}",
      nlev,
      atm_path.size(),
      spectral_propmat_path.size());
  ARTS_USER_ERROR_IF(freq_index < 0 or freq_index >= nf,
                     "freq_index must be in [0, {}) for a freq_grid of {} frequencies, got {}",
                     nf,
                     nf,
                     freq_index);
  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_propmat_path, [nf](const PropmatVector& v) { return static_cast<Index>(v.size()) != nf; }),
      "Every spectral_propmat_path level must have freq_grid.size() = {} propagation matrices",
      nf);
  for (Index l = 0; l < nlev - 1; l++)
    ARTS_USER_ERROR_IF(not(ray_path[l].altitude() > ray_path[l + 1].altitude()),
                       "The ray_path altitudes must decrease strictly from the first point (top of the atmosphere) "
                       "to the last (surface)");
  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_propmat_path,
                   [freq_index](const PropmatVector& v) { return v[freq_index].is_polarized(); }),
      "The gas extinction of RT3 and RT4 is scalar: the gas propagation matrices in spectral_propmat_path must not "
      "be polarized (only A may be non-zero) at frequency index {}",
      freq_index);

  const Index nlay = nlev - 1;
  path_layers layers{.frequency      = freq_grid[freq_index],
                     .height         = Vector(nlev),
                     .temperature    = Vector(nlev),
                     .gas_extinction = Vector(nlay)};
  for (Index i = 0; i < nlev; i++) {
    layers.height[i]      = ray_path[i].altitude();
    layers.temperature[i] = atm_path[i].temperature;
  }
  for (Index i = 0; i < nlay; i++)
    layers.gas_extinction[i] =
        std::midpoint(spectral_propmat_path[i][freq_index].A(), spectral_propmat_path[i + 1][freq_index].A());
  return layers;
}
}  // namespace polradtran

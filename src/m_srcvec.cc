#include <workspace.h>

#include <ranges>

#include "matpack_mdspan_helpers_grid_t.h"
#include "rtepack.h"

void spectral_rad_srcvec_pathFromPropmat(SourceVector&                 spectral_rad_srcvec_path,
                                         const ArrayOfSpectralPropmat& spectral_propmat_path,
                                         const ArrayOfAscendingGrid&   freq_grid_path,
                                         const ArrayOfAtmPoint&        atm_path,
                                         const JacobianTargets&        jac_targets) try {
  ARTS_TIME_REPORT

  const Index it = jac_targets.target_position(AtmKey::t);

  const Vector ts{std::from_range, atm_path | stdv::transform(&AtmPoint::temperature)};

  spectral_rad_srcvec_path.init(spectral_propmat_path, freq_grid_path, ts, it);
}
ARTS_METHOD_ERROR_CATCH

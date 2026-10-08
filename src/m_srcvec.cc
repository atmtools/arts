#include <workspace.h>

#include <ranges>

#include <array_algo.h>
#include <arts_omp.h>

#include "matpack_mdspan_helpers_check.h"
#include "matpack_mdspan_helpers_grid_t.h"
#include "rtepack.h"

void spectral_rad_srcvec_pathFromPropmat(SourceVector&               spectral_rad_srcvec_path,
                                         const ArrayOfPropmatVector& spectral_propmat_path,
                                         const ArrayOfStokvecVector& spectral_nlte_srcvec_path,
                                         const ArrayOfPropmatMatrix& spectral_propmat_jac_path,
                                         const ArrayOfStokvecMatrix& spectral_nlte_srcvec_jac_path,
                                         const ArrayOfAscendingGrid& freq_grid_path,
                                         const ArrayOfAtmPoint&      atm_path,
                                         const JacobianTargets&      jac_targets) try {
  ARTS_TIME_REPORT

  const Index it = jac_targets.target_position(AtmKey::t);

  const Vector ts{std::from_range, atm_path | stdv::transform(&AtmPoint::temperature)};

  spectral_rad_srcvec_path.init(spectral_propmat_path,
                                spectral_propmat_jac_path,
                                spectral_nlte_srcvec_path,
                                spectral_nlte_srcvec_jac_path,
                                freq_grid_path,
                                ts,
                                it);
}
ARTS_METHOD_ERROR_CATCH

void spectral_rad_srcvec_pathZero(SourceVector&               spectral_rad_srcvec_path,
                                  const ArrayOfPropmatVector& spectral_propmat_path,
                                  const ArrayOfPropmatMatrix& spectral_propmat_jac_path,
                                  const ArrayOfAscendingGrid& freq_grid_path) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(not arr::same_size(spectral_propmat_path, spectral_propmat_jac_path, freq_grid_path),
                     "All input must have the same size.");

  const Size np = spectral_propmat_jac_path.size();

  if (np == 0) {
    spectral_rad_srcvec_path.J.resize(spectral_rad_srcvec_path.J.nrows(), 0);
    spectral_rad_srcvec_path.dJ.resize(spectral_rad_srcvec_path.dJ.npages(), 0, spectral_rad_srcvec_path.dJ.ncols());
    return;
  }

  const Size nq = spectral_propmat_jac_path.front().nrows();
  const Size nf = spectral_propmat_jac_path.front().ncols();

  spectral_rad_srcvec_path.J.resize(nf, np);
  spectral_rad_srcvec_path.dJ.resize(nf, np, nq);

  ARTS_USER_ERROR_IF(not all_same_shape({nq, nf}, spectral_propmat_jac_path),
                     "All derivative parameters must have same shape ({}, {}).",
                     nq,
                     nf);

  ARTS_USER_ERROR_IF(not all_same_shape({nf}, spectral_propmat_path, freq_grid_path),
                     "All forward parameters must have same shape ({}).",
                     nf);

#pragma omp parallel for collapse(2) if (!arts_omp_in_parallel())
  for (Size i = 0; i < np; i++) {
    for (Size j = 0; j < nf; j++) {
      spectral_rad_srcvec_path.J[j, i]  = Stokvec{0.0};
      spectral_rad_srcvec_path.dJ[j, i] = Stokvec{0.0};
    }
  }
}
ARTS_METHOD_ERROR_CATCH

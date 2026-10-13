#include "rtepack_mueller_matrix.h"

#include <configtypes.h>

/* Here rather than in rtepack.h, so that the note appears once per build and not
   for every file that includes rtepack. */
#if __cplusplus > 202302L
#warning \
    "rtepack: C++26 makes <cmath> and <complex> constexpr (P1383).  These rtepack functions can then be constexpr in their headers: tran (its constructors and evaluations), exp, sqrt and logK of a propmat, level_emission (with planck), rayleigh_scattering, power_reflectance(rv, rh), specular_reflected_direction, fresnel_reflectance_specular, fresnel_reflectance_nonspecular, specular_radiance, nonspecular_radiance, meridional_basis, scattering_plane_rotation, lab_frame_phase_matrix and jones(const muelmat &)."
#endif

namespace rtepack {
void forward_cumulative_transmission(Array<muelmat_vector> &Pi, const Array<muelmat_vector> &T) {
  const Size N = T.size();

  Pi.resize(N);

  if (N == 0) return;

  const Index nv = T.front().size();
  for (auto &p : Pi) p.resize(nv);

  Pi.front() = T.front();

  for (Size i = 1; i < N; i++) {
    for (auto &&[Pi1, Pi0, T] : stdv::zip(Pi[i], Pi[i - 1], T[i])) { Pi1 = Pi0 * T; }
  }
}

Array<muelmat_vector> forward_cumulative_transmission(const Array<muelmat_vector> &T) {
  Array<muelmat_vector> out;
  forward_cumulative_transmission(out, T);
  return out;
}
}  // namespace rtepack

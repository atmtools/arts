#include "rtepack_polarization.h"

#include <cmath>

namespace rtepack {
polarization_basis meridional_basis(const Vector3& k, const Vector3& h_pole) {
  Vector3       h      = cross(k, Vector3{0.0, 0.0, 1.0});
  const Numeric norm_h = hypot(h);
  if (norm_h < 1e-12) {
    h = h_pole;
  } else {
    h /= norm_h;
  }
  return {cross(h, k), h};
}

Vector4 scattering_plane_rotation(const polarization_basis& in, const polarization_basis& out) {
  const Vector3 k_in = cross(in.first, in.second), k_out = cross(out.first, out.second);
  Vector3       N      = cross(k_in, k_out);
  const Numeric norm_N = hypot(N);
  if (norm_N < 1e-12) {
    N = in.second;
  } else {
    N /= norm_N;
  }

  const auto [c_in, s_in]   = detail::rotation_angle(in, {cross(N, k_in), N});
  const auto [c_out, s_out] = detail::rotation_angle({cross(N, k_out), N}, out);
  return {c_in, s_in, c_out, s_out};
}

muelmat lab_frame_phase_matrix(const compact_planar_muelmat& F,
                               const polarization_basis&     in,
                               const polarization_basis&     out) {
  const auto [c_in, s_in, c_out, s_out] = scattering_plane_rotation(in, out);
  return rotated(F, c_in, s_in, c_out, s_out);
}
}  // namespace rtepack

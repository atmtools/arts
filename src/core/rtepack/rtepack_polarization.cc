#include "rtepack_polarization.h"

#include <cmath>

namespace rtepack {
polarization_basis meridional_basis(const Vector3& k, const Vector3& h_pole) {
  // k x z = (k_y, -k_x, 0) gives the meridian to rounding, except along +-z, where it is 0 / 0 and, within the
  // rounding of k, resolves none
  Vector3       h      = cross(k, Vector3{0.0, 0.0, 1.0});
  const Numeric norm_h = hypot(h);
  if (not(norm_h > detail::parallel)) {
    h = h_pole;
  } else {
    h /= norm_h;
  }
  return {cross(h, k), h};
}

Vector4 scattering_plane_rotation(const polarization_basis& in, const polarization_basis& out) {
  // The rounding of the cross product, about epsilon / |N| in direction, also tilts N out of the plane across k_in,
  // which would shorten the basis vectors of the scattering plane.  Projected back across k_in, the rounding only
  // turns the plane, which reaches the phase matrix through the part of F that depends on the plane alone, vanishing
  // towards forward and backward; so N gives the plane until the rays are parallel in double precision.
  const Vector3 k_in = cross(in.first, in.second), k_out = cross(out.first, out.second);
  Vector3       N       = cross(k_in, k_out);
  N                    -= dot(N, k_in) * k_in;
  const Numeric norm_N  = hypot(N);
  if (not(norm_N > detail::parallel)) {
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

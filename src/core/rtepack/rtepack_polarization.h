#pragma once

#include <debug.h>
#include <matpack.h>

#include <limits>
#include <utility>

#include "rtepack_compact_planar_mueller_matrix.h"
#include "rtepack_mueller_matrix.h"

namespace rtepack {
/** The polarization basis {v, h} of a ray
 *
 *  Orthonormal, perpendicular to the ray's propagation direction k and
 *  right-handed with it, v x h = k.  The Stokes vector of a field
 *  E = E_v v + E_h h in it is I = |E_v|^2 + |E_h|^2, Q = |E_v|^2 - |E_h|^2,
 *  U = 2 Re(E_v E_h*) and V = 2 Im(E_v E_h*).
 */
using polarization_basis = std::pair<Vector3, Vector3>;

namespace detail {
/** The largest cross product of two parallel directions
 *
 *  Unit vectors accurate to rounding, and the bases computed from them, give
 *  cross products accurate to a few epsilon, which this bounds.  At or below
 *  it a cross product resolves no direction, and the two are parallel in
 *  double precision.
 */
constexpr Numeric parallel = 16 * std::numeric_limits<Numeric>::epsilon();

//! {cos 2 psi, sin 2 psi} of the rotation from one basis of a ray to another, cos psi = v_to . v_from and sin psi = v_to . h_from
constexpr std::pair<Numeric, Numeric> rotation_angle(const polarization_basis& from, const polarization_basis& to) {
  const Numeric c = dot(to.first, from.first), s = dot(to.first, from.second);
  return {c * c - s * s, 2.0 * s * c};
}
}  // namespace detail

/** The meridional basis of a ray propagating along the unit vector k
 *
 *  With z up, h = k x z / |k x z| and v = h x k.  v lies in the plane of k
 *  and z with v . z >= 0, the "projected up" direction of the Zeeman
 *  angles, so this is the basis of ARTS's propagation matrices and of
 *  VDISORT, RT3 and RT4.  Along +-z (|k x z| at most detail::parallel) the
 *  plane is that of the horizontal unit vector h_pole: h = h_pole, the limit
 *  along the meridian with that h.  Any other ray has the h of its own
 *  meridian, to rounding.
 *
 *  @param k      Unit propagation direction
 *  @param h_pole Horizontal unit vector, the h of a ray along +-z
 *  @return       {v, h}
 */
polarization_basis meridional_basis(const Vector3& k, const Vector3& h_pole);

/** The Stokes rotation from one polarization basis of a ray to another
 *
 *  stokes_rotation(cos 2 psi, sin 2 psi) with cos psi = v_to . v_from and
 *  sin psi = v_to . h_from.  Throws unless both bases belong to the same
 *  propagation direction, v_from x h_from = v_to x h_to.
 *
 *  @param from The basis of the input Stokes vector
 *  @param to   The basis of the output Stokes vector
 *  @return     The Mueller matrix taking the Stokes vector from one to the other
 */
constexpr muelmat basis_rotation(const polarization_basis& from, const polarization_basis& to) {
  const Vector3 k_from = cross(from.first, from.second), k_to = cross(to.first, to.second);
  const Vector3 dk = k_from - k_to;
  ARTS_USER_ERROR_IF(dot(dk, dk) > 1e-18,
                     "The polarization bases belong to different propagation directions, {} and {}",
                     k_from,
                     k_to);
  const auto [c2, s2] = detail::rotation_angle(from, to);
  return stokes_rotation(c2, s2);
}

/** The rotations of a scattering-plane matrix to the bases of two rays
 *
 *  The scattering-plane basis of a ray along k is (N x k, N) with
 *  N = k_in x k_out / |k_in x k_out|, and a scattering-plane matrix F has
 *  Q = I_par - I_perp in it.  Its phase matrix between the polarization
 *  bases in and out is
 *
 *    rotated(F, cos2psi_in, sin2psi_in, cos2psi_out, sin2psi_out)
 *
 *  i.e. basis_rotation from in to the scattering plane, F, and
 *  basis_rotation from the scattering plane to out.  Exactly forward or
 *  backward (|k_in x k_out| at most detail::parallel) the scattering plane
 *  is undefined, and that of k_in and h_in is used.  A regular F there
 *  (F12 = F34 = 0 and F22 = +-F33) gives the same phase matrix for every
 *  plane.  Near them the cross product gives the plane to about
 *  epsilon / |k_in x k_out|, which reaches the phase matrix only through the
 *  part of F that depends on the plane, vanishing there for a physical F.
 *
 *  @param in  The polarization basis of the incident ray
 *  @param out The polarization basis of the scattered ray
 *  @return    [cos2psi_in, sin2psi_in, cos2psi_out, sin2psi_out]
 */
Vector4 scattering_plane_rotation(const polarization_basis& in, const polarization_basis& out);

/** The phase matrix between two rays of a scattering-plane matrix
 *
 *  F at the scattering angle between in and out, rotated as by
 *  scattering_plane_rotation.
 *
 *  @param F   The scattering-plane matrix at cos(Theta) = k_in . k_out
 *  @param in  The polarization basis of the incident ray
 *  @param out The polarization basis of the scattered ray
 *  @return    The phase matrix from the Stokes vector in in to that in out
 */
muelmat lab_frame_phase_matrix(const compact_planar_muelmat& F,
                               const polarization_basis&     in,
                               const polarization_basis&     out);
}  // namespace rtepack

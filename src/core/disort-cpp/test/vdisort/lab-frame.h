#pragma once

/* The laboratory-frame phase matrix of randomly oriented particles with a
   plane of symmetry, from vector geometry.  Shared by the VDISORT
   comparisons with RT4 (vdisort-rt4-comparison.cpp) and RT3
   (vdisort-rt3-comparison.cpp).  Neither solver's rotation formulas are
   used.

   Frame: z points up.  A ray propagating at zenith cosine mu (> 0 upward)
   and azimuth phi has the direction n = (sin cos phi, sin sin phi, mu).
   Its meridional basis is e_v = (mu cos phi, mu sin phi, -sin) in the plane
   of n and z and e_h = (-sin phi, cos phi, 0), with e_v x e_h = n.  Stokes:
   I = I_v + I_h, Q = I_v - I_h, U = 2 Re(E_v E_h*).  The basis
   h = n x z / |n x z|, v = h x n of the RT3 wrapper and of VDISORT's
   polarized tests is (-e_v, -e_h), which gives the same Stokes vector. */

#include <matpack.h>
#include <rtepack_mueller_matrix.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

namespace vdisort_test {
//! The scattering matrix in the scattering-plane basis (e_par, e_perp), Q = I_par - I_perp
struct tro_elements {
  Numeric F11{}, F12{}, F22{}, F33{}, F34{}, F44{};
};

//! F as a function of cos(Theta)
using tro_matrix = std::function<tro_elements(Numeric)>;

using vec3 = std::array<Numeric, 3>;

inline Numeric dot(const vec3& a, const vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

inline vec3 cross(const vec3& a, const vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

inline vec3 direction(Numeric cos_theta, Numeric phi) {
  const Numeric st = std::sqrt(std::max(0.0, 1 - cos_theta * cos_theta));
  return {st * std::cos(phi), st * std::sin(phi), cos_theta};
}

/* The Stokes rotation from the meridional basis (e_v, e_h) to the
   scattering-plane basis (e_par = N x n, e_perp = N) of the ray
   n = direction(cos_theta, phi).  Both bases are right-handed with n.  With
   e_par = cos(a) e_v + sin(a) e_h:
     Q' = cos(2a) Q + sin(2a) U,  U' = -sin(2a) Q + cos(2a) U. */
inline rtepack::muelmat to_scattering_plane(const vec3& N, Numeric cos_theta, Numeric phi) {
  const Numeric st  = std::sqrt(std::max(0.0, 1 - cos_theta * cos_theta));
  const vec3    n   = direction(cos_theta, phi);
  const vec3    ev  = {cos_theta * std::cos(phi), cos_theta * std::sin(phi), -st};
  const vec3    eh  = {-std::sin(phi), std::cos(phi), 0.0};
  const vec3    par = cross(N, n);
  const Numeric c = dot(par, ev), s = dot(par, eh);
  const Numeric c2 = c * c - s * s, s2 = 2 * s * c;
  return {1, 0, 0, 0, 0, c2, s2, 0, 0, -s2, c2, 0, 0, 0, 0, 1};
}

inline rtepack::muelmat transposed(const rtepack::muelmat& m) {
  rtepack::muelmat t{0.0};
  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 4; j++) t[i, j] = m[j, i];
  return t;
}

/* The lab-frame phase matrix Z(out <- in) = L_out^T F(Theta) L_in for
   in = (ci, phii) and out = (co, phio), in the meridional basis.  Exactly
   forward or backward (|n_in x n_out| < 1e-12) the scattering plane is
   undefined.  There the plane through n_in and e_h(n_in) is used.  For
   |mu| < 1 a forward pair has the same and a backward pair the opposite
   azimuth, both rotations are then by 0 or pi, and Z = F, as in RT3
   (ROTATE_PHASE_MATRIX for SIN_SCAT = 0).  At mu = +-1 the meridional
   basis depends on the azimuth label, and Z is F rotated between the two
   bases.  For a regular F (F12 = F34 = 0 and |F22| = |F33| there) Z does
   not depend on the choice of plane. */
inline rtepack::muelmat lab_frame(const tro_matrix& F, Numeric ci, Numeric phii, Numeric co, Numeric phio) {
  const vec3    ni   = direction(ci, phii);
  const vec3    no   = direction(co, phio);
  vec3          N    = cross(ni, no);
  const Numeric norm = std::sqrt(dot(N, N));
  if (norm < 1e-12)
    N = {-std::sin(phii), std::cos(phii), 0.0};
  else
    for (auto& x : N) x /= norm;
  const auto             f = F(std::clamp(dot(ni, no), -1.0, 1.0));
  const rtepack::muelmat Fm{f.F11, f.F12, 0, 0, f.F12, f.F22, 0, 0, 0, 0, f.F33, f.F34, 0, 0, -f.F34, f.F44};
  return transposed(to_scattering_plane(N, co, phio)) * Fm * to_scattering_plane(N, ci, phii);
}
}  // namespace vdisort_test

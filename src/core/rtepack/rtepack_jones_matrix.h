#pragma once

#include <matpack.h>

#include <complex>

#include "rtepack_mueller_matrix.h"

namespace rtepack {
/** A Jones matrix
 *
 *  The linear map of the field of an incident ray to that of an outgoing ray,
 *  [E_v', E_h'] = J [E_v, E_h], with each field's components in its ray's
 *  polarization basis (v, h).  The rows are the outgoing components and the
 *  columns the incident ones.  It acts on fields, so it is neither the
 *  coherency matrix of a field nor an amplitude scattering matrix, whose
 *  element order and prefactors differ.
 *
 *  to_muelmat(J) (rtepack_multitype.h) is its Mueller matrix, and jones{M} the
 *  Jones matrix of a non-depolarizing Mueller matrix M.
 */
struct jones final : matpack::cdata_t<Complex, 2, 2> {
  //! a times the identity
  explicit constexpr jones(Complex a) noexcept : matpack::cdata_t<Complex, 2, 2>{a, Complex{}, Complex{}, a} {}

  //! The identity
  constexpr jones() noexcept : jones(Complex{1.0}) {}

  constexpr jones(Complex vv, Complex vh, Complex hv, Complex hh) noexcept
      : matpack::cdata_t<Complex, 2, 2>{vv, vh, hv, hh} {}

  /** The Jones matrix of a non-depolarizing Mueller matrix
   *
   *  M determines J up to a global phase, which is chosen to make the largest
   *  element of J real and positive, so to_muelmat(jones{M}) is M.  Throws if M is
   *  not the Mueller matrix of any Jones matrix, e.g. if it depolarizes.
   */
  explicit jones(const muelmat &M);

  //! The identity
  static constexpr jones id() { return jones{}; }

  //! Independent amplitudes v of the vertical and h of the horizontal component, e.g. Fresnel reflection
  static constexpr jones diagonal(Complex v, Complex h) { return {v, Complex{}, Complex{}, h}; }

  /** The change to the basis of the same ray turned by psi about its propagation direction
   *
   *  v' = cos(psi) v + sin(psi) h and h' = cos(psi) h - sin(psi) v, so
   *  to_muelmat(rotation(cos(psi), sin(psi))) is stokes_rotation(cos(2 psi), sin(2 psi)).
   */
  static constexpr jones rotation(Numeric cos_psi, Numeric sin_psi) { return {cos_psi, sin_psi, -sin_psi, cos_psi}; }

  constexpr jones &operator*=(Complex a) {
    for (auto &x : data) x *= a;
    return *this;
  }
};

//! The composition: first b, then a
constexpr jones operator*(const jones &a, const jones &b) {
  return {a[0, 0] * b[0, 0] + a[0, 1] * b[1, 0],
          a[0, 0] * b[0, 1] + a[0, 1] * b[1, 1],
          a[1, 0] * b[0, 0] + a[1, 1] * b[1, 0],
          a[1, 0] * b[0, 1] + a[1, 1] * b[1, 1]};
}

constexpr jones operator*(Complex a, jones b) { return b *= a; }
constexpr jones operator*(jones a, Complex b) { return a *= b; }
}  // namespace rtepack

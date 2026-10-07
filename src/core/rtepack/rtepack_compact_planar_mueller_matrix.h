#pragma once

#include <matpack.h>
#include <xml.h>

#include "rtepack_common.h"
#include "rtepack_mueller_matrix.h"

namespace rtepack {

/** The Mueller matrix of scattering by randomly oriented particles with a plane of symmetry, in compact form.

    In the scattering-plane basis (Q = I_par - I_perp, the reference plane
    being the plane of the incident and scattered directions) the scattering
    matrix of such particles has six independent elements,

        F = [[F11, F12,    0,   0],
             [F12, F22,    0,   0],
             [  0,   0,  F33, F34],
             [  0,   0, -F34, F44]],

    stored in ARTS's order [F11, F12, F22, F33, F34, F44].  Other codes use
    other orders (Evans' RT3: F11, F12, F33, F34, F22, F44), and the sign of
    F34 depends on the sign convention of V; convert by name, not by index.

    F is in the scattering-plane basis, not in a laboratory (meridional)
    basis: a laboratory-frame phase matrix is L_out^T F L_in, with L the
    Stokes rotations from the meridional bases of the two directions into the
    scattering plane.  Products with a muelmat therefore give a muelmat, and
    so does a product of two of these (it is not closed under products).
*/
struct compact_planar_muelmat final : Vector6 {
  constexpr compact_planar_muelmat() : Vector6{0., 0., 0., 0., 0., 0.} {}

  constexpr compact_planar_muelmat(Numeric f11, Numeric f12, Numeric f22, Numeric f33, Numeric f34, Numeric f44)
      : Vector6{f11, f12, f22, f33, f34, f44} {}

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  explicit constexpr compact_planar_muelmat(const T &input) noexcept : compact_planar_muelmat{} {
    assert(stdr::size(input) == 6);
    stdr::copy(input, data.begin());
  }

  [[nodiscard]] constexpr Numeric F11() const { return data[0]; }
  [[nodiscard]] constexpr Numeric F12() const { return data[1]; }
  [[nodiscard]] constexpr Numeric F22() const { return data[2]; }
  [[nodiscard]] constexpr Numeric F33() const { return data[3]; }
  [[nodiscard]] constexpr Numeric F34() const { return data[4]; }
  [[nodiscard]] constexpr Numeric F44() const { return data[5]; }

  [[nodiscard]] constexpr Numeric &F11() { return data[0]; }
  [[nodiscard]] constexpr Numeric &F12() { return data[1]; }
  [[nodiscard]] constexpr Numeric &F22() { return data[2]; }
  [[nodiscard]] constexpr Numeric &F33() { return data[3]; }
  [[nodiscard]] constexpr Numeric &F34() { return data[4]; }
  [[nodiscard]] constexpr Numeric &F44() { return data[5]; }

  //! The full 4 x 4 matrix, in the scattering-plane basis
  [[nodiscard]] constexpr muelmat expand() const {
    return {F11(), F12(), 0, 0, F12(), F22(), 0, 0, 0, 0, F33(), F34(), 0, 0, -F34(), F44()};
  }

  constexpr compact_planar_muelmat &operator*=(Numeric x) {
    for (auto &v : data) v *= x;
    return *this;
  }

  constexpr compact_planar_muelmat &operator+=(const compact_planar_muelmat &x) {
    for (Size i = 0; i < 6; i++) data[i] += x.data[i];
    return *this;
  }
};

constexpr compact_planar_muelmat operator*(compact_planar_muelmat a, Numeric x) { return a *= x; }
constexpr compact_planar_muelmat operator*(Numeric x, compact_planar_muelmat a) { return a *= x; }
constexpr compact_planar_muelmat operator+(compact_planar_muelmat a, const compact_planar_muelmat &b) { return a += b; }

//! A * F, a muelmat
constexpr muelmat operator*(const muelmat &a, const compact_planar_muelmat &f) { return a * f.expand(); }

//! F * A, a muelmat
constexpr muelmat operator*(const compact_planar_muelmat &f, const muelmat &a) { return f.expand() * a; }

//! F * G, a muelmat: the product is not of the compact form
constexpr muelmat operator*(const compact_planar_muelmat &f, const compact_planar_muelmat &g) {
  return f.expand() * g.expand();
}

using compact_planar_muelmat_vector            = matpack::data_t<compact_planar_muelmat, 1>;
using compact_planar_muelmat_vector_view       = matpack::view_t<compact_planar_muelmat, 1>;
using compact_planar_muelmat_vector_const_view = matpack::view_t<const compact_planar_muelmat, 1>;
}  // namespace rtepack

template <> struct std::formatter<rtepack::compact_planar_muelmat> : std::formatter<Vector6> {};

template <> struct xml_io_stream<rtepack::compact_planar_muelmat>
    : xml_io_stream_inherit<Vector6, rtepack::compact_planar_muelmat> {};

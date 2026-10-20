#pragma once

#include <matpack.h>
#include <xml.h>

#include "rtepack_common.h"
#include "rtepack_propagation_matrix.h"

namespace rtepack {

/** The propagation matrix of a medium with a plane of symmetry through the propagation direction, in compact form.

    In the basis of that plane (Q = I_par - I_perp, the reference plane being
    the plane of symmetry) the propagation matrix of such a medium has three
    independent elements,

        K = [[A, B,  0, 0],
             [B, A,  0, 0],
             [0, 0,  A, W],
             [0, 0, -W, A]],

    the A, B and W of a propmat, stored in that order.  For azimuthally
    randomly oriented particles they are K11, K12 and K34 in the laboratory
    (meridional) basis; for totally randomly oriented particles B = W = 0.
*/
struct compact_planar_propmat final : Vector3 {
  constexpr compact_planar_propmat() : Vector3{0., 0., 0.} {}

  explicit constexpr compact_planar_propmat(Numeric a) : Vector3{a, 0., 0.} {}

  constexpr compact_planar_propmat(Numeric a, Numeric b, Numeric w) : Vector3{a, b, w} {}

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  explicit constexpr compact_planar_propmat(const T &input) noexcept : compact_planar_propmat{} {
    assert(stdr::size(input) == 3);
    stdr::copy(input, data.begin());
  }

  [[nodiscard]] constexpr Numeric A() const { return data[0]; }
  [[nodiscard]] constexpr Numeric B() const { return data[1]; }
  [[nodiscard]] constexpr Numeric W() const { return data[2]; }

  [[nodiscard]] constexpr Numeric &A() { return data[0]; }
  [[nodiscard]] constexpr Numeric &B() { return data[1]; }
  [[nodiscard]] constexpr Numeric &W() { return data[2]; }

  //! The full propagation matrix, in the basis of the plane of symmetry
  [[nodiscard]] constexpr propmat expand() const { return {A(), B(), 0, 0, 0, 0, W()}; }

  constexpr compact_planar_propmat &operator*=(Numeric x) {
    for (auto &v : data) v *= x;
    return *this;
  }

  constexpr compact_planar_propmat &operator/=(Numeric x) {
    for (auto &v : data) v /= x;
    return *this;
  }

  constexpr compact_planar_propmat &operator+=(const compact_planar_propmat &x) {
    for (Size i = 0; i < 3; i++) data[i] += x.data[i];
    return *this;
  }

  constexpr compact_planar_propmat &operator-=(const compact_planar_propmat &x) {
    for (Size i = 0; i < 3; i++) data[i] -= x.data[i];
    return *this;
  }
};

constexpr compact_planar_propmat operator*(compact_planar_propmat a, Numeric x) { return a *= x; }
constexpr compact_planar_propmat operator*(Numeric x, compact_planar_propmat a) { return a *= x; }
constexpr compact_planar_propmat operator/(compact_planar_propmat a, Numeric x) { return a /= x; }
constexpr compact_planar_propmat operator+(compact_planar_propmat a, const compact_planar_propmat &b) { return a += b; }
constexpr compact_planar_propmat operator-(compact_planar_propmat a, const compact_planar_propmat &b) { return a -= b; }

using compact_planar_propmat_vector             = matpack::data_t<compact_planar_propmat, 1>;
using compact_planar_propmat_vector_view        = matpack::view_t<compact_planar_propmat, 1>;
using compact_planar_propmat_vector_const_view  = matpack::view_t<const compact_planar_propmat, 1>;
using compact_planar_propmat_matrix             = matpack::data_t<compact_planar_propmat, 2>;
using compact_planar_propmat_matrix_view        = matpack::view_t<compact_planar_propmat, 2>;
using compact_planar_propmat_matrix_const_view  = matpack::view_t<const compact_planar_propmat, 2>;
using compact_planar_propmat_tensor3            = matpack::data_t<compact_planar_propmat, 3>;
using compact_planar_propmat_tensor3_view       = matpack::view_t<compact_planar_propmat, 3>;
using compact_planar_propmat_tensor3_const_view = matpack::view_t<const compact_planar_propmat, 3>;
}  // namespace rtepack

template <> struct std::formatter<rtepack::compact_planar_propmat> : std::formatter<Vector3> {};

template <> struct xml_io_stream<rtepack::compact_planar_propmat>
    : xml_io_stream_inherit<Vector3, rtepack::compact_planar_propmat> {};

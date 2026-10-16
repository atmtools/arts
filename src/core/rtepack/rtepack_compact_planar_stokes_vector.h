#pragma once

#include <matpack.h>
#include <xml.h>

#include "rtepack_common.h"
#include "rtepack_stokes_vector.h"

namespace rtepack {

/** The Stokes vector of a medium with a plane of symmetry through the propagation direction, in compact form.

    In the basis of that plane (Q = I_par - I_perp, the reference plane being
    the plane of symmetry) U = V = 0, so only

        s = [I, Q, 0, 0]

    the I and Q of a stokvec are stored, in that order.  For the absorption
    vector of azimuthally randomly oriented particles they are a1 and a2 in
    the laboratory (meridional) basis; for totally randomly oriented particles
    Q = 0.
*/
struct compact_planar_stokvec final : Vector2 {
  constexpr compact_planar_stokvec() : Vector2{0., 0.} {}

  explicit constexpr compact_planar_stokvec(Numeric i) : Vector2{i, 0.} {}

  constexpr compact_planar_stokvec(Numeric i, Numeric q) : Vector2{i, q} {}

  template <stdr::forward_range T> requires std::is_convertible_v<stdr::range_value_t<T>, Numeric>
  explicit constexpr compact_planar_stokvec(const T &input) noexcept : compact_planar_stokvec{} {
    assert(stdr::size(input) == 2);
    stdr::copy(input, data.begin());
  }

  [[nodiscard]] constexpr Numeric I() const { return data[0]; }
  [[nodiscard]] constexpr Numeric Q() const { return data[1]; }

  [[nodiscard]] constexpr Numeric &I() { return data[0]; }
  [[nodiscard]] constexpr Numeric &Q() { return data[1]; }

  //! The full Stokes vector, in the basis of the plane of symmetry
  [[nodiscard]] constexpr stokvec expand() const { return {I(), Q(), 0, 0}; }

  constexpr compact_planar_stokvec &operator*=(Numeric x) {
    for (auto &v : data) v *= x;
    return *this;
  }

  constexpr compact_planar_stokvec &operator/=(Numeric x) {
    for (auto &v : data) v /= x;
    return *this;
  }

  constexpr compact_planar_stokvec &operator+=(const compact_planar_stokvec &x) {
    for (Size i = 0; i < 2; i++) data[i] += x.data[i];
    return *this;
  }

  constexpr compact_planar_stokvec &operator-=(const compact_planar_stokvec &x) {
    for (Size i = 0; i < 2; i++) data[i] -= x.data[i];
    return *this;
  }
};

constexpr compact_planar_stokvec operator*(compact_planar_stokvec a, Numeric x) { return a *= x; }
constexpr compact_planar_stokvec operator*(Numeric x, compact_planar_stokvec a) { return a *= x; }
constexpr compact_planar_stokvec operator/(compact_planar_stokvec a, Numeric x) { return a /= x; }
constexpr compact_planar_stokvec operator+(compact_planar_stokvec a, const compact_planar_stokvec &b) { return a += b; }
constexpr compact_planar_stokvec operator-(compact_planar_stokvec a, const compact_planar_stokvec &b) { return a -= b; }

using compact_planar_stokvec_vector             = matpack::data_t<compact_planar_stokvec, 1>;
using compact_planar_stokvec_vector_view        = matpack::view_t<compact_planar_stokvec, 1>;
using compact_planar_stokvec_vector_const_view  = matpack::view_t<const compact_planar_stokvec, 1>;
using compact_planar_stokvec_matrix             = matpack::data_t<compact_planar_stokvec, 2>;
using compact_planar_stokvec_matrix_view        = matpack::view_t<compact_planar_stokvec, 2>;
using compact_planar_stokvec_matrix_const_view  = matpack::view_t<const compact_planar_stokvec, 2>;
using compact_planar_stokvec_tensor3            = matpack::data_t<compact_planar_stokvec, 3>;
using compact_planar_stokvec_tensor3_view       = matpack::view_t<compact_planar_stokvec, 3>;
using compact_planar_stokvec_tensor3_const_view = matpack::view_t<const compact_planar_stokvec, 3>;
}  // namespace rtepack

template <> struct std::formatter<rtepack::compact_planar_stokvec> : std::formatter<Vector2> {};

template <> struct xml_io_stream<rtepack::compact_planar_stokvec>
    : xml_io_stream_inherit<Vector2, rtepack::compact_planar_stokvec> {};

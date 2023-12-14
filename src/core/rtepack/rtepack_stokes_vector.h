#pragma once

#include "rtepack_concepts.h"

namespace rtepack {
struct stokvec final : vec4 {
  [[nodiscard]] constexpr stokvec(Numeric a = 0.0,
                                  Numeric b = 0.0,
                                  Numeric c = 0.0,
                                  Numeric d = 0.0)
      : vec4{a, b, c, d} {}

  [[nodiscard]] constexpr Numeric I() const { return data[0]; }
  [[nodiscard]] constexpr Numeric Q() const { return data[1]; }
  [[nodiscard]] constexpr Numeric U() const { return data[2]; }
  [[nodiscard]] constexpr Numeric V() const { return data[3]; }

  [[nodiscard]] constexpr Numeric &I() { return data[0]; }
  [[nodiscard]] constexpr Numeric &Q() { return data[1]; }
  [[nodiscard]] constexpr Numeric &U() { return data[2]; }
  [[nodiscard]] constexpr Numeric &V() { return data[3]; }
};

//! Addition of two stokvec vectors
constexpr auto operator+(const stokvec &a, const stokvec &b) {
  stokvec c{a};
  c += b;
  return c;
}

//! Subtraction between two stokvec vectors
constexpr auto operator-(const stokvec &a, const stokvec &b) {
  stokvec c{a};
  c -= b;
  return c;
}

//! Scaling a stokvec vector
constexpr auto operator*(const Numeric &a, const stokvec &b) {
  stokvec c{b};
  c *= a;
  return c;
}

//! Scaling a stokvec vector
constexpr auto operator*(const stokvec &a, const Numeric &b) {
  stokvec c{a};
  c *= b;
  return c;
}

//! Take the average of two stokvec vectors
constexpr auto avg(const stokvec &a, const stokvec &b) {
  return 0.5 * a + 0.5 * b;
}

using stokvec_vector = matpack::matpack_data<stokvec, 1>;
using stokvec_vector_view = matpack::matpack_view<stokvec, 1, false, false>;
using stokvec_vector_const_view =
    matpack::matpack_view<stokvec, 1, true, false>;

using stokvec_matrix = matpack::matpack_data<stokvec, 2>;
using stokvec_matrix_view = matpack::matpack_view<stokvec, 2, false, false>;
using stokvec_matrix_const_view =
    matpack::matpack_view<stokvec, 2, true, false>;

using stokvec_tensor3 = matpack::matpack_data<stokvec, 3>;
using stokvec_tensor3_view = matpack::matpack_view<stokvec, 3, false, false>;
using stokvec_tensor3_const_view =
    matpack::matpack_view<stokvec, 3, true, false>;
}  // namespace rtepack

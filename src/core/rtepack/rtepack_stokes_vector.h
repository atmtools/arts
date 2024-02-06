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

  constexpr stokvec& operator += (const stokvec &b) {
    I() += b.I();
    Q() += b.Q();
    U() += b.U();
    V() += b.V();
    return *this;
  }

  constexpr stokvec& operator -= (const stokvec &b) {
    I() -= b.I();
    Q() -= b.Q();
    U() -= b.U();
    V() -= b.V();
    return *this;
  }
};

//! Addition of two stokvec vectors
constexpr stokvec operator+(stokvec a, const stokvec &b) {
  a += b;
  return a;
}

//! Subtraction between two stokvec vectors
constexpr stokvec operator-(stokvec a, const stokvec &b) {
  a -= b;
  return a;
}

//! Scaling a stokvec vector
constexpr stokvec operator*(const Numeric &a, stokvec b) {
  b *= a;
  return b;
}

//! Scaling a stokvec vector
constexpr stokvec operator*(stokvec a, const Numeric &b) {
  a *= b;
  return a;
}

//! Take the average of two stokvec vectors
constexpr stokvec avg(const stokvec &a, const stokvec &b) {
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

using stokvec_tensor4 = matpack::matpack_data<stokvec, 4>;
using stokvec_tensor4_view = matpack::matpack_view<stokvec, 4, false, false>;
using stokvec_tensor4_const_view =
    matpack::matpack_view<stokvec, 4, true, false>;
}  // namespace rtepack

#pragma once

#include <enums.h>

#include <utility>

#include "rtepack_concepts.h"

namespace rtepack {
struct stokvec final : vec4 {
  [[nodiscard]] constexpr stokvec(Numeric a = 0.0,
                                  Numeric b = 0.0,
                                  Numeric c = 0.0,
                                  Numeric d = 0.0)
      : vec4{a, b, c, d} {}

  constexpr stokvec(std::array<Numeric, 4> data) noexcept : vec4{data} {}

  constexpr stokvec &operator=(Numeric a) {
    data = {a, 0., 0., 0.};
    return *this;
  }

  [[nodiscard]] constexpr Numeric I() const { return data[0]; }
  [[nodiscard]] constexpr Numeric Q() const { return data[1]; }
  [[nodiscard]] constexpr Numeric U() const { return data[2]; }
  [[nodiscard]] constexpr Numeric V() const { return data[3]; }

  [[nodiscard]] constexpr Numeric &I() { return data[0]; }
  [[nodiscard]] constexpr Numeric &Q() { return data[1]; }
  [[nodiscard]] constexpr Numeric &U() { return data[2]; }
  [[nodiscard]] constexpr Numeric &V() { return data[3]; }

  constexpr stokvec &operator+=(const stokvec &b) {
    I() += b.I();
    Q() += b.Q();
    U() += b.U();
    V() += b.V();
    return *this;
  }

  constexpr stokvec &operator-=(const stokvec &b) {
    I() -= b.I();
    Q() -= b.Q();
    U() -= b.U();
    V() -= b.V();
    return *this;
  }

  [[nodiscard]] constexpr bool is_zero() const {
    return I() == 0.0 && Q() == 0.0 && U() == 0.0 && V() == 0.0;
  }
};

constexpr stokvec to_stokvec(PolarizationChoice p) {
  using enum PolarizationChoice;
  switch (p) {
    case I:
      return {1.0, 0.0, 0.0, 0.0};
    case Q:
      return {0.0, 1.0, 0.0, 0.0};
    case U:
      return {0.0, 0.0, 1.0, 0.0};
    case V:
      return {0.0, 0.0, 0.0, 1.0};
    case Iv:
      return {1.0, 1.0, 0.0, 0.0};
    case Ih:
      return {1.0, -1.0, 0.0, 0.0};
    case Ip45:
      return {1.0, 0.0, 1.0, 0.0};
    case Im45:
      return {1.0, 0.0, -1.0, 0.0};
    case Ilhc:
      return {1.0, 0.0, 0.0, -1.0};
    case Irhc:
      return {1.0, 0.0, 0.0, 1.0};
  }
  std::unreachable();
}

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

/** Convertion methods for stokes vectors
 * 
  * @param type Spectral Radiance Unit Type
 * @return A function that takes a stokes vector and a frequency and returns a converted stokes vector
 */
std::function<stokvec(const stokvec, const Numeric)> unit_converter(
    const SpectralRadianceUnitType type, const Numeric n = 1.0);
std::function<stokvec(const stokvec, const stokvec, const Numeric)>
dunit_converter(const SpectralRadianceUnitType type, const Numeric n = 1.0);

using stokvec_vector      = matpack::matpack_data<stokvec, 1>;
using stokvec_vector_view = matpack::matpack_view<stokvec, 1, false, false>;
using stokvec_vector_const_view =
    matpack::matpack_view<stokvec, 1, true, false>;

using stokvec_matrix      = matpack::matpack_data<stokvec, 2>;
using stokvec_matrix_view = matpack::matpack_view<stokvec, 2, false, false>;
using stokvec_matrix_const_view =
    matpack::matpack_view<stokvec, 2, true, false>;

using stokvec_tensor3      = matpack::matpack_data<stokvec, 3>;
using stokvec_tensor3_view = matpack::matpack_view<stokvec, 3, false, false>;
using stokvec_tensor3_const_view =
    matpack::matpack_view<stokvec, 3, true, false>;

using stokvec_tensor4      = matpack::matpack_data<stokvec, 4>;
using stokvec_tensor4_view = matpack::matpack_view<stokvec, 4, false, false>;
using stokvec_tensor4_const_view =
    matpack::matpack_view<stokvec, 4, true, false>;

using stokvec_tensor5      = matpack::matpack_data<stokvec, 5>;
using stokvec_tensor5_view = matpack::matpack_view<stokvec, 5, false, false>;
using stokvec_tensor5_const_view =
    matpack::matpack_view<stokvec, 5, true, false>;

using stokvec_tensor6      = matpack::matpack_data<stokvec, 6>;
using stokvec_tensor6_view = matpack::matpack_view<stokvec, 6, false, false>;
using stokvec_tensor6_const_view =
    matpack::matpack_view<stokvec, 6, true, false>;
}  // namespace rtepack

template <>

struct std::formatter<rtepack::stokvec> {
  std::formatter<rtepack::vec4> fmt;

  [[nodiscard]] constexpr auto &inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto &inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const rtepack::stokvec &v,
                              FmtContext &ctx) const {
    return fmt.format(v, ctx);
  }
};

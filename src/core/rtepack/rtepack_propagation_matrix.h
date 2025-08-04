#pragma once

#include <compare.h>
#include <xml.h>

#include "rtepack_concepts.h"

namespace rtepack {
struct propmat final : vec7 {
  constexpr propmat(Numeric a = 0.0,
                    Numeric b = 0.0,
                    Numeric c = 0.0,
                    Numeric d = 0.0,
                    Numeric u = 0.0,
                    Numeric v = 0.0,
                    Numeric w = 0.0)
      : vec7{a, b, c, d, u, v, w} {}

  constexpr propmat(std::array<Numeric, 7> data) noexcept : vec7{data} {}

  [[nodiscard]] constexpr decltype(auto) A() const { return data[0]; }
  [[nodiscard]] constexpr decltype(auto) B() const { return data[1]; }
  [[nodiscard]] constexpr decltype(auto) C() const { return data[2]; }
  [[nodiscard]] constexpr decltype(auto) D() const { return data[3]; }
  [[nodiscard]] constexpr decltype(auto) U() const { return data[4]; }
  [[nodiscard]] constexpr decltype(auto) V() const { return data[5]; }
  [[nodiscard]] constexpr decltype(auto) W() const { return data[6]; }

  [[nodiscard]] constexpr decltype(auto) A() { return data[0]; }
  [[nodiscard]] constexpr decltype(auto) B() { return data[1]; }
  [[nodiscard]] constexpr decltype(auto) C() { return data[2]; }
  [[nodiscard]] constexpr decltype(auto) D() { return data[3]; }
  [[nodiscard]] constexpr decltype(auto) U() { return data[4]; }
  [[nodiscard]] constexpr decltype(auto) V() { return data[5]; }
  [[nodiscard]] constexpr decltype(auto) W() { return data[6]; }

  //! Check if the matrix is purely rotational
  [[nodiscard]] constexpr bool is_rotational() const { return A() == 0.0; }

  //! Check if the matrix is polarized
  [[nodiscard]] bool is_polarized() const;

  constexpr auto operator<=>(const propmat &pm) const { return A() <=> pm.A(); }
};

//! Addition of two propmat matrixes
constexpr propmat operator+(propmat a, const propmat &b) {
  a += b;
  return a;
}

//! Subtraction between two propmat matrixes
constexpr propmat operator-(propmat a, const propmat &b) {
  a -= b;
  return a;
}

//! Scaling a propmat matrix
constexpr propmat operator*(propmat a, const Numeric &b) {
  a *= b;
  return a;
}

//! Scaling a propmat matrix
constexpr propmat operator*(const Numeric &a, propmat b) { return b * a; }

//! Scaling a propmat matrix
constexpr propmat operator/(propmat a, const propmat &b) {
  a /= b;
  return a;
}

//! Take the average of two propmat matrixes
constexpr propmat avg(propmat a, const propmat &b) {
  a += b;
  a *= 0.5;
  return a;
}

using propmat_vector            = matpack::data_t<propmat, 1>;
using propmat_vector_view       = matpack::view_t<propmat, 1>;
using propmat_vector_const_view = matpack::view_t<const propmat, 1>;

using propmat_matrix            = matpack::data_t<propmat, 2>;
using propmat_matrix_view       = matpack::view_t<propmat, 2>;
using propmat_matrix_const_view = matpack::view_t<const propmat, 2>;

using propmat_tensor3            = matpack::data_t<propmat, 3>;
using propmat_tensor3_view       = matpack::view_t<propmat, 3>;
using propmat_tensor3_const_view = matpack::view_t<const propmat, 3>;

propmat_vector operator*(Numeric x, const propmat_vector_const_view &y);
}  // namespace rtepack

template <>

struct std::formatter<rtepack::propmat> {
  std::formatter<rtepack::vec7> fmt;

  [[nodiscard]] constexpr auto &inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto &inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const rtepack::propmat &v,
                              FmtContext &ctx) const {
    return fmt.format(v, ctx);
  }
};

template <>
struct xml_io_stream<rtepack::propmat> {
  static constexpr std::string_view type_name = "Propmat"sv;

  static void write(std::ostream &os,
                    const rtepack::propmat &x,
                    bofstream *pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream &is,
                   rtepack::propmat &x,
                   bifstream *pbifs = nullptr);
  static void put(std::span<const rtepack::propmat> x, bofstream *);
  static void get(std::span<rtepack::propmat> x, bifstream *pbifs);
  static void parse(std::span<rtepack::propmat> x, std::istream &);
};

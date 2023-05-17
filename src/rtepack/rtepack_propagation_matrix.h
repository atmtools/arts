#pragma once

#include "debug.h"
#include "matpack_data.h"
#include "matpack_lazy.h"
#include "rtepack_concepts.h"

#include <type_traits>

namespace rtepack {
using namespace matpack::lazy;

template <typename T>
concept lazy_propmat =
    constexpr_vec_data_like<T> and std::remove_cvref_t<T>::size() == 7;

template <typename T>
concept propmat_convertible =
    matpack::column_keeper<T> and matpack::row_keeper<T> and
    matpack::rank<T>() == 2 and matpack::mdvalue_type_compatible<T, Numeric>;

struct propmat final : vec7 {
  constexpr propmat(Numeric a = 0.0, Numeric b = 0.0, Numeric c = 0.0,
                    Numeric d = 0.0, Numeric u = 0.0, Numeric v = 0.0,
                    Numeric w = 0.0)
      : vec7{a, b, c, d, u, v, w} {}

  constexpr propmat(lazy_propmat auto &&a)
    requires(std::is_rvalue_reference_v<decltype(a)>)
      : vec7(std::move(a)) {}

  template <propmat_convertible T> propmat(const T &t) {
    ARTS_USER_ERROR_IF(matpack::column_size(t) != 4,
                       "The matrix must have 4 columns.")
    ARTS_USER_ERROR_IF(matpack::row_size(t) != 4,
                       "The matrix must have 4 rows.")

    ARTS_USER_ERROR_IF(
        matpack::mdvalue(t, std::array<Index, 2>{0, 0}) not_eq
        matpack::mdvalue(t, std::array<Index, 2>{1, 1}) or
            matpack::mdvalue(t, std::array<Index, 2>{0, 0}) not_eq
            matpack::mdvalue(t, std::array<Index, 2>{2, 2}) or
                matpack::mdvalue(t, std::array<Index, 2>{0, 0}) not_eq
                matpack::mdvalue(t, std::array<Index, 2>{3, 3}),
        "Must have same diagonal value")

    ARTS_USER_ERROR_IF(matpack::mdvalue(t, std::array<Index, 2>{0, 1}) not_eq
                       matpack::mdvalue(t, std::array<Index, 2>{1, 0}),
                       "Must be 10 - 01 symmetric")

    ARTS_USER_ERROR_IF(matpack::mdvalue(t, std::array<Index, 2>{0, 2}) not_eq
                       matpack::mdvalue(t, std::array<Index, 2>{2, 0}),
                       "Must be 20 - 02 symmetric")

    ARTS_USER_ERROR_IF(matpack::mdvalue(t, std::array<Index, 2>{0, 3}) not_eq
                       matpack::mdvalue(t, std::array<Index, 2>{3, 0}),
                       "Must be 30 - 03 symmetric")

    ARTS_USER_ERROR_IF(matpack::mdvalue(t, std::array<Index, 2>{1, 2}) not_eq
                      -matpack::mdvalue(t, std::array<Index, 2>{2, 1}),
                       "Must be 12 - 21 asymmetric")

    ARTS_USER_ERROR_IF(matpack::mdvalue(t, std::array<Index, 2>{1, 3}) not_eq
                      -matpack::mdvalue(t, std::array<Index, 2>{3, 1}),
                       "Must be 13 - 31 asymmetric")

    ARTS_USER_ERROR_IF(matpack::mdvalue(t, std::array<Index, 2>{2, 3}) not_eq
                      -matpack::mdvalue(t, std::array<Index, 2>{3, 2}),
                       "Must be 23 - 32 asymmetric")

    data[0] = matpack::mdvalue(t, std::array<Index, 2>{0, 0});
    data[1] = matpack::mdvalue(t, std::array<Index, 2>{0, 1});
    data[2] = matpack::mdvalue(t, std::array<Index, 2>{0, 2});
    data[3] = matpack::mdvalue(t, std::array<Index, 2>{0, 3});
    data[4] = matpack::mdvalue(t, std::array<Index, 2>{1, 2});
    data[5] = matpack::mdvalue(t, std::array<Index, 2>{1, 3});
    data[6] = matpack::mdvalue(t, std::array<Index, 2>{2, 3});
  }

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

  [[nodiscard]] constexpr bool is_rotational() const { return A() == 0.0; }

  constexpr auto operator<=>(const propmat &pm) const { return A() <=> pm.A(); }
};

//! Addition of two propmat matrixes
constexpr auto operator+(const propmat &a, const propmat &b) {
  return constexpr_vec_data{a} + constexpr_vec_data{b};
}

//! Addition of a propmat matrix and a lazy one
constexpr auto operator+(const propmat &a, const lazy_propmat auto &b) {
  return constexpr_vec_data{a} + b;
}

//! Addition of a propmat matrix and a lazy one
constexpr auto operator+(const lazy_propmat auto &a, const propmat &b) {
  return a + constexpr_vec_data{b};
}

//! Subtraction between two propmat matrixes
constexpr auto operator-(const propmat &a, const propmat &b) {
  return constexpr_vec_data{a} - constexpr_vec_data{b};
}

//! Subtraction between a propmat matrix and a lazy one
constexpr auto operator-(const propmat &a, const lazy_propmat auto &b) {
  return constexpr_vec_data{a} - b;
}

//! Subtraction between a propmat matrix and a lazy one
constexpr auto operator-(const lazy_propmat auto &a, const propmat &b) {
  return a - constexpr_vec_data{b};
}

//! Scaling a propmat matrix
constexpr auto operator*(const Numeric &a, const propmat &b) {
  return vscl{a, constexpr_vec_data{b}};
}

//! Scaling a propmat matrix
constexpr auto operator*(const propmat &a, const Numeric &b) {
  return vscl{b, constexpr_vec_data{a}};
}

//! Scaling a propmat matrix
constexpr propmat operator/(const propmat &a, const propmat &b) {
  auto c = a;
  c /= b;
  return c;
}

//! Take the average of two propmat matrixes
constexpr auto avg(const propmat &a, const propmat &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a propmat matrix and a lazy propmat matrix
constexpr auto avg(const propmat &a, const lazy_propmat auto &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a propmat matrix and a lazy propmat matrix
constexpr auto avg(const lazy_propmat auto &a, const propmat &b) {
  return 0.5 * a + 0.5 * b;
}

using propmat_vector = matpack::matpack_data<propmat, 1>;
using propmat_vector_view = matpack::matpack_view<propmat, 1, false, false>;
using propmat_vector_const_view =
    matpack::matpack_view<propmat, 1, true, false>;

using propmat_matrix = matpack::matpack_data<propmat, 2>;
using propmat_matrix_view = matpack::matpack_view<propmat, 2, false, false>;
using propmat_matrix_const_view =
    matpack::matpack_view<propmat, 2, true, false>;

using propmat_tensor3 = matpack::matpack_data<propmat, 3>;
using propmat_tensor3_view = matpack::matpack_view<propmat, 3, false, false>;
using propmat_tensor3_const_view =
    matpack::matpack_view<propmat, 3, true, false>;
} // namespace rtepack

#pragma once

#include "matpack_data.h"
#include "matpack_lazy.h"
#include "rtepack_concepts.h"

#include <type_traits>

namespace rtepack {
using namespace matpack::lazy;

template <typename T>
concept lazy_propmat =
    constexpr_vec_data_like<T> and std::remove_cvref_t<T>::size() == 7;

struct propmat final : vec7 {
  constexpr propmat(Numeric a = 0.0, Numeric b = 0.0, Numeric c = 0.0,
                    Numeric d = 0.0, Numeric u = 0.0, Numeric v = 0.0,
                    Numeric w = 0.0)
      : vec7{a, b, c, d, u, v, w} {}

  constexpr propmat(lazy_propmat auto &&a)
    requires(std::is_rvalue_reference_v<decltype(a)>)
      : vec7(std::move(a)) {}

  [[nodiscard]] constexpr Numeric A() const { return this->operator[](0); }
  [[nodiscard]] constexpr Numeric B() const { return this->operator[](1); }
  [[nodiscard]] constexpr Numeric C() const { return this->operator[](2); }
  [[nodiscard]] constexpr Numeric D() const { return this->operator[](3); }
  [[nodiscard]] constexpr Numeric U() const { return this->operator[](4); }
  [[nodiscard]] constexpr Numeric V() const { return this->operator[](5); }
  [[nodiscard]] constexpr Numeric W() const { return this->operator[](6); }

  [[nodiscard]] constexpr Numeric &A() { return this->operator[](0); }
  [[nodiscard]] constexpr Numeric &B() { return this->operator[](1); }
  [[nodiscard]] constexpr Numeric &C() { return this->operator[](2); }
  [[nodiscard]] constexpr Numeric &D() { return this->operator[](3); }
  [[nodiscard]] constexpr Numeric &U() { return this->operator[](4); }
  [[nodiscard]] constexpr Numeric &V() { return this->operator[](5); }
  [[nodiscard]] constexpr Numeric &W() { return this->operator[](6); }

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
using propmat_vector_const_view = matpack::matpack_view<propmat, 1, true, false>;

using propmat_matrix = matpack::matpack_data<propmat, 2>;
using propmat_matrix_view = matpack::matpack_view<propmat, 2, false, false>;
using propmat_matrix_const_view = matpack::matpack_view<propmat, 2, true, false>;

using propmat_tensor3 = matpack::matpack_data<propmat, 3>;
using propmat_tensor3_view = matpack::matpack_view<propmat, 3, false, false>;
using propmat_tensor3_const_view = matpack::matpack_view<propmat, 3, true, false>;
} // namespace rtepack



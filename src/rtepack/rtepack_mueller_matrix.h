#pragma once

#include "matpack_lazy.h"
#include "rtepack_concepts.h"

#include <type_traits>

namespace rtepack {
using namespace matpack::lazy;

template <typename T>
concept lazy_muelmat =
    constexpr_smat_data_like<T> and std::remove_cvref_t<T>::size() == 4;

struct muelmat final : mat44 {
  constexpr muelmat(Numeric tau = 1.0) noexcept
      : mat44{tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau} {}

  constexpr muelmat(lazy_muelmat auto &&a)
    requires(std::is_rvalue_reference_v<decltype(a)>)
      : mat44(static_cast<mat44>(std::move(a))) {}

  constexpr muelmat(Numeric a, Numeric b, Numeric c, Numeric d, Numeric e,
                    Numeric f, Numeric g, Numeric h, Numeric i, Numeric j,
                    Numeric k, Numeric l, Numeric m, Numeric n, Numeric o,
                    Numeric p) noexcept
      : mat44{a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p} {}
};

//! Addition between muelmat matrices
constexpr auto operator+(const muelmat &a, const muelmat &b) {
  return constexpr_smat_data{a} + constexpr_smat_data{b};
}

//! Addition between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator+(const lazy_muelmat auto &a, const muelmat &b) {
  return a + constexpr_smat_data{b};
}

//! Addition between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator+(const muelmat &a, const lazy_muelmat auto &b) {
  return constexpr_smat_data{a} + b;
}

//! Subtraction between muelmat matrices
constexpr auto operator-(const muelmat &a, const muelmat &b) {
  return constexpr_smat_data{a} - constexpr_smat_data{b};
}

//! Subtraction between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator-(const lazy_muelmat auto &a, const muelmat &b) {
  return a - constexpr_smat_data{b};
}

//! Subtraction between muelmat and lazy types (lazy @ muelmat)
constexpr auto operator-(const muelmat &a, const lazy_muelmat auto &b) {
  return constexpr_smat_data{a} - b;
}

//! Scaling a muelmat matrix
constexpr auto operator*(const Numeric &a, const muelmat &b) {
  return smscl{a, constexpr_smat_data{b}};
}

//! Scaling a muelmat matrix
constexpr auto operator*(const muelmat &a, const Numeric &b) {
  return smscl{b, constexpr_smat_data{a}};
}

//! Take the average of two muelmat matrices
constexpr auto avg(const muelmat &a, const muelmat &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a muelmat matrix and a non-muelmat matrix
constexpr auto avg(const muelmat &a, const lazy_muelmat auto &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a muelmat matrix and a non-muelmat matrix
constexpr auto avg(const lazy_muelmat auto &a, const muelmat &b) {
  return 0.5 * a + 0.5 * b;
}

using muelmat_vector = matpack::matpack_data<muelmat, 1>;
using muelmat_vector_view = matpack::matpack_view<muelmat, 1, false, false>;
using muelmat_vector_const_view = matpack::matpack_view<muelmat, 1, true, false>;

using muelmat_matrix = matpack::matpack_data<muelmat, 2>;
using muelmat_matrix_view = matpack::matpack_view<muelmat, 2, false, false>;
using muelmat_matrix_const_view = matpack::matpack_view<muelmat, 2, true, false>;

using muelmat_tensor3 = matpack::matpack_data<muelmat, 3>;
using muelmat_tensor3_view = matpack::matpack_view<muelmat, 3, false, false>;
using muelmat_tensor3_const_view = matpack::matpack_view<muelmat, 3, true, false>;
} // namespace rtepack

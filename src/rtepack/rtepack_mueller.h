#pragma once

#include "matpack_lazy.h"
#include "rtepack_concepts.h"

#include <type_traits>

namespace rtepack {
using namespace matpack::lazy;

template <typename T>
concept lazy_mueller =
    constexpr_smat_data_like<T> and std::remove_cvref_t<T>::size() == 4;

struct mueller final : mat44 {
  constexpr mueller(Numeric tau = 1.0) noexcept
      : mat44{tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau, 0, 0, 0, 0, tau} {}

  constexpr mueller(lazy_mueller auto &&a)
    requires(std::is_rvalue_reference_v<decltype(a)>)
      : mat44(static_cast<mat44>(std::move(a))) {}

  constexpr mueller(Numeric a, Numeric b, Numeric c, Numeric d, Numeric e,
                    Numeric f, Numeric g, Numeric h, Numeric i, Numeric j,
                    Numeric k, Numeric l, Numeric m, Numeric n, Numeric o,
                    Numeric p) noexcept
      : mat44{a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p} {}
};

//! Addition between mueller matrices
constexpr auto operator+(const mueller &a, const mueller &b) {
  return constexpr_smat_data{a} + constexpr_smat_data{b};
}

//! Addition between mueller and lazy types (lazy @ mueller)
constexpr auto operator+(const lazy_mueller auto &a, const mueller &b) {
  return a + constexpr_smat_data{b};
}

//! Addition between mueller and lazy types (lazy @ mueller)
constexpr auto operator+(const mueller &a, const lazy_mueller auto &b) {
  return constexpr_smat_data{a} + b;
}

//! Subtraction between mueller matrices
constexpr auto operator-(const mueller &a, const mueller &b) {
  return constexpr_smat_data{a} - constexpr_smat_data{b};
}

//! Subtraction between mueller and lazy types (lazy @ mueller)
constexpr auto operator-(const lazy_mueller auto &a, const mueller &b) {
  return a - constexpr_smat_data{b};
}

//! Subtraction between mueller and lazy types (lazy @ mueller)
constexpr auto operator-(const mueller &a, const lazy_mueller auto &b) {
  return constexpr_smat_data{a} - b;
}

//! Scaling a mueller matrix
constexpr auto operator*(const Numeric &a, const mueller &b) {
  return smscl{a, constexpr_smat_data{b}};
}

//! Scaling a mueller matrix
constexpr auto operator*(const mueller &a, const Numeric &b) {
  return smscl{b, constexpr_smat_data{a}};
}

//! Take the average of two mueller matrices
constexpr auto avg(const mueller &a, const mueller &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a mueller matrix and a non-mueller matrix
constexpr auto avg(const mueller &a, const lazy_mueller auto &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a mueller matrix and a non-mueller matrix
constexpr auto avg(const lazy_mueller auto &a, const mueller &b) {
  return 0.5 * a + 0.5 * b;
}

using mueller_view = matpack::matpack_view<mueller, 1, false, false>;
using const_mueller_view = matpack::matpack_view<mueller, 1, true, false>;
} // namespace rtepack

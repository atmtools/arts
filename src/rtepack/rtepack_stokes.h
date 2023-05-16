#pragma once

#include "matpack_lazy.h"
#include "rtepack_concepts.h"

#include <type_traits>

namespace rtepack {
using namespace matpack::lazy;

template <typename T>
concept lazy_stokes =
    constexpr_vec_data_like<T> and std::remove_cvref_t<T>::size() == 4;

struct stokes final : vec4 {
  constexpr stokes(Numeric a = 0.0, Numeric b = 0.0, Numeric c = 0.0,
                   Numeric d = 0.0)
      : vec4{a, b, c, d} {}

  constexpr stokes(lazy_stokes auto &&a)
    requires(std::is_rvalue_reference_v<decltype(a)>)
      : vec4(std::move(a)) {}
};

//! Addition of two stokes vectors
constexpr auto operator+(const stokes &a, const stokes &b) {
  return constexpr_vec_data{a} + constexpr_vec_data{b};
}

//! Addition of a stokes vector and a lazy one
constexpr auto operator+(const stokes &a, const lazy_stokes auto &b) {
  return constexpr_vec_data{a} + b;
}

//! Addition of a stokes vector and a lazy one
constexpr auto operator+(const lazy_stokes auto &a, const stokes &b) {
  return a + constexpr_vec_data{b};
}

//! Subtraction between two stokes vectors
constexpr auto operator-(const stokes &a, const stokes &b) {
  return constexpr_vec_data{a} - constexpr_vec_data{b};
}

//! Subtraction between a stokes vector and a lazy one
constexpr auto operator-(const stokes &a, const lazy_stokes auto &b) {
  return constexpr_vec_data{a} - b;
}

//! Subtraction between a stokes vector and a lazy one
constexpr auto operator-(const lazy_stokes auto &a, const stokes &b) {
  return a - constexpr_vec_data{b};
}

//! Scaling a stokes vector
constexpr auto operator*(const Numeric &a, const stokes &b) {
  return vscl{a, constexpr_vec_data{b}};
}

//! Scaling a stokes vector
constexpr auto operator*(const stokes &a, const Numeric &b) {
  return vscl{b, constexpr_vec_data{a}};
}

//! Take the average of two stokes vectors
constexpr auto avg(const stokes &a, const stokes &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a stokes vector and a non-stokes vector
constexpr auto avg(const stokes &a, const lazy_stokes auto &b) {
  return 0.5 * a + 0.5 * b;
}

//! Take the average of a stokes vector and a non-stokes vector
constexpr auto avg(const lazy_stokes auto &a, const stokes &b) {
  return 0.5 * a + 0.5 * b;
}

using stokes_view = matpack::matpack_view<stokes, 1, false, false>;
using const_stokes_view = matpack::matpack_view<stokes, 1, true, false>;
} // namespace rtepack

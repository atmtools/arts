#pragma once

#include "matpack_mdspan.h"

#pragma GCC diagnostic push
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#else
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
#include <Eigen/Dense>
#pragma GCC diagnostic pop

namespace matpack::eigen {
template <md::any_exhaustive<2> U> auto mat(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic,
                                  Eigen::Dynamic, Eigen::RowMajor>>(
      x.data_handle(), x.extent(0), x.extent(1));
}

template <md::any_strided<2> U> auto mat(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic,
                                  Eigen::Dynamic, Eigen::RowMajor>,
                    Eigen::Unaligned,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
      x.unsafe_data_handle(), x.extent(0), x.extent(1),
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{x.stride(0), x.stride(1)});
}

template <md::any_exhaustive<1> U> auto col_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, 1, Eigen::Dynamic>>(
      x.data_handle(), 1, x.extent(0));
}

template <md::any_exhaustive<1> U> auto row_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic, 1>>(
      x.data_handle(), x.extent(0), 1);
}

template <md::any_strided<1> U> auto col_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, 1, Eigen::Dynamic>,
                    Eigen::Unaligned,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
      x.unsafe_data_handle(), 1, x.extent(0),
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
}

template <md::any_strided<1> U> auto row_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic, 1>,
                    Eigen::Unaligned,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
     x.unsafe_data_handle(), x.extent(0), 1,
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
}

template <md::matpack_matrix U> auto as_eigen(const U &x) { return mat(x); }

template <md::matpack_vector U> auto as_eigen(const U &x) { return row_vec(x); }
} // namespace eigen::matpack

namespace matpack::md {
template <matpack_matrix MAT, matpack_vector VEC>
auto operator*(const MAT &A, const VEC&x) {
  return eigen::mat(A) * eigen::row_vec(x);
}

template <matpack_matrix MAT1, matpack_matrix MAT2>
auto operator*(const MAT1 &A, const MAT2 &B) {
  return eigen::mat(A) * eigen::mat(B);
}

template <arithmetic CONST, matpack_matrix_or_vector MATVEC>
auto operator*(const CONST& a, const MATVEC &b) {
  return a * eigen::as_eigen(b);
}
}

auto operator-(matpack::md::matpack_matrix_or_vector auto &&x, matpack::md::matpack_matrix_or_vector auto &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) -
         matpack::eigen::as_eigen(std::forward<decltype(y)>(y));
}

template <typename Derived>
auto operator*(Eigen::MatrixBase<Derived> &&A, matpack::md::matpack_matrix_or_vector auto &&x) {
  return std::forward<decltype(A)>(A) *
         matpack::eigen::as_eigen(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator+(Eigen::MatrixBase<Derived> &&A, matpack::md::matpack_matrix_or_vector auto &&x) {
  return std::forward<decltype(A)>(A) +
         matpack::eigen::as_eigen(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator-(Eigen::MatrixBase<Derived> &&A, matpack::md::matpack_matrix_or_vector auto &&x) {
  return std::forward<decltype(A)>(A) -
         matpack::eigen::as_eigen(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator*(matpack::md::matpack_matrix_or_vector auto &&A, Eigen::MatrixBase<Derived> &&x) {
  return matpack::eigen::as_eigen(std::forward<decltype(A)>(A)) *
         std::forward<decltype(x)>(x);
}

template <typename Derived>
auto operator+(matpack::md::matpack_matrix_or_vector auto &&x, Eigen::MatrixBase<Derived> &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) +
         std::forward<decltype(y)>(y);
}

template <typename Derived>
auto operator-(matpack::md::matpack_matrix_or_vector auto &&x, Eigen::MatrixBase<Derived> &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) -
         std::forward<decltype(y)>(y);
}

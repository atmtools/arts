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

namespace matpack {
namespace eigen {
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
} // namespace eigen

auto operator*(md::matpack_matrix auto &&A, md::matpack_vector auto &&x) {
  return matpack::eigen::mat(std::forward<decltype(A)>(A)) *
         matpack::eigen::row_vec(std::forward<decltype(x)>(x));
}

auto operator*(md::matpack_matrix auto &&A, md::matpack_matrix auto &&B) {
  return matpack::eigen::mat(std::forward<decltype(A)>(A)) *
         matpack::eigen::mat(std::forward<decltype(B)>(B));
}

auto operator*(md::arithmetic auto &&a, md::matpack_matrix_or_vector auto &&b) {
  return std::forward<decltype(a)>(a) *
         matpack::eigen::as_eigen(std::forward<decltype(b)>(b));
}

auto operator+(md::matpack_matrix_or_vector auto &&x, md::matpack_matrix_or_vector auto &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) +
         matpack::eigen::as_eigen(std::forward<decltype(y)>(y));
}

auto operator-(md::matpack_matrix_or_vector auto &&x, md::matpack_matrix_or_vector auto &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) -
         matpack::eigen::as_eigen(std::forward<decltype(y)>(y));
}

template <typename Derived>
auto operator*(Eigen::MatrixBase<Derived> &&A, md::matpack_matrix_or_vector auto &&x) {
  return std::forward<decltype(A)>(A) *
         matpack::eigen::as_eigen(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator+(Eigen::MatrixBase<Derived> &&A, md::matpack_matrix_or_vector auto &&x) {
  return std::forward<decltype(A)>(A) +
         matpack::eigen::as_eigen(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator-(Eigen::MatrixBase<Derived> &&A, md::matpack_matrix_or_vector auto &&x) {
  return std::forward<decltype(A)>(A) -
         matpack::eigen::as_eigen(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator*(md::matpack_matrix_or_vector auto &&A, Eigen::MatrixBase<Derived> &&x) {
  return matpack::eigen::as_eigen(std::forward<decltype(A)>(A)) *
         std::forward<decltype(x)>(x);
}

template <typename Derived>
auto operator+(md::matpack_matrix_or_vector auto &&x, Eigen::MatrixBase<Derived> &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) +
         std::forward<decltype(y)>(y);
}

template <typename Derived>
auto operator-(md::matpack_matrix_or_vector auto &&x, Eigen::MatrixBase<Derived> &&y) {
  return matpack::eigen::as_eigen(std::forward<decltype(x)>(x)) -
         std::forward<decltype(y)>(y);
}
} // namespace matpack

#pragma once

#include <type_traits>

#include "matpackI.h"
#include "matpack_complex.h"
#include "matpack_concepts.h"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-user-provided-dtor"
#else
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif

#include <Eigen/Dense>

#pragma GCC diagnostic pop

namespace matpack::eigen {

//! Setup for eigen mapping
template <typename internal_type, bool constant_type = true>
struct eigen_type {
  using stride_type = Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>;
  using matrix_type = Eigen::
      Matrix<internal_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using constant_matrix_type = std::add_const_t<matrix_type>;
  using mutable_matrix_type = std::remove_const_t<matrix_type>;
  using map = Eigen::Map<std::conditional_t<constant_type,
                                            constant_matrix_type,
                                            mutable_matrix_type>,
                         Eigen::Unaligned,
                         stride_type>;
};

//! An actual map to an eigen type
template <typename internal_type, bool constant_type = true>
using eigen_map = typename eigen_type<internal_type, constant_type>::map;

//! The stride of the eigen type
template <typename internal_type, bool constant_type = true>
using eigen_stride =
    typename eigen_type<internal_type, constant_type>::stride_type;

//! Map the input to a non-owning const-correct Eigen Map representing a column vector
auto col_vec(matpack::vector auto&& x) {
  using internal_type =
      std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;
  constexpr bool real_type = std::is_same_v<internal_type, Numeric>;
  using matrix_type_real =
      std::conditional_t<real_type, Eigen::RowVectorXd, Eigen::RowVectorXcd>;
  using matrix_type = std::
      conditional_t<constant_type, const matrix_type_real, matrix_type_real>;
  using stride_type = Eigen::InnerStride<Eigen::Dynamic>;
  using map_type = Eigen::Map<matrix_type, Eigen::Unaligned, stride_type>;

  static_assert(real_type or std::is_same_v<internal_type, Complex>,
                "Only for Complex and Numeric");

  return map_type(
      x.get_c_array() + x.selem(), x.nelem(), stride_type{x.delem()});
}

//! Map the input to a non-owning const-correct Eigen Map representing a column vector
auto row_vec(matpack::vector auto&& x) {
  using internal_type =
      std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;
  constexpr bool real_type = std::is_same_v<internal_type, Numeric>;
  using matrix_type_real =
      std::conditional_t<real_type, Eigen::VectorXd, Eigen::VectorXcd>;
  using matrix_type = std::
      conditional_t<constant_type, const matrix_type_real, matrix_type_real>;
  using stride_type = Eigen::InnerStride<Eigen::Dynamic>;
  using map_type = Eigen::Map<matrix_type, Eigen::Unaligned, stride_type>;

  static_assert(real_type or std::is_same_v<internal_type, Complex>,
                "Only for Complex and Numeric");

  return map_type(
      x.get_c_array() + x.selem(), x.nelem(), stride_type{x.delem()});
}

//! Map the input to a non-owning const-correct Eigen Map representing a row matrix
auto mat(matpack::matrix auto&& x) {
  using internal_type =
      std::remove_cvref_t<std::remove_pointer_t<decltype(x.get_c_array())>>;
  constexpr bool constant_type = std::is_const_v<decltype(x)>;

  using stride_type = eigen_stride<internal_type, constant_type>;
  using matrix_map = eigen_map<internal_type, constant_type>;

  return matrix_map(x.get_c_array() + x.selem(),
                    x.nrows(),
                    x.ncols(),
                    stride_type(x.drows(), x.dcols()));
}

//! Test if the type represents a standard layout std::vector or std::array
template <typename T>
concept standard_vector = requires(T a) {
  a.data();
  { a.size() } -> std::integral;
  { a[0] } -> complex_or_real;
};

//! Map the input to a non-owning const Eigen Map representing a column vector
auto col_vec(standard_vector auto&& x) {
  using internal_type =
      std::remove_cvref_t<std::remove_pointer_t<decltype(x.data())>>;

  using stride_type = eigen_stride<internal_type>;
  using matrix_map = eigen_map<internal_type>;

  return matrix_map(x.data(), 1, x.size(), stride_type(1, 1));
}

//! Map the input to a non-owning const Eigen Map representing a row vector
auto row_vec(standard_vector auto&& x) {
  using internal_type =
      std::remove_cvref_t<std::remove_pointer_t<decltype(x.data())>>;

  using stride_type = eigen_stride<internal_type>;
  using matrix_map = eigen_map<internal_type>;

  return matrix_map(x.data(), x.size(), 1, stride_type(1, 1));
}

//! A generic concept we might want to move out of here
template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;
}  // namespace matpack::eigen

auto operator*(matpack::matrix auto&& A, matpack::vector auto&& x) {
  return matpack::eigen::mat(std::forward<decltype(A)>(A)) *
         matpack::eigen::row_vec(std::forward<decltype(x)>(x));
}

auto operator*(matpack::matrix auto&& A, matpack::matrix auto&& B) {
  return matpack::eigen::mat(std::forward<decltype(A)>(A)) *
         matpack::eigen::mat(std::forward<decltype(B)>(B));
}

template <typename Derived>
auto operator*(Eigen::MatrixBase<Derived>&& A, matpack::vector auto&& x) {
  return std::forward<decltype(A)>(A) *
         matpack::eigen::row_vec(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator*(Eigen::MatrixBase<Derived>&& A, matpack::matrix auto&& B) {
  return std::forward<decltype(A)>(A) *
         matpack::eigen::mat(std::forward<decltype(B)>(B));
}

template <typename Derived>
auto operator*(matpack::matrix auto&& A, Eigen::MatrixBase<Derived>&& x) {
  return matpack::eigen::mat(std::forward<decltype(A)>(A)) *
         std::forward<decltype(x)>(x);
}

auto operator+(matpack::vector auto&& x, matpack::vector auto&& y) {
  return matpack::eigen::row_vec(std::forward<decltype(x)>(x)) +
         matpack::eigen::row_vec(std::forward<decltype(y)>(y));
}

auto operator+(matpack::matrix auto&& x, matpack::matrix auto&& y) {
  return matpack::eigen::row_vec(std::forward<decltype(x)>(x)) +
         matpack::eigen::row_vec(std::forward<decltype(y)>(y));
}

template <typename Derived>
auto operator+(Eigen::MatrixBase<Derived>&& A, matpack::vector auto&& x) {
  return std::forward<decltype(A)>(A) +
         matpack::eigen::row_vec(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator+(Eigen::MatrixBase<Derived>&& A, matpack::matrix auto&& B) {
  return std::forward<decltype(A)>(A) +
         matpack::eigen::mat(std::forward<decltype(B)>(B));
}

template <typename Derived>
auto operator+(matpack::vector auto&& x, Eigen::MatrixBase<Derived>&& y) {
  return matpack::eigen::row_vec(std::forward<decltype(x)>(x)) +
         std::forward<decltype(y)>(y);
}

template <typename Derived>
auto operator+(matpack::matrix auto&& x, Eigen::MatrixBase<Derived>&& y) {
  return matpack::eigen::mat(std::forward<decltype(x)>(x)) +
         std::forward<decltype(y)>(y);
}

auto operator-(matpack::vector auto&& x, matpack::vector auto&& y) {
  return matpack::eigen::row_vec(std::forward<decltype(x)>(x)) -
         matpack::eigen::row_vec(std::forward<decltype(y)>(y));
}

auto operator-(matpack::matrix auto&& x, matpack::matrix auto&& y) {
  return matpack::eigen::row_vec(std::forward<decltype(x)>(x)) -
         matpack::eigen::row_vec(std::forward<decltype(y)>(y));
}

template <typename Derived>
auto operator-(Eigen::MatrixBase<Derived>&& A, matpack::vector auto&& x) {
  return std::forward<decltype(A)>(A) -
         matpack::eigen::row_vec(std::forward<decltype(x)>(x));
}

template <typename Derived>
auto operator-(Eigen::MatrixBase<Derived>&& A, matpack::matrix auto&& B) {
  return std::forward<decltype(A)>(A) -
         matpack::eigen::mat(std::forward<decltype(B)>(B));
}

template <typename Derived>
auto operator-(matpack::vector auto&& x, Eigen::MatrixBase<Derived>&& y) {
  return matpack::eigen::row_vec(std::forward<decltype(x)>(x)) -
         std::forward<decltype(y)>(y);
}

auto operator*(matpack::eigen::arithmetic auto&& a, matpack::vector auto&& b) {
  return std::forward<decltype(a)>(a) *
         matpack::eigen::row_vec(std::forward<decltype(b)>(b));
}

auto operator*(matpack::eigen::arithmetic auto&& a, matpack::matrix auto&& b) {
  return std::forward<decltype(a)>(a) *
         matpack::eigen::mat(std::forward<decltype(b)>(b));
}

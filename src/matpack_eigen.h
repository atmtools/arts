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

//! Can perform basic vector operations
template <typename T, typename U>
concept vector_linear_algebraic =
    not vector<T> and not matrix<T> and vector<U> and requires(T a, U b) {
  a * row_vec(b);
  a + row_vec(b);
  a - row_vec(b);
  row_vec(b) + a;
  row_vec(b) - a;
};

//! Can perform basic matrix operations
template <typename T, typename U>
concept matrix_linear_algebraic =
    not vector<T> and not matrix<T> and matrix<U> and requires(T a, U b) {
  a * mat(b);
  a + mat(b);
  a - mat(b);
  mat(b) + a;
  mat(b) - a;
};

/*! Concept that the second type is a matpack type and the former an eigen type,
    or at least something that supports the intended operations */
template <typename T, typename U>
concept linear_algebraic =
    matrix_linear_algebraic<T, U> or vector_linear_algebraic<T, U>;
}  // namespace matpack::eigen

auto operator*(matpack::matrix auto&& A, matpack::vector auto&& x) {
  using namespace matpack::eigen;
  return mat(std::forward<decltype(A)>(A)) *
         row_vec(std::forward<decltype(x)>(x));
}

auto operator*(matpack::matrix auto&& A, matpack::matrix auto&& B) {
  using namespace matpack::eigen;
  return mat(std::forward<decltype(A)>(A)) * mat(std::forward<decltype(B)>(B));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Vector>
auto operator*(T&& A, matpack::vector auto&& x) {
  using namespace matpack::eigen;
  return std::forward<T>(A) * row_vec(std::forward<decltype(x)>(x));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Matrix>
auto operator*(T&& A, matpack::matrix auto&& B) {
  using namespace matpack::eigen;
  return std::forward<T>(A) * mat(std::forward<decltype(B)>(B));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Matrix>
auto operator*(matpack::matrix auto&& A, T&& x) {
  using namespace matpack::eigen;
  return mat(std::forward<decltype(A)>(A)) * std::forward<T>(x);
}

auto operator+(matpack::vector auto&& x, matpack::vector auto&& y) {
  using namespace matpack::eigen;
  return row_vec(std::forward<decltype(x)>(x)) +
         row_vec(std::forward<decltype(y)>(y));
}

auto operator+(matpack::matrix auto&& x, matpack::matrix auto&& y) {
  using namespace matpack::eigen;
  return row_vec(std::forward<decltype(x)>(x)) +
         row_vec(std::forward<decltype(y)>(y));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Vector>
auto operator+(T&& A, matpack::vector auto&& x) {
  using namespace matpack::eigen;
  return std::forward<T>(A) + row_vec(std::forward<decltype(x)>(x));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Matrix>
auto operator+(T&& A, matpack::matrix auto&& B) {
  using namespace matpack::eigen;
  return std::forward<T>(A) + mat(std::forward<decltype(B)>(B));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Vector>
auto operator+(matpack::vector auto&& x, T&& y) {
  using namespace matpack::eigen;
  return row_vec(std::forward<decltype(x)>(x)) + std::forward<T>(y);
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Matrix>
auto operator+(matpack::matrix auto&& x, T&& y) {
  using namespace matpack::eigen;
  return mat(std::forward<decltype(x)>(x)) + std::forward<T>(y);
}

auto operator-(matpack::vector auto&& x, matpack::vector auto&& y) {
  using namespace matpack::eigen;
  return row_vec(std::forward<decltype(x)>(x)) -
         row_vec(std::forward<decltype(y)>(y));
}

auto operator-(matpack::matrix auto&& x, matpack::matrix auto&& y) {
  using namespace matpack::eigen;
  return row_vec(std::forward<decltype(x)>(x)) -
         row_vec(std::forward<decltype(y)>(y));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Vector>
auto operator-(T&& A, matpack::vector auto&& x) {
  using namespace matpack::eigen;
  return std::forward<T>(A) - row_vec(std::forward<decltype(x)>(x));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Matrix>
auto operator-(T&& A, matpack::matrix auto&& B) {
  using namespace matpack::eigen;
  return std::forward<T>(A) - mat(std::forward<decltype(B)>(B));
}

template <typename T>
requires matpack::eigen::linear_algebraic<T, Vector>
auto operator-(matpack::vector auto&& x, T&& y) {
  using namespace matpack::eigen;
  return row_vec(std::forward<decltype(x)>(x)) - std::forward<T>(y);
}

auto operator*(matpack::eigen::arithmetic auto&& a, matpack::vector auto&& b) {
  using namespace matpack::eigen;
  return std::forward<decltype(a)>(a) * row_vec(std::forward<decltype(b)>(b));
}

auto operator*(matpack::eigen::arithmetic auto&& a, matpack::matrix auto&& b) {
  using namespace matpack::eigen;
  return std::forward<decltype(a)>(a) * mat(std::forward<decltype(b)>(b));
}

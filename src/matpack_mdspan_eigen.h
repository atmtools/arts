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
/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U An exhaustive matpack matrix type
 * @param x The value
 * @return An appropriate eigen matrix map
 */
template <md::any_exhaustive<2> U> auto mat(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic,
                                  Eigen::Dynamic, Eigen::RowMajor>>(
      x.data_handle(), x.extent(0), x.extent(1));
}

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U A strided matpack matrix type
 * @param x The value
 * @return An appropriate eigen matrix map
 */
template <md::any_strided<2> U> auto mat(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic,
                                  Eigen::Dynamic, Eigen::RowMajor>,
                    Eigen::Unaligned,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
      x.unsafe_data_handle(), x.extent(0), x.extent(1),
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{x.stride(0), x.stride(1)});
}

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U An exhaustive matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as column vector
 */
template <md::any_exhaustive<1> U> auto col_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, 1, Eigen::Dynamic>>(
      x.data_handle(), 1, x.extent(0));
}

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U An exhaustive matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as row vector
 */
template <md::any_exhaustive<1> U> auto row_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic, 1>>(
      x.data_handle(), x.extent(0), 1);
}

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U A strided matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as column vector
 */
template <md::any_strided<1> U> auto col_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, 1, Eigen::Dynamic>,
                    Eigen::Unaligned,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
      x.unsafe_data_handle(), 1, x.extent(0),
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
}

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U A strided matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as row vector
 */
template <md::any_strided<1> U> auto row_vec(const U &x) {
  return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic, 1>,
                    Eigen::Unaligned,
                    Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
     x.unsafe_data_handle(), x.extent(0), 1,
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
}

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U A strided matpack matrix type
 * @param x The value
 * @return An appropriate eigen matrix map
 */
template <md::matpack_matrix U> auto as_eigen(const U &x) { return mat(x); }

/** Treat the input as an Eigen matrix type 
 * 
 * @tparam U A strided matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as row vector
 */
template <md::matpack_vector U> auto as_eigen(const U &x) { return row_vec(x); }
} // namespace matpack::eigen

namespace matpack::md {
//! Perform the math and return a lazy object; A * x
template <matpack_matrix MAT, matpack_vector VEC>
auto operator*(const MAT &A, const VEC&x) {
  return eigen::as_eigen(A) * eigen::as_eigen(x);
}

//! Perform the math and return a lazy object; A * B
template <matpack_matrix MAT1, matpack_matrix MAT2>
auto operator*(const MAT1 &A, const MAT2 &B) {
  return eigen::as_eigen(A) * eigen::as_eigen(B);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, matpack_matrix_or_vector MATVEC>
auto operator*(const CONST& a, const MATVEC &b) {
  return a * eigen::as_eigen(b);
}

/*  FIXME: THIS FAILS TO COMPILE FOR UNKNOWN REASONS
//! Perform the math and return a lazy object; x + y
auto operator+(matpack_matrix_or_vector auto &&x, matpack_matrix_or_vector auto &&y) {
  return eigen::as_eigen(std::forward<decltype(x)>(x)) +
         eigen::as_eigen(std::forward<decltype(y)>(y));
}
*/

//! Perform the math and return a lazy object; x - y
auto operator-(matpack_matrix_or_vector auto &&x, matpack_matrix_or_vector auto &&y) {
  return eigen::as_eigen(std::forward<decltype(x)>(x)) -
         eigen::as_eigen(std::forward<decltype(y)>(y));
}

//! Perform the math and return a lazy object; A * x
template <typename Derived, matpack_matrix_or_vector MATVEC>
auto operator*(Eigen::MatrixBase<Derived> &&A, const MATVEC &x) {
  return std::forward<Eigen::MatrixBase<Derived>>(A) * eigen::as_eigen(x);
}

//! Perform the math and return a lazy object; A * x
template <typename Derived, matpack_matrix_or_vector MATVEC>
auto operator*(const MATVEC &A, Eigen::MatrixBase<Derived> &&x) {
  return eigen::as_eigen(A) * std::forward<Eigen::MatrixBase<Derived>>(x);
}

//! Perform the math and return a lazy object; a + b
template <typename Derived, matpack_matrix_or_vector MATVEC>
auto operator+(Eigen::MatrixBase<Derived> &&a, const MATVEC &b) {
  return std::forward<Eigen::MatrixBase<Derived>>(a) + matpack::eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a + b
template <typename Derived, matpack_matrix_or_vector MATVEC>
auto operator+(const MATVEC &a, Eigen::MatrixBase<Derived> &&b) {
  return matpack::eigen::as_eigen(a) + std::forward<Eigen::MatrixBase<Derived>>(b);
}

//! Perform the math and return a lazy object; a - b
template <typename Derived, matpack_matrix_or_vector MATVEC>
auto operator-(Eigen::MatrixBase<Derived> &&a, const MATVEC &b) {
  return std::forward<Eigen::MatrixBase<Derived>>(a) - matpack::eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a - b
template <typename Derived, matpack_matrix_or_vector MATVEC>
auto operator-(const MATVEC &a, Eigen::MatrixBase<Derived> &&b) {
  return matpack::eigen::as_eigen(a) - std::forward<Eigen::MatrixBase<Derived>>(b);
}
}  // namespace matpack::md

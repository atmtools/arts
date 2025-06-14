#pragma once

#include "matpack_mdspan_common.h"
#ifndef _MSC_VER
#pragma GCC diagnostic push
#if defined(__clang__)
#pragma GCC diagnostic ignored "-Wdeprecated-copy-with-dtor"
#else
#pragma GCC diagnostic ignored "-Wclass-memaccess"
#endif
#endif

#include <Eigen/Dense>

#ifndef _MSC_VER
#pragma GCC diagnostic pop
#endif

#include "matpack_mdspan.h"

namespace matpack::eigen {
#define EIGEN_STRIDED_MAT(U, T)                                          \
  Eigen::Map<                                                            \
      Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, \
      Eigen::Unaligned,                                                  \
      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>                     \
      mat(T);
EIGEN_STRIDED_MAT(Numeric, StridedMatrixView &)
EIGEN_STRIDED_MAT(Numeric, StridedConstMatrixView &)
EIGEN_STRIDED_MAT(Numeric, const StridedMatrixView &)
EIGEN_STRIDED_MAT(Numeric, const StridedConstMatrixView &)
EIGEN_STRIDED_MAT(Complex, StridedComplexMatrixView &)
EIGEN_STRIDED_MAT(Complex, StridedConstComplexMatrixView &)
EIGEN_STRIDED_MAT(Complex, const StridedComplexMatrixView &)
EIGEN_STRIDED_MAT(Complex, const StridedConstComplexMatrixView &)
#undef EIGEN_STRIDED_MAT

#define EIGEN_MAT(U, T)                                                  \
  Eigen::Map<                                                            \
      Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> \
      mat(T);
EIGEN_MAT(Numeric, Matrix &)
EIGEN_MAT(Numeric, MatrixView &)
EIGEN_MAT(Numeric, ConstMatrixView &)
EIGEN_MAT(Numeric, const Matrix &)
EIGEN_MAT(Numeric, const MatrixView &)
EIGEN_MAT(Numeric, const ConstMatrixView &)
EIGEN_MAT(Complex, ComplexMatrix &)
EIGEN_MAT(Complex, ComplexMatrixView &)
EIGEN_MAT(Complex, ConstComplexMatrixView &)
EIGEN_MAT(Complex, const ComplexMatrix &)
EIGEN_MAT(Complex, const ComplexMatrixView &)
EIGEN_MAT(Complex, const ConstComplexMatrixView &)
#undef EIGEN_MAT

#define EIGEN_STRIDED_VEC(U, T)                             \
  Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>,           \
             Eigen::Unaligned,                              \
             Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> \
      col_vec(T);                                           \
  Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>,           \
             Eigen::Unaligned,                              \
             Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> \
      row_vec(T);
EIGEN_STRIDED_VEC(Numeric, StridedVectorView &)
EIGEN_STRIDED_VEC(Numeric, StridedConstVectorView &)
EIGEN_STRIDED_VEC(Numeric, const StridedVectorView &)
EIGEN_STRIDED_VEC(Numeric, const StridedConstVectorView &)
EIGEN_STRIDED_VEC(Complex, StridedComplexVectorView &)
EIGEN_STRIDED_VEC(Complex, StridedConstComplexVectorView &)
EIGEN_STRIDED_VEC(Complex, const StridedComplexVectorView &)
EIGEN_STRIDED_VEC(Complex, const StridedConstComplexVectorView &)
#undef EIGEN_STRIDED_VEC

#define EIGEN_VEC(U, T)                                       \
  Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>> col_vec(T); \
  Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>> row_vec(T);
EIGEN_VEC(Numeric, Vector &)
EIGEN_VEC(Numeric, VectorView &)
EIGEN_VEC(Numeric, ConstVectorView &)
EIGEN_VEC(Numeric, const Vector &)
EIGEN_VEC(Numeric, const VectorView &)
EIGEN_VEC(Numeric, const ConstVectorView &)
EIGEN_VEC(Complex, ComplexVector &)
EIGEN_VEC(Complex, ComplexVectorView &)
EIGEN_VEC(Complex, ConstComplexVectorView &)
EIGEN_VEC(Complex, const ComplexVector &)
EIGEN_VEC(Complex, const ComplexVectorView &)
EIGEN_VEC(Complex, const ConstComplexVectorView &)
#undef EIGEN_VEC

template <dyn_ranked_md<2> T>
constexpr auto as_eigen(T &&x) {
  return mat(std::forward<T>(x));
}

template <dyn_ranked_md<1> T>
constexpr auto as_eigen(T &&x) {
  return row_vec(std::forward<T>(x));
}

template <typename T>
concept eigenable = dyn_ranked_md<T, 2> or dyn_ranked_md<T, 1>;
}  // namespace matpack::eigen

namespace matpack {
//! Perform the math and return a lazy object; A * x
template <dyn_ranked_md<2> MAT, dyn_ranked_md<1> VEC>
constexpr auto operator*(const MAT &A, const VEC &x) {
  return eigen::as_eigen(A) * eigen::as_eigen(x);
}

//! Perform the math and return a lazy object; A * B
template <dyn_ranked_md<2> MATONE, dyn_ranked_md<2> MATTWO>
constexpr auto operator*(const MATONE &A, const MATTWO &B) {
  return eigen::as_eigen(A) * eigen::as_eigen(B);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, dyn_ranked_md<2> MAT>
constexpr auto operator*(const std::complex<CONST> &a, const MAT &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, dyn_ranked_md<1> VEC>
constexpr auto operator*(const std::complex<CONST> &a, const VEC &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, dyn_ranked_md<2> MAT>
constexpr auto operator*(const MAT &b, const std::complex<CONST> &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, dyn_ranked_md<1> VEC>
constexpr auto operator*(const VEC &b, const std::complex<CONST> &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, dyn_ranked_md<2> MAT>
constexpr auto operator*(const CONST &a, const MAT &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, dyn_ranked_md<1> VEC>
constexpr auto operator*(const CONST &a, const VEC &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, dyn_ranked_md<2> MAT>
constexpr auto operator*(const MAT &b, const CONST &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, dyn_ranked_md<1> VEC>
constexpr auto operator*(const VEC &b, const CONST &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; x + y
template <dyn_ranked_md<2> ONE, dyn_ranked_md<2> TWO>
constexpr auto operator+(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <dyn_ranked_md<1> ONE, dyn_ranked_md<1> TWO>
constexpr auto operator+(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x = y
template <dyn_ranked_md<2> ONE, dyn_ranked_md<2> TWO>
constexpr auto operator-(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <dyn_ranked_md<1> ONE, dyn_ranked_md<1> TWO>
constexpr auto operator-(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x * y
template <typename Derived, dyn_ranked_md<1> VEC>
constexpr auto operator*(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x * eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x * y
template <typename Derived, dyn_ranked_md<2> MAT>
constexpr auto operator*(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x * eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <typename Derived, dyn_ranked_md<1> VEC>
constexpr auto operator+(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <typename Derived, dyn_ranked_md<2> MAT>
constexpr auto operator+(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <typename Derived, dyn_ranked_md<1> VEC>
constexpr auto operator-(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <typename Derived, dyn_ranked_md<2> MAT>
constexpr auto operator-(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; y * x
template <typename Derived, dyn_ranked_md<1> VEC>
constexpr auto operator*(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) * x;
}

//! Perform the math and return a lazy object; y * x
template <typename Derived, dyn_ranked_md<2> MAT>
constexpr auto operator*(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) * x;
}

//! Perform the math and return a lazy object; y + x
template <typename Derived, dyn_ranked_md<1> VEC>
constexpr auto operator+(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) + x;
}

//! Perform the math and return a lazy object; y + x
template <typename Derived, dyn_ranked_md<2> MAT>
constexpr auto operator+(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) + x;
}

//! Perform the math and return a lazy object; y - x
template <typename Derived, dyn_ranked_md<1> VEC>
constexpr auto operator-(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) - x;
}

//! Perform the math and return a lazy object; y - x
template <typename Derived, dyn_ranked_md<2> MAT>
constexpr auto operator-(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) - x;
}
}  // namespace matpack

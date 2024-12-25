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
template <ranked_md<2> T>
constexpr auto mat(T &&x) {
  using U = value_type<T>;
  if constexpr (any_strided_view<T>) {
    return Eigen::Map<
        Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>,
        Eigen::Unaligned,
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        const_cast<U *>(x.data_handle()),
        x.extent(0),
        x.extent(1),
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{x.stride(0),
                                                      x.stride(1)});
  } else {
    return Eigen::Map<
        Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        const_cast<U *>(x.data_handle()), x.extent(0), x.extent(1));
  }
}

template <ranked_md<1> T>
constexpr auto col_vec(T &&x) {
  using U = value_type<T>;
  if constexpr (any_strided_view<T>) {
    return Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>,
                      Eigen::Unaligned,
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        const_cast<U *>(x.data_handle()),
        1,
        x.extent(0),
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
  } else {
    return Eigen::Map<Eigen::Matrix<U, 1, Eigen::Dynamic>>(
        const_cast<U *>(x.data_handle()), 1, x.extent(0));
  }
}

template <ranked_md<1> T>
constexpr auto row_vec(T &&x) {
  using U = value_type<T>;
  if constexpr (any_strided_view<T>) {
    return Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>,
                      Eigen::Unaligned,
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        const_cast<U *>(x.data_handle()),
        x.extent(0),
        1,
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
  } else {
    return Eigen::Map<Eigen::Matrix<U, Eigen::Dynamic, 1>>(
        const_cast<U *>(x.data_handle()), x.extent(0), 1);
  }
}

template <ranked_md<2> T>
constexpr auto as_eigen(T &&x) {
  return mat(std::forward<T>(x));
}

template <ranked_md<1> T>
constexpr auto as_eigen(T &&x) {
  return row_vec(std::forward<T>(x));
}

template <typename T>
concept eigenable = ranked_md<T, 2> or ranked_md<T, 1>;
}  // namespace matpack::eigen

namespace matpack {
//! Perform the math and return a lazy object; A * x
template <ranked_md<2> MAT, ranked_md<1> VEC>
constexpr auto operator*(const MAT &A, const VEC &x) {
  return eigen::as_eigen(A) * eigen::as_eigen(x);
}

//! Perform the math and return a lazy object; A * B
template <ranked_md<2> MATONE, ranked_md<2> MATTWO>
constexpr auto operator*(const MATONE &A, const MATTWO &B) {
  return eigen::as_eigen(A) * eigen::as_eigen(B);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, ranked_md<2> MAT>
constexpr auto operator*(const std::complex<CONST> &a, const MAT &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, ranked_md<1> VEC>
constexpr auto operator*(const std::complex<CONST> &a, const VEC &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, ranked_md<2> MAT>
constexpr auto operator*(const MAT &b, const std::complex<CONST> &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, ranked_md<1> VEC>
constexpr auto operator*(const VEC &b, const std::complex<CONST> &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, ranked_md<2> MAT>
constexpr auto operator*(const CONST &a, const MAT &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, ranked_md<1> VEC>
constexpr auto operator*(const CONST &a, const VEC &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, ranked_md<2> MAT>
constexpr auto operator*(const MAT &b, const CONST &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, ranked_md<1> VEC>
constexpr auto operator*(const VEC &b, const CONST &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; x + y
template <ranked_md<2> ONE, ranked_md<2> TWO>
constexpr auto operator+(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <ranked_md<1> ONE, ranked_md<1> TWO>
constexpr auto operator+(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x = y
template <ranked_md<2> ONE, ranked_md<2> TWO>
constexpr auto operator-(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <ranked_md<1> ONE, ranked_md<1> TWO>
constexpr auto operator-(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x * y
template <typename Derived, ranked_md<1> VEC>
constexpr auto operator*(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x * eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x * y
template <typename Derived, ranked_md<2> MAT>
constexpr auto operator*(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x * eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <typename Derived, ranked_md<1> VEC>
constexpr auto operator+(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <typename Derived, ranked_md<2> MAT>
constexpr auto operator+(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <typename Derived, ranked_md<1> VEC>
constexpr auto operator-(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <typename Derived, ranked_md<2> MAT>
constexpr auto operator-(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; y * x
template <typename Derived, ranked_md<1> VEC>
constexpr auto operator*(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) * x;
}

//! Perform the math and return a lazy object; y * x
template <typename Derived, ranked_md<2> MAT>
constexpr auto operator*(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) * x;
}

//! Perform the math and return a lazy object; y + x
template <typename Derived, ranked_md<1> VEC>
constexpr auto operator+(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) + x;
}

//! Perform the math and return a lazy object; y + x
template <typename Derived, ranked_md<2> MAT>
constexpr auto operator+(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) + x;
}

//! Perform the math and return a lazy object; y - x
template <typename Derived, ranked_md<1> VEC>
constexpr auto operator-(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) - x;
}

//! Perform the math and return a lazy object; y - x
template <typename Derived, ranked_md<2> MAT>
constexpr auto operator-(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) - x;
}
}  // namespace matpack

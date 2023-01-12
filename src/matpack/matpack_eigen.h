#pragma once

#include "matpack_concepts.h"
#include "matpack_view.h"

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
 * @tparam U A matpack matrix view type
 * @param x The value
 * @return An appropriate eigen matrix map
 */
template <strict_sized_matpack_type<2> U> constexpr auto mat(const U &x) {
  if constexpr (is_always_exhaustive_v<U>) {
    return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic,
                                    Eigen::Dynamic, Eigen::RowMajor>>(
        x.data_handle(), x.extent(0), x.extent(1));
  } else {
    return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic,
                                    Eigen::Dynamic, Eigen::RowMajor>,
                      Eigen::Unaligned,
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        x.unsafe_data_handle(), x.extent(0), x.extent(1),
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{x.stride(0),
                                                      x.stride(1)});
  }
}

/** Treat the input as an Eigen matrix type
 *
 * @tparam U An matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as column vector
 */
template <strict_sized_matpack_type<1> U> constexpr auto col_vec(const U &x) {
  if constexpr (is_always_exhaustive_v<U>) {
    return Eigen::Map<Eigen::Matrix<typename U::value_type, 1, Eigen::Dynamic>>(
        x.data_handle(), 1, x.extent(0));
  } else {
    return Eigen::Map<Eigen::Matrix<typename U::value_type, 1, Eigen::Dynamic>,
                      Eigen::Unaligned,
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        x.unsafe_data_handle(), 1, x.extent(0),
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
  }
}

/** Treat the input as an Eigen matrix type
 *
 * @tparam U A matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as row vector
 */
template <strict_sized_matpack_type<1> U> constexpr auto row_vec(const U &x) {
  if constexpr (is_always_exhaustive_v<U>) {
    return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic, 1>>(
        x.data_handle(), x.extent(0), 1);
  } else {
    return Eigen::Map<Eigen::Matrix<typename U::value_type, Eigen::Dynamic, 1>,
                      Eigen::Unaligned,
                      Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>>(
        x.unsafe_data_handle(), x.extent(0), 1,
        Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>{1, x.stride(0)});
  }
}

/** Treat the input as an Eigen matrix type
 *
 * @tparam U A matpack matrix type
 * @param x The value
 * @return An appropriate eigen matrix map as matrix
 */
template <strict_sized_matpack_type<2> U> constexpr auto as_eigen(const U &x) {
  return mat(x);
}

/** Treat the input as an Eigen matrix type
 *
 * @tparam U A matpack vector type
 * @param x The value
 * @return An appropriate eigen matrix map as row vector
 */
template <strict_sized_matpack_type<1> U> constexpr auto as_eigen(const U &x) {
  return row_vec(x);
}
} // namespace matpack::eigen

namespace matpack {
//! Compute dot-product of x and y
template <strict_sized_matpack_type<1> VECONE,
          strict_sized_matpack_type<1> VECTWO>
constexpr auto operator*(const VECONE &x, const VECTWO &y) {
  return eigen::as_eigen(x).dot(eigen::as_eigen(y));
}

//! Perform the math and return a lazy object; A * x
template <strict_sized_matpack_type<2> MAT, strict_sized_matpack_type<1> VEC>
constexpr auto operator*(const MAT &A, const VEC &x) {
  return eigen::as_eigen(A) * eigen::as_eigen(x);
}

//! Perform the math and return a lazy object; A * B
template <strict_sized_matpack_type<2> MATONE,
          strict_sized_matpack_type<2> MATTWO>
constexpr auto operator*(const MATONE &A, const MATTWO &B) {
  return eigen::as_eigen(A) * eigen::as_eigen(B);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, strict_sized_matpack_type<2> MAT>
constexpr auto operator*(const CONST &a, const MAT &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; a * b
template <arithmetic CONST, strict_sized_matpack_type<1> VEC>
constexpr auto operator*(const CONST &a, const VEC &b) {
  return a * eigen::as_eigen(b);
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, strict_sized_matpack_type<2> MAT>
constexpr auto operator*(const MAT &b, const CONST &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; b * a
template <arithmetic CONST, strict_sized_matpack_type<1> VEC>
constexpr auto operator*(const VEC &b, const CONST &a) {
  return eigen::as_eigen(b) * a;
}

//! Perform the math and return a lazy object; x + y
template <strict_sized_matpack_type<2> ONE, strict_sized_matpack_type<2> TWO>
constexpr auto operator+(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <strict_sized_matpack_type<1> ONE, strict_sized_matpack_type<1> TWO>
constexpr auto operator+(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x = y
template <strict_sized_matpack_type<2> ONE, strict_sized_matpack_type<2> TWO>
constexpr auto operator-(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <strict_sized_matpack_type<1> ONE, strict_sized_matpack_type<1> TWO>
constexpr auto operator-(const ONE &x, const TWO &y) {
  return eigen::as_eigen(x) - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x * y
template <typename Derived, strict_sized_matpack_type<1> VEC>
constexpr auto operator*(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x * eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x * y
template <typename Derived, strict_sized_matpack_type<2> MAT>
constexpr auto operator*(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x * eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <typename Derived, strict_sized_matpack_type<1> VEC>
constexpr auto operator+(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x + y
template <typename Derived, strict_sized_matpack_type<2> MAT>
constexpr auto operator+(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x + eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <typename Derived, strict_sized_matpack_type<1> VEC>
constexpr auto operator-(const Eigen::MatrixBase<Derived> &x, const VEC &y) {
  return x - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; x - y
template <typename Derived, strict_sized_matpack_type<2> MAT>
constexpr auto operator-(const Eigen::MatrixBase<Derived> &x, const MAT &y) {
  return x - eigen::as_eigen(y);
}

//! Perform the math and return a lazy object; y * x
template <typename Derived, strict_sized_matpack_type<1> VEC>
constexpr auto operator*(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) * x;
}

//! Perform the math and return a lazy object; y * x
template <typename Derived, strict_sized_matpack_type<2> MAT>
constexpr auto operator*(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) * x;
}

//! Perform the math and return a lazy object; y + x
template <typename Derived, strict_sized_matpack_type<1> VEC>
constexpr auto operator+(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) + x;
}

//! Perform the math and return a lazy object; y + x
template <typename Derived, strict_sized_matpack_type<2> MAT>
constexpr auto operator+(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) + x;
}

//! Perform the math and return a lazy object; y - x
template <typename Derived, strict_sized_matpack_type<1> VEC>
constexpr auto operator-(const VEC &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) - x;
}

//! Perform the math and return a lazy object; y - x
template <typename Derived, strict_sized_matpack_type<2> MAT>
constexpr auto operator-(const MAT &y, const Eigen::MatrixBase<Derived> &x) {
  return eigen::as_eigen(y) - x;
}
} // namespace matpack

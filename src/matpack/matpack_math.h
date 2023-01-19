#pragma once

#include <limits>
#include <nonstd.h>

#include "debug.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"

#include <algorithm>
#include <numeric>

/** Return a transposed view of this matrix type
 * 
 * @param x Any matpack type that is 2-dimensional
 * @return A transpose view of that matpack type
 */
template <matpack::strict_size_matpack_type<2> MAT>
auto transpose(const MAT &x)
    -> matpack::matpack_view<matpack::matpack_value_type<MAT>, 2,
                             MAT::is_const(), true> {
  return matpack::strided_mdspan<matpack::matpack_value_type<MAT>, 2>{
      x.unsafe_data_handle(),
      {std::array{x.extent(1), x.extent(0)},
       std::array{x.stride(1), x.stride(0)}}};
}

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult_fast(ExhaustiveMatrixView A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstMatrixView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult_slow(MatrixView A, const ConstMatrixView& B, const ConstMatrixView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult_fast(ExhaustiveComplexMatrixView A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexMatrixView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult_slow(ComplexMatrixView A, const ConstComplexMatrixView& B, const ConstComplexMatrixView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any vector
 */
void mult_fast(ExhaustiveVectorView A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstVectorView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any vector
 */
void mult_slow(VectorView A, const ConstMatrixView& B, const ConstVectorView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any vector
 */
void mult_fast(ExhaustiveComplexVectorView A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexVectorView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any vector
 */
void mult_slow(ComplexVectorView A, const ConstComplexMatrixView& B, const ConstComplexVectorView& C);

/** Computes the 3-dim cross-product of B and C
 * 
 * @param A May not point at the same data as B or C
 * @param B Any vector
 * @param C Any vector
 */
void cross3(VectorView A, const ConstVectorView& B, const ConstVectorView& C);

/** Computes the 3-dim cross-product of B and C
 * 
 * @param A May not point at the same data as B or C
 * @param B Any vector
 * @param C Any vector
 */
void cross3(ComplexVectorView A, const ConstComplexVectorView& B, const ConstComplexVectorView& C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
template <matpack::strict_size_matpack_type<2> OUTMAT,
          matpack::strict_size_matpack_type<2> LMAT,
          matpack::strict_size_matpack_type<2> RMAT>
void mult(OUTMAT &&A, const LMAT &B, const RMAT &C)
  requires(matpack::same_value_type_v<OUTMAT, LMAT, RMAT>)
{
  if constexpr (matpack::is_always_exhaustive_v<OUTMAT> and
                matpack::is_always_exhaustive_v<LMAT> and
                matpack::is_always_exhaustive_v<RMAT>)
    mult_fast(std::forward<OUTMAT>(A), B, C);
  else
    mult_slow(std::forward<OUTMAT>(A), B, C);
}

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any vector
 */
template <matpack::strict_size_matpack_type<1> LHSVEC,
          matpack::strict_size_matpack_type<2> RHSMAT,
          matpack::strict_size_matpack_type<1> RHSVEC>
void mult(LHSVEC &&A, const RHSMAT &B, const RHSVEC &C)
  requires(matpack::same_value_type_v<LHSVEC, RHSMAT, RHSVEC>)
{
  if constexpr (matpack::is_always_exhaustive_v<LHSVEC> and
                matpack::is_always_exhaustive_v<RHSMAT> and
                matpack::is_always_exhaustive_v<RHSVEC>)
    mult_fast(std::forward<LHSVEC>(A), B, C);
  else
    mult_slow(std::forward<LHSVEC>(A), B, C);
}

//! Copies the diagonal of matrix view A and returns it
Vector diagonal(const ConstMatrixView& A);

/** Returns a uniform grid
 * 
 * @param x0 Initial value
 * @param N Number of elemets
 * @param dx Spacing between elements
 * @return Vector {x0, x0+dx, ... x0+(N-1)*dx}
 */
Vector uniform_grid(Numeric x0, Index N, Numeric dx);

/** Returns a uniform grid
 * 
 * @param x0 Initial value
 * @param N Number of elemets
 * @param dx Spacing between elements
 * @return Vector {x0, x0+dx, ... x0+(N-1)*dx}
 */
ComplexVector uniform_grid(Complex x0, Index N, Complex dx);

/** Sets out = f(in) for each element in in and out
 *
 * The two matpack types must have the same initial size
 *
 * @param[out] out The output, may be the same as in
 * @param[in] f The functional transformation
 * @param[in] in The input
 */
template <matpack::any_matpack_type OUT, matpack::any_matpack_type IN>
void transform(OUT &&out,
               matpack::matpack_value_type<OUT> (&f)(
                   matpack::matpack_value_type<IN>),
               const IN &in)
  requires(not std::remove_cvref_t<OUT>::is_const())
{
  ARTS_ASSERT(out.size() == in.size())
  std::transform(in.elem_begin(), in.elem_end(), out.elem_begin(),
                 f);
}

/** Returns the minimal value of the object
 * 
 * @param in Any matpack type
 * @return The minimal value
 */
template <matpack::any_matpack_type IN>
constexpr auto min(const IN &in) {
  ARTS_ASSERT(in.size() > 0)
  using T = matpack::matpack_value_type<IN>;
  //FIXME: CLANG bug means we want to use std::reduce rather than min_element
  return std::reduce(in.elem_begin(), in.elem_end(), std::numeric_limits<T>::max(), [](auto a, auto b){return a < b ? a : b;});
}

/** Returns the maximum value of the object
 * 
 * @param in Any matpack type
 * @return The maximum value
 */
template <matpack::any_matpack_type IN>
constexpr auto max(const IN &in) {
  ARTS_ASSERT(in.size() > 0)
  using T = matpack::matpack_value_type<IN>;
  //FIXME: CLANG bug means we want to use std::reduce rather than max_element
  return std::reduce(in.elem_begin(), in.elem_end(), std::numeric_limits<T>::lowest(), [](auto a, auto b){return a > b ? a : b;});
}

/** Returns the minimum and maximum value of the object
 * 
 * @param in Any matpack type
 * @return [min(in), max(in)]
 */
template <matpack::any_matpack_type IN>
constexpr auto minmax(const IN &in)
    -> std::pair<matpack::matpack_value_type<IN>,
                 matpack::matpack_value_type<IN>> {
  // FIXME: CLANG bug means we cant use minmax_element so we need to be less efficient
  return {min(in), max(in)};
}

/** Returns the sum of all elements of this object
 * 
 * @param in Any matpack type
 * @return The sum
 */
template <matpack::any_matpack_type IN>
constexpr auto sum(const IN &in) {
  return std::reduce(in.elem_begin(), in.elem_end());
}

/** Returns the mean of all elements of this object
 * 
 * @param in Any matpack type
 * @return The mean
 */
template <matpack::any_matpack_type IN> constexpr auto mean(const IN &in) {
  return sum(in) / static_cast<decltype(sum(in))>(in.size());
}

/** Returns the mean of all elements of this object except for NaN
 * 
 * @param in Any matpack type
 * @return The mean
 */
template <matpack::any_matpack_type IN> constexpr auto nanmean(const IN &in) {
  using T = matpack::matpack_value_type<IN>;
  using pt = std::pair<Index, T>;
  const auto [count, sum] = std::transform_reduce(
      in.elem_begin(), in.elem_end(), pt(0, 0),
      [](pt a, pt b) {
        return pt{a.first + b.first, a.second + b.second};
      },
      [](T a) {
        using namespace std;
        using namespace nonstd;
        return isnan(a) ? pt(0, 0) : pt(1, a);
      });
  return sum / static_cast<T>(count);
}

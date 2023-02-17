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
template <matpack::strict_rank_matpack_type<2> MAT>
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
void mult(MatrixView A, const ConstMatrixView &B, const ConstMatrixView &C);

/** Makes A = B * C
 * 
 * @param[out] A May not point at the same data as B or C
 * @param B Any matrix
 * @param C Any matrix
 */
void mult(ComplexMatrixView A, const ConstComplexMatrixView &B, const ConstComplexMatrixView &C);

/** Makes y = M * x
 * 
 * @param[out] y May not point at the same data as M or x
 * @param M Any matrix
 * @param x Any vector
 */
void mult(VectorView y, const ConstMatrixView &M, const ConstVectorView &x);

/** Makes y = M * x
 * 
 * @param[out] y May not point at the same data as M or x
 * @param M Any matrix
 * @param x Any vector
 */
void mult(ComplexVectorView y, const ConstComplexMatrixView &M, const ConstComplexVectorView &x);

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
        return nonstd::isnan(a) ? pt(0, 0) : pt(1, a);
      });
  return sum / static_cast<T>(count);
}

namespace matpack {
//! Compute dot-product of x and y
template <strict_rank_matpack_type<1> VECONE,
          strict_rank_matpack_type<1> VECTWO>
constexpr auto operator*(const VECONE &x, const VECTWO &y) {
  ARTS_ASSERT(x.size() == y.size(), x.size(), " vs ", y.size())
  using T = std::remove_cvref_t<decltype(x[0] * y[0])>;
  return std::transform_reduce(x.elem_begin(), x.elem_end(), y.elem_begin(), T{0});
}
}  // namespace matpack

/** Reverses the input in place using elementwise iteration
 *
 * Note that this will not allocate any new memory
 *
 * Calls std::reverse(inout.elem_begin(), inout.elem_end())
 *
 * @param[in] in Any mutable view of matpack data
 */
template <matpack::any_matpack_type INOUT>
constexpr void reverse_inplace(INOUT &&inout) {
  std::reverse(inout.elem_begin(), inout.elem_end());
}

/** Returns the reverse of the input
 *
 * Note that this will allocate new memory before calling reverse_inplace()
 * to reverse the allocated data
 *
 * @param[in] in Any input view of data
 * @return Allocated memory of the reverse of the input
 */
template <matpack::any_matpack_type IN>
constexpr matpack::matpack_data<matpack::matpack_value_type<IN>,
                                matpack::rank<IN>()>
reverse(const IN &in) {
  matpack::matpack_data<matpack::matpack_value_type<IN>, matpack::rank<IN>()>
      out(in);
  reverse_inplace(out);
  return out;
}

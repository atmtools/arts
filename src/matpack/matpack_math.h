#pragma once

#include <nonstd.h>

#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"

#include <algorithm>
#include <numeric>

template <matpack::strict_sized_matpack_type<2> MAT>
auto transpose(const MAT &x)
    -> matpack::matpack_view<matpack::matpack_value_type<MAT>, 2,
                             MAT::is_const(), true> {
  return matpack::strided_mdspan<Numeric, 2>{
      x.unsafe_data_handle(),
      {std::array{x.extent(1), x.extent(0)},
       std::array{x.stride(1), x.stride(0)}}};
}

void mult_fast(ExhaustiveMatrixView&& A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstMatrixView& C);

void mult_slow(MatrixView&& A, const ConstMatrixView& B, const ConstMatrixView& C);

void mult_fast(ExhaustiveComplexMatrixView&& A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexMatrixView& C);

void mult_slow(ComplexMatrixView&& A, const ConstComplexMatrixView& B, const ConstComplexMatrixView& C);

void mult_fast(ExhaustiveVectorView&& A, const ExhaustiveConstMatrixView& B, const ExhaustiveConstVectorView& C);

void mult_slow(VectorView&& A, const ConstMatrixView& B, const ConstVectorView& C);

void mult_fast(ExhaustiveComplexVectorView&& A, const ExhaustiveConstComplexMatrixView& B, const ExhaustiveConstComplexVectorView& C);

void mult_slow(ComplexVectorView&& A, const ConstComplexMatrixView& B, const ConstComplexVectorView& C);

template <matpack::strict_sized_matpack_type<2> OUTMAT,
          matpack::strict_sized_matpack_type<2> LMAT,
          matpack::strict_sized_matpack_type<2> RMAT>
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

template <matpack::strict_sized_matpack_type<1> LHSVEC,
          matpack::strict_sized_matpack_type<2> RHSMAT,
          matpack::strict_sized_matpack_type<1> RHSVEC>
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

Vector uniform_grid(Numeric x0, Index N, Numeric dx);

ComplexVector uniform_grid(Complex x0, Index N, Complex dx);

template <matpack::any_matpack_type OUT, matpack::any_matpack_type IN>
void transform(OUT &&out,
               matpack::matpack_value_type<OUT> (&my_func)(
                   matpack::matpack_value_type<IN>),
               const IN &in)
  requires(not std::remove_cvref_t<OUT>::is_const())
{
  std::transform(in.elem_begin(), in.elem_end(), out.elem_begin(),
                 my_func);
}

template <matpack::any_matpack_type IN>
constexpr auto min_element(const IN &in)  {
  return std::min_element(in.elem_begin(), in.elem_end());
}

template <matpack::any_matpack_type IN>
constexpr auto min(const IN &in) -> matpack::matpack_value_type<IN> {
  return *min_element(in);
}

template <matpack::any_matpack_type IN>
constexpr auto max_element(const IN &in)  {
  return std::max_element(in.elem_begin(), in.elem_end());
}

template <matpack::any_matpack_type IN>
constexpr auto max(const IN &in) -> matpack::matpack_value_type<IN> {
  return *max_element(in);
}

template <matpack::any_matpack_type IN>
constexpr auto minmax_element(const IN &in)  {
  return std::minmax_element(in.elem_begin(), in.elem_end());
}

template <matpack::any_matpack_type IN>
constexpr auto minmax(const IN &in) -> std::pair<matpack::matpack_value_type<IN>, matpack::matpack_value_type<IN>> {
  auto [mn, mx] = minmax_element(in);
  return {*mn, *mx};
}

template <matpack::any_matpack_type IN>
constexpr auto sum(const IN &in) {
  return std::reduce(in.elem_begin(), in.elem_end());
}

template <matpack::any_matpack_type IN>
constexpr auto has_nan(const IN &in) {
  return std::any_of(in.elem_begin(), in.elem_end(), [](auto& a){return nonstd::isnan(a);});
}

template <matpack::any_matpack_type IN>
constexpr auto is_decreasing(const IN &in) {
  return not has_nan(in) and in.elem_end() == std::adjacent_find(in.elem_begin(), in.elem_end(), [](auto& a, auto& b){return a <= b;});
}

template <matpack::any_matpack_type IN>
constexpr auto is_increasing(const IN &in) {
  return not has_nan(in) and in.elem_end() == std::adjacent_find(in.elem_begin(), in.elem_end(), [](auto& a, auto& b){return a >= b;});
}

template <matpack::any_matpack_type IN>
constexpr auto is_sorted(const IN &in) {
  return not has_nan(in) and std::is_sorted(in.elem_begin(), in.elem_end());
}

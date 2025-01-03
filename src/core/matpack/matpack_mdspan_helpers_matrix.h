#pragma once

#include "matpack_mdspan.h"

namespace matpack {
template <ranked_md<2> MAT>
constexpr auto transpose(MAT &&A) {
  using T       = element_type<MAT>;
  using rettype = matpack::strided_view_t<T, 2>;
  mapping_type<rettype> map{std::array{A.extent(1), A.extent(0)},
                            std::array{A.stride(1), A.stride(0)}};
  return rettype{mdstrided_t<element_type<MAT>, 2>{
      const_cast<T *>(std::forward<MAT>(A).data_handle()), map}};
}

template <mut_ranked_md<2> MAT>
constexpr MAT &inplace_transpose(MAT &x)
  requires(not MAT::is_const)
{
  assert(x.nrows() == x.ncols());
  for (Index i = 0; i < x.nrows(); ++i) {
    for (Index j = 0; j < i; ++j) {
      std::swap(x[i, j], x[j, i]);
    }
  }
  return x;
}

template <ranked_md<2> MAT>
constexpr strided_view_t<element_type<MAT>, 1> diagonal(MAT &&A) {
  assert(A.nrows() == A.ncols());
  return strided_view_t<element_type<MAT>, 1>{mdstrided_t<element_type<MAT>, 1>{
      A.data_handle(), {std::array{A.extent(0)}, std::array{A.extent(0) + 1}}}};
}
}  // namespace matpack

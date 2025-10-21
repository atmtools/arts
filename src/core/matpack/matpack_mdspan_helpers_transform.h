#pragma once

#include <nonstd.h>

#include <algorithm>
#include <type_traits>

#include "matpack_mdspan.h"

namespace matpack {
/** Reverse the matpack type elementwise in-place */
template <mut_any_md Self>
constexpr auto reverse_inplace(Self&& self) {
  return std::ranges::reverse(elemwise_range(self));
}

/** Reverse the matpack type by copy
 * 
 * Wraps and calls reverse_inplace(), returning the copied md-type
 */
template <any_md Self>
constexpr auto reverse(const Self& self) {
  using T    = value_type<Self>;
  using md_t = std::
      conditional_t<any_cdata<Self>, std::remove_cvref_t<Self>, data_t<T, 1>>;

  md_t out{self};
  reverse_inplace(out);
  return out;
}

/** Projects a on b */
template <any_md Self>
constexpr auto proj(const Self& a, Self b) {
  b *= dot(a, b) / dot(b, b);
  return b;
}
}  // namespace matpack

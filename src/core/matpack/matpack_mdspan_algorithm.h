#pragma once

#include <compare.h>

#include <algorithm>
#include <vector>

#include "matpack_mdspan_common.h"
#include "matpack_mdspan_common_sizes.h"
#include "matpack_mdspan_common_types.h"

namespace matpack {
template <ranked_md<1> T>
Range sorted_range(const T& container,
                   const value_type<T>& lower_bound,
                   const value_type<T>& upper_bound) {
  auto low = std::ranges::lower_bound(container, lower_bound);
  auto upp = std::ranges::upper_bound(low, container.end(), upper_bound);
  return {static_cast<Index>(std::distance(container.begin(), low)),
          static_cast<Index>(std::distance(low, upp))};
}

//! Divide N elements into n ranges for OpenMP parallelization, returning the
//! ranges
std::vector<Range> omp_offset_count(const Index N, const Index n);

template <any_md T>
data_t<value_type<T>, 1 + rank<T>()> create(const std::vector<T>& vec) {
  assert(not vec.empty());
  assert(not stdr::all_of(vec, Cmp::eq(vec[0].shape()), &T::shape));

  constexpr Size new_rank = 1 + rank<T>();

  std::array<Index, new_rank> new_sizes;
  new_sizes[0] = static_cast<Index>(vec.size());
  for (Size i = 0; i < rank<T>(); i++) {
    new_sizes[i + 1] = static_cast<Index>(vec[0].extent(i));
  }

  data_t<value_type<T>, new_rank> out(new_sizes);
  for (Size i = 0; i < vec.size(); i++) out[i] = vec[i];
  return out;
}
}  // namespace matpack

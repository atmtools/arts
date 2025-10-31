#pragma once

#include <algorithm>

#include "matpack_mdspan_common.h"

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
}  // namespace matpack

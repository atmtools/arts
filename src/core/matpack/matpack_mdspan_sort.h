#pragma once

#include <configtypes.h>

#include <algorithm>
#include <functional>
#include <iterator>
#include <ranges>

namespace matpack {
template <std::random_access_iterator Iter,
          typename Compare = std::ranges::less,
          typename Proj    = std::identity>
constexpr void sort(Iter first, Iter last, Compare comp = {}, Proj proj = {}) {
  using std::iter_swap;

  const auto dist = last - first;
  if (dist < 2) return;

  if (dist < 16) {
    for (auto i = first; i != last; ++i) {
      auto min  = i;
      auto minv = std::invoke(proj, *min);

      for (auto j = i + 1; j != last; ++j) {
        auto&& v = std::invoke(proj, *j);
        if (comp(v, minv)) {
          min  = j;
          minv = std::move(v);
        }
      }

      if (min != i) iter_swap(i, min);
    }
    return;
  }

  auto pivot_pos = first + dist / 2;
  iter_swap(pivot_pos, last - 1);
  const auto pivot_val = std::invoke(proj, *(last - 1));
  auto store           = first;
  for (auto it = first; it != last - 1; ++it) {
    if (comp(std::invoke(proj, *it), pivot_val)) {
      iter_swap(it, store);
      ++store;
    }
  }
  iter_swap(store, last - 1);
  sort(first, store, comp, proj);
  sort(store + 1, last, comp, proj);
}

template <std::ranges::random_access_range R,
          typename Compare = std::ranges::less,
          typename Proj    = std::identity>
constexpr void sort(R&& r, Compare comp = {}, Proj proj = {}) {
  sort(r.begin(), r.end(), comp, proj);
}
}  // namespace matpack
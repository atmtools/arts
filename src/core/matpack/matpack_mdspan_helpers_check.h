#pragma once

#include <compare.h>
#include <nonstd.h>

#include <utility>
#include <vector>

#include "matpack_mdspan.h"
#include "matpack_mdspan_helpers_grid_t.h"

namespace matpack {
template <Size N>
struct integer_helper {
  std::array<Index, N> shape;
  constexpr integer_helper(std::array<Index, N> shape) : shape(shape) {}
  template <integral... ind>
  constexpr integer_helper(ind&&... i)
    requires(sizeof...(ind) == N)
      : shape({static_cast<Index>(i)...}) {}
};

template <integral... ind>
integer_helper(ind&&...) -> integer_helper<sizeof...(ind)>;

template <Size N>
integer_helper(std::array<Index, N>) -> integer_helper<N>;

template <Size N, ranked_md<N> First, ranked_md<N>... Rest>
constexpr bool same_shape(const integer_helper<N>& a,
                          const First& b,
                          const Rest&... c) {
  return (a.shape == b.shape()) and ((a.shape == c.shape()) and ...);
}

template <Size N, ranked_md<N> Orig, ranked_md<N> First, ranked_md<N>... Rest>
constexpr bool same_shape(const Orig& a, const First& b, const Rest&... c) {
  return same_shape(integer_helper<N>{a.shape()}, b, c...);
}

template <Size N, ranked_md<N> B, ranked_md<N>... C>
constexpr bool all_same_shape(const integer_helper<N>& a,
                              const std::vector<B>& b,
                              const std::vector<C>&... c) {
  const auto t = [shape = a.shape](auto& x) { return x.shape() == shape; };
  return stdr::all_of(b, t) and (stdr::all_of(c, t) and ...);
}

template <Size N, ranked_md<N> B, ranked_md<N>... C>
constexpr bool all_same_shape(const ranked_md<N> auto& a,
                              const std::vector<B>& b,
                              const std::vector<C>&... c) {
  return all_same_shape(integer_helper<N>{a.shape()}, b, c...);
}

constexpr bool is_increasing(const exact_md<Numeric, 1> auto& a) {
  return AscendingGrid::is_sorted(a);
}

constexpr bool is_decreasing(const exact_md<Numeric, 1> auto& a) {
  return DescendingGrid::is_sorted(a);
}
}  // namespace matpack

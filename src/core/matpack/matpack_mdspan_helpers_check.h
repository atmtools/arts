#pragma once

#include <compare.h>
#include <nonstd.h>

#include <concepts>
#include <vector>

#include "matpack_mdspan_common_types.h"
#include "matpack_mdspan_helpers_grid_t.h"

namespace matpack {
template <typename T>
concept any_md_compat = any_md<T> or requires(T a) { a.shape(); };

template <typename T, Size N>
concept ranked_md_compat = ranked_md<T, N> or rank<T>() == N;

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

template <any_md_compat First, ranked_md_compat<rank<First>()>... Rest>
constexpr bool same_shape(const integer_helper<rank<First>()>& sz,
                          const First& b,
                          const Rest&... c) {
  return (sz.shape == b.shape()) and ((sz.shape == c.shape()) and ...);
}

template <any_md_compat Orig,
          ranked_md_compat<rank<Orig>()> First,
          ranked_md_compat<rank<Orig>()>... Rest>
constexpr bool same_shape(const Orig& a, const First& b, const Rest&... c) {
  return same_shape(integer_helper<rank<Orig>()>{a.shape()}, b, c...);
}

template <any_md_compat B, ranked_md_compat<rank<B>()>... C>
constexpr bool all_same_shape(const integer_helper<rank<B>()>& sz,
                              const std::vector<B>& b,
                              const std::vector<C>&... c) {
  const auto t = [shape = sz.shape](auto& x) { return x.shape() == shape; };
  return stdr::all_of(b, t) and (stdr::all_of(c, t) and ...);
}

template <any_md_compat A,
          ranked_md_compat<rank<A>()> B,
          ranked_md_compat<rank<A>()>... C>
constexpr bool all_same_shape(const A& a,
                              const std::vector<B>& b,
                              const std::vector<C>&... c) {
  return all_same_shape(integer_helper<rank<A>()>{a.shape()}, b, c...);
}

constexpr bool is_increasing(const exact_md<Numeric, 1> auto& a) {
  return AscendingGrid::is_sorted(a);
}

constexpr bool is_decreasing(const exact_md<Numeric, 1> auto& a) {
  return DescendingGrid::is_sorted(a);
}
}  // namespace matpack

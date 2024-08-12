#pragma once

#include <matpack.h>

constexpr bool all_same_size(const auto& x, const auto&... xs) {
  return ((x.size() == xs.size()) and ...);
}

template <Size N, typename... Ts>
constexpr bool all_same_shape(const std::array<Index, N>& s,
                              const std::vector<Ts>&... xs) {
  return (std::ranges::all_of(xs, Cmp::eq(s), &Ts::shape) and ...);
}

template <typename T, typename... Ts>
constexpr bool all_same_shape(const T& s, const std::vector<Ts>&... xs) {
  return all_same_shape(s.shape(), xs...);
}

template <typename T, typename... Ts>
constexpr bool same_shape(const T& s, const Ts&... xs) {
  return ((s.shape() == xs.shape()) and ...);
}

template <Size N, typename... Ts>
constexpr bool same_shape(const std::array<Index, N>& s, const Ts&... xs) {
  return ((s == xs.shape()) and ...);
}

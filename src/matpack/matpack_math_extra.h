#pragma once

/*
This file only exist to allow standard library
interactions to work as matpack interactions

It is intended to be included sparingly and
all functionality should mimic the matpack
functionality, but not all matpack functionality
*/

#include "matpack_concepts.h"

#include <algorithm>
#include <concepts>
#include <numeric>

template <typename T>
concept standard_iterable = not matpack::any_matpack_type<T> and matpack::rank<T>() == 1 and requires (T a) {
  { a.begin() } -> std::random_access_iterator;
  { a.end() } -> std::random_access_iterator;
};

/** Find minimum of x by reduction */
constexpr auto min(standard_iterable auto &&x) {
  return std::reduce(
      x.begin(), x.end(),
      std::numeric_limits<std::remove_cvref_t<decltype(x[0])>>::max(),
      [](auto a, auto b) { return a < b ? a : b; });
}

/** Find miximum of x by reduction */
constexpr auto max(standard_iterable auto &&x) {
  return std::reduce(
      x.begin(), x.end(),
      std::numeric_limits<std::remove_cvref_t<decltype(x[0])>>::lowest(),
      [](auto a, auto b) { return a > b ? a : b; });
}

/** Find minimum and maximum of x by reduction */
constexpr auto minmax(standard_iterable auto&& x) {
  return std::pair{min(x), max(x)};
}

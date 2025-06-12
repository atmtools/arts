#pragma once

#include <configtypes.h>
#include <debug.h>

#include <array>
#include <limits>
#include <tuple>
#include <type_traits>
#include <utility>

#include "matpack_mdspan.h"

namespace matpack {
namespace {
template <char c, std::array cs, Size N = cs.size()>
consteval Size position_first() {
  for (Size i = 0; i < N; ++i) {
    if (cs[i] == c) {
      return i;
    }
  }
  return std::numeric_limits<Size>::max();
}

template <char c, std::array cs, Size N = cs.size()>
consteval std::array<char, N - 1> drop_first() {
  std::array<char, N - 1> result{};
  bool done = false;
  for (Size i = 0; i < N; ++i) {
    if (cs[i] == c and not done) {
      done = true;
    } else {
      result[i - done] = cs[i];
    }
  }
  return result;
}

template <char c, std::array cs, Size N = cs.size()>
consteval Size count_char() {
  Size n = 0;
  for (Size i = 0; i < N; ++i) {
    n += cs[i] == c;
  }
  return n;
}

template <char c,
          std::array cs,
          typename T,
          Size N = cs.size(),
          Size n = count_char<c, cs>()>
constexpr decltype(auto) reduce_mdrank(T&& arr, Index i [[maybe_unused]]) {
  if constexpr (n == 0 or N == 0) {
    return std::forward<T>(arr);
  } else if constexpr (n == 1) {
    if constexpr (N == 1 or cs[0] == c) {
      return arr[i];
    } else {
      return sub<position_first<c, cs>(), N>(arr, i);
    }
  } else {
    return reduce_mdrank<c, drop_first<c, cs>()>(
        sub<position_first<c, cs>(), N>(arr, i), i);
  }
}

template <char c, std::array cs>
consteval auto reduce_charrank() {
  if constexpr (count_char<c, cs>() == 0) {
    return cs;
  } else {
    return reduce_charrank<c, drop_first<c, cs>()>();
  }
}

template <std::array cs, Size N = cs.size()>
consteval bool empty() {
  return N == 0;
}

template <std::array cf,
          std::array... cs,
          Size N = cf.size(),
          Size M = sizeof...(cs)>
consteval char find_first_char() {
  if constexpr (N == 0 and M == 0) {
    return '\0';
  } else if constexpr (N != 0) {
    return cf.front();
  } else if constexpr (M != 0) {
    return find_first_char<cs...>();
  } else {
    return '\0';
  }
}

template <char c, std::array cf, std::array... cs, Size N = cf.size()>
constexpr Size mddimsize(const auto& xf, const auto&... xs) {
  if constexpr (N > 0 and find_first_char<cf>() == c) {
    constexpr Size dim = position_first<c, cf>();
    return static_cast<Size>(std::get<dim>(mdshape(xf)));
  } else if constexpr (sizeof...(cs) == 0) {
    return std::numeric_limits<Size>::max();
  } else {
    return mddimsize<c, cs...>(xs...);
  }
}

template <std::array... cs>
constexpr bool good_sizes(const auto&... xs) {
  constexpr char first_char = find_first_char<cs...>();
  const Size n              = mddimsize<first_char, cs...>(xs...);
  return ((mddimsize<first_char, cs>(xs) == std::numeric_limits<Size>::max() or
           n == mddimsize<first_char, cs>(xs)) and
          ...);
}
}  // namespace

template <typename T, std::array... cs>
constexpr T einsum_reduce(const auto&... xs) {
  assert((good_sizes<cs...>(xs...)));

  if constexpr ((empty<cs>() and ...)) {
    return (static_cast<T>(xs) * ...);
  } else {
    constexpr char first_char = find_first_char<cs...>();
    const Size n              = mddimsize<first_char, cs...>(xs...);

    T sum{};
    for (Size i = 0; i < n; ++i) {
      sum += einsum_reduce<T, reduce_charrank<first_char, cs>()...>(
          reduce_mdrank<first_char, cs>(xs, i)...);
    }
    return sum;
  }
}

namespace {
template <std::array cr, std::array cs>
consteval bool copy_chars() {
  if constexpr (cr.size() != cs.size()) {
    return false;
  } else {
    return cr == cs;
  }
}

template <std::array... c>
consteval bool copy_chars()
  requires(sizeof...(c) != 2)
{
  return false;
}

template <std::array A,
          std::array B,
          std::array C,
          Size NA = A.size(),
          Size NB = B.size(),
          Size NC = C.size()>
consteval bool multiply_chars() {
  if constexpr (NA == 2 and NB == 2 and NC == 2) {
    return A[0] == B[0] and A[1] == C[1] and B[1] == C[0];
  } else if (NA == 1 and NB == 2 and NC == 1) {
    return A[0] == B[0] and B[1] == C[0];
  } else {
    return false;
  }
}

template <std::array... c>
consteval bool multiply_chars()
  requires(sizeof...(c) != 3)
{
  return false;
}

template <typename T>
constexpr void transform_reduce(T& xr,
                                const any_md auto& x1,
                                const any_md auto& x2) {
  xr = std::transform_reduce(
      x1.elem_begin(), x1.elem_end(), x2.elem_begin(), T{});
}

constexpr void copy_arrs(auto&& xr, const auto& xs) { xr = xs; }

template <std::array cr, std::array... cs>
constexpr void einsum_arr(auto&& xr, const auto&... xs) {
  assert((good_sizes<cr, cs...>(xr, xs...)));
  assert((good_sizes<cr, cs...>(xr, xs...)));
  if constexpr (empty<cr>()) {
    xr = einsum_reduce<std::remove_cvref_t<decltype(xr)>, cs...>(xs...);
  }

  else if constexpr (copy_chars<cr, cs...>()) {
    copy_arrs(xr, xs...);
  }

  else if constexpr (multiply_chars<cr, cs...>()) {
    mult(xr, xs...);
  }

  else {
    constexpr char first_char = find_first_char<cr>();
    const Size n              = mddimsize<first_char, cr>(xr);
    for (Size i = 0; i < n; i++) {
      einsum_arr<reduce_charrank<first_char, cr>(),
                 reduce_charrank<first_char, cs>()...>(
          reduce_mdrank<first_char, cr>(xr, i),
          reduce_mdrank<first_char, cs>(xs, i)...);
    }
  }
}

template <std::array cr, std::array... cs>
consteval bool einsum_arr_optpath() {
  if constexpr (empty<cr>()) {
    return false;
  }

  else if constexpr (copy_chars<cr, cs...>()) {
    return true;
  }

  else if constexpr (multiply_chars<cr, cs...>()) {
    return true;
  }

  else {
    constexpr char first_char = find_first_char<cr>();
    return einsum_arr_optpath<reduce_charrank<first_char, cr>(),
                              reduce_charrank<first_char, cs>()...>();
  }
}

template <Size N>
struct string_literal {
  constexpr string_literal(const char (&in)[N]) { std::copy_n(in, N, str); }
  char str[N];

  constexpr std::array<char, N - 1> to_array() const {
    std::array<char, N - 1> out{};
    std::copy_n(str, N - 1, out.begin());
    return out;
  }
};
}  // namespace

template <string_literal sr,
          string_literal... si,
          Size N = sr.to_array().size()>
constexpr void einsum(auto&& xr, const auto&... xi)
  requires(sizeof...(si) == sizeof...(xi))
{
  einsum_arr<sr.to_array(), si.to_array()...>(std::forward<decltype(xr)>(xr),
                                              xi...);
}

template <typename T,
          string_literal... s,
          Size N = std::get<0>(std::tuple{s...}).to_array().size()>
constexpr T einsum(std::array<Index, N> sz, const auto&... xi)
  requires(sizeof...(s) == sizeof...(xi) + 1 and 1 != sizeof...(s))
{
  if constexpr (N == 0) {
    T out{};
    einsum<s...>(out, xi...);
    return out;
  } else {
    T out(sz);
    einsum<s...>(out, xi...);
    return out;
  }
}

template <string_literal... s>
consteval bool einsum_optpath() {
  return einsum_arr_optpath<s.to_array()...>();
}
}  // namespace matpack

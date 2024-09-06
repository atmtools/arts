#pragma once

#include <matpack_concepts.h>
#include <matpack_iter.h>
#include <matpack_math.h>

#include <limits>
#include <source_location>
#include <tuple>
#include <type_traits>
#include <utility>

#include "debug.h"

namespace matpack {
namespace detail {
template <char c, std::array cs, std::size_t N = cs.size()>
consteval std::size_t position_first() {
  for (std::size_t i = 0; i < N; ++i) {
    if (cs[i] == c) {
      return i;
    }
  }
  return std::numeric_limits<std::size_t>::max();
}

template <char c, std::array cs, std::size_t N = cs.size()>
consteval std::array<char, N - 1> drop_first() {
  std::array<char, N - 1> result{};
  bool done = false;
  for (std::size_t i = 0; i < N; ++i) {
    if (cs[i] == c and not done) {
      done = true;
    } else {
      result[i - done] = cs[i];
    }
  }
  return result;
}

template <char c, std::array cs, std::size_t N = cs.size()>
consteval std::size_t count_char() {
  std::size_t n = 0;
  for (std::size_t i = 0; i < N; ++i) {
    n += cs[i] == c;
  }
  return n;
}

template <char c,
          std::array cs,
          std::size_t N = cs.size(),
          std::size_t n = count_char<c, cs>()>
constexpr decltype(auto) reduce_mdrank(auto&& arr, Index i [[maybe_unused]]) {
  if constexpr (n == 0 or N == 0) {
    return arr;
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

template <std::array cs, std::size_t N = cs.size()>
consteval bool empty() {
  return N == 0;
}

template <std::array cf,
          std::array... cs,
          std::size_t N = cf.size(),
          std::size_t M = sizeof...(cs)>
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

template <char c, std::array cf, std::array... cs, std::size_t N = cf.size()>
constexpr std::size_t mddimsize(const auto& xf, const auto&... xs) {
  if constexpr (N > 0 and find_first_char<cf>() == c) {
    constexpr std::size_t dim = position_first<c, cf>();
    return static_cast<std::size_t>(std::get<dim>(mdshape(xf)));
  } else if constexpr (sizeof...(cs) != 0) {
    return mddimsize<c, cs...>(xs...);
  } else {
    return std::numeric_limits<std::size_t>::max();
  }
}

constexpr std::string shape(const any_matpack_type auto& x) {
  std::ostringstream os;
  std::string_view sep = "(";
  for (auto c : x.shape()) {
    os << std::exchange(sep, ", ") << c;
  }
  os << ')';
  return os.str();
}

constexpr std::string shape(const auto&) { return "()"; }

template <std::array... cs>
std::string error_msg(const auto&... xs) {
  std::ostringstream os;
  os << "operands could not be broadcast together with shapes ";
  ((os << shape(xs) << ' '), ...);
  return os.str();
}

template <std::array... cs>
constexpr bool good_sizes(const auto&... xs) {
  constexpr char first_char = detail::find_first_char<cs...>();
  const std::size_t n = detail::mddimsize<first_char, cs...>(xs...);
  return ((detail::mddimsize<first_char, cs>(xs) ==
               std::numeric_limits<std::size_t>::max() or
           n == detail::mddimsize<first_char, cs>(xs)) and
          ...);
}
}  // namespace detail

template <typename T, std::array... cs>
constexpr T einsum_reduce(const auto&... xs) {
  ARTS_ASSERT((detail::good_sizes<cs...>(xs...)),
              "einsum_reduce: ",
              detail::error_msg<cs...>(xs...))

  if constexpr ((detail::empty<cs>() and ...)) {
    return static_cast<T>((xs * ...));
  } else {
    constexpr char first_char = detail::find_first_char<cs...>();
    const std::size_t n = detail::mddimsize<first_char, cs...>(xs...);

    T sum{};
    for (std::size_t i = 0; i < n; ++i) {
      sum += einsum_reduce<T, detail::reduce_charrank<first_char, cs>()...>(
          detail::reduce_mdrank<first_char, cs>(xs, i)...);
    }
    return sum;
  }
}

namespace detail {
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
          std::size_t NA = A.size(),
          std::size_t NB = B.size(),
          std::size_t NC = C.size()>
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
                                const any_matpack_type auto& x1,
                                const any_matpack_type auto& x2) {
  xr = std::transform_reduce(
      x1.elem_begin(), x1.elem_end(), x2.elem_begin(), T{});
}

constexpr void copy_arrs(auto&& xr, const auto& xs) { xr = xs; }

template <std::array cr, std::array... cs>
constexpr void einsum_arr(auto&& xr, const auto&... xs) {
  assert((detail::good_sizes<cr, cs...>(xr, xs...)));
  ARTS_ASSERT((detail::good_sizes<cr, cs...>(xr, xs...)),
              "einsum: ",
              detail::error_msg<cr, cs...>(xr, xs...))
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
    const std::size_t n = mddimsize<first_char, cr>(xr);
    for (std::size_t i = 0; i < n; i++) {
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

template <std::size_t N>
struct string_literal {
  constexpr string_literal(const char (&in)[N]) { std::copy_n(in, N, str); }
  char str[N];

  constexpr std::array<char, N - 1> to_array() const {
    std::array<char, N - 1> out{};
    std::copy_n(str, N - 1, out.begin());
    return out;
  }
};
}  // namespace detail

template <detail::string_literal sr,
          detail::string_literal... si,
          std::size_t N = sr.to_array().size()>
constexpr void einsum(auto&& xr, const auto&... xi)
  requires(sizeof...(si) == sizeof...(xi))
{
  detail::einsum_arr<sr.to_array(), si.to_array()...>(
      std::forward<decltype(xr)>(xr), xi...);
}

template <typename T,
          detail::string_literal... s,
          std::size_t N = std::get<0>(std::tuple{s...}).to_array().size()>
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

template <detail::string_literal... s>
consteval bool einsum_optpath() {
  return detail::einsum_arr_optpath<s.to_array()...>();
}
}  // namespace matpack

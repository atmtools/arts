#pragma once

#include <nonstd.h>

#include <algorithm>
#include <array>
#include <type_traits>
#include <utility>

#include "compare.h"
#include "configtypes.h"
#include "matpack_concepts.h"
#include "matpack_iter.h"

namespace matpack {
struct par {};

template <char... dim>
struct einsum_id {
  static constexpr Size N = sizeof...(dim);
  static constexpr std::array<char, N> dims{dim...};

  template <char... other>
  using extend_einsum_id = einsum_id<dim..., other...>;

  template <char... other>
  static consteval auto extend(einsum_id<other...>) {
    return extend_einsum_id<other...>{};
  }

  template <char c>
  static consteval Size shape_pos() {
    return static_cast<Size>(
        std::distance(dims.begin(), std::ranges::find(dims, c)));
  }

  template <char c>
  static constexpr Size pos = shape_pos<c>();

  template <char c, class T, char n, char... rem>
  static consteval auto exclude_impl() {
    if constexpr (n == c) {
      return typename T::template extend_einsum_id<rem...>{};
    } else if constexpr (sizeof...(rem) == 0) {
      return T{};
    } else {
      return exclude_impl<c,
                          typename T::template extend_einsum_id<n>,
                          rem...>();
    }
  }

  template <char c>
  static consteval auto exclude() {
    if constexpr (std::ranges::find(dims, c) == dims.end()) {
      return einsum_id<dim...>{};
    } else {
      return exclude_impl<c, einsum_id<>, dim...>();
    }
  }

  template <char c>
  using inner = std::decay_t<decltype(exclude<c>())>;

  //! Selects the p-th element of the array x if c is to be selected
  template <char c, typename T>
  static constexpr decltype(auto) select_if(T& x, Index p) {
    if constexpr (pos<c> < N) {
      if constexpr (N == 1) {
        return x[p];
      } else {
        return sub<pos<c>, N>(x, p);
      }
    } else {
      return x;
    }
  }

  //! Selects the size of the dimension c if it exists, otherwise it is 0
  template <char c, typename T>
  static constexpr Size dimsize_if(const T& x) {
    if constexpr (pos<c> < N) {
      return static_cast<Size>(std::get<pos<c>>(mdshape(x)));
    } else {
      return 0;
    }
  }

  //! What is the leftmost character, if any
  static consteval char next_char() {
    if constexpr (N > 0) {
      return dims.front();
    } else {
      return ' ';
    }
  }

  //! Built-in tests
  static consteval bool no_duplicate() {
    std::array<char, N> x{dim...};
    std::ranges::sort(x);
    return std::ranges::adjacent_find(x) == x.end();
  }
  static consteval bool no_space() {
    return std::ranges::find(dims, ' ') == dims.end();
  }
  static_assert(
      no_duplicate(),
      "duplicate dimensions in einsum_id, not supported by einsum (yet)");
  static_assert(no_space(),
                "space character is reserved for internal use only");
};

template <typename T>
concept einsumable = T::no_duplicate() and T::no_space();

template <einsumable out, einsumable... in>
struct einsum_dims {
  static constexpr Size N = 1 + sizeof...(in);
  static constexpr bool remchar = out::N > 0 or ((in::N > 0) or ...);
  static constexpr bool reduce = out::N == 0;

  static consteval char find_next_char() {
    if constexpr (remchar) {
      constexpr std::array<char, N> all_chars{out::next_char(),
                                              (in::next_char())...};
      for (auto c : all_chars) {
        if (c != ' ') {
          return c;
        }
      }
      std::unreachable();
    } else {
      return ' ';
    }
  }

  static constexpr char next = find_next_char();

  template <typename... Ts>
  static constexpr Size inner_dimsize(Ts&&... x) {
    Size sz = 0;
    for (Size s : {in::template dimsize_if<next>(x)...}) {
      if (s != 0) {
        if (sz == 0) {
          sz = s;
        } else if (sz != s) {
          throw std::runtime_error("einsum_dims: dimension mismatch");
        }
      }
    }
    return sz;
  }

  using inner = einsum_dims<typename out::template inner<next>,
                            typename in::template inner<next>...>;

  template <typename T, typename... Ts>
  static T reduce_einsum(const Ts&... x)
    requires(sizeof...(Ts) == sizeof...(in))
  {
    if constexpr (next == ' ') {
      return (x * ...);
    } else {
      const Size n = inner_dimsize(x...);
      T x0 = T{};
      for (Size i = 0; i < n; i++) {
        x0 += inner::template reduce_einsum<T>(
            in::template select_if<next>(x, i)...);
      }
      return x0;
    }
  }

  template <typename T, typename... Ts>
  static T reduce_einsum(par, const Ts&... x)
    requires(sizeof...(Ts) == sizeof...(in))
  {
    if constexpr (next == ' ') {
      return (x * ...);
    } else {
      const Size n = inner_dimsize(x...);
      T x0 = T{};
#pragma omp parallel for
      for (Size i = 0; i < n; i++) {
        x0 += inner::template reduce_einsum<T>(
            in::template select_if<next>(x, i)...);
      }
      return x0;
    }
  }

  template <typename T, typename... Ts>
  static void einsum(T&& x0, const Ts&... x)
    requires(sizeof...(Ts) == sizeof...(in))
  {
    if constexpr (reduce) {
      x0 = reduce_einsum<std::decay_t<T>>(x...);
    } else {
      const Size n = out::template dimsize_if<next>(x0);
      for (Size i = 0; i < n; i++) {
        inner::einsum(out::template select_if<next>(x0, i),
                      in::template select_if<next>(x, i)...);
      }
    }
  }

  template <typename T, typename... Ts>
  static void einsum(par, T&& x0, const Ts&... x)
    requires(sizeof...(Ts) == sizeof...(in))
  {
    if constexpr (reduce) {
      x0 = reduce_einsum<std::decay_t<T>>(par{}, x...);
    } else {
      const Size n = out::template dimsize_if<next>(x0);
#pragma omp parallel for
      for (Size i = 0; i < n; i++) {
        inner::einsum(out::template select_if<next>(x0, i),
                      in::template select_if<next>(x, i)...);
      }
    }
  }
};

template <std::array s>
struct einsum_id_converter {
  static consteval Size sz() { return s.size(); }

  static consteval std::array<char, sz() - 1> red() {
    std::array<char, sz() - 1> x;
    for (Size i = 0; i < sz() - 1; i++) {
      x[i] = s[i + 1];
    }
    return x;
  }

  static consteval auto convert() {
    if constexpr (sz() == 0) {
      return einsum_id<>{};
    } else if constexpr (sz() == 1) {
      return einsum_id<s[0]>{};
    } else {
      return einsum_id<s[0]>::extend(einsum_id_converter<red()>::convert());
    }
  }

  using converted = decltype(convert());
};

template <std::array... s>
constexpr void einsum(auto&&... x) {
  einsum_dims<typename einsum_id_converter<s>::converted...>::einsum(x...);
}
}  // namespace matpack

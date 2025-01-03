#pragma once

#include <configtypes.h>
#include <format_tags.h>

#include <array>
#include <concepts>
#include <cstdlib>
#include <experimental/mdspan>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>

#include "matpack_mdspan_common_sizes.h"
#include "matpack_mdspan_common_types.h"

namespace matpack {
namespace stdx = std::experimental;
namespace stdr = std::ranges;
namespace stdv = std::ranges::views;

using Joker                 = stdx::full_extent_t;
constexpr inline auto joker = stdx::full_extent;

struct [[nodiscard]] StridedRange {
  Index offset;
  Index nelem;
  Index stride;

  template <integral Offset, integral Nelem, integral Stride>
  constexpr StridedRange(Offset i0, Nelem n, Stride d)
      : offset(static_cast<Index>(i0)),
        nelem(static_cast<Index>(n)),
        stride(static_cast<Index>(d)) {}

  [[nodiscard]] constexpr stdx::strided_slice<Index, Index, Index>
  to_strided_slice() const {
    return {
        .offset = offset, .extent = 1 + (nelem - 1) * stride, .stride = stride};
  }
};

struct [[nodiscard]] Range {
  Index offset;
  Index nelem;

  template <integral Offset, integral Nelem>
  constexpr Range(Offset i0, Nelem n)
      : offset(static_cast<Index>(i0)), nelem(static_cast<Index>(n)) {}

  [[nodiscard]] constexpr stdx::
      strided_slice<Index, Index, std::integral_constant<Index, 1>>
      to_strided_slice() const {
    return {.offset = offset, .extent = nelem, .stride = {}};
  }
};

template <typename T>
concept exhaustive_access_operator =
    integral<T> or std::is_same_v<std::remove_cvref_t<T>, Joker> or
    std::is_same_v<std::remove_cvref_t<T>, Range>;

template <typename T>
concept strided_access_operator =
    std::is_same_v<std::remove_cvref_t<T>, StridedRange>;

template <typename T>
concept access_operator =
    exhaustive_access_operator<T> or strided_access_operator<T>;

template <access_operator L, access_operator... Acc>
consteval bool left_exhaustive() {
  if constexpr (strided_access_operator<L>) {
    return false;
  } else if constexpr (exhaustive_access_operator<L> and
                       (std::same_as<std::remove_cvref_t<Acc>, Joker> and
                        ...)) {
    return true;
  } else {
    return integral<L> and left_exhaustive<Acc...>();
  }
}

template <class T, Size N, bool is_exhaustive, access_operator... Acc>
  requires(not std::is_reference_v<T>)
struct left_mdsel {
  static constexpr bool access_exhaustive =
      is_exhaustive and left_exhaustive<Acc...>();

  static constexpr Size reduces_rank = (integral<Acc> + ...);

  static_assert(N >= reduces_rank, "Rank reduction is too large");

  using type = std::conditional_t<
      N - reduces_rank == 0,
      T&,
      std::conditional_t<access_exhaustive,
                         view_t<T, N - reduces_rank>,
                         strided_view_t<T, N - reduces_rank>>>;
};

template <class T, Size N, bool is_exhaustive, access_operator... Acc>
using left_mdsel_t = typename left_mdsel<T, N, is_exhaustive, Acc...>::type;

template <Size N>
consteval std::array<Joker, N> jokers() {
  std::array<Joker, N> out;
  out.fill(joker);
  return out;
}

template <typename Acc>
constexpr decltype(auto) to_base(Acc&& i) {
  if constexpr (std::is_same_v<std::remove_cvref_t<Acc>, StridedRange> or
                std::is_same_v<std::remove_cvref_t<Acc>, Range>) {
    return std::forward<Acc>(i).to_strided_slice();
  } else {
    return std::forward<Acc>(i);
  }
}

template <Size N, access_operator... Acc>
constexpr decltype(auto) acc(Acc&&... i)
  requires(sizeof...(Acc) <= N)
{
  if constexpr (N == sizeof...(Acc)) {
    return std::forward_as_tuple(std::forward<Acc>(i)...);
  } else {
    return std::tuple_cat(std::forward_as_tuple(std::forward<Acc>(i)...),
                          jokers<N - sizeof...(Acc)>());
  }
}

template <Size N, typename Self, access_operator... Acc>
constexpr decltype(auto) left_sub(Self&& s, Acc&&... i) {
  return std::apply(
      [&s]<access_operator... AccT>(AccT&&... x) {
        if constexpr (any_md<Self>) {
          return stdx::submdspan(s.base_md(),
                                 to_base(std::forward<AccT>(x))...);
        } else {
          return stdx::submdspan(s, to_base(std::forward<AccT>(x))...);
        }
      },
      acc<N>(std::forward<Acc>(i)...));
}

//! Create a tuple to access an index surrounded by jokers
template <Size M, Size N, typename T>
[[nodiscard]] constexpr auto tup(T&& i)
  requires(M < N)
{
  constexpr Size l = M;
  constexpr Size r = N - M - 1;
  static_assert(l < N, "Left must be less than N");
  static_assert(r < N, "Right must be less than N");
  constexpr bool has_left  = l > 0;
  constexpr bool has_right = r > 0;

  if constexpr (has_left and has_right) {
    return std::tuple_cat(
        jokers<l>(), std::forward_as_tuple(std::forward<T>(i)), jokers<r>());
  } else if constexpr (has_left) {
    return std::tuple_cat(jokers<l>(),
                          std::forward_as_tuple(std::forward<T>(i)));
  } else if constexpr (has_right) {
    return std::tuple_cat(std::forward_as_tuple(to_base(std::forward<T>(i))),
                          jokers<r>());
  } else {
    static_assert(M == 0, "M must be 0");
    static_assert(N == 1, "N must be 1");
    return std::forward_as_tuple(std::forward<T>(i));
  }
}

template <Size M, Size N, typename Self, access_operator Acc>
constexpr decltype(auto) sub(Self&& s, Acc&& i)
  requires(M < N)
{
  if constexpr (std::same_as<std::remove_cvref_t<Acc>, Joker>) {
    return std::forward<Self>(s);
  } else {
    return std::apply(
        [&s]<access_operator... AccT>(AccT&&... x) {
          if constexpr (any_md<Self>) {
            return stdx::submdspan(s.base_md(),
                                   to_base(std::forward<AccT>(x))...);
          } else {
            return stdx::submdspan(s, to_base(std::forward<AccT>(x))...);
          }
        },
        tup<M, N>(i));
  }
}

//! Get a positional value
template <rankable T, Size dim = rank<T>()>
constexpr auto mdvalue(const T& v, const std::array<Index, dim>& pos) {
  if constexpr (any_md<T>) {
    return v.base_md()[pos];
  } else if constexpr (dim == 1) {
    if constexpr (has_index_access<T>)
      return v[pos[0]];
    else
      return *(std::begin(v) + pos[0]);
  } else {
    return std::apply(
        [&v](auto... inds) {
          if constexpr (requires { v(inds...); }) {
            return v(inds...);
          } else {
            return v[inds...];
          }
        },
        pos);
  }
}

//! Test that you can convert the underlying value type of an object to a
//! specific type
template <typename U, typename T>
concept mdvalue_type_compatible = requires(U a) {
  { mdvalue(a, mdshape(a)) } -> std::convertible_to<T>;
};

template <Size N>
constexpr Size mdsize(const std::array<Index, N>& shape) {
  Size size = 1;
  for (Size i = 0; i < N; i++) size *= shape[i];
  return size;
}

template <Size N>
constexpr std::array<Index, N> mdpos(const std::array<Index, N>& shape,
                                     Index i) {
  std::array<Index, N> pos{};

  Index n = mdsize<N>(shape);
  for (Size j = 0; j < N; j++) {
    n              /= shape[j];
    const auto res  = std::div(i, n);
    pos[j]          = res.quot;
    i               = res.rem;
  }

  return pos;
}

template <rankable T, Size dim = rank<T>()>
constexpr auto mdvalue(const T& v,
                       const std::array<Index, dim>& shape,
                       const Index i) {
  return mdvalue(v, mdpos(shape, i));
}
}  // namespace matpack

using StridedRange           = matpack::StridedRange;
using Range                  = matpack::Range;
using Joker                  = matpack::Joker;
inline constexpr Joker joker = matpack::joker;

template <>
struct std::formatter<StridedRange> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const StridedRange& v, FmtContext& ctx) const {
    return tags.format(
        ctx, "sr["sv, v.offset, tags.sep(), v.nelem, tags.sep(), v.stride, ']');
  }
};

template <>
struct std::formatter<Range> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Range& v, FmtContext& ctx) const {
    return tags.format(ctx, "r["sv, v.offset, tags.sep(), v.nelem, ']');
  }
};

template <>
struct std::formatter<Joker> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const Joker&, FmtContext& ctx) const {
    return tags.format(ctx, "joker"sv);
  }
};

#pragma once

#include <debug.h>

#include <algorithm>
#include <concepts>

#include "matpack_mdspan.h"

struct Ascending {
  consteval static std::string_view name() { return "ascending"sv; }
  constexpr bool operator()(Numeric a, Numeric b) const { return a <= b; }
};

struct Descending {
  consteval static std::string_view name() { return "descending"sv; }
  constexpr bool operator()(Numeric a, Numeric b) const { return a >= b; }
};

namespace matpack {
template <typename T>
concept sorting_t = std::same_as<T, Ascending> or std::same_as<T, Descending>;

template <sorting_t Compare>
class grid_t {
  Vector x;

 public:
  [[nodiscard]] constexpr const Vector& vec() const { return x; }

  using value_type = Numeric;

  static constexpr bool is_sorted(const exact_md<Numeric, 1> auto& x) {
    return stdr::is_sorted(x, Compare{});
  }

  static constexpr void assert_sorted(const exact_md<Numeric, 1> auto& x) {
    if (not is_sorted(x))
      throw std::runtime_error(std::format(
          "Expected {} sort.  Input list:\n{:B,}", Compare::name(), x));
  }

  grid_t(Index N) : x(N) {
    constexpr Numeric mi        = std::numeric_limits<Numeric>::lowest();
    constexpr Numeric ma        = std::numeric_limits<Numeric>::max();
    constexpr bool is_ascending = Compare{}(mi, ma);
    constexpr Numeric val0      = is_ascending ? mi : ma;
    constexpr Numeric val1      = is_ascending ? ma : mi;

    Numeric val{val0};
    for (Index i = 0; i < N; i++) {
      x[i] = val;
      val  = std::nextafter(val, val1);
    }

    // Safety check
    assert_sorted(x);
  }

  grid_t(std::initializer_list<Numeric> il) : x(il) { assert_sorted(x); }

  template <typename... Ts>
  explicit constexpr grid_t(Ts&&... ts)
    requires std::constructible_from<Vector, Ts...>
      : x(Vector{std::forward<Ts>(ts)...}) {
    assert_sorted(x);
  }

  template <class Iterator, class Func>
  constexpr grid_t(const Iterator& begin, const Iterator& end, Func fun)
      : x(std::distance(begin, end)) {
    std::transform(begin, end, x.begin(), fun);
    assert_sorted(x);
  }

  constexpr grid_t(Vector&& in) : x(std::move(in)) { assert_sorted(x); }

  constexpr grid_t& operator=(Vector&& in) {
    assert_sorted(in);
    x = std::move(in);
    return *this;
  }

  template <typename T>
  constexpr grid_t& operator=(const T& in) {
    assert_sorted(in);
    x = in;
    return *this;
  }

  constexpr grid_t()                         = default;
  constexpr grid_t(grid_t&&)                 = default;
  constexpr grid_t(const grid_t&)            = default;
  constexpr grid_t& operator=(grid_t&&)      = default;
  constexpr grid_t& operator=(const grid_t&) = default;

  constexpr operator Vector() && { return Vector{std::move(x)}; }
  constexpr operator const Vector&() const { return vec(); }
  constexpr operator ConstVectorView() const { return vec(); }
  constexpr operator StridedConstVectorView() const { return vec(); }
  constexpr operator std::span<const Numeric>() const { return vec(); }

  template <access_operator Op>
  [[nodiscard]] constexpr auto operator[](const Op& op) const {
    return x[op];
  }

  [[nodiscard]] constexpr auto size() const { return x.size(); }
  [[nodiscard]] constexpr auto begin() const { return x.begin(); }
  [[nodiscard]] constexpr auto end() const { return x.end(); }
  [[nodiscard]] constexpr auto empty() const { return x.empty(); }
  [[nodiscard]] constexpr Numeric front() const { return x.front(); }
  [[nodiscard]] constexpr Numeric back() const { return x.back(); }

  [[nodiscard]] constexpr std::array<Index, 1> shape() const {
    return x.shape();
  }
};

template <sorting_t C1, sorting_t C2>
[[nodiscard]] constexpr bool operator==(const grid_t<C1>& x,
                                        const grid_t<C2>& y) {
  return x.vec() == y.vec();
}

template <sorting_t C1, sorting_t C2>
[[nodiscard]] constexpr bool operator!=(const grid_t<C1>& x,
                                        const grid_t<C2>& y) {
  return x.vec() != y.vec();
}

template <sorting_t Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator==(const grid_t<Compare>& x, const md& y) {
  return x.vec() == y;
}

template <sorting_t Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator!=(const grid_t<Compare>& x, const md& y) {
  return x.vec() != y;
}

template <sorting_t Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator==(const md& y, const grid_t<Compare>& x) {
  return y == x.vec();
}

template <sorting_t Compare, exact_md<Numeric, 1> md>
[[nodiscard]] constexpr bool operator!=(const md& y, const grid_t<Compare>& x) {
  return y != x.vec();
}
}  // namespace matpack

using AscendingGrid        = matpack::grid_t<Ascending>;
using DescendingGrid       = matpack::grid_t<Descending>;
using ArrayOfAscendingGrid = std::vector<AscendingGrid>;

template <matpack::sorting_t Compare>
struct std::formatter<matpack::grid_t<Compare>> {
  std::formatter<matpack::strided_view_t<const Numeric, 1>> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::grid_t<Compare>& v,
                              FmtContext& ctx) const {
    return fmt.format(v, ctx);
  }
};

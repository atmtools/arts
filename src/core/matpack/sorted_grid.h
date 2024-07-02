#pragma once

#include <array.h>

#include <algorithm>
#include <concepts>
#include <functional>
#include <ranges>
#include <stdexcept>

#include "configtypes.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"

namespace matpack {
template <class Compare>
class grid;

template <class Compare>
class grid_view {
  grid<Compare>& x;

 public:
  [[nodiscard]] Index size() const { return x.x.size(); }

  grid_view(grid<Compare>& x_) : x(x_) {}

  [[nodiscard]] const Vector& vec() const { return x.x; }

  grid_view& operator=(const grid<Compare>& y) {
    if (&y.x != &x.x) ExhaustiveVectorView{x.x} = y.x;
    return *this;
  }

  grid_view& operator=(const ConstVectorView& y) {
    *this = grid<Compare>{y};
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const grid_view& y) {
    return os << y.x;
  }
};

template <class Compare>
class grid {
  Vector x;

  friend class grid_view<Compare>;

 public:
  [[nodiscard]] const Vector& vec() const { return x; }

  using value_type = Numeric;

  static void assert_sorted(const Vector& x) {
    ARTS_USER_ERROR_IF(
        not std::ranges::is_sorted(x, Compare{}), "Wrong sorting:\n", x);
  }

  //! Unsafe constructor
  grid(Index N) : x(N) {}

  grid(std::initializer_list<Numeric> il) : x(il) { assert_sorted(x); }

  template <typename... Ts>
  explicit grid(Ts&&... ts)
    requires std::constructible_from<Vector, Ts...>
      : x(Vector{std::forward<Ts>(ts)...}) {
    assert_sorted(x);
  }

  grid(Vector&& in) : x(std::move(in)) { assert_sorted(x); }

  grid& operator=(Vector&& in) {
    x = std::move(in);
    assert_sorted(x);
    return *this;
  }

  template <typename T>
  grid& operator=(T&& in) {
    x = std::forward<T>(in);
    assert_sorted(x);
    return *this;
  }

  grid() = default;
  grid(grid&&) = default;
  grid(const grid&) = default;
  grid& operator=(grid&&) = default;
  grid& operator=(const grid&) = default;

  operator const Vector&() const { return x; }
  operator ConstVectorView() const { return x; }
  operator ExhaustiveConstVectorView() const { return x; }

  template <access_operator Op>
  [[nodiscard]] constexpr auto operator()(const Op& op) const {
    return x(op);
  }

  template <access_operator Op>
  [[nodiscard]] constexpr auto operator[](const Op& op) const {
    return x[op];
  }

  [[nodiscard]] constexpr auto nelem() const { return x.nelem(); }
  [[nodiscard]] constexpr auto size() const { return x.size(); }
  [[nodiscard]] constexpr auto begin() const { return x.begin(); }
  [[nodiscard]] constexpr auto end() const { return x.end(); }
  [[nodiscard]] constexpr auto empty() const { return x.empty(); }
  [[nodiscard]] constexpr Numeric front() const { return x.front(); }
  [[nodiscard]] constexpr Numeric back() const { return x.back(); }
  void clear() { x.clear(); }
  void reserve(Size n) { x.reserve(n); }
  void erase(auto... it) { x.erase(it...); }
  constexpr Numeric& emplace_back(Numeric v) {
    return x.size() == 0 ? x.emplace_back(v)
           : Compare{}(v, x.back())
               ? throw std::runtime_error("grid::emplace_back: not sorted")
               : x.emplace_back(v);
  }

  void unsafe_resize(const Index n) { x.resize(n); }
  [[nodiscard]] constexpr auto unsafe_begin() { return x.begin(); }
  [[nodiscard]] constexpr auto unsafe_end() { return x.end(); }
  void unsafe_push_back(Numeric v) { x.push_back(v); }
  Numeric& unsafe_emplace_back(Numeric v) { return x.emplace_back(v); }

  //! Return a pointer to this object's allocated data
  [[nodiscard]] constexpr auto data_handle() const { return x.data_handle(); }

  auto operator<=>(const grid&) const = default;

  friend std::ostream& operator<<(std::ostream& os, const grid& g) {
    return os << g.x;
  }
};
}  // namespace matpack

using AscendingGrid = matpack::grid<std::less_equal<>>;
using AscendingGridView = matpack::grid_view<std::less_equal<>>;
using DescendingGrid = matpack::grid<std::greater_equal<>>;
using ArrayOfAscendingGrid = Array<AscendingGrid>;

namespace matpack {
template <class T>
std::ostream& operator<<(std::ostream& os, const Array<grid<T>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Array<Array<grid<T>>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}

template <class T>
std::ostream& operator<<(std::ostream& os,
                         const Array<Array<Array<grid<T>>>>& a) {
  for (auto& x : a) os << x << '\n';
  return os;
}
}  // namespace matpack

template <class Compare>
struct std::formatter<matpack::grid<Compare>> {
  std::formatter<matpack::matpack_view<Numeric, 1, true, false>> fmt;

  [[nodiscard]] constexpr auto& inner_fmt() { return fmt.inner_fmt(); }
  [[nodiscard]] constexpr auto& inner_fmt() const { return fmt.inner_fmt(); }

  template <typename... Ts>
  constexpr void make_compat(std::formatter<Ts>&... xs) const {
    inner_fmt().make_compat(xs...);
  }

  template <typename U>
  constexpr void compat(const std::formatter<U>& x) {
    x.make_compat(*this);
  }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return fmt.parse(ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::grid<Compare>& v,
                              FmtContext& ctx) const {
    return fmt.format(v, ctx);
  }
};

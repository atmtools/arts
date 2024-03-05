#pragma once

#include <array.h>

#include <concepts>
#include <functional>
#include <ranges>

#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"

namespace matpack {
template <class Compare>
class grid {
  Vector x;

 public:
  using value_type = Numeric;

  void assert_sorted() const {
    ARTS_USER_ERROR_IF(
        not std::ranges::is_sorted(x, Compare{}), "Wrong sorting:\n", x);
  }

  grid(std::initializer_list<Numeric> il) : x(il) { assert_sorted(); }

  template <typename... Ts>
  explicit grid(Ts&&... ts)
    requires std::constructible_from<Vector, Ts...>
      : x(Vector{std::forward<Ts>(ts)...}) {
    assert_sorted();
  }

  grid(Vector&& in) : x(std::move(in)) {
    assert_sorted();
  }

  grid& operator=(Vector&& in) {
    x = std::move(in);
    assert_sorted();
    return *this;
  }

  template <typename T>
  grid& operator=(T&& in) {
    x = std::forward<T>(in);
    assert_sorted();
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
  void erase(auto ... it) { x.erase(it...); }

  void unsafe_resize(const Index n) { x.resize(n); }
  [[nodiscard]] constexpr auto unsafe_begin() { return x.begin(); }
  [[nodiscard]] constexpr auto unsafe_end() { return x.end(); }
  void unsafe_push_back(Numeric v) { x.push_back(v); }
  Numeric& unsafe_emplace_back(Numeric v) {
    return x.emplace_back(v);
  }

  //! Return a pointer to this object's allocated data
  [[nodiscard]] constexpr auto data_handle() const { return x.data_handle(); }

  auto operator<=>(const grid&) const = default;

  friend std::ostream& operator<<(std::ostream& os, const grid& g) {
    return os << g.x;
  }
};
}  // namespace matpack

using AscendingGrid = matpack::grid<std::less_equal<>>;
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

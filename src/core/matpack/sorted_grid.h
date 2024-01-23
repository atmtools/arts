#pragma once

#include <ranges>

#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_view.h"

namespace matpack {
template <bool ascending>
class grid {
  Vector x;

  void assert_sorted() const {
    if constexpr (ascending) {
      ARTS_USER_ERROR_IF(not std::ranges::is_sorted(
                             x, [](auto& a, auto& b) { return a <= b; }),
                         "Must be sorted in strictly ascending order, got:\n",
                         x);
    } else {
      ARTS_USER_ERROR_IF(not std::ranges::is_sorted(
                             x, [](auto& a, auto& b) { return a >= b; }),
                         "Must be sorted in strictly descending order, got:\n",
                         x);
    }
  }

 public:
  grid(Vector in) : x(std::move(in)) { assert_sorted(); }

  grid& operator=(Vector&& in) {
    x = std::move(in);
    assert_sorted();
    return *this;
  }

  grid& operator=(const Vector& in) {
    x = in;
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
  auto operator()(const Op& op) const {
    return op(x);
  }
};
}  // namespace matpack

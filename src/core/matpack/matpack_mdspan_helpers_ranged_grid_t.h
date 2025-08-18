#pragma once

#include <limits>

#include "lagrange_interp.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_helpers_grid_t.h"

namespace matpack {
template <sorting_t Compare,
          Numeric lower_limit,
          Numeric upper_limit,
          bool include_lower,
          bool include_upper>
  requires(lower_limit < upper_limit)
class ranged_grid_t {
  Vector grid;

 public:
  using lag_sorting = std::conditional_t<std::same_as<Compare, Ascending>,
                                         lagrange_interp::ascending_grid_t,
                                         lagrange_interp::descending_grid_t>;

  static constexpr bool is_valid(Numeric value) {
    bool lower = include_lower ? value >= lower_limit : value > lower_limit;
    bool upper = include_upper ? value <= upper_limit : value < upper_limit;
    return lower and upper;
  }

  static void assert_ranged(const exact_md<Numeric, 1> auto& grid) {
    if (grid.empty()) return;

    ARTS_USER_ERROR_IF(
        not grid_t<Compare>::is_sorted(grid) or not is_valid(grid.front()) or
            not is_valid(grid.back()),
        "Wrong range. Expected range [{}, {}) of {} sorting, got:\n{:B,}",
        lower_limit,
        upper_limit,
        Compare::name(),
        grid);
  }

  constexpr ranged_grid_t()                                    = default;
  constexpr ranged_grid_t(const ranged_grid_t&)                = default;
  constexpr ranged_grid_t(ranged_grid_t&&) noexcept            = default;
  constexpr ranged_grid_t& operator=(const ranged_grid_t&)     = default;
  constexpr ranged_grid_t& operator=(ranged_grid_t&&) noexcept = default;

  constexpr ranged_grid_t(Vector&& x) {
    assert_ranged(x);
    grid = std::move(x);
  }

  constexpr ranged_grid_t(const Vector& x) {
    assert_ranged(x);
    grid = x;
  }

  constexpr ranged_grid_t& operator=(Vector&& x) {
    assert_ranged(x);
    grid = std::move(x);
    return *this;
  }

  constexpr ranged_grid_t& operator=(const Vector& x) {
    assert_ranged(x);
    grid = x;
    return *this;
  }

  [[nodiscard]] const Vector& vec() const { return grid; }
  constexpr operator Vector() && { return std::move(grid); }
  constexpr operator const Vector&() const { return grid; }
  constexpr operator ConstVectorView() const { return grid; }
  constexpr operator StridedConstVectorView() const { return grid; }
  constexpr operator std::span<const Numeric>() const { return grid; }

  template <access_operator Op>
  constexpr decltype(auto) operator[](Op&& x) const {
    return grid[std::forward<Op>(x)];
  }

  template <Size N,
            lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(Numeric x) const {
    return lagrange_interp::lag_t<N, transform>{grid, x, lag_sorting{}};
  }
  template <lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(Numeric x, Size N) const {
    return lagrange_interp::lag_t<-1, transform>{grid, x, N, lag_sorting{}};
  }

  template <Size N,
            lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(std::span<const Numeric> x,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") const {
    return lagrange_interp::make_lags<N, transform>(
        grid, x, extrapolation_limit, info);
  }

  template <lagrange_interp::transformer transform = lagrange_interp::identity>
  [[nodiscard]] auto lag(std::span<const Numeric> x,
                         Size N,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") const {
    return lagrange_interp::make_lags<>(grid, x, N, extrapolation_limit, info);
  }
};
}  // namespace matpack

using LongitudeGrid =
    matpack::ranged_grid_t<Ascending, -180.0, 180.0, true, false>;
using LatitudeGrid = matpack::ranged_grid_t<Ascending, -90.0, 90.0, true, true>;
using AltitudeGrid =
    matpack::ranged_grid_t<Ascending,
                           -std::numeric_limits<Numeric>::infinity(),
                           std::numeric_limits<Numeric>::infinity(),
                           false,
                           false>;
using ZenithGrid  = matpack::ranged_grid_t<Descending, 0.0, 180.0, true, true>;
using AzimuthGrid = matpack::ranged_grid_t<Descending, 0.0, 360.0, true, false>;
using FrequencyGrid =
    matpack::ranged_grid_t<Ascending,
                           0.0,
                           std::numeric_limits<Numeric>::infinity(),
                           false,
                           false>;

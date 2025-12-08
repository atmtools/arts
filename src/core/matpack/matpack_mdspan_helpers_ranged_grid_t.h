#pragma once

#include <algorithm>
#include <limits>
#include <stdexcept>

#include "lagrange_interp.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_helpers_grid_t.h"
#include "matpack_mdspan_view_t.h"
#include "xml_io_base.h"
#include "xml_io_stream.h"
#include "xml_io_stream_matpack_mdspan.h"

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
  using lag_sorting_t = std::conditional_t<std::same_as<Compare, Ascending>,
                                           lagrange_interp::ascending_grid_t,
                                           lagrange_interp::descending_grid_t>;

  static constexpr bool is_valid(Numeric value) {
    bool lower = include_lower ? value >= lower_limit : value > lower_limit;
    bool upper = include_upper ? value <= upper_limit : value < upper_limit;
    return lower and upper;
  }

  static void assert_ranged(const exact_md<Numeric, 1> auto& grid) {
    if (grid.empty()) return;

    if (not grid_t<Compare>::is_sorted(grid) or not is_valid(grid.front()) or
        not is_valid(grid.back())) {
      throw std::runtime_error(std::format(
          "Wrong range. Expected range {}{}, {}{} of {} sort, got:\n{:B,}",
          include_lower ? '[' : '(',
          lower_limit,
          upper_limit,
          include_upper ? ']' : ')',
          Compare::name(),
          grid));
    }
  }

  static void assert_ranged(Numeric x) { assert_ranged(ConstVectorView{x}); }

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
  [[nodiscard]] Vector&& rvec() && noexcept { return std::move(grid); }
  constexpr operator Vector() && { return std::move(*this).rvec(); }
  constexpr operator const Vector&() const { return grid; }
  constexpr operator ConstVectorView() const { return grid; }
  constexpr operator StridedConstVectorView() const { return grid; }
  constexpr operator std::span<const Numeric, std::dynamic_extent>() const {
    return grid;
  }

  [[nodiscard]] constexpr Size size() const { return grid.size(); }
  [[nodiscard]] constexpr bool empty() const { return grid.empty(); }

  [[nodiscard]] constexpr auto begin() const { return grid.begin(); }
  [[nodiscard]] constexpr auto end() const { return grid.end(); }

  [[nodiscard]] constexpr Numeric front() const { return grid.front(); }
  [[nodiscard]] constexpr Numeric back() const { return grid.back(); }

  [[nodiscard]] constexpr auto operator<=>(const ranged_grid_t&) const =
      default;

  template <access_operator Op>
  constexpr decltype(auto) operator[](Op&& x) const {
    return grid[std::forward<Op>(x)];
  }

  template <Size N, lagrange_interp::transformer transform>
  [[nodiscard]] auto lag(Numeric x) const {
    return lagrange_interp::lag_t<N, transform>{grid, x, lag_sorting_t{}};
  }

  template <lagrange_interp::transformer transform>
  [[nodiscard]] auto lag(Numeric x, Size N) const {
    return lagrange_interp::lag_t<-1, transform>{grid, x, N, lag_sorting_t{}};
  }

  template <Size N, lagrange_interp::transformer transform>
  [[nodiscard]] auto lag(std::span<const Numeric> x,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") const {
    return lagrange_interp::make_lags<N, transform>(
        grid, x, extrapolation_limit, info);
  }

  template <lagrange_interp::transformer transform>
  [[nodiscard]] auto lag(std::span<const Numeric> x,
                         Size N,
                         Numeric extrapolation_limit = 0.5,
                         const char* info            = "UNNAMED") const {
    return lagrange_interp::make_lags<transform>(
        grid, x, N, extrapolation_limit, info);
  }
};
}  // namespace matpack

//! Guaranteed [-180, 180)
using LonGrid = matpack::ranged_grid_t<Ascending, -180.0, 180.0, true, false>;

//! Guaranteed [-90, 90]
using LatGrid = matpack::ranged_grid_t<Ascending, -90.0, 90.0, true, true>;

//! Guaranteed [0, 180]
using ZenGrid = matpack::ranged_grid_t<Ascending, 0.0, 180.0, true, true>;

//! Guaranteed [0, 360)
using AziGrid = matpack::ranged_grid_t<Ascending, 0.0, 360.0, true, false>;

template <>
struct xml_io_stream_name<LonGrid> {
  static constexpr std::string_view name = "LonGrid"sv;
};

template <>
struct xml_io_stream_name<LatGrid> {
  static constexpr std::string_view name = "LatGrid"sv;
};

template <>
struct xml_io_stream_name<ZenGrid> {
  static constexpr std::string_view name = "ZenGrid"sv;
};

template <>
struct xml_io_stream_name<AziGrid> {
  static constexpr std::string_view name = "AziGrid"sv;
};

template <matpack::sorting_t Compare,
          Numeric lower_limit,
          Numeric upper_limit,
          bool include_lower,
          bool include_upper>
struct xml_io_stream<matpack::ranged_grid_t<Compare,
                                            lower_limit,
                                            upper_limit,
                                            include_lower,
                                            include_upper>> {
  using T = matpack::ranged_grid_t<Compare,
                                   lower_limit,
                                   upper_limit,
                                   include_lower,
                                   include_upper>;

  static constexpr std::string_view type_name = xml_io_stream_name_v<T>;

  static void write(std::ostream& os,
                    const T& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    xml_write_to_stream(os, x.vec(), pbofs, name);
  }

  static void read(std::istream& is, T& x, bifstream* pbifs = nullptr) {
    Vector data;
    xml_read_from_stream(is, data, pbifs);
    x = std::move(data);
  }
};

template <matpack::sorting_t Compare,
          Numeric lower_limit,
          Numeric upper_limit,
          bool include_lower,
          bool include_upper>
struct std::formatter<matpack::ranged_grid_t<Compare,
                                             lower_limit,
                                             upper_limit,
                                             include_lower,
                                             include_upper>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::ranged_grid_t<Compare,
                                                           lower_limit,
                                                           upper_limit,
                                                           include_lower,
                                                           include_upper>& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.vec());
  }
};

#pragma once

#include <xml.h>

#include <concepts>
#include <limits>
#include <stdexcept>

#include "lagrange_interp.h"
#include "matpack_mdspan_data_t.h"
#include "matpack_mdspan_helpers_grid_t.h"
#include "matpack_mdspan_view_t.h"
#include "xml_io_stream_matpack_mdspan.h"

namespace matpack {
//! A limit of a ranged grid is excluded, included, or infinite (no limit)
enum class lim : char { excl, incl, inf };

template <Index ll_, Index ul_, lim il_ = lim::incl, lim iu_ = lim::incl> class ranged_grid_t
    : public grid_t<std::conditional_t<ll_ <= ul_, ascending_t, descending_t>> {
  static_assert(ll_ != ul_ or il_ == lim::inf or iu_ == lim::inf, "Lower and upper limits cannot be equal");

  static constexpr bool ascending = ll_ <= ul_;

  static consteval Numeric limit(Index x, lim l, bool upper) {
    if (l == lim::inf)
      return upper ? std::numeric_limits<Numeric>::infinity() : -std::numeric_limits<Numeric>::infinity();
    return static_cast<Numeric>(x);
  }

  static constexpr bool    include_lower = (ascending ? il_ : iu_) != lim::excl;
  static constexpr bool    include_upper = (ascending ? iu_ : il_) != lim::excl;
  static constexpr Numeric lower_limit   = ascending ? limit(ll_, il_, false) : limit(ul_, iu_, false);
  static constexpr Numeric upper_limit   = ascending ? limit(ul_, iu_, true) : limit(ll_, il_, true);

  using sorter_t = grid_t<std::conditional_t<ascending, ascending_t, descending_t>>;

  sorter_t&       grid() noexcept { return *this; }
  const sorter_t& grid() const noexcept { return *this; }

 public:
  using lag_sorting_t =
      std::conditional_t<ascending, lagrange_interp::ascending_grid_t, lagrange_interp::descending_grid_t>;

  static constexpr bool is_valid(Numeric value) {
    bool lower = include_lower ? value >= lower_limit : value > lower_limit;
    bool upper = include_upper ? value <= upper_limit : value < upper_limit;
    return lower and upper;
  }

  static void assert_ranged(const exact_md<Numeric, 1> auto& grid) {
    if (grid.empty()) return;

    if (not is_valid(grid.front()) or not is_valid(grid.back())) {
      throw std::runtime_error(std::format("Wrong range. Expected range {}{}, {}{}, got:\n{:B,}",
                                           include_lower ? '[' : '(',
                                           lower_limit,
                                           upper_limit,
                                           include_upper ? ']' : ')',
                                           grid));
    }
  }

  static void assert_ranged(Numeric x) { assert_ranged(ConstVectorView{x}); }

  constexpr ranged_grid_t()                                    = default;
  constexpr ranged_grid_t(const ranged_grid_t&)                = default;
  constexpr ranged_grid_t(ranged_grid_t&&) noexcept            = default;
  constexpr ranged_grid_t& operator=(const ranged_grid_t&)     = default;
  constexpr ranged_grid_t& operator=(ranged_grid_t&&) noexcept = default;

  constexpr ranged_grid_t(Vector&& v) {
    assert_ranged(v);
    grid() = std::move(v);
  }

  constexpr ranged_grid_t(const Vector& v) {
    assert_ranged(v);
    grid() = v;
  }

  constexpr ranged_grid_t& operator=(Vector&& v) {
    assert_ranged(v);
    grid() = std::move(v);
    return *this;
  }

  constexpr ranged_grid_t& operator=(const Vector& v) {
    assert_ranged(v);
    grid() = v;
    return *this;
  }

  //! The grid is already sorted, so only its range is checked
  template <std::derived_from<sorter_t> G> constexpr ranged_grid_t(const G& g) {
    assert_ranged(g.vec());
    grid() = static_cast<const sorter_t&>(g);
  }

  //! The grid is already sorted, so only its range is checked
  template <std::derived_from<sorter_t> G> constexpr ranged_grid_t(G&& g) {
    assert_ranged(g.vec());
    grid() = static_cast<sorter_t&&>(g);
  }

  //! The grid is already sorted, so only its range is checked
  template <std::derived_from<sorter_t> G> constexpr ranged_grid_t& operator=(const G& g) {
    assert_ranged(g.vec());
    grid() = static_cast<const sorter_t&>(g);
    return *this;
  }

  //! The grid is already sorted, so only its range is checked
  template <std::derived_from<sorter_t> G> constexpr ranged_grid_t& operator=(G&& g) {
    assert_ranged(g.vec());
    grid() = static_cast<sorter_t&&>(g);
    return *this;
  }

  template <Size N, lagrange_interp::grid_transformer transform> [[nodiscard]] auto lag(Numeric x) const {
    return lagrange_interp::lag_t<N, transform>{*this, x, lag_sorting_t{}};
  }

  template <lagrange_interp::grid_transformer transform> [[nodiscard]] auto lag(Numeric x, Size N) const {
    return lagrange_interp::lag_t<-1, transform>{*this, x, N, lag_sorting_t{}};
  }

  template <Size N, lagrange_interp::grid_transformer transform>
  [[nodiscard]] auto lag(std::span<const Numeric> x,
                         Numeric                  extrapolation_limit = 0.5,
                         const char*              info                = "UNNAMED") const {
    return lagrange_interp::make_lags<N, transform>(*this, x, extrapolation_limit, info);
  }

  template <lagrange_interp::grid_transformer transform> [[nodiscard]] auto lag(std::span<const Numeric> x,
                                                                                Size                     N,
                                                                                Numeric     extrapolation_limit = 0.5,
                                                                                const char* info = "UNNAMED") const {
    return lagrange_interp::make_lags<transform>(*this, x, N, extrapolation_limit, info);
  }
};
}  // namespace matpack

//! Guaranteed [-180, 180)
using LonGrid = matpack::ranged_grid_t<-180, 180, matpack::lim::incl, matpack::lim::excl>;

//! Guaranteed [-90, 90]
using LatGrid = matpack::ranged_grid_t<-90, 90>;

//! Guaranteed [0, 180]
using ZenGrid = matpack::ranged_grid_t<0, 180>;

//! Guaranteed [0, 360)
using AziGrid = matpack::ranged_grid_t<0, 360, matpack::lim::incl, matpack::lim::excl>;

//! Altitude, unbounded
using AltGrid = matpack::ranged_grid_t<0, 0, matpack::lim::inf, matpack::lim::inf>;

//! Frequency, guaranteed (0, inf]
using FreqGrid = matpack::ranged_grid_t<0, 0, matpack::lim::excl, matpack::lim::inf>;

template <> struct xml_io_stream_name<LonGrid> {
  static constexpr std::string_view name = "LonGrid"sv;
};

template <> struct xml_io_stream_name<LatGrid> {
  static constexpr std::string_view name = "LatGrid"sv;
};

template <> struct xml_io_stream_name<ZenGrid> {
  static constexpr std::string_view name = "ZenGrid"sv;
};

template <> struct xml_io_stream_name<AziGrid> {
  static constexpr std::string_view name = "AziGrid"sv;
};

template <> struct xml_io_stream_name<AltGrid> {
  static constexpr std::string_view name = "AltGrid"sv;
};

template <> struct xml_io_stream_name<FreqGrid> {
  static constexpr std::string_view name = "FreqGrid"sv;
};

template <Index lower_limit, Index upper_limit, matpack::lim include_lower, matpack::lim include_upper>
struct xml_io_stream<matpack::ranged_grid_t<lower_limit, upper_limit, include_lower, include_upper>> {
  using T = matpack::ranged_grid_t<lower_limit, upper_limit, include_lower, include_upper>;

  static constexpr std::string_view type_name = xml_io_stream_name_v<T>;

  static void write(std::ostream& os, const T& x, bofstream* pbofs = nullptr, std::string_view name = ""sv) {
    xml_write_to_stream(os, x.vec(), pbofs, name);
  }

  static void read(std::istream& is, T& x, bifstream* pbifs = nullptr) {
    Vector data;
    xml_read_from_stream(is, data, pbifs);
    x = std::move(data);
  }
};

template <Index lower_limit, Index upper_limit, matpack::lim include_lower, matpack::lim include_upper>
struct std::formatter<matpack::ranged_grid_t<lower_limit, upper_limit, include_lower, include_upper>> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const matpack::ranged_grid_t<lower_limit, upper_limit, include_lower, include_upper>& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.vec());
  }
};

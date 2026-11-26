#pragma once

#include <enumsInterpolationExtrapolation.h>
#include <field_concepts.h>
#include <matpack.h>
#include <operators.h>

#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "subsurf_field.h"

//! A property of the spectral subsurface field, a type distinct from SubsurfacePropertyTag
struct SubsurfaceSpectralPropertyTag {
  String name;

  auto operator<=>(const SubsurfaceSpectralPropertyTag &x) const = default;
};

namespace std {
template <> struct hash<SubsurfaceSpectralPropertyTag> {
  static std::size_t operator()(const SubsurfaceSpectralPropertyTag &pp) noexcept {
    return std::hash<String>{}(pp.name);
  }
};
}  // namespace std

namespace Subsurface {
//! A function of frequency [Hz], altitude relative to the local surface [m], latitude [deg], and longitude [deg]
using SpectralFunctionalData = NumericQuaternaryOperator;
using SpectralFieldData      = std::variant<SpectralGeodeticField4, Numeric, SpectralFunctionalData>;

/** One parameter of the subsurface as a function of frequency and position
 *
 * The altitude is relative to the local surface, as for SubsurfaceField.
 *
 * Gridded data is interpolated linearly in frequency, altitude, latitude and
 * longitude, and frequency is the slowest dimension of its flat view.
 * Outside a grid, the extrapolation of that side of that dimension applies:
 * None throws, Zero gives zero, Nearest uses the edge of the grid, and Linear
 * extrapolates linearly.  With log_interpolation, the logarithm of the data is
 * interpolated instead, in all dimensions.
 *
 * See field_concepts.h for the flat weights.
 */
struct SpectralData {
  SpectralFieldData          data{SpectralFunctionalData{}};
  InterpolationExtrapolation frq_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation frq_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation alt_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation alt_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation lat_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lat_low{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_upp{InterpolationExtrapolation::None};
  InterpolationExtrapolation lon_low{InterpolationExtrapolation::None};
  //! Interpolate the logarithm of the data, which must then be positive
  bool log_interpolation{false};

  SpectralData();
  SpectralData(const SpectralData &);
  SpectralData(SpectralData &&) noexcept;
  SpectralData &operator=(const SpectralData &);
  SpectralData &operator=(SpectralData &&) noexcept;

  void adjust_interpolation_extrapolation();

  explicit SpectralData(Numeric x);
  SpectralData &operator=(Numeric x);

  explicit SpectralData(SpectralGeodeticField4 x);
  SpectralData &operator=(SpectralGeodeticField4 x);

  explicit SpectralData(SpectralFunctionalData x);
  SpectralData &operator=(SpectralFunctionalData x);

  [[nodiscard]] String data_type() const;

  //! The parameter at a frequency and position
  [[nodiscard]] Numeric at(const Numeric frq, const Numeric alt, const Numeric lat, const Numeric lon) const;

  //! The parameter at each frequency of freq_grid at a position
  [[nodiscard]] Vector at(const AscendingGrid &freq_grid,
                          const Numeric        alt,
                          const Numeric        lat,
                          const Numeric        lon) const;

  //! Flat weights at a frequency and position, into weights, which is cleared first to reuse its allocation; zero weights are left out
  void flat_weight(std::vector<std::pair<Index, Numeric>> &weights,
                   const Numeric                           frq,
                   const Numeric                           alt,
                   const Numeric                           lat,
                   const Numeric                           lon) const;

  [[nodiscard]] ConstVectorView flat_view() const;

  [[nodiscard]] VectorView flat_view();

  [[nodiscard]] bool ok() const;
};

//! The parameters of a spectral subsurface field at one position, each as values on a frequency grid
struct SpectralPoint {
  std::unordered_map<SubsurfaceSpectralPropertyTag, Vector> props;

  Vector       &operator[](const SubsurfaceSpectralPropertyTag &key);
  const Vector &operator[](const SubsurfaceSpectralPropertyTag &key) const;

  [[nodiscard]] bool contains(const SubsurfaceSpectralPropertyTag &key) const;

  [[nodiscard]] std::vector<SubsurfaceSpectralPropertyTag> keys() const;

  [[nodiscard]] Size size() const;
};

//! Parameters of the subsurface as functions of frequency and position, for the methods that use them
struct SpectralField {
  std::unordered_map<SubsurfaceSpectralPropertyTag, SpectralData> props;

  SpectralData       &operator[](const SubsurfaceSpectralPropertyTag &key);
  const SpectralData &operator[](const SubsurfaceSpectralPropertyTag &key) const;

  [[nodiscard]] bool contains(const SubsurfaceSpectralPropertyTag &key) const;

  [[nodiscard]] std::vector<SubsurfaceSpectralPropertyTag> keys() const;

  [[nodiscard]] Size size() const;

  //! The parameters at each frequency of freq_grid at a position
  [[nodiscard]] SpectralPoint at(const AscendingGrid &freq_grid,
                                 const Numeric        alt,
                                 const Numeric        lat,
                                 const Numeric        lon) const;
};
}  // namespace Subsurface

using SubsurfaceSpectralPoint = Subsurface::SpectralPoint;
using SubsurfaceSpectralData  = Subsurface::SpectralData;
using SubsurfaceSpectralField = Subsurface::SpectralField;

template <> struct std::formatter<SubsurfaceSpectralPropertyTag> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SubsurfaceSpectralPropertyTag &v, FmtContext &ctx) const {
    const std::string_view quote = tags.quote();
    return tags.format(ctx, quote, v.name, quote);
  }
};

template <> struct std::formatter<SubsurfaceSpectralData> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SubsurfaceSpectralData &v, FmtContext &ctx) const {
    const std::string_view sep = tags.sep();
    tags.add_if_bracket(ctx, "["sv);
    tags.format(ctx, v.data, sep);
    tags.add_if_bracket(ctx, "["sv);
    tags.format(ctx,
                v.frq_upp,
                sep,
                v.frq_low,
                sep,
                v.alt_upp,
                sep,
                v.alt_low,
                sep,
                v.lat_upp,
                sep,
                v.lat_low,
                sep,
                v.lon_upp,
                sep,
                v.lon_low);
    tags.add_if_bracket(ctx, "]]"sv);
    return ctx.out();
  }
};

template <> struct std::formatter<SubsurfaceSpectralField> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SubsurfaceSpectralField &v, FmtContext &ctx) const {
    if (tags.short_str) return tags.format(ctx, R"("Prop": )"sv, v.props.size());
    return tags.format(ctx, R"("Prop": )"sv, v.props);
  }
};

template <> struct std::formatter<SubsurfaceSpectralPoint> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SubsurfaceSpectralPoint &v, FmtContext &ctx) const {
    if (tags.short_str) return tags.format(ctx, R"("Prop": )"sv, v.props.size());
    return tags.format(ctx, R"("Prop": )"sv, v.props);
  }
};

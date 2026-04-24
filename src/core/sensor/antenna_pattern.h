#pragma once

#include <enumsAntennaType.h>
#include <matpack.h>

namespace sensor {
//! A 2D angular antenna pattern on local zenith and azimuth offsets.
struct AntennaPattern;

//! 2D gridded field of antenna weights on local zenith and azimuth offsets.
using AntennaPatternGriddedField =
    matpack::gridded_data_t<Numeric, AscendingGrid, AscendingGrid>;

//! Concept for selecting antenna patterns.
template <typename T>
concept AntennaPatternSelection = std::derived_from<T, AntennaPattern>;

struct AntennaPattern {
  AntennaType type{AntennaType::PencilBeam};

  Numeric sigma_zenith{0.0};
  Numeric sigma_azimuth{0.0};
  bool lookup_uses_angular_grids{false};

  AscendingGrid lookup_zenith_grid{};
  AscendingGrid lookup_azimuth_grid{};
  Matrix lookup_response{};

  AntennaPattern() = default;

  static AntennaPattern pencil();
  static AntennaPattern gaussian(Numeric zenith_std, Numeric azimuth_std);
  static AntennaPattern gaussian_fwhm(Numeric zenith_fwhm,
                                      Numeric azimuth_fwhm);
  static AntennaPattern lookup(const AscendingGrid& zenith_grid,
                               const AscendingGrid& azimuth_grid,
                               const Matrix& response_lookup);
  static AntennaPattern lookup(const ZenGrid& zenith_grid,
                               const AziGrid& azimuth_grid,
                               const Matrix& response_lookup);

  void set_pencil_beam();
  void set_gaussian(Numeric zenith_std, Numeric azimuth_std);
  void set_gaussian_fwhm(Numeric zenith_fwhm, Numeric azimuth_fwhm);
  void set_lookup(const AscendingGrid& zenith_grid,
                  const AscendingGrid& azimuth_grid,
                  ConstMatrixView response_lookup);
  void set_lookup(const ZenGrid& zenith_grid,
                  const AziGrid& azimuth_grid,
                  ConstMatrixView response_lookup);

  [[nodiscard]] Numeric operator()(Numeric delta_zenith,
                                   Numeric delta_azimuth) const;

  [[nodiscard]] AntennaPatternGriddedField response(
      const AscendingGrid& zenith_grid,
      const AscendingGrid& azimuth_grid) const;
  [[nodiscard]] AntennaPatternGriddedField normalized_response(
      const AscendingGrid& zenith_grid,
      const AscendingGrid& azimuth_grid) const;

  [[nodiscard]] AntennaPatternGriddedField raw_sensor(
      const AscendingGrid& dzen_grid, const AscendingGrid& dazi_grid) const;
  [[nodiscard]] AntennaPatternGriddedField normalized_raw_sensor(
      const AscendingGrid& dzen_grid, const AscendingGrid& dazi_grid) const;
};
}  // namespace sensor

template <>
struct std::formatter<sensor::AntennaPattern> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::AntennaPattern& v,
                              FmtContext& ctx) const {
    return tags.format(
        ctx,
        "type="sv,
        v.type,
        tags.sep(),
        "sigma_zenith="sv,
        v.sigma_zenith,
        tags.sep(),
        "sigma_azimuth="sv,
        v.sigma_azimuth,
        tags.sep(),
        "lookup_shape="sv,
        std::array{v.lookup_response.nrows(), v.lookup_response.ncols()});
  }
};
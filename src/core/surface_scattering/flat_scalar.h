#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <enumsInterpolationExtrapolation.h>
#include <format_tags.h>
#include <matpack.h>
#include <surf.h>
#include <xml_io_stream.h>

#include "surface_scattering_properties.h"

namespace surface_scattering {

/** Flat-scalar (specular) surface scattering model.
 *
 * Mirrors the workspace method spectral_surf_reflFlatScalar: stores a scalar
 * reflectivity R(f) in [0, 1] (clamped) as a spectral field on an arbitrary,
 * sorted frequency grid (SortedGriddedField1).
 *
 *   BRDF_specular(I->I) = Muelmat{R}        (all angles)
 *   emissivity_specular(I) = Muelmat{1 - R} (all angles)
 *
 * The diffuse tensors stay zero.
 */
struct FlatScalarSurfaceScatterer {
  /// Scalar reflectivity spectrum on an arbitrary sorted frequency grid.
  /// The single grid dimension must be in Hz (ascending order).
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  SortedGriddedField1 reflectivity_spectrum{};

  /// Interpolation and extrapolation method for frequency.
  /// Controls how values outside the grid domain are handled.
  InterpolationExtrapolation interp_extrapolation{
      InterpolationExtrapolation::Nearest};

  FlatScalarSurfaceScatterer() = default;
  explicit FlatScalarSurfaceScatterer(SortedGriddedField1 spectrum_);

  FlatScalarSurfaceScatterer(const FlatScalarSurfaceScatterer&)            = default;
  FlatScalarSurfaceScatterer(FlatScalarSurfaceScatterer&&) noexcept        = default;
  FlatScalarSurfaceScatterer& operator=(const FlatScalarSurfaceScatterer&) = default;
  FlatScalarSurfaceScatterer& operator=(FlatScalarSurfaceScatterer&&) noexcept = default;

  [[nodiscard]] SurfaceScatteringModelProperties
  get_surface_scattering_model_properties(const SurfacePoint& surf_point,
                                          Numeric lat,
                                          Numeric lon,
                                          const Vector& f_grid,
                                          const Vector& za_inc_grid,
                                          const Vector& aa_inc_grid,
                                          const Vector& za_scat_grid,
                                          const Vector& aa_scat_grid) const;

  [[nodiscard]] const SortedGriddedField1& get_reflectivity_spectrum() const {
    return reflectivity_spectrum;
  }
  void set_reflectivity_spectrum(const SortedGriddedField1& s) {
    reflectivity_spectrum = s;
  }

  friend std::ostream& operator<<(std::ostream& os, const FlatScalarSurfaceScatterer& s);
};

/** Flat-scalar surface scattering model with spatially-varying reflectivity.
 *
 * Extends :class:`FlatScalarSurfaceScatterer` by storing the reflectivity as
 * a 3-dimensional field over (latitude [deg], longitude [deg], frequency [Hz])
 * using a SortedGriddedField3.  At runtime the field is bilinearly
 * interpolated in the geographic dimensions and linearly interpolated onto the
 * simulation f_grid.
 */
struct FlatScalarSurfaceScattererField {
  /// Scalar reflectivity on a sorted (lat [deg], lon [deg], freq [Hz]) grid.
  /// All three grid dimensions must be in ascending order.
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  SortedGriddedField3 reflectivity_field{};

  /// Interpolation and extrapolation method for frequency grid.
  InterpolationExtrapolation interp_extrapolation{
      InterpolationExtrapolation::Nearest};

  FlatScalarSurfaceScattererField() = default;
  explicit FlatScalarSurfaceScattererField(SortedGriddedField3 field_);

  FlatScalarSurfaceScattererField(const FlatScalarSurfaceScattererField&)            = default;
  FlatScalarSurfaceScattererField(FlatScalarSurfaceScattererField&&) noexcept        = default;
  FlatScalarSurfaceScattererField& operator=(const FlatScalarSurfaceScattererField&) = default;
  FlatScalarSurfaceScattererField& operator=(FlatScalarSurfaceScattererField&&) noexcept = default;

  [[nodiscard]] SurfaceScatteringModelProperties
  get_surface_scattering_model_properties(const SurfacePoint& surf_point,
                                          Numeric lat,
                                          Numeric lon,
                                          const Vector& f_grid,
                                          const Vector& za_inc_grid,
                                          const Vector& aa_inc_grid,
                                          const Vector& za_scat_grid,
                                          const Vector& aa_scat_grid) const;

  [[nodiscard]] const SortedGriddedField3& get_reflectivity_field() const {
    return reflectivity_field;
  }
  void set_reflectivity_field(const SortedGriddedField3& f) {
    reflectivity_field = f;
  }

  friend std::ostream& operator<<(std::ostream& os, const FlatScalarSurfaceScattererField& s);
};

}  // namespace surface_scattering

template <>
struct std::formatter<surface_scattering::FlatScalarSurfaceScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const surface_scattering::FlatScalarSurfaceScatterer& v,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "FlatScalarSurfaceScatterer"sv);
    }
    return tags.format(ctx,
                       "FlatScalarSurfaceScatterer"sv,
                       ": "sv,
                       v.reflectivity_spectrum,
                       ", "sv,
                       v.interp_extrapolation);
  }
};

template <>
struct std::formatter<surface_scattering::FlatScalarSurfaceScattererField> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const surface_scattering::FlatScalarSurfaceScattererField& v,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "FlatScalarSurfaceScattererField"sv);
    }
    return tags.format(ctx,
                       "FlatScalarSurfaceScattererField"sv,
                       ": "sv,
                       v.reflectivity_field,
                       ", "sv,
                       v.interp_extrapolation);
  }
};

template <>
struct xml_io_stream<surface_scattering::FlatScalarSurfaceScatterer> {
  static constexpr std::string_view type_name = "FlatScalarSurfaceScatterer";

  static void write(std::ostream& os,
                    const surface_scattering::FlatScalarSurfaceScatterer& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   surface_scattering::FlatScalarSurfaceScatterer& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<surface_scattering::FlatScalarSurfaceScattererField> {
  static constexpr std::string_view type_name = "FlatScalarSurfaceScattererField";

  static void write(std::ostream& os,
                    const surface_scattering::FlatScalarSurfaceScattererField& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   surface_scattering::FlatScalarSurfaceScattererField& x,
                   bifstream* pbifs = nullptr);
};

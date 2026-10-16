#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <enumsInterpolationExtrapolation.h>
#include <format_tags.h>
#include <matpack.h>
#include <surf.h>
#include <xml_io_stream.h>

#include "scattering_internal.h"
#include "surface_scattering_properties.h"

namespace surface_scattering {

/** Fresnel (specular) surface scattering model.
 *
 * Mirrors the workspace method spectral_surf_reflFlatRealFresnel:
 * stores a real refractive index of the surface as a spectral field on an
 * arbitrary, sorted frequency grid (SortedGriddedField1) plus the refractive
 * index n1 of the medium the radiation propagates in (default 1.0).
 *
 *   BRDF_specular(I->I) = fresnel_reflectance(Rv, Rh)(za_inc)
 *   emissivity_specular(I)   = I4 - fresnel_reflectance(Rv, Rh)(za)
 *
 * The reflectance depends only on the incidence zenith angle (in degrees);
 * the emissivity tensor has no incidence axis, so it is evaluated at the
 * scattering zenith angle, which (for the specular consumer) is the mirror
 * of the outgoing direction and equals the physical incidence angle.
 * The diffuse tensors stay zero.
 */
struct FresnelSurfaceScatterer {
  /// Surface (medium 2) refractive index spectrum on an arbitrary sorted
  /// frequency grid.  The single grid dimension must be in Hz (ascending).
  /// Values must be > 0 (checked when the properties are computed).
  SortedGriddedField1 refractive_index_spectrum{};

  /// Refractive index of medium 1 (the medium the radiation propagates in).
  /// Must be > 0.
  Numeric n1{1.0};

  /// Interpolation and extrapolation method for frequency.
  /// Linear: unlimited linear extrapolation beyond the stored grid.
  /// Nearest: values outside the stored grid evaluate to the edge value.
  /// None: frequencies outside the stored grid are a user error.
  /// Zero: values outside the stored grid are 0 (a user error here, since the
  ///       refractive index must stay positive).
  InterpolationExtrapolation interp_extrapolation{
      InterpolationExtrapolation::Nearest};

  FresnelSurfaceScatterer() = default;
  explicit FresnelSurfaceScatterer(SortedGriddedField1 spectrum_, Numeric n1_ = 1.0);

  FresnelSurfaceScatterer(const FresnelSurfaceScatterer&)            = default;
  FresnelSurfaceScatterer(FresnelSurfaceScatterer&&) noexcept        = default;
  FresnelSurfaceScatterer& operator=(const FresnelSurfaceScatterer&) = default;
  FresnelSurfaceScatterer& operator=(FresnelSurfaceScatterer&&) noexcept = default;

  void check_n1() const;

  [[nodiscard]] SurfaceScatteringModelProperties
  get_surface_scattering_model_properties(const SurfacePoint& surf_point,
                                          Numeric lat,
                                          Numeric lon,
                                          const Vector& f_grid,
                                          const Vector& za_inc_grid,
                                          const Vector& aa_inc_grid,
                                          const Vector& za_scat_grid,
                                          const Vector& aa_scat_grid) const;

  [[nodiscard]] const SortedGriddedField1& get_refractive_index_spectrum() const {
    return refractive_index_spectrum;
  }
  void set_refractive_index_spectrum(const SortedGriddedField1& s) {
    refractive_index_spectrum = s;
  }

  friend std::ostream& operator<<(std::ostream& os, const FresnelSurfaceScatterer& s);
};

/** Fresnel surface scattering model with spatially-varying refractive index.
 *
 * Extends :class:`FresnelSurfaceScatterer` by storing the surface refractive
 * index as a 3-dimensional field over (latitude [deg], longitude [deg],
 * frequency [Hz]) using a SortedGriddedField3.  At runtime the field is
 * bilinearly interpolated in the geographic dimensions and linearly
 * interpolated onto the simulation f_grid.
 */
struct FresnelSurfaceScattererField {
  /// Surface refractive index on a sorted (lat [deg], lon [deg], freq [Hz])
  /// grid.  All three grid dimensions must be in ascending order.
  /// Values must be > 0 (checked when the properties are computed).
  SortedGriddedField3 refractive_index_field{};

  /// Refractive index of medium 1 (the medium the radiation propagates in).
  /// Must be > 0.
  Numeric n1{1.0};

  /// Interpolation and extrapolation method for frequency grid.
  /// Linear: unlimited linear extrapolation; Nearest: edge-value clamp;
  /// None: user error outside the stored grid; Zero: 0 outside (user error
  /// here, the refractive index must stay positive).
  InterpolationExtrapolation interp_extrapolation{
      InterpolationExtrapolation::Nearest};

  FresnelSurfaceScattererField() = default;
  explicit FresnelSurfaceScattererField(SortedGriddedField3 field_, Numeric n1_ = 1.0);

  void check_n1() const;

  FresnelSurfaceScattererField(const FresnelSurfaceScattererField&)            = default;
  FresnelSurfaceScattererField(FresnelSurfaceScattererField&&) noexcept        = default;
  FresnelSurfaceScattererField& operator=(const FresnelSurfaceScattererField&) = default;
  FresnelSurfaceScattererField& operator=(FresnelSurfaceScattererField&&) noexcept = default;

  [[nodiscard]] SurfaceScatteringModelProperties
  get_surface_scattering_model_properties(const SurfacePoint& surf_point,
                                          Numeric lat,
                                          Numeric lon,
                                          const Vector& f_grid,
                                          const Vector& za_inc_grid,
                                          const Vector& aa_inc_grid,
                                          const Vector& za_scat_grid,
                                          const Vector& aa_scat_grid) const;

  [[nodiscard]] const SortedGriddedField3& get_refractive_index_field() const {
    return refractive_index_field;
  }
  void set_refractive_index_field(const SortedGriddedField3& f) {
    validate_longitude_grid(f.grid<1>());
    refractive_index_field = f;
  }

  friend std::ostream& operator<<(std::ostream& os, const FresnelSurfaceScattererField& s);
};

}  // namespace surface_scattering

template <>
struct std::formatter<surface_scattering::FresnelSurfaceScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const surface_scattering::FresnelSurfaceScatterer& v,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "FresnelSurfaceScatterer"sv);
    }
    return tags.format(ctx,
                       "FresnelSurfaceScatterer"sv,
                       ": "sv,
                       v.n1,
                       ", "sv,
                       v.refractive_index_spectrum,
                       ", "sv,
                       v.interp_extrapolation);
  }
};

template <>
struct std::formatter<surface_scattering::FresnelSurfaceScattererField> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const surface_scattering::FresnelSurfaceScattererField& v,
                              FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "FresnelSurfaceScattererField"sv);
    }
    return tags.format(ctx,
                       "FresnelSurfaceScattererField"sv,
                       ": "sv,
                       v.n1,
                       ", "sv,
                       v.refractive_index_field,
                       ", "sv,
                       v.interp_extrapolation);
  }
};

template <>
struct xml_io_stream<surface_scattering::FresnelSurfaceScatterer> {
  static constexpr std::string_view type_name = "FresnelSurfaceScatterer";

  static void write(std::ostream& os,
                    const surface_scattering::FresnelSurfaceScatterer& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   surface_scattering::FresnelSurfaceScatterer& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<surface_scattering::FresnelSurfaceScattererField> {
  static constexpr std::string_view type_name = "FresnelSurfaceScattererField";

  static void write(std::ostream& os,
                    const surface_scattering::FresnelSurfaceScattererField& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   surface_scattering::FresnelSurfaceScattererField& x,
                   bifstream* pbifs = nullptr);
};

#pragma once

#include <arts_constants.h>
#include <configtypes.h>
#include <debug.h>
#include <format_tags.h>
#include <matpack.h>
#include <surf.h>
#include <xml_io_stream.h>

#include "surface_scattering_properties.h"

namespace surface_scattering {

/** Lambertian surface scattering model.
 *
 * Implements a Lambertian (isotropic) BRDF.  The reflectivity is stored as a
 * spectral field on an arbitrary, sorted frequency grid
 * (SortedGriddedField1).  At runtime the field is linearly interpolated onto
 * the simulation's f_grid, decoupling the stored spectral resolution from the
 * simulation grid.
 *
 *   BRDF(I->I) = r(f) / pi
 *   emissivity(I) = 1 - r(f)
 *
 * where r(f) is the reflectivity interpolated to frequency f.
 *
 * The SurfacePropertyTag names the surface property this model represents,
 * providing a semantic key for future lookup from SurfacePoint.
 */
struct LambertianSurfaceScatterer {
  /// Reflectivity spectrum on an arbitrary sorted frequency grid.
  /// The single grid dimension must be in Hz (ascending order).
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  SortedGriddedField1 reflectivity_spectrum{};

  LambertianSurfaceScatterer() = default;
  LambertianSurfaceScatterer(SortedGriddedField1 spectrum_);

  LambertianSurfaceScatterer(const LambertianSurfaceScatterer&)            = default;
  LambertianSurfaceScatterer(LambertianSurfaceScatterer&&) noexcept        = default;
  LambertianSurfaceScatterer& operator=(const LambertianSurfaceScatterer&) = default;
  LambertianSurfaceScatterer& operator=(LambertianSurfaceScatterer&&) noexcept = default;

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

  friend std::ostream& operator<<(std::ostream& os,
                                   const LambertianSurfaceScatterer& s);
};

/** Lambertian surface scattering model with spatially-varying reflectivity.
 *
 * Extends :class:`LambertianSurfaceScatterer` by storing the reflectivity as
 * a 3-dimensional field over (latitude [deg], longitude [deg], frequency [Hz])
 * using a SortedGriddedField3.  At runtime the field is bilinearly interpolated
 * in the geographic dimensions and linearly interpolated onto the simulation
 * f_grid.
 *
 *   BRDF(I->I) = r(lat, lon, f) / pi
 *   emissivity(I) = 1 - r(lat, lon, f)
 *
 * Values are clamped to [0, 1] after interpolation.
 * Extrapolation beyond the stored grids is permitted (values are still clamped).
 */
struct LambertianSurfaceScattererField {
  /// Reflectivity on a sorted (lat [deg], lon [deg], freq [Hz]) grid.
  /// All three grid dimensions must be in ascending order.
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  SortedGriddedField3 reflectivity_field{};

  LambertianSurfaceScattererField() = default;
  LambertianSurfaceScattererField(SortedGriddedField3 field_);

  LambertianSurfaceScattererField(const LambertianSurfaceScattererField&)            = default;
  LambertianSurfaceScattererField(LambertianSurfaceScattererField&&) noexcept        = default;
  LambertianSurfaceScattererField& operator=(const LambertianSurfaceScattererField&) = default;
  LambertianSurfaceScattererField& operator=(LambertianSurfaceScattererField&&) noexcept = default;

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

  friend std::ostream& operator<<(std::ostream& os,
                                   const LambertianSurfaceScattererField& s);
};

}  // namespace surface_scattering

template <>
struct std::formatter<surface_scattering::LambertianSurfaceScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(
      const surface_scattering::LambertianSurfaceScatterer& v,
      FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "LambertianSurfaceScatterer"sv);
    }
    return tags.format(ctx,
                       "LambertianSurfaceScatterer"sv,
                       ": "sv,
                       v.reflectivity_spectrum);
  }
};

template <>
struct std::formatter<surface_scattering::LambertianSurfaceScattererField> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(
      const surface_scattering::LambertianSurfaceScattererField& v,
      FmtContext& ctx) const {
    if (tags.names) {
      return tags.format(ctx, "LambertianSurfaceScattererField"sv);
    }
    return tags.format(ctx,
                       "LambertianSurfaceScattererField"sv,
                       ": "sv,
                       v.reflectivity_field);
  }
};

template <>
struct xml_io_stream<surface_scattering::LambertianSurfaceScatterer> {
  static constexpr std::string_view type_name = "LambertianSurfaceScatterer";

  static void write(std::ostream& os,
                    const surface_scattering::LambertianSurfaceScatterer& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   surface_scattering::LambertianSurfaceScatterer& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<surface_scattering::LambertianSurfaceScattererField> {
  static constexpr std::string_view type_name = "LambertianSurfaceScattererField";

  static void write(std::ostream& os,
                    const surface_scattering::LambertianSurfaceScattererField& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   surface_scattering::LambertianSurfaceScattererField& x,
                   bifstream* pbifs = nullptr);
};

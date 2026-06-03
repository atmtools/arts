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
  /// Tag identifying the surface property (e.g., "albedo")
  SurfacePropertyTag reflectivity_tag{};
  /// Reflectivity spectrum on an arbitrary sorted frequency grid.
  /// The single grid dimension must be in Hz (ascending order).
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  SortedGriddedField1 reflectivity_spectrum{};

  LambertianSurfaceScatterer() = default;
  LambertianSurfaceScatterer(SurfacePropertyTag tag,
                             SortedGriddedField1 spectrum_);

  LambertianSurfaceScatterer(const LambertianSurfaceScatterer&)            = default;
  LambertianSurfaceScatterer(LambertianSurfaceScatterer&&) noexcept        = default;
  LambertianSurfaceScatterer& operator=(const LambertianSurfaceScatterer&) = default;
  LambertianSurfaceScatterer& operator=(LambertianSurfaceScatterer&&) noexcept = default;

  [[nodiscard]] SurfaceScatteringModelProperties
  get_surface_scattering_model_properties(const SurfacePoint& surf_point,
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
                       v.reflectivity_tag.name,
                       ": "sv,
                       v.reflectivity_spectrum);
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

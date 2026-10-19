#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <surf.h>
#include <xml_io_stream.h>

#include "surface_scattering_properties.h"
#include "surface_spectral_field.h"

namespace surface_scattering {

/** Lambertian surface scattering model.
 *
 * Implements a Lambertian (isotropic) BRDF.  The reflectivity is stored as
 * surface spectral data over frequency, latitude, and longitude.  At runtime it
 * is interpolated onto the simulation's f_grid at the surface position,
 * decoupling the stored spectral resolution from the simulation grid.
 *
 *   BRDF(I->I) = r(f, lat, lon)
 *   emissivity(I) = 1 - r(f, lat, lon)
 *
 * Values are clamped to [0, 1] after interpolation.
 *
 * Units convention: the BRDF is stored as the *dimensionless scattering
 * kernel* r(f), not as a physical BRDF in [1/sr] (the physical Lambertian
 * BRDF is r/pi).  Consistent use therefore requires the angular quadrature
 * weights of *spectral_radSurfaceScatteringDiffuse* to carry the projected
 * solid angle factor cos(za_in) dOmega -- with pure sin(za) dza daa weights a
 * uniform incoming radiance field would reflect with effective reflectance
 * 2*pi*r.  *spectral_radSurfaceScatteringDiffuseDirect* applies the
 * projected-area factor cos(theta_inc) at the actual surface normal
 * internally, so a single sun beam scatters as r * cos(theta_inc) * I_beam.
 *
 * The SurfacePropertyTag names the surface property this model represents,
 * providing a semantic key for future lookup from SurfacePoint.
 */
struct LambertianSurfaceScatterer {
  /// Reflectivity over frequency [Hz], latitude [deg], and longitude [deg].
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  Data reflectivity{};

  LambertianSurfaceScatterer() = default;
  explicit LambertianSurfaceScatterer(Data reflectivity_);

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
                       "LambertianSurfaceScatterer"sv,
                       ": "sv,
                       v.reflectivity);
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

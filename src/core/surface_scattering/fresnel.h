#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <surf.h>
#include <xml_io_stream.h>

#include "surface_scattering_properties.h"
#include "surface_spectral_field.h"

namespace surface_scattering {

/** Fresnel (specular) surface scattering model.
 *
 * Mirrors the workspace method spectral_surf_reflFlatRealFresnel:
 * stores a real refractive index of the surface as surface spectral data over
 * frequency, latitude, and longitude, plus the refractive index n1 of the
 * medium the radiation propagates in (default 1.0).
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
  /// Surface (medium 2) refractive index over frequency [Hz], latitude [deg],
  /// and longitude [deg].  Values must be > 0 (checked when the properties are
  /// computed).
  Data refractive_index{};

  /// Refractive index of medium 1 (the medium the radiation propagates in).
  /// Must be > 0.
  Numeric n1{1.0};

  FresnelSurfaceScatterer() = default;
  explicit FresnelSurfaceScatterer(Data refractive_index_, Numeric n1_ = 1.0);

  FresnelSurfaceScatterer(const FresnelSurfaceScatterer&)                = default;
  FresnelSurfaceScatterer(FresnelSurfaceScatterer&&) noexcept            = default;
  FresnelSurfaceScatterer& operator=(const FresnelSurfaceScatterer&)     = default;
  FresnelSurfaceScatterer& operator=(FresnelSurfaceScatterer&&) noexcept = default;

  void check_n1() const;

  [[nodiscard]] SurfaceScatteringModelProperties get_surface_scattering_model_properties(
      const SurfacePoint& surf_point,
      Numeric             lat,
      Numeric             lon,
      const Vector&       f_grid,
      const Vector&       za_inc_grid,
      const Vector&       aa_inc_grid,
      const Vector&       za_scat_grid,
      const Vector&       aa_scat_grid) const;

  friend std::ostream& operator<<(std::ostream& os, const FresnelSurfaceScatterer& s);
};

}  // namespace surface_scattering

template <> struct std::formatter<surface_scattering::FresnelSurfaceScatterer> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const surface_scattering::FresnelSurfaceScatterer& v, FmtContext& ctx) const {
    if (tags.names) { return tags.format(ctx, "FresnelSurfaceScatterer"sv); }
    return tags.format(ctx, "FresnelSurfaceScatterer"sv, ": "sv, v.n1, ", "sv, v.refractive_index);
  }
};

template <> struct xml_io_stream<surface_scattering::FresnelSurfaceScatterer> {
  static constexpr std::string_view type_name = "FresnelSurfaceScatterer";

  static void write(std::ostream&                                      os,
                    const surface_scattering::FresnelSurfaceScatterer& x,
                    bofstream*                                         pbofs = nullptr,
                    std::string_view                                   name  = ""sv);

  static void read(std::istream& is, surface_scattering::FresnelSurfaceScatterer& x, bifstream* pbifs = nullptr);
};

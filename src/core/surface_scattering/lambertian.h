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
 * Implements a Lambertian (isotropic) BRDF:
 *   BRDF(I->I) = reflectivity[f] / pi
 *   emissivity(I) = 1 - reflectivity[f]
 *
 * The reflectivity vector is carried directly in the struct. The
 * SurfacePropertyTag names the surface property this model represents,
 * providing a semantic key for future lookup from SurfacePoint.
 */
struct LambertianSurfaceScatterer {
  /// Tag identifying the surface property (e.g., "albedo")
  SurfacePropertyTag reflectivity_tag{};
  /// Reflectivity values over f_grid (must match the f_grid passed to get_*)
  Vector reflectivity{};

  LambertianSurfaceScatterer() = default;
  LambertianSurfaceScatterer(SurfacePropertyTag tag, Vector reflectivity_);

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

  [[nodiscard]] const Vector& get_reflectivity() const { return reflectivity; }
  void set_reflectivity(const Vector& r) { reflectivity = r; }

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
    return tags.format(ctx, v.reflectivity_tag.name, ": "sv, v.reflectivity);
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


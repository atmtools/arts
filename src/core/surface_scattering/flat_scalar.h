#pragma once

#include <format_tags.h>
#include <matpack.h>
#include <surf.h>
#include <xml_io_stream.h>

#include "surface_scattering_properties.h"
#include "surface_spectral_field.h"

namespace surface_scattering {

/** Flat-scalar (specular) surface scattering model.
 *
 * Mirrors the workspace method spectral_surf_reflFlatScalar: stores a scalar
 * reflectivity R(f, lat, lon) in [0, 1] (clamped) as surface spectral data.
 *
 *   BRDF_specular(I->I) = Muelmat{R}        (all angles)
 *   emissivity_specular(I) = Muelmat{1 - R} (all angles)
 *
 * The diffuse tensors stay zero.
 */
struct FlatScalarSurfaceScatterer {
  /// Scalar reflectivity over frequency [Hz], latitude [deg], and longitude [deg].
  /// Values are expected in [0, 1]; out-of-range values are clamped.
  Data reflectivity{};

  FlatScalarSurfaceScatterer() = default;
  explicit FlatScalarSurfaceScatterer(Data reflectivity_);

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

  friend std::ostream& operator<<(std::ostream& os, const FlatScalarSurfaceScatterer& s);
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
    return tags.format(ctx, "FlatScalarSurfaceScatterer"sv, ": "sv, v.reflectivity);
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

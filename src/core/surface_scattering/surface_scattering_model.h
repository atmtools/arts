#pragma once

#include <configtypes.h>
#include <debug.h>
#include <format_tags.h>
#include <matpack.h>
#include <xml_io_base.h>
#include <xml_io_stream.h>
#include <xml_io_stream_variant.h>

#include <cstdint>
#include <format>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>

#include "flat_scalar.h"
#include "fresnel.h"
#include "lambertian.h"
#include "surface_scattering_properties.h"

namespace surface_scattering {

/// Variant type holding any concrete surface scattering model
using SurfaceScatteringModel =
    std::variant<LambertianSurfaceScatterer,
                 LambertianSurfaceScattererField,
                 FresnelSurfaceScatterer,
                 FresnelSurfaceScattererField,
                 FlatScalarSurfaceScatterer,
                 FlatScalarSurfaceScattererField>;

}  // namespace surface_scattering

using SurfaceScatteringModel = surface_scattering::SurfaceScatteringModel;

/// Pull LambertianSurfaceScatterer into the global namespace (mirrors
/// the pattern for HenyeyGreensteinScatterer in scattering_species.h)
using LambertianSurfaceScatterer = surface_scattering::LambertianSurfaceScatterer;
using LambertianSurfaceScattererField = surface_scattering::LambertianSurfaceScattererField;
using FresnelSurfaceScatterer = surface_scattering::FresnelSurfaceScatterer;
using FresnelSurfaceScattererField = surface_scattering::FresnelSurfaceScattererField;
using FlatScalarSurfaceScatterer = surface_scattering::FlatScalarSurfaceScatterer;
using FlatScalarSurfaceScattererField = surface_scattering::FlatScalarSurfaceScattererField;

/** Named map of surface scattering models.
 *
 * Models are stored by name for later individual lookup, and their bulk
 * surface scattering properties are accumulated over the models with
 * per-model weights derived from the *SurfacePropertyTag* masks of the
 * SurfacePoint (iteration follows std::map key order).  See Weighting for
 * how the masks become model weights.
 */
struct MapOfSurfaceScatteringModel {
  std::map<std::string, surface_scattering::SurfaceScatteringModel> models;

  /// Insert (or replace) a named model
  void add(const std::string& name,
           const surface_scattering::SurfaceScatteringModel& model);

  /// Accumulate bulk surface scattering properties from all stored models
  [[nodiscard]] surface_scattering::SurfaceScatteringModelProperties
  get_surface_scattering_model_properties(const SurfacePoint& surf_point,
                                         Numeric lat,
                                         Numeric lon,
                                         const Vector& f_grid,
                                         const Vector& za_inc_grid,
                                         const Vector& aa_inc_grid,
                                         const Vector& za_scat_grid,
                                         const Vector& aa_scat_grid) const;

  // Weighting options for combining multiple models.  The raw per-model
  // weights are the SurfacePoint mask values under the model's name key
  // (missing key -> 0); masks must be non-negative.
  enum class Weighting : std::uint8_t  {
    /// Winner takes all: the model(s) with the largest mask value get weight
    /// 1 (split equally on ties), all others 0.  All-zero masks give every
    /// model weight 1/N.
    Maximum,
    /// Mask-weighted average: weights are the raw mask values normalized to
    /// sum 1.  All-zero masks give all-zero weights (no scattering).
    Average,
  };
  Weighting weighting_option = Weighting::Maximum;

private:
  [[nodiscard]] Vector maximum_weighting(const SurfacePoint& surf_point) const;
  [[nodiscard]] Vector average_weighting(const SurfacePoint& surf_point) const;
  [[nodiscard]] Vector get_raw_weighting(const SurfacePoint& surf_point) const;


};

template <>
struct std::formatter<MapOfSurfaceScatteringModel> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const MapOfSurfaceScatteringModel& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    tags.add_if_bracket(ctx, "{");
    bool first = true;
    for (const auto& [name, model] : v.models) {
      if (!first) tags.format(ctx, sep);
      first = false;
      tags.format(ctx, name);
    }
    tags.add_if_bracket(ctx, "}");
    return ctx.out();
  }
};

//! XML type name for the SurfaceScatteringModel variant
template <>
struct xml_io_stream_name<SurfaceScatteringModel> {
  static constexpr std::string_view name = "SurfaceScatteringModel"sv;
};

//! XML type name for MapOfSurfaceScatteringModel
template <>
struct xml_io_stream_name<MapOfSurfaceScatteringModel> {
  static constexpr std::string_view name = "MapOfSurfaceScatteringModel"sv;
};

template <>
struct xml_io_stream<MapOfSurfaceScatteringModel> {
  static constexpr std::string_view type_name = "MapOfSurfaceScatteringModel";

  static void write(std::ostream& os,
                    const MapOfSurfaceScatteringModel& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   MapOfSurfaceScatteringModel& x,
                   bifstream* pbifs = nullptr);
};


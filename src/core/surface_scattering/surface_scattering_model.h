#pragma once

#include <configtypes.h>
#include <debug.h>
#include <format_tags.h>
#include <matpack.h>
#include <xml_io_base.h>
#include <xml_io_stream.h>
#include <xml_io_stream_variant.h>

#include <format>
#include <map>
#include <stdexcept>
#include <string>
#include <variant>

#include "lambertian.h"
#include "surface_scattering_properties.h"

namespace surface_scattering {

/// Variant type holding any concrete surface scattering model
using SurfaceScatteringModel = std::variant<LambertianSurfaceScatterer>;

}  // namespace surface_scattering

using SurfaceScatteringModel = surface_scattering::SurfaceScatteringModel;

/// Pull LambertianSurfaceScatterer into the global namespace (mirrors
/// the pattern for HenyeyGreensteinScatterer in scattering_species.h)
using LambertianSurfaceScatterer = surface_scattering::LambertianSurfaceScatterer;

/** Named map of surface scattering models.
 *
 * Models are stored by name for later individual lookup, and their bulk
 * surface scattering properties are accumulated in insertion order when
 * get_bulk_surface_scattering_properties() is called.
 * Mirrors ArrayOfScatteringSpecies but uses a named std::map instead of a
 * plain vector, matching the plan for MapOfSurfaceScatteringModel.
 */
struct MapOfSurfaceScatteringModel {
  std::map<std::string, surface_scattering::SurfaceScatteringModel> models;

  /// Insert (or replace) a named model
  void add(const std::string& name,
           const surface_scattering::SurfaceScatteringModel& model);

  /// Accumulate bulk surface scattering properties from all stored models
  [[nodiscard]] surface_scattering::BulkSurfaceScatteringProperties
  get_bulk_surface_scattering_properties(const SurfacePoint& surf_point,
                                         const Vector& f_grid,
                                         const Vector& za_inc_grid,
                                         const Vector& aa_inc_grid,
                                         const Vector& za_scat_grid,
                                         const Vector& aa_scat_grid) const;
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
    tags.add_if_bracket(ctx, '{');
    bool first = true;
    for (const auto& [name, model] : v.models) {
      if (!first) tags.format(ctx, sep);
      first = false;
      tags.format(ctx, name);
    }
    tags.add_if_bracket(ctx, '}');
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


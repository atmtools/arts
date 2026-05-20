#include "surface_scattering_model.h"

#include <xml_io_base.h>
#include <xml_io_stream_core.h>
#include <xml_io_stream_variant.h>

void MapOfSurfaceScatteringModel::add(
    const std::string& name,
    const surface_scattering::SurfaceScatteringModel& model) {
  models[name] = model;
}

surface_scattering::SurfaceScatteringModelProperties
MapOfSurfaceScatteringModel::get_surface_scattering_model_properties(
    const SurfacePoint& surf_point,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  if (models.empty()) {
    const Index nf  = f_grid.size();
    const Index nzs = za_scat_grid.size();
    return {std::nullopt, Tensor3(nf, nzs, 4, 0.0)};
  }

  const auto visitor =
      [&](const auto& model) -> surface_scattering::SurfaceScatteringModelProperties {
    if constexpr (requires {
                    model.get_surface_scattering_model_properties(
                        surf_point,
                        f_grid,
                        za_inc_grid,
                        aa_inc_grid,
                        za_scat_grid,
                        aa_scat_grid);
                  }) {
      return model.get_surface_scattering_model_properties(
          surf_point, f_grid, za_inc_grid, aa_inc_grid, za_scat_grid, aa_scat_grid);
    } else {
      throw std::runtime_error(std::format(
          "Method not implemented for surface scattering model:\n{:N}", model));
    }
    std::unreachable();
  };

  auto it  = models.begin();
  auto bsp = std::visit(visitor, it->second);
  for (++it; it != models.end(); ++it) {
    bsp += std::visit(visitor, it->second);
  }
  return bsp;
}

void xml_io_stream<MapOfSurfaceScatteringModel>::write(
    std::ostream& os,
    const MapOfSurfaceScatteringModel& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name, "nelem", x.models.size());
  tag.write_to_stream(os);

  for (const auto& [key, model] : x.models) {
    // Write the entry key as a String
    xml_write_to_stream(os, String{key}, pbofs);
    // Write the variant model
    xml_write_to_stream(os, model, pbofs);
  }

  tag.write_to_end_stream(os);
}

void xml_io_stream<MapOfSurfaceScatteringModel>::read(
    std::istream& is,
    MapOfSurfaceScatteringModel& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  Size nelem = 0;
  tag.get_attribute_value("nelem", nelem);
  x.models.clear();

  for (Size i = 0; i < nelem; ++i) {
    String key;
    xml_read_from_stream(is, key, pbifs);
    surface_scattering::SurfaceScatteringModel model{
        surface_scattering::LambertianSurfaceScatterer{}};
    xml_read_from_stream(is, model, pbifs);
    x.models[key] = std::move(model);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}


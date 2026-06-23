#include "surface_scattering_model.h"

#include <xml_io_base.h>
#include <xml_io_stream_core.h>
#include <xml_io_stream_variant.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include "rtepack.h"

void MapOfSurfaceScatteringModel::add(
    const std::string& name,
    const surface_scattering::SurfaceScatteringModel& model) {
  models[name] = model;
}

Vector MapOfSurfaceScatteringModel::get_raw_weighting(
    const SurfacePoint& surf_point) const {
  const Index N_models = models.size();
  Vector weights(N_models, 0.);
  //loop over MapOfSurfaceScatteringModel
  Index i = 0;
  for (const auto& [key, model] : models) {
    //now we have to check if in Surface point is a variable with the same name
    // as the key in the map
    if (surf_point.contains(SurfacePropertyTag{key})) {
      //if it is we have to get the value of this variable and compare it with the current maximum
      weights[i] = surf_point[SurfacePropertyTag{key}];
    } else {
      weights[i] = 0.;
    }
    i++;
  }
  return weights;
}

Vector MapOfSurfaceScatteringModel::maximum_weighting(
    const SurfacePoint& surf_point) const {
  Vector weights = get_raw_weighting(surf_point);
  // Since we now have the weights, we now set every weight except the maximum to zero and
  // the maximum to 1
  Numeric max_weight = *std::max_element(weights.begin(), weights.end());
  for (auto& w : weights) {
    if (w < max_weight) {
      w = 0.;
    } else {
      w = 1.;
    }
  }

  // Now we check if the sum of the weights is 1, if not we normalize the weights
  // This can happen if more then one value has the same value as the maximum
  Numeric sum_weights = std::accumulate(weights.begin(), weights.end(), 0.);
  if (sum_weights > 1) {
    for (auto& w : weights) {
      w /= sum_weights;
    }
  }

  return weights;
}

Vector MapOfSurfaceScatteringModel::average_weighting(
    const SurfacePoint& surf_point) const {
  Vector weights = get_raw_weighting(surf_point);
  // Since we now have the weights, we now set every weight to the average of the weights
  Numeric sum_weights = std::accumulate(weights.begin(), weights.end(), 0.);

  if (sum_weights > 0) {
    for (auto& w : weights) {
      w /= sum_weights;
    }
  }

  return weights;
}

surface_scattering::SurfaceScatteringModelProperties
MapOfSurfaceScatteringModel::get_surface_scattering_model_properties(
    const SurfacePoint& surf_point,
    Numeric lat,
    Numeric lon,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid
    ) const
{
  if (models.empty()) {
    const Index nf  = f_grid.size();
    const Index nzs = za_scat_grid.size();
    const Index nas = aa_scat_grid.size();
    const Index nzi = za_inc_grid.size();
    const Index nai = aa_inc_grid.size();
    return {.brdf_matrix=MuelmatTensor5(nf, nzi, nai, nzs, nas,0.0), .emissivity_vector=MuelmatTensor3(nf, nzs, nas, 0.0)};
  }

  const auto visitor = [&](const auto& model)
      -> surface_scattering::SurfaceScatteringModelProperties {
    if constexpr (requires {
                    model.get_surface_scattering_model_properties(surf_point,
                                                                  lat,
                                                                  lon,
                                                                  f_grid,
                                                                  za_inc_grid,
                                                                  aa_inc_grid,
                                                                  za_scat_grid,
                                                                  aa_scat_grid);
                  }) {
      return model.get_surface_scattering_model_properties(surf_point,
                                                           lat,
                                                           lon,
                                                           f_grid,
                                                           za_inc_grid,
                                                           aa_inc_grid,
                                                           za_scat_grid,
                                                           aa_scat_grid);
    } else {
      throw std::runtime_error(std::format(
          "Method not implemented for surface scattering model:\n{:N}", model));
    }
    std::unreachable();
  };

  // Now we need the weighting according to weighting_option
  Vector weights;

  std::cout << "weighting option: " << (weighting_option == Weighting::Maximum ? "Maximum" : "Average") << "\n";
  if (weighting_option == Weighting::Maximum) {
    weights = maximum_weighting(surf_point);
  } else if (weighting_option == Weighting::Average) {
    weights = average_weighting(surf_point);
  } else {
    throw std::runtime_error("Invalid weighting option");
  }

  //Now we have to sum if we have more than one model, but we have to weight the models according to the weights
  surface_scattering::SurfaceScatteringModelProperties bsp;
  Index i = 0;
  for (const auto& [key, model] : models) {
    surface_scattering::SurfaceScatteringModelProperties model_props =
        std::visit(visitor, model);

    model_props *= weights[i];
    bsp += model_props;
    i++;
  }
  return bsp;
}

surface_scattering::SurfaceScatteringModelProperties&
surface_scattering::SurfaceScatteringModelProperties::operator*=(Numeric scalar) {  
  brdf_matrix *= scalar;
  emissivity_vector *= scalar;
  return *this;
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


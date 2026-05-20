#include "lambertian.h"

#include <xml_io_base.h>

namespace surface_scattering {

LambertianSurfaceScatterer::LambertianSurfaceScatterer(SurfacePropertyTag tag,
                                                       Vector reflectivity_)
    : reflectivity_tag(std::move(tag)),
      reflectivity(std::move(reflectivity_)) {}

SurfaceScatteringModelProperties
LambertianSurfaceScatterer::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  ARTS_USER_ERROR_IF(
      reflectivity.size() != static_cast<Size>(f_grid.size()),
      "Reflectivity vector size ({}) must match f_grid size ({}).",
      reflectivity.size(),
      f_grid.size());

  const Index nf   = f_grid.size();
  const Index nzi  = za_inc_grid.size();
  const Index nai  = aa_inc_grid.size();
  const Index nzs  = za_scat_grid.size();
  const Index nas  = aa_scat_grid.size();

  Tensor7 brdf(nf, nzi, nai, nzs, nas, 4, 4, 0.0);
  Tensor3 emissivity(nf, nzs, 4, 0.0);

  for (Index f = 0; f < nf; ++f) {
    const Numeric r        = reflectivity[f];
    const Numeric brdf_val = r / Constant::pi;
    for (Index zi = 0; zi < nzi; ++zi)
      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as)
            brdf[f, zi, ai, zs, as, 0, 0] = brdf_val;
    for (Index zs = 0; zs < nzs; ++zs)
      emissivity[f, zs, 0] = 1.0 - r;
  }

  return SurfaceScatteringModelProperties{
      .brdf_matrix       = std::move(brdf),
      .emissivity_vector = std::move(emissivity),
  };
}

std::ostream& operator<<(std::ostream& os,
                         const LambertianSurfaceScatterer& s) {
  return os << "LambertianSurfaceScatterer(" << s.reflectivity_tag.name << ")";
}

SurfaceScatteringModelProperties& SurfaceScatteringModelProperties::operator+=(
    const SurfaceScatteringModelProperties& other) {
  if (brdf_matrix.has_value()) {
    ARTS_USER_ERROR_IF(
        !other.brdf_matrix.has_value(),
        "BRDF matrix missing in calculation of bulk surface scattering properties.");
    *brdf_matrix += *other.brdf_matrix;
  } else if (other.brdf_matrix.has_value()) {
    brdf_matrix = other.brdf_matrix;
  }
  emissivity_vector += other.emissivity_vector;
  return *this;
}

}  // namespace surface_scattering

void xml_io_stream<surface_scattering::LambertianSurfaceScatterer>::write(
    std::ostream& os,
    const surface_scattering::LambertianSurfaceScatterer& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.reflectivity_tag.name, pbofs);
  xml_write_to_stream(os, x.reflectivity, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<surface_scattering::LambertianSurfaceScatterer>::read(
    std::istream& is,
    surface_scattering::LambertianSurfaceScatterer& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity_tag.name, pbifs);
  xml_read_from_stream(is, x.reflectivity, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}


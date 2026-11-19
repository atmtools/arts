#include "flat_scalar.h"

#include <xml_io_base.h>

#include <algorithm>

#include "rtepack.h"

namespace surface_scattering {

namespace {

//! Compute the flat-scalar (specular) BRDF and complementary emissivity from a
//! per-frequency (possibly clamped) reflectivity vector.
void flat_scalar_tensors(const auto&     r_data,
                         Index           nf,
                         Index           nzi,
                         Index           nai,
                         Index           nzs,
                         Index           nas,
                         MuelmatTensor5& brdf_specular,
                         MuelmatTensor3& emissivity_specular) {
  Muelmat brdf;
  Muelmat emissivity;

  for (Index f = 0; f < nf; ++f) {
    const Numeric r = std::clamp(r_data[f], Numeric{0}, Numeric{1});
    brdf            = Muelmat{r};
    emissivity      = Muelmat{1.0 - r};

    for (Index zi = 0; zi < nzi; ++zi)
      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as) brdf_specular[f, zi, ai, zs, as] = brdf;
    for (Index zs = 0; zs < nzs; ++zs)
      for (Index as = 0; as < nas; ++as) emissivity_specular[f, zs, as] = emissivity;
  }
}

}  // namespace

FlatScalarSurfaceScatterer::FlatScalarSurfaceScatterer(Data reflectivity_) : reflectivity(std::move(reflectivity_)) {}

SurfaceScatteringModelProperties FlatScalarSurfaceScatterer::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric       lat,
    Numeric       lon,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  const Vector r_data = reflectivity.at(f_grid, lat, lon);

  MuelmatTensor5 brdf_specular(f_grid.size(),
                               za_inc_grid.size(),
                               aa_inc_grid.size(),
                               za_scat_grid.size(),
                               aa_scat_grid.size(),
                               rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_specular(f_grid.size(), za_scat_grid.size(), aa_scat_grid.size(), rtepack::muelmat{0.0});

  flat_scalar_tensors(r_data,
                      f_grid.size(),
                      za_inc_grid.size(),
                      aa_inc_grid.size(),
                      za_scat_grid.size(),
                      aa_scat_grid.size(),
                      brdf_specular,
                      emissivity_specular);

  // Purely specular model: the diffuse components stay missing (std::nullopt)
  return SurfaceScatteringModelProperties{
      .brdf_matrix_diffuse        = std::nullopt,
      .emissivity_vector_diffuse  = std::nullopt,
      .brdf_matrix_specular       = std::move(brdf_specular),
      .emissivity_vector_specular = std::move(emissivity_specular),
  };
}

std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const FlatScalarSurfaceScatterer& s) {
  return os << "FlatScalarSurfaceScatterer";
}

}  // namespace surface_scattering

using namespace surface_scattering;

void xml_io_stream<FlatScalarSurfaceScatterer>::write(std::ostream&                     os,
                                                      const FlatScalarSurfaceScatterer& x,
                                                      bofstream*                        pbofs,
                                                      std::string_view                  name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.reflectivity, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<FlatScalarSurfaceScatterer>::read(std::istream&               is,
                                                     FlatScalarSurfaceScatterer& x,
                                                     bifstream*                  pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

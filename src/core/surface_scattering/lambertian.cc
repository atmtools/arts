#include "lambertian.h"

#include <xml_io_base.h>

#include <algorithm>

#include "rtepack.h"

namespace surface_scattering {

namespace {

/** Compute Lambertian BRDF and emissivity from a per-frequency reflectivity vector.
 *
 * @param r_data  Reflectivity values (one per frequency); may be raw interpolated
 *                output (not yet clamped).
 * @param nf      Number of frequencies (== r_data.size() == f_grid.size()).
 * @param nzi     Number of incoming zenith angles.
 * @param nai     Number of incoming azimuth angles.
 * @param nzs     Number of scattering zenith angles.
 * @param nas     Number of scattering azimuth angles.
 */
SurfaceScatteringModelProperties lambertian_properties(
    const auto& r_data, Index nf, Index nzi, Index nai, Index nzs, Index nas) {
  MuelmatTensor5 brdf(nf, nzi, nai, nzs, nas, rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity(nf, nzs, nas, rtepack::muelmat{0.0});

  // Only the [0,0] element is non-zero; it is set per frequency below
  Muelmat isotropic_brdf{0.0};
  Muelmat isotropic_emissivity{0.0};

  for (Index f = 0; f < nf; ++f) {
    const Numeric r = std::clamp(r_data[f], Numeric{0}, Numeric{1});

    isotropic_brdf[0, 0]       = r;
    isotropic_emissivity[0, 0] = 1.0 - r;

    for (Index zi = 0; zi < nzi; ++zi)
      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as) brdf[f, zi, ai, zs, as] = isotropic_brdf;
    for (Index zs = 0; zs < nzs; ++zs)
      for (Index as = 0; as < nas; ++as) emissivity[f, zs, as] = isotropic_emissivity;
  }

  // Lambertian scattering has no specular component; those stay missing
  // (std::nullopt) so consumers read them as zero
  return SurfaceScatteringModelProperties{
      .brdf_matrix_diffuse        = std::move(brdf),
      .emissivity_vector_diffuse  = std::move(emissivity),
      .brdf_matrix_specular       = std::nullopt,
      .emissivity_vector_specular = std::nullopt,
  };
}

}  // namespace

LambertianSurfaceScatterer::LambertianSurfaceScatterer(Data reflectivity_) : reflectivity(std::move(reflectivity_)) {}

SurfaceScatteringModelProperties LambertianSurfaceScatterer::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric       lat,
    Numeric       lon,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  return lambertian_properties(reflectivity.at(f_grid, lat, lon),
                               f_grid.size(),
                               za_inc_grid.size(),
                               aa_inc_grid.size(),
                               za_scat_grid.size(),
                               aa_scat_grid.size());
}

std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const LambertianSurfaceScatterer& s) {
  return os << "LambertianSurfaceScatterer";
}

}  // namespace surface_scattering

void xml_io_stream<surface_scattering::LambertianSurfaceScatterer>::write(
    std::ostream&                                         os,
    const surface_scattering::LambertianSurfaceScatterer& x,
    bofstream*                                            pbofs,
    std::string_view                                      name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.reflectivity, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<surface_scattering::LambertianSurfaceScatterer>::read(
    std::istream& is, surface_scattering::LambertianSurfaceScatterer& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

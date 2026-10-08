#include "lambertian.h"
#include "rtepack.h"

#include <lagrange_interp.h>
#include <xml_io_base.h>
#include <cmath>

#include "scattering_internal.h"

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
    const auto& r_data,
    Index nf, Index nzi, Index nai, Index nzs, Index nas) {
  MuelmatTensor5 brdf(nf, nzi, nai, nzs, nas, rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity(nf, nzs, nas, rtepack::muelmat{0.0});

  // Lambertian scattering has no specular component; keep the tensors sized
  // (but zero) so consumers can index them uniformly
  MuelmatTensor5 brdf_specular(nf, nzi, nai, nzs, nas, rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_specular(nf, nzs, nas, rtepack::muelmat{0.0});

  Muelmat isotropic_brdf;
  Muelmat isotropic_emissivity;

  for (Index f = 0; f < nf; ++f) {
    const Numeric r        = std::clamp(r_data[f], Numeric{0}, Numeric{1});

    isotropic_brdf[0, 0] = r;
    isotropic_emissivity[0,0] = 1.0 - r;

    for (Index zi = 0; zi < nzi; ++zi)
      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as)
            brdf[f, zi, ai, zs, as] = isotropic_brdf;
    for (Index zs = 0; zs < nzs; ++zs)
      for (Index as = 0; as < nas; ++as)
        emissivity[f, zs, as] = isotropic_emissivity;
  }

  return SurfaceScatteringModelProperties{
      .brdf_matrix_diffuse       = std::move(brdf),
      .emissivity_vector_diffuse = std::move(emissivity),
      .brdf_matrix_specular      = std::move(brdf_specular),
      .emissivity_vector_specular = std::move(emissivity_specular),
  };
}

}  // namespace

LambertianSurfaceScatterer::LambertianSurfaceScatterer(SortedGriddedField1 spectrum_)
    : reflectivity_spectrum(std::move(spectrum_)) {}

SurfaceScatteringModelProperties
LambertianSurfaceScatterer::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric /*lat*/,
    Numeric /*lon*/,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  ARTS_USER_ERROR_IF(
      !reflectivity_spectrum.ok(),
      "reflectivity_spectrum is not valid (grid size does not match data size).");
  ARTS_USER_ERROR_IF(
      reflectivity_spectrum.grid<0>().empty(),
      "reflectivity_spectrum frequency grid is empty.");

  // Linearly interpolate the stored spectral reflectivity onto f_grid.
  // Extrapolation beyond the stored grid is permitted (extrapolation_limit =
  // max) so that simulations whose f_grid slightly exceeds the stored range
  // are handled gracefully; values are clamped to [0, 1] afterwards.
  using id = lagrange_interp::grid_identity;

  // Frequency extrapolation limit based on member setting
  const Numeric extrap_limit = frequency_extrap_limit(interp_extrapolation);

  const auto f_lag = lagrange_interp::make_lags<1, id>(
      reflectivity_spectrum.grid<0>(),
      f_grid,
      extrap_limit,
      "Reflectivity frequency grid");
  const auto r_data = lagrange_interp::reinterp(reflectivity_spectrum.data, f_lag);

  return lambertian_properties(r_data,
                               f_grid.size(),
                               za_inc_grid.size(),
                               aa_inc_grid.size(),
                               za_scat_grid.size(),
                               aa_scat_grid.size());
}

std::ostream& operator<<(std::ostream& os,
                         [[maybe_unused]] const LambertianSurfaceScatterer& s) {
  return os << "LambertianSurfaceScatterer";
}

LambertianSurfaceScattererField::LambertianSurfaceScattererField(
    SortedGriddedField3 field_)
    : reflectivity_field(std::move(field_)) {
  validate_longitude_grid(reflectivity_field.grid<1>());
}

SurfaceScatteringModelProperties
LambertianSurfaceScattererField::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric lat,
    Numeric lon,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  ARTS_USER_ERROR_IF(
      !reflectivity_field.ok(),
      "reflectivity_field is not valid (grid size does not match data size).");
  ARTS_USER_ERROR_IF(
      reflectivity_field.grid<0>().empty(),
      "reflectivity_field latitude grid is empty.");
  ARTS_USER_ERROR_IF(
      reflectivity_field.grid<1>().empty(),
      "reflectivity_field longitude grid is empty.");
  ARTS_USER_ERROR_IF(
      reflectivity_field.grid<2>().empty(),
      "reflectivity_field frequency grid is empty.");

  using id = lagrange_interp::grid_identity;

  // Single-point spatial lags with cyclic interpolation
  // Latitude: non-cyclic (use identity), as poles are not continuous
  const auto lat_lag = reflectivity_field.grid<0>().lag<1, id>(lat);
  
  // Frequency extrapolation limit based on member setting
  const Numeric extrap_limit = frequency_extrap_limit(interp_extrapolation);
  
  // Multi-point frequency lag using member-controlled extrapolation.
  const auto freq_lag = reflectivity_field.grid<2>().lag<1, id>(
      f_grid,
      extrap_limit,
      "Reflectivity frequency grid");

  // Interpolate: for each target frequency, fix spatial position and
  // linearly interpolate across (lat, lon, freq).
  const Index nf = f_grid.size();
  Vector r_data(nf);
  
  // Longitude: use [-180, 180] cycler (validated at construction time)
  const auto lon_lag = reflectivity_field.grid<1>().lag<1, lon_cycler>(lon);
  for (Index f = 0; f < nf; ++f) {
    r_data[f] = lagrange_interp::interp(
        reflectivity_field.data, lat_lag, lon_lag, freq_lag[f]);
  }

  return lambertian_properties(r_data,
                               nf,
                               za_inc_grid.size(),
                               aa_inc_grid.size(),
                               za_scat_grid.size(),
                               aa_scat_grid.size());
}

std::ostream& operator<<(std::ostream& os,
                         [[maybe_unused]] const LambertianSurfaceScattererField& s) {
  return os << "LambertianSurfaceScattererField";
}

SurfaceScatteringModelProperties& SurfaceScatteringModelProperties::operator+=(
    const SurfaceScatteringModelProperties& other) {
  brdf_matrix_diffuse += other.brdf_matrix_diffuse;
  emissivity_vector_diffuse += other.emissivity_vector_diffuse;
  brdf_matrix_specular += other.brdf_matrix_specular;
  emissivity_vector_specular += other.emissivity_vector_specular;
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

  xml_write_to_stream(os, x.reflectivity_spectrum, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<surface_scattering::LambertianSurfaceScatterer>::read(
    std::istream& is,
    surface_scattering::LambertianSurfaceScatterer& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity_spectrum, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<surface_scattering::LambertianSurfaceScattererField>::write(
    std::ostream& os,
    const surface_scattering::LambertianSurfaceScattererField& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.reflectivity_field, pbofs);
  xml_write_to_stream(os, x.interp_extrapolation, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<surface_scattering::LambertianSurfaceScattererField>::read(
    std::istream& is,
    surface_scattering::LambertianSurfaceScattererField& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity_field, pbifs);
  xml_read_from_stream(is, x.interp_extrapolation, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

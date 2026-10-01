#include "flat_scalar.h"
#include "rtepack.h"

#include <lagrange_interp.h>
#include <xml_io_base.h>
#include <cmath>
#include <limits>

#include "scattering_internal.h"

namespace surface_scattering {

namespace {

//! Compute the flat-scalar (specular) BRDF and complementary emissivity from a
//! per-frequency (possibly clamped) reflectivity vector.
void flat_scalar_tensors(
    const auto& r_data,
    Index nf,
    Index nzi,
    Index nai,
    Index nzs,
    Index nas,
    MuelmatTensor5& brdf_specular,
    MuelmatTensor3& emissivity_specular) {
  Muelmat brdf;
  Muelmat emissivity;

  for (Index f = 0; f < nf; ++f) {
    const Numeric r = std::clamp(r_data[f], Numeric{0}, Numeric{1});
    brdf     = Muelmat{r};
    emissivity = Muelmat{1.0 - r};

    for (Index zi = 0; zi < nzi; ++zi)
      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as)
            brdf_specular[f, zi, ai, zs, as] = brdf;
    for (Index zs = 0; zs < nzs; ++zs)
      for (Index as = 0; as < nas; ++as)
        emissivity_specular[f, zs, as] = emissivity;
  }
}

}  // namespace

FlatScalarSurfaceScatterer::FlatScalarSurfaceScatterer(SortedGriddedField1 spectrum_)
    : reflectivity_spectrum(std::move(spectrum_)) {}

SurfaceScatteringModelProperties
FlatScalarSurfaceScatterer::get_surface_scattering_model_properties(
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

  using id = lagrange_interp::grid_identity;
  const Numeric extrap_limit = frequency_extrap_limit(interp_extrapolation);

  const auto r_lag = lagrange_interp::make_lags<1, id>(
      reflectivity_spectrum.grid<0>(),
      f_grid,
      extrap_limit,
      "Reflectivity frequency grid");
  const auto r_data = lagrange_interp::reinterp(reflectivity_spectrum.data, r_lag);

  MuelmatTensor5 brdf_specular(f_grid.size(),
                               za_inc_grid.size(),
                               aa_inc_grid.size(),
                               za_scat_grid.size(),
                               aa_scat_grid.size(),
                               rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_specular(f_grid.size(),
                                     za_scat_grid.size(),
                                     aa_scat_grid.size(),
                                     rtepack::muelmat{0.0});
  MuelmatTensor5 brdf_diffuse(f_grid.size(),
                              za_inc_grid.size(),
                              aa_inc_grid.size(),
                              za_scat_grid.size(),
                              aa_scat_grid.size(),
                              rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_diffuse(f_grid.size(),
                                    za_scat_grid.size(),
                                    aa_scat_grid.size(),
                                    rtepack::muelmat{0.0});

  flat_scalar_tensors(
      r_data,
      f_grid.size(),
      za_inc_grid.size(),
      aa_inc_grid.size(),
      za_scat_grid.size(),
      aa_scat_grid.size(),
      brdf_specular,
      emissivity_specular);

  return SurfaceScatteringModelProperties{
      .brdf_matrix_diffuse        = std::move(brdf_diffuse),
      .emissivity_vector_diffuse  = std::move(emissivity_diffuse),
      .brdf_matrix_specular       = std::move(brdf_specular),
      .emissivity_vector_specular = std::move(emissivity_specular),
  };
}

std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const FlatScalarSurfaceScatterer& s) {
  return os << "FlatScalarSurfaceScatterer";
}

FlatScalarSurfaceScattererField::FlatScalarSurfaceScattererField(SortedGriddedField3 field_)
    : reflectivity_field(std::move(field_)) {
  validate_longitude_grid(reflectivity_field.grid<1>());
}

SurfaceScatteringModelProperties
FlatScalarSurfaceScattererField::get_surface_scattering_model_properties(
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
  const Numeric extrap_limit = frequency_extrap_limit(interp_extrapolation);

  const auto lat_lag  = reflectivity_field.grid<0>().lag<1, id>(lat);
  const auto lon_lag  = reflectivity_field.grid<1>().lag<1, lon_cycler>(lon);
  const auto freq_lag = reflectivity_field.grid<2>().lag<1, id>(
      f_grid, extrap_limit, "Reflectivity frequency grid");

  const Index nf = f_grid.size();
  Vector        r_data(nf);
  for (Index f = 0; f < nf; ++f) {
    r_data[f] = lagrange_interp::interp(reflectivity_field.data, lat_lag, lon_lag, freq_lag[f]);
  }

  MuelmatTensor5 brdf_specular(f_grid.size(),
                               za_inc_grid.size(),
                               aa_inc_grid.size(),
                               za_scat_grid.size(),
                               aa_scat_grid.size(),
                               rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_specular(f_grid.size(),
                                     za_scat_grid.size(),
                                     aa_scat_grid.size(),
                                     rtepack::muelmat{0.0});
  MuelmatTensor5 brdf_diffuse(f_grid.size(),
                              za_inc_grid.size(),
                              aa_inc_grid.size(),
                              za_scat_grid.size(),
                              aa_scat_grid.size(),
                              rtepack::muelmat{0.0});
  MuelmatTensor3 emissivity_diffuse(f_grid.size(),
                                    za_scat_grid.size(),
                                    aa_scat_grid.size(),
                                    rtepack::muelmat{0.0});

  flat_scalar_tensors(
      r_data,
      nf,
      za_inc_grid.size(),
      aa_inc_grid.size(),
      za_scat_grid.size(),
      aa_scat_grid.size(),
      brdf_specular,
      emissivity_specular);

  return SurfaceScatteringModelProperties{
      .brdf_matrix_diffuse        = std::move(brdf_diffuse),
      .emissivity_vector_diffuse  = std::move(emissivity_diffuse),
      .brdf_matrix_specular       = std::move(brdf_specular),
      .emissivity_vector_specular = std::move(emissivity_specular),
  };
}

std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const FlatScalarSurfaceScattererField& s) {
  return os << "FlatScalarSurfaceScattererField";
}

}  // namespace surface_scattering

using namespace surface_scattering;

void xml_io_stream<FlatScalarSurfaceScatterer>::write(
    std::ostream& os,
    const FlatScalarSurfaceScatterer& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.reflectivity_spectrum, pbofs);
  xml_write_to_stream(os, x.interp_extrapolation, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<FlatScalarSurfaceScatterer>::read(
    std::istream& is,
    FlatScalarSurfaceScatterer& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity_spectrum, pbifs);
  xml_read_from_stream(is, x.interp_extrapolation, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<FlatScalarSurfaceScattererField>::write(
    std::ostream& os,
    const FlatScalarSurfaceScattererField& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.reflectivity_field, pbofs);
  xml_write_to_stream(os, x.interp_extrapolation, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<FlatScalarSurfaceScattererField>::read(
    std::istream& is,
    FlatScalarSurfaceScattererField& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.reflectivity_field, pbifs);
  xml_read_from_stream(is, x.interp_extrapolation, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

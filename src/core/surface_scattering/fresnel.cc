#include "fresnel.h"
#include "rtepack.h"

#include <lagrange_interp.h>
#include <physics_funcs.h>
#include <xml_io_base.h>
#include <cmath>
#include <limits>

#include "scattering_internal.h"

namespace {

//! Compute the Fresnel (specular) BRDF and complementary emissivity from a
//! per-frequency real refractive-index pair (n1 constant, n2 varying).
//!
//! Uses the *pair* overload of fresnel() (physics_funcs.h) which carries the
//! total-internal-reflection guard (|n1 sinθ1 / n2| > 1 -> {1,1}); the
//! out-parameter overload lacks that guard.
//!
//! The BRDF matrix at a frequency depends only on the incoming zenith angle
//! (degrees); the emissivity tensor (no incidence axis) is the energy
//! complement I4 - M evaluated at each scattering zenith angle.
void fresnel_tensors(const Numeric n1,
                     const auto& n2_data,
                     const Vector& za_inc_grid,
                     const Vector& za_scat_grid,
                     Index nf,
                     Index nzi,
                     Index nai,
                     Index nzs,
                     Index nas,
                     MuelmatTensor5& brdf_specular,
                     MuelmatTensor3& emissivity_specular) {
  for (Index f = 0; f < nf; ++f) {
    const Numeric n2 = n2_data[f];

    for (Index zi = 0; zi < nzi; ++zi) {
      const auto [Rv, Rh] = fresnel(Complex{n1}, Complex{n2}, za_inc_grid[zi]);
      const Muelmat M     = rtepack::fresnel_reflectance(Rv, Rh);

      for (Index ai = 0; ai < nai; ++ai)
        for (Index zs = 0; zs < nzs; ++zs)
          for (Index as = 0; as < nas; ++as)
            brdf_specular[f, zi, ai, zs, as] = M;
    }

    for (Index zs = 0; zs < nzs; ++zs) {
      for (Index as = 0; as < nas; ++as) {
        const auto [Rv, Rh] = fresnel(Complex{n1}, Complex{n2}, za_scat_grid[zs]);
        const Muelmat M     = rtepack::fresnel_reflectance(Rv, Rh);
        Muelmat E{};
        for (Index r = 0; r < 4; ++r)
          for (Index c = 0; c < 4; ++c)
            E[r, c] = (r == c ? 1.0 : 0.0) - M[r, c];
        emissivity_specular[f, zs, as] = E;
      }
    }
  }
}

}  // namespace

using namespace surface_scattering;

void FresnelSurfaceScatterer::check_n1() const {
  ARTS_USER_ERROR_IF(
      not (n1 > 0.0),
      "FresnelSurfaceScatterer n1 must be positive, but is {}.", n1);
}

FresnelSurfaceScatterer::FresnelSurfaceScatterer(SortedGriddedField1 spectrum_, Numeric n1_)
    : refractive_index_spectrum(std::move(spectrum_)), n1{n1_} {
  check_n1();
}

SurfaceScatteringModelProperties
FresnelSurfaceScatterer::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric /*lat*/,
    Numeric /*lon*/,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  check_n1();
  ARTS_USER_ERROR_IF(
      !refractive_index_spectrum.ok(),
      "refractive_index_spectrum is not valid (grid size does not match data size).");
  ARTS_USER_ERROR_IF(
      refractive_index_spectrum.grid<0>().empty(),
      "refractive_index_spectrum frequency grid is empty.");

  using id = lagrange_interp::grid_identity;
  const Numeric extrap_limit = frequency_extrap_limit(interp_extrapolation);

  const auto n2_lag = lagrange_interp::make_lags<1, id>(
      refractive_index_spectrum.grid<0>(),
      f_grid,
      extrap_limit,
      "Refractive index frequency grid");
  const auto n2_data = lagrange_interp::reinterp(refractive_index_spectrum.data, n2_lag);

  const Size nf_size = f_grid.size();
  for (Index f = 0; f < Index(nf_size); ++f) {
    ARTS_USER_ERROR_IF(not (n2_data[f] > 0.0),
                       "Refractive index must be positive, but {} (f={}).",
                       n2_data[f],
                       f_grid[f]);
  }

  MuelmatTensor5 brdf_specular(nf_size,
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

  fresnel_tensors(
      n1,
      n2_data,
      za_inc_grid,
      za_scat_grid,
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

std::ostream& operator<<(std::ostream& os, [[maybe_unused]] const FresnelSurfaceScatterer& s) {
  return os << std::format("FresnelSurfaceScatterer n1={}", s.n1);
}

void FresnelSurfaceScattererField::check_n1() const {
  ARTS_USER_ERROR_IF(
      not (n1 > 0.0),
      "FresnelSurfaceScattererField n1 must be positive, but is {}.", n1);
}

FresnelSurfaceScattererField::FresnelSurfaceScattererField(SortedGriddedField3 field_, Numeric n1_)
    : refractive_index_field(std::move(field_)), n1{n1_} {
  validate_longitude_grid(refractive_index_field.grid<1>());
  check_n1();
}

SurfaceScatteringModelProperties
FresnelSurfaceScattererField::get_surface_scattering_model_properties(
    const SurfacePoint& /*surf_point*/,
    Numeric lat,
    Numeric lon,
    const Vector& f_grid,
    const Vector& za_inc_grid,
    const Vector& aa_inc_grid,
    const Vector& za_scat_grid,
    const Vector& aa_scat_grid) const {
  check_n1();
  ARTS_USER_ERROR_IF(
      !refractive_index_field.ok(),
      "refractive_index_field is not valid (grid size does not match data size).");
  ARTS_USER_ERROR_IF(
      refractive_index_field.grid<0>().empty(),
      "refractive_index_field latitude grid is empty.");
  ARTS_USER_ERROR_IF(
      refractive_index_field.grid<1>().empty(),
      "refractive_index_field longitude grid is empty.");
  ARTS_USER_ERROR_IF(
      refractive_index_field.grid<2>().empty(),
      "refractive_index_field frequency grid is empty.");

  using id = lagrange_interp::grid_identity;
  const Numeric extrap_limit = frequency_extrap_limit(interp_extrapolation);

  const auto lat_lag  = refractive_index_field.grid<0>().lag<1, id>(lat);
  const auto lon_lag  = refractive_index_field.grid<1>().lag<1, lon_cycler>(lon);
  const auto freq_lag = refractive_index_field.grid<2>().lag<1, id>(
      f_grid, extrap_limit, "Refractive index frequency grid");

  const Index nf = f_grid.size();
  Vector       n2_data(nf);
  for (Index f = 0; f < nf; ++f) {
    n2_data[f] = lagrange_interp::interp(refractive_index_field.data, lat_lag, lon_lag, freq_lag[f]);
  }
  for (Index f = 0; f < nf; ++f) {
    ARTS_USER_ERROR_IF(not (n2_data[f] > 0.0),
                       "Refractive index must be positive, but {} (f={}).",
                       n2_data[f],
                       f_grid[f]);
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

  fresnel_tensors(
      n1,
      n2_data,
      za_inc_grid,
      za_scat_grid,
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

std::ostream&
operator<<(std::ostream& os, [[maybe_unused]] const FresnelSurfaceScattererField& s) {
  return os << std::format("FresnelSurfaceScattererField n1={}", s.n1);
}

void xml_io_stream<FresnelSurfaceScatterer>::write(
    std::ostream& os,
    const FresnelSurfaceScatterer& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.n1, pbofs, "n1"sv);
  xml_write_to_stream(os, x.refractive_index_spectrum, pbofs);
  xml_write_to_stream(os, x.interp_extrapolation, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<FresnelSurfaceScatterer>::read(
    std::istream& is,
    FresnelSurfaceScatterer& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.n1, pbifs);
  xml_read_from_stream(is, x.refractive_index_spectrum, pbifs);
  xml_read_from_stream(is, x.interp_extrapolation, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<FresnelSurfaceScattererField>::write(
    std::ostream& os,
    const FresnelSurfaceScattererField& x,
    bofstream* pbofs,
    std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.n1, pbofs, "n1"sv);
  xml_write_to_stream(os, x.refractive_index_field, pbofs);
  xml_write_to_stream(os, x.interp_extrapolation, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<FresnelSurfaceScattererField>::read(
    std::istream& is,
    FresnelSurfaceScattererField& x,
    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.n1, pbifs);
  xml_read_from_stream(is, x.refractive_index_field, pbifs);
  xml_read_from_stream(is, x.interp_extrapolation, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

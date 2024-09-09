#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <format>
#include <numeric>
#include <vector>

#include "arts_constants.h"
#include "configtypes.h"
#include "debug.h"
#include "mh_checks.h"

////////////////////////////////////////////////////////////////////////
// Core Disort
////////////////////////////////////////////////////////////////////////

struct disort_sizes {
  Index nv;
  Index np;
  Index nquad;
  Index nleg;
  Index nfou;
  Index nsource;
  Index nbdrf;
};

disort_sizes disort_get_and_check_sizes(
    const Matrix& disort_optical_thicknesses,
    const Matrix& disort_single_scattering_albedo,
    const Matrix& disort_fractional_scattering,
    const Tensor3& disort_legendre_coefficients,
    const Tensor3& disort_source_polynomial,
    const Tensor3& disort_positive_boundary_condition,
    const Tensor3& disort_negative_boundary_condition,
    const MatrixOfDisortBDRF&
        disort_bidirectional_reflectance_distribution_functions,
    const Vector& disort_solar_zenith_angle,
    const Vector& disort_solar_azimuth_angle,
    const Vector& disort_solar_source) {
  ARTS_USER_ERROR_IF(not same_shape(disort_optical_thicknesses,
                                    disort_single_scattering_albedo,
                                    disort_fractional_scattering),
                     R"(Bad shapes: {:B,} != {:B,} != {:B,}

disort_optical_thicknesses.shape()      = {:B,},
disort_single_scattering_albedo.shape() = {:B,},
disort_fractional_scattering.shape()    = {:B,}
)",
                     disort_optical_thicknesses.shape(),
                     disort_single_scattering_albedo.shape(),
                     disort_fractional_scattering.shape(),
                     disort_optical_thicknesses.shape(),
                     disort_single_scattering_albedo.shape(),
                     disort_fractional_scattering.shape())

  const Index nv = disort_optical_thicknesses.nrows();
  const Index np = disort_optical_thicknesses.ncols();

  ARTS_USER_ERROR_IF(not same_shape(std::array{nv},
                                    disort_solar_zenith_angle,
                                    disort_solar_azimuth_angle,
                                    disort_solar_source),
                     R"(Bad shapes: [{}] != {:B,} != {:B,} != {:B,}

disort_solar_zenith_angle.shape()  = {:B,},
disort_solar_azimuth_angle.shape() = {:B,},
disort_solar_source.shape()        = {:B,}
)",
                     nv,
                     disort_solar_zenith_angle.shape(),
                     disort_solar_azimuth_angle.shape(),
                     disort_solar_source.shape(),
                     disort_solar_zenith_angle.shape(),
                     disort_solar_azimuth_angle.shape(),
                     disort_solar_source.shape())

  const auto nleg = disort_legendre_coefficients.ncols();
  ARTS_USER_ERROR_IF(disort_legendre_coefficients.npages() != nv or
                         disort_legendre_coefficients.nrows() != np,
                     R"(Bad shapes: [{}, {}, {}] != {:B,}

disort_legendre_coefficients.shape() = {:B,}
)",
                     nv,
                     np,
                     nleg,
                     disort_legendre_coefficients.shape(),
                     disort_legendre_coefficients.shape())

  const Index nfou       = disort_negative_boundary_condition.nrows();
  const Index nquad_half = disort_negative_boundary_condition.ncols();

  ARTS_USER_ERROR_IF(not same_shape(std::array{nv, nfou, nquad_half},
                                    disort_negative_boundary_condition,
                                    disort_positive_boundary_condition),
                     R"(Bad shapes: [{}, {}, {}] != {:B,} != {:B,}

disort_negative_boundary_condition.shape() = {:B,},
disort_positive_boundary_condition.shape() = {:B,}
)",
                     nv,
                     nfou,
                     nquad_half,
                     disort_negative_boundary_condition.shape(),
                     disort_positive_boundary_condition.shape(),
                     disort_negative_boundary_condition.shape(),
                     disort_positive_boundary_condition.shape())

  const Index nsource = disort_source_polynomial.ncols();

  ARTS_USER_ERROR_IF(disort_source_polynomial.npages() != nv or
                         disort_source_polynomial.nrows() != np,
                     R"(Bad shapes: [{}, {}, {}] != {:B,}

disort_source_polynomial.shape() = {:B,}
)",
                     nv,
                     np,
                     nsource,
                     disort_source_polynomial.shape(),
                     disort_source_polynomial.shape())

  const auto nbdrf =
      disort_bidirectional_reflectance_distribution_functions.ncols();

  ARTS_USER_ERROR_IF(
      disort_bidirectional_reflectance_distribution_functions.nrows() != nv,

      R"(Bad shapes: [{} {}] != {:B,}

disort_bidirectional_reflectance_distribution_functions.shape() = {:B,}
)",
      nv,
      nbdrf,
      disort_bidirectional_reflectance_distribution_functions.shape(),
      disort_bidirectional_reflectance_distribution_functions.shape())

  return {.nv      = nv,
          .np      = np,
          .nquad   = nquad_half * 2,
          .nleg    = nleg,
          .nfou    = nfou,
          .nsource = nsource,
          .nbdrf   = nbdrf};
}

void disort_spectral_radiance_fieldCalc(
    Tensor4& disort_spectral_radiance_field,
    Vector& disort_quadrature_angles,
    Vector& disort_quadrature_weights,
    const Matrix& disort_optical_thicknesses,
    const Matrix& disort_single_scattering_albedo,
    const Matrix& disort_fractional_scattering,
    const Tensor3& disort_legendre_coefficients,
    const Tensor3& disort_source_polynomial,
    const Tensor3& disort_positive_boundary_condition,
    const Tensor3& disort_negative_boundary_condition,
    const MatrixOfDisortBDRF&
        disort_bidirectional_reflectance_distribution_functions,
    const Vector& disort_solar_zenith_angle,
    const Vector& disort_solar_azimuth_angle,
    const Vector& disort_solar_source,
    const Vector& phis) {
  using Conversion::acosd, Conversion::cosd, Conversion::deg2rad;

  const auto sizes = disort_get_and_check_sizes(
      disort_optical_thicknesses,
      disort_single_scattering_albedo,
      disort_fractional_scattering,
      disort_legendre_coefficients,
      disort_source_polynomial,
      disort_positive_boundary_condition,
      disort_negative_boundary_condition,
      disort_bidirectional_reflectance_distribution_functions,
      disort_solar_zenith_angle,
      disort_solar_azimuth_angle,
      disort_solar_source);

  const Index nv      = sizes.nv;
  const Index np      = sizes.np;
  const Index nquad   = sizes.nquad;
  const Index nleg    = sizes.nleg;
  const Index nfou    = sizes.nfou;
  const Index nsource = sizes.nsource;
  const Index nbdrf   = sizes.nbdrf;

  disort_spectral_radiance_field.resize(nv, np, phis.size(), nquad);

  const Index nleg_reduced = nleg;
  disort::main_data dis(np, nquad, nleg, nfou, nsource, nleg_reduced, nbdrf);

  //! Supplementary outputs
  disort_quadrature_weights = dis.weights();
  disort_quadrature_angles.resize(nquad);
  std::transform(dis.mu().begin(),
                 dis.mu().end(),
                 disort_quadrature_angles.begin(),
                 [](const Numeric& mu) { return acosd(mu); });

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      for (Index i = 0; i < nbdrf; i++) {
        dis.brdf_modes()[i] =
            disort_bidirectional_reflectance_distribution_functions(iv, i);
      }

      dis.solar_zenith()        = cosd(disort_solar_zenith_angle[iv]);
      dis.beam_azimuth()        = deg2rad(disort_solar_azimuth_angle[iv]);
      dis.tau()                 = disort_optical_thicknesses[iv];
      dis.omega()               = disort_single_scattering_albedo[iv];
      dis.f()                   = disort_fractional_scattering[iv];
      dis.all_legendre_coeffs() = disort_legendre_coefficients[iv];
      dis.positive_boundary()   = disort_positive_boundary_condition[iv];
      dis.negative_boundary()   = disort_negative_boundary_condition[iv];
      dis.source_poly()         = disort_source_polynomial[iv];

      dis.update_all(disort_solar_source[iv]);

      dis.gridded_u(disort_spectral_radiance_field[iv], phis);
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in disort:\n{}", error);
}

void disort_spectral_flux_fieldCalc(
    Tensor3& disort_spectral_flux_field,
    const Matrix& disort_optical_thicknesses,
    const Matrix& disort_single_scattering_albedo,
    const Matrix& disort_fractional_scattering,
    const Tensor3& disort_legendre_coefficients,
    const Tensor3& disort_source_polynomial,
    const Tensor3& disort_positive_boundary_condition,
    const Tensor3& disort_negative_boundary_condition,
    const MatrixOfDisortBDRF&
        disort_bidirectional_reflectance_distribution_functions,
    const Vector& disort_solar_zenith_angle,
    const Vector& disort_solar_azimuth_angle,
    const Vector& disort_solar_source) {
  using Conversion::acosd, Conversion::cosd, Conversion::deg2rad;

  const auto sizes = disort_get_and_check_sizes(
      disort_optical_thicknesses,
      disort_single_scattering_albedo,
      disort_fractional_scattering,
      disort_legendre_coefficients,
      disort_source_polynomial,
      disort_positive_boundary_condition,
      disort_negative_boundary_condition,
      disort_bidirectional_reflectance_distribution_functions,
      disort_solar_zenith_angle,
      disort_solar_azimuth_angle,
      disort_solar_source);

  const Index nv      = sizes.nv;
  const Index np      = sizes.np;
  const Index nquad   = sizes.nquad;
  const Index nleg    = sizes.nleg;
  const Index nfou    = sizes.nfou;
  const Index nsource = sizes.nsource;
  const Index nbdrf   = sizes.nbdrf;

  disort_spectral_flux_field.resize(nv, 3, np);

  const Index nleg_reduced = nleg;
  disort::main_data dis(np, nquad, nleg, nfou, nsource, nleg_reduced, nbdrf);

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      for (Index i = 0; i < nbdrf; i++) {
        dis.brdf_modes()[i] =
            disort_bidirectional_reflectance_distribution_functions(iv, i);
      }

      dis.solar_zenith()        = cosd(disort_solar_zenith_angle[iv]);
      dis.beam_azimuth()        = deg2rad(disort_solar_azimuth_angle[iv]);
      dis.tau()                 = disort_optical_thicknesses[iv];
      dis.omega()               = disort_single_scattering_albedo[iv];
      dis.f()                   = disort_fractional_scattering[iv];
      dis.all_legendre_coeffs() = disort_legendre_coefficients[iv];
      dis.positive_boundary()   = disort_positive_boundary_condition[iv];
      dis.negative_boundary()   = disort_negative_boundary_condition[iv];
      dis.source_poly()         = disort_source_polynomial[iv];

      dis.update_all(disort_solar_source[iv]);

      dis.gridded_flux(disort_spectral_flux_field(iv, 0, joker),
                       disort_spectral_flux_field(iv, 1, joker),
                       disort_spectral_flux_field(iv, 2, joker));
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in disort:\n{}", error);
}

////////////////////////////////////////////////////////////////////////
// Integrate functions
////////////////////////////////////////////////////////////////////////

void spectral_radianceIntegrateDisort(
    StokvecVector& /*spectral_radiance*/,
    const Tensor4& /*disort_spectral_radiance_field*/,
    const Vector& /*disort_quadrature_angles*/,
    const Vector& /*disort_quadrature_weights*/) {
  ARTS_USER_ERROR("Not implemented")
}

void SpectralFluxDisort(Matrix& spectral_flux_field_up,
                        Matrix& spectral_flux_field_down,
                        const Tensor3& disort_spectral_flux_field) {
  ARTS_USER_ERROR_IF(disort_spectral_flux_field.nrows() != 3,
                     "Must have shape (*, 3, *), but got shape {:B,}",
                     disort_spectral_flux_field.shape())

  spectral_flux_field_up    = disort_spectral_flux_field(joker, 0, joker);
  spectral_flux_field_down  = disort_spectral_flux_field(joker, 1, joker);
  spectral_flux_field_down += disort_spectral_flux_field(joker, 2, joker);
}

////////////////////////////////////////////////////////////////////////
// Disort Agenda wrappers
////////////////////////////////////////////////////////////////////////

void disort_spectral_radiance_fieldFromAgenda(
    const Workspace& ws,
    Tensor4& disort_spectral_radiance_field,
    Vector& disort_quadrature_angles,
    Vector& disort_quadrature_weights,
    const Agenda& disort_settings_agenda,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension,
    const Index& disort_legendre_polynomial_dimension,
    const Vector& phis) {
  Matrix disort_optical_thicknesses;
  Matrix disort_single_scattering_albedo;
  Matrix disort_fractional_scattering;
  Tensor3 disort_legendre_coefficients;
  Tensor3 disort_source_polynomial;
  Tensor3 disort_positive_boundary_condition;
  Tensor3 disort_negative_boundary_condition;
  MatrixOfDisortBDRF disort_bidirectional_reflectance_distribution_functions;
  Vector disort_solar_zenith_angle;
  Vector disort_solar_azimuth_angle;
  Vector disort_solar_source;

  disort_settings_agendaExecute(
      ws,
      disort_optical_thicknesses,
      disort_single_scattering_albedo,
      disort_fractional_scattering,
      disort_legendre_coefficients,
      disort_source_polynomial,
      disort_positive_boundary_condition,
      disort_negative_boundary_condition,
      disort_bidirectional_reflectance_distribution_functions,
      disort_solar_zenith_angle,
      disort_solar_azimuth_angle,
      disort_solar_source,
      disort_quadrature_dimension,
      disort_fourier_mode_dimension,
      disort_legendre_polynomial_dimension,
      disort_settings_agenda);

  disort_spectral_radiance_fieldCalc(
      disort_spectral_radiance_field,
      disort_quadrature_angles,
      disort_quadrature_weights,
      disort_optical_thicknesses,
      disort_single_scattering_albedo,
      disort_fractional_scattering,
      disort_legendre_coefficients,
      disort_source_polynomial,
      disort_positive_boundary_condition,
      disort_negative_boundary_condition,
      disort_bidirectional_reflectance_distribution_functions,
      disort_solar_zenith_angle,
      disort_solar_azimuth_angle,
      disort_solar_source,
      phis);
}

void disort_spectral_flux_fieldFromAgenda(
    const Workspace& ws,
    Tensor3& disort_spectral_flux_field,
    const Agenda& disort_settings_agenda,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension,
    const Index& disort_legendre_polynomial_dimension) {
  Matrix disort_optical_thicknesses;
  Matrix disort_single_scattering_albedo;
  Matrix disort_fractional_scattering;
  Tensor3 disort_legendre_coefficients;
  Tensor3 disort_source_polynomial;
  Tensor3 disort_positive_boundary_condition;
  Tensor3 disort_negative_boundary_condition;
  MatrixOfDisortBDRF disort_bidirectional_reflectance_distribution_functions;
  Vector disort_solar_zenith_angle;
  Vector disort_solar_azimuth_angle;
  Vector disort_solar_source;

  disort_settings_agendaExecute(
      ws,
      disort_optical_thicknesses,
      disort_single_scattering_albedo,
      disort_fractional_scattering,
      disort_legendre_coefficients,
      disort_source_polynomial,
      disort_positive_boundary_condition,
      disort_negative_boundary_condition,
      disort_bidirectional_reflectance_distribution_functions,
      disort_solar_zenith_angle,
      disort_solar_azimuth_angle,
      disort_solar_source,
      disort_quadrature_dimension,
      disort_fourier_mode_dimension,
      disort_legendre_polynomial_dimension,
      disort_settings_agenda);

  disort_spectral_flux_fieldCalc(
      disort_spectral_flux_field,
      disort_optical_thicknesses,
      disort_single_scattering_albedo,
      disort_fractional_scattering,
      disort_legendre_coefficients,
      disort_source_polynomial,
      disort_positive_boundary_condition,
      disort_negative_boundary_condition,
      disort_bidirectional_reflectance_distribution_functions,
      disort_solar_zenith_angle,
      disort_solar_azimuth_angle,
      disort_solar_source);
}
#include <arts_conversions.h>
#include <arts_omp.h>
#include <disort.h>
#include <workspace.h>

#include <algorithm>
#include <numeric>
#include <vector>

#include "arts_constants.h"
#include "auto_wsm.h"
#include "debug.h"
#include "matpack_data.h"
#include "matpack_view.h"
#include "mh_checks.h"
#include "mystring.h"
#include "path_point.h"
#include "physics_funcs.h"
#include "sorted_grid.h"
#include "surf.h"

void ray_pathGeometricUplooking(ArrayOfPropagationPathPoint& ray_path,
                                const AtmField& atmospheric_field,
                                const SurfaceField& surface_field,
                                const Numeric& latitude,
                                const Numeric& longitude,
                                const Numeric& max_step) {
  ray_pathGeometric(
      ray_path,
      atmospheric_field,
      surface_field,
      {surface_field.single_value(SurfaceKey::h, latitude, longitude),
       latitude,
       longitude},
      {0, 0},
      max_step,
      1.0,
      true,
      false,
      true,
      true,
      false);
}

void disort_solar_sourceTurnOff(Vector& disort_solar_source,
                                Vector& disort_solar_zenith_angle,
                                Vector& disort_solar_azimuth_angle,
                                const AscendingGrid& frequency_grid) {
  const auto n = frequency_grid.size();

  disort_solar_source.resize(n);
  disort_solar_zenith_angle.resize(n);
  disort_solar_azimuth_angle.resize(n);

  disort_solar_source        = 0.0;
  disort_solar_zenith_angle  = 0.0;
  disort_solar_azimuth_angle = 0.0;
}

void disort_source_polynomialTurnOff(
    Tensor3& disort_source_polynomial,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& frequency_grid) {
  disort_source_polynomial.resize(
      frequency_grid.size(), ray_path.size() - 1, 0);
}

void disort_source_polynomialLinearInTau(
    Tensor3& disort_source_polynomial,
    const Matrix& disort_optical_thicknesses,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const AscendingGrid& frequency_grid) {
  const auto nv = frequency_grid.size();
  const auto N  = ray_path_atmospheric_point.size();

  disort_source_polynomial.resize(nv, N - 1, 2);

  ARTS_USER_ERROR_IF(
      not same_shape(std::array<Index, 2>{nv, static_cast<Index>(N) - 1},
                     disort_optical_thicknesses),
      "Incorrect shape for optical thicknesses.");

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    const Numeric& f = frequency_grid[iv];

    for (Size i = 0; i < N - 1; i++) {
      const Numeric& t0 = ray_path_atmospheric_point[i + 0].temperature;
      const Numeric& t1 = ray_path_atmospheric_point[i + 1].temperature;

      const Numeric y0 = planck(f, t0);
      const Numeric y1 = planck(f, t1);

      const Numeric x0 = i == 0 ? 0.0 : disort_optical_thicknesses(iv, i - 1);
      const Numeric x1 = disort_optical_thicknesses(iv, i);

      const Numeric b                    = (y1 - y0) / (x1 - x0);
      disort_source_polynomial(iv, i, 0) = y0 - b * x0;
      disort_source_polynomial(iv, i, 1) = b;
    }
  }
}

void disort_negative_boundary_conditionTurnOff(
    Tensor3& disort_negative_boundary_condition,
    const AscendingGrid& frequency_grid,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: ",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: ",
                     disort_fourier_mode_dimension)

  disort_negative_boundary_condition.resize(frequency_grid.size(),
                                            disort_fourier_mode_dimension,
                                            disort_quadrature_dimension / 2);
  disort_negative_boundary_condition = 0.0;
}

void disort_positive_boundary_conditionTurnOff(
    Tensor3& disort_positive_boundary_condition,
    const AscendingGrid& frequency_grid,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: ",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: ",
                     disort_fourier_mode_dimension)

  disort_positive_boundary_condition.resize(frequency_grid.size(),
                                            disort_fourier_mode_dimension,
                                            disort_quadrature_dimension / 2);
  disort_positive_boundary_condition = 0.0;
}

void disort_positive_boundary_conditionCosmicBackgroundRadiation(
    Tensor3& disort_positive_boundary_condition,
    const AscendingGrid& frequency_grid,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: ",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: ",
                     disort_fourier_mode_dimension)

  const auto nv = frequency_grid.size();
  disort_positive_boundary_condition.resize(
      nv, disort_fourier_mode_dimension, disort_quadrature_dimension / 2);
  disort_positive_boundary_condition = 0.0;

  for (Index iv = 0; iv < nv; iv++) {
    disort_positive_boundary_condition(iv, 0, joker) = planck(
        frequency_grid[iv], Constant::cosmic_microwave_background_temperature);
  }
}

void disort_negative_boundary_conditionSurfaceTemperature(
    Tensor3& disort_negative_boundary_condition,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: ",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: ",
                     disort_fourier_mode_dimension)

  const auto nv = frequency_grid.size();
  disort_negative_boundary_condition.resize(
      nv, disort_fourier_mode_dimension, disort_quadrature_dimension / 2);
  disort_negative_boundary_condition = 0.0;

  const Numeric T = surface_field.single_value(
      SurfaceKey::t, ray_path_point.latitude(), ray_path_point.longitude());

  for (Index iv = 0; iv < nv; iv++) {
    disort_negative_boundary_condition(iv, 0, joker) =
        planck(frequency_grid[iv], T);
  }
}

void disort_legendre_coefficientsTurnOff(
    Tensor3& disort_legendre_coefficients,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& frequency_grid,
    const Index& disort_legendre_polynomial_dimension) {
  ARTS_USER_ERROR_IF(disort_legendre_polynomial_dimension < 0,
                     "Legendre polynomial must be positive: ",
                     disort_legendre_polynomial_dimension)

  const Size N   = ray_path.size();
  const Index nv = frequency_grid.size();

  disort_legendre_coefficients.resize(
      nv, N - 1, disort_legendre_polynomial_dimension);

  disort_legendre_coefficients                  = 0.0;
  disort_legendre_coefficients(joker, joker, 0) = 1.0;
}

void disort_fractional_scatteringTurnOff(
    Matrix& disort_fractional_scattering,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& frequency_grid) {
  const Size N   = ray_path.size();
  const Index nv = frequency_grid.size();

  disort_fractional_scattering.resize(nv, N - 1);
  disort_fractional_scattering = 0.0;
}

void disort_single_scattering_albedoTurnOff(
    Matrix& disort_single_scattering_albedo,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& frequency_grid) {
  const Size N   = ray_path.size();
  const Index nv = frequency_grid.size();

  disort_single_scattering_albedo.resize(nv, N - 1);
  disort_single_scattering_albedo = 0.0;
}

void disort_optical_thicknessesFromPath(
    Matrix& disort_optical_thicknesses,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfPropmatVector& ray_path_propagation_matrix) {
  const Size N = ray_path.size();
  ARTS_USER_ERROR_IF(N < 1, "Must have a non-empty ray path.")

  ARTS_USER_ERROR_IF(
      not all_same_size(ray_path, ray_path_propagation_matrix),
      "Propagation matrices and frequency grids must have the same size.")

  const Index nv = ray_path_propagation_matrix.front().size();

  ARTS_USER_ERROR_IF(
      not all_same_shape(std::array{nv}, ray_path_propagation_matrix),
      "Propagation matrices and frequency grids must have the same shape.")

  // No polarization allowed
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(ray_path_propagation_matrix,
                          [](const PropmatVector& pms) {
                            return std::ranges::any_of(
                                pms, Cmp::eq(true), &Propmat::is_polarized);
                          }),
      "No implementation for polarized propagation matrices.");

  // Altitude is increasing
  ARTS_USER_ERROR_IF(
      std::ranges::is_sorted(
          ray_path,
          [](const PropagationPathPoint& a, const PropagationPathPoint& b) {
            return a.altitude() > b.altitude();
          }),
      "Ray path points must be sorted by increasing altitude.");

  const Vector r = [n = N - 1, &ray_path]() {
    Vector out(n);
    for (Size i = 0; i < n; i++) {
      out[i] = ray_path[i + 1].altitude() - ray_path[i].altitude();
    }

    return out;
  }();

  disort_optical_thicknesses.resize(nv, N - 1);

  for (Index iv = 0; iv < nv; iv++) {
    for (Size i = 0; i < N - 1; i++) {
      disort_optical_thicknesses(iv, i) =
          r[i] * std::midpoint(ray_path_propagation_matrix[i + 1][iv].A(),
                               ray_path_propagation_matrix[i + 0][iv].A());
      if (i > 0) {
        disort_optical_thicknesses(iv, i) +=
            disort_optical_thicknesses(iv, i - 1);
      }
    }
  }
}

void disort_bidirectional_reflectance_distribution_functionsTurnOff(
    MatrixOfDisortBDRF& disort_bidirectional_reflectance_distribution_functions,
    const AscendingGrid& frequency_grid) {
  const auto nv = frequency_grid.size();
  disort_bidirectional_reflectance_distribution_functions.resize(nv, 0);
}

void disort_bidirectional_reflectance_distribution_functionsLambertianConstant(
    MatrixOfDisortBDRF& disort_bidirectional_reflectance_distribution_functions,
    const AscendingGrid& frequency_grid,
    const Numeric& value) {
  const auto nv = frequency_grid.size();

  disort_bidirectional_reflectance_distribution_functions.resize(nv, 1);

  const auto f =
      DisortBDRF{[value](ExhaustiveMatrixView x,
                         const ExhaustiveConstVectorView&,
                         const ExhaustiveConstVectorView&) { x = value; }};

  disort_bidirectional_reflectance_distribution_functions = f;
}

void disort_spectral_radiance_fieldCalc(
    Tensor3& disort_spectral_radiance_field,
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
    const Vector& disort_solar_source) {
  ARTS_USER_ERROR_IF(not same_shape(disort_optical_thicknesses,
                                    disort_single_scattering_albedo,
                                    disort_fractional_scattering),
                     std::format(R"(Bad shapes: {:B,} != {:B,} != {:B,}

disort_optical_thicknesses.shape()      = {:B,},
disort_single_scattering_albedo.shape() = {:B,},
disort_fractional_scattering.shape()    = {:B,}
)",
                                 disort_optical_thicknesses.shape(),
                                 disort_single_scattering_albedo.shape(),
                                 disort_fractional_scattering.shape(),
                                 disort_optical_thicknesses.shape(),
                                 disort_single_scattering_albedo.shape(),
                                 disort_fractional_scattering.shape()))

  const Index nv = disort_optical_thicknesses.nrows();
  const Index np = disort_optical_thicknesses.ncols();

  ARTS_USER_ERROR_IF(not same_shape(std::array{nv},
                                    disort_solar_zenith_angle,
                                    disort_solar_azimuth_angle,
                                    disort_solar_source),
                     std::format(R"(Bad shapes: [{}] != {:B,} != {:B,} != {:B,}

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
                                 disort_solar_source.shape()))

  const auto nleg = disort_legendre_coefficients.ncols();
  ARTS_USER_ERROR_IF(disort_legendre_coefficients.npages() != nv or
                         disort_legendre_coefficients.nrows() != np,
                     std::format(R"(Bad shapes: [{}, {}, {}] != {:B,}

disort_legendre_coefficients.shape() = {:B,}
)",
                                 nv,
                                 np,
                                 nleg,
                                 disort_legendre_coefficients.shape(),
                                 disort_legendre_coefficients.shape()))

  const Index nfou       = disort_negative_boundary_condition.nrows();
  const Index nquad_half = disort_negative_boundary_condition.ncols();
  const Index nquad      = nquad_half * 2;

  ARTS_USER_ERROR_IF(not same_shape(std::array{nv, nfou, nquad_half},
                                    disort_negative_boundary_condition,
                                    disort_positive_boundary_condition),
                     std::format(R"(Bad shapes: [{}, {}, {}] != {:B,} != {:B,}

disort_negative_boundary_condition.shape() = {:B,},
disort_positive_boundary_condition.shape() = {:B,}
)",
                                 nv,
                                 nfou,
                                 nquad_half,
                                 disort_negative_boundary_condition.shape(),
                                 disort_positive_boundary_condition.shape(),
                                 disort_negative_boundary_condition.shape(),
                                 disort_positive_boundary_condition.shape()))

  const Index nsource = disort_source_polynomial.ncols();

  ARTS_USER_ERROR_IF(disort_source_polynomial.npages() != nv or
                         disort_source_polynomial.nrows() != np,
                     std::format(R"(Bad shapes: [{}, {}, {}] != {:B,}

disort_source_polynomial.shape() = {:B,}
)",
                                 nv,
                                 np,
                                 nsource,
                                 disort_source_polynomial.shape(),
                                 disort_source_polynomial.shape()))

  const auto nbdrf =
      disort_bidirectional_reflectance_distribution_functions.ncols();

  ARTS_USER_ERROR_IF(
      disort_bidirectional_reflectance_distribution_functions.nrows() != nv,
      std::format(
          R"(Bad shapes: [{} {}] != {:B,}

disort_bidirectional_reflectance_distribution_functions.shape() = {:B,}
)",
          nv,
          nbdrf,
          disort_bidirectional_reflectance_distribution_functions.shape(),
          disort_bidirectional_reflectance_distribution_functions.shape()))

  disort_spectral_radiance_field.resize(nv, np, nquad);

  const Index nleg_reduced = nleg;
  disort::main_data dis(np, nquad, nleg, nfou, nsource, nleg_reduced, nbdrf);

  //! Supplementary outputs
  disort_quadrature_weights = dis.weights();
  disort_quadrature_angles.resize(nquad);
  std::transform(dis.mu().begin(),
                 dis.mu().end(),
                 disort_quadrature_angles.begin(),
                 [](const Numeric& mu) { return Conversion::acosd(mu); });

  String error;

#pragma omp parallel for if (not arts_omp_in_parallel()) firstprivate(dis)
  for (Index iv = 0; iv < nv; iv++) {
    try {
      for (Index i = 0; i < nbdrf; i++) {
        dis.brdf_modes()[i] =
            disort_bidirectional_reflectance_distribution_functions(iv, i);
      }

      dis.solar_zenith() = Conversion::cosd(disort_solar_zenith_angle[iv]);
      dis.beam_azimuth() = Conversion::deg2rad(disort_solar_azimuth_angle[iv]);
      dis.tau()          = disort_optical_thicknesses[iv];
      dis.omega()        = disort_single_scattering_albedo[iv];
      dis.f()            = disort_fractional_scattering[iv];
      dis.all_legendre_coeffs() = disort_legendre_coefficients[iv];
      dis.positive_boundary()   = disort_positive_boundary_condition[iv];
      dis.negative_boundary()   = disort_negative_boundary_condition[iv];
      dis.source_poly()         = disort_source_polynomial[iv];

      dis.update_all(disort_solar_source[iv]);

      dis.gridded_u(disort_spectral_radiance_field[iv].reshape_as(np, 1, nquad),
                    {0.0});
    } catch (const std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(error.size(), "Error occurred in disort:\n", error);
}

void spectral_radianceIntegrateDisort(
    StokvecVector& /*spectral_radiance*/,
    const Tensor3& /*disort_spectral_radiance_field*/,
    const Vector& /*disort_quadrature_angles*/,
    const Vector& /*disort_quadrature_weights*/) {
  ARTS_USER_ERROR("Not implemented")
}

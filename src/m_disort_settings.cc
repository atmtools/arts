#include <arts_omp.h>
#include <disort.h>
#include <matpack.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>

#include "mh_checks.h"

////////////////////////////////////////////////////////////////////////////////
// Disort solar source
////////////////////////////////////////////////////////////////////////////////

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


////////////////////////////////////////////////////////////////////////////////
// Disort source polynomial
////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////
// Disort negative boundary condition
////////////////////////////////////////////////////////////////////////////////

void disort_negative_boundary_conditionTurnOff(
    Tensor3& disort_negative_boundary_condition,
    const AscendingGrid& frequency_grid,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: {}",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: {}",
                     disort_fourier_mode_dimension)

  disort_negative_boundary_condition.resize(frequency_grid.size(),
                                            disort_fourier_mode_dimension,
                                            disort_quadrature_dimension / 2);
  disort_negative_boundary_condition = 0.0;
}

void disort_negative_boundary_conditionSurfaceTemperature(
    Tensor3& disort_negative_boundary_condition,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: {}",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: {}",
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

////////////////////////////////////////////////////////////////////////////////
// Disort positive boundary condition
////////////////////////////////////////////////////////////////////////////////

void disort_positive_boundary_conditionTurnOff(
    Tensor3& disort_positive_boundary_condition,
    const AscendingGrid& frequency_grid,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension) {
  ARTS_USER_ERROR_IF(disort_quadrature_dimension < 0,
                     "Quadratures must be positive: {}",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: {}",
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
                     "Quadratures must be positive: {}",
                     disort_quadrature_dimension)
  ARTS_USER_ERROR_IF(disort_fourier_mode_dimension < 0,
                     "Fourier modes must be positive: {}",
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

////////////////////////////////////////////////////////////////////////////////
// Disort Legendre coefficients
////////////////////////////////////////////////////////////////////////////////

void disort_legendre_coefficientsTurnOff(
    Tensor3& disort_legendre_coefficients,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& frequency_grid,
    const Index& disort_legendre_polynomial_dimension) {
  ARTS_USER_ERROR_IF(disort_legendre_polynomial_dimension < 0,
                     "Legendre polynomial must be positive: {}",
                     disort_legendre_polynomial_dimension)

  const Size N   = ray_path.size();
  const Index nv = frequency_grid.size();

  disort_legendre_coefficients.resize(
      nv, N - 1, disort_legendre_polynomial_dimension);

  disort_legendre_coefficients                  = 0.0;
  disort_legendre_coefficients(joker, joker, 0) = 1.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort fractional scattering
////////////////////////////////////////////////////////////////////////////////

void disort_fractional_scatteringTurnOff(
    Matrix& disort_fractional_scattering,
    const ArrayOfPropagationPathPoint& ray_path,
    const AscendingGrid& frequency_grid) {
  const Size N   = ray_path.size();
  const Index nv = frequency_grid.size();

  disort_fractional_scattering.resize(nv, N - 1);
  disort_fractional_scattering = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort optical thicknesses
////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////
// Disort BRDF / BDRF
////////////////////////////////////////////////////////////////////////////////

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

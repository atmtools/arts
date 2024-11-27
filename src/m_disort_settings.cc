#include <arts_omp.h>
#include <disort.h>
#include <matpack.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <workspace.h>

#include <numeric>

#include "atm.h"
#include "debug.h"
#include "enumsSurfaceKey.h"
#include "mh_checks.h"
#include "sorted_grid.h"
#include "sun_methods.h"

////////////////////////////////////////////////////////////////////////////////
// Disort settings initialization
////////////////////////////////////////////////////////////////////////////////

void disort_settingsInit(DisortSettings& disort_settings,
                         const AscendingGrid& frequency_grid,
                         const ArrayOfPropagationPathPoint& ray_path,
                         const Index& quadrature_dimension,
                         const Index& legendre_polynomial_dimension,
                         const Index& fourier_mode_dimension) {
  disort_settings   = DisortSettings();
  const Index nfreq = frequency_grid.size();
  const Index nlay  = ray_path.size() - 1;
  disort_settings.resize(quadrature_dimension,
                         legendre_polynomial_dimension,
                         fourier_mode_dimension,
                         nfreq,
                         nlay);
}

////////////////////////////////////////////////////////////////////////////////
// Disort solar source
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSun(DisortSettings& disort_settings) {
  disort_settings.solar_source        = 0.0;
  disort_settings.solar_zenith_angle  = 0.0;
  disort_settings.solar_azimuth_angle = 0.0;
}

void disort_settingsSetSun(DisortSettings& disort_settings,
                           const AscendingGrid& frequency_grid,
                           const SurfaceField& surface_field,
                           const Sun& sun,
                           const PropagationPathPoint& ray_path_point) {
  const Numeric h =
      surface_field.single_value(SurfaceKey::h, sun.latitude, sun.longitude);

  const Vector3 sun_pos{sun.distance - h, sun.latitude, sun.longitude};
  const Vector2 los =
      geometric_los(ray_path_point.pos, sun_pos, surface_field.ellipsoid);

  const Numeric sin2_alpha =
      sun.sin_alpha_squared(ray_path_point.pos, surface_field.ellipsoid);

  const Index nv = frequency_grid.nelem();

  ARTS_USER_ERROR_IF(
      disort_settings.solar_source.size() != nv or sun.spectrum.nrows() != nv,
      R"(Solar spectrum not agreeing with frequency grids:

frequency_grid.nelem():              {}
disort_settings.solar_source.size(): {}
sun.spectrum.nrows():                {}
)",
      nv,
      disort_settings.solar_source.size(),
      sun.spectrum.nrows())

  for (Index iv = 0; iv < nv; iv++) {
    disort_settings.solar_source[iv] = sun.spectrum(iv, 0) * sin2_alpha;
  }

  disort_settings.solar_zenith_angle  = los[0];
  disort_settings.solar_azimuth_angle = los[1];
}

////////////////////////////////////////////////////////////////////////////////
// Disort source polynomial
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoLayerThermalEmission(DisortSettings& disort_settings) {
  disort_settings.source_polynomial.resize(
      disort_settings.nfreq, disort_settings.nlay, 0);
}

void disort_settingsLayerThermalEmissionLinearInTau(
    DisortSettings& disort_settings,
    const ArrayOfAtmPoint& ray_path_atmospheric_point,
    const AscendingGrid& frequency_grid) {
  const auto nv = frequency_grid.size();
  const auto N  = ray_path_atmospheric_point.size();

  disort_settings.source_polynomial.resize(nv, N - 1, 2);

  ARTS_USER_ERROR_IF(
      not same_shape(std::array<Index, 2>{nv, static_cast<Index>(N) - 1},
                     disort_settings.optical_thicknesses),
      "Incorrect shape: {:B,} vs {:B,}",
      std::array<Index, 2>{nv, static_cast<Index>(N) - 1},
      disort_settings.optical_thicknesses.shape());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Index iv = 0; iv < nv; iv++) {
    const Numeric& f = frequency_grid[iv];

    for (Size i = 0; i < N - 1; i++) {
      const Numeric& t0 = ray_path_atmospheric_point[i + 0].temperature;
      const Numeric& t1 = ray_path_atmospheric_point[i + 1].temperature;

      const Numeric y0 = planck(f, t0);
      const Numeric y1 = planck(f, t1);

      const Numeric x0 =
          i == 0 ? 0.0 : disort_settings.optical_thicknesses(iv, i - 1);
      const Numeric x1 = disort_settings.optical_thicknesses(iv, i);

      const Numeric b                             = (y1 - y0) / (x1 - x0);
      disort_settings.source_polynomial(iv, i, 0) = y0 - b * x0;
      disort_settings.source_polynomial(iv, i, 1) = b;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort positive boundary condition
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSurfaceEmission(DisortSettings& disort_settings) {
  disort_settings.positive_boundary_condition = 0.0;
}

void disort_settingsSurfaceEmissionByTemperature(
    DisortSettings& disort_settings,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field) {
  const auto nv = frequency_grid.size();

  const Numeric T = surface_field.single_value(
      SurfaceKey::t, ray_path_point.latitude(), ray_path_point.longitude());

  disort_settings.positive_boundary_condition = 0.0;

  ARTS_USER_ERROR_IF(
      nv != disort_settings.positive_boundary_condition.npages(),
      "Frequency grid size does not match the positive boundary condition size: {} vs {}",
      nv,
      disort_settings.positive_boundary_condition.npages())

  ARTS_USER_ERROR_IF(
      disort_settings.positive_boundary_condition.nrows() < 1,
      "Must have at least one fourier mode to use the positive boundary condition.")

  for (Index iv = 0; iv < nv; iv++) {
    disort_settings.positive_boundary_condition(iv, 0, joker) =
        planck(frequency_grid[iv], T);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort negative boundary condition
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSpaceEmission(DisortSettings& disort_settings) {
  disort_settings.negative_boundary_condition = 0.0;
}

void disort_settingsCosmicMicrowaveBackgroundRadiation(
    DisortSettings& disort_settings, const AscendingGrid& frequency_grid) {
  const Index nv = frequency_grid.size();

  disort_settings.negative_boundary_condition = 0.0;

  ARTS_USER_ERROR_IF(
      nv != disort_settings.negative_boundary_condition.npages(),
      "Frequency grid size does not match the negative boundary condition size: {} vs {}",
      nv,
      disort_settings.negative_boundary_condition.npages())

  ARTS_USER_ERROR_IF(
      disort_settings.negative_boundary_condition.nrows() < 1,
      "Must have at leaat one fourier mode to use the negative boundary condition.")

  for (Index iv = 0; iv < nv; iv++) {
    disort_settings.negative_boundary_condition(iv, 0, joker) = planck(
        frequency_grid[iv], Constant::cosmic_microwave_background_temperature);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort Legendre coefficients
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoLegendre(DisortSettings& disort_settings) {
  ARTS_USER_ERROR_IF(
      disort_settings.legendre_coefficients.nrows() < 1,
      "Must have at least one Legendre mode to use the Legendre coefficients.")

  disort_settings.legendre_coefficients                  = 0.0;
  disort_settings.legendre_coefficients(joker, joker, 0) = 1.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort fractional scattering
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoFractionalScattering(DisortSettings& disort_settings) {
  disort_settings.fractional_scattering = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort optical thicknesses
////////////////////////////////////////////////////////////////////////////////

void disort_settingsOpticalThicknessFromPath(
    DisortSettings& disort_settings,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfPropmatVector& ray_path_propagation_matrix) {
  const Index N  = disort_settings.nlay;
  const Index nv = disort_settings.nfreq;

  ARTS_USER_ERROR_IF(ray_path.size() != ray_path_propagation_matrix.size() or
                         ray_path.size() != static_cast<Size>(N + 1),
                     "Wrong path size.")

  if (N == 0) return;

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

  const Vector r = [n = N, &ray_path]() {
    Vector out(n);
    for (Index i = 0; i < n; i++) {
      out[i] = ray_path[i].altitude() - ray_path[i + 1].altitude();
    }

    return out;
  }();

  ARTS_USER_ERROR_IF(std::ranges::any_of(r, Cmp::le(0.0)),
                     R"(Atmospheric layer thickness must be positive.

Values:   {:B,}

Ray path points must be sorted by decreasing altitude.

ray_path: {:B,}
)",
                     r,
                     ray_path);

  for (Index iv = 0; iv < nv; iv++) {
    for (Index i = 0; i < N; i++) {
      disort_settings.optical_thicknesses(iv, i) =
          r[i] * std::midpoint(ray_path_propagation_matrix[i + 1][iv].A(),
                               ray_path_propagation_matrix[i + 0][iv].A());
      if (i > 0) {
        disort_settings.optical_thicknesses(iv, i) +=
            disort_settings.optical_thicknesses(iv, i - 1);

        ARTS_USER_ERROR_IF(disort_settings.optical_thicknesses(iv, i) <=
                               disort_settings.optical_thicknesses(iv, i - 1),
                           R"(
Not strictly increasing optical thicknesses between layers.

Check *ray_path_propagation_matrix* contain zeroes or negative values for A().

Old value:            {}
New value:            {}
Old layer:            {}
New layer:            {}
Frequency grid point: {}
)",
                           disort_settings.optical_thicknesses(iv, i - 1),
                           disort_settings.optical_thicknesses(iv, i),
                           i - 1,
                           i,
                           iv);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort BRDF / BDRF
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSurfaceScattering(DisortSettings& disort_settings) {
  disort_settings.bidirectional_reflectance_distribution_functions.resize(
      disort_settings.nfreq, 0);
}

void disort_settingsSurfaceLambertian(DisortSettings& disort_settings,
                                      const Numeric& value) {
  disort_settings.bidirectional_reflectance_distribution_functions.resize(
      disort_settings.nfreq, 1);

  const auto f =
      DisortBDRF{[value](ExhaustiveMatrixView x,
                         const ExhaustiveConstVectorView&,
                         const ExhaustiveConstVectorView&) { x = value; }};

  disort_settings.bidirectional_reflectance_distribution_functions = f;
}

////////////////////////////////////////////////////////////////////////////////
// Disort Single scattering albedo
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSingleScatteringAlbedo(DisortSettings& disort_settings) {
  disort_settings.single_scattering_albedo = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort Scattering Species Interface
////////////////////////////////////////////////////////////////////////////////

void disort_settingsLegendreCoefficientsFromPath(
    DisortSettings& disort_settings,
    const ArrayOfMatrix&
        ray_path_phase_matrix_scattering_totally_random_orientation_spectral) try {
  const Size N  = disort_settings.nlay;
  const Index F = disort_settings.nfreq;
  const Index L = disort_settings.legendre_polynomial_dimension;

  ARTS_USER_ERROR_IF(
      (N + 1) !=
          ray_path_phase_matrix_scattering_totally_random_orientation_spectral
              .size(),
      R"(The number of levels in ray_path_phase_matrix_scattering_totally_random_orientation_spectral must be one more than the number of layers in disort_settings.

  disort_settings.nlay + 1:                                                    {}
  ray_path_phase_matrix_scattering_totally_random_orientation_spectral.size(): {}
)",
      N,
      ray_path_phase_matrix_scattering_totally_random_orientation_spectral
          .size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(
          std::array{F, L},
          ray_path_phase_matrix_scattering_totally_random_orientation_spectral),
      "The shape of ray_path_phase_matrix_scattering_totally_random_orientation_spectral must be {:B,}, at least one is not",
      std::array{F, L});

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{F, static_cast<Index>(N), L},
                     disort_settings.legendre_coefficients),
      R"(The shape of disort_settings.legendre_coefficients must be {:B,}, but is {:B,})",
      std::array{F, static_cast<Index>(N), L},
      disort_settings.legendre_coefficients.shape());

  Vector invfac(L);
  for (Index j = 0; j < L; j++) {
    invfac[j] = 1.0 / std::sqrt(2 * j + 1);
  }

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(3)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      for (Index j = 0; j < L; j++) {
        disort_settings.legendre_coefficients(iv, i, j) =
            invfac[j] *
            std::midpoint(
                ray_path_phase_matrix_scattering_totally_random_orientation_spectral
                    [i](iv, j),
                ray_path_phase_matrix_scattering_totally_random_orientation_spectral
                    [i + 1](iv, j));
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      // Disort wants the first value to be 1.0, so we normalize
      if (std::isnormal(disort_settings.legendre_coefficients(iv, i, 0))) {
        disort_settings.legendre_coefficients(iv, i, joker) /=
            disort_settings.legendre_coefficients(iv, i, 0);

        //! WARNING: Numerical instabiliy occurs for large j-values in invfac cf what scattering_species does
        for (auto& v : disort_settings.legendre_coefficients(iv, i, joker)) {
          v = std::clamp(v, -1.0, 1.0);
        }
      } else {
        disort_settings.legendre_coefficients(iv, i, joker) = 0.0;
        disort_settings.legendre_coefficients(iv, i, 0)     = 1.0;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void disort_settingsSingleScatteringAlbedoFromPath(
    DisortSettings& disort_settings,
    const ArrayOfPropmatVector&
        ray_path_propagation_matrix_scattering_totally_random_orientation_spectral,
    const ArrayOfStokvecVector&
        ray_path_absorption_vector_scattering_totally_random_orientation_spectral) try {
  const Size N  = disort_settings.nlay;
  const Index F = disort_settings.nfreq;

  ARTS_USER_ERROR_IF(
      (N + 1) !=
          ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
              .size(),
      R"(The number of levels in ray_path_propagation_matrix_scattering_totally_random_orientation_spectral must be one more than the number of layers in disort_settings.

  disort_settings.nlay:                                                              {}
  ray_path_propagation_matrix_scattering_totally_random_orientation_spectral.size(): {}
)",
      N,
      ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
          .size());

  ARTS_USER_ERROR_IF(
      (N + 1) !=
          ray_path_absorption_vector_scattering_totally_random_orientation_spectral
              .size(),
      R"(The number of levels in ray_path_absorption_vector_scattering_totally_random_orientation_spectral must be one more than the number of layers in disort_settings.

  disort_settings.nlay:                                                             {}
  ray_path_absorption_vector_scattering_totally_random_orientation_spectral.size(): {}
)",
      N,
      ray_path_absorption_vector_scattering_totally_random_orientation_spectral
          .size());

  ARTS_USER_ERROR_IF(
      not all_same_shape(
          std::array{F},
          ray_path_absorption_vector_scattering_totally_random_orientation_spectral,
          ray_path_propagation_matrix_scattering_totally_random_orientation_spectral),
      "The shape of ray_path_propagation_matrix_scattering_totally_random_orientation_spectral and ray_path_absorption_vector_scattering_totally_random_orientation_spectral must be {:B,}, at least one is not",
      std::array{F});

  ARTS_USER_ERROR_IF(
      not same_shape(std::array{F, static_cast<Index>(N)},
                     disort_settings.single_scattering_albedo),
      R"(The shape of disort_settings.single_scattering_albedo must be {:B,}, but is {:B,})",
      std::array{F, static_cast<Index>(N)},
      disort_settings.single_scattering_albedo.shape());

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      const Numeric x =
          1.0 -
          std::midpoint(
              (inv(ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
                       [i][iv]) *
               ray_path_absorption_vector_scattering_totally_random_orientation_spectral
                   [i][iv])
                  .I(),
              (inv(ray_path_propagation_matrix_scattering_totally_random_orientation_spectral
                       [i][iv]) *
               ray_path_absorption_vector_scattering_totally_random_orientation_spectral
                   [i][iv])
                  .I());
      disort_settings.single_scattering_albedo(iv, i) =
          std::isnormal(x) ? x : 0.0;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

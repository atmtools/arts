#include <array_algo.h>
#include <arts_omp.h>
#include <disort.h>
#include <legendre.h>
#include <matpack.h>
#include <path_point.h>
#include <physics_funcs.h>
#include <rtepack.h>
#include <sun_methods.h>
#include <workspace.h>

#include <exception>
#include <numeric>

////////////////////////////////////////////////////////////////////////////////
// Disort settings initialization
////////////////////////////////////////////////////////////////////////////////

void disort_settingsInit(DisortSettings& disort_settings,
                         const AscendingGrid& freq_grid,
                         const ArrayOfPropagationPathPoint& ray_path,
                         const Index& quadrature_dimension,
                         const Index& legendre_polynomial_dimension,
                         const Index& fourier_mode_dimension) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(quadrature_dimension % 2,
                     "Quadrature dimension ({}) must be even",
                     quadrature_dimension);

  ARTS_USER_ERROR_IF(
      ray_path.size() < 2,
      "Must have at least one layer (two levels) to use disort solver");

  disort_settings.resize(quadrature_dimension,
                         legendre_polynomial_dimension,
                         fourier_mode_dimension,
                         freq_grid,
                         DescendingGrid{ray_path.begin(),
                                        ray_path.end(),
                                        [](auto& v) { return v.altitude(); }});
}

////////////////////////////////////////////////////////////////////////////////
// Disort solar source
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSun(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.solar_source        = 0.0;
  disort_settings.solar_zenith_angle  = 0.0;
  disort_settings.solar_azimuth_angle = 0.0;
}

void disort_settingsSetSun(DisortSettings& disort_settings,
                           const AscendingGrid& freq_grid,
                           const SurfaceField& surface_field,
                           const Sun& sun,
                           const PropagationPathPoint& ray_path_point) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surface_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surface_field.ellipsoid)

  const Numeric h =
      surface_field.single_value(SurfaceKey::h, sun.latitude, sun.longitude);

  const Vector3 sun_pos{sun.distance - h, sun.latitude, sun.longitude};
  const Vector2 los =
      geometric_los(ray_path_point.pos, sun_pos, surface_field.ellipsoid);

  const Numeric sin2_alpha =
      sun.sin_alpha_squared(ray_path_point.pos, surface_field.ellipsoid);

  const Size nv = freq_grid.size();

  ARTS_USER_ERROR_IF(disort_settings.solar_source.size() != nv or
                         static_cast<Size>(sun.spectrum.nrows()) != nv,
                     R"(Solar spectrum not agreeing with frequency grids:

freq_grid.size():               {}
disort_settings.solar_source.size(): {}
sun.spectrum.nrows():                {}
)",
                     nv,
                     disort_settings.solar_source.size(),
                     sun.spectrum.nrows())

  for (Size iv = 0; iv < nv; iv++) {
    disort_settings.solar_source[iv] = sun.spectrum[iv, 0] * sin2_alpha;
  }

  disort_settings.solar_zenith_angle  = los[0];
  disort_settings.solar_azimuth_angle = los[1];
}

////////////////////////////////////////////////////////////////////////////////
// Disort source polynomial
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoLayerThermalEmission(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.source_polynomial.resize(
      disort_settings.frequency_count(), disort_settings.layer_count(), 0);
}

namespace {
template <typename T>
void disort_settingsLayerThermalEmissionLinearInTauImpl(
    DisortSettings& disort_settings,
    const T& ray_path_points,
    const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  const Size nv = freq_grid.size();
  const Size N  = ray_path_points.size();

  disort_settings.source_polynomial.resize(nv, N - 1, 2);

  ARTS_USER_ERROR_IF(
      not same_shape<2>({static_cast<Index>(nv), static_cast<Index>(N) - 1},
                        disort_settings.optical_thicknesses),
      "Incorrect shape: [{}, {}] vs {:B,}",
      nv,
      N - 1,
      disort_settings.optical_thicknesses.shape());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const Numeric& f = freq_grid[iv];

    for (Size i = 0; i < N - 1; i++) {
      const Numeric& t0 = ray_path_points[i + 0].temperature;
      const Numeric& t1 = ray_path_points[i + 1].temperature;

      const Numeric y0 = planck(f, t0);
      const Numeric y1 = planck(f, t1);

      const Numeric x0 =
          i == 0 ? 0.0 : disort_settings.optical_thicknesses[iv, i - 1];
      const Numeric x1 = disort_settings.optical_thicknesses[iv, i];

      const Numeric b                             = (y1 - y0) / (x1 - x0);
      disort_settings.source_polynomial[iv, i, 0] = y0 - b * x0;
      disort_settings.source_polynomial[iv, i, 1] = b;
    }
  }
}
}  // namespace

void disort_settingsLayerThermalEmissionLinearInTau(
    DisortSettings& disort_settings,
    const ArrayOfAtmPoint& ray_path_atm_point,
    const AscendingGrid& freq_grid) {
  disort_settingsLayerThermalEmissionLinearInTauImpl(
      disort_settings, ray_path_atm_point, freq_grid);
}

void disort_settingsSubsurfaceLayerThermalEmissionLinearInTau(
    DisortSettings& disort_settings,
    const ArrayOfSubsurfacePoint& subsurface_profile,
    const AscendingGrid& freq_grid) {
  disort_settingsLayerThermalEmissionLinearInTauImpl(
      disort_settings, subsurface_profile, freq_grid);
}

void disort_settingsLayerNonThermalEmissionLinearInTau(
    DisortSettings& disort_settings,
    const ArrayOfAtmPoint& ray_path_atm_point,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfStokvecVector&
        ray_path_propagation_matrix_source_vector_nonlte,
    const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  const Size nv = freq_grid.size();
  const Size N  = ray_path_atm_point.size();

  disort_settings.source_polynomial.resize(nv, N - 1, 2);

  ARTS_USER_ERROR_IF(
      not same_shape<2>({static_cast<Index>(nv), static_cast<Index>(N) - 1},
                        disort_settings.optical_thicknesses),
      "Incorrect shape: [{}, {}] vs {:B,}",
      nv,
      N - 1,
      disort_settings.optical_thicknesses.shape());

  ARTS_USER_ERROR_IF(
      not arr::same_size(ray_path_atm_point,
                         ray_path_propagation_matrix,
                         ray_path_propagation_matrix_source_vector_nonlte),
      R"(Not same size:

ray_path_atm_point.size():   {}
ray_path_propagation_matrix.size():  {}
ray_path_source_vector_nonlte.size(): {}
)",
      ray_path_atm_point.size(),
      ray_path_propagation_matrix.size(),
      ray_path_propagation_matrix_source_vector_nonlte.size());

  ARTS_USER_ERROR_IF(not arr::elemwise_same_size(
                         ray_path_propagation_matrix,
                         ray_path_propagation_matrix_source_vector_nonlte),
                     R"(Not same size:

ray_path_propagation_matrix.size():   {}
ray_path_source_vector_nonlte.size(): {}
)",
                     ray_path_propagation_matrix.size(),
                     ray_path_propagation_matrix_source_vector_nonlte.size());

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size iv = 0; iv < nv; iv++) {
    const Numeric& f = freq_grid[iv];

    for (Size i = 0; i < N - 1; i++) {
      const Numeric& t0 = ray_path_atm_point[i + 0].temperature;
      const Numeric& t1 = ray_path_atm_point[i + 1].temperature;

      const Muelmat invK0 = inv(ray_path_propagation_matrix[i + 0][iv]);
      const Muelmat invK1 = inv(ray_path_propagation_matrix[i + 1][iv]);

      const Stokvec& S0 =
          ray_path_propagation_matrix_source_vector_nonlte[i + 0][iv];
      const Stokvec& S1 =
          ray_path_propagation_matrix_source_vector_nonlte[i + 1][iv];

      Numeric y0 = planck(f, t0) + (invK0 * S0).I();
      Numeric y1 = planck(f, t1) + (invK1 * S1).I();

      // Numerically, these may be slightly negative even though they can mathematically only be exactly 0.
      if (y0 < 0) y0 = 0.0;
      if (y1 < 0) y1 = 0.0;

      const Numeric x0 =
          i == 0 ? 0.0 : disort_settings.optical_thicknesses[iv, i - 1];
      const Numeric x1 = disort_settings.optical_thicknesses[iv, i];

      const Numeric b                             = (y1 - y0) / (x1 - x0);
      disort_settings.source_polynomial[iv, i, 0] = y0 - b * x0;
      disort_settings.source_polynomial[iv, i, 1] = b;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort positive boundary condition - from "surface", "subsurface" or "below"
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSurfaceEmission(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.positive_boundary_condition = 0.0;
}

void disort_settingsSurfaceEmissionByTemperature(
    DisortSettings& disort_settings,
    const AscendingGrid& freq_grid,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field) {
  ARTS_TIME_REPORT

  const auto nv = freq_grid.size();

  const Numeric T = surface_field.single_value(
      SurfaceKey::t, ray_path_point.latitude(), ray_path_point.longitude());

  auto& limit = disort_settings.positive_boundary_condition = 0.0;

  ARTS_USER_ERROR_IF(
      static_cast<Index>(nv) != limit.npages(),
      "Frequency grid size does not match the positive boundary condition size: {} vs {}",
      nv,
      limit.npages())

  ARTS_USER_ERROR_IF(
      limit.nrows() < 1,
      "Must have at least one fourier mode to use the positive boundary condition.")

  for (Size iv = 0; iv < nv; iv++) {
    limit[iv, 0, joker] = planck(freq_grid[iv], T);
  }
}

void disort_settingsSubsurfaceEmissionByTemperature(
    DisortSettings& disort_settings,
    const AscendingGrid& freq_grid,
    const ArrayOfSubsurfacePoint& subsurface_profile) {
  ARTS_TIME_REPORT

  const auto nv = freq_grid.size();

  auto& limit = disort_settings.positive_boundary_condition = 0.0;

  ARTS_USER_ERROR_IF(subsurface_profile.size() < 2, "Need at least two points")

  ARTS_USER_ERROR_IF(
      static_cast<Index>(nv) != limit.npages(),
      "Frequency grid size does not match the positive boundary condition size: {} vs {}",
      nv,
      limit.npages())

  ARTS_USER_ERROR_IF(
      limit.nrows() < 1,
      "Must have at least one fourier mode to use the positive boundary condition.")

  const Numeric Tbot = subsurface_profile.back().temperature;

  for (Size iv = 0; iv < nv; iv++) {
    limit[iv, 0, joker] = planck(freq_grid[iv], Tbot);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort negative boundary condition - from "space", "atmosphere" or "above"
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSpaceEmission(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.negative_boundary_condition = 0.0;
}

void disort_settingsCosmicMicrowaveBackgroundRadiation(
    DisortSettings& disort_settings, const AscendingGrid& freq_grid) {
  ARTS_TIME_REPORT

  const Index nv = freq_grid.size();

  disort_settings.negative_boundary_condition = 0.0;

  ARTS_USER_ERROR_IF(
      nv != disort_settings.negative_boundary_condition.npages(),
      "Frequency grid size does not match the negative boundary condition size: {} vs {}",
      nv,
      disort_settings.negative_boundary_condition.npages())

  ARTS_USER_ERROR_IF(
      disort_settings.negative_boundary_condition.nrows() < 1,
      "Must have at least one fourier mode to use the negative boundary condition.")

  for (Index iv = 0; iv < nv; iv++) {
    disort_settings.negative_boundary_condition[iv, 0, joker] = planck(
        freq_grid[iv], Constant::cosmic_microwave_background_temperature);
  }
}

void disort_settingsDownwellingObserver(
    const Workspace& ws,
    DisortSettings& disort_settings,
    const AscendingGrid& freq_grid,
    const ArrayOfPropagationPathPoint& ray_path,
    const AtmField& atm_field,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const Agenda& spectral_radiance_observer_agenda,
    const Stokvec& pol) {
  ARTS_TIME_REPORT

  const auto& ray_path_point = ray_path.front();

  auto& limit = disort_settings.negative_boundary_condition = 0;

  const Index nv = freq_grid.size();
  const Index N  = disort_settings.quadrature_dimension / 2;

  ARTS_USER_ERROR_IF(
      nv != limit.npages(),
      "Frequency grid size does not match the boundary condition size: {} vs {}",
      nv,
      limit.npages())

  ARTS_USER_ERROR_IF(limit.nrows() < 1, "Must have at least one fourier mode.")

  ARTS_USER_ERROR_IF(limit.ncols() != N,
                     "Must have at least one quadrature dimension.")

  Vector mu(N);
  Vector W(N);
  Legendre::PositiveDoubleGaussLegendre(mu, W);

  StokvecVector spectral_radiance;
  StokvecMatrix spectral_radiance_jacobian;
  ArrayOfPropagationPathPoint ray_path_up;
  const JacobianTargets jacobian_targets{};

  String error{};

#pragma omp parallel for if (not arts_omp_in_parallel()) \
    firstprivate(spectral_radiance, spectral_radiance_jacobian, ray_path_up)
  for (Index i = 0; i < N; i++) {
    try {
      spectral_radiance_observer_agendaExecute(
          ws,
          spectral_radiance,
          spectral_radiance_jacobian,
          ray_path_up,
          freq_grid,
          jacobian_targets,
          ray_path_point.pos,
          {Conversion::acosd(mu[i]), ray_path_point.azimuth()},
          atm_field,
          surface_field,
          subsurface_field,
          spectral_radiance_observer_agenda);

      for (Index iv = 0; iv < nv; iv++) {
        limit[iv, 0, i] = dot(spectral_radiance[iv], pol);
      }
    } catch (std::exception& e) {
#pragma omp critical
      error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), error)
}

////////////////////////////////////////////////////////////////////////////////
// Disort Legendre coefficients
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoLegendre(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      disort_settings.legendre_coefficients.nrows() < 1,
      "Must have at least one Legendre mode to use the Legendre coefficients.")

  disort_settings.legendre_coefficients                  = 0.0;
  disort_settings.legendre_coefficients[joker, joker, 0] = 1.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort fractional scattering
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoFractionalScattering(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.fractional_scattering = 0.0;
}

////////////////////////////////////////////////////////////////////////////////
// Disort optical thicknesses
////////////////////////////////////////////////////////////////////////////////

void disort_settingsOpticalThicknessFromPath(
    DisortSettings& disort_settings,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const Numeric& min_optical_depth) {
  ARTS_TIME_REPORT

  const Index N  = disort_settings.layer_count();
  const Index nv = disort_settings.frequency_count();

  ARTS_USER_ERROR_IF(ray_path.size() != ray_path_propagation_matrix.size() or
                         ray_path.size() != static_cast<Size>(N + 1),
                     "Wrong path size.")

  if (N == 0) return;

  ARTS_USER_ERROR_IF(
      not all_same_shape<1>({nv}, ray_path_propagation_matrix),
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
      disort_settings.optical_thicknesses[iv, i] = std::max(
          r[i] * std::midpoint(ray_path_propagation_matrix[i + 1][iv].A(),
                               ray_path_propagation_matrix[i + 0][iv].A()),
          min_optical_depth);
      if (i > 0) {
        disort_settings.optical_thicknesses[iv, i] +=
            disort_settings.optical_thicknesses[iv, i - 1];

        ARTS_USER_ERROR_IF((disort_settings.optical_thicknesses[iv, i] <=
                            disort_settings.optical_thicknesses[iv, i - 1]),
                           R"(
Not strictly increasing optical thicknesses between layers.

Check *ray_path_propagation_matrix* contain zeroes or negative values for A().

Value:                {}
Frequency grid index: {}
)",
                           disort_settings.optical_thicknesses[iv, i - 1],
                           iv);
      }
    }
  }
}

void disort_settingsSubsurfaceScalarAbsorption(
    DisortSettings& disort_settings,
    const ArrayOfPropagationPathPoint& ray_path,
    const ArrayOfSubsurfacePoint& subsurface_profile,
    const Numeric& min_optical_depth) {
  ARTS_TIME_REPORT

  const Index N = disort_settings.layer_count();

  ARTS_USER_ERROR_IF(ray_path.size() != subsurface_profile.size() or
                         subsurface_profile.size() != static_cast<Size>(N + 1),
                     "Wrong path size.")

  if (N == 0) return;

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

  const SubsurfacePropertyTag absorption_tag{"scalar absorption"};

  ARTS_USER_ERROR_IF(
      not stdr::all_of(subsurface_profile, Cmp::contains(absorption_tag)),
      R"(Missing '{}' in some or all of the subsurface profile)",
      absorption_tag);

  for (Index i = 0; i < N; i++) {
    disort_settings.optical_thicknesses[joker, i] = std::max(
        r[i] * std::midpoint(subsurface_profile[i + 1][absorption_tag],
                             subsurface_profile[i + 0][absorption_tag]),
        min_optical_depth);
    if (i > 0) {
      disort_settings.optical_thicknesses[joker, i] +=
          disort_settings.optical_thicknesses[joker, i - 1];
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort BRDF / BDRF
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSurfaceScattering(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.bidirectional_reflectance_distribution_functions.resize(
      disort_settings.frequency_count(), 0);
}

void disort_settingsSurfaceLambertian(DisortSettings& disort_settings,
                                      const Vector& vec) {
  ARTS_TIME_REPORT

  disort_settings.bidirectional_reflectance_distribution_functions.resize(
      disort_settings.frequency_count(), 1);

  for (Index iv = 0; iv < disort_settings.frequency_count(); iv++) {
    disort_settings.bidirectional_reflectance_distribution_functions[iv, 0] =
        DisortBDRF{[value = vec[iv]](MatrixView x,
                                     const ConstVectorView&,
                                     const ConstVectorView&) { x = value; }};
  }
}

void disort_settingsSurfaceLambertian(DisortSettings& disort_settings,
                                      const Numeric& value) {
  ARTS_TIME_REPORT

  disort_settings.bidirectional_reflectance_distribution_functions.resize(
      disort_settings.frequency_count(), 1);

  const auto f = DisortBDRF{[value](MatrixView x,
                                    const ConstVectorView&,
                                    const ConstVectorView&) { x = value; }};

  disort_settings.bidirectional_reflectance_distribution_functions = f;
}

////////////////////////////////////////////////////////////////////////////////
// Disort Single scattering albedo
////////////////////////////////////////////////////////////////////////////////

void disort_settingsNoSingleScatteringAlbedo(DisortSettings& disort_settings) {
  ARTS_TIME_REPORT

  disort_settings.single_scattering_albedo = 0.0;
}

void disort_settingsSubsurfaceScalarSingleScatteringAlbedo(
    DisortSettings& disort_settings,
    const ArrayOfSubsurfacePoint& subsurface_profile) {
  ARTS_TIME_REPORT

  const Index N = disort_settings.layer_count();

  if (N == 0) return;

  const SubsurfacePropertyTag ssa_tag{"scalar ssa"};

  ARTS_USER_ERROR_IF(
      not stdr::all_of(subsurface_profile, Cmp::contains(ssa_tag)),
      R"(Missing '{}' in some or all of the subsurface profile)",
      ssa_tag);

  for (Index i = 0; i < N; i++) {
    disort_settings.single_scattering_albedo[joker, i] = std::midpoint(
        subsurface_profile[i + 1][ssa_tag], subsurface_profile[i + 0][ssa_tag]);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Disort Scattering Species Interface
////////////////////////////////////////////////////////////////////////////////

void disort_settingsLegendreCoefficientsFromPath(
    DisortSettings& disort_settings,
    const ArrayOfSpecmatMatrix& ray_path_phase_matrix_scattering_spectral) try {
  ARTS_TIME_REPORT

  const Size N  = disort_settings.layer_count();
  const Index F = disort_settings.frequency_count();
  const Index L = disort_settings.legendre_polynomial_dimension;

  ARTS_USER_ERROR_IF(
      (N + 1) != ray_path_phase_matrix_scattering_spectral.size(),
      R"(The number of levels in ray_path_phase_matrix_scattering_spectral must be one more than the number of layers in disort_settings.

  disort_settings.layer_count() + 1:                         {}
  ray_path_phase_matrix_scattering_spectral.size(): {}
)",
      N,
      ray_path_phase_matrix_scattering_spectral.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape<2>({static_cast<Index>(F), static_cast<Index>(L)},
                            ray_path_phase_matrix_scattering_spectral),
      "The shape of ray_path_phase_matrix_scattering_spectral must be {:B,}, at least one is not",
      std::array{F, L});

  ARTS_USER_ERROR_IF(
      not same_shape<3>({F, N, L}, disort_settings.legendre_coefficients),
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
        disort_settings.legendre_coefficients[iv, i, j] =
            invfac[j] *
            std::midpoint(
                ray_path_phase_matrix_scattering_spectral[i][iv, j][0, 0]
                    .real(),
                ray_path_phase_matrix_scattering_spectral[i + 1][iv, j][0, 0]
                    .real());
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      // Disort wants the first value to be 1.0, so we normalize
      if (std::isnormal(disort_settings.legendre_coefficients[iv, i, 0])) {
        disort_settings.legendre_coefficients[iv, i, joker] /=
            disort_settings.legendre_coefficients[iv, i, 0];

        //! WARNING: Numerical instabiliy occurs for large j-values in invfac cf what scattering_species does
        for (auto& v : disort_settings.legendre_coefficients[iv, i, joker]) {
          v = std::clamp(v, -1.0, 1.0);
        }
      } else {
        disort_settings.legendre_coefficients[iv, i, joker] = 0.0;
        disort_settings.legendre_coefficients[iv, i, 0]     = 1.0;
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void disort_settingsSingleScatteringAlbedoFromPath(
    DisortSettings& disort_settings,
    const ArrayOfPropmatVector& ray_path_propagation_matrix,
    const ArrayOfPropmatVector& ray_path_propagation_matrix_scattering,
    const ArrayOfStokvecVector& ray_path_absorption_vector_scattering) try {
  ARTS_TIME_REPORT

  const Size N  = disort_settings.layer_count();
  const Index F = disort_settings.frequency_count();

  ARTS_USER_ERROR_IF(
      (N + 1) != ray_path_propagation_matrix.size(),
      R"(The number of levels in ray_path_propagation_matrix must be one more than the number of layers in disort_settings.

  disort_settings.layer_count():                          {}
  ray_path_propagation_matrix.size(): {}
)",
      N,
      ray_path_propagation_matrix_scattering.size());

  ARTS_USER_ERROR_IF(
      (N + 1) != ray_path_propagation_matrix_scattering.size(),
      R"(The number of levels in ray_path_propagation_matrix_scattering must be one more than the number of layers in disort_settings.

  disort_settings.layer_count():                          {}
  ray_path_propagation_matrix_scattering.size(): {}
)",
      N,
      ray_path_propagation_matrix_scattering.size());

  ARTS_USER_ERROR_IF(
      (N + 1) != ray_path_absorption_vector_scattering.size(),
      R"(The number of levels in ray_path_absorption_vector_scattering must be one more than the number of layers in disort_settings.

  disort_settings.layer_count():                         {}
  ray_path_absorption_vector_scattering.size(): {}
)",
      N,
      ray_path_absorption_vector_scattering.size());

  ARTS_USER_ERROR_IF(
      not all_same_shape<1>({F},
                            ray_path_absorption_vector_scattering,
                            ray_path_propagation_matrix_scattering),
      "The shape of ray_path_propagation_matrix_scattering and ray_path_absorption_vector_scattering must be {:B,}, at least one is not",
      std::array{F});

  ARTS_USER_ERROR_IF(
      not same_shape<2>({F, N}, disort_settings.single_scattering_albedo),
      R"(The shape of disort_settings.single_scattering_albedo must be {:B,}, but is {:B,})",
      std::array{F, static_cast<Index>(N)},
      disort_settings.single_scattering_albedo.shape());

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
  for (Size i = 0; i < N; i++) {
    for (Index iv = 0; iv < F; iv++) {
      const Numeric ext_upper = ray_path_propagation_matrix[i][iv][0];
      const Numeric abs_scat_upper =
          ray_path_absorption_vector_scattering[i][iv][0];
      const Numeric ext_scat_upper =
          ray_path_propagation_matrix_scattering[i][iv][0];

      const Numeric ext_lower = ray_path_propagation_matrix[i + 1][iv][0];
      const Numeric abs_scat_lower =
          ray_path_absorption_vector_scattering[i + 1][iv][0];
      const Numeric ext_scat_lower =
          ray_path_propagation_matrix_scattering[i + 1][iv][0];

      const Numeric ext      = std::midpoint(ext_upper, ext_lower);
      const Numeric abs_scat = std::midpoint(abs_scat_upper, abs_scat_lower);
      const Numeric ext_scat = std::midpoint(ext_scat_upper, ext_scat_lower);

      // const Numeric abs = ext + abs_scat;

      const Numeric x = (ext_scat - abs_scat) / ext;
      disort_settings.single_scattering_albedo[iv, i] =
          std::isnormal(x) ? x : 0.0;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

namespace {
Agenda disort_settings_agendaSetup(
    const disort_settings_agenda_setup_layer_emission_type&
        layer_emission_setting,
    const disort_settings_agenda_setup_scattering_type& scattering_setting,
    const disort_settings_agenda_setup_space_type& space_setting,
    const disort_settings_agenda_setup_sun_type& sun_setting,
    const disort_settings_agenda_setup_surface_type& surface_setting,
    const Vector& surface_lambertian_value,
    const Numeric& min_optical_depth) {
  ARTS_TIME_REPORT

  AgendaCreator agenda("disort_settings_agenda");

  agenda.add("jacobian_targetsOff");

  // Clearsky absorption
  agenda.add("ray_path_atm_pointFromPath");
  agenda.add("freq_grid_pathFromPath");
  agenda.add("ray_path_propagation_matrixFromPath");

  agenda.add("disort_settingsInit");

  switch (scattering_setting) {
    using enum disort_settings_agenda_setup_scattering_type;
    case ScatteringSpecies:
      agenda.add("legendre_degreeFromDisortSettings");
      agenda.add("ray_path_propagation_matrix_scatteringFromSpectralAgenda");
      agenda.add("ray_path_propagation_matrixAddScattering");

      agenda.add("disort_settingsNoFractionalScattering");
      agenda.add("disort_settingsLegendreCoefficientsFromPath");
      agenda.add("disort_settingsSingleScatteringAlbedoFromPath");
      break;
    case None:
      agenda.add("disort_settingsNoFractionalScattering");
      agenda.add("disort_settingsNoLegendre");
      agenda.add("disort_settingsNoSingleScatteringAlbedo");
      break;
  }

  // We have both scattering and clearsky absorption, so we can set the optical thickness
  agenda.add("disort_settingsOpticalThicknessFromPath",
             SetWsv("min_optical_depth", min_optical_depth));

  // Since we have the optical thickness, we can set the thermal emission
  switch (layer_emission_setting) {
    using enum disort_settings_agenda_setup_layer_emission_type;
    case LinearInTau:
      agenda.add("disort_settingsLayerThermalEmissionLinearInTau");
      break;
    case LinearInTauNonLTE:
      agenda.add("disort_settingsLayerNonThermalEmissionLinearInTau");
      break;
    case None: agenda.add("disort_settingsNoLayerThermalEmission"); break;
  }

  switch (space_setting) {
    using enum disort_settings_agenda_setup_space_type;
    case None: agenda.add("disort_settingsNoSpaceEmission"); break;
    case CosmicMicrowaveBackgroundRadiation:
      agenda.add("disort_settingsCosmicMicrowaveBackgroundRadiation");
      break;
  }

  switch (surface_setting) {
    using enum disort_settings_agenda_setup_surface_type;
    case None:
      agenda.add("disort_settingsNoSurfaceEmission");
      agenda.add("disort_settingsNoSurfaceScattering");
      break;
    case Thermal:
      agenda.add("ray_path_pointLowestFromPath");
      agenda.add("disort_settingsSurfaceEmissionByTemperature");
      agenda.add("disort_settingsNoSurfaceScattering");
      break;
    case ThermalLambertian:
      agenda.add("ray_path_pointLowestFromPath");
      agenda.add("disort_settingsSurfaceEmissionByTemperature");
      agenda.add("disort_settingsSurfaceLambertian",
                 SetWsv{"value", surface_lambertian_value});
      break;
    case Lambertian:
      agenda.add("disort_settingsNoSurfaceEmission");
      agenda.add("disort_settingsSurfaceLambertian",
                 SetWsv{"value", surface_lambertian_value});
      break;
  }

  switch (sun_setting) {
    using enum disort_settings_agenda_setup_sun_type;
    case None: agenda.add("disort_settingsNoSun"); break;
    case Sun:
      agenda.add("ray_path_pointHighestFromPath");
      agenda.add("disort_settingsSetSun");
      break;
  }

  return std::move(agenda).finalize(false);
}

Agenda disort_settings_agendaSubsurfaceSetup(
    const disort_settings_agenda_setup_sun_type& sun_setting,
    const Numeric& min_optical_depth,
    const Index& fading_bottom) {
  ARTS_TIME_REPORT

  AgendaCreator agenda("disort_settings_agenda");

  agenda.add("jacobian_targetsOff");

  // Clearsky absorption
  agenda.add("subsurface_profileFromPath");

  agenda.add("disort_settingsInit");
  agenda.add("disort_settingsSubsurfaceScalarAbsorption",
             SetWsv("min_optical_depth", min_optical_depth));
  agenda.add("disort_settingsSubsurfaceScalarSingleScatteringAlbedo");

  if (fading_bottom) {
    agenda.add("disort_settingsNoSurfaceEmission");
  } else {
    agenda.add("disort_settingsSubsurfaceEmissionByTemperature");
  }
  agenda.add("disort_settingsSubsurfaceLayerThermalEmissionLinearInTau");
  agenda.add("disort_settingsNoLegendre");

  switch (sun_setting) {
    using enum disort_settings_agenda_setup_sun_type;
    case None: agenda.add("disort_settingsNoSun"); break;
    case Sun:
      agenda.add("ray_path_pointHighestFromPath");
      agenda.add("disort_settingsSetSun");
      break;
  }

  return std::move(agenda).finalize(false);
}
}  // namespace

void disort_settings_agendaSetup(Agenda& disort_settings_agenda,
                                 const String& layer_emission_setting,
                                 const String& scattering_setting,
                                 const String& space_setting,
                                 const String& sun_setting,
                                 const String& surface_setting,
                                 const Vector& surface_lambertian_value,
                                 const Numeric& min_optical_depth) {
  ARTS_TIME_REPORT

  disort_settings_agenda = disort_settings_agendaSetup(
      to<disort_settings_agenda_setup_layer_emission_type>(
          layer_emission_setting),
      to<disort_settings_agenda_setup_scattering_type>(scattering_setting),
      to<disort_settings_agenda_setup_space_type>(space_setting),
      to<disort_settings_agenda_setup_sun_type>(sun_setting),
      to<disort_settings_agenda_setup_surface_type>(surface_setting),
      surface_lambertian_value,
      min_optical_depth);
}

void disort_settings_agendaSubsurfaceSetup(Agenda& disort_settings_agenda,
                                           const String& sun_setting,
                                           const Numeric& min_optical_depth,
                                           const Index& fading_bottom) {
  ARTS_TIME_REPORT

  disort_settings_agenda = disort_settings_agendaSubsurfaceSetup(
      to<disort_settings_agenda_setup_sun_type>(sun_setting),
      min_optical_depth,
      fading_bottom);
}

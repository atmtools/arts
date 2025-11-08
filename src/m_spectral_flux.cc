#include <array_algo.h>
#include <arts_omp.h>
#include <workspace.h>

namespace {
void ray_path_spectral_radianceStepByStepEmissionForwardOnly(
    ArrayOfStokvecVector& ray_path_spectral_radiance,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix,
    const ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const StokvecVector& spectral_radiance_background) try {
  ARTS_TIME_REPORT

  ray_path_spectral_radiance.resize(ray_path_transmission_matrix.size());
  arr::elemwise_resize(spectral_radiance_background.size(),
                       ray_path_spectral_radiance);

  ray_path_spectral_radiance.back() = spectral_radiance_background;

  two_level_linear_emission_step_by_step_full(
      ray_path_spectral_radiance,
      ray_path_transmission_matrix,
      ray_path_spectral_radiance_source);
}
ARTS_METHOD_ERROR_CATCH

void ray_path_spectral_radianceClearskyEmission(
    const Workspace& ws,
    ArrayOfStokvecVector& ray_path_spectral_radiance,
    const AtmField& atm_field,
    const AscendingGrid& frequency_grid,
    const Agenda& propagation_matrix_agenda,
    const ArrayOfPropagationPathPoint& ray_path,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field) try {
  ARTS_TIME_REPORT

  PropagationPathPoint ray_path_point;
  ray_path_pointBackground(ray_path_point, ray_path);
  StokvecVector spectral_radiance_background;
  StokvecMatrix spectral_radiance_background_jacobian;
  spectral_radiance_backgroundAgendasAtEndOfPath(
      ws,
      spectral_radiance_background,
      spectral_radiance_background_jacobian,
      frequency_grid,
      {},
      ray_path_point,
      surface_field,
      subsurface_field,
      spectral_radiance_space_agenda,
      spectral_radiance_surface_agenda);
  ArrayOfAtmPoint ray_path_atm_point;
  ray_path_atm_pointFromPath(ray_path_atm_point, ray_path, atm_field);
  ArrayOfAscendingGrid ray_path_frequency_grid;
  ArrayOfVector3 ray_path_frequency_wind_shift_jacobian;
  ray_path_frequency_gridFromPath(ray_path_frequency_grid,
                                  ray_path_frequency_wind_shift_jacobian,
                                  frequency_grid,
                                  ray_path,
                                  ray_path_atm_point);
  ArrayOfPropmatVector ray_path_propagation_matrix;
  ArrayOfStokvecVector ray_path_propagation_matrix_source_vector_nonlte;
  ArrayOfPropmatMatrix ray_path_propagation_matrix_jacobian;
  ArrayOfStokvecMatrix
      ray_path_propagation_matrix_source_vector_nonlte_jacobian;
  ray_path_propagation_matrixFromPath(
      ws,
      ray_path_propagation_matrix,
      ray_path_propagation_matrix_source_vector_nonlte,
      ray_path_propagation_matrix_jacobian,
      ray_path_propagation_matrix_source_vector_nonlte_jacobian,
      propagation_matrix_agenda,
      ray_path_frequency_grid,
      ray_path_frequency_wind_shift_jacobian,
      {},
      ray_path,
      ray_path_atm_point);
  ArrayOfMuelmatVector ray_path_transmission_matrix;
  ArrayOfMuelmatTensor3 ray_path_transmission_matrix_jacobian;
  ray_path_transmission_matrixFromPath(ray_path_transmission_matrix,
                                       ray_path_transmission_matrix_jacobian,
                                       ray_path_propagation_matrix,
                                       ray_path_propagation_matrix_jacobian,
                                       ray_path,
                                       ray_path_atm_point,
                                       surface_field,
                                       {},
                                       0);
  ArrayOfStokvecVector ray_path_spectral_radiance_source;
  ArrayOfStokvecMatrix ray_path_spectral_radiance_source_jacobian;
  ray_path_spectral_radiance_sourceFromPropmat(
      ray_path_spectral_radiance_source,
      ray_path_spectral_radiance_source_jacobian,
      ray_path_propagation_matrix,
      ray_path_propagation_matrix_source_vector_nonlte,
      ray_path_propagation_matrix_jacobian,
      ray_path_propagation_matrix_source_vector_nonlte_jacobian,
      ray_path_frequency_grid,
      ray_path_atm_point,
      {});
  ray_path_spectral_radianceStepByStepEmissionForwardOnly(
      ray_path_spectral_radiance,
      ray_path_transmission_matrix,
      ray_path_spectral_radiance_source,
      spectral_radiance_background);
}
ARTS_METHOD_ERROR_CATCH
}  // namespace

void spectral_flux_profileFromPathField(
    const Workspace& ws,
    Matrix& spectral_flux_profile,
    const ArrayOfArrayOfPropagationPathPoint& ray_path_field,
    const AtmField& atm_field,
    const Agenda& propagation_matrix_agenda,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& altitude_grid) try {
  ARTS_TIME_REPORT

  const Size N = ray_path_field.size();
  const Size M = altitude_grid.size();
  const Size K = frequency_grid.size();

  spectral_flux_profile.resize(M, K);
  spectral_flux_profile = 0.0;

  ArrayOfArrayOfStokvecVector ray_path_spectral_radiance_field(N);

  String error{};
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size n = 0; n < N; n++) {
    try {
      ray_path_spectral_radianceClearskyEmission(
          ws,
          ray_path_spectral_radiance_field[n],
          atm_field,
          frequency_grid,
          propagation_matrix_agenda,
          ray_path_field[n],
          spectral_radiance_space_agenda,
          spectral_radiance_surface_agenda,
          surface_field,
          subsurface_field);
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  if (not error.empty()) throw std::runtime_error(error);

  ARTS_USER_ERROR_IF(
      not arr::each_same_size(ray_path_spectral_radiance_field, ray_path_field),
      "Not all ray paths have the same altitude count")

  for (auto& sradf : ray_path_spectral_radiance_field) {
    ARTS_USER_ERROR_IF(
        stdr::any_of(sradf,
                     [K](const StokvecVector& v) { return v.size() != K; }),
        "Not all ray path spectral radiances in the field have the same frequency count")
  }

  struct Zenith {
    Size outer;
    Size inner;
    Numeric angle;
  };

  std::vector<Zenith> zenith_angles(2 * M);  // Up and down
  for (Size m = 0; m < M; m++) {
    zenith_angles.resize(0);
    const Numeric alt = altitude_grid[m];
    VectorView t      = spectral_flux_profile[m];

    for (Size i = 0; i < ray_path_field.size(); i++) {
      for (Size j = 0; j < ray_path_field[i].size(); j++) {
        if (alt == ray_path_field[i][j].altitude()) {
          zenith_angles.emplace_back(i, j, ray_path_field[i][j].zenith());
        }
      }
    }

    ARTS_USER_ERROR_IF(
        zenith_angles.size() == 0, "No ray paths intersects altitude {} m", alt)

    stdr::sort(zenith_angles, {}, &Zenith::angle);

    // Integrate
    using Constant::pi;
    using Conversion::cosd;
    for (Size i = 0; i < zenith_angles.size() - 1; i++) {
      const auto& z0 = zenith_angles[i];
      const auto& z1 = zenith_angles[i + 1];

      const auto& y0 = ray_path_spectral_radiance_field[z0.outer][z0.inner];
      const auto& y1 = ray_path_spectral_radiance_field[z1.outer][z1.inner];

      const Numeric w = 0.5 * pi * (cosd(z0.angle) - cosd(z1.angle));
      for (Size k = 0; k < K; k++) t[k] += w * (y0[k][0] + y1[k][0]);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void flux_profileIntegrate(Vector& flux_profile,
                           const Matrix& spectral_flux_profile,
                           const AscendingGrid& frequency_grid) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(static_cast<Index>(frequency_grid.size()) !=
                         spectral_flux_profile.extent(1),
                     "Frequency grid and spectral flux profile size mismatch")

  const Size K = spectral_flux_profile.extent(0);

  flux_profile.resize(K);
  flux_profile = 0.0;

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size k = 0; k < K; k++) {
    auto s  = spectral_flux_profile[k];
    auto& y = flux_profile[k];
    for (Size i = 0; i < frequency_grid.size() - 1; i++) {
      const auto w  = 0.5 * (frequency_grid[i + 1] - frequency_grid[i]);
      y            += w * (s[i] + s[i + 1]);
    }
  }
}

void nlte_line_flux_profileIntegrate(
    QuantumIdentifierVectorMap& nlte_line_flux_profile,
    const Matrix& spectral_flux_profile,
    const AbsorptionBands& abs_bands,
    const ArrayOfAtmPoint& ray_path_atm_point,
    const AscendingGrid& frequency_grid) {
  ARTS_TIME_REPORT

  const Size K = spectral_flux_profile.extent(0);
  const Size M = spectral_flux_profile.extent(1);

  ARTS_USER_ERROR_IF(frequency_grid.size() != M,
                     "Frequency grid and spectral flux profile size mismatch")

  ARTS_USER_ERROR_IF(
      ray_path_atm_point.size() != K,
      "Atmospheric point and spectral flux profile size mismatch");

  nlte_line_flux_profile.clear();
  for (const auto& [key, band] : abs_bands) {
    ARTS_USER_ERROR_IF(band.size() != 1, "Only one line per band is supported");

    ARTS_USER_ERROR_IF(not is_voigt(band.lineshape),
                       "Only one line per band is supported");

    Matrix weighted_spectral_flux_profile(spectral_flux_profile.shape(), 0.0);

    for (Size k = 0; k < K; k++) {
      lbl::compute_voigt(weighted_spectral_flux_profile[k],
                         band.lines.front(),
                         frequency_grid,
                         ray_path_atm_point[k],
                         key.isot.mass);
    }

    weighted_spectral_flux_profile *= spectral_flux_profile;

    flux_profileIntegrate(nlte_line_flux_profile[key],
                          weighted_spectral_flux_profile,
                          frequency_grid);
  }
}

void spectral_flux_profileFromSpectralRadianceField(
    Matrix& spectral_flux_profile,
    const GriddedSpectralField6& spectral_radiance_field,
    const Stokvec& pol) {
  ARTS_TIME_REPORT

  using Constant::pi;
  using Conversion::cosd;
  using Conversion::deg2rad;

  const Size NALT = spectral_radiance_field.grid<0>().size();
  const Size NLAT = spectral_radiance_field.grid<1>().size();
  const Size NLON = spectral_radiance_field.grid<2>().size();
  const Size NZA  = spectral_radiance_field.grid<3>().size();
  const Size NAA  = spectral_radiance_field.grid<4>().size();
  const Size NFRE = spectral_radiance_field.grid<5>().size();

  const auto& za = spectral_radiance_field.grid<3>();
  const auto& aa = spectral_radiance_field.grid<4>();

  ARTS_USER_ERROR_IF(not spectral_radiance_field.ok(),
                     R"(Spectral radiance field is not OK (shape is wrong):

spectral_radiance_field.data.shape(): {:B,}
Should be:                            [NALT, NLAT, NLON, NZA, NAA, NFRE]

NALT: {}
NLAT: {}
NLON: {}
NZA:  {}
NAA:  {}
NFRE: {}
)",
                     spectral_radiance_field.data.shape(),
                     NALT,
                     NLAT,
                     NLON,
                     NZA,
                     NAA,
                     NFRE);
  ARTS_USER_ERROR_IF(NLAT != 1, "Only for one latitude point");
  ARTS_USER_ERROR_IF(NLON != 1, "Only for one longitude point");
  ARTS_USER_ERROR_IF(NZA < 2, "Must have more than one zenith angle.");
  ARTS_USER_ERROR_IF(NAA != 1, "Only for one azimuth angle.")

  spectral_flux_profile.resize(NALT, NFRE);
  spectral_flux_profile = 0.0;

  if (NAA == 1) {
    const auto&& s = spectral_radiance_field.data.view_as(NALT, NZA, NFRE);

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
    for (Size i = 0; i < NALT; i++) {
      for (Size j = 0; j < NFRE; j++) {
        for (Size iza0 = 0; iza0 < NZA - 1; iza0++) {
          const Size iza1   = iza0 + 1;
          const Numeric wza = 0.5 * pi * (cosd(za[iza0]) - cosd(za[iza1]));

          spectral_flux_profile[i, j] +=
              wza * (dot(s[i, iza0, j], pol) + dot(s[i, iza1, j], pol));
        }
      }
    }
  } else {
    const auto&& s = spectral_radiance_field.data.view_as(NALT, NZA, NAA, NFRE);

#pragma omp parallel for if (not arts_omp_in_parallel()) collapse(2)
    for (Size i = 0; i < NALT; i++) {
      for (Size j = 0; j < NFRE; j++) {
        for (Size iza0 = 0; iza0 < NZA - 1; iza0++) {
          const Size iza1   = iza0 + 1;
          const Numeric wza = 0.25 * (cosd(za[iza0]) - cosd(za[iza1]));

          for (Size iaa0 = 0; iaa0 < NAA - 1; iaa0++) {
            const Size iaa1 = iaa0 + 1;
            const Numeric w = wza * deg2rad(aa[iaa0] - aa[iaa1]);

            spectral_flux_profile[i, j] +=
                w *
                (dot(s[i, iza0, iaa0, j], pol) + dot(s[i, iza1, iaa0, j], pol) +
                 dot(s[i, iza0, iaa1, j], pol) + dot(s[i, iza1, iaa1, j], pol));
          }
        }
      }
    }
  }
}

#include <arts_omp.h>
#include <workspace.h>

void ray_path_spectral_radianceStepByStepEmissionForwardOnly(
    ArrayOfStokvecVector& ray_path_spectral_radiance,
    const ArrayOfMuelmatVector& ray_path_transmission_matrix,
    const ArrayOfStokvecVector& ray_path_spectral_radiance_source,
    const StokvecVector& spectral_radiance_background) try {
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
    const AtmField& atmospheric_field,
    const AscendingGrid& frequency_grid,
    const Agenda& propagation_matrix_agenda,
    const ArrayOfPropagationPathPoint& ray_path,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda,
    const SurfaceField& surface_field) try {
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
      spectral_radiance_space_agenda,
      spectral_radiance_surface_agenda);
  ArrayOfAtmPoint ray_path_atmospheric_point;
  ray_path_atmospheric_pointFromPath(
      ray_path_atmospheric_point, ray_path, atmospheric_field);
  ArrayOfAscendingGrid ray_path_frequency_grid;
  ArrayOfVector3 ray_path_frequency_grid_wind_shift_jacobian;
  ray_path_frequency_gridFromPath(ray_path_frequency_grid,
                                  ray_path_frequency_grid_wind_shift_jacobian,
                                  frequency_grid,
                                  ray_path,
                                  ray_path_atmospheric_point);
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
      ray_path_frequency_grid_wind_shift_jacobian,
      {},
      ray_path,
      ray_path_atmospheric_point);
  ArrayOfMuelmatVector ray_path_transmission_matrix;
  ArrayOfMuelmatTensor3 ray_path_transmission_matrix_jacobian;
  ray_path_transmission_matrixFromPath(ray_path_transmission_matrix,
                                       ray_path_transmission_matrix_jacobian,
                                       ray_path_propagation_matrix,
                                       ray_path_propagation_matrix_jacobian,
                                       ray_path,
                                       ray_path_atmospheric_point,
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
      ray_path_atmospheric_point,
      {});
  ray_path_spectral_radianceStepByStepEmissionForwardOnly(
      ray_path_spectral_radiance,
      ray_path_transmission_matrix,
      ray_path_spectral_radiance_source,
      spectral_radiance_background);
}
ARTS_METHOD_ERROR_CATCH

void spectral_flux_profileFromPathField(
    const Workspace& ws,
    Matrix& spectral_flux_profile,
    const ArrayOfArrayOfPropagationPathPoint& ray_path_field,
    const AtmField& atmospheric_field,
    const Agenda& propagation_matrix_agenda,
    const Agenda& spectral_radiance_space_agenda,
    const Agenda& spectral_radiance_surface_agenda,
    const SurfaceField& surface_field,
    const AscendingGrid& frequency_grid,
    const AscendingGrid& altitude_grid) try {
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
          atmospheric_field,
          frequency_grid,
          propagation_matrix_agenda,
          ray_path_field[n],
          spectral_radiance_space_agenda,
          spectral_radiance_surface_agenda,
          surface_field);
    } catch (std::exception& e) {
#pragma omp critical
      if (error.empty()) error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(), "{}", error)

  for (Size n = 0; n < N; n++) {
    ARTS_USER_ERROR_IF(
        ray_path_spectral_radiance_field[n].size() != ray_path_field[n].size(),
        "Not all ray paths have the same altitude count")
  }

  ARTS_USER_ERROR_IF(
      stdr::any_of(ray_path_spectral_radiance_field | stdv::join,
                   [K](const StokvecVector& v) { return v.size() != K; }),
      "Not all ray path spectral radiances in the field have the same frequency count")

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

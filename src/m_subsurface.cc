#include <arts_omp.h>
#include <workspace.h>

#include <algorithm>
#include <exception>
#include <ranges>

void subsurface_profileFromPath(ArrayOfSubsurfacePoint& subsurface_profile,
                                const SubsurfaceField& subsurface_field,
                                const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  const auto to_vec3 = stdv::transform([](const auto& p) { return p.pos; });

  const auto extract = stdv::transform(
      [&](Vector3 pos) { return subsurface_field.at(pos[0], pos[1], pos[2]); });

  subsurface_profile =
      ray_path | to_vec3 | extract | stdr::to<ArrayOfSubsurfacePoint>();
}

void ray_pathFromPointAndDepth(ArrayOfPropagationPathPoint& ray_path,
                               const PropagationPathPoint& ray_path_point,
                               const DescendingGrid& depth_profile) {
  ARTS_TIME_REPORT

  const auto to_ppp = stdv::transform([ray_path_point](const auto& d) {
    auto v    = ray_path_point;
    v.pos[0] += d;
    return v;
  });

  ray_path = depth_profile | to_ppp | stdr::to<ArrayOfPropagationPathPoint>();
}

void spectral_radianceSubsurfaceDisortEmissionWithJacobian(
    const Workspace& ws,
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field_,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension,
    const Index& disort_legendre_polynomial_dimension,
    const Agenda& disort_settings_agenda,
    const Agenda& disort_settings_downwelling_wrapper_agenda,
    const DescendingGrid& depth_profile) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surface_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surface_field.ellipsoid)

  DisortSettings disort_settings                = {};
  ArrayOfPropagationPathPoint ray_path          = {};
  DisortRadiance disort_spectral_radiance_field = {};
  ZenithGriddedField1 disort_quadrature         = {};
  Vector model_state_vector                     = {};
  Vector model_state_vector_perturbation        = {};
  StokvecVector spectral_radiance2              = {};
  AtmField atm_field                            = atm_field_;
  SurfaceField surf_field                       = surface_field;
  SubsurfaceField subsurf_field                 = subsurface_field;
  const AzimuthGrid azimuth_grid = Vector{ray_path_point.azimuth()};

  spectral_radiance_jacobian.resize(jacobian_targets.x_size(),
                                    freq_grid.size());

  spectral_radianceSubsurfaceDisortEmission(
      ws,
      spectral_radiance,
      disort_settings,
      ray_path,
      disort_spectral_radiance_field,
      disort_quadrature,
      atm_field,
      disort_fourier_mode_dimension,
      disort_legendre_polynomial_dimension,
      disort_quadrature_dimension,
      disort_settings_agenda,
      disort_settings_downwelling_wrapper_agenda,
      freq_grid,
      ray_path_point,
      subsurface_field,
      surface_field,
      depth_profile,
      azimuth_grid);

  model_state_vectorPerturbations(model_state_vector_perturbation,
                                  jacobian_targets);
  model_state_vectorInit(model_state_vector, jacobian_targets);
  model_state_vectorFromAtmosphere(
      model_state_vector, atm_field, jacobian_targets);
  model_state_vectorFromSurface(
      model_state_vector, surface_field, jacobian_targets);
  model_state_vectorFromSubsurface(
      model_state_vector, subsurface_field, jacobian_targets);

  ARTS_USER_ERROR_IF(
      model_state_vector.size() != jacobian_targets.x_size() or
          model_state_vector.size() != model_state_vector_perturbation.size(),
      "Model state vector, model state vector perturbations, "
      "and Jacobian targets size do not match: {} vs {} vs {}",
      model_state_vector.size(),
      model_state_vector_perturbation.size(),
      jacobian_targets.x_size());

  String error{};

#pragma omp parallel for if (not arts_omp_in_parallel() and             \
                                 static_cast<Size>(                     \
                                         disort_quadrature_dimension) < \
                                         model_state_vector.size())     \
    firstprivate(model_state_vector,                                    \
                     atm_field,                                         \
                     surf_field,                                        \
                     subsurf_field,                                     \
                     disort_settings,                                   \
                     ray_path,                                          \
                     disort_spectral_radiance_field,                    \
                     disort_quadrature,                                 \
                     spectral_radiance2)
  for (Size i = 0; i < model_state_vector.size(); i++) {
    try {
      const Numeric orig     = model_state_vector[i];
      model_state_vector[i] += model_state_vector_perturbation[i];

      surface_fieldFromModelState(
          surf_field, model_state_vector, jacobian_targets);
      subsurface_fieldFromModelState(
          subsurf_field, model_state_vector, jacobian_targets);
      atm_fieldFromModelState(atm_field, model_state_vector, jacobian_targets);

      model_state_vector[i] = orig;

      spectral_radianceSubsurfaceDisortEmission(
          ws,
          spectral_radiance2,
          disort_settings,
          ray_path,
          disort_spectral_radiance_field,
          disort_quadrature,
          atm_field,
          disort_fourier_mode_dimension,
          disort_legendre_polynomial_dimension,
          disort_quadrature_dimension,
          disort_settings_agenda,
          disort_settings_downwelling_wrapper_agenda,
          freq_grid,
          ray_path_point,
          subsurf_field,
          surf_field,
          depth_profile,
          azimuth_grid);

      std::transform(spectral_radiance2.begin(),
                     spectral_radiance2.end(),
                     spectral_radiance.begin(),
                     spectral_radiance_jacobian[i].begin(),
                     [d = model_state_vector_perturbation[i]](
                         const Stokvec& a, const Stokvec& b) -> Stokvec {
                       return (a - b) / d;
                     });
    } catch (std::exception& e) {
#pragma omp critical
      error = e.what();
    }
  }

  ARTS_USER_ERROR_IF(not error.empty(),
                     "Error occurred in disort-spectral with jacobian:\n{}",
                     error)
}

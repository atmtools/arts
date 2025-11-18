#include <arts_omp.h>
#include <workspace.h>

#include <algorithm>
#include <exception>
#include <ranges>

void subsurf_profileFromPath(ArrayOfSubsurfacePoint& subsurf_profile,
                             const SubsurfaceField& subsurf_field,
                             const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  const auto to_vec3 = stdv::transform([](const auto& p) { return p.pos; });

  const auto extract = stdv::transform(
      [&](Vector3 pos) { return subsurf_field.at(pos[0], pos[1], pos[2]); });

  subsurf_profile =
      ray_path | to_vec3 | extract | stdr::to<ArrayOfSubsurfacePoint>();
}

void ray_pathFromPointAndDepth(ArrayOfPropagationPathPoint& ray_path,
                               const PropagationPathPoint& ray_point,
                               const DescendingGrid& depth_profile) {
  ARTS_TIME_REPORT

  const auto to_ppp = stdv::transform([ray_point](const auto& d) {
    auto v    = ray_point;
    v.pos[0] += d;
    return v;
  });

  ray_path = depth_profile | to_ppp | stdr::to<ArrayOfPropagationPathPoint>();
}

void spectral_radSubsurfaceDisortEmissionWithJacobian(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field_,
    const SurfaceField& surf_field_,
    const SubsurfaceField& subsurf_field_,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_point,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension,
    const Index& disort_legendre_polynomial_dimension,
    const Agenda& disort_settings_agenda,
    const Agenda& disort_settings_downwelling_wrapper_agenda,
    const DescendingGrid& depth_profile) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field_.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field_.ellipsoid)

  DisortSettings disort_settings           = {};
  ArrayOfPropagationPathPoint ray_path     = {};
  DisortRadiance disort_spectral_rad_field = {};
  ZenGriddedField1 disort_quadrature       = {};
  Vector model_state_vec                   = {};
  Vector model_state_vec_perturbation      = {};
  StokvecVector spectral_rad2              = {};
  AtmField atm_field                       = atm_field_;
  SurfaceField surf_field                  = surf_field_;
  SubsurfaceField subsurf_field            = subsurf_field_;
  const AziGrid azi_grid                   = Vector{ray_point.azimuth()};

  spectral_rad_jac.resize(jac_targets.x_size(), freq_grid.size());

  spectral_radSubsurfaceDisortEmission(
      ws,
      spectral_rad,
      disort_settings,
      ray_path,
      disort_spectral_rad_field,
      disort_quadrature,
      atm_field,
      disort_fourier_mode_dimension,
      disort_legendre_polynomial_dimension,
      disort_quadrature_dimension,
      disort_settings_agenda,
      disort_settings_downwelling_wrapper_agenda,
      freq_grid,
      ray_point,
      subsurf_field,
      surf_field,
      depth_profile,
      azi_grid);

  model_state_vecPerturbations(model_state_vec_perturbation, jac_targets);
  model_state_vecInit(model_state_vec, jac_targets);
  model_state_vecFromAtmosphere(model_state_vec, atm_field, jac_targets);
  model_state_vecFromSurface(model_state_vec, surf_field, jac_targets);
  model_state_vecFromSubsurface(model_state_vec, subsurf_field, jac_targets);

  ARTS_USER_ERROR_IF(
      model_state_vec.size() != jac_targets.x_size() or
          model_state_vec.size() != model_state_vec_perturbation.size(),
      "Model state vector, model state vector perturbations, "
      "and Jacobian targets size do not match: {} vs {} vs {}",
      model_state_vec.size(),
      model_state_vec_perturbation.size(),
      jac_targets.x_size());

  String error{};

#pragma omp parallel for if (not arts_omp_in_parallel() and             \
                                 static_cast<Size>(                     \
                                         disort_quadrature_dimension) < \
                                         model_state_vec.size())        \
    firstprivate(model_state_vec,                                       \
                     atm_field,                                         \
                     surf_field,                                        \
                     subsurf_field,                                     \
                     disort_settings,                                   \
                     ray_path,                                          \
                     disort_spectral_rad_field,                         \
                     disort_quadrature,                                 \
                     spectral_rad2)
  for (Size i = 0; i < model_state_vec.size(); i++) {
    try {
      const Numeric orig  = model_state_vec[i];
      model_state_vec[i] += model_state_vec_perturbation[i];

      surf_fieldFromModelState(surf_field, model_state_vec, jac_targets);
      subsurf_fieldFromModelState(subsurf_field, model_state_vec, jac_targets);
      atm_fieldFromModelState(atm_field, model_state_vec, jac_targets);

      model_state_vec[i] = orig;

      spectral_radSubsurfaceDisortEmission(
          ws,
          spectral_rad2,
          disort_settings,
          ray_path,
          disort_spectral_rad_field,
          disort_quadrature,
          atm_field,
          disort_fourier_mode_dimension,
          disort_legendre_polynomial_dimension,
          disort_quadrature_dimension,
          disort_settings_agenda,
          disort_settings_downwelling_wrapper_agenda,
          freq_grid,
          ray_point,
          subsurf_field,
          surf_field,
          depth_profile,
          azi_grid);

      std::transform(spectral_rad2.begin(),
                     spectral_rad2.end(),
                     spectral_rad.begin(),
                     spectral_rad_jac[i].begin(),
                     [d = model_state_vec_perturbation[i]](
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

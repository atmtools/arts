#include <workspace.h>

#include <ranges>

void subsurface_profileFromPath(ArrayOfSubsurfacePoint subsurface_profile,
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

void spectral_radianceScalarSubsurface(
    const Workspace& ws,
    const AscendingGrid& frequency_grid,
    const PropagationPathPoint& ray_path_point,
    const Index& disort_quadrature_dimension,
    const Index& disort_fourier_mode_dimension,
    const Index& disort_legendre_polynomial_dimension,
    const Agenda& disort_settings_agenda,
    const DescendingGrid& depth_profile,
    const AscendingGrid& phis) {
  ARTS_TIME_REPORT

  ArrayOfPropagationPathPoint ray_path;
  ray_pathFromPointAndDepth(ray_path, ray_path_point, depth_profile);

  Tensor4 disort_spectral_radiance_field;
  Vector disort_quadrature_angles, disort_quadrature_weights;
  disort_spectral_radiance_fieldFromAgenda(ws,
                                           disort_spectral_radiance_field,
                                           disort_quadrature_angles,
                                           disort_quadrature_weights,
                                           disort_fourier_mode_dimension,
                                           disort_legendre_polynomial_dimension,
                                           disort_quadrature_dimension,
                                           disort_settings_agenda,
                                           frequency_grid,
                                           ray_path,
                                           phis);
}

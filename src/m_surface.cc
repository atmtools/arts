#include <arts_omp.h>
#include <workspace.h>

namespace {
Vector2 specular_losNormal(const Vector2& normal,
                           const Vector2& los,
                           const Vector3& pos,
                           const Vector2& ell) {
  const auto [ecef_pos, ecef_los] = path::geodetic_poslos2ecef(pos, los, ell);
  const auto [_, ecef_normal] = path::geodetic_poslos2ecef(pos, normal, ell);

  // Specular direction is 2(dn*di)dn-di, where dn is the normal vector
  const Numeric fac       = 2 * dot(ecef_normal, ecef_los);
  const Vector3 ecef_spec = {fac * ecef_normal[0] - ecef_los[0],
                             fac * ecef_normal[1] - ecef_los[1],
                             fac * ecef_normal[2] - ecef_los[2]};

  return path::ecef2geodetic_poslos(ecef_pos, normalized(ecef_spec), ell)
      .second;
}
}  // namespace

void surface_reflectanceFlatRealFresnel(
    MuelmatVector& surface_reflectance,
    MuelmatMatrix& surface_reflectance_jacobian,
    const AscendingGrid& frequency_grid,
    const SurfaceField& surface_field,
    const PropagationPathPoint& ray_path_point,
    const JacobianTargets& jacobian_targets) try {
  ARTS_TIME_REPORT

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag refraction_target{"scalar refractive index"};

  ARTS_USER_ERROR_IF(not surface_field.contains(refraction_target),
                     R"--(Missing key property tag for method.

Tag "scalar refractive index" not in the surface field.

surface_field:
{}
)--",
                     surface_field);

  const auto& pos                  = ray_path_point.pos;
  const auto& los                  = ray_path_point.los;
  const auto& ell                  = surface_field.ellipsoid;
  const Numeric lat                = pos[1];
  const Numeric lon                = pos[2];
  const SurfacePoint surface_point = surface_field.at(lat, lon);
  const Numeric n1                 = ray_path_point.nreal;
  const Numeric n2                 = surface_point[refraction_target];
  const Size nf                    = frequency_grid.size();

  surface_reflectance.resize(nf);
  surface_reflectance_jacobian.resize(jacobian_targets.target_count(), nf);

  const auto [Rv, Rh] = fresnel(
      n1,
      n2,
      std::acos(dot(
          path::geodetic_poslos2ecef(pos, los, ell).second,
          path::geodetic_poslos2ecef(pos, surface_point.normal, ell).second)));

  const Muelmat R              = rtepack::fresnel_reflectance(Rv, Rh);
  surface_reflectance          = R;
  surface_reflectance_jacobian = 0.0;

  for (auto& target : jacobian_targets.surf) {
    if (target.type == refraction_target) {
      const auto [Rv2, Rh2] =
          fresnel(n1,
                  n2 + 1e-3,
                  std::acos(dot(
                      path::geodetic_poslos2ecef(pos, los, ell).second,
                      path::geodetic_poslos2ecef(pos, surface_point.normal, ell)
                          .second)));

      const Muelmat dR = 1000. * (rtepack::fresnel_reflectance(Rv2, Rh2) - R);

      surface_reflectance_jacobian[target.target_pos] = dR;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void surface_reflectanceFlatScalar(
    MuelmatVector& surface_reflectance,
    MuelmatMatrix& surface_reflectance_jacobian,
    const AscendingGrid& frequency_grid,
    const SurfaceField& surface_field,
    const PropagationPathPoint& ray_path_point,
    const JacobianTargets& jacobian_targets) try {
  ARTS_TIME_REPORT

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag reflectance_target{"flat scalar reflectance"};

  ARTS_USER_ERROR_IF(not surface_field.contains(reflectance_target),
                     R"--(Missing key property tag for method.

Tag "flat scalar reflectance" not in the surface field.

surface_field:
{}
)--",
                     surface_field);

  const Numeric lat                = ray_path_point.pos[1];
  const Numeric lon                = ray_path_point.pos[2];
  const SurfacePoint surface_point = surface_field.at(lat, lon);
  const Numeric R                  = surface_point[reflectance_target];
  const Size nf                    = frequency_grid.size();

  ARTS_USER_ERROR_IF(
      R < 0.0 or R > 1.0,
      "Flat scalar reflectance must be between 0 and 1, but is {}.",
      R)

  surface_reflectance.resize(nf);
  surface_reflectance_jacobian.resize(jacobian_targets.target_count(), nf);

  surface_reflectance          = R;
  surface_reflectance_jacobian = 0.0;

  for (auto& target : jacobian_targets.surf) {
    if (target.type == reflectance_target) {
      surface_reflectance_jacobian[target.target_pos] = 1.0;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceSurfaceReflectance(
    const Workspace& ws,
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& frequency_grid,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const Agenda& spectral_radiance_observer_agenda,
    const Agenda& spectral_radiance_closed_surface_agenda,
    const Agenda& surface_reflectance_agenda) try {
  ARTS_TIME_REPORT

  const Size NF                    = frequency_grid.size();
  const Size NX                    = jacobian_targets.x_size();
  const Numeric lat                = ray_path_point.pos[1];
  const Numeric lon                = ray_path_point.pos[2];
  const SurfacePoint surface_point = surface_field.at(lat, lon);

  MuelmatVector surface_reflectance;
  MuelmatMatrix surface_reflectance_jacobian;
  surface_reflectance_agendaExecute(ws,
                                    surface_reflectance,
                                    surface_reflectance_jacobian,
                                    frequency_grid,
                                    surface_field,
                                    ray_path_point,
                                    jacobian_targets,
                                    surface_reflectance_agenda);

  // Get the direction of the incoming radiation
  const Vector2 los = specular_losNormal(surface_point.normal,
                                         ray_path_point.los,
                                         ray_path_point.pos,
                                         surface_field.ellipsoid);

  ArrayOfPropagationPathPoint ray_path;
  StokvecVector spectral_radiance_surface;
  StokvecMatrix spectral_radiance_jacobian_surface;

  spectral_radiance_observer_agendaExecute(ws,
                                           spectral_radiance,
                                           spectral_radiance_jacobian,
                                           ray_path,
                                           frequency_grid,
                                           jacobian_targets,
                                           ray_path_point.pos,
                                           los,
                                           atmospheric_field,
                                           surface_field,
                                           subsurface_field,
                                           spectral_radiance_observer_agenda);

  spectral_radiance_surface_agendaExecute(
      ws,
      spectral_radiance_surface,
      spectral_radiance_jacobian_surface,
      frequency_grid,
      jacobian_targets,
      ray_path_point,
      surface_field,
      subsurface_field,
      spectral_radiance_closed_surface_agenda);

#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size j = 0; j < NX; j++) {
    for (Size i = 0; i < NF; i++) {
      spectral_radiance_jacobian[j, i] =
          rtepack::reflection(spectral_radiance_jacobian[j, i],
                              surface_reflectance[i],
                              spectral_radiance_jacobian_surface[j, i]);
    }
  }

  for (auto& target : jacobian_targets.surf) {
    const SurfaceData& data = surface_field[target.type];
    const auto ws           = data.flat_weights(lat, lon);

#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size i = 0; i < NF; i++) {
      for (const auto& [j, w] : ws) {
        spectral_radiance_jacobian[j + target.x_start, i] +=
            w * rtepack::dreflection(
                    spectral_radiance[i],
                    surface_reflectance_jacobian[target.target_pos, i],
                    spectral_radiance_surface[i]);
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < NF; i++) {
    spectral_radiance[i] = rtepack::reflection(spectral_radiance[i],
                                               surface_reflectance[i],
                                               spectral_radiance_surface[i]);
  }
}
ARTS_METHOD_ERROR_CATCH

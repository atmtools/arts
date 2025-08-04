#include <arts_omp.h>
#include <workspace.h>

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

void spectral_radianceFlatScalarReflectance(
    const Workspace& ws,
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& frequency_grid,
    const AtmField& atmospheric_field,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const Agenda& spectral_radiance_observer_agenda) try {
  ARTS_TIME_REPORT

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag reflectance_target{"flat scalar reflectance"};

  ARTS_USER_ERROR_IF(not surface_field.contains(reflectance_target) or
                         not surface_field.contains(SurfaceKey::t),
                     R"--(Missing key property tag for method.

Tag "flat scalar reflectance" not in the surface field.

surface_field:
{}
)--",
                     surface_field);

  const SurfacePoint surface_point =
      surface_field.at(ray_path_point.pos[1], ray_path_point.pos[2]);

  // Get the direction of the incoming radiation
  const Vector2 los = specular_losNormal(surface_point.normal,
                                         ray_path_point.los,
                                         ray_path_point.pos,
                                         surface_field.ellipsoid);

  const Size NF = frequency_grid.size();
  const Size NX = jacobian_targets.x_size();

  ArrayOfPropagationPathPoint ray_path;
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

  const Numeric R = surface_point[reflectance_target];
  const Numeric T = surface_point[SurfaceKey::t];

  ARTS_USER_ERROR_IF(
      R < 0.0 or R > 1.0,
      "Flat scalar reflectance must be between 0 and 1, but is {}.",
      R)

#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size j = 0; j < NX; j++) {
    for (Size i = 0; i < NF; i++) {
      spectral_radiance_jacobian[j, i] =
          rtepack::flat_scalar_reflection_pure_reflect(
              spectral_radiance_jacobian[j, i], R);
    }
  }

  for (auto& target : jacobian_targets.surf) {
    if (target.type == SurfaceKey::t) {
      const auto& data = surface_field[target.type];

      const auto ws =
          data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);

#pragma omp parallel for if (not arts_omp_in_parallel())
      for (Size i = 0; i < NF; i++) {
        const Numeric dBdt = dplanck_dt(frequency_grid[i], T);
        for (const auto& w : ws) {
          spectral_radiance_jacobian[w.first + target.x_start, i] +=
              w.second * rtepack::dflat_scalar_reflection_db(R, dBdt);
        }
      }
    } else if (target.type == reflectance_target) {
      const auto& data = surface_field[target.type];

      const auto ws =
          data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);

#pragma omp parallel for if (not arts_omp_in_parallel())
      for (Size i = 0; i < NF; i++) {
        const Numeric B = planck(frequency_grid[i], T);
        for (const auto& w : ws) {
          spectral_radiance_jacobian[w.first + target.x_start, i] +=
              w.second *
              rtepack::dflat_scalar_reflection_dr(spectral_radiance[i], 1.0, B);
        }
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < NF; i++) {
    const Numeric B = planck(frequency_grid[i], T);
    spectral_radiance[i] =
        rtepack::flat_scalar_reflection(spectral_radiance[i], R, B);
  }
}
ARTS_METHOD_ERROR_CATCH

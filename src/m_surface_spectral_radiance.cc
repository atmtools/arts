#include <arts_omp.h>
#include <nonstd.h>
#include <workspace.h>

#include "path_point.h"
#include "physics_funcs.h"

namespace {
constexpr bool is_close_enough(Numeric a, Numeric b, Numeric lim) {
  return nonstd::abs(a - b) <= lim;
}

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
    const Agenda& spectral_radiance_observer_agenda,
    const Agenda& spectral_radiance_subsurface_agenda) try {
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

  const SurfacePoint surface_point =
      surface_field.at(ray_path_point.pos[1], ray_path_point.pos[2]);

  const Numeric R = surface_point[reflectance_target];
  const Numeric z = surface_point.elevation;

  ARTS_USER_ERROR_IF(
      not is_close_enough(z, ray_path_point.pos[0], 1e-2),
      R"(The surface elevation {} does not match the ray path point altitude {}.
)",
      z,
      ray_path_point.pos[0])

  ARTS_USER_ERROR_IF(
      R < 0.0 or R > 1.0,
      "Flat scalar reflectance must be between 0 and 1, but is {}.",
      R)

  // Get the direction of the incoming radiation
  const Vector2 los = specular_losNormal(surface_point.normal,
                                         ray_path_point.los,
                                         ray_path_point.pos,
                                         surface_field.ellipsoid);

  ArrayOfPropagationPathPoint ray_path;
  StokvecVector spectral_radiance_subsurface;
  StokvecMatrix spectral_radiance_jacobian_subsurface;

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

  spectral_radiance_surface_agendaExecute(ws,
                                          spectral_radiance_subsurface,
                                          spectral_radiance_jacobian_subsurface,
                                          frequency_grid,
                                          jacobian_targets,
                                          ray_path_point,
                                          surface_field,
                                          subsurface_field,
                                          spectral_radiance_subsurface_agenda);

  const Size NF = frequency_grid.size();
  const Size NX = jacobian_targets.x_size();

#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size j = 0; j < NX; j++) {
    for (Size i = 0; i < NF; i++) {
      spectral_radiance_jacobian[j, i] = rtepack::flat_scalar_reflection(
          spectral_radiance_jacobian[j, i],
          R,
          spectral_radiance_jacobian_subsurface[j, i]);
    }
  }

  for (auto& target : jacobian_targets.surf) {
    if (target.type == reflectance_target) {
      const auto& data = surface_field[target.type];

      const auto ws =
          data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);

#pragma omp parallel for if (not arts_omp_in_parallel())
      for (Size i = 0; i < NF; i++) {
        for (const auto& w : ws) {
          spectral_radiance_jacobian[w.first + target.x_start, i] +=
              w.second *
              rtepack::dflat_scalar_reflection_dr(
                  spectral_radiance[i], spectral_radiance_subsurface[i]);
        }
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < NF; i++) {
    spectral_radiance[i] = rtepack::flat_scalar_reflection(
        spectral_radiance[i], R, spectral_radiance_subsurface[i]);
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radianceFlatFresnelReflectance(
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
    const Agenda& spectral_radiance_subsurface_agenda) try {
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

  const SurfacePoint surface_point =
      surface_field.at(ray_path_point.pos[1], ray_path_point.pos[2]);

  const Numeric n1 = ray_path_point.nreal;
  const Numeric n2 = surface_point[refraction_target];
  const Numeric z  = surface_point.elevation;

  ARTS_USER_ERROR_IF(
      not is_close_enough(z, ray_path_point.pos[0], 1e-2),
      R"(The surface elevation {} does not match the ray path point altitude {}.
)",
      z,
      ray_path_point.pos[0])

  // Get the direction of the incoming radiation
  const Vector2 los = specular_losNormal(surface_point.normal,
                                         ray_path_point.los,
                                         ray_path_point.pos,
                                         surface_field.ellipsoid);

  ArrayOfPropagationPathPoint ray_path;
  StokvecVector spectral_radiance_subsurface;
  StokvecMatrix spectral_radiance_jacobian_subsurface;

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

  spectral_radiance_surface_agendaExecute(ws,
                                          spectral_radiance_subsurface,
                                          spectral_radiance_jacobian_subsurface,
                                          frequency_grid,
                                          jacobian_targets,
                                          ray_path_point,
                                          surface_field,
                                          subsurface_field,
                                          spectral_radiance_subsurface_agenda);

  const Size NF = frequency_grid.size();
  const Size NX = jacobian_targets.x_size();

  const auto [Rv, Rh] = fresnel(
      n1,
      n2,
      std::acos(dot(
          path::geodetic_poslos2ecef(
              ray_path_point.pos, los, surface_field.ellipsoid)
              .second,
          path::geodetic_poslos2ecef(
              ray_path_point.pos, surface_point.normal, surface_field.ellipsoid)
              .second)));

  const Muelmat R = rtepack::fresnel_reflectance(Rv, Rh);

#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size j = 0; j < NX; j++) {
    for (Size i = 0; i < NF; i++) {
      spectral_radiance_jacobian[j, i] =
          rtepack::reflection(spectral_radiance_jacobian[j, i],
                              R,
                              spectral_radiance_jacobian_subsurface[j, i]);
    }
  }

  for (auto& target : jacobian_targets.surf) {
    if (target.type == refraction_target) {
      //! FIXME: This should be done properly, but I am not going to spend time getting that done (RL)
      const auto [Rv2, Rh2] = fresnel(
          n1,
          n2 + 1e-3,
          std::acos(dot(path::geodetic_poslos2ecef(
                            ray_path_point.pos, los, surface_field.ellipsoid)
                            .second,
                        path::geodetic_poslos2ecef(ray_path_point.pos,
                                                   surface_point.normal,
                                                   surface_field.ellipsoid)
                            .second)));
      const Muelmat dR = 1000. * (rtepack::fresnel_reflectance(Rv2, Rh2) - R);

      const auto& data = surface_field[target.type];

      const auto ws =
          data.flat_weights(ray_path_point.pos[1], ray_path_point.pos[2]);

#pragma omp parallel for if (not arts_omp_in_parallel())
      for (Size i = 0; i < NF; i++) {
        for (const auto& w : ws) {
          spectral_radiance_jacobian[w.first + target.x_start, i] +=
              w.second *
              rtepack::dreflection_dn2(
                  spectral_radiance[i], dR, spectral_radiance_subsurface[i]);
        }
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < NF; i++) {
    spectral_radiance[i] = rtepack::reflection(
        spectral_radiance[i], R, spectral_radiance_subsurface[i]);
  }
}
ARTS_METHOD_ERROR_CATCH

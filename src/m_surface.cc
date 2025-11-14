#include <arts_omp.h>
#include <geodetic.h>
#include <workspace.h>

namespace {
Vector2 specular_losNormal(const Vector2& normal,
                           const Vector2& los,
                           const Vector3& pos,
                           const Vector2& ell) {
  const auto [ecef_pos, ecef_los] = geodetic_los2ecef(pos, los, ell);
  const auto [_, ecef_normal]     = geodetic_los2ecef(pos, normal, ell);

  // Specular direction is 2(dn*di)dn-di, where dn is the normal vector
  const Numeric fac       = 2 * dot(ecef_normal, ecef_los);
  const Vector3 ecef_spec = {fac * ecef_normal[0] - ecef_los[0],
                             fac * ecef_normal[1] - ecef_los[1],
                             fac * ecef_normal[2] - ecef_los[2]};

  return ecef2geodetic_los(ecef_pos, normalized(ecef_spec), ell).second;
}
}  // namespace

void spectral_surf_reflFlatRealFresnel(MuelmatVector& spectral_surf_refl,
                                       MuelmatMatrix& spectral_surf_refl_jac,
                                       const AscendingGrid& freq_grid,
                                       const SurfaceField& surf_field,
                                       const PropagationPathPoint& ray_point,
                                       const JacobianTargets& jac_targets) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag refraction_target{"scalar refractive index"};

  ARTS_USER_ERROR_IF(not surf_field.contains(refraction_target),
                     R"--(Missing key property tag for method.

Tag "scalar refractive index" not in the surface field.

surf_field:
{}
)--",
                     surf_field);

  const auto& pos               = ray_point.pos;
  const auto& los               = ray_point.los;
  const auto& ell               = surf_field.ellipsoid;
  const Numeric lat             = pos[1];
  const Numeric lon             = pos[2];
  const SurfacePoint surf_point = surf_field.at(lat, lon);
  const Numeric n1              = ray_point.nreal;
  const Numeric n2              = surf_point[refraction_target];
  const Size nf                 = freq_grid.size();

  spectral_surf_refl.resize(nf);
  spectral_surf_refl_jac.resize(jac_targets.target_count(), nf);

  const auto [Rv, Rh] = fresnel(
      n1,
      n2,
      std::acos(dot(geodetic_los2ecef(pos, los, ell).second,
                    geodetic_los2ecef(pos, surf_point.normal, ell).second)));

  const Muelmat R        = rtepack::fresnel_reflectance(Rv, Rh);
  spectral_surf_refl     = R;
  spectral_surf_refl_jac = 0.0;

  for (auto& target : jac_targets.surf) {
    if (target.type == refraction_target) {
      const auto [Rv2, Rh2] =
          fresnel(n1,
                  n2 + 1e-3,
                  std::acos(dot(
                      geodetic_los2ecef(pos, los, ell).second,
                      geodetic_los2ecef(pos, surf_point.normal, ell).second)));

      const Muelmat dR = 1000. * (rtepack::fresnel_reflectance(Rv2, Rh2) - R);

      spectral_surf_refl_jac[target.target_pos] = dR;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_surf_reflFlatScalar(MuelmatVector& spectral_surf_refl,
                                  MuelmatMatrix& spectral_surf_refl_jac,
                                  const AscendingGrid& freq_grid,
                                  const SurfaceField& surf_field,
                                  const PropagationPathPoint& ray_point,
                                  const JacobianTargets& jac_targets) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag reflectance_target{"flat scalar reflectance"};

  ARTS_USER_ERROR_IF(not surf_field.contains(reflectance_target),
                     R"--(Missing key property tag for method.

Tag "flat scalar reflectance" not in the surface field.

surf_field:
{}
)--",
                     surf_field);

  const Numeric lat             = ray_point.pos[1];
  const Numeric lon             = ray_point.pos[2];
  const SurfacePoint surf_point = surf_field.at(lat, lon);
  const Numeric R               = surf_point[reflectance_target];
  const Size nf                 = freq_grid.size();

  ARTS_USER_ERROR_IF(
      R < 0.0 or R > 1.0,
      "Flat scalar reflectance must be between 0 and 1, but is {}.",
      R)

  spectral_surf_refl.resize(nf);
  spectral_surf_refl_jac.resize(jac_targets.target_count(), nf);

  spectral_surf_refl     = R;
  spectral_surf_refl_jac = 0.0;

  for (auto& target : jac_targets.surf) {
    if (target.type == reflectance_target) {
      spectral_surf_refl_jac[target.target_pos] = 1.0;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSurfaceReflectance(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_point,
    const Agenda& spectral_rad_observer_agenda,
    const Agenda& spectral_rad_closed_surface_agenda,
    const Agenda& spectral_surf_refl_agenda) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      surf_field.bad_ellipsoid(),
      "Surface field not properly set up - bad reference ellipsoid: {:B,}",
      surf_field.ellipsoid)

  const Size NF                 = freq_grid.size();
  const Size NX                 = jac_targets.x_size();
  const Numeric lat             = ray_point.pos[1];
  const Numeric lon             = ray_point.pos[2];
  const SurfacePoint surf_point = surf_field.at(lat, lon);

  MuelmatVector spectral_surf_refl;
  MuelmatMatrix spectral_surf_refl_jac;
  spectral_surf_refl_agendaExecute(ws,
                                   spectral_surf_refl,
                                   spectral_surf_refl_jac,
                                   freq_grid,
                                   surf_field,
                                   ray_point,
                                   jac_targets,
                                   spectral_surf_refl_agenda);

  // Get the direction of the incoming radiation
  const Vector2 los = specular_losNormal(
      surf_point.normal, ray_point.los, ray_point.pos, surf_field.ellipsoid);

  ArrayOfPropagationPathPoint ray_path;
  StokvecVector spectral_rad_surface;
  StokvecMatrix spectral_rad_jac_surface;

  spectral_rad_observer_agendaExecute(ws,
                                      spectral_rad,
                                      spectral_rad_jac,
                                      ray_path,
                                      freq_grid,
                                      jac_targets,
                                      ray_point.pos,
                                      los,
                                      atm_field,
                                      surf_field,
                                      subsurf_field,
                                      spectral_rad_observer_agenda);

  spectral_rad_surface_agendaExecute(ws,
                                     spectral_rad_surface,
                                     spectral_rad_jac_surface,
                                     freq_grid,
                                     jac_targets,
                                     ray_point,
                                     surf_field,
                                     subsurf_field,
                                     spectral_rad_closed_surface_agenda);

#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size j = 0; j < NX; j++) {
    for (Size i = 0; i < NF; i++) {
      spectral_rad_jac[j, i] =
          rtepack::reflection(spectral_rad_jac[j, i],
                              spectral_surf_refl[i],
                              spectral_rad_jac_surface[j, i]);
    }
  }

  for (auto& target : jac_targets.surf) {
    const SurfaceData& data = surf_field[target.type];
    const auto ws           = data.flat_weights(lat, lon);

#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size i = 0; i < NF; i++) {
      for (const auto& [j, w] : ws) {
        spectral_rad_jac[j + target.x_start, i] +=
            w *
            rtepack::dreflection(spectral_rad[i],
                                 spectral_surf_refl_jac[target.target_pos, i],
                                 spectral_rad_surface[i]);
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < NF; i++) {
    spectral_rad[i] = rtepack::reflection(
        spectral_rad[i], spectral_surf_refl[i], spectral_rad_surface[i]);
  }
}
ARTS_METHOD_ERROR_CATCH

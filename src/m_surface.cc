#include <arts_omp.h>
#include <geodetic.h>
#include <sun.h>
#include <sun_methods.h>
#include <workspace.h>
#include "rtepack.h"

#include <algorithm>
#include <cmath>

namespace {
Vector2 specular_losNormal(const Vector2& normal, const Vector2& los, const Vector3& pos, const Vector2& ell) {
  const auto [ecef_pos, ecef_los] = geodetic_los2ecef(pos, los, ell);
  const auto [_, ecef_normal]     = geodetic_los2ecef(pos, normal, ell);

  // Specular direction is 2(dn*di)dn-di, where dn is the normal vector
  const Numeric fac       = 2 * dot(ecef_normal, ecef_los);
  const Vector3 ecef_spec = {
      fac * ecef_normal[0] - ecef_los[0], fac * ecef_normal[1] - ecef_los[1], fac * ecef_normal[2] - ecef_los[2]};

  return ecef2geodetic_los(ecef_pos, normalized(ecef_spec), ell).second;
}
}  // namespace

void spectral_surf_reflFlatRealFresnel(MuelmatVector&              spectral_surf_refl,
                                       MuelmatMatrix&              spectral_surf_refl_jac,
                                       const AscendingGrid&        freq_grid,
                                       const SurfaceField&         surf_field,
                                       const PropagationPathPoint& ray_point,
                                       const JacobianTargets&      jac_targets) try {
  ARTS_TIME_REPORT

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag refraction_target{"scalar refractive index"};

  ARTS_USER_ERROR_IF(not surf_field.contains(refraction_target),
                     R"--(Missing key property tag for method.

Tag "scalar refractive index" not in the surface field.

surf_field:
{}
)--",
                     surf_field);

  const auto&        pos        = ray_point.pos;
  const auto&        los        = ray_point.los;
  const auto&        ell        = surf_field.ellipsoid;
  const Numeric      lat        = pos[1];
  const Numeric      lon        = pos[2];
  const SurfacePoint surf_point = surf_field.at(lat, lon);
  const Numeric      n1         = ray_point.nreal;
  const Numeric      n2         = surf_point[refraction_target];
  const Size         nf         = freq_grid.size();

  spectral_surf_refl.resize(nf);
  spectral_surf_refl_jac.resize(jac_targets.target_count(), nf);

  const auto [Rv, Rh] = fresnel(
      n1,
      n2,
      std::acos(dot(geodetic_los2ecef(pos, los, ell).second, geodetic_los2ecef(pos, surf_point.normal, ell).second)));

  const Muelmat R        = rtepack::fresnel_reflectance(Rv, Rh);
  spectral_surf_refl     = R;
  spectral_surf_refl_jac = Muelmat{0.0};

  for (auto& target : jac_targets.surf) {
    if (target.type == refraction_target) {
      const auto [Rv2, Rh2] = fresnel(n1,
                                      n2 + 1e-3,
                                      std::acos(dot(geodetic_los2ecef(pos, los, ell).second,
                                                    geodetic_los2ecef(pos, surf_point.normal, ell).second)));

      const Muelmat dR = 1000. * (rtepack::fresnel_reflectance(Rv2, Rh2) - R);

      spectral_surf_refl_jac[target.target_pos] = dR;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_surf_reflFlatScalar(MuelmatVector&              spectral_surf_refl,
                                  MuelmatMatrix&              spectral_surf_refl_jac,
                                  const AscendingGrid&        freq_grid,
                                  const SurfaceField&         surf_field,
                                  const PropagationPathPoint& ray_point,
                                  const JacobianTargets&      jac_targets) try {
  ARTS_TIME_REPORT

  //! NOTE: Feel free to change the name and style of this key, it is unique to this method
  const SurfacePropertyTag reflectance_target{"flat scalar reflectance"};

  ARTS_USER_ERROR_IF(not surf_field.contains(reflectance_target),
                     R"--(Missing key property tag for method.

Tag "flat scalar reflectance" not in the surface field.

surf_field:
{}
)--",
                     surf_field);

  const Numeric      lat        = ray_point.pos[1];
  const Numeric      lon        = ray_point.pos[2];
  const SurfacePoint surf_point = surf_field.at(lat, lon);
  const Numeric      R          = surf_point[reflectance_target];
  const Size         nf         = freq_grid.size();

  ARTS_USER_ERROR_IF(R < 0.0 or R > 1.0, "Flat scalar reflectance must be between 0 and 1, but is {}.", R)

  spectral_surf_refl.resize(nf);
  spectral_surf_refl_jac.resize(jac_targets.target_count(), nf);

  // An isotropic mirror, Rv = -Rh = sqrt(R): a power reflectance has no phase between them
  spectral_surf_refl     = rtepack::power_reflectance(R);
  spectral_surf_refl_jac = Muelmat{0.0};

  for (auto& target : jac_targets.surf) {
    if (target.type == reflectance_target) {
      spectral_surf_refl_jac[target.target_pos] = rtepack::power_reflectance(1.0);
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSurfaceReflectance(const Workspace&            ws,
                                    StokvecVector&              spectral_rad,
                                    StokvecMatrix&              spectral_rad_jac,
                                    const AscendingGrid&        freq_grid,
                                    const AtmField&             atm_field,
                                    const SurfaceField&         surf_field,
                                    const SubsurfaceField&      subsurf_field,
                                    const JacobianTargets&      jac_targets,
                                    const PropagationPathPoint& ray_point,
                                    const Agenda&               spectral_rad_observer_agenda,
                                    const Agenda&               spectral_rad_closed_surface_agenda,
                                    const Agenda&               spectral_surf_refl_agenda) try {
  ARTS_TIME_REPORT

  const Size         NF         = freq_grid.size();
  const Size         NX         = jac_targets.x_size();
  const Numeric      lat        = ray_point.pos[1];
  const Numeric      lon        = ray_point.pos[2];
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
  const Vector2 los = specular_losNormal(surf_point.normal, ray_point.los, ray_point.pos, surf_field.ellipsoid);

  ArrayOfPropagationPathPoint ray_path;
  StokvecVector               spectral_rad_surface;
  StokvecMatrix               spectral_rad_jac_surface;

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
          rtepack::reflection(spectral_rad_jac[j, i], spectral_surf_refl[i], spectral_rad_jac_surface[j, i]);
    }
  }

  for (auto& target : jac_targets.surf) {
    const SurfaceData& data = surf_field[target.type];
    const auto         ws   = data.flat_weights(lat, lon);

#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size i = 0; i < NF; i++) {
      for (const auto& [j, w] : ws) {
        spectral_rad_jac[j + target.x_start, i] +=
            w * rtepack::dreflection(
                    spectral_rad[i], spectral_surf_refl_jac[target.target_pos, i], spectral_rad_surface[i]);
      }
    }
  }

#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i = 0; i < NF; i++) {
    spectral_rad[i] = rtepack::reflection(spectral_rad[i], spectral_surf_refl[i], spectral_rad_surface[i]);
  }
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSurfaceScatteringDiffuse(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const MapOfSurfaceScatteringModel& surface_models,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_point,
    const ZenGrid& zen_grid,
    const AziGrid& az_grid,
    const Vector&  zen_grid_weights,
    const Vector&   az_grid_weights,
    const Agenda& spectral_rad_incoming_agenda,
    const Agenda& spectral_rad_closed_surface_agenda) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(surf_field.bad_ellipsoid(),
                     "Surface field not properly set up - bad reference ellipsoid: {:B,}",
                     surf_field.ellipsoid)

  const Size nf = freq_grid.size();
  const Size nq = jac_targets.x_size();

  spectral_rad.resize(nf);
  spectral_rad = 0.0;

  spectral_rad_jac.resize(nq, nf);
  spectral_rad_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  const Vector       za_out     = {ray_point.los[0]};
  const Vector       aa_out     = {ray_point.los[1]};
  const SurfacePoint surf_point = surf_field.at(ray_point.pos[1], ray_point.pos[2]);

  // get the emissivity matrix and BRDF matrix, which are members of surface_props
  const surface_scattering::SurfaceScatteringModelProperties surface_props =
      surface_models.get_surface_scattering_model_properties(
          surf_point,
                                                             ray_point.pos[1],
                                                             ray_point.pos[2],
                                                             freq_grid,
                                                             zen_grid,
                                                             az_grid,
                                                             za_out,
                                                             aa_out);

  // get the subsurface emission
  StokvecVector spectral_rad_surface;
  StokvecMatrix spectral_rad_jac_surface;
  spectral_rad_surface_agendaExecute(ws,
                                     spectral_rad_surface,
                                     spectral_rad_jac_surface,
                                     freq_grid,
                                     jac_targets,
                                     ray_point,
                                     surf_field,
                                     subsurf_field,
                                     spectral_rad_closed_surface_agenda);

  // The surface normal is stored with the outward direction at za = 180, so its
  // ECEF image points inward; an incoming direction is visible when it opposes
  // that inward vector, i.e. when the dot product with it is negative.  This is
  // the same horizon visibility test as in spectral_radSurfaceScatteringDirect,
  // applied per quadrature direction: directions below the horizon defined by
  // the actual surface normal contribute hard-zero (no error), emission remains.
  const auto [__, ecef_normal] = geodetic_los2ecef(ray_point.pos, surf_point.normal, surf_field.ellipsoid);

  // get the incoming radiation
  StokvecTensor3 spectral_rad_incoming(zen_grid.size(), az_grid.size(), freq_grid.size());
  StokvecTensor4 spectral_rad_incoming_jac(zen_grid.size(), az_grid.size(), jac_targets.x_size(), freq_grid.size());
  for (Size j = 0; j < zen_grid.size(); j++) {
    for (Size k = 0; k < az_grid.size(); k++) {
      const Vector2 los_incoming = {zen_grid[j], az_grid[k]};

      const auto [_, ecef_los_in] = geodetic_los2ecef(ray_point.pos, los_incoming, surf_field.ellipsoid);
      if (dot(ecef_los_in, ecef_normal) >= 0.0) continue;  // Below the tilted horizon: zero contribution

      StokvecVector spectral_rad_incoming_temp;
      StokvecMatrix spectral_rad_incoming_jac_temp;

      spectral_rad_incoming_agendaExecute(ws,
                                          spectral_rad_incoming_temp,
                                          spectral_rad_incoming_jac_temp,
                                          freq_grid,
                                          jac_targets,
                                          ray_point.pos,
                                          los_incoming,
                                          atm_field,
                                          surf_field,
                                          subsurf_field,
                                          spectral_rad_incoming_agenda);

      spectral_rad_incoming[j, k, joker] = spectral_rad_incoming_temp;
      spectral_rad_incoming_jac[j,k,joker,joker] =
          spectral_rad_incoming_jac_temp;
    }
  }

  // Calculate scattered radiation
  StokvecVector  spectral_rad_scattered(freq_grid.size());
  StokvecMatrix  spectral_rad_scattered_jac(jac_targets.x_size(),freq_grid.size());

  // integrate over the incoming directions to get the scattered upward radiation
  for (Size i_za  = 0; i_za < zen_grid.size(); i_za ++) {
    for (Size i_aa = 0; i_aa < az_grid.size(); i_aa ++) {

      // calculate scattered upward radiation
      for (Size i_f = 0; i_f < freq_grid.size(); i_f++) {
        const Muelmat R  = surface_props.brdf_matrix_diffuse[i_f, i_za, i_aa, 0, 0];
        spectral_rad_scattered[i_f] += R * spectral_rad_incoming[i_za, i_aa, i_f] *
                           zen_grid_weights[i_za] * az_grid_weights[i_aa];

      }

      //Calculate scattered upward radiation jacobian 
      //For now, there is no jacobian for the surface scattering model!!!
      for (Size i_jac = 0; i_jac < jac_targets.x_size(); i_jac++) {
        for (Size i_f = 0; i_f < freq_grid.size(); i_f++) {
          const Muelmat R  = surface_props.brdf_matrix_diffuse[i_f, i_za, i_aa, 0, 0];
          spectral_rad_scattered_jac[i_jac, i_f] += R * spectral_rad_incoming_jac[i_za, i_aa, i_jac, i_f] *
                           zen_grid_weights[i_za] * az_grid_weights[i_aa];
        }
      }
    }
  }

  // Calculate upward emission
  for (Size i_f = 0; i_f < freq_grid.size(); i_f++) {
    spectral_rad_surface[i_f] = surface_props.emissivity_vector_diffuse[i_f, 0, 0] * spectral_rad_surface[i_f];
  }

  //Calculate jacobian for subsurface emission
  //For now, there is no jacobian for the surface scattering model!!!
  StokvecMatrix spectral_rad_jac_subsurface(jac_targets.x_size(), freq_grid.size());
  for (Size i_jac = 0; i_jac < jac_targets.x_size(); i_jac++) {
    for (Size i_f = 0; i_f < freq_grid.size(); i_f++) {
      spectral_rad_jac_subsurface[i_jac, i_f] = surface_props.emissivity_vector_diffuse[i_f, 0, 0] * spectral_rad_jac_surface[i_jac, i_f];
    }
  }



  spectral_rad += spectral_rad_scattered;
  spectral_rad += spectral_rad_surface;
  spectral_rad_jac += spectral_rad_scattered_jac;
  spectral_rad_jac += spectral_rad_jac_subsurface;

}
ARTS_METHOD_ERROR_CATCH

void spectral_radSurfaceScatteringSpecular(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const MapOfSurfaceScatteringModel& surface_models,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_point,
    const Agenda& spectral_rad_incoming_agenda,
    const Agenda& spectral_rad_closed_surface_agenda) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(surf_field.bad_ellipsoid(),
                     "Surface field not properly set up - bad reference ellipsoid: {:B,}",
                     surf_field.ellipsoid)

  const Size nf = freq_grid.size();
  const Size nq = jac_targets.x_size();

  spectral_rad.resize(nf);
  spectral_rad = 0.0;

  spectral_rad_jac.resize(nq, nf);
  spectral_rad_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  const SurfacePoint surf_point = surf_field.at(ray_point.pos[1], ray_point.pos[2]);

  // Incoming direction is the mirror reflection of the outgoing direction
  // about the local surface normal (see spectral_radSurfaceReflectance)
  const Vector2 los_in = specular_losNormal(surf_point.normal, ray_point.los, ray_point.pos, surf_field.ellipsoid);

  // get the emissivity vector and BRDF matrix at the exact (single)
  // incident and outgoing directions
  const Vector za_in  {los_in[0]};
  const Vector aa_in  {los_in[1]};
  const Vector za_out {ray_point.los[0]};
  const Vector aa_out {ray_point.los[1]};

  const surface_scattering::SurfaceScatteringModelProperties surface_props =
      surface_models.get_surface_scattering_model_properties(
          surf_point,
          ray_point.pos[1],
          ray_point.pos[2],
          freq_grid,
          za_in,
          aa_in,
          za_out,
          aa_out);

  // get the subsurface emission
  StokvecVector spectral_rad_surface;
  StokvecMatrix spectral_rad_jac_surface;
  spectral_rad_surface_agendaExecute(ws,
                                     spectral_rad_surface,
                                     spectral_rad_jac_surface,
                                     freq_grid,
                                     jac_targets,
                                     ray_point,
                                     surf_field,
                                     subsurf_field,
                                     spectral_rad_closed_surface_agenda);

  // get the incoming radiation from the single specular direction
  StokvecVector spectral_rad_incoming;
  StokvecMatrix spectral_rad_incoming_jac;
  spectral_rad_incoming_agendaExecute(ws,
                                      spectral_rad_incoming,
                                      spectral_rad_incoming_jac,
                                      freq_grid,
                                      jac_targets,
                                      ray_point.pos,
                                      los_in,
                                      atm_field,
                                      surf_field,
                                      subsurf_field,
                                      spectral_rad_incoming_agenda);

  // Calculate reflected radiation (no angular integration: single direction)
  StokvecVector spectral_rad_reflected(nf);
  spectral_rad_reflected = 0.0;
  StokvecMatrix spectral_rad_reflected_jac(nq, nf);
  spectral_rad_reflected_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  // For now, there is no jacobian for the surface scattering model!!!
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i_f = 0; i_f < nf; i_f++) {
    spectral_rad_reflected[i_f] +=
        surface_props.brdf_matrix_specular[i_f, 0, 0, 0, 0] * spectral_rad_incoming[i_f];
  }
#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size i_jac = 0; i_jac < nq; i_jac++) {
    for (Size i_f = 0; i_f < nf; i_f++) {
      spectral_rad_reflected_jac[i_jac, i_f] =
          surface_props.brdf_matrix_specular[i_f, 0, 0, 0, 0] * spectral_rad_incoming_jac[i_jac, i_f];
    }
  }

  // Calculate upward emission with the specular emissivity
  // For now, there is no jacobian for the surface scattering model!!!
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i_f = 0; i_f < nf; i_f++) {
    spectral_rad_surface[i_f] = surface_props.emissivity_vector_specular[i_f, 0, 0] * spectral_rad_surface[i_f];
  }
#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size i_jac = 0; i_jac < nq; i_jac++) {
    for (Size i_f = 0; i_f < nf; i_f++) {
      spectral_rad_jac[i_jac, i_f] =
          surface_props.emissivity_vector_specular[i_f, 0, 0] * spectral_rad_jac_surface[i_jac, i_f];
    }
  }

  spectral_rad += spectral_rad_reflected;
  spectral_rad += spectral_rad_surface;
  spectral_rad_jac += spectral_rad_reflected_jac;
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSurfaceScatteringSpecularDirect(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const MapOfSurfaceScatteringModel& surface_models,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_point,
    const ArrayOfSun& suns,
    const Agenda& ray_path_observer_agenda,
    const Agenda& spectral_rad_incoming_agenda,
    const Agenda& spectral_rad_closed_surface_agenda,
    const Numeric& angle_cut,
    const Index& refinement) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(surf_field.bad_ellipsoid(),
                     "Surface field not properly set up - bad reference ellipsoid: {:B,}",
                     surf_field.ellipsoid)

  const Size nf = freq_grid.size();
  const Size nq = jac_targets.x_size();

  spectral_rad.resize(nf);
  spectral_rad = 0.0;

  spectral_rad_jac.resize(nq, nf);
  spectral_rad_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  const SurfacePoint surf_point = surf_field.at(ray_point.pos[1], ray_point.pos[2]);

  // Glint direction: the mirror reflection of the outgoing (ray) direction about
  // the local surface normal.  It is stored as the looking direction towards the
  // beam source that mirrors into the ray direction, i.e. directly comparable
  // with sun_geometric_los and the obs_los of spectral_rad_incoming_agenda.
  const Vector2 los_spec =
      specular_losNormal(surf_point.normal, ray_point.los, ray_point.pos, surf_field.ellipsoid);

  // Single outgoing (ray) direction; the incoming directions are the sun beams
  const Vector za_out {ray_point.los[0]};
  const Vector aa_out {ray_point.los[1]};

  // The specular emissivity depends only on the outgoing direction, so it is
  // evaluated once here with a dummy incidence for the emission term
  const surface_scattering::SurfaceScatteringModelProperties emission_props =
      surface_models.get_surface_scattering_model_properties(
          surf_point,
          ray_point.pos[1],
          ray_point.pos[2],
          freq_grid,
          Vector{0.0},
          Vector{0.0},
          za_out,
          aa_out);

  // get the subsurface emission
  StokvecVector spectral_rad_surface;
  StokvecMatrix spectral_rad_jac_surface;
  spectral_rad_surface_agendaExecute(ws,
                                     spectral_rad_surface,
                                     spectral_rad_jac_surface,
                                     freq_grid,
                                     jac_targets,
                                     ray_point,
                                     surf_field,
                                     subsurf_field,
                                     spectral_rad_closed_surface_agenda);

  // ECEF image of the glint direction, used for the refractive disc re-test
  const auto [_, ecef_los_spec] = geodetic_los2ecef(ray_point.pos, los_spec, surf_field.ellipsoid);

  StokvecVector spectral_rad_scattered(nf);
  spectral_rad_scattered = 0.0;
  StokvecMatrix spectral_rad_scattered_jac(nq, nf);
  spectral_rad_scattered_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  // One delta-function beam per sun, contributions summed; no suns means no
  // scattered term at all
  for (const auto& sun : suns) {
    // Solar disc test on the glint direction: the sun contributes only if the
    // specular direction falls inside the disc, i.e. beta <= alpha with alpha
    // the angular radius of the sun at the point.  This gate subsumes the
    // horizon gate: a sun on the far side of the (possibly tilted) surface
    // normal can never satisfy beta <= alpha.
    const auto [beta, hit] = hit_sun(sun, ray_point.pos, los_spec, surf_field.ellipsoid);
    if (not hit) continue;

    // Refraction-aware beam direction via the observer agenda
    const Vector2 los_in =
        sun_refractive_los(ws, sun, ray_point.pos, surf_field, ray_path_observer_agenda, angle_cut, refinement);

    // Disc re-test on the refracted LOS against the glint direction
    const Numeric alpha = std::asin(std::sqrt(sun.sin_alpha_squared(ray_point.pos, surf_field.ellipsoid)));
    const auto [__, ecef_los_in] = geodetic_los2ecef(ray_point.pos, los_in, surf_field.ellipsoid);
    const Numeric beta_ref = std::acos(std::clamp(dot(ecef_los_spec, ecef_los_in), -1.0, 1.0));
    if (beta_ref > alpha) continue;

    const Vector za_in {los_in[0]};
    const Vector aa_in {los_in[1]};

    const surface_scattering::SurfaceScatteringModelProperties sun_props =
        surface_models.get_surface_scattering_model_properties(
            surf_point,
            ray_point.pos[1],
            ray_point.pos[2],
            freq_grid,
            za_in,
            aa_in,
            za_out,
            aa_out);

    // get the incoming radiation from this single beam direction
    StokvecVector spectral_rad_incoming;
    StokvecMatrix spectral_rad_incoming_jac;
    spectral_rad_incoming_agendaExecute(ws,
                                        spectral_rad_incoming,
                                        spectral_rad_incoming_jac,
                                        freq_grid,
                                        jac_targets,
                                        ray_point.pos,
                                        los_in,
                                        atm_field,
                                        surf_field,
                                        subsurf_field,
                                        spectral_rad_incoming_agenda);

    // No quadrature weights: the beam radiance is delta-weighted; radiance is
    // conserved by specular reflection so BRDF * I_in is the exact normalisation
    // For now, there is no jacobian for the surface scattering model!!!
#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size i_f = 0; i_f < nf; i_f++) {
      spectral_rad_scattered[i_f] += sun_props.brdf_matrix_specular[i_f, 0, 0, 0, 0] * spectral_rad_incoming[i_f];
    }
#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
    for (Size i_jac = 0; i_jac < nq; i_jac++) {
      for (Size i_f = 0; i_f < nf; i_f++) {
        spectral_rad_scattered_jac[i_jac, i_f] +=
            sun_props.brdf_matrix_specular[i_f, 0, 0, 0, 0] * spectral_rad_incoming_jac[i_jac, i_f];
      }
    }
  }

  // Calculate upward emission with the specular emissivity
  // For now, there is no jacobian for the surface scattering model!!!
  StokvecMatrix spectral_rad_jac_subsurface(nq, nf);
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i_f = 0; i_f < nf; i_f++) {
    spectral_rad_surface[i_f] = emission_props.emissivity_vector_specular[i_f, 0, 0] * spectral_rad_surface[i_f];
  }
#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size i_jac = 0; i_jac < nq; i_jac++) {
    for (Size i_f = 0; i_f < nf; i_f++) {
      spectral_rad_jac_subsurface[i_jac, i_f] =
          emission_props.emissivity_vector_specular[i_f, 0, 0] * spectral_rad_jac_surface[i_jac, i_f];
    }
  }

  spectral_rad += spectral_rad_scattered;
  spectral_rad += spectral_rad_surface;
  spectral_rad_jac += spectral_rad_scattered_jac;
  spectral_rad_jac += spectral_rad_jac_subsurface;
}
ARTS_METHOD_ERROR_CATCH

void spectral_radSurfaceScatteringDirect(
    const Workspace& ws,
    StokvecVector& spectral_rad,
    StokvecMatrix& spectral_rad_jac,
    const AscendingGrid& freq_grid,
    const AtmField& atm_field,
    const SurfaceField& surf_field,
    const SubsurfaceField& subsurf_field,
    const MapOfSurfaceScatteringModel& surface_models,
    const JacobianTargets& jac_targets,
    const PropagationPathPoint& ray_point,
    const ArrayOfSun& suns,
    const Agenda& ray_path_observer_agenda,
    const Agenda& spectral_rad_incoming_agenda,
    const Agenda& spectral_rad_closed_surface_agenda,
    const Numeric& angle_cut,
    const Index& refinement) try {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(surf_field.bad_ellipsoid(),
                     "Surface field not properly set up - bad reference ellipsoid: {:B,}",
                     surf_field.ellipsoid)

  const Size nf = freq_grid.size();
  const Size nq = jac_targets.x_size();

  spectral_rad.resize(nf);
  spectral_rad = 0.0;

  spectral_rad_jac.resize(nq, nf);
  spectral_rad_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  const SurfacePoint surf_point = surf_field.at(ray_point.pos[1], ray_point.pos[2]);

  // Single outgoing (ray) direction; the incoming directions are the sun beams
  const Vector za_out {ray_point.los[0]};
  const Vector aa_out {ray_point.los[1]};

  // The diffuse emissivity depends only on the outgoing direction, so it is
  // evaluated once here with a dummy incidence for the emission term
  const surface_scattering::SurfaceScatteringModelProperties emission_props =
      surface_models.get_surface_scattering_model_properties(
          surf_point,
          ray_point.pos[1],
          ray_point.pos[2],
          freq_grid,
          Vector{0.0},
          Vector{0.0},
          za_out,
          aa_out);

  // get the subsurface emission
  StokvecVector spectral_rad_surface;
  StokvecMatrix spectral_rad_jac_surface;
  spectral_rad_surface_agendaExecute(ws,
                                     spectral_rad_surface,
                                     spectral_rad_jac_surface,
                                     freq_grid,
                                     jac_targets,
                                     ray_point,
                                     surf_field,
                                     subsurf_field,
                                     spectral_rad_closed_surface_agenda);

  // A beam is visible only above the horizon defined by the actual surface
  // normal; below it the direct term is hard-zero (no error), emission remains.
  // The surface normal is stored with the outward direction at za = 180, so its
  // ECEF image points inward; a beam is visible when it opposes that inward
  // vector, i.e. when the dot product with it is negative.
  const auto [__, ecef_normal] = geodetic_los2ecef(ray_point.pos, surf_point.normal, surf_field.ellipsoid);
  const auto above_horizon     = [&](const Vector2& los) {
    const auto [_, ecef_los_in] = geodetic_los2ecef(ray_point.pos, los, surf_field.ellipsoid);
    return dot(ecef_los_in, ecef_normal) < 0.0;
  };

  StokvecVector spectral_rad_scattered(nf);
  spectral_rad_scattered = 0.0;
  StokvecMatrix spectral_rad_scattered_jac(nq, nf);
  spectral_rad_scattered_jac = Stokvec{0.0, 0.0, 0.0, 0.0};

  // One delta-function beam per sun, contributions summed; no suns means no
  // scattered term at all
  for (const auto& sun : suns) {
    // Gate on the geometric LOS first to avoid the search for invisible suns
    Vector2 los = sun_geometric_los(sun, ray_point.pos, surf_field);
    if (not above_horizon(los)) continue;

    // Refraction-aware beam direction via the observer agenda
    los = sun_refractive_los(ws, sun, ray_point.pos, surf_field, ray_path_observer_agenda, angle_cut, refinement);
    if (not above_horizon(los)) continue;

    const Vector za_in {los[0]};
    const Vector aa_in {los[1]};

    const surface_scattering::SurfaceScatteringModelProperties sun_props =
        surface_models.get_surface_scattering_model_properties(
            surf_point,
            ray_point.pos[1],
            ray_point.pos[2],
            freq_grid,
            za_in,
            aa_in,
            za_out,
            aa_out);

    // get the incoming radiation from this single beam direction
    StokvecVector spectral_rad_incoming;
    StokvecMatrix spectral_rad_incoming_jac;
    spectral_rad_incoming_agendaExecute(ws,
                                        spectral_rad_incoming,
                                        spectral_rad_incoming_jac,
                                        freq_grid,
                                        jac_targets,
                                        ray_point.pos,
                                        los,
                                        atm_field,
                                        surf_field,
                                        subsurf_field,
                                        spectral_rad_incoming_agenda);

    // No quadrature weights: the beam radiance is delta-weighted
    // For now, there is no jacobian for the surface scattering model!!!
#pragma omp parallel for if (not arts_omp_in_parallel())
    for (Size i_f = 0; i_f < nf; i_f++) {
      spectral_rad_scattered[i_f] += sun_props.brdf_matrix_diffuse[i_f, 0, 0, 0, 0] * spectral_rad_incoming[i_f];
    }
#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
    for (Size i_jac = 0; i_jac < nq; i_jac++) {
      for (Size i_f = 0; i_f < nf; i_f++) {
        spectral_rad_scattered_jac[i_jac, i_f] +=
            sun_props.brdf_matrix_diffuse[i_f, 0, 0, 0, 0] * spectral_rad_incoming_jac[i_jac, i_f];
      }
    }
  }

  // Calculate upward emission with the diffuse emissivity
  // For now, there is no jacobian for the surface scattering model!!!
  StokvecMatrix spectral_rad_jac_subsurface(nq, nf);
#pragma omp parallel for if (not arts_omp_in_parallel())
  for (Size i_f = 0; i_f < nf; i_f++) {
    spectral_rad_surface[i_f] = emission_props.emissivity_vector_diffuse[i_f, 0, 0] * spectral_rad_surface[i_f];
  }
#pragma omp parallel for collapse(2) if (not arts_omp_in_parallel())
  for (Size i_jac = 0; i_jac < nq; i_jac++) {
    for (Size i_f = 0; i_f < nf; i_f++) {
      spectral_rad_jac_subsurface[i_jac, i_f] =
          emission_props.emissivity_vector_diffuse[i_f, 0, 0] * spectral_rad_jac_surface[i_jac, i_f];
    }
  }

  spectral_rad += spectral_rad_scattered;
  spectral_rad += spectral_rad_surface;
  spectral_rad_jac += spectral_rad_scattered_jac;
  spectral_rad_jac += spectral_rad_jac_subsurface;
}
ARTS_METHOD_ERROR_CATCH


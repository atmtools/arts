#include <workspace.h>

#include <ranges>
#include "path_point.h"

void subsurface_profileFromPath(ArrayOfSubsurfacePoint subsurface_profile,
                                const SubsurfaceField& subsurface_field,
                                const ArrayOfPropagationPathPoint& ray_path) {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      not subsurface_field.contains(SubsurfaceKey::t) or
          not subsurface_field.contains(SubsurfacePropertyTag{"scalar tau"}) or
          not subsurface_field.contains(SubsurfacePropertyTag{"scalar omega"}),
      R"(Sub-surface field does not contain all required components.

It must contain: `temperature`, `scalar tau`, `scalar omega`
)")

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

  ray_path = depth_profile | stdv::transform([ray_path_point](const auto& d) {
               auto v   = ray_path_point;
               v.pos[0] = d;
               return v;
             }) |
             stdr::to<ArrayOfPropagationPathPoint>();
}

void spectral_radianceScalarSubsurface(
    StokvecVector& spectral_radiance,
    StokvecMatrix& spectral_radiance_jacobian,
    const AscendingGrid& frequency_grid,
    const JacobianTargets& jacobian_targets,
    const PropagationPathPoint& ray_path_point,
    const SurfaceField& surface_field,
    const SubsurfaceField& subsurface_field,
    const DescendingGrid& depth_profile) {
  ARTS_TIME_REPORT

  ArrayOfPropagationPathPoint ray_path;
  ray_pathFromPointAndDepth(ray_path, ray_path_point, depth_profile);
}

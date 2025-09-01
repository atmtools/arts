#include <workspace.h>

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

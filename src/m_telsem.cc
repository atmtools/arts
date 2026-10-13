#include <arts_conversions.h>
#include <geodetic.h>
#include <workspace.h>

namespace {
Numeric incidence_angle(const SurfacePoint&         surf_point,
                        const PropagationPathPoint& ray_point,
                        const Vector2&              ellipsoid) {
  const auto los    = geodetic_los2ecef(ray_point.pos, ray_point.los, ellipsoid).second;
  const auto normal = geodetic_los2ecef(ray_point.pos, surf_point.normal, ellipsoid).second;
  return Conversion::rad2deg(std::acos(std::clamp(std::abs(dot(los, normal)), Numeric{0}, Numeric{1})));
}
}  // namespace

void telsem_atlasReadAscii(TelsemAtlas& telsem_atlas, const String& filename, const Index& month) {
  ARTS_TIME_REPORT
  telsem_read_ascii(filename, telsem_atlas, month);
}

void spectral_surf_reflTelsem(MuelmatVector&              spectral_surf_refl,
                              MuelmatMatrix&              spectral_surf_refl_jac,
                              const AscendingGrid&        freq_grid,
                              const SurfaceField&         surf_field,
                              const PropagationPathPoint& ray_point,
                              const JacobianTargets&      jac_targets,
                              const TelsemAtlas&          telsem_atlas,
                              const Numeric&              max_distance) try {
  ARTS_TIME_REPORT

  const Numeric      lat   = ray_point.pos[1];
  const Numeric      lon   = ray_point.pos[2];
  const SurfacePoint point = surf_field.at(lat, lon);
  const Numeric      angle = incidence_angle(point, ray_point, surf_field.ellipsoid);

  spectral_surf_refl.resize(freq_grid.size());
  spectral_surf_refl_jac.resize(jac_targets.target_count(), freq_grid.size());
  spectral_surf_refl_jac = Muelmat{0.0};

  for (Size iv = 0; iv < freq_grid.size(); ++iv) {
    const auto emissivity  = telsem_atlas.emissivity(lat, lon, angle, freq_grid[iv], max_distance);
    spectral_surf_refl[iv] = rtepack::power_reflectance(std::clamp(1 - emissivity[0], Numeric{0}, Numeric{1}),
                                                        std::clamp(1 - emissivity[1], Numeric{0}, Numeric{1}));
  }
}
ARTS_METHOD_ERROR_CATCH

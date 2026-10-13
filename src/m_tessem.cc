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

Muelmat reflectance(const TessemNN& neth,
                    const TessemNN& netv,
                    Numeric         frequency,
                    Numeric         angle,
                    Numeric         wind_speed,
                    Numeric         temperature,
                    Numeric         salinity) {
  const Vector input{frequency, angle, wind_speed, temperature, salinity};
  const Vector eh = tessem_emissivity(neth, input);
  const Vector ev = tessem_emissivity(netv, input);
  ARTS_USER_ERROR_IF(eh.empty() or ev.empty(), "TESSEM neural networks must have at least one output")
  return rtepack::power_reflectance(std::clamp(1 - ev[0], Numeric{0}, Numeric{1}),
                                    std::clamp(1 - eh[0], Numeric{0}, Numeric{1}));
}
}  // namespace

void tessem_nnReadAscii(TessemNN& tessem_nn, const String& filename) {
  ARTS_TIME_REPORT
  tessem_read_ascii(filename, tessem_nn);
}

void spectral_surf_reflTessem(MuelmatVector&              spectral_surf_refl,
                              MuelmatMatrix&              spectral_surf_refl_jac,
                              const AscendingGrid&        freq_grid,
                              const SurfaceField&         surf_field,
                              const PropagationPathPoint& ray_point,
                              const JacobianTargets&      jac_targets,
                              const TessemNN&             tessem_neth,
                              const TessemNN&             tessem_netv) try {
  ARTS_TIME_REPORT

  constexpr auto           temperature_key = SurfaceKey::t;
  const SurfacePropertyTag wind_key{"wind speed"};
  const SurfacePropertyTag salinity_key{"salinity"};

  ARTS_USER_ERROR_IF(not surf_field.contains(temperature_key), "TESSEM requires surface temperature")
  ARTS_USER_ERROR_IF(not surf_field.contains(wind_key), "TESSEM requires surf_field property \"wind speed\" [m/s]")
  ARTS_USER_ERROR_IF(not surf_field.contains(salinity_key), "TESSEM requires surf_field property \"salinity\" [kg/kg]")

  const Numeric      lat         = ray_point.pos[1];
  const Numeric      lon         = ray_point.pos[2];
  const SurfacePoint point       = surf_field.at(lat, lon);
  const Numeric      angle       = incidence_angle(point, ray_point, surf_field.ellipsoid);
  const Numeric      wind        = point[wind_key];
  const Numeric      temperature = point[temperature_key];
  const Numeric      salinity    = point[salinity_key];

  spectral_surf_refl.resize(freq_grid.size());
  spectral_surf_refl_jac.resize(jac_targets.target_count(), freq_grid.size());
  spectral_surf_refl_jac = Muelmat{0.0};

  for (Size iv = 0; iv < freq_grid.size(); ++iv) {
    spectral_surf_refl[iv] = reflectance(tessem_neth, tessem_netv, freq_grid[iv], angle, wind, temperature, salinity);
  }

  for (const auto& target : jac_targets.surf) {
    const bool is_temperature = target.type == temperature_key;
    const bool is_wind        = target.type == wind_key;
    const bool is_salinity    = target.type == salinity_key;
    if (not(is_temperature or is_wind or is_salinity)) continue;
    ARTS_USER_ERROR_IF(target.d == 0, "TESSEM Jacobian perturbation must be nonzero for target {}", target.type)

    for (Size iv = 0; iv < freq_grid.size(); ++iv) {
      const Muelmat perturbed                       = reflectance(tessem_neth,
                                                                  tessem_netv,
                                                                  freq_grid[iv],
                                                                  angle,
                                                                  wind + (is_wind ? target.d : 0),
                                                                  temperature + (is_temperature ? target.d : 0),
                                                                  salinity + (is_salinity ? target.d : 0));
      spectral_surf_refl_jac[target.target_pos, iv] = (perturbed - spectral_surf_refl[iv]) / target.d;
    }
  }
}
ARTS_METHOD_ERROR_CATCH

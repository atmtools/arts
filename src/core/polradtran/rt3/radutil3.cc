#include "radutil3.h"

#include <debug.h>
#include <radutil.h>

#include <array>
#include <type_traits>
#include <variant>

namespace polradtran::rt3 {
void ground_surface(const surface&  ground,
                    Index           src_code,
                    ConstVectorView mu_values,
                    ConstVectorView quad_weights,
                    Numeric         frequency,
                    Numeric         ground_temp,
                    Tensor5View     surf_reflect,
                    Tensor3View     gnd_radiance,
                    Tensor3View     direct_reflect) {
  const Index nmode   = gnd_radiance.extent(0);
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = gnd_radiance.extent(2);
  ARTS_USER_ERROR_IF(nmode < 1 or nstokes < 1 or nstokes > 4 or quad_weights.extent(0) != nummu or
                         surf_reflect.shape() != (std::array<Index, 5>{nmode, nummu, nstokes, nummu, nstokes}) or
                         gnd_radiance.shape() != (std::array<Index, 3>{nmode, nummu, nstokes}) or
                         direct_reflect.shape() != gnd_radiance.shape(),
                     "The ground with {} mu_values needs quad_weights [nummu], surf_reflect [aziorder + 1, nummu, "
                     "nstokes, nummu, nstokes] and gnd_radiance and direct_reflect [aziorder + 1, nummu, nstokes] "
                     "with 1 to 4 Stokes components; got {}, {:B,}, {:B,} and {:B,}",
                     nummu,
                     quad_weights.extent(0),
                     surf_reflect.shape(),
                     gnd_radiance.shape(),
                     direct_reflect.shape());
  ARTS_USER_ERROR_IF(src_code < 0 or src_code > 3, "The ground needs SRC_CODE 0 to 3, got {}", src_code);
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "The ground needs a positive frequency, got {} Hz", frequency);

  const bool solar   = src_code == 1 or src_code == 3;
  const bool thermal = src_code == 2 or src_code == 3;

  // The ground routines make the whole surface layer, of which
  // SURF_REFLECT is REFLECT(..., 2)
  Tensor5 reflect(2, nummu, nstokes, nummu, nstokes), trans(2, nummu, nstokes, nummu, nstokes);
  Tensor3 source(2, nummu, nstokes);

  std::visit(
      [&](const auto& g) {
        using T = std::remove_cvref_t<decltype(g)>;
        if constexpr (std::is_same_v<T, lambertian_surface>) {
          for (Index mode = 0; mode < nmode; mode++) {
            // For a Lambertian surface
            lambert_surface_layer(mode, mu_values, quad_weights, g.albedo, reflect, trans, source);
            surf_reflect[mode] = reflect[1];
            // The radiance from the ground is thermal (without the solar
            // source and direct flux) ...
            lambert_radiance(mode, thermal ? 2 : 0, g.albedo, ground_temp, frequency, 0.0, gnd_radiance[mode]);
            // ... and reflected direct, for a unit direct flux (the
            // solar source alone)
            lambert_radiance(mode, 1, g.albedo, ground_temp, frequency, 1.0, direct_reflect[mode]);
          }
        } else {
          static_assert(std::is_same_v<T, fresnel_surface>);
          ARTS_USER_ERROR_IF(solar,
                             "RT3 cannot reflect the direct beam from a Fresnel ground (its reflection is specular); "
                             "the solar source needs a Lambertian ground");
          // For a Fresnel surface, the same in every mode
          fresnel_surface_layer(mu_values, g.refractive_index, reflect, trans, source);
          for (Index mode = 0; mode < nmode; mode++) {
            surf_reflect[mode] = reflect[1];
            // The radiance from the ground is thermal
            fresnel_radiance(mode, mu_values, g.refractive_index, ground_temp, frequency, gnd_radiance[mode]);
          }
          direct_reflect = 0.0;
        }
      },
      ground);
}
}  // namespace polradtran::rt3

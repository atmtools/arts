#include "radutil4.h"

#include <arts_conversions.h>
#include <debug.h>
#include <physics_funcs.h>
#include <rtepack.h>

#include <array>
#include <type_traits>
#include <variant>

namespace rt4 {
void lambert_surface_layer(Index           mode,
                           ConstVectorView mu_values,
                           ConstVectorView quad_weights,
                           Numeric         ground_albedo,
                           Tensor5View     reflect,
                           Tensor5View     trans,
                           Tensor3View     source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(quad_weights.extent(0) != nummu or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "LAMBERT_SURFACE with {} mu_values needs quad_weights [nummu], reflect and trans [2, nummu, "
                     "nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {}, {:B,}, {:B,} and {:B,}",
                     nummu,
                     quad_weights.extent(0),
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));
  // The Lambertian ground reflects the flux equally in all direction
  // and completely unpolarizes the radiation
  if (mode == 0) {
    for (Index j2 = 0; j2 < nummu; j2++)
      reflect[1, j2, 0, joker, 0] = 2.0 * ground_albedo * mu_values[j2] * quad_weights[j2];
  }
}

void lambert_radiance(Numeric ground_albedo, Numeric ground_temp, Numeric frequency, MatrixView radiance) {
  ARTS_USER_ERROR_IF(
      not(ground_temp >= 0.0), "LAMBERT_RADIANCE needs a ground temperature >= 0 K, got {} K", ground_temp);

  // Thermal radiation going up
  radiance = 0.0;

  const Numeric thermal = (1.0 - ground_albedo) * planck(frequency, ground_temp);
  radiance[joker, 0]    = thermal;
}

void fresnel_surface_layer(
    ConstVectorView mu_values, Complex index, Tensor5View reflect, Tensor5View trans, Tensor3View source) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = source.extent(2);
  ARTS_USER_ERROR_IF(nstokes > 4 or reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "FRESNEL_SURFACE with {} mu_values needs at most 4 Stokes components and reflect and trans "
                     "[2, nummu, nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {:B,}, {:B,} and {:B,}",
                     nummu,
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));

  // REFLECT(I, J, K, J, 2) = R(I, K): R1 and R2 in the [I, Q] block, R3 and
  // R4 in the [U, V] block
  const Range stokes{0, nstokes};
  for (Index j = 0; j < nummu; j++) {
    const auto [rv, rh]            = fresnel(1.0, index, Conversion::acosd(mu_values[j]));
    const Muelmat r                = rtepack::fresnel_reflectance(rv, rh);
    reflect[1, j, joker, j, joker] = transpose(r.view()[stokes, stokes]);
  }
}

void fresnel_radiance(
    ConstVectorView mu_values, Complex index, Numeric ground_temp, Numeric frequency, MatrixView radiance) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = radiance.extent(1);
  ARTS_USER_ERROR_IF(nstokes > 4 or radiance.extent(0) != nummu,
                     "FRESNEL_RADIANCE with {} mu_values needs radiance [nummu, nstokes <= 4]; got {:B,}",
                     nummu,
                     radiance.shape());
  ARTS_USER_ERROR_IF(
      not(ground_temp >= 0.0), "FRESNEL_RADIANCE needs a ground temperature >= 0 K, got {} K", ground_temp);

  // Thermal radiation going up: [(1 - R1) B, -R2 B, 0, 0]
  radiance = 0.0;
  const Range   stokes{0, nstokes};
  const Stokvec planck_ground{planck(frequency, ground_temp)};
  for (Index j = 0; j < nummu; j++) {
    const auto [rv, rh] = fresnel(1.0, index, Conversion::acosd(mu_values[j]));
    const Stokvec e     = (Muelmat::id() - rtepack::fresnel_reflectance(rv, rh)) * planck_ground;
    radiance[j]         = e.view()[stokes];
  }
}

void specular_surface_layer(ConstMatrixView ground_reflec, Tensor5View reflect, Tensor5View trans, Tensor3View source) {
  const Index nummu   = reflect.extent(1);
  const Index nstokes = ground_reflec.extent(0);
  ARTS_USER_ERROR_IF(ground_reflec.shape() != (std::array<Index, 2>{nstokes, nstokes}) or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "SPECULAR_SURFACE needs a square ground_reflec [nstokes, nstokes], reflect and trans [2, nummu, "
                     "nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {:B,}, {:B,}, {:B,} and {:B,}",
                     ground_reflec.shape(),
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));

  // REFLECT(S1, J, S2, J, 2) = GROUND_REFLEC(S2, S1) = R(S1, S2)
  for (Index j = 0; j < nummu; j++) reflect[1, j, joker, j, joker] = transpose(ground_reflec);
}

void specular_radiance(ConstMatrixView ground_reflec, Numeric ground_temp, Numeric frequency, MatrixView radiance) {
  const Index nstokes = ground_reflec.extent(0);
  ARTS_USER_ERROR_IF(
      nstokes > 4 or ground_reflec.shape() != (std::array<Index, 2>{nstokes, nstokes}) or radiance.extent(1) != nstokes,
      "SPECULAR_RADIANCE needs a square ground_reflec [nstokes, nstokes], nstokes <= 4, and radiance "
      "[nummu, nstokes]; got {:B,} and {:B,}",
      ground_reflec.shape(),
      radiance.shape());
  ARTS_USER_ERROR_IF(
      not(ground_temp >= 0.0), "SPECULAR_RADIANCE needs a ground temperature >= 0 K, got {} K", ground_temp);

  // Thermal radiation going up: (1 - R) B
  const Range stokes{0, nstokes};
  Muelmat     r{0.0};
  r.view()[stokes, stokes] = ground_reflec;
  const Stokvec e          = (Muelmat::id() - r) * Stokvec{planck(frequency, ground_temp)};

  radiance = 0.0;
  for (Index s = 0; s < nstokes; s++) radiance[joker, s] = e[s];
}

void external_surface_layer(ConstTensor4View surf_reflect, Tensor5View reflect, Tensor5View trans, Tensor3View source) {
  const Index nummu   = surf_reflect.extent(0);
  const Index nstokes = surf_reflect.extent(1);
  ARTS_USER_ERROR_IF(surf_reflect.shape() != (std::array<Index, 4>{nummu, nstokes, nummu, nstokes}) or
                         reflect.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         trans.shape() != (std::array<Index, 5>{2, nummu, nstokes, nummu, nstokes}) or
                         source.shape() != (std::array<Index, 3>{2, nummu, nstokes}),
                     "EXTERNAL_SURFACE needs surf_reflect [nummu, nstokes, nummu, nstokes], reflect and trans [2, "
                     "nummu, nstokes, nummu, nstokes] and source [2, nummu, nstokes]; got {:B,}, {:B,}, {:B,} and "
                     "{:B,}",
                     surf_reflect.shape(),
                     reflect.shape(),
                     trans.shape(),
                     source.shape());

  reflect = 0.0;
  source  = 0.0;
  for (Index h = 0; h < 2; h++) identity(trans[h].view_as(nummu * nstokes, nummu * nstokes));

  // REFLECT(I1, J1, I2, J2, 2) = SURF_REFL(I1, J1, I2, J2)
  reflect[1] = surf_reflect;
}

void thermal_radiance(Numeric temperature, Numeric albedo, Numeric frequency, Tensor3View radiance) {
  ARTS_USER_ERROR_IF(
      radiance.extent(0) != 2, "THERMAL_RADIANCE needs radiance [2, nummu, nstokes]; got {:B,}", radiance.shape());
  ARTS_USER_ERROR_IF(not(temperature >= 0.0), "THERMAL_RADIANCE needs a temperature >= 0 K, got {} K", temperature);

  radiance                  = 0.0;
  const Numeric thermal     = (1.0 - albedo) * planck(frequency, temperature);
  radiance[joker, joker, 0] = thermal;
}

void ground_surface(const surface&  ground,
                    ConstVectorView mu_values,
                    ConstVectorView quad_weights,
                    Numeric         frequency,
                    Numeric         ground_temp,
                    Tensor4View     surf_reflect,
                    MatrixView      gnd_radiance) {
  const Index nummu   = mu_values.extent(0);
  const Index nstokes = gnd_radiance.extent(1);
  ARTS_USER_ERROR_IF(quad_weights.extent(0) != nummu or
                         surf_reflect.shape() != (std::array<Index, 4>{nummu, nstokes, nummu, nstokes}) or
                         gnd_radiance.shape() != (std::array<Index, 2>{nummu, nstokes}),
                     "The ground with {} mu_values needs quad_weights [nummu], surf_reflect [nummu, nstokes, nummu, "
                     "nstokes] and gnd_radiance [nummu, nstokes]; got {}, {:B,} and {:B,}",
                     nummu,
                     quad_weights.extent(0),
                     surf_reflect.shape(),
                     gnd_radiance.shape());
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "The ground needs a positive frequency, got {} Hz", frequency);

  // The ground routines make the whole surface layer, of which
  // SURF_REFLECT is REFLECT(..., 2)
  Tensor5 reflect(2, nummu, nstokes, nummu, nstokes), trans(2, nummu, nstokes, nummu, nstokes);
  Tensor3 source(2, nummu, nstokes);

  std::visit(
      [&](const auto& g) {
        using T = std::remove_cvref_t<decltype(g)>;
        if constexpr (std::is_same_v<T, lambertian_surface>) {
          // For a Lambertian surface
          lambert_surface_layer(0, mu_values, quad_weights, g.albedo, reflect, trans, source);
          // The radiance from the ground is thermal and reflected direct
          lambert_radiance(g.albedo, ground_temp, frequency, gnd_radiance);
          surf_reflect = reflect[1];
        } else if constexpr (std::is_same_v<T, fresnel_surface>) {
          // For a Fresnel surface
          fresnel_surface_layer(mu_values, g.refractive_index, reflect, trans, source);
          // The radiance from the ground is thermal
          fresnel_radiance(mu_values, g.refractive_index, ground_temp, frequency, gnd_radiance);
          surf_reflect = reflect[1];
        } else if constexpr (std::is_same_v<T, specular_surface>) {
          // For a Specular surface
          ARTS_USER_ERROR_IF(g.reflectivity.shape() != (std::array<Index, 2>{nstokes, nstokes}),
                             "specular_surface reflectivity must be [{}, {}] (nstokes = {}), got {:B,}",
                             nstokes,
                             nstokes,
                             nstokes,
                             g.reflectivity.shape());
          specular_surface_layer(g.reflectivity, reflect, trans, source);
          // The radiance from the ground is thermal and reflected direct
          specular_radiance(g.reflectivity, ground_temp, frequency, gnd_radiance);
          surf_reflect = reflect[1];
        } else {
          // The reflection and emission given: SURF_REFLECT(out s, out mu,
          // in s, in mu) is [in mu, in s, out mu, out s]
          static_assert(std::is_same_v<T, discrete_surface>);
          ARTS_USER_ERROR_IF(g.reflection.shape() != (std::array<Index, 4>{nummu, nummu, nstokes, nstokes}) or
                                 g.emission.shape() != (std::array<Index, 2>{nummu, nstokes}),
                             "discrete_surface must have reflection [{}, {}, {}, {}] and emission [{}, {}] "
                             "(nmu_total = {}, nstokes = {}); got {:B,} and {:B,}",
                             nummu,
                             nummu,
                             nstokes,
                             nstokes,
                             nummu,
                             nstokes,
                             nummu,
                             nstokes,
                             g.reflection.shape(),
                             g.emission.shape());
          gnd_radiance = g.emission;
          for (Index io = 0; io < nummu; io++)
            for (Index ii = 0; ii < nummu; ii++) surf_reflect[ii, joker, io, joker] = transpose(g.reflection[io, ii]);
        }
      },
      ground);
}
}  // namespace rt4

#include "rt3_arts.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>
#include <polradtran_arts.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <utility>

namespace polradtran::rt3 {
namespace {
scattering_set no_scattering(Index degree) {
  scattering_set s{.extinction = 0.0, .scattering = 0.0, .legendre = CompactPlanarMuelmatVector(degree + 1)};
  s.legendre[0].F11() = 1.0;
  return s;
}
}  // namespace

scattering_set scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                                 const AtmPoint&                 atm_point,
                                 Numeric                         frequency,
                                 Index                           degree,
                                 Numeric                         normalisation_tolerance) {
  ARTS_USER_ERROR_IF(degree < 0, "The Legendre degree must be >= 0, got {}", degree);
  ARTS_USER_ERROR_IF(
      not(normalisation_tolerance >= 0.0), "normalisation_tolerance must be >= 0, got {}", normalisation_tolerance);
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "frequency must be positive, got {} Hz", frequency);

  if (scattering_species.species.empty()) return no_scattering(degree);

  const auto bulk =
      scattering_species.get_bulk_scattering_properties_tro_spectral(atm_point, Vector{frequency}, degree);
  const Numeric extinction = bulk.extinction_matrix[0].A();
  const Numeric scattering = extinction - bulk.absorption_vector[0][0];

  // c_l = a_l sqrt((2 l + 1) / 4 pi); each a_l is the scattering-plane Mueller matrix
  CompactPlanarMuelmatVector c(degree + 1);
  for (Index l = 0; l <= degree; l++) {
    const auto&   a = (*bulk.phase_matrix)[0, l];
    const Numeric y = std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * Constant::pi));
    c[l]            = y * CompactPlanarMuelmat{
                              a[0, 0].real(), a[0, 1].real(), a[1, 1].real(), a[2, 2].real(), a[2, 3].real(), a[3, 3].real()};
  }

  // The scattering coefficient implied by the phase matrix, 2 pi int F11 dx
  const Numeric phase_integral = 4.0 * Constant::pi * c[0].F11();
  ARTS_USER_ERROR_IF(not std::isinf(normalisation_tolerance) and
                         not(std::abs(phase_integral - scattering) <= normalisation_tolerance * extinction),
                     "The scattering coefficient from the phase matrix, 2 pi int F11 dcos(Theta) = {} per m, and the "
                     "extinction minus the absorption, {} per m, must agree to normalisation_tolerance times the "
                     "extinction, {} * {} per m (RT3 normalises the phase function and takes the albedo from the "
                     "latter)",
                     phase_integral,
                     scattering,
                     normalisation_tolerance,
                     extinction);

  if (phase_integral == 0.0) {
    auto s       = no_scattering(degree);
    s.extinction = extinction;
    return s;
  }

  const Numeric c0 = c[0].F11();
  for (auto& cl : c) cl /= c0;
  return {.extinction = extinction, .scattering = scattering, .legendre = std::move(c)};
}

problem problem_from_path(const ArrayOfPropagationPathPoint& ray_path,
                          const ArrayOfAtmPoint&             atm_path,
                          const ArrayOfPropmatVector&        spectral_propmat_path,
                          const AscendingGrid&               freq_grid,
                          Index                              freq_index,
                          const ArrayOfScatteringSpecies&    scattering_species,
                          const path_settings&               settings,
                          const surface&                     ground,
                          Numeric                            surface_temperature,
                          Numeric                            sky_temperature) {
  path_layers layers = layers_from_path(ray_path, atm_path, spectral_propmat_path, freq_grid, freq_index);
  const Index nlev   = layers.height.size();

  const Index nmu_total = settings.nmu + static_cast<Index>(settings.extra_mu.size());
  Index       degree    = settings.legendre_degree;
  if (degree < 0) {
    degree = max_legendre_degree(settings.nmu, settings.quad);
    if (settings.delta_m) degree = std::max(degree, 2 * nmu_total);
  }

  const Index nlay = nlev - 1;
  problem     p{.nstokes                = settings.nstokes,
                .nmu                    = settings.nmu,
                .quad                   = settings.quad,
                .extra_mu               = settings.extra_mu,
                .aziorder               = settings.aziorder,
                .max_delta_tau          = settings.max_delta_tau,
                .delta_m                = settings.delta_m,
                .direct_flux            = 0.0,
                .direct_mu              = 1.0,
                .thermal                = true,
                .frequency              = layers.frequency,
                .height                 = std::move(layers.height),
                .temperature            = std::move(layers.temperature),
                .gas_extinction         = std::move(layers.gas_extinction),
                .scattering_sets        = {},
                .layer_scattering_index = ArrayOfIndex(nlay, -1),
                .sky_temperature        = sky_temperature,
                .surface_temperature    = surface_temperature,
                .ground                 = ground};

  std::vector<scattering_set> level;
  level.reserve(nlev);
  for (const auto& atm : atm_path)
    level.push_back(scattering_optics(scattering_species, atm, p.frequency, degree, settings.normalisation_tolerance));

  for (Index l = 0; l < nlay; l++) {
    const auto& a = level[l];
    const auto& b = level[l + 1];

    scattering_set s{.extinction = std::midpoint(a.extinction, b.extinction),
                     .scattering = std::midpoint(a.scattering, b.scattering),
                     .legendre   = a.legendre};
    if (s.extinction == 0.0) continue;
    if (a.scattering + b.scattering != 0.0) {
      const Numeric wa = a.scattering / (a.scattering + b.scattering),
                    wb = b.scattering / (a.scattering + b.scattering);
      for (Index i = 0; i <= degree; i++) s.legendre[i] = wa * s.legendre[i] + wb * b.legendre[i];
    }
    p.layer_scattering_index[l] = static_cast<Index>(p.scattering_sets.size());
    p.scattering_sets.push_back(std::move(s));
  }
  return p;
}
}  // namespace polradtran::rt3

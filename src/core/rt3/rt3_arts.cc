#include "rt3_arts.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <debug.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <numeric>
#include <utility>

namespace rt3 {
namespace {
//! ARTS's compact TRO element of each RT3 column (F11, F12, F33, F34, F22, F44)
constexpr std::array<Index, 6> arts_element{0, 1, 3, 4, 2, 5};

scattering_set no_scattering(Index degree) {
  scattering_set s{.extinction = 0.0, .scattering = 0.0, .legendre = Matrix(degree + 1, 6, 0.0)};
  s.legendre[0, 0] = 1.0;
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

  // c_l = a_l sqrt((2 l + 1) / 4 pi), in ARTS's element order; each a_l is the scattering-plane Mueller matrix
  Matrix c(degree + 1, 6);
  for (Index l = 0; l <= degree; l++) {
    const auto&   a = (*bulk.phase_matrix)[0, l];
    const Numeric y = std::sqrt(static_cast<Numeric>(2 * l + 1) / (4.0 * Constant::pi));
    c[l, 0]         = y * a[0, 0].real();
    c[l, 1]         = y * a[0, 1].real();
    c[l, 2]         = y * a[1, 1].real();
    c[l, 3]         = y * a[2, 2].real();
    c[l, 4]         = y * a[2, 3].real();
    c[l, 5]         = y * a[3, 3].real();
  }

  // The scattering coefficient implied by the phase matrix, 2 pi int F11 dx
  const Numeric phase_integral = 4.0 * Constant::pi * c[0, 0];
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

  scattering_set s{.extinction = extinction, .scattering = scattering, .legendre = Matrix(degree + 1, 6)};
  for (Index l = 0; l <= degree; l++)
    for (Index k = 0; k < 6; k++) s.legendre[l, k] = c[l, arts_element[k]] / c[0, 0];
  return s;
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
  const Index nlev = static_cast<Index>(ray_path.size());
  const Index nf   = static_cast<Index>(freq_grid.size());

  ARTS_USER_ERROR_IF(nlev < 2, "ray_path needs at least 2 points (1 layer), got {}", nlev);
  ARTS_USER_ERROR_IF(
      static_cast<Index>(atm_path.size()) != nlev or static_cast<Index>(spectral_propmat_path.size()) != nlev,
      "ray_path, atm_path and spectral_propmat_path must have one entry per level; they have {}, {} "
      "and {}",
      nlev,
      atm_path.size(),
      spectral_propmat_path.size());
  ARTS_USER_ERROR_IF(freq_index < 0 or freq_index >= nf,
                     "freq_index must be in [0, {}) for a freq_grid of {} frequencies, got {}",
                     nf,
                     nf,
                     freq_index);
  ARTS_USER_ERROR_IF(
      stdr::any_of(spectral_propmat_path, [nf](const PropmatVector& v) { return static_cast<Index>(v.size()) != nf; }),
      "Every spectral_propmat_path level must have freq_grid.size() = {} propagation matrices",
      nf);
  for (Index l = 0; l < nlev - 1; l++)
    ARTS_USER_ERROR_IF(not(ray_path[l].altitude() > ray_path[l + 1].altitude()),
                       "The ray_path altitudes must decrease strictly from the first point (top of the atmosphere) "
                       "to the last (surface)");
  ARTS_USER_ERROR_IF(stdr::any_of(spectral_propmat_path,
                                  [freq_index](const PropmatVector& v) { return v[freq_index].is_polarized(); }),
                     "RT3's gas extinction is scalar: the gas propagation matrices in spectral_propmat_path must not "
                     "be polarized (only A may be non-zero) at frequency index {}",
                     freq_index);

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
                .frequency              = freq_grid[freq_index],
                .height                 = Vector(nlev),
                .temperature            = Vector(nlev),
                .gas_extinction         = Vector(nlay),
                .scattering_sets        = {},
                .layer_scattering_index = ArrayOfIndex(nlay, -1),
                .sky_temperature        = sky_temperature,
                .surface_temperature    = surface_temperature,
                .ground                 = ground};

  for (Index l = 0; l < nlev; l++) {
    p.height[l]      = ray_path[l].altitude();
    p.temperature[l] = atm_path[l].temperature;
  }
  for (Index l = 0; l < nlay; l++)
    p.gas_extinction[l] =
        std::midpoint(spectral_propmat_path[l][freq_index].A(), spectral_propmat_path[l + 1][freq_index].A());

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
      s.legendre *= a.scattering / (a.scattering + b.scattering);
      for (Index i = 0; i <= degree; i++)
        for (Index k = 0; k < 6; k++)
          s.legendre[i, k] += b.scattering / (a.scattering + b.scattering) * b.legendre[i, k];
    }
    p.layer_scattering_index[l] = static_cast<Index>(p.scattering_sets.size());
    p.scattering_sets.push_back(std::move(s));
  }
  return p;
}
}  // namespace rt3

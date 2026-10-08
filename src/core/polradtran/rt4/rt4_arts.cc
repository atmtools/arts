#include "rt4_arts.h"

#include <arts_conversions.h>
#include <debug.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>
#include <utility>
#include <vector>

namespace polradtran::rt4 {
namespace {
bool all_zero(const layer_optics& o) {
  const auto zero = [](Numeric x) { return x == 0.0; };
  return stdr::all_of(o.extinction | by_elem, zero) and stdr::all_of(o.absorption | by_elem, zero) and
         stdr::all_of(o.phase | by_elem, zero);
}
}  // namespace

layer_optics scattering_optics(const ArrayOfScatteringSpecies& scattering_species,
                               const AtmPoint&                 atm_point,
                               Numeric                         frequency,
                               const Vector&                   mu,
                               Index                           nstokes) {
  const Index ns = nstokes;
  const Index n  = static_cast<Index>(mu.size());

  ARTS_USER_ERROR_IF(ns != 1 and ns != 2, "RT4 supports nstokes 1 ([I]) or 2 ([I, Q]), got {}", ns);
  ARTS_USER_ERROR_IF(n < 1, "RT4 needs at least one stream per hemisphere, got an empty mu");
  ARTS_USER_ERROR_IF(stdr::any_of(mu, [](Numeric x) { return not(x > 0.0 and x <= 1.0); }),
                     "RT4 stream cosines mu must be in (0, 1]");
  ARTS_USER_ERROR_IF(not(frequency > 0.0), "frequency must be positive, got {} Hz", frequency);

  layer_optics o{.extinction = Tensor4(2, n, ns, ns, 0.0),
                 .absorption = Tensor3(2, n, ns, 0.0),
                 .phase      = Tensor6(2, 2, n, n, ns, ns, 0.0)};
  if (scattering_species.species.empty()) return o;

  // The propagation zenith angles [deg] of the streams, ascending as ARTS's
  // angular grids must be: RT4's (h, mu_i) is entry stream[h * n + i].
  // Downward streams propagate toward the surface.
  std::vector<std::pair<Numeric, Index>> sorted(2 * n);
  for (Index i = 0; i < n; i++) {
    const Numeric theta  = Conversion::rad2deg(std::acos(mu[i]));
    sorted[down * n + i] = {180.0 - theta, down * n + i};
    sorted[up * n + i]   = {theta, up * n + i};
  }
  stdr::sort(sorted);
  Vector       za(2 * n);
  ArrayOfIndex stream(2 * n);
  for (Index j = 0; j < 2 * n; j++) {
    za[j]                    = sorted[j].first;
    stream[sorted[j].second] = j;
  }
  ARTS_USER_ERROR_IF(not ZenGrid::is_sorted(za), "RT4 stream cosines mu must be distinct");

  const auto bulk =
      scattering_species.get_bulk_scattering_properties_aro_fourier(atm_point, Vector{frequency}, za, za, 0);
  ARTS_USER_ERROR_IF(not bulk.phase_matrix.has_value(),
                     "RT4 needs the phase matrix of every scattering species; the bulk scattering properties have "
                     "none");

  const auto& ext = bulk.extinction_matrix;
  const auto& abs = bulk.absorption_vector;
  const auto& pha = *bulk.phase_matrix;
  ARTS_USER_ERROR_IF(ext.extent(0) != 1 or abs.extent(0) != 1 or pha.extent(0) != 1,
                     "The bulk scattering properties must be at a single temperature; they have {}, {} and {} "
                     "temperatures for the extinction, absorption and phase matrix",
                     ext.extent(0),
                     abs.extent(0),
                     pha.extent(0));

  for (Index h = 0; h < 2; h++) {
    for (Index i = 0; i < n; i++) {
      const Index d = stream[h * n + i];
      for (Index s = 0; s < ns; s++) {
        o.extinction[h, i, s, s] = ext[0, 0, d, 0];
        o.absorption[h, i, s]    = abs[0, 0, d, s];
      }
      if (ns > 1) o.extinction[h, i, 0, 1] = o.extinction[h, i, 1, 0] = ext[0, 0, d, 1];
    }
  }

  // The azimuthal mean C_0 [t, f, za_inc, za_scat, m = 0, cosine, 4 * row + col]
  for (Index ho = 0; ho < 2; ho++) {
    for (Index hi = 0; hi < 2; hi++) {
      for (Index io = 0; io < n; io++) {
        for (Index ii = 0; ii < n; ii++) {
          const Index out = stream[ho * n + io], in = stream[hi * n + ii];
          for (Index so = 0; so < ns; so++)
            for (Index si = 0; si < ns; si++) o.phase[ho, hi, io, ii, so, si] = pha[0, 0, in, out, 0, 0, 4 * so + si];
        }
      }
    }
  }
  return o;
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
                     "RT4's gas extinction is scalar: the gas propagation matrices in spectral_propmat_path must not "
                     "be polarized (only A may be non-zero) at frequency index {}",
                     freq_index);

  const Index nlay = nlev - 1;
  problem     p{.nstokes                 = settings.nstokes,
                .nmu                     = settings.nmu,
                .quad                    = settings.quad,
                .extra_mu                = settings.extra_mu,
                .max_delta_tau           = settings.max_delta_tau,
                .normalisation_tolerance = settings.normalisation_tolerance,
                .frequency               = freq_grid[freq_index],
                .height                  = Vector(nlev),
                .temperature             = Vector(nlev),
                .gas_extinction          = Vector(nlay),
                .optics                  = {},
                .layer_optics_index      = ArrayOfIndex(nlay, -1),
                .sky_temperature         = sky_temperature,
                .surface_temperature     = surface_temperature,
                .ground                  = ground};

  for (Index l = 0; l < nlev; l++) {
    p.height[l]      = ray_path[l].altitude();
    p.temperature[l] = atm_path[l].temperature;
  }
  for (Index l = 0; l < nlay; l++)
    p.gas_extinction[l] =
        std::midpoint(spectral_propmat_path[l][freq_index].A(), spectral_propmat_path[l + 1][freq_index].A());

  const auto q = get_quadrature(settings.nmu, settings.quad);
  Vector     mu(settings.nmu + static_cast<Index>(settings.extra_mu.size()));
  std::ranges::copy(q.mu, mu.begin());
  std::ranges::copy(settings.extra_mu, mu.begin() + settings.nmu);

  std::vector<layer_optics> level;
  level.reserve(nlev);
  for (const auto& atm : atm_path)
    level.push_back(scattering_optics(scattering_species, atm, p.frequency, mu, settings.nstokes));

  for (Index l = 0; l < nlay; l++) {
    layer_optics o  = level[l];
    o.extinction   += level[l + 1].extinction;
    o.absorption   += level[l + 1].absorption;
    o.phase        += level[l + 1].phase;
    o.extinction   *= 0.5;
    o.absorption   *= 0.5;
    o.phase        *= 0.5;
    if (all_zero(o)) continue;
    p.layer_optics_index[l] = static_cast<Index>(p.optics.size());
    p.optics.push_back(std::move(o));
  }
  return p;
}
}  // namespace polradtran::rt4

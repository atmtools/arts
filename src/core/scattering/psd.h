#ifndef ARTS_CORE_SCATTERING_PSD_H_
#define ARTS_CORE_SCATTERING_PSD_H_

#include <matpack.h>
#include <optional>

#include "atm.h"
#include "properties.h"
#include "utils.h"

namespace scattering {

enum class SizeParameter {
    Mass,
    MaximumDiameter,
    VolumeEqDiameter
};


/*** Single-moment modified gamma distribution
 *
 * Implements a modified gamma distribution with a single free moment.
 * Currently this moment is expected to be the mass density.
 */
struct MGDSingleMoment {
  ScatteringSpeciesProperty moment;
  Numeric n_alpha;
  Numeric n_b;
  Numeric mu;
  Numeric gamma;
  Numeric t_min;
  Numeric t_max;
  bool picky;

  MGDSingleMoment() = default;

  MGDSingleMoment(ScatteringSpeciesProperty moment_,
                  Numeric n_alpha_,
                  Numeric n_b_,
                  Numeric mu_,
                  Numeric gamma_,
                  Numeric t_min_,
                  Numeric t_max_,
                  bool picky_)
      : moment(moment_),
        n_alpha(n_alpha_),
        n_b(n_b_),
        mu(mu_),
        gamma(gamma_),
        t_min(t_min_),
        t_max(t_max_),
    picky(picky_)
  {}

  MGDSingleMoment(ScatteringSpeciesProperty moment_,
                  std::string name,
                  Numeric t_min_,
                  Numeric t_max_,
                  bool picky_)
      : moment(moment_), t_min(t_min_), t_max(t_max_), picky(picky_) {
    if (name == "Abel12") {
      n_alpha = 0.22;
      n_b = 2.2;
      mu = 0.0;
      gamma = 1.0;
    } else if (name == "Wang16") {
      // Wang 16 parameters converted to SI units
      n_alpha = 14.764;
      n_b = 1.49;
      mu = 0.0;
      gamma = 1.0;
    } else if (name == "Field19") {
      n_alpha = 7.9e9;
      n_b = -2.58;
      mu = 0.0;
      gamma = 1.0;
    } else {
      std::ostringstream os;
      os << "The PSD configuration '" << name << "' is currently not supported."
         << " Supported config names are 'Abel12', 'Wang16', 'Field19'.";
      throw std::runtime_error(os.str());
    }
  }

  /** Evaluate PSD at given atmospheric point.
   *
   * @param point The atmospheric point.
   * @param particle_size A vector containing the particles sizes at which to
   * evaluate the PSD.
   * @param scat_species_a The a parameter of the mass-size relationship of
   * the particle data.
   * @param scat_species_b The b parameter of the mass-size relationship of
   * the particle data.
   *
   */
  Vector evaluate(const AtmPoint& point,
                  const Vector& particle_sizes,
                  const Numeric& scat_species_a,
                  const Numeric& scat_species_b) const;
};


/*** Binned PSD
 *
 * The BinnedPSD class represents a particle size distribution using a fixed number of particle counts
 * over a sequence of size bins. Particles with sizes outside the size and temperature range are set
 * to zero.
 */
struct BinnedPSD {
  Vector bins;
  Vector counts;
  Numeric t_min = 0.0;
  Numeric t_max = 350.0;

  BinnedPSD() = default;

  BinnedPSD(Vector bins_,
            Vector counts_,
            Numeric t_min_ = 0.0,
            Numeric t_max_ = 350.0):
  bins(bins_), counts(counts_), t_min(t_min_), t_max(t_max_) {
    if (bins_.size() != (counts_.size() + 1)) {
        ARTS_USER_ERROR("The bin vector must have exactly one element more than the counts vector.");
    }
    if (!std::is_sorted(bins.begin(), bins.end())) {
        ARTS_USER_ERROR("The bins vector must be strictly increasing.");
    }
  }

  Vector evaluate(const AtmPoint& point,
                  const Vector& particle_sizes,
                  const Numeric& scat_species_a,
                  const Numeric& scat_species_b) const {

    Index n_parts = particle_sizes.size();
    Vector pnd = Vector(n_parts);
    for (Index ind = 0; ind < n_parts; ++ind) {
      if ((point.temperature < t_min) || (t_max < point.temperature)) {
        pnd[ind] = 0.0;
      } else {
        Index bin_ind = digitize(bins, particle_sizes[ind]);
        if (bin_ind < 0) {
          pnd[ind] = 0.0;
        } else if (bin_ind >= bins.size()) {
          pnd[ind] = 0.0;
        } else {
          pnd[ind] = counts[bin_ind];
        }
      }
    }
    return pnd;
  }
};

}
#endif  // ARTS_CORE_PSD_H_

#ifndef ARTS_CORE_SCATTERING_PSD_H_
#define ARTS_CORE_SCATTERING_PSD_H_

#include <matpack.h>
#include <optional>

#include "atm.h"
#include "enums.h"
#include "properties.h"
#include "utils.h"

namespace scattering {


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
                  bool picky_);

  MGDSingleMoment(ScatteringSpeciesProperty moment_,
                  std::string name,
                  Numeric t_min_,
                  Numeric t_max_,
                  bool picky_);

  static constexpr SizeParameter get_size_parameter() {
    return SizeParameter::DVeq;
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

  SizeParameter size_parameter = SizeParameter::DVeq;
  Vector bins;
  Vector counts;
  Numeric t_min = 0.0;
  Numeric t_max = 350.0;

  BinnedPSD() = default;

  BinnedPSD(SizeParameter size_parameter_,
            Vector bins_,
            Vector counts_,
            Numeric t_min_ = 0.0,
            Numeric t_max_ = 350.0);

  static constexpr SizeParameter get_size_parameter() {
    return SizeParameter::Mass;
  }

  Vector evaluate(const AtmPoint& point,
                  const Vector& particle_sizes,
                  const Numeric& /*scat_species_a*/,
                  const Numeric& /*scat_species_b*/) const;
};

}
#endif  // ARTS_CORE_PSD_H_

#ifndef ARTS_CORE_SCATTERING_PSD_H_
#define ARTS_CORE_SCATTERING_PSD_H_

#include <math_funcs.h>
#include <matpack_math.h>

#include <optional>

#include "atm.h"
#include "properties.h"

using Parameter = std::variant<Numeric, ScatteringSpeciesProperty>;

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
      : n_alpha(n_alpha_),
        n_b(n_b_),
        mu(mu_),
        gamma(gamma_),
        t_min(t_min_),
        t_max(t_max_) {}

  MGDSingleMoment(ScatteringSpeciesProperty moment_,
                  std::string name,
                  Numeric t_min_,
                  Numeric t_max_,
                  bool picky_)
      : t_min(t_min_), t_max(t_max_), picky(picky_) {
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

  /** Evaluate PSD a atmospheric point.
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
                  const Numeric& scat_species_b) const {
    if (!point.has(moment)) {
      std::ostringstream os;
      os << "The PSD requires water content from '" << moment
         << "' but it is not part of the AtmPoint.";
      throw std::runtime_error(os.str());
    }

    Numeric water_content = point[moment];
    Numeric t = point[Atm::Key::t];

    // Outside of [t_min,tmax]?
    if (t < t_min || t > t_max) {
      if (picky) {
        std::ostringstream os;
        os << "Method called with a temperature of " << t << " K.\n"
           << "This is outside the specified allowed range: [ max(0.," << t_min
           << "), " << t_max << " ]";
        throw std::runtime_error(os.str());
      }
    }

    // Negative wc?
    Numeric psd_weight = 1.0;
    if (water_content < 0) {
      psd_weight = -1.0;
      water_content *= -1.0;
    }

    auto nsi = particle_sizes.size();

    // Calculate PSD
    // Calculate lambda for modified gamma distribution from mass density
    Numeric k = (scat_species_b + mu + 1 - gamma) / gamma;
    Numeric expo = 1.0 / (n_b - k - 1);
    Numeric denom = scat_species_a * n_alpha * tgamma(k + 1);
    Numeric lam = pow(water_content * gamma / denom, expo);
    Numeric n_0 = n_alpha * pow(lam, n_b);
    Vector psd_1p(nsi);
    Matrix jac_data(4, nsi);

    psd_1p = 0.0;
    mgd_with_derivatives(psd_1p,
                         jac_data,
                         particle_sizes,
                         n_0,
                         mu,
                         lam,
                         gamma,
                         false,   // n_0 jacobian
                         false,   // mu jacobian
                         false,   // lambda jacobian
                         false);  // gamma jacobian
    return psd_1p;
  }
};

#endif  // ARTS_CORE_PSD_H_

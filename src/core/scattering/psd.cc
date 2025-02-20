#include "psd.h"

namespace scattering {

Vector MGDSingleMoment::evaluate(const AtmPoint& point,
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
  Numeric t = point[AtmKey::t];

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
  if (water_content < 0) {
    water_content *= -1.0;
  }

  auto nsi = particle_sizes.size();
  Vector psd(nsi);
  psd = 0.0;
  if (water_content == 0.0) return psd;

  // Calculate PSD
  // Calculate lambda for modified gamma distribution from mass density
  Numeric k = (scat_species_b + mu + 1 - gamma) / gamma;
  Numeric expo = 1.0 / (n_b - k - 1);
  Numeric denom = scat_species_a * n_alpha * tgamma(k + 1);
  Numeric lam = pow(water_content * gamma / denom, expo);
  Numeric n_0 = n_alpha * pow(lam, n_b);
  Matrix jac_data(4, nsi);

  mgd_with_derivatives(psd,
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
  return psd;
}

}

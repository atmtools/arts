#include "psd.h"

namespace scattering {
MGDSingleMoment::MGDSingleMoment(ScatteringSpeciesProperty moment_,
                                 Numeric n_alpha_,
                                 Numeric n_b_,
                                 Numeric mu_,
                                 Numeric gamma_,
                                 Numeric t_min_,
                                 Numeric t_max_,
                                 bool picky_)
    : moment(std::move(moment_)),
      n_alpha(n_alpha_),
      n_b(n_b_),
      mu(mu_),
      gamma(gamma_),
      t_min(t_min_),
      t_max(t_max_),
      picky(picky_) {}

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
  Numeric t             = point[AtmKey::t];

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
  Numeric k     = (scat_species_b + mu + 1 - gamma) / gamma;
  Numeric expo  = 1.0 / (n_b - k - 1);
  Numeric denom = scat_species_a * n_alpha * tgamma(k + 1);
  Numeric lam   = pow(water_content * gamma / denom, expo);
  Numeric n_0   = n_alpha * pow(lam, n_b);
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

MGDSingleMoment::MGDSingleMoment(ScatteringSpeciesProperty moment_,
                                 std::string name,
                                 Numeric t_min_,
                                 Numeric t_max_,
                                 bool picky_)
    : moment(std::move(moment_)), t_min(t_min_), t_max(t_max_), picky(picky_) {
  if (name == "Abel12") {
    n_alpha = 0.22;
    n_b     = 2.2;
    mu      = 0.0;
    gamma   = 1.0;
  } else if (name == "Wang16") {
    // Wang 16 parameters converted to SI units
    n_alpha = 14.764;
    n_b     = 1.49;
    mu      = 0.0;
    gamma   = 1.0;
  } else if (name == "Field19") {
    n_alpha = 7.9e9;
    n_b     = -2.58;
    mu      = 0.0;
    gamma   = 1.0;
  } else {
    std::ostringstream os;
    os << "The PSD configuration '" << name << "' is currently not supported."
       << " Supported config names are 'Abel12', 'Wang16', 'Field19'.";
    throw std::runtime_error(os.str());
  }
}

BinnedPSD::BinnedPSD(SizeParameter size_parameter_,
                     Vector bins_,
                     Vector counts_,
                     Numeric t_min_,
                     Numeric t_max_)
    : size_parameter(size_parameter_),
      bins(bins_),
      counts(counts_),
      t_min(t_min_),
      t_max(t_max_) {
  if (bins_.size() != (counts_.size() + 1)) {
    ARTS_USER_ERROR(
        "The bin vector must have exactly one element more than the counts vector.");
  }
  if (!std::is_sorted(bins.begin(), bins.end())) {
    ARTS_USER_ERROR("The bins vector must be strictly increasing.");
  }
}

Vector BinnedPSD::evaluate(const AtmPoint& point,
                           const Vector& particle_sizes,
                           const Numeric& /*scat_species_a*/,
                           const Numeric& /*scat_species_b*/) const {
  Index n_parts = particle_sizes.size();
  Vector pnd    = Vector(n_parts);
  for (Index ind = 0; ind < n_parts; ++ind) {
    if ((point.temperature < t_min) || (t_max < point.temperature)) {
      pnd[ind] = 0.0;
    } else {
      Index bin_ind = digitize(bins, particle_sizes[ind]);
      Index n_bins  = bins.size();
      if (bin_ind < 0) {
        pnd[ind] = 0.0;
      } else if (bin_ind >= n_bins) {
        pnd[ind] = 0.0;
      } else {
        pnd[ind] = counts[bin_ind];
      }
    }
  }
  return pnd;
}
}  // namespace scattering

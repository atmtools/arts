#include "scattering_species.h"

namespace Scattering {

HenyeyGreenstein::HenyeyGreenstein(const Numeric& g_) : g(g_){};

Numeric HenyeyGreenstein::evaluate_phase_function(const Numeric& theta) {
  Numeric g2 = g * g;
  return 0.25 * Constant::inv_pi * (1.0 - g2) *
         std::pow(1.0 + g2 - 2.0 * g * std::cos(theta), -3.0 / 2.0);
}

Vector HenyeyGreenstein::evaluate_phase_function(const Vector& theta) {
  Vector results(theta.size());
  std::transform(theta.begin(),
                 theta.end(),
                 results.begin(),
                 [this](const Numeric& theta_) {
                   return evaluate_phase_function(theta_);
                 });
  return results;
}

std::ostream& operator<<(std::ostream& os,
                         const Scattering::HenyeyGreenstein& scatterer) {
  return os << "HenyeyGreenstein(g = " << scatterer.g << ")";
}

std::ostream& operator<<(std::ostream& os,
                         const Scattering::Species& /*species*/) {
  os << "A scattering species." << std::endl;
  return os;
}

}  // namespace Scattering

#include "lbl_lineshape_linemixing.h"

namespace lbl::linemixing {
Numeric species_data::Q(const Rational J,
                        const Numeric T,
                        const Numeric T0,
                        const Numeric energy) const {
  return std::exp(- beta(T0, T) * energy / (Constant::k * T)) * scaling(T0, T) / pow(J * (J+1), lambda(T0, T));
}

Numeric species_data::Omega(const Numeric T,
                            const Numeric T0,
                            const Numeric other_mass,
                            const Numeric energy_x,
                            const Numeric energy_xm2) const {
  using Constant::h;
  using Constant::k;
  using Constant::pi;
  using Constant::m_u;
  using Constant::h_bar;
  using Math::pow2;
  
  // Constants for the expression
  constexpr Numeric fac = 8 * k / (m_u * pi);
  
  const Numeric wnnm2 = (energy_x - energy_xm2) / h_bar;
  
  const Numeric inv_eff_mass = 1 / mass + 1 / other_mass;
  const Numeric v_bar_pow2 = fac*T*inv_eff_mass;
  const Numeric tauc_pow2 = pow2(collisional_distance(T0, T)) / v_bar_pow2;
  
  return 1.0 / pow2(1 + pow2(wnnm2) * tauc_pow2 / 24.0);
}

std::ostream& operator<<(std::ostream&os , const isot_map& self) {
  for (auto& [k, v]: self) {
    os << k << ":\n";
    for (auto& [k2, v2]: v) {
      os << "  " << k2 << ":\n";
      os << "    scaling:              " << v2.scaling << "\n";
      os << "    beta:                 " << v2.beta << "\n";
      os << "    lambda:               " << v2.lambda << "\n";
      os << "    collisional_distance: " << v2.collisional_distance << "\n";
      os << "    mass:                 " << v2.mass << "\n";
    }
  }
  return os;
}
}  // namespace lbl::linemixing
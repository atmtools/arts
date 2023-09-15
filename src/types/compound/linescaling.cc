#include "linescaling.h"
#include "partfun.h"

Numeric single_partition_function(const Numeric& T,
                                  const Species::IsotopeRecord& ir) {
  return PartitionFunctions::Q(T, ir);
}

Numeric dsingle_partition_function_dT(const Numeric& T,
                                      const Species::IsotopeRecord& ir) {
  return PartitionFunctions::dQdT(T, ir);
}

Numeric stimulated_emission(Numeric T, Numeric F0) {
  using namespace Constant;
  static constexpr Numeric c1 = -h / k;
  return std::exp(c1 * F0 / T);
}

Numeric dstimulated_emissiondT(Numeric T, Numeric F0) {
  using namespace Constant;
  using namespace Math;
  static constexpr Numeric c1 = -h / k;
  return -F0 * c1 * std::exp(F0 * c1 / T) / pow2(T);
}

Numeric dstimulated_emissiondF0(Numeric T, Numeric F0) {
  using namespace Constant;
  static constexpr Numeric c1 = -h / k;
  return c1 * std::exp(F0 * c1 / T) / T;
}

Numeric stimulated_relative_emission(const Numeric F0, const Numeric T0, const Numeric T) noexcept {
  return std::expm1(-Conversion::hz2joule(F0) / Conversion::kelvin2joule(T)) / std::expm1(-Conversion::hz2joule(F0) / Conversion::kelvin2joule(T0));
}

Numeric dstimulated_relative_emission_dT(const Numeric F0, const Numeric T0, const Numeric T) noexcept {
  return Conversion::hz2joule(F0) * std::exp(-Conversion::hz2joule(F0) / Conversion::kelvin2joule(T)) /
  (T * Conversion::kelvin2joule(T) * std::expm1(-Conversion::hz2joule(F0) / Conversion::kelvin2joule(T0)));
}

Numeric dstimulated_relative_emission_dF0(const Numeric F0, const Numeric T0, const Numeric T) noexcept {
  const Numeric exp0m1 = std::expm1(-Conversion::hz2joule(F0) / Conversion::kelvin2joule(T0));
  const Numeric exptm1 = std::expm1(-Conversion::hz2joule(F0) / Conversion::kelvin2joule(T));
  return Constant::h * (T * exptm1 * (1 + exp0m1) - T0 * exp0m1 * (1 + exptm1)) / (Constant::k * T * T0 * Math::pow2(exp0m1));
}

Numeric stimulated_relative_emission(const Numeric& gamma,
                                     const Numeric& gamma_ref) {
  return (1. - gamma) / (1. - gamma_ref);
}

Numeric dstimulated_relative_emission_dT(const Numeric& gamma,
                                         const Numeric& gamma_ref,
                                         const Numeric& F0,
                                         const Numeric& T) {
  static constexpr Numeric c = -Constant::h / Constant::k;

  return c * F0 * gamma / (T * T * (1. - gamma_ref));
}

Numeric dstimulated_relative_emission_dF0(const Numeric& gamma,
                                          const Numeric& gamma_ref,
                                          const Numeric& T,
                                          const Numeric& T0) {
  static constexpr Numeric c = -Constant::h / Constant::k;

  const Numeric g0 = 1 - gamma_ref;
  const Numeric g = 1 - gamma;

  return c * (g * gamma_ref / (T0 * g0 * g0) - gamma / (T * g0));
}

// Ratio of boltzman emission at T and T0
Numeric boltzman_ratio(const Numeric& T, const Numeric& T0, const Numeric& E0) {
  static constexpr Numeric c = 1 / Constant::k;

  return exp(E0 * c * (T - T0) / (T * T0));
}

Numeric dboltzman_ratio_dT(const Numeric& boltzmann_ratio,
                           const Numeric& T,
                           const Numeric& E0) {
  static constexpr Numeric c = 1 / Constant::k;

  return E0 * c * boltzmann_ratio / (T * T);
}

// Boltzmann factor at T
Numeric boltzman_factor(Numeric T, Numeric E0)
{
  return std::exp(- E0 / (Constant::k*T));
}

// Boltzmann factor at T
Numeric dboltzman_factordT(Numeric T, Numeric E0) {
  using namespace Constant;
  using namespace Math;
  static constexpr Numeric c1 = -1 / k;
  return -E0 * c1 * std::exp(E0 * c1 / T) / pow2(T);
}

// Boltzmann factor at T
Numeric dboltzman_factordE0(Numeric T, Numeric E0) {
  using namespace Constant;
  static constexpr Numeric c1 = -1 / k;
  return c1 * std::exp(E0 * c1 / T) / T;
}

Numeric absorption_nlte_ratio(const Numeric& gamma,
                              const Numeric& r_upp,
                              const Numeric& r_low) noexcept {
  return (r_low - r_upp * gamma) / (1 - gamma);
}

Numeric dabsorption_nlte_rate_dT(const Numeric& gamma,
                                 const Numeric& T,
                                 const Numeric& F0,
                                 const Numeric& El,
                                 const Numeric& Eu,
                                 const Numeric& r_upp,
                                 const Numeric& r_low) {
  static constexpr Numeric c = 1 / Constant::k;

  ARTS_USER_ERROR_IF (El < 0 or Eu < 0,
    "It is considered undefined behavior to NLTE and "
    "temperature Jacobian without defining all "
    "vibrational energy states")

  const Numeric x = 1 / (T * (gamma - 1));
  const Numeric hf = F0 * Constant::h;

  return x * x * c *
         ((gamma - 1) * (El * r_low - Eu * gamma * r_upp) -
          hf * gamma * (r_low - r_upp));
}

Numeric dabsorption_nlte_rate_dF0(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& r_upp,
                                  const Numeric& r_low) {
  static constexpr Numeric c = -Constant::h / Constant::k;

  return c * gamma * (r_low - r_upp) / (T * Math::pow2(gamma - 1));
}

Numeric dabsorption_nlte_rate_dTl(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tl,
                                  const Numeric& El,
                                  const Numeric& r_low) {
  const Numeric x = 1 / (Constant::k * T);
  const Numeric y = 1 / Tl;

  return El * x * y * y * T * r_low / (gamma - 1);
}

Numeric dabsorption_nlte_rate_dTu(const Numeric& gamma,
                                  const Numeric& T,
                                  const Numeric& Tu,
                                  const Numeric& Eu,
                                  const Numeric& r_upp) {
  const Numeric x = 1 / (Constant::k * T);
  const Numeric y = 1 / Tu;

  return Eu * x * y * y * T * gamma * r_upp / (gamma - 1);
}

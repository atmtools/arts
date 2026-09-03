#include "common.h"

#include <arts_constants.h>
#include <debug.h>
#include <legendre.h>

#include <algorithm>
#include <limits>
#include <ranges>
#include <stdexcept>

namespace disort_common {
namespace {
Numeric checked_surface_mu(const Numeric mu) {
  const Numeric value = std::abs(mu);
  if (not(value > 0.0 and value <= 1.0)) throw std::domain_error("Surface direction cosine must satisfy 0 < |mu| <= 1");
  return value;
}

Numeric unit_clamp(const Numeric value) { return std::clamp(value, -1.0, 1.0); }

Numeric shadow_eta(const Numeric mu, const Numeric slope_variance) {
  const Numeric sine = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  if (sine == 0.0) return 0.0;
  const Numeric cotangent = mu / sine;
  const Numeric root      = std::sqrt(slope_variance);
  return 0.5 * (root / (std::sqrt(Constant::pi) * cotangent) * std::exp(-cotangent * cotangent / slope_variance) -
                std::erfc(cotangent / root));
}
}  // namespace

Numeric lambertian_brdf(const Numeric albedo) {
  if (not std::isfinite(albedo) or albedo < 0.0 or albedo > 1.0)
    throw std::domain_error("Lambertian albedo must be finite and in [0, 1]");
  return albedo * Constant::inv_pi;
}

Numeric hapke_brdf(const Numeric outgoing_mu,
                   const Numeric incoming_mu,
                   const Numeric relative_azimuth,
                   const Numeric opposition_amplitude,
                   const Numeric opposition_width,
                   const Numeric single_scattering_albedo) {
  const Numeric mu  = checked_surface_mu(outgoing_mu);
  const Numeric mup = checked_surface_mu(incoming_mu);
  if (not std::isfinite(opposition_amplitude) or not std::isfinite(opposition_width) or
      not std::isfinite(single_scattering_albedo) or opposition_width < 0.0 or single_scattering_albedo < 0.0 or
      single_scattering_albedo > 1.0)
    throw std::domain_error("Invalid Hapke BRDF parameters");

  const Numeric sin_mu    = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  const Numeric sin_mup   = std::sqrt(std::max(0.0, 1.0 - mup * mup));
  const Numeric cos_alpha = unit_clamp(mu * mup - sin_mu * sin_mup * std::cos(relative_azimuth));
  const Numeric alpha     = std::acos(cos_alpha);
  const Numeric phase     = 1.0 + 0.5 * cos_alpha;
  const Numeric tangent   = std::tan(0.5 * alpha);
  const Numeric opposition =
      opposition_width == 0.0 ? 0.0 : opposition_amplitude * opposition_width / (opposition_width + tangent);
  const Numeric gamma = std::sqrt(1.0 - single_scattering_albedo);
  const Numeric h0    = (1.0 + 2.0 * mup) / (1.0 + 2.0 * mup * gamma);
  const Numeric h     = (1.0 + 2.0 * mu) / (1.0 + 2.0 * mu * gamma);
  return single_scattering_albedo / (4.0 * Constant::pi * (mu + mup)) * ((1.0 + opposition) * phase + h0 * h - 1.0);
}

Numeric rpv_brdf(const Numeric outgoing_mu,
                 const Numeric incoming_mu,
                 const Numeric relative_azimuth,
                 const Numeric rho0,
                 const Numeric kappa,
                 const Numeric asymmetry,
                 const Numeric hotspot) {
  const Numeric mu_r = checked_surface_mu(outgoing_mu);
  const Numeric mu_i = checked_surface_mu(incoming_mu);
  if (not std::isfinite(rho0) or not std::isfinite(kappa) or not std::isfinite(asymmetry) or
      not std::isfinite(hotspot) or rho0 < 0.0 or not(asymmetry > -1.0 and asymmetry < 1.0) or hotspot < 0.0)
    throw std::domain_error("Invalid RPV BRDF parameters");

  const Numeric sin_i     = std::sqrt(std::max(0.0, 1.0 - mu_i * mu_i));
  const Numeric sin_r     = std::sqrt(std::max(0.0, 1.0 - mu_r * mu_r));
  const Numeric tan_i     = sin_i / mu_i;
  const Numeric tan_r     = sin_r / mu_r;
  const Numeric cosine    = std::cos(relative_azimuth);
  const Numeric cos_alpha = unit_clamp(mu_i * mu_r - sin_i * sin_r * cosine);
  const Numeric distance  = std::sqrt(std::max(0.0, tan_i * tan_i + tan_r * tan_r + 2.0 * tan_i * tan_r * cosine));
  const Numeric phase =
      (1.0 - asymmetry * asymmetry) / std::pow(1.0 + asymmetry * asymmetry + 2.0 * asymmetry * cos_alpha, 1.5);
  return rho0 * std::pow(mu_i * mu_r * (mu_i + mu_r), kappa - 1.0) * phase * (1.0 + (1.0 - hotspot) / (1.0 + distance));
}

Numeric ross_li_brdf(const Numeric outgoing_mu,
                     const Numeric incoming_mu,
                     const Numeric relative_azimuth,
                     const Numeric isotropic,
                     const Numeric volumetric,
                     const Numeric geometric,
                     const Numeric hotspot_angle) {
  const Numeric mu_r = checked_surface_mu(outgoing_mu);
  const Numeric mu_i = checked_surface_mu(incoming_mu);
  if (not std::isfinite(isotropic) or not std::isfinite(volumetric) or not std::isfinite(geometric) or
      not std::isfinite(hotspot_angle) or hotspot_angle <= 0.0)
    throw std::domain_error("Invalid Ross-Li BRDF parameters");

  const Numeric sin_i          = std::sqrt(std::max(0.0, 1.0 - mu_i * mu_i));
  const Numeric sin_r          = std::sqrt(std::max(0.0, 1.0 - mu_r * mu_r));
  const Numeric tan_i          = sin_i / mu_i;
  const Numeric tan_r          = sin_r / mu_r;
  const Numeric cosine         = std::cos(relative_azimuth);
  const Numeric sine           = std::sin(relative_azimuth);
  const Numeric cos_alpha      = unit_clamp(mu_i * mu_r - sin_i * sin_r * cosine);
  const Numeric sin_alpha      = std::sqrt(std::max(0.0, 1.0 - cos_alpha * cos_alpha));
  const Numeric alpha          = std::acos(cos_alpha);
  const Numeric hotspot_factor = 1.0 + 1.0 / (1.0 + alpha / hotspot_angle);
  const Numeric volume         = 4.0 / (3.0 * Constant::pi * (mu_i + mu_r)) *
                                     ((Constant::pi / 2.0 - alpha) * cos_alpha + sin_alpha) * hotspot_factor -
                                 1.0 / 3.0;

  constexpr Numeric height_to_base = 2.0;
  const Numeric     distance_sq    = tan_i * tan_i + tan_r * tan_r + 2.0 * tan_i * tan_r * cosine;
  const Numeric     cos_t = height_to_base * mu_i * mu_r / (mu_i + mu_r) *
                            std::sqrt(std::max(0.0, distance_sq + tan_i * tan_i * tan_r * tan_r * sine * sine));
  const Numeric     t     = cos_t >= -1.0 and cos_t <= 1.0 ? std::acos(cos_t) : 0.0;
  const Numeric     geometric_kernel =
      (mu_i + mu_r) / (Constant::pi * mu_i * mu_r) * (t - std::sin(t) * std::cos(t) - Constant::pi) +
      (1.0 + cos_alpha) / (2.0 * mu_i * mu_r);
  return std::max(0.0, isotropic + geometric * geometric_kernel + volumetric * volume);
}

void check_surface_weight(const Numeric weight) {
  if (not std::isfinite(weight) or weight < 0.0)
    throw std::domain_error("Surface-mixture weights must be finite and nonnegative");
}

void check_surface_fraction(const Numeric fraction) {
  if (not std::isfinite(fraction) or fraction < 0.0 or fraction > 1.0)
    throw std::domain_error("Surface-mixture fraction must be finite and in [0, 1]");
}

fresnel_amplitudes dielectric_fresnel_amplitudes(const Numeric incident_mu, const Complex refractive_index) {
  const Numeric mu = checked_surface_mu(incident_mu);
  if (not std::isfinite(refractive_index.real()) or not std::isfinite(refractive_index.imag()) or
      refractive_index.real() <= 0.0)
    throw std::domain_error("The Fresnel refractive index must be finite with positive real part");
  const Complex transmitted_sine = std::sqrt(std::max(0.0, 1.0 - mu * mu)) / refractive_index;
  Complex       transmitted_mu   = std::sqrt(Complex{1.0, 0.0} - transmitted_sine * transmitted_sine);
  if (transmitted_mu.real() < 0.0 or (transmitted_mu.real() == 0.0 and transmitted_mu.imag() < 0.0))
    transmitted_mu = -transmitted_mu;
  return {
      .vertical   = (refractive_index * mu - transmitted_mu) / (refractive_index * mu + transmitted_mu),
      .horizontal = (mu - refractive_index * transmitted_mu) / (mu + refractive_index * transmitted_mu),
  };
}

cox_munk_optics cox_munk_reflection(const Numeric outgoing_mu,
                                    const Numeric incoming_mu,
                                    const Numeric relative_azimuth,
                                    const Numeric wind_speed,
                                    const Complex refractive_index,
                                    const bool    shadowing) {
  const Numeric mu_r = checked_surface_mu(outgoing_mu);
  const Numeric mu_i = checked_surface_mu(incoming_mu);
  if (not std::isfinite(wind_speed) or wind_speed < 0.0)
    throw std::domain_error("Cox-Munk wind speed must be finite and nonnegative");

  const Numeric sin_i       = std::sqrt(std::max(0.0, 1.0 - mu_i * mu_i));
  const Numeric sin_r       = std::sqrt(std::max(0.0, 1.0 - mu_r * mu_r));
  const Numeric cos_theta   = std::clamp(-mu_i * mu_r + sin_i * sin_r * std::cos(relative_azimuth), -1.0, 1.0);
  const Numeric denominator = 2.0 * (1.0 - cos_theta);
  if (denominator == 0.0) return {};
  const Numeric normal_mu_sq = (mu_i + mu_r) * (mu_i + mu_r) / denominator;
  if (normal_mu_sq <= 0.0) return {};

  const Numeric facet_mu       = std::sqrt(std::max(0.0, 0.5 * (1.0 - cos_theta)));
  const auto    amplitudes     = dielectric_fresnel_amplitudes(facet_mu, refractive_index);
  const Numeric slope_variance = 0.003 + 0.00512 * wind_speed;
  const Numeric probability =
      std::exp(-(1.0 - normal_mu_sq) / (slope_variance * normal_mu_sq)) / (Constant::pi * slope_variance);
  Numeric factor = probability / (4.0 * mu_i * mu_r * normal_mu_sq * normal_mu_sq);
  if (shadowing) factor /= 1.0 + shadow_eta(mu_i, slope_variance) + shadow_eta(mu_r, slope_variance);
  return {.factor = factor, .amplitudes = amplitudes};
}

Numeric phi1_neg(const Numeric x) {
  if (x == 0.0) return 1.0;
  return -std::expm1(-x) / x;
}

Numeric phi2(const Numeric x) {
  if (std::abs(x) >= 1.0) return (std::expm1(x) - x) / (x * x);

  Numeric term = 0.5;
  Numeric sum  = term;
  for (Index k = 1; k <= 16; ++k) {
    term *= x / static_cast<Numeric>(k + 2);
    sum  += term;
  }
  return sum;
}

void barycentric_weights(Vector& weights, const ConstVectorView& nodes) {
  const Index n = static_cast<Index>(nodes.size());
  weights.resize(n);
  Vector  log_weights(n);
  Numeric largest = -std::numeric_limits<Numeric>::infinity();
  for (Index i = 0; i < n; ++i) {
    Numeric sign   = 1.0;
    log_weights[i] = 0.0;
    for (Index j = 0; j < n; ++j) {
      if (i == j) continue;
      const Numeric difference  = nodes[i] - nodes[j];
      log_weights[i]           -= std::log(std::abs(difference));
      if (difference < 0.0) sign = -sign;
    }
    weights[i] = sign;
    largest    = std::max(largest, log_weights[i]);
  }
  for (Index i = 0; i < n; ++i) weights[i] *= std::exp(log_weights[i] - largest);
}

void initialize_streams(VectorView mu, VectorView inv_mu, VectorView weights) {
  ARTS_USER_ERROR_IF(mu.size() != 2 * weights.size() or inv_mu.size() != mu.size(),
                     "Double-Gauss arrays have sizes mu={}, inv_mu={}, weights={}; expected [2N, 2N, N]",
                     mu.size(),
                     inv_mu.size(),
                     weights.size());
  const Index n = static_cast<Index>(weights.size());
  Legendre::PositiveDoubleGaussLegendre(mu[Range{0, n}], weights);
  std::transform(mu.begin(), mu.begin() + n, mu.begin() + n, [](const Numeric x) { return -x; });
  std::transform(mu.begin(), mu.end(), inv_mu.begin(), [](const Numeric x) { return 1.0 / x; });
}

void check_layer_input_sizes(const Index n_layers, const AscendingGrid& tau, const ConstVectorView& omega) {
  ARTS_USER_ERROR_IF(n_layers <= 0, "A DISORT solver requires at least one layer");
  ARTS_USER_ERROR_IF(
      static_cast<Index>(tau.size()) != n_layers, "tau_arr has size {}, expected {}", tau.size(), n_layers);
  ARTS_USER_ERROR_IF(
      static_cast<Index>(omega.size()) != n_layers, "omega_arr has size {}, expected {}", omega.size(), n_layers);
}

void check_layer_input_values(const AscendingGrid&   tau,
                              const ConstVectorView& omega,
                              const Numeric          mu0,
                              const Numeric          phi0) {
  ARTS_USER_ERROR_IF(
      tau.empty() or tau.front() <= 0.0, "tau_arr must be nonempty and strictly positive, got {:B,}", tau);
  ARTS_USER_ERROR_IF(stdr::any_of(omega, [](const Numeric x) { return x < 0.0 or x > 1.0; }),
                     "omega_arr must be in [0, 1], got {:B,}",
                     omega);
  ARTS_USER_ERROR_IF(mu0 < 0.0 or mu0 > 1.0, "mu0 must be in [0, 1], got {}", mu0);
  ARTS_USER_ERROR_IF(phi0 < 0.0 or phi0 >= Constant::two_pi, "phi0 must be in [0, 2*pi), got {}", phi0);
}

Index layer_index(const AscendingGrid& tau, const Numeric value) {
  ARTS_USER_ERROR_IF(tau.empty(), "Cannot locate an optical depth in an empty layer grid");
  ARTS_USER_ERROR_IF(value < 0.0 or value > tau.back(), "tau ({}) must be in [0, {}]", value, tau.back());
  const Index layer = std::distance(tau.begin(), stdr::lower_bound(tau, value));
  return std::min(layer, static_cast<Index>(tau.size()) - 1);
}

Numeric user_angle_centered_first_moment(const Numeric center,
                                         const Numeric lower,
                                         const Numeric upper,
                                         const Numeric observation,
                                         const Numeric abs_mu,
                                         const bool    downward) {
  if (upper <= lower) return 0.0;
  const Numeric near        = downward ? upper : lower;
  const Numeric distance    = downward ? observation - near : near - observation;
  const Numeric z           = (upper - lower) / abs_mu;
  const Numeric attenuation = std::exp(-distance / abs_mu);
  const Numeric l0          = -std::expm1(-z);
  Numeric       l1;
  if (std::abs(z) < 1.0e-3) {
    const Numeric z2 = z * z;
    l1               = z2 * (0.5 + z * (-1.0 / 3.0 + z * (1.0 / 8.0 + z * (-1.0 / 30.0 + z / 144.0))));
  } else {
    l1 = 1.0 - (1.0 + z) * std::exp(-z);
  }
  const Numeric direction = downward ? -1.0 : 1.0;
  return attenuation * ((near - center) * l0 + direction * abs_mu * l1);
}

Numeric downward_tms_kernel(const Numeric mu,
                            const Numeric mu0,
                            const Numeric thickness,
                            const Numeric value_at_lower_end,
                            const Numeric value_at_upper_end) {
  const Numeric y = 1.0 / mu - 1.0 / mu0;
  const Numeric w = thickness * y;
  if (std::abs(w) < 1.0) return thickness / mu * phi1_neg(w) * value_at_lower_end;
  return (value_at_lower_end - value_at_upper_end) / (mu * y);
}

Numeric ims_chi(const Numeric tau, const Numeric mu, const Numeric scaled_mu0) {
  const Numeric x = 1.0 / mu - 1.0 / scaled_mu0;
  const Numeric z = -tau * x;
  if (std::abs(z) < 1.0) return tau * tau * std::exp(-tau / scaled_mu0) * phi2(z) / (mu * scaled_mu0);
  return ((tau - 1.0 / x) * std::exp(-tau / scaled_mu0) + std::exp(-tau / mu) / x) / (mu * scaled_mu0 * x);
}

std::array<Numeric, 2> centered_pair_integrals(const Numeric kappa,
                                               const Numeric center,
                                               const Numeric lower,
                                               const Numeric upper,
                                               const Numeric observation,
                                               const Numeric abs_mu,
                                               const bool    downward) {
  const Numeric max_s = std::max(std::abs(lower - center), std::abs(upper - center));
  if (std::abs(kappa) * max_s < centered_integral_limit) {
    return {user_angle_exponential_integral<Numeric>(0.0, center, lower, upper, observation, abs_mu, downward),
            user_angle_centered_first_moment(center, lower, upper, observation, abs_mu, downward)};
  }

  const Numeric plus =
      user_angle_exponential_integral<Numeric>(kappa, center, lower, upper, observation, abs_mu, downward);
  const Numeric minus =
      user_angle_exponential_integral<Numeric>(-kappa, center, lower, upper, observation, abs_mu, downward);
  return {0.5 * (plus + minus), (plus - minus) / (2.0 * kappa)};
}

}  // namespace disort_common

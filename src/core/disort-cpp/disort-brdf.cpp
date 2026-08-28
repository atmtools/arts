#include "disort-brdf.h"

#include <legendre.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace disort::brdf {
namespace {
Numeric checked_mu(const Numeric mu) {
  const Numeric value = std::abs(mu);
  if (not(value > 0.0 and value <= 1.0)) throw std::domain_error("BRDF direction cosine must satisfy 0 < |mu| <= 1");
  return value;
}

Numeric unit_clamp(const Numeric value) { return std::clamp(value, -1.0, 1.0); }

Numeric shadow_eta(const Numeric mu, const Numeric slope_variance) {
  const Numeric sine = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  if (sine == 0.0) return 0.0;
  const Numeric cotangent     = mu / sine;
  const Numeric root_variance = std::sqrt(slope_variance);
  const Numeric first =
      root_variance / (std::sqrt(Constant::pi) * cotangent) * std::exp(-cotangent * cotangent / slope_variance);
  return 0.5 * (first - std::erfc(cotangent / root_variance));
}
}  // namespace

Numeric Hapke::operator()(const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) const {
  const Numeric mu  = checked_mu(outgoing_mu);
  const Numeric mup = checked_mu(incoming_mu);
  if (not(single_scattering_albedo >= 0.0 and single_scattering_albedo <= 1.0) or opposition_width < 0.0)
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

Numeric CoxMunk::operator()(const Numeric outgoing_mu,
                            const Numeric incoming_mu,
                            const Numeric relative_azimuth) const {
  const Numeric mu_r = checked_mu(outgoing_mu);
  const Numeric mu_i = checked_mu(incoming_mu);
  if (wind_speed < 0.0 or refractive_index <= 0.0) throw std::domain_error("Invalid Cox-Munk BRDF parameters");

  const Numeric sin_i       = std::sqrt(std::max(0.0, 1.0 - mu_i * mu_i));
  const Numeric sin_r       = std::sqrt(std::max(0.0, 1.0 - mu_r * mu_r));
  const Numeric cos_theta   = unit_clamp(-mu_i * mu_r + sin_i * sin_r * std::cos(relative_azimuth));
  const Numeric denominator = 2.0 * (1.0 - cos_theta);
  if (denominator == 0.0) return 0.0;
  const Numeric mu_n_sq = (mu_i + mu_r) * (mu_i + mu_r) / denominator;
  if (mu_n_sq <= 0.0) return 0.0;

  const Numeric slope_variance = 0.003 + 0.00512 * wind_speed;
  const Numeric probability = std::exp(-(1.0 - mu_n_sq) / (slope_variance * mu_n_sq)) / (Constant::pi * slope_variance);
  const Numeric sin_li      = std::sqrt(std::max(0.0, 1.0 - 0.5 * (1.0 - cos_theta)));
  const Numeric cos_li      = std::sqrt(std::max(0.0, 0.5 * (1.0 - cos_theta)));
  const Numeric sin_lt      = sin_li / refractive_index;
  Numeric       fresnel     = 1.0;
  if (sin_lt <= 1.0) {
    const Numeric cos_lt = std::sqrt(std::max(0.0, 1.0 - sin_lt * sin_lt));
    const Numeric rs     = (cos_li - refractive_index * cos_lt) / (cos_li + refractive_index * cos_lt);
    const Numeric rp     = (refractive_index * cos_li - cos_lt) / (cos_lt + refractive_index * cos_li);
    fresnel              = 0.5 * (rs * rs + rp * rp);
  }
  Numeric result = probability * fresnel / (4.0 * mu_i * mu_r * mu_n_sq * mu_n_sq);
  if (shadowing) result /= shadow_eta(mu_i, slope_variance) + shadow_eta(mu_r, slope_variance) + 1.0;
  return result;
}

Numeric RPV::operator()(const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) const {
  const Numeric mu_r = checked_mu(outgoing_mu);
  const Numeric mu_i = checked_mu(incoming_mu);
  if (rho0 < 0.0 or not(asymmetry > -1.0 and asymmetry < 1.0) or hotspot < 0.0)
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

Numeric RossLi::operator()(const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) const {
  const Numeric mu_r = checked_mu(outgoing_mu);
  const Numeric mu_i = checked_mu(incoming_mu);
  if (hotspot_angle <= 0.0) throw std::domain_error("Ross-Li hotspot angle must be positive");

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
  const Numeric     geometric =
      (mu_i + mu_r) / (Constant::pi * mu_i * mu_r) * (t - std::sin(t) * std::cos(t) - Constant::pi) +
      (1.0 + cos_alpha) / (2.0 * mu_i * mu_r);
  return std::max(0.0, isotropic + geometric * this->geometric + volumetric * volume);
}

std::vector<BDRF> fourier_modes(RawFunction raw, const Index number_of_modes, const Index azimuth_quadrature_points) {
  if (not raw) throw std::invalid_argument("A raw BRDF function is required");
  if (number_of_modes < 0) throw std::invalid_argument("The number of BRDF Fourier modes cannot be negative");
  if (azimuth_quadrature_points <= 0)
    throw std::invalid_argument("The BRDF azimuth quadrature must contain at least one point");

  struct Quadrature {
    Vector nodes;
    Vector weights;
  };
  auto quadrature =
      std::make_shared<Quadrature>(Quadrature{Vector(azimuth_quadrature_points), Vector(azimuth_quadrature_points)});
  Legendre::PositiveDoubleGaussLegendre(quadrature->nodes, quadrature->weights);
  std::vector<BDRF> result;
  result.reserve(static_cast<std::size_t>(number_of_modes));
  for (Index mode = 0; mode < number_of_modes; ++mode) {
    result.emplace_back(BDRF{
        [raw, quadrature, mode](MatrixView output, const ConstVectorView& outgoing, const ConstVectorView& incoming) {
          const Numeric factor = mode == 0 ? 1.0 : 2.0;
          for (Size i = 0; i < outgoing.size(); ++i)
            for (Size j = 0; j < incoming.size(); ++j) {
              Numeric coefficient = 0.0;
              for (Size k = 0; k < quadrature->nodes.size(); ++k) {
                const Numeric azimuth  = Constant::pi * quadrature->nodes[k];
                coefficient           += quadrature->weights[k] * raw(outgoing[i], std::abs(incoming[j]), azimuth) *
                                         std::cos(static_cast<Numeric>(mode) * azimuth);
              }
              output[i, j] = factor * coefficient;
            }
        }});
  }
  return result;
}

}  // namespace disort::brdf

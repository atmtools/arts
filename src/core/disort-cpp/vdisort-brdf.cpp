#include "vdisort-brdf.h"

#include <arts_constants.h>
#include <legendre.h>
#include <rtepack_surface.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <memory>
#include <stdexcept>

namespace vdisort::brdf {
namespace {
Numeric checked_mu(const Numeric mu) {
  const Numeric value = std::abs(mu);
  if (not(value > 0.0 and value <= 1.0)) throw std::domain_error("BPrDF direction cosine must satisfy 0 < |mu| <= 1");
  return value;
}

Numeric shadow_eta(const Numeric mu, const Numeric slope_variance) {
  const Numeric sine = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  if (sine == 0.0) return 0.0;
  const Numeric cotangent = mu / sine;
  const Numeric root      = std::sqrt(slope_variance);
  return 0.5 * (root / (std::sqrt(Constant::pi) * cotangent) * std::exp(-cotangent * cotangent / slope_variance) -
                std::erfc(cotangent / root));
}

Vector3 direction(const Numeric mu, const Numeric azimuth) {
  const Numeric sine = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  return {sine * std::cos(azimuth), sine * std::sin(azimuth), mu};
}

void fill_combined(rtepack::muelmat&       cosine,
                   rtepack::muelmat&       sine,
                   const rtepack::muelmat& ordinary_cosine,
                   const rtepack::muelmat& ordinary_sine,
                   const Index             mode) {
  cosine = rtepack::muelmat{0.0};
  sine   = rtepack::muelmat{0.0};
  if (mode == 0) {
    for (Index i = 0; i < 2; ++i)
      for (Index j = 0; j < 2; ++j) cosine[i, j] = ordinary_cosine[i, j];
    for (Index i = 2; i < 4; ++i)
      for (Index j = 2; j < 4; ++j) sine[i, j] = ordinary_cosine[i, j];
    return;
  }
  for (Index i = 0; i < 4; ++i)
    for (Index j = 0; j < 4; ++j) {
      const bool same_block = (i < 2) == (j < 2);
      if (same_block) {
        cosine[i, j] = ordinary_cosine[i, j];
        sine[i, j]   = ordinary_cosine[i, j];
      } else if (i < 2) {
        cosine[i, j] = -ordinary_sine[i, j];
        sine[i, j]   = ordinary_sine[i, j];
      } else {
        cosine[i, j] = ordinary_sine[i, j];
        sine[i, j]   = -ordinary_sine[i, j];
      }
    }
}
}  // namespace

rtepack::muelmat CoxMunk::operator()(const Numeric outgoing_mu,
                                     const Numeric incoming_mu,
                                     const Numeric relative_azimuth) const {
  const Numeric mu_r = checked_mu(outgoing_mu);
  const Numeric mu_i = checked_mu(incoming_mu);
  if (wind_speed < 0.0 or refractive_index <= 0.0) throw std::domain_error("Invalid polarized Cox-Munk parameters");

  const Vector3 incident = direction(-mu_i, 0.0);
  const Vector3 outgoing = direction(mu_r, relative_azimuth);
  Vector3       normal   = outgoing - incident;
  const Numeric norm     = hypot(normal);
  if (norm == 0.0) return rtepack::muelmat{0.0};
  normal /= norm;
  if (normal[2] <= 0.0) return rtepack::muelmat{0.0};

  const Numeric facet_mu = -dot(incident, normal);
  if (facet_mu <= 0.0) return rtepack::muelmat{0.0};
  const Numeric transmitted_sine = std::sqrt(std::max(0.0, 1.0 - facet_mu * facet_mu)) / refractive_index;
  if (transmitted_sine > 1.0) return rtepack::muelmat{0.0};
  const Numeric transmitted_mu = std::sqrt(std::max(0.0, 1.0 - transmitted_sine * transmitted_sine));
  const Numeric rv = (refractive_index * facet_mu - transmitted_mu) / (refractive_index * facet_mu + transmitted_mu);
  const Numeric rh = (facet_mu - refractive_index * transmitted_mu) / (facet_mu + refractive_index * transmitted_mu);

  const Numeric slope_variance = 0.003 + 0.00512 * wind_speed;
  const Numeric normal_mu_sq   = normal[2] * normal[2];
  const Numeric probability =
      std::exp(-(1.0 - normal_mu_sq) / (slope_variance * normal_mu_sq)) / (Constant::pi * slope_variance);
  Numeric factor = probability / (4.0 * mu_i * mu_r * normal_mu_sq * normal_mu_sq);
  if (shadowing) factor /= 1.0 + shadow_eta(mu_i, slope_variance) + shadow_eta(mu_r, slope_variance);
  return factor * rtepack::fresnel_reflectance_nonspecular(rv, rh, incident, outgoing, normal);
}

rtepack::muelmat Fresnel::operator()(const Numeric incident_mu) const {
  const Numeric mu = checked_mu(incident_mu);
  if (refractive_index <= 0.0) throw std::domain_error("The Fresnel refractive index must be positive");
  const Numeric transmitted_sine = std::sqrt(std::max(0.0, 1.0 - mu * mu)) / refractive_index;
  const Complex transmitted_mu   = std::sqrt(Complex{1.0 - transmitted_sine * transmitted_sine, 0.0});
  const Complex rv               = (refractive_index * mu - transmitted_mu) / (refractive_index * mu + transmitted_mu);
  const Complex rh               = (mu - refractive_index * transmitted_mu) / (mu + refractive_index * transmitted_mu);
  return rtepack::fresnel_reflectance(rv, rh);
}

std::vector<BDRF> fourier_modes(RawFunction raw, const Index number_of_modes, const Index azimuth_quadrature_points) {
  if (not raw) throw std::invalid_argument("A raw polarized BPrDF is required");
  if (number_of_modes < 0) throw std::invalid_argument("The number of BPrDF Fourier modes cannot be negative");
  if (azimuth_quadrature_points <= 0)
    throw std::invalid_argument("The BPrDF azimuth quadrature must contain at least one point");

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
    const auto evaluate = [raw, quadrature, mode](const Index                  alpha,
                                                  rtepack::muelmat_matrix_view output,
                                                  const ConstVectorView&       outgoing,
                                                  const ConstVectorView&       incoming) {
      const Numeric factor = 2.0 * (mode == 0 ? 1.0 : 2.0);
      for (Index i = 0; i < static_cast<Index>(outgoing.size()); ++i)
        for (Index j = 0; j < static_cast<Index>(incoming.size()); ++j) {
          rtepack::muelmat ordinary_cosine{0.0}, ordinary_sine{0.0};
          for (Index k = 0; k < static_cast<Index>(quadrature->nodes.size()); ++k) {
            const Numeric azimuth  = Constant::two_pi * quadrature->nodes[k];
            const auto    matrix   = raw(outgoing[i], std::abs(incoming[j]), azimuth);
            ordinary_cosine       += quadrature->weights[k] * std::cos(static_cast<Numeric>(mode) * azimuth) * matrix;
            ordinary_sine         += quadrature->weights[k] * std::sin(static_cast<Numeric>(mode) * azimuth) * matrix;
          }
          ordinary_cosine *= factor;
          ordinary_sine   *= factor;
          rtepack::muelmat combined_cosine, combined_sine;
          fill_combined(combined_cosine, combined_sine, ordinary_cosine, ordinary_sine, mode);
          output[i, j] = alpha == cosine_mode ? combined_cosine : combined_sine;
        }
    };
    const auto evaluate_beam = [raw, quadrature, mode](const Index                  alpha,
                                                       rtepack::muelmat_matrix_view output,
                                                       const ConstVectorView&       outgoing,
                                                       const ConstVectorView&       incoming) {
      const Numeric factor = 2.0 * (mode == 0 ? 1.0 : 2.0);
      for (Index i = 0; i < static_cast<Index>(outgoing.size()); ++i)
        for (Index j = 0; j < static_cast<Index>(incoming.size()); ++j) {
          rtepack::muelmat ordinary_cosine{0.0}, ordinary_sine{0.0};
          for (Index k = 0; k < static_cast<Index>(quadrature->nodes.size()); ++k) {
            const Numeric azimuth  = Constant::two_pi * quadrature->nodes[k];
            const auto    matrix   = raw(outgoing[i], std::abs(incoming[j]), azimuth);
            ordinary_cosine       += quadrature->weights[k] * std::cos(static_cast<Numeric>(mode) * azimuth) * matrix;
            ordinary_sine         += quadrature->weights[k] * std::sin(static_cast<Numeric>(mode) * azimuth) * matrix;
          }
          ordinary_cosine *= factor;
          // VDISORT reconstructs modes with phi_in - phi_out, whereas the
          // raw BPrDF is parameterized by phi_out - phi_in.
          ordinary_sine *= -factor;
          rtepack::muelmat combined{0.0};
          for (Index row = 0; row < 4; ++row)
            for (Index column = 0; column < 4; ++column) {
              if (alpha == cosine_mode)
                combined[row, column] = row < 2 ? ordinary_cosine[row, column] : ordinary_sine[row, column];
              else
                combined[row, column] = row < 2 ? ordinary_sine[row, column] : ordinary_cosine[row, column];
            }
          output[i, j] = combined;
        }
    };
    result.push_back(BDRF{
        .cosine = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view out,
                                          const ConstVectorView&       mu_out,
                                          const ConstVectorView& mu_in) { evaluate(cosine_mode, out, mu_out, mu_in); }},
        .sine   = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view out,
                                          const ConstVectorView&       mu_out,
                                          const ConstVectorView& mu_in) { evaluate(sine_mode, out, mu_out, mu_in); }},
        .beam_cosine = BDRF::func_t{[evaluate_beam](rtepack::muelmat_matrix_view out,
                                                    const ConstVectorView&       mu_out,
                                                    const ConstVectorView&       mu_in) {
          evaluate_beam(cosine_mode, out, mu_out, mu_in);
        }},
        .beam_sine   = BDRF::func_t{
            [evaluate_beam](rtepack::muelmat_matrix_view out,
                            const ConstVectorView&       mu_out,
                            const ConstVectorView&       mu_in) { evaluate_beam(sine_mode, out, mu_out, mu_in); }}});
  }
  return result;
}

std::vector<BDRF> cox_munk_fourier_modes(const Numeric wind_speed,
                                         const Numeric refractive_index,
                                         const bool    shadowing,
                                         const Index   number_of_modes,
                                         const Index   azimuth_quadrature_points) {
  return fourier_modes(CoxMunk{wind_speed, refractive_index, shadowing}, number_of_modes, azimuth_quadrature_points);
}

std::vector<BDRF> fresnel_fourier_modes(const Numeric refractive_index, const Index number_of_modes) {
  if (number_of_modes < 0) throw std::invalid_argument("The number of Fresnel Fourier modes cannot be negative");
  const Fresnel fresnel{refractive_index};
  // Validate eagerly, including when zero modes are requested.
  static_cast<void>(fresnel(1.0));

  std::vector<BDRF> result;
  result.reserve(static_cast<std::size_t>(number_of_modes));
  for (Index mode = 0; mode < number_of_modes; ++mode) {
    const auto evaluate = [fresnel, mode](const Index                  alpha,
                                          rtepack::muelmat_matrix_view output,
                                          const ConstVectorView&       outgoing,
                                          const ConstVectorView&       incoming) {
      output = rtepack::muelmat{0.0};
      if (alpha == sine_mode) return;

      Vector nodes(incoming.size()), weights(incoming.size());
      Legendre::PositiveDoubleGaussLegendre(nodes, weights);
      const Numeric azimuth_factor = mode == 0 ? 1.0 : 2.0;
      for (Index i = 0; i < static_cast<Index>(outgoing.size()); ++i)
        for (Index j = 0; j < static_cast<Index>(incoming.size()); ++j) {
          const Numeric mu_in = std::abs(incoming[j]);
          if (std::abs(outgoing[i] - mu_in) > 64.0 * std::numeric_limits<Numeric>::epsilon()) continue;
          output[i, j] = azimuth_factor / (Constant::pi * weights[j] * mu_in) * fresnel(mu_in);
        }
    };
    const auto unsupported_beam = [](rtepack::muelmat_matrix_view, const ConstVectorView&, const ConstVectorView&) {
      throw std::runtime_error(
          "Ideal Fresnel reflection is a directional delta distribution and cannot reflect a direct beam into the "
          "native VDISORT quadrature; use a finite-width surface model such as Cox-Munk");
    };
    result.push_back(BDRF{
        .cosine = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view out,
                                          const ConstVectorView&       mu_out,
                                          const ConstVectorView& mu_in) { evaluate(cosine_mode, out, mu_out, mu_in); }},
        .sine   = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view out,
                                          const ConstVectorView&       mu_out,
                                          const ConstVectorView& mu_in) { evaluate(sine_mode, out, mu_out, mu_in); }},
        .beam_cosine = BDRF::func_t{unsupported_beam},
        .beam_sine   = BDRF::func_t{unsupported_beam}});
  }
  return result;
}

}  // namespace vdisort::brdf

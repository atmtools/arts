#include "vdisort-brdf.h"

#include <arts_constants.h>
#include <common.h>
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
  const auto reflection = disort_common::cox_munk_reflection(
      outgoing_mu, incoming_mu, relative_azimuth, wind_speed, refractive_index, shadowing);
  if (reflection.factor == 0.0) return rtepack::muelmat{0.0};

  const Numeric mu_r     = std::abs(outgoing_mu);
  const Numeric mu_i     = std::abs(incoming_mu);
  const Vector3 incident = direction(-mu_i, 0.0);
  const Vector3 outgoing = direction(mu_r, relative_azimuth);
  Vector3       normal   = outgoing - incident;
  const Numeric norm     = hypot(normal);
  if (norm == 0.0) return rtepack::muelmat{0.0};
  normal /= norm;
  return reflection.factor *
         rtepack::fresnel_reflectance_nonspecular(
             reflection.amplitudes.vertical, reflection.amplitudes.horizontal, incident, outgoing, normal);
}

rtepack::muelmat Fresnel::operator()(const Numeric incident_mu) const {
  const auto amplitudes = disort_common::dielectric_fresnel_amplitudes(incident_mu, refractive_index);
  return rtepack::fresnel_reflectance(amplitudes.vertical, amplitudes.horizontal);
}

std::vector<BDRF> lambertian_fourier_modes(const Numeric albedo, const Index number_of_modes) {
  const Numeric coefficient = 2.0 * disort_common::lambertian_brdf(albedo);
  if (number_of_modes < 0) throw std::invalid_argument("The number of Lambertian Fourier modes cannot be negative");

  std::vector<BDRF> result;
  result.reserve(static_cast<std::size_t>(number_of_modes));
  for (Index mode = 0; mode < number_of_modes; ++mode) {
    const auto evaluate = [coefficient = mode == 0 ? coefficient : 0.0](
                              rtepack::muelmat_matrix_view output, const ConstVectorView&, const ConstVectorView&) {
      output = rtepack::muelmat{0.0};
      if (coefficient != 0.0)
        for (Index i = 0; i < output.nrows(); ++i)
          for (Index j = 0; j < output.ncols(); ++j) output[i, j][0, 0] = coefficient;
    };
    const auto zero = [](rtepack::muelmat_matrix_view output, const ConstVectorView&, const ConstVectorView&) {
      output = rtepack::muelmat{0.0};
    };
    result.push_back(BDRF{.cosine      = BDRF::func_t{evaluate},
                          .sine        = BDRF::func_t{zero},
                          .beam_cosine = BDRF::func_t{evaluate},
                          .beam_sine   = BDRF::func_t{zero}});
  }
  return result;
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

std::vector<BDRF> combine_fourier_modes(std::vector<BDRF> first,
                                        const Numeric     first_weight,
                                        std::vector<BDRF> second,
                                        const Numeric     second_weight) {
  disort_common::check_surface_weight(first_weight);
  disort_common::check_surface_weight(second_weight);

  struct Components {
    std::vector<BDRF> first;
    Numeric           first_weight;
    std::vector<BDRF> second;
    Numeric           second_weight;
  };
  auto components =
      std::make_shared<Components>(Components{std::move(first), first_weight, std::move(second), second_weight});
  const auto number_of_modes = std::max(components->first.size(), components->second.size());

  std::vector<BDRF> result;
  result.reserve(number_of_modes);
  for (std::size_t mode = 0; mode < number_of_modes; ++mode) {
    const auto evaluate = [components, mode](const bool                   beam,
                                             const Index                  alpha,
                                             rtepack::muelmat_matrix_view output,
                                             const ConstVectorView&       outgoing,
                                             const ConstVectorView&       incoming) {
      output = rtepack::muelmat{0.0};
      rtepack::muelmat_matrix scratch(output.nrows(), output.ncols(), rtepack::muelmat{0.0});
      const auto              add = [&](const std::vector<BDRF>& modes, const Numeric weight) {
        if (weight == 0.0 or mode >= modes.size()) return;
        if (beam)
          modes[mode].beam(alpha, scratch, outgoing, incoming);
        else
          modes[mode](alpha, scratch, outgoing, incoming);
        for (Index i = 0; i < output.nrows(); ++i)
          for (Index j = 0; j < output.ncols(); ++j) output[i, j] += weight * scratch[i, j];
      };
      add(components->first, components->first_weight);
      add(components->second, components->second_weight);
    };
    result.push_back(BDRF{.cosine      = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view output,
                                                                 const ConstVectorView&       outgoing,
                                                                 const ConstVectorView&       incoming) {
                            evaluate(false, cosine_mode, output, outgoing, incoming);
                          }},
                          .sine        = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view output,
                                                                 const ConstVectorView&       outgoing,
                                                                 const ConstVectorView&       incoming) {
                            evaluate(false, sine_mode, output, outgoing, incoming);
                          }},
                          .beam_cosine = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view output,
                                                                 const ConstVectorView&       outgoing,
                                                                 const ConstVectorView&       incoming) {
                            evaluate(true, cosine_mode, output, outgoing, incoming);
                          }},
                          .beam_sine   = BDRF::func_t{[evaluate](rtepack::muelmat_matrix_view output,
                                                                 const ConstVectorView&       outgoing,
                                                                 const ConstVectorView&       incoming) {
                            evaluate(true, sine_mode, output, outgoing, incoming);
                          }}});
  }
  return result;
}

std::vector<BDRF> fresnel_lambertian_fourier_modes(const Numeric fresnel_fraction,
                                                   const Numeric lambertian_albedo,
                                                   const Numeric refractive_index,
                                                   const Index   number_of_modes) {
  disort_common::check_surface_fraction(fresnel_fraction);
  return combine_fourier_modes(fresnel_fourier_modes(refractive_index, number_of_modes),
                               fresnel_fraction,
                               lambertian_fourier_modes(lambertian_albedo, number_of_modes),
                               1.0 - fresnel_fraction);
}

std::vector<BDRF> cox_munk_lambertian_fourier_modes(const Numeric cox_munk_fraction,
                                                    const Numeric lambertian_albedo,
                                                    const Numeric wind_speed,
                                                    const Numeric refractive_index,
                                                    const bool    shadowing,
                                                    const Index   number_of_modes,
                                                    const Index   azimuth_quadrature_points) {
  disort_common::check_surface_fraction(cox_munk_fraction);
  return combine_fourier_modes(
      cox_munk_fourier_modes(wind_speed, refractive_index, shadowing, number_of_modes, azimuth_quadrature_points),
      cox_munk_fraction,
      lambertian_fourier_modes(lambertian_albedo, number_of_modes),
      1.0 - cox_munk_fraction);
}

}  // namespace vdisort::brdf

#include "disort-brdf.h"

#include <common.h>
#include <legendre.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace disort::brdf {

Numeric Hapke::operator()(const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) const {
  return disort_common::hapke_brdf(
      outgoing_mu, incoming_mu, relative_azimuth, opposition_amplitude, opposition_width, single_scattering_albedo);
}

Numeric CoxMunk::operator()(const Numeric outgoing_mu,
                            const Numeric incoming_mu,
                            const Numeric relative_azimuth) const {
  const auto reflection = disort_common::cox_munk_reflection(
      outgoing_mu, incoming_mu, relative_azimuth, wind_speed, refractive_index, shadowing);
  return 0.5 * reflection.factor *
         (std::norm(reflection.amplitudes.vertical) + std::norm(reflection.amplitudes.horizontal));
}

Numeric RPV::operator()(const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) const {
  return disort_common::rpv_brdf(outgoing_mu, incoming_mu, relative_azimuth, rho0, kappa, asymmetry, hotspot);
}

Numeric RossLi::operator()(const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) const {
  return disort_common::ross_li_brdf(
      outgoing_mu, incoming_mu, relative_azimuth, isotropic, volumetric, geometric, hotspot_angle);
}

std::vector<BDRF> lambertian_fourier_modes(const Numeric albedo, const Index number_of_modes) {
  const Numeric coefficient = Constant::pi * disort_common::lambertian_brdf(albedo);
  if (number_of_modes < 0) throw std::invalid_argument("The number of Lambertian Fourier modes cannot be negative");

  std::vector<BDRF> result;
  result.reserve(static_cast<std::size_t>(number_of_modes));
  for (Index mode = 0; mode < number_of_modes; ++mode) {
    result.emplace_back(
        BDRF{[coefficient = mode == 0 ? coefficient : 0.0](
                 MatrixView output, const ConstVectorView&, const ConstVectorView&) { output = coefficient; }});
  }
  return result;
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
          const Numeric factor = Constant::pi * (mode == 0 ? 1.0 : 2.0);
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
    result.emplace_back(
        BDRF{[components, mode](MatrixView output, const ConstVectorView& outgoing, const ConstVectorView& incoming) {
          output = 0.0;
          Matrix     scratch(output.nrows(), output.ncols());
          const auto add = [&](const std::vector<BDRF>& modes, const Numeric weight) {
            if (weight == 0.0 or mode >= modes.size()) return;
            modes[mode](scratch, outgoing, incoming);
            for (Index i = 0; i < output.nrows(); ++i)
              for (Index j = 0; j < output.ncols(); ++j) output[i, j] += weight * scratch[i, j];
          };
          add(components->first, components->first_weight);
          add(components->second, components->second_weight);
        }});
  }
  return result;
}

std::vector<BDRF> cox_munk_lambertian_fourier_modes(const Numeric cox_munk_fraction,
                                                    const Numeric lambertian_albedo,
                                                    const Numeric wind_speed,
                                                    const Complex refractive_index,
                                                    const bool    shadowing,
                                                    const Index   number_of_modes,
                                                    const Index   azimuth_quadrature_points) {
  disort_common::check_surface_fraction(cox_munk_fraction);
  return combine_fourier_modes(
      fourier_modes(CoxMunk{wind_speed, refractive_index, shadowing}, number_of_modes, azimuth_quadrature_points),
      cox_munk_fraction,
      lambertian_fourier_modes(lambertian_albedo, number_of_modes),
      1.0 - cox_munk_fraction);
}

}  // namespace disort::brdf

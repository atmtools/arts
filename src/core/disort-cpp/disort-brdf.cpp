#include "disort-brdf.h"

#include <common.h>
#include <legendre.h>

#include <algorithm>
#include <cmath>
#include <format>
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

namespace {
//! Zero reflection of number_of_modes modes into mu_out from mu_in, for the beam mu0 (NaN: none)
surface_reflection zero_reflection(const Index            number_of_modes,
                                   const ConstVectorView& mu_out,
                                   const ConstVectorView& mu_in,
                                   const Numeric          mu0) {
  if (number_of_modes < 0) throw std::invalid_argument("The number of BRDF Fourier modes cannot be negative");
  return {number_of_modes, static_cast<Index>(mu_out.size()), static_cast<Index>(mu_in.size()), mu0};
}
}  // namespace

surface_reflection lambertian_fourier_modes(const Numeric          albedo,
                                            const Index            number_of_modes,
                                            const ConstVectorView& mu_out,
                                            const ConstVectorView& mu_in,
                                            const Numeric          mu0) {
  const Numeric coefficient = Constant::pi * disort_common::lambertian_brdf(albedo);
  auto          out         = zero_reflection(number_of_modes, mu_out, mu_in, mu0);
  if (number_of_modes > 0) {
    out.brdf()[0] = coefficient;
    if (not std::isnan(mu0)) out.beam_brdf()[0] = coefficient;
  }
  return out;
}

surface_reflection fourier_modes(const RawFunction&     raw,
                                 const Index            number_of_modes,
                                 const ConstVectorView& mu_out,
                                 const ConstVectorView& mu_in,
                                 const Numeric          mu0,
                                 const Index            azimuth_quadrature_points) {
  if (not raw) throw std::invalid_argument("A raw BRDF function is required");
  if (azimuth_quadrature_points <= 0)
    throw std::invalid_argument("The BRDF azimuth quadrature must contain at least one point");

  auto   out = zero_reflection(number_of_modes, mu_out, mu_in, mu0);
  Vector nodes(azimuth_quadrature_points), weights(azimuth_quadrature_points);
  Legendre::PositiveDoubleGaussLegendre(nodes, weights);

  // cos(m pi x_k) at the nodes x_k on [0, 1]
  Matrix cosine(number_of_modes, azimuth_quadrature_points);
  for (Index m = 0; m < number_of_modes; ++m)
    for (Index k = 0; k < azimuth_quadrature_points; ++k)
      cosine[m, k] = std::cos(static_cast<Numeric>(m) * (Constant::pi * nodes[k]));

  // R_m(o, i) = pi (2 - delta_m0) sum_k w_k rho(o, i, pi x_k) cos(m pi x_k), rho sampled once for all modes
  Vector     rho(azimuth_quadrature_points);
  const auto modes = [&](auto&& coefficients, const Numeric outgoing, const Numeric incoming) {
    for (Index k = 0; k < azimuth_quadrature_points; ++k)
      rho[k] = raw(outgoing, std::abs(incoming), Constant::pi * nodes[k]);
    for (Index m = 0; m < number_of_modes; ++m) {
      Numeric coefficient = 0.0;
      for (Index k = 0; k < azimuth_quadrature_points; ++k) coefficient += weights[k] * rho[k] * cosine[m, k];
      coefficients[m] = Constant::pi * (m == 0 ? 1.0 : 2.0) * coefficient;
    }
  };
  for (Size i = 0; i < mu_out.size(); ++i) {
    for (Size j = 0; j < mu_in.size(); ++j) modes(out.brdf()[joker, i, j], mu_out[i], mu_in[j]);
    if (not std::isnan(mu0)) modes(out.beam_brdf()[joker, i], mu_out[i], mu0);
  }
  return out;
}

surface_reflection combine_fourier_modes(const surface_reflection& first,
                                         const Numeric             first_weight,
                                         const surface_reflection& second,
                                         const Numeric             second_weight) {
  disort_common::check_surface_weight(first_weight);
  disort_common::check_surface_weight(second_weight);
  if (first.directions() != second.directions() or first.streams() != second.streams())
    throw std::invalid_argument(std::format(
        "Surface reflections can only be combined on the same directions, got {:B,} and {:B,} (modes, out, in)",
        first.brdf().shape(),
        second.brdf().shape()));
  if (not(first.beam_mu() == second.beam_mu() or (std::isnan(first.beam_mu()) and std::isnan(second.beam_mu()))))
    throw std::invalid_argument(std::format(
        "Surface reflections can only be combined for the same beam, got the beam cosines {} and {} (NaN: none)",
        first.beam_mu(),
        second.beam_mu()));

  const Index        nmodes = std::max(first.modes(), second.modes());
  surface_reflection out(nmodes, first.directions(), first.streams(), first.beam_mu());
  const auto         add = [&out](const surface_reflection& x, const Numeric weight) {
    if (weight == 0.0) return;
    for (Index m = 0; m < x.modes(); ++m) {
      for (Index i = 0; i < x.directions(); ++i) {
        for (Index j = 0; j < x.streams(); ++j) out.brdf()[m, i, j] += weight * x.brdf()[m, i, j];
        out.beam_brdf()[m, i] += weight * x.beam_brdf()[m, i];
      }
    }
  };
  add(first, first_weight);
  add(second, second_weight);
  return out;
}

surface_reflection cox_munk_lambertian_fourier_modes(const Numeric          cox_munk_fraction,
                                                     const Numeric          lambertian_albedo,
                                                     const Numeric          wind_speed,
                                                     const Complex          refractive_index,
                                                     const bool             shadowing,
                                                     const Index            number_of_modes,
                                                     const ConstVectorView& mu_out,
                                                     const ConstVectorView& mu_in,
                                                     const Numeric          mu0,
                                                     const Index            azimuth_quadrature_points) {
  disort_common::check_surface_fraction(cox_munk_fraction);
  return combine_fourier_modes(fourier_modes(CoxMunk{wind_speed, refractive_index, shadowing},
                                             number_of_modes,
                                             mu_out,
                                             mu_in,
                                             mu0,
                                             azimuth_quadrature_points),
                               cox_munk_fraction,
                               lambertian_fourier_modes(lambertian_albedo, number_of_modes, mu_out, mu_in, mu0),
                               1.0 - cox_munk_fraction);
}

}  // namespace disort::brdf

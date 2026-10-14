#include "vdisort-brdf.h"

#include <arts_constants.h>
#include <common.h>
#include <legendre.h>
#include <rtepack_surface.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <format>
#include <limits>
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

namespace {
//! Zero reflection of number_of_modes modes into mu_out from mu_in, for the beam mu0 (NaN: none), with a (zero)
//! specular part if specular
surface_reflection zero_reflection(const Index            number_of_modes,
                                   const ConstVectorView& mu_out,
                                   const ConstVectorView& mu_in,
                                   const Numeric          mu0,
                                   const bool             specular = false) {
  if (number_of_modes < 0) throw std::invalid_argument("The number of BRDF Fourier modes cannot be negative");
  return {number_of_modes, static_cast<Index>(mu_out.size()), static_cast<Index>(mu_in.size()), mu0, specular};
}
}  // namespace

surface_reflection lambertian_fourier_modes(const Numeric          albedo,
                                            const Index            number_of_modes,
                                            const ConstVectorView& mu_out,
                                            const ConstVectorView& mu_in,
                                            const Numeric          mu0) {
  const Numeric coefficient = 2.0 * disort_common::lambertian_brdf(albedo);
  auto          out         = zero_reflection(number_of_modes, mu_out, mu_in, mu0);
  if (number_of_modes > 0 and coefficient != 0.0) {
    for (Size i = 0; i < mu_out.size(); ++i) {
      for (Size j = 0; j < mu_in.size(); ++j) out.brdf()[cosine_mode, 0, i, j][0, 0] = coefficient;
      if (not std::isnan(mu0)) out.beam_brdf()[cosine_mode, 0, i][0, 0] = coefficient;
    }
  }
  return out;
}

surface_reflection fourier_modes(const RawFunction&     raw,
                                 const Index            number_of_modes,
                                 const ConstVectorView& mu_out,
                                 const ConstVectorView& mu_in,
                                 const Numeric          mu0,
                                 const Index            azimuth_quadrature_points) {
  if (not raw) throw std::invalid_argument("A raw polarized BPrDF is required");
  if (azimuth_quadrature_points <= 0)
    throw std::invalid_argument("The BPrDF azimuth quadrature must contain at least one point");

  auto   out = zero_reflection(number_of_modes, mu_out, mu_in, mu0);
  Vector nodes(azimuth_quadrature_points), weights(azimuth_quadrature_points);
  Legendre::PositiveDoubleGaussLegendre(nodes, weights);

  // cos(m phi_k) and sin(m phi_k) at the nodes phi_k on [0, 2 pi]
  Matrix cosine(number_of_modes, azimuth_quadrature_points), sine(number_of_modes, azimuth_quadrature_points);
  for (Index m = 0; m < number_of_modes; ++m)
    for (Index k = 0; k < azimuth_quadrature_points; ++k) {
      cosine[m, k] = std::cos(static_cast<Numeric>(m) * (Constant::two_pi * nodes[k]));
      sine[m, k]   = std::sin(static_cast<Numeric>(m) * (Constant::two_pi * nodes[k]));
    }

  // The BPrDF at the nodes, sampled once for all modes
  rtepack::muelmat_vector samples(azimuth_quadrature_points);
  const auto              sample = [&](const Numeric outgoing, const Numeric incoming) {
    for (Index k = 0; k < azimuth_quadrature_points; ++k)
      samples[k] = raw(outgoing, std::abs(incoming), Constant::two_pi * nodes[k]);
  };

  // The ordinary cosine and sine coefficients (1 / 2 pi) int R cos(m phi) dphi of mode m, by the rule on [0, 2 pi],
  // with the factor 2 (2 - delta_m0); the sine one also with sine_sign
  const auto ordinary = [&](const Index m, const Numeric sine_sign) {
    std::pair<rtepack::muelmat, rtepack::muelmat> cs{rtepack::muelmat{0.0}, rtepack::muelmat{0.0}};
    for (Index k = 0; k < azimuth_quadrature_points; ++k) {
      cs.first  += weights[k] * cosine[m, k] * samples[k];
      cs.second += weights[k] * sine[m, k] * samples[k];
    }
    const Numeric factor  = 2.0 * (m == 0 ? 1.0 : 2.0);
    cs.first             *= factor;
    cs.second            *= sine_sign * factor;
    return cs;
  };

  for (Size i = 0; i < mu_out.size(); ++i) {
    for (Size j = 0; j < mu_in.size(); ++j) {
      sample(mu_out[i], mu_in[j]);
      for (Index m = 0; m < number_of_modes; ++m) {
        const auto [ordinary_cosine, ordinary_sine] = ordinary(m, 1.0);
        fill_combined(
            out.brdf()[cosine_mode, m, i, j], out.brdf()[sine_mode, m, i, j], ordinary_cosine, ordinary_sine, m);
      }
    }
    if (std::isnan(mu0)) continue;

    // The beam has only cosine azimuth terms; VDISORT reconstructs its modes with phi_in - phi_out, whereas the
    // raw BPrDF is parameterized by phi_out - phi_in
    sample(mu_out[i], mu0);
    for (Index m = 0; m < number_of_modes; ++m) {
      const auto [ordinary_cosine, ordinary_sine] = ordinary(m, -1.0);
      for (Index row = 0; row < 4; ++row)
        for (Index column = 0; column < 4; ++column) {
          out.beam_brdf()[cosine_mode, m, i][row, column] =
              row < 2 ? ordinary_cosine[row, column] : ordinary_sine[row, column];
          out.beam_brdf()[sine_mode, m, i][row, column] =
              row < 2 ? ordinary_sine[row, column] : ordinary_cosine[row, column];
        }
    }
  }
  return out;
}

surface_reflection depolarizing_fourier_modes(const ScalarRawFunction& raw,
                                              const Index              number_of_modes,
                                              const ConstVectorView&   mu_out,
                                              const ConstVectorView&   mu_in,
                                              const Numeric            mu0,
                                              const Index              azimuth_quadrature_points) {
  if (not raw) throw std::invalid_argument("A raw scalar BRDF function is required");
  return fourier_modes(
      [&raw](const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) {
        rtepack::muelmat result{0.0};
        result[0, 0] = raw(outgoing_mu, incoming_mu, relative_azimuth);
        return result;
      },
      number_of_modes,
      mu_out,
      mu_in,
      mu0,
      azimuth_quadrature_points);
}

surface_reflection hapke_fourier_modes(const Numeric          opposition_amplitude,
                                       const Numeric          opposition_width,
                                       const Numeric          single_scattering_albedo,
                                       const Index            number_of_modes,
                                       const ConstVectorView& mu_out,
                                       const ConstVectorView& mu_in,
                                       const Numeric          mu0,
                                       const Index            azimuth_quadrature_points) {
  return depolarizing_fourier_modes(
      [=](const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) {
        return disort_common::hapke_brdf(outgoing_mu,
                                         incoming_mu,
                                         relative_azimuth,
                                         opposition_amplitude,
                                         opposition_width,
                                         single_scattering_albedo);
      },
      number_of_modes,
      mu_out,
      mu_in,
      mu0,
      azimuth_quadrature_points);
}

surface_reflection rpv_fourier_modes(const Numeric          rho0,
                                     const Numeric          kappa,
                                     const Numeric          asymmetry,
                                     const Numeric          hotspot,
                                     const Index            number_of_modes,
                                     const ConstVectorView& mu_out,
                                     const ConstVectorView& mu_in,
                                     const Numeric          mu0,
                                     const Index            azimuth_quadrature_points) {
  return depolarizing_fourier_modes(
      [=](const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) {
        return disort_common::rpv_brdf(outgoing_mu, incoming_mu, relative_azimuth, rho0, kappa, asymmetry, hotspot);
      },
      number_of_modes,
      mu_out,
      mu_in,
      mu0,
      azimuth_quadrature_points);
}

surface_reflection ross_li_fourier_modes(const Numeric          isotropic,
                                         const Numeric          volumetric,
                                         const Numeric          geometric,
                                         const Numeric          hotspot_angle,
                                         const Index            number_of_modes,
                                         const ConstVectorView& mu_out,
                                         const ConstVectorView& mu_in,
                                         const Numeric          mu0,
                                         const Index            azimuth_quadrature_points) {
  return depolarizing_fourier_modes(
      [=](const Numeric outgoing_mu, const Numeric incoming_mu, const Numeric relative_azimuth) {
        return disort_common::ross_li_brdf(
            outgoing_mu, incoming_mu, relative_azimuth, isotropic, volumetric, geometric, hotspot_angle);
      },
      number_of_modes,
      mu_out,
      mu_in,
      mu0,
      azimuth_quadrature_points);
}

surface_reflection cox_munk_fourier_modes(const Numeric          wind_speed,
                                          const Complex          refractive_index,
                                          const bool             shadowing,
                                          const Index            number_of_modes,
                                          const ConstVectorView& mu_out,
                                          const ConstVectorView& mu_in,
                                          const Numeric          mu0,
                                          const Index            azimuth_quadrature_points) {
  return fourier_modes(
      CoxMunk{wind_speed, refractive_index, shadowing}, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points);
}

surface_reflection fresnel_fourier_modes(const Complex          refractive_index,
                                         const Index            number_of_modes,
                                         const ConstVectorView& mu_out,
                                         const ConstVectorView& mu_in,
                                         const Numeric          mu0) {
  const Fresnel fresnel{refractive_index};
  // Validate eagerly, including when there are no directions.
  static_cast<void>(fresnel(1.0));

  // A flat surface reflects specularly only: no reflection kernel, and the Fresnel matrix as the specular part (which
  // cannot reflect a beam, so mu0 must be NaN)
  auto out = zero_reflection(number_of_modes, mu_out, mu_in, mu0, true);
  for (Size i = 0; i < mu_out.size(); ++i) out.specular()[i] = fresnel(mu_out[i]);
  return out;
}

surface_reflection combine_fourier_modes(const surface_reflection& first,
                                         const Numeric             first_weight,
                                         const surface_reflection& second,
                                         const Numeric             second_weight) {
  disort_common::check_surface_weight(first_weight);
  disort_common::check_surface_weight(second_weight);
  const Index n_out = first.directions(), n_in = first.streams();
  if (n_out != second.directions() or n_in != second.streams())
    throw std::invalid_argument(std::format(
        "Surface reflections can only be combined on the same directions, got {:B,} and {:B,} (2, modes, out, in)",
        first.brdf().shape(),
        second.brdf().shape()));
  if (not(first.beam_mu() == second.beam_mu() or (std::isnan(first.beam_mu()) and std::isnan(second.beam_mu()))))
    throw std::invalid_argument(std::format(
        "Surface reflections can only be combined for the same beam, got the beam cosines {} and {} (NaN: none)",
        first.beam_mu(),
        second.beam_mu()));

  // The weighted sum of the specular parts, or none if no weighted component has one
  const auto has = [](const surface_reflection& x, const Numeric weight) { return weight != 0.0 and x.has_specular(); };
  const bool specular = has(first, first_weight) or has(second, second_weight);

  surface_reflection out(std::max(first.modes(), second.modes()), n_out, n_in, first.beam_mu(), specular);

  const auto add = [&out](const surface_reflection& x, const Numeric weight) {
    if (weight == 0.0) return;
    for (Index alpha = 0; alpha < 2; ++alpha)
      for (Index m = 0; m < x.modes(); ++m)
        for (Index i = 0; i < x.directions(); ++i) {
          for (Index j = 0; j < x.streams(); ++j) out.brdf()[alpha, m, i, j] += weight * x.brdf()[alpha, m, i, j];
          out.beam_brdf()[alpha, m, i] += weight * x.beam_brdf()[alpha, m, i];
        }
  };
  add(first, first_weight);
  add(second, second_weight);

  for (Index i = 0; i < n_out and specular; ++i) {
    if (has(first, first_weight)) out.specular()[i] += first_weight * first.specular()[i];
    if (has(second, second_weight)) out.specular()[i] += second_weight * second.specular()[i];
  }
  return out;
}

surface_reflection fresnel_lambertian_fourier_modes(const Numeric          fresnel_fraction,
                                                    const Numeric          lambertian_albedo,
                                                    const Complex          refractive_index,
                                                    const Index            number_of_modes,
                                                    const ConstVectorView& mu_out,
                                                    const ConstVectorView& mu_in,
                                                    const Numeric          mu0) {
  disort_common::check_surface_fraction(fresnel_fraction);
  return combine_fourier_modes(fresnel_fourier_modes(refractive_index, number_of_modes, mu_out, mu_in, mu0),
                               fresnel_fraction,
                               lambertian_fourier_modes(lambertian_albedo, number_of_modes, mu_out, mu_in, mu0),
                               1.0 - fresnel_fraction);
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
  return combine_fourier_modes(
      cox_munk_fourier_modes(
          wind_speed, refractive_index, shadowing, number_of_modes, mu_out, mu_in, mu0, azimuth_quadrature_points),
      cox_munk_fraction,
      lambertian_fourier_modes(lambertian_albedo, number_of_modes, mu_out, mu_in, mu0),
      1.0 - cox_munk_fraction);
}

}  // namespace vdisort::brdf

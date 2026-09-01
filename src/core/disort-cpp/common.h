#pragma once

#include <matpack.h>

#include <array>
#include <cmath>
#include <complex>
#include <concepts>
#include <functional>
#include <limits>
#include <type_traits>
#include <utility>

namespace disort_common {

inline constexpr Numeric conservative_omega_tolerance = 1.0e-8;
inline constexpr Numeric conservative_kappa_limit     = 1.0e-4;
inline constexpr Numeric centered_integral_limit      = 1.0e-5;

struct flux_values {
  Numeric up{};
  Numeric down_diffuse{};
  Numeric down_direct{};
  Numeric dfdt{};
};

enum class ims_convention : bool {
  /** DISORT 4.0.99 INTCOR: subtract IMS within the 10-degree aureole. */
  disort,
  /** Pythonic-DISORT NT correction: add IMS in every downward direction. */
  pythonic_disort,
};

struct diffuse_moments {
  Numeric upward{};
  Numeric downward{};
  Numeric mean_intensity{};
};

struct centered_pair_basis {
  Numeric cosh{};
  Numeric s_sinhc{};
  Numeric kappa_sinh{};
};

/** Fresnel amplitude coefficients in the local vertical/horizontal basis. */
struct fresnel_amplitudes {
  Complex vertical{};
  Complex horizontal{};
};

/** Shared scalar part of a Cox-Munk microfacet reflection.
 *
 * `factor` contains the slope probability, projected-area denominator, and
 * optional shadowing.  Multiplying it by a Fresnel Mueller matrix gives the
 * polarized BPrDF; multiplying it by the mean squared Fresnel amplitudes gives
 * the scalar BRDF.
 */
struct cox_munk_optics {
  Numeric            factor{};
  fresnel_amplitudes amplitudes{};
};

/** Return the physical Lambertian BRDF A/pi after validating its albedo. */
[[nodiscard]] Numeric lambertian_brdf(Numeric albedo);

/** Evaluate the shared scalar Hapke BRDF used by both transport cores. */
[[nodiscard]] Numeric hapke_brdf(Numeric outgoing_mu,
                                 Numeric incoming_mu,
                                 Numeric relative_azimuth,
                                 Numeric opposition_amplitude,
                                 Numeric opposition_width,
                                 Numeric single_scattering_albedo);

/** Evaluate the shared scalar RPV BRDF used by both transport cores. */
[[nodiscard]] Numeric rpv_brdf(Numeric outgoing_mu,
                               Numeric incoming_mu,
                               Numeric relative_azimuth,
                               Numeric rho0,
                               Numeric kappa,
                               Numeric asymmetry,
                               Numeric hotspot);

/** Evaluate the shared scalar Ross-Li BRDF used by both transport cores. */
[[nodiscard]] Numeric ross_li_brdf(Numeric outgoing_mu,
                                   Numeric incoming_mu,
                                   Numeric relative_azimuth,
                                   Numeric isotropic,
                                   Numeric volumetric,
                                   Numeric geometric,
                                   Numeric hotspot_angle);

/** Validate a nonnegative finite surface-mixture weight. */
void check_surface_weight(Numeric weight);

/** Validate a finite surface-mixture fraction in [0, 1]. */
void check_surface_fraction(Numeric fraction);

/** Compute dielectric Fresnel amplitudes, including total internal reflection. */
[[nodiscard]] fresnel_amplitudes dielectric_fresnel_amplitudes(Numeric incident_mu, Complex refractive_index);

/** Compute the shared optical and microfacet factor of Cox-Munk reflection. */
[[nodiscard]] cox_munk_optics cox_munk_reflection(Numeric outgoing_mu,
                                                  Numeric incoming_mu,
                                                  Numeric relative_azimuth,
                                                  Numeric wind_speed,
                                                  Complex refractive_index,
                                                  bool    shadowing);

/** Evaluate sinh(x)/x with a cancellation-free series near zero. */
[[nodiscard]] inline Numeric sinhc(const Numeric x) {
  if (std::abs(x) >= 1.0e-4) return std::sinh(x) / x;
  const Numeric x2 = x * x;
  return 1.0 + x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 / 5040.0));
}

/** Evaluate (1 - exp(-x))/x accurately, including at x = 0. */
[[nodiscard]] Numeric phi1_neg(Numeric x);

/** Evaluate (exp(x) - 1 - x)/x^2 accurately, including at x = 0. */
[[nodiscard]] Numeric phi2(Numeric x);

/** Compute scale-normalized first-form barycentric interpolation weights. */
void barycentric_weights(Vector& weights, const ConstVectorView& nodes);

/** Initialize the positive/negative double-Gauss streams and inverse cosines. */
void initialize_streams(VectorView mu, VectorView inv_mu, VectorView weights);

/** Validate that layer-indexed optical-depth and albedo arrays have consistent sizes. */
void check_layer_input_sizes(Index n_layers, const AscendingGrid& tau, const ConstVectorView& omega);

/** Validate optical depths, single-scattering albedos, and beam angles. */
void check_layer_input_values(const AscendingGrid& tau, const ConstVectorView& omega, Numeric mu0, Numeric phi0);

/** Return the layer containing an optical depth, assigning interfaces to the layer below. */
[[nodiscard]] Index layer_index(const AscendingGrid& tau, Numeric value);

/** Integrate the linear member of a centered conservative-mode pair along a user ray. */
[[nodiscard]] Numeric user_angle_centered_first_moment(
    Numeric center, Numeric lower, Numeric upper, Numeric observation, Numeric abs_mu, bool downward);

/** Interpolate the layer-endpoint TMS source and integrate it along a downward ray. */
[[nodiscard]] Numeric downward_tms_kernel(
    Numeric mu, Numeric mu0, Numeric thickness, Numeric value_at_lower_end, Numeric value_at_upper_end);

/** Evaluate the depth-dependent IMS attenuation kernel. */
[[nodiscard]] Numeric ims_chi(Numeric tau, Numeric mu, Numeric scaled_mu0);

/** Form the centered hyperbolic basis for eigenvalues +kappa and -kappa. */
[[nodiscard]] inline centered_pair_basis centered_pair(const Numeric kappa, const Numeric distance) {
  const Numeric z = kappa * distance;
  return {.cosh = std::cosh(z), .s_sinhc = distance * sinhc(z), .kappa_sinh = kappa * std::sinh(z)};
}

/** Integrate both centered conservative-pair basis functions along a user ray. */
[[nodiscard]] std::array<Numeric, 2> centered_pair_integrals(
    Numeric kappa, Numeric center, Numeric lower, Numeric upper, Numeric observation, Numeric abs_mu, bool downward);

/** Select the centered representation in the conservative or small-eigenvalue regime. */
[[nodiscard]] constexpr bool use_centered_pair(const Numeric omega, const Numeric kappa) {
  return omega >= 1.0 - conservative_omega_tolerance or kappa <= conservative_kappa_limit;
}

/** Evaluate (exp(z) - 1)/z for real or complex z without cancellation near zero. */
template <typename Scalar> requires(std::same_as<Scalar, Numeric> or std::same_as<Scalar, Complex>)
[[nodiscard]] Scalar exprel(const Scalar z) {
  if (z == Scalar{0.0}) return Scalar{1.0};

  if constexpr (std::same_as<Scalar, Numeric>) {
    return std::expm1(z) / z;
  } else {
    if (std::abs(z) > 0.5) return (std::exp(z) - Scalar{1.0}) / z;

    Scalar          term      = 1.0;
    Scalar          sum       = term;
    constexpr Index max_terms = 4 * std::numeric_limits<Numeric>::max_digits10;
    for (Index n = 1; n <= max_terms; ++n) {
      term *= z / static_cast<Numeric>(n + 1);
      sum  += term;
      if (std::abs(term) <= std::numeric_limits<Numeric>::epsilon() * std::abs(sum)) break;
    }
    return sum;
  }
}

/** Analytically integrate one exponential eigenmode along a signed user ray. */
template <typename Scalar> requires(std::same_as<Scalar, Numeric> or std::same_as<Scalar, Complex>)
[[nodiscard]] Scalar user_angle_exponential_integral(const Scalar  k,
                                                     const Numeric reference,
                                                     const Numeric lower,
                                                     const Numeric upper,
                                                     const Numeric observation,
                                                     const Numeric abs_mu,
                                                     const bool    downward) {
  if (upper <= lower) return Scalar{0.0};
  const auto exponent = [&](const Numeric optical_depth) {
    const Numeric distance = downward ? observation - optical_depth : optical_depth - observation;
    return k * (optical_depth - reference) - distance / abs_mu;
  };
  const Scalar  e0     = exponent(lower);
  const Scalar  e1     = exponent(upper);
  const Scalar  z      = e1 - e0;
  const Numeric z_real = [&] {
    if constexpr (std::same_as<Scalar, Numeric>)
      return z;
    else
      return z.real();
  }();
  const Scalar average = z_real > 0.0 ? std::exp(e1) * exprel(-z) : std::exp(e0) * exprel(z);
  return (upper - lower) / abs_mu * average;
}

/** Interpolate scalar or Stokes-valued samples with precomputed barycentric weights. */
template <typename Values> [[nodiscard]] auto barycentric_interpolate(const ConstVectorView& nodes,
                                                                      const ConstVectorView& weights,
                                                                      const Values&          values,
                                                                      const Numeric          x) {
  using Value = std::remove_cvref_t<decltype(values[0])>;
  Value   numerator{};
  Numeric denominator = 0.0;
  for (Index i = 0; i < static_cast<Index>(nodes.size()); ++i) {
    if (x == nodes[i]) return Value{values[i]};
    const Numeric term  = weights[i] / (x - nodes[i]);
    numerator          += term * values[i];
    denominator        += term;
  }
  return numerator / denominator;
}

/** Transform a value/derivative pair through a centered hyperbolic basis. */
template <typename Scalar> [[nodiscard]] std::array<Scalar, 2> centered_pair_columns(const Scalar&              x,
                                                                                     const Scalar&              r,
                                                                                     const centered_pair_basis& basis) {
  return {basis.cosh * x + basis.kappa_sinh * r, basis.s_sinhc * x + basis.cosh * r};
}

/** Transform centered-pair coefficients to value and derivative amplitudes. */
template <typename Scalar>
[[nodiscard]] std::array<Scalar, 2> centered_pair_amplitudes(const Scalar&              c0,
                                                             const Scalar&              c1,
                                                             const centered_pair_basis& basis) {
  if constexpr (std::same_as<Scalar, Numeric>) {
    return {std::fma(c1, basis.s_sinhc, c0 * basis.cosh), std::fma(c0, basis.kappa_sinh, c1 * basis.cosh)};
  } else {
    return {c0 * basis.cosh + c1 * basis.s_sinhc, c0 * basis.kappa_sinh + c1 * basis.cosh};
  }
}

/** Integrate upward, downward, and mean diffuse intensity over double-Gauss streams. */
template <typename Getter> [[nodiscard]] diffuse_moments integrate_diffuse(const ConstVectorView& positive_mu,
                                                                           const ConstVectorView& weights,
                                                                           Getter&&               intensity) {
  diffuse_moments result;
  const Index     n = static_cast<Index>(weights.size());
  for (Index i = 0; i < n; ++i) {
    const Numeric upward    = std::invoke(intensity, i);
    const Numeric downward  = std::invoke(intensity, n + i);
    result.upward          += weights[i] * positive_mu[i] * upward;
    result.downward        += weights[i] * positive_mu[i] * downward;
    result.mean_intensity  += 0.5 * weights[i] * (upward + downward);
  }
  return result;
}

/** Evaluate a scalar or vector-valued polynomial with Horner's method. */
template <typename Getter>
[[nodiscard]] auto horner_polynomial(const Index count, const Numeric x, Getter&& coefficient) {
  using Value = std::remove_cvref_t<std::invoke_result_t<Getter, Index>>;
  Value value{};
  for (Index i = count - 1; i >= 0; --i) {
    if constexpr (std::same_as<Value, Numeric>)
      value = std::fma(value, x, std::invoke(coefficient, i));
    else
      value = value * x + std::invoke(coefficient, i);
  }
  return value;
}

/** Return the attenuated direct-beam irradiance normal to a horizontal surface. */
[[nodiscard]] inline Numeric direct_beam_flux(const Numeric intensity, const Numeric mu0, const Numeric tau) {
  return mu0 * intensity * std::exp(-tau / mu0);
}

/** Return the attenuated direct-beam radiance normalization. */
[[nodiscard]] inline Numeric direct_beam_radiance(const Numeric intensity, const Numeric mu0, const Numeric tau) {
  return intensity * std::exp(-tau / mu0);
}

}  // namespace disort_common

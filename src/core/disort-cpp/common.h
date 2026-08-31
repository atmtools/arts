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

[[nodiscard]] inline Numeric sinhc(const Numeric x) {
  if (std::abs(x) >= 1.0e-4) return std::sinh(x) / x;
  const Numeric x2 = x * x;
  return 1.0 + x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 / 5040.0));
}
[[nodiscard]] Numeric phi1_neg(Numeric x);
[[nodiscard]] Numeric phi2(Numeric x);

void barycentric_weights(Vector& weights, const ConstVectorView& nodes);

void initialize_streams(VectorView mu, VectorView inv_mu, VectorView weights);

void check_layer_input_sizes(Index n_layers, const AscendingGrid& tau, const ConstVectorView& omega);
void check_layer_input_values(const AscendingGrid& tau, const ConstVectorView& omega, Numeric mu0, Numeric phi0);
[[nodiscard]] Index layer_index(const AscendingGrid& tau, Numeric value);

[[nodiscard]] Numeric user_angle_centered_first_moment(
    Numeric center, Numeric lower, Numeric upper, Numeric observation, Numeric abs_mu, bool downward);

[[nodiscard]] Numeric downward_tms_kernel(
    Numeric mu, Numeric mu0, Numeric thickness, Numeric value_at_lower_end, Numeric value_at_upper_end);
[[nodiscard]] Numeric ims_chi(Numeric tau, Numeric mu, Numeric scaled_mu0);

[[nodiscard]] inline centered_pair_basis centered_pair(const Numeric kappa, const Numeric distance) {
  const Numeric z = kappa * distance;
  return {.cosh = std::cosh(z), .s_sinhc = distance * sinhc(z), .kappa_sinh = kappa * std::sinh(z)};
}

[[nodiscard]] std::array<Numeric, 2> centered_pair_integrals(
    Numeric kappa, Numeric center, Numeric lower, Numeric upper, Numeric observation, Numeric abs_mu, bool downward);

[[nodiscard]] constexpr bool use_centered_pair(const Numeric omega, const Numeric kappa) {
  return omega >= 1.0 - conservative_omega_tolerance or kappa <= conservative_kappa_limit;
}

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

template <typename Scalar> [[nodiscard]] std::array<Scalar, 2> centered_pair_columns(const Scalar&              x,
                                                                                     const Scalar&              r,
                                                                                     const centered_pair_basis& basis) {
  return {basis.cosh * x + basis.kappa_sinh * r, basis.s_sinhc * x + basis.cosh * r};
}

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

[[nodiscard]] inline Numeric direct_beam_flux(const Numeric intensity, const Numeric mu0, const Numeric tau) {
  return mu0 * intensity * std::exp(-tau / mu0);
}

[[nodiscard]] inline Numeric direct_beam_radiance(const Numeric intensity, const Numeric mu0, const Numeric tau) {
  return intensity * std::exp(-tau / mu0);
}

}  // namespace disort_common

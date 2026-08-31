#include "vdisort.h"

#include <arts_constants.h>
#include <debug.h>
#include <legendre.h>
#include <rtepack_multitype.h>
#include <time_report.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <ranges>

namespace {
[[nodiscard]] bool polarized_source(const rtepack::stokvec& stokes) { return stokes.I() > 0.0; }

rtepack::stokvec to_stokvec(const ConstVectorView block) {
  ARTS_USER_ERROR_IF(
      block.size() != vdisort::stokes_dimension, "A Stokes vector must have four components, got {}", block.size());
  return {block[0], block[1], block[2], block[3]};
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

Numeric barycentric_interpolate(const ConstVectorView& nodes,
                                const ConstVectorView& weights,
                                const ConstVectorView& values,
                                const Numeric          x) {
  Numeric numerator   = 0.0;
  Numeric denominator = 0.0;
  for (Index i = 0; i < static_cast<Index>(nodes.size()); ++i) {
    if (x == nodes[i]) return values[i];
    const Numeric term  = weights[i] / (x - nodes[i]);
    numerator          += term * values[i];
    denominator        += term;
  }
  return numerator / denominator;
}

Complex exprel(const Complex z) {
  if (z == Complex{0.0}) return 1.0;

  // The direct quotient is well conditioned away from zero.  Around zero,
  // sum exprel(z) = sum_n z^n/(n+1)! to a floating-point stopping criterion.
  if (std::abs(z) > 0.5) return (std::exp(z) - 1.0) / z;

  Complex         term      = 1.0;
  Complex         sum       = term;
  constexpr Index max_terms = 4 * std::numeric_limits<Numeric>::max_digits10;
  for (Index n = 1; n <= max_terms; ++n) {
    term *= z / static_cast<Numeric>(n + 1);
    sum  += term;
    if (std::abs(term) <= std::numeric_limits<Numeric>::epsilon() * std::abs(sum)) break;
  }
  return sum;
}

Complex user_angle_exponential_integral(const Complex k,
                                        const Numeric reference,
                                        const Numeric lower,
                                        const Numeric upper,
                                        const Numeric observation,
                                        const Numeric abs_mu,
                                        const bool    downward) {
  if (upper <= lower) return 0.0;
  const auto exponent = [&](const Numeric optical_depth) {
    const Numeric distance = downward ? observation - optical_depth : optical_depth - observation;
    return k * (optical_depth - reference) - distance / abs_mu;
  };
  const Complex e0      = exponent(lower);
  const Complex e1      = exponent(upper);
  const Complex z       = e1 - e0;
  const Complex average = z.real() > 0.0 ? std::exp(e1) * exprel(-z) : std::exp(e0) * exprel(z);
  return (upper - lower) / abs_mu * average;
}

void exponential_moments(Vector& moments, const Numeric z) {
  const Index count = static_cast<Index>(moments.size());
  if (count == 0) return;
  if (z == 0.0) {
    for (Index degree = 0; degree < count; ++degree) moments[degree] = 1.0 / static_cast<Numeric>(degree + 1);
    return;
  }

  if (std::abs(z) <= 1.0) {
    constexpr Index max_terms = 4 * std::numeric_limits<Numeric>::max_digits10;
    for (Index degree = 0; degree < count; ++degree) {
      Numeric term = 1.0 / static_cast<Numeric>(degree + 1);
      Numeric sum  = term;
      for (Index n = 1; n <= max_terms; ++n) {
        term *= z / static_cast<Numeric>(n) * static_cast<Numeric>(n + degree) / static_cast<Numeric>(n + degree + 1);
        sum  += term;
        if (std::abs(term) <= std::numeric_limits<Numeric>::epsilon() * std::abs(sum)) break;
      }
      moments[degree] = sum;
    }
    return;
  }

  const Numeric exponential = std::exp(z);
  moments[0]                = -std::expm1(z) / -z;
  for (Index degree = 1; degree < count; ++degree)
    moments[degree] = (exponential - static_cast<Numeric>(degree) * moments[degree - 1]) / z;
}

rtepack::stokvec user_angle_polynomial_integral(const rtepack::stokvec_vector& coefficients,
                                                const Numeric                  source_scale,
                                                const Numeric                  source_offset,
                                                const Numeric                  lower,
                                                const Numeric                  upper,
                                                const Numeric                  observation,
                                                const Numeric                  abs_mu,
                                                const bool                     downward) {
  rtepack::stokvec result{};
  const Index      count = static_cast<Index>(coefficients.size());
  if (count == 0 or upper <= lower) return result;

  const Numeric width    = upper - lower;
  const Numeric endpoint = downward ? upper : lower;
  const Numeric distance = downward ? observation - upper : lower - observation;
  const Numeric base     = std::fma(source_scale, endpoint, source_offset);
  const Numeric step     = (downward ? -1.0 : 1.0) * source_scale * width;
  const Numeric z        = -width / abs_mu;

  Vector moments(count);
  exponential_moments(moments, z);
  for (Index power = 0; power < count; ++power) {
    Numeric choose = 1.0;
    for (Index degree = 0; degree <= power; ++degree) {
      const Numeric factor  = choose * std::pow(base, power - degree) * std::pow(step, degree) * moments[degree];
      result               += factor * coefficients[power];
      if (degree < power) choose *= static_cast<Numeric>(power - degree) / static_cast<Numeric>(degree + 1);
    }
  }
  result *= width / abs_mu * std::exp(-distance / abs_mu);
  return result;
}

rtepack::stokvec_tensor3 to_boundary_data(const Tensor4& in) {
  if (in.size() == 0) return {};
  const auto [nalpha, nfourier, nstream, nstokes] = in.shape();
  ARTS_USER_ERROR_IF(
      nstokes != vdisort::stokes_dimension, "The boundary Stokes dimension must be 4, got {:B,}", in.shape());
  rtepack::stokvec_tensor3 out(nalpha, nfourier, nstream);
  for (Index a = 0; a < nalpha; ++a)
    for (Index m = 0; m < nfourier; ++m)
      for (Index i = 0; i < nstream; ++i) out[a, m, i] = to_stokvec(in[a, m, i]);
  return out;
}

rtepack::stokvec_matrix to_source_poly_data(const Tensor3& in) {
  const auto [nlayers, ncoeffs, nstokes] = in.shape();
  ARTS_USER_ERROR_IF(
      nstokes != vdisort::stokes_dimension, "The source Stokes dimension must be 4, got {:B,}", in.shape());
  rtepack::stokvec_matrix out(nlayers, ncoeffs);
  for (Index l = 0; l < nlayers; ++l)
    for (Index p = 0; p < ncoeffs; ++p) out[l, p] = to_stokvec(in[l, p]);
  return out;
}

void fill_combined(rtepack::muelmat&       out_cos,
                   rtepack::muelmat&       out_sin,
                   const rtepack::muelmat& ordinary_cos,
                   const rtepack::muelmat& ordinary_sin,
                   const Index             m) {
  out_cos = rtepack::muelmat{0.0};
  out_sin = rtepack::muelmat{0.0};

  if (m == 0) {
    // Paper Eq. (82): the two independent m=0 systems carry [I,Q] and [U,V].
    for (Index i = 0; i < 2; ++i)
      for (Index j = 0; j < 2; ++j) out_cos[i, j] = ordinary_cos[i, j];
    for (Index i = 2; i < 4; ++i)
      for (Index j = 2; j < 4; ++j) out_sin[i, j] = ordinary_cos[i, j];
    return;
  }

  // Paper Eq. (81).
  out_cos[0, 0] = ordinary_cos[0, 0];
  out_cos[0, 1] = ordinary_cos[0, 1];
  out_cos[0, 2] = -ordinary_sin[0, 2];
  out_cos[0, 3] = -ordinary_sin[0, 3];
  out_cos[1, 0] = ordinary_cos[1, 0];
  out_cos[1, 1] = ordinary_cos[1, 1];
  out_cos[1, 2] = -ordinary_sin[1, 2];
  out_cos[1, 3] = -ordinary_sin[1, 3];
  out_cos[2, 0] = ordinary_sin[2, 0];
  out_cos[2, 1] = ordinary_sin[2, 1];
  out_cos[2, 2] = ordinary_cos[2, 2];
  out_cos[2, 3] = ordinary_cos[2, 3];
  out_cos[3, 0] = ordinary_sin[3, 0];
  out_cos[3, 1] = ordinary_sin[3, 1];
  out_cos[3, 2] = ordinary_cos[3, 2];
  out_cos[3, 3] = ordinary_cos[3, 3];

  out_sin[0, 0] = ordinary_cos[0, 0];
  out_sin[0, 1] = ordinary_cos[0, 1];
  out_sin[0, 2] = ordinary_sin[0, 2];
  out_sin[0, 3] = ordinary_sin[0, 3];
  out_sin[1, 0] = ordinary_cos[1, 0];
  out_sin[1, 1] = ordinary_cos[1, 1];
  out_sin[1, 2] = ordinary_sin[1, 2];
  out_sin[1, 3] = ordinary_sin[1, 3];
  out_sin[2, 0] = -ordinary_sin[2, 0];
  out_sin[2, 1] = -ordinary_sin[2, 1];
  out_sin[2, 2] = ordinary_cos[2, 2];
  out_sin[2, 3] = ordinary_cos[2, 3];
  out_sin[3, 0] = -ordinary_sin[3, 0];
  out_sin[3, 1] = -ordinary_sin[3, 1];
  out_sin[3, 2] = ordinary_cos[3, 2];
  out_sin[3, 3] = ordinary_cos[3, 3];
}

rtepack::muelmat to_muelmat(const ConstMatrixView block) {
  assert(block.nrows() == vdisort::stokes_dimension and block.ncols() == vdisort::stokes_dimension);
  rtepack::muelmat out{0.0};
  for (Index i = 0; i < vdisort::stokes_dimension; ++i)
    for (Index j = 0; j < vdisort::stokes_dimension; ++j) out[i, j] = block[i, j];
  return out;
}

vdisort::phase_matrix_data to_phase_matrix_data(const Tensor7& in) {
  if (in.size() == 0) return {};
  const auto [nalpha, nfourier, nlayers, nout, nin, ns1, ns2] = in.shape();
  ARTS_USER_ERROR_IF(ns1 != vdisort::stokes_dimension or ns2 != vdisort::stokes_dimension,
                     "The last two phase-matrix dimensions must both be 4, got {:B,}",
                     in.shape());
  vdisort::phase_matrix_data out(nalpha, nfourier, nlayers, nout, nin, rtepack::muelmat{0.0});
  for (Index a = 0; a < nalpha; ++a)
    for (Index m = 0; m < nfourier; ++m)
      for (Index l = 0; l < nlayers; ++l)
        for (Index i = 0; i < nout; ++i)
          for (Index j = 0; j < nin; ++j) out[a, m, l, i, j] = to_muelmat(in[a, m, l, i, j]);
  return out;
}

vdisort::beam_phase_matrix_data to_beam_phase_matrix_data(const Tensor6& in) {
  if (in.size() == 0) return {};
  const auto [nalpha, nfourier, nlayers, nout, ns1, ns2] = in.shape();
  ARTS_USER_ERROR_IF(ns1 != vdisort::stokes_dimension or ns2 != vdisort::stokes_dimension,
                     "The last two beam phase-matrix dimensions must both be 4, got {:B,}",
                     in.shape());
  vdisort::beam_phase_matrix_data out(nalpha, nfourier, nlayers, nout, rtepack::muelmat{0.0});
  for (Index a = 0; a < nalpha; ++a)
    for (Index m = 0; m < nfourier; ++m)
      for (Index l = 0; l < nlayers; ++l)
        for (Index i = 0; i < nout; ++i) out[a, m, l, i] = to_muelmat(in[a, m, l, i]);
  return out;
}
}  // namespace

namespace vdisort {
namespace {
Numeric correction_phi2(const Numeric x) {
  if (std::abs(x) >= 1.0) return (std::expm1(x) - x) / (x * x);
  Numeric term = 0.5;
  Numeric sum  = term;
  for (Index k = 1; k <= 20; ++k) {
    term *= x / static_cast<Numeric>(k + 2);
    sum  += term;
  }
  return sum;
}

Numeric correction_downward_kernel(const Numeric mu,
                                   const Numeric mu0,
                                   const Numeric thickness,
                                   const Numeric lower_value,
                                   const Numeric upper_value) {
  const Numeric y = 1.0 / mu - 1.0 / mu0;
  const Numeric w = thickness * y;
  if (std::abs(w) < 1.0) return thickness / mu * (w == 0.0 ? 1.0 : -std::expm1(-w) / w) * lower_value;
  return (lower_value - upper_value) / (mu * y);
}

Numeric correction_ims_chi(const Numeric tau, const Numeric mu, const Numeric scaled_mu0) {
  const Numeric x = 1.0 / mu - 1.0 / scaled_mu0;
  const Numeric z = -tau * x;
  if (std::abs(z) < 1.0) return tau * tau * std::exp(-tau / scaled_mu0) * correction_phi2(z) / (mu * scaled_mu0);
  return ((tau - 1.0 / x) * std::exp(-tau / scaled_mu0) + std::exp(-tau / mu) / x) / (mu * scaled_mu0 * x);
}
}  // namespace

delta_m_correction_cache::delta_m_correction_cache(AscendingGrid                       physical_tau,
                                                   Vector                              omega,
                                                   Vector                              fraction,
                                                   const Numeric                       mu0,
                                                   const Numeric                       phi0,
                                                   rtepack::stokvec                    beam,
                                                   Vector                              user_mu,
                                                   Vector                              phi,
                                                   lab_phase_function                  original_phase,
                                                   lab_phase_function                  transport_phase,
                                                   lab_phase_function                  removed_phase,
                                                   const Index                         intermediate_mu,
                                                   const Index                         intermediate_phi,
                                                   lab_phase_pair_convolution_function removed_pair_convolution)
    : physical_tau_(std::move(physical_tau)),
      omega_(std::move(omega)),
      fraction_(std::move(fraction)),
      user_mu_(std::move(user_mu)),
      phi_(std::move(phi)),
      mu0_(mu0),
      phi0_(phi0),
      beam_(beam) {
  ARTS_TIME_REPORT

  const Index nlayers = static_cast<Index>(physical_tau_.size());
  const Index nuser   = static_cast<Index>(user_mu_.size());
  const Index nphi    = static_cast<Index>(phi_.size());
  ARTS_USER_ERROR_IF(omega_.size() != physical_tau_.size() or fraction_.size() != physical_tau_.size(),
                     "Polarized delta-M arrays must all contain {} layers, got omega={} and fraction={}",
                     nlayers,
                     omega_.size(),
                     fraction_.size());
  ARTS_USER_ERROR_IF(mu0_ <= 0.0 or mu0_ > 1.0, "The delta-M beam cosine must be in (0, 1], got {}", mu0_);
  ARTS_USER_ERROR_IF(intermediate_mu <= 0 or intermediate_phi <= 0,
                     "IMS quadrature sizes must be positive, got {} x {}",
                     intermediate_mu,
                     intermediate_phi);
  ARTS_USER_ERROR_IF(std::ranges::any_of(user_mu_, [](const Numeric mu) { return mu == 0.0 or std::abs(mu) > 1.0; }),
                     "Delta-M user cosines must be nonzero and in [-1, 1], got {:B,}",
                     user_mu_);

  scale_.resize(nlayers);
  Vector  scaled_tau(nlayers);
  Numeric physical_top = 0.0, scaled_top = 0.0;
  for (Index layer = 0; layer < nlayers; ++layer) {
    ARTS_USER_ERROR_IF(omega_[layer] < 0.0 or omega_[layer] > 1.0 or fraction_[layer] < 0.0 or fraction_[layer] >= 1.0,
                       "Invalid delta-M omega/fraction in layer {}: {}, {}",
                       layer,
                       omega_[layer],
                       fraction_[layer]);
    scale_[layer]      = 1.0 - omega_[layer] * fraction_[layer];
    scaled_top        += scale_[layer] * (physical_tau_[layer] - physical_top);
    scaled_tau[layer]  = scaled_top;
    physical_top       = physical_tau_[layer];
  }
  scaled_tau_ = AscendingGrid{std::move(scaled_tau)};

  tms_operator_.resize(nlayers, nphi, nuser);
  for (Index layer = 0; layer < nlayers; ++layer) {
    const Numeric transport_omega = omega_[layer] * (1.0 - fraction_[layer]) / scale_[layer];
    for (Index p = 0; p < nphi; ++p)
      for (Index u = 0; u < nuser; ++u) {
        const Numeric beam_factor = user_mu_[u] > 0.0 ? mu0_ / (mu0_ + user_mu_[u]) : 1.0;
        tms_operator_[layer, p, u] =
            transport_omega * Constant::inv_pi * 0.25 * beam_factor *
            (original_phase(layer, user_mu_[u], phi_[p], -mu0_, phi0_) / (1.0 - fraction_[layer]) -
             transport_phase(layer, user_mu_[u], phi_[p], -mu0_, phi0_));
      }
  }

  Vector mid_mu(intermediate_mu), mid_weight(intermediate_mu);
  Legendre::GaussLegendre(mid_mu, mid_weight);
  ims_operator_.resize(nlayers + 1, nphi, nuser);
  ims_operator_ = rtepack::muelmat{0.0};
  ims_scalar_.resize(nlayers + 1);
  ims_mu0_.resize(nlayers + 1);
  Numeric cumulative_peak = 0.0;
  Vector  layer_peak_weight(nlayers, 0.0);
  physical_top = 0.0;
  for (Index boundary = 1; boundary <= nlayers; ++boundary) {
    const Index layer           = boundary - 1;
    layer_peak_weight[layer]    = omega_[layer] * fraction_[layer] * (physical_tau_[layer] - physical_top);
    cumulative_peak            += layer_peak_weight[layer];
    physical_top                = physical_tau_[layer];
    const Numeric average_peak  = cumulative_peak / physical_tau_[layer];
    ims_mu0_[boundary]          = mu0_ / (1.0 - average_peak);
    ims_scalar_[boundary]       = Constant::inv_pi * 0.25 * average_peak * average_peak / (1.0 - average_peak);
    if (cumulative_peak == 0.0) continue;

    const auto average_removed =
        [&](const Numeric out_mu, const Numeric out_phi, const Numeric in_mu, const Numeric in_phi) {
          rtepack::muelmat value{0.0};
          for (Index l = 0; l < boundary; ++l)
            if (layer_peak_weight[l] != 0.0)
              value += layer_peak_weight[l] * removed_phase(l, out_mu, out_phi, in_mu, in_phi);
          return value / cumulative_peak;
        };

    for (Index p = 0; p < nphi; ++p)
      for (Index u = 0; u < nuser; ++u) {
        const rtepack::muelmat direct = average_removed(user_mu_[u], phi_[p], -mu0_, phi0_);
        rtepack::muelmat       convolution{0.0};
        if (removed_pair_convolution) {
          for (Index first = 0; first < boundary; ++first)
            for (Index second = 0; second < boundary; ++second)
              if (layer_peak_weight[first] != 0.0 and layer_peak_weight[second] != 0.0)
                convolution += layer_peak_weight[first] * layer_peak_weight[second] *
                               removed_pair_convolution(first, second, user_mu_[u], phi_[p], -mu0_, phi0_);
          convolution /= cumulative_peak * cumulative_peak;
        } else {
          for (Index q = 0; q < intermediate_mu; ++q)
            for (Index a = 0; a < intermediate_phi; ++a) {
              const Numeric mid_phi =
                  Constant::two_pi * (static_cast<Numeric>(a) + 0.5) / static_cast<Numeric>(intermediate_phi);
              convolution += mid_weight[q] * (average_removed(user_mu_[u], phi_[p], mid_mu[q], mid_phi) *
                                              average_removed(mid_mu[q], mid_phi, -mu0_, phi0_));
            }
          convolution *= 0.5 / static_cast<Numeric>(intermediate_phi);
        }

        ims_operator_[boundary, p, u] = 2.0 * direct - convolution;
      }
  }
  ims_operator_[0] = ims_operator_[1];
  ims_scalar_[0]   = ims_scalar_[1];
  ims_mu0_[0]      = ims_mu0_[1];
}

rtepack::stokvec_vector delta_m_correction_cache::evaluate(const Numeric tau, const Index phi_index) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0.0 or tau > physical_tau_.back(),
                     "Physical correction depth must be in [0, {}], got {}",
                     physical_tau_.back(),
                     tau);
  ARTS_USER_ERROR_IF(phi_index < 0 or phi_index >= static_cast<Index>(phi_.size()),
                     "Correction azimuth index {} is outside [0, {})",
                     phi_index,
                     phi_.size());
  const Index layer = std::min<Index>(
      std::distance(physical_tau_.begin(), std::ranges::lower_bound(physical_tau_, tau)), physical_tau_.size() - 1);
  const Numeric physical_layer_top = layer == 0 ? 0.0 : physical_tau_[layer - 1];
  const Numeric scaled_layer_top   = layer == 0 ? 0.0 : scaled_tau_[layer - 1];
  const Numeric scaled_tau         = scaled_layer_top + scale_[layer] * (tau - physical_layer_top);

  rtepack::stokvec_vector out(user_mu_.size(), rtepack::stokvec{});
  for (Index u = 0; u < static_cast<Index>(user_mu_.size()); ++u) {
    const Numeric mu = user_mu_[u], abs_mu = std::abs(mu);
    if (mu > 0.0) {
      const Numeric at_observation  = std::exp(-scaled_tau / mu0_);
      const Numeric at_bottom       = std::exp((scaled_tau - scaled_tau_[layer]) / abs_mu - scaled_tau_[layer] / mu0_);
      out[u]                       += (at_observation - at_bottom) * (tms_operator_[layer, phi_index, u] * beam_);
      for (Index l = layer + 1; l < static_cast<Index>(physical_tau_.size()); ++l) {
        const Numeric top = scaled_tau_[l - 1], bottom = scaled_tau_[l];
        const Numeric a  = std::exp((scaled_tau - top) / abs_mu - top / mu0_);
        const Numeric b  = std::exp((scaled_tau - bottom) / abs_mu - bottom / mu0_);
        out[u]          += (a - b) * (tms_operator_[l, phi_index, u] * beam_);
      }
    } else {
      const Numeric at_observation = std::exp(-scaled_tau / mu0_);
      const Numeric at_top         = std::exp((scaled_layer_top - scaled_tau) / abs_mu - scaled_layer_top / mu0_);
      out[u] += correction_downward_kernel(abs_mu, mu0_, scaled_tau - scaled_layer_top, at_observation, at_top) *
                (tms_operator_[layer, phi_index, u] * beam_);
      for (Index l = 0; l < layer; ++l) {
        const Numeric top = l == 0 ? 0.0 : scaled_tau_[l - 1], bottom = scaled_tau_[l];
        const Numeric a = std::exp((bottom - scaled_tau) / abs_mu - bottom / mu0_);
        const Numeric b = std::exp((top - scaled_tau) / abs_mu - top / mu0_);
        out[u] +=
            correction_downward_kernel(abs_mu, mu0_, bottom - top, a, b) * (tms_operator_[l, phi_index, u] * beam_);
      }

      const Numeric beam_theta = std::acos(-mu0_);
      const Numeric ray_theta  = std::acos(mu);
      if (std::abs(beam_theta - ray_theta) <= Constant::pi / 18.0) {
        const Numeric layer_bottom        = physical_tau_[layer];
        const Numeric weight              = (tau - physical_layer_top) / (layer_bottom - physical_layer_top);
        const auto    boundary_correction = [&](const Index boundary) {
          return ims_scalar_[boundary] * correction_ims_chi(tau, abs_mu, ims_mu0_[boundary]) *
                 (ims_operator_[boundary, phi_index, u] * beam_);
        };
        out[u] -= (1.0 - weight) * boundary_correction(layer) + weight * boundary_correction(layer + 1);
      }
    }
  }
  return out;
}

phase_matrix_fourier_coefficients phase_matrix_fourier_split(const rtepack::specmat_matrix_const_view& phase_matrix) {
  phase_matrix_fourier_coefficients out{
      .cosine = rtepack::muelmat_matrix(phase_matrix.nrows(), phase_matrix.ncols(), rtepack::muelmat{0.0}),
      .sine   = rtepack::muelmat_matrix(phase_matrix.nrows(), phase_matrix.ncols(), rtepack::muelmat{0.0})};
  for (Index frequency = 0; frequency < phase_matrix.nrows(); ++frequency) {
    for (Index coefficient = 0; coefficient < phase_matrix.ncols(); ++coefficient) {
      for (Index i = 0; i < stokes_dimension; ++i) {
        for (Index j = 0; j < stokes_dimension; ++j) {
          const Complex value                      = phase_matrix[frequency, coefficient][i, j];
          out.cosine[frequency, coefficient][i, j] = value.real();
          out.sine[frequency, coefficient][i, j]   = -value.imag();
        }
      }
    }
  }
  return out;
}

void BDRF::operator()(const Index                  alpha,
                      rtepack::muelmat_matrix_view out,
                      const ConstVectorView&       mu_out,
                      const ConstVectorView&       mu_in) const {
  ARTS_USER_ERROR_IF(alpha != cosine_mode and alpha != sine_mode,
                     "The VDISORT BRDF mode must be 0 (cosine) or 1 (sine), got {}",
                     alpha);
  (alpha == cosine_mode ? cosine : sine)(out, mu_out, mu_in);
}

phase_matrix_data combine_phase_matrices(const rtepack::muelmat_tensor4& cosine, const rtepack::muelmat_tensor4& sine) {
  ARTS_USER_ERROR_IF(cosine.shape() != sine.shape(),
                     "Cosine and sine phase matrices have different shapes: {:B,} and {:B,}",
                     cosine.shape(),
                     sine.shape());
  const auto [nfourier, nlayers, nout, nin] = cosine.shape();
  phase_matrix_data out(2, nfourier, nlayers, nout, nin, rtepack::muelmat{0.0});
  for (Index m = 0; m < nfourier; ++m)
    for (Index l = 0; l < nlayers; ++l)
      for (Index i = 0; i < nout; ++i)
        for (Index j = 0; j < nin; ++j)
          fill_combined(
              out[cosine_mode, m, l, i, j], out[sine_mode, m, l, i, j], cosine[m, l, i, j], sine[m, l, i, j], m);
  return out;
}

Tensor7 combine_phase_matrices(const Tensor6& cosine, const Tensor6& sine) {
  ARTS_USER_ERROR_IF(cosine.shape() != sine.shape(),
                     "Cosine and sine phase matrices have different shapes: {:B,} and {:B,}",
                     cosine.shape(),
                     sine.shape());
  const auto [nfourier, nlayers, nout, nin, ns1, ns2] = cosine.shape();
  ARTS_USER_ERROR_IF(ns1 != stokes_dimension or ns2 != stokes_dimension,
                     "The last two phase-matrix dimensions must both be 4, got {:B,}",
                     cosine.shape());
  rtepack::muelmat_tensor4 c(nfourier, nlayers, nout, nin), s(nfourier, nlayers, nout, nin);
  for (Index m = 0; m < nfourier; ++m)
    for (Index l = 0; l < nlayers; ++l)
      for (Index i = 0; i < nout; ++i)
        for (Index j = 0; j < nin; ++j) {
          c[m, l, i, j] = to_muelmat(cosine[m, l, i, j]);
          s[m, l, i, j] = to_muelmat(sine[m, l, i, j]);
        }
  const auto combined = combine_phase_matrices(c, s);
  Tensor7    out(2, nfourier, nlayers, nout, nin, stokes_dimension, stokes_dimension, 0.0);
  for (Index a = 0; a < 2; ++a)
    for (Index m = 0; m < nfourier; ++m)
      for (Index l = 0; l < nlayers; ++l)
        for (Index i = 0; i < nout; ++i)
          for (Index j = 0; j < nin; ++j)
            for (Index so = 0; so < stokes_dimension; ++so)
              for (Index si = 0; si < stokes_dimension; ++si)
                out[a, m, l, i, j, so, si] = combined[a, m, l, i, j][so, si];
  return out;
}

beam_phase_matrix_data combine_beam_phase_matrices(const rtepack::muelmat_tensor3& cosine,
                                                   const rtepack::muelmat_tensor3& sine) {
  ARTS_USER_ERROR_IF(cosine.shape() != sine.shape(),
                     "Cosine and sine beam phase matrices have different shapes: {:B,} and {:B,}",
                     cosine.shape(),
                     sine.shape());
  const auto [nfourier, nlayers, nout] = cosine.shape();
  beam_phase_matrix_data out(2, nfourier, nlayers, nout, rtepack::muelmat{0.0});
  for (Index m = 0; m < nfourier; ++m)
    for (Index l = 0; l < nlayers; ++l)
      for (Index i = 0; i < nout; ++i)
        fill_combined(out[cosine_mode, m, l, i], out[sine_mode, m, l, i], cosine[m, l, i], sine[m, l, i], m);
  return out;
}

Tensor6 combine_beam_phase_matrices(const Tensor5& cosine, const Tensor5& sine) {
  ARTS_USER_ERROR_IF(cosine.shape() != sine.shape(),
                     "Cosine and sine beam phase matrices have different shapes: {:B,} and {:B,}",
                     cosine.shape(),
                     sine.shape());
  const auto [nfourier, nlayers, nout, ns1, ns2] = cosine.shape();
  ARTS_USER_ERROR_IF(ns1 != stokes_dimension or ns2 != stokes_dimension,
                     "The last two beam phase-matrix dimensions must both be 4, got {:B,}",
                     cosine.shape());
  rtepack::muelmat_tensor3 c(nfourier, nlayers, nout), s(nfourier, nlayers, nout);
  for (Index m = 0; m < nfourier; ++m)
    for (Index l = 0; l < nlayers; ++l)
      for (Index i = 0; i < nout; ++i) {
        c[m, l, i] = to_muelmat(cosine[m, l, i]);
        s[m, l, i] = to_muelmat(sine[m, l, i]);
      }
  const auto combined = combine_beam_phase_matrices(c, s);
  Tensor6    out(2, nfourier, nlayers, nout, stokes_dimension, stokes_dimension, 0.0);
  for (Index a = 0; a < 2; ++a)
    for (Index m = 0; m < nfourier; ++m)
      for (Index l = 0; l < nlayers; ++l)
        for (Index i = 0; i < nout; ++i)
          for (Index so = 0; so < stokes_dimension; ++so)
            for (Index si = 0; si < stokes_dimension; ++si) out[a, m, l, i, so, si] = combined[a, m, l, i][so, si];
  return out;
}

main_data::main_data(
    const Index NLayers_, const Index NQuad_, const Index NFourier_, const Index Nscoeffs_, const Index NBDRF_)
    : NLayers(NLayers_),
      NQuad(NQuad_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(Nscoeffs_),
      NBDRF(NBDRF_),
      NState(stokes_dimension * NQuad),
      NHalfState(stokes_dimension * N),
      has_source_poly(Nscoeffs > 0),
      tau_arr(NLayers),
      omega_arr(NLayers),
      source_poly_coeffs(NLayers, Nscoeffs),
      source_coordinate_scale(NLayers, 1.0),
      source_coordinate_offset(NLayers, 0.0),
      phase_matrix(2, NFourier, NLayers, NQuad, NQuad, rtepack::muelmat{0.0}),
      boundary_up(2, NFourier, N),
      boundary_down(2, NFourier, N),
      brdf_fourier_modes(NBDRF),
      beam_stokes{},
      beam_phase_matrix(2, NFourier, NLayers, NQuad, rtepack::muelmat{0.0}),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      half_range_barycentric_weights(N),
      scaled_source_poly_coeffs(NLayers, Nscoeffs),
      G_collect(2, NFourier, NLayers, NState, NState),
      K_collect(2, NFourier, NLayers, NState),
      GC_collect(2, NFourier, NLayers, NState),
      B_collect(2, NFourier, NLayers, NQuad),
      source_collect(2, NFourier, NLayers, NQuad, Nscoeffs),
      um(NLayers, 2, NFourier, NQuad),
      top_anchored(static_cast<std::size_t>(2 * NFourier * NLayers * NState), 0) {
  ARTS_USER_ERROR_IF(NQuad <= 0 or NQuad % 2 != 0, "NQuad must be a positive even number, got {}", NQuad);
  Legendre::PositiveDoubleGaussLegendre(mu_arr[Range{0, N}], W);
  std::transform(mu_arr.begin(), mu_arr.begin() + N, mu_arr.begin() + N, [](const Numeric x) { return -x; });
  std::transform(mu_arr.begin(), mu_arr.end(), inv_mu_arr.begin(), [](const Numeric x) { return 1.0 / x; });
  barycentric_weights(half_range_barycentric_weights, mu_arr[Range{0, N}]);
}

main_data::main_data(Index             NQuad_,
                     Index             NFourier_,
                     AscendingGrid     tau_arr_,
                     Vector            omega_arr_,
                     Tensor7           phase_matrix_,
                     Tensor4           boundary_up_,
                     Tensor4           boundary_down_,
                     Tensor3           source_poly_coeffs_,
                     std::vector<BDRF> brdf_fourier_modes_,
                     Numeric           mu0_,
                     Vector            beam_stokes_,
                     Numeric           phi0_,
                     Tensor6           beam_phase_matrix_,
                     Vector            source_coordinate_scale_,
                     Vector            source_coordinate_offset_)
    : main_data(NQuad_,
                NFourier_,
                std::move(tau_arr_),
                std::move(omega_arr_),
                to_phase_matrix_data(phase_matrix_),
                to_boundary_data(boundary_up_),
                to_boundary_data(boundary_down_),
                to_source_poly_data(source_poly_coeffs_),
                std::move(brdf_fourier_modes_),
                mu0_,
                to_stokvec(beam_stokes_),
                phi0_,
                to_beam_phase_matrix_data(beam_phase_matrix_),
                std::move(source_coordinate_scale_),
                std::move(source_coordinate_offset_)) {}

main_data::main_data(Index                    NQuad_,
                     Index                    NFourier_,
                     AscendingGrid            tau_arr_,
                     Vector                   omega_arr_,
                     phase_matrix_data        phase_matrix_,
                     rtepack::stokvec_tensor3 boundary_up_,
                     rtepack::stokvec_tensor3 boundary_down_,
                     rtepack::stokvec_matrix  source_poly_coeffs_,
                     std::vector<BDRF>        brdf_fourier_modes_,
                     Numeric                  mu0_,
                     rtepack::stokvec         beam_stokes_,
                     Numeric                  phi0_,
                     beam_phase_matrix_data   beam_phase_matrix_,
                     Vector                   source_coordinate_scale_,
                     Vector                   source_coordinate_offset_)
    : main_data(static_cast<Index>(tau_arr_.size()),
                NQuad_,
                NFourier_,
                source_poly_coeffs_.size() == 0 ? 0 : source_poly_coeffs_.shape()[1],
                static_cast<Index>(brdf_fourier_modes_.size())) {
  tau_arr            = std::move(tau_arr_);
  omega_arr          = std::move(omega_arr_);
  phase_matrix       = std::move(phase_matrix_);
  boundary_up        = std::move(boundary_up_);
  boundary_down      = std::move(boundary_down_);
  source_poly_coeffs = std::move(source_poly_coeffs_);
  if (source_coordinate_scale_.size() != 0) source_coordinate_scale = std::move(source_coordinate_scale_);
  if (source_coordinate_offset_.size() != 0) source_coordinate_offset = std::move(source_coordinate_offset_);
  brdf_fourier_modes = std::move(brdf_fourier_modes_);
  mu0                = mu0_;
  beam_stokes        = std::move(beam_stokes_);
  phi0               = phi0_;
  has_beam_source    = polarized_source(beam_stokes);
  if (beam_phase_matrix_.size() != 0) beam_phase_matrix = std::move(beam_phase_matrix_);

  check_input_size();
  update_all();
}

Index main_data::state_index(const Index stream, const Index stokes) const {
  return stokes_dimension * stream + stokes;
}

Index main_data::anchor_index(const Index alpha, const Index m, const Index layer, const Index eigen) const {
  return (((alpha * NFourier + m) * NLayers + layer) * NState + eigen);
}

Numeric main_data::layer_top(const Index layer) const { return layer == 0 ? 0.0 : tau_arr[layer - 1]; }

Complex main_data::homogeneous(const Index   alpha,
                               const Index   m,
                               const Index   layer,
                               const Index   state,
                               const Index   eigen,
                               const Numeric tau) const {
  const Numeric anchor = top_anchored[anchor_index(alpha, m, layer, eigen)] ? layer_top(layer) : tau_arr[layer];
  return G_collect[alpha, m, layer, state, eigen] * std::exp(K_collect[alpha, m, layer, eigen] * (tau - anchor));
}

Numeric main_data::particular(
    const Index alpha, const Index m, const Index layer, const Index state, const Numeric tau) const {
  Numeric value = 0.0;
  // VDISORT CHANGE: permit a layer-local affine source coordinate.  The
  // scalar-limit adapter uses this to retain DISORT's original-tau source
  // convention while its delta-M transport operator runs on scaled tau.
  const Numeric source_tau = std::fma(source_coordinate_scale[layer], tau, source_coordinate_offset[layer]);
  const Index   stream     = state / stokes_dimension;
  const Index   stokes     = state % stokes_dimension;
  for (Index p = Nscoeffs - 1; p >= 0; --p)
    value = std::fma(value, source_tau, source_collect[alpha, m, layer, stream, p][stokes]);
  if (has_beam_source) value += B_collect[alpha, m, layer, stream][stokes] * std::exp(-tau / mu0);
  return value;
}

void main_data::check_input_size() const {
  ARTS_USER_ERROR_IF(
      static_cast<Index>(tau_arr.size()) != NLayers, "tau_arr has size {}, expected {}", tau_arr.size(), NLayers);
  ARTS_USER_ERROR_IF(
      static_cast<Index>(omega_arr.size()) != NLayers, "omega_arr has size {}, expected {}", omega_arr.size(), NLayers);
  ARTS_USER_ERROR_IF((source_poly_coeffs.shape() != std::array<Index, 2>{NLayers, Nscoeffs}),
                     "source_poly_coeffs has shape {:B,}, expected [{}, {}] Stokes blocks",
                     source_poly_coeffs.shape(),
                     NLayers,
                     Nscoeffs);
  ARTS_USER_ERROR_IF(static_cast<Index>(source_coordinate_scale.size()) != NLayers or
                         static_cast<Index>(source_coordinate_offset.size()) != NLayers,
                     "Source-coordinate arrays have sizes {} and {}, expected {}",
                     source_coordinate_scale.size(),
                     source_coordinate_offset.size(),
                     NLayers);
  ARTS_USER_ERROR_IF((phase_matrix.shape() != std::array<Index, 5>{2, NFourier, NLayers, NQuad, NQuad}),
                     "phase_matrix has shape {:B,}, expected [2, {}, {}, {}, {}] Mueller blocks",
                     phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     NQuad,
                     NQuad);
  ARTS_USER_ERROR_IF((boundary_up.shape() != std::array<Index, 3>{2, NFourier, N}),
                     "boundary_up has shape {:B,}, expected [2, {}, {}] Stokes blocks",
                     boundary_up.shape(),
                     NFourier,
                     N);
  ARTS_USER_ERROR_IF((boundary_down.shape() != std::array<Index, 3>{2, NFourier, N}),
                     "boundary_down has shape {:B,}, expected [2, {}, {}] Stokes blocks",
                     boundary_down.shape(),
                     NFourier,
                     N);
  ARTS_USER_ERROR_IF((beam_phase_matrix.shape() != std::array<Index, 4>{2, NFourier, NLayers, NQuad}),
                     "beam_phase_matrix has shape {:B,}, expected [2, {}, {}, {}] Mueller blocks",
                     beam_phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     NQuad);
  ARTS_USER_ERROR_IF(static_cast<Index>(brdf_fourier_modes.size()) != NBDRF,
                     "There are {} BRDF modes, expected {}",
                     brdf_fourier_modes.size(),
                     NBDRF);
}

void main_data::check_input_value() const {
  ARTS_USER_ERROR_IF(NLayers <= 0, "VDISORT requires at least one layer");
  ARTS_USER_ERROR_IF(NFourier <= 0, "VDISORT requires at least one Fourier mode");
  ARTS_USER_ERROR_IF(tau_arr.front() <= 0.0, "tau_arr must be strictly positive, got {:B,}", tau_arr);
  ARTS_USER_ERROR_IF(std::ranges::any_of(omega_arr, [](const Numeric x) { return x < 0.0 or x > 1.0; }),
                     "omega_arr must be in [0, 1], got {:B,}",
                     omega_arr);
  ARTS_USER_ERROR_IF(mu0 < 0.0 or mu0 > 1.0, "mu0 must be in [0, 1], got {}", mu0);
  ARTS_USER_ERROR_IF(phi0 < 0.0 or phi0 >= Constant::two_pi, "phi0 must be in [0, 2*pi), got {}", phi0);
  ARTS_USER_ERROR_IF(beam_stokes[0] < 0.0, "The beam I component must be non-negative, got {}", beam_stokes[0]);
  const Numeric polarized_norm = std::hypot(beam_stokes[1], beam_stokes[2], beam_stokes[3]);
  ARTS_USER_ERROR_IF(polarized_norm > beam_stokes[0] * (1.0 + 1e-12),
                     "The beam Stokes vector is non-physical: sqrt(Q^2+U^2+V^2)={} > I={}",
                     polarized_norm,
                     beam_stokes[0]);
  ARTS_USER_ERROR_IF(has_beam_source and mu0 == 0.0, "A direct beam requires mu0 > 0");
  ARTS_USER_ERROR_IF(NBDRF > NFourier, "There are {} BRDF modes but only {} Fourier modes", NBDRF, NFourier);
}

Index main_data::tau_index(const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0.0 or tau > tau_arr.back(), "tau ({}) must be in [0, {}]", tau, tau_arr.back());
  const Index layer = std::distance(tau_arr.begin(), std::ranges::lower_bound(tau_arr, tau));
  return std::min(layer, NLayers - 1);
}

void main_data::diagonalize() {
  ARTS_TIME_REPORT

  // VDISORT CHANGE BEGIN: solve the real, generally non-symmetric 8N x 8N
  // eigenproblem with complex arithmetic (paper Eqs. 85-90).
  Matrix                       A_real(NState, NState);
  ComplexMatrix                A(NState, NState);
  ComplexMatrix                eigenvectors(NState, NState);
  ComplexVector                eigenvalues(NState);
  Vector                       rhs(NState);
  std::vector<Index>           order(static_cast<std::size_t>(NState));
  complex_diagonalize_workdata eigen_work(NState);
  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      for (Index l = 0; l < NLayers; ++l) {
        A_real = 0.0;
        for (Index i = 0; i < NQuad; ++i) {
          for (Index so = 0; so < stokes_dimension; ++so) {
            const Index row = state_index(i, so);
            for (Index j = 0; j < NQuad; ++j) {
              const Numeric factor = -0.5 * omega_arr[l] * W[j % N] * inv_mu_arr[i];
              for (Index si = 0; si < stokes_dimension; ++si) {
                A_real[row, state_index(j, si)] += factor * phase_matrix[alpha, m, l, i, j][so, si];
              }
            }
            A_real[row, row] += inv_mu_arr[i];
          }
        }

        for (Index i = 0; i < NState; ++i)
          for (Index j = 0; j < NState; ++j) A[i, j] = A_real[i, j];

        ::diagonalize(eigenvectors, eigenvalues, A, eigen_work);

        std::iota(order.begin(), order.end(), Index{0});
        std::ranges::sort(order, [&](const Index a, const Index b) {
          if (eigenvalues[a].real() != eigenvalues[b].real()) return eigenvalues[a].real() < eigenvalues[b].real();
          return eigenvalues[a].imag() < eigenvalues[b].imag();
        });

        for (Index e = 0; e < NState; ++e) {
          const Index old                            = order[static_cast<std::size_t>(e)];
          K_collect[alpha, m, l, e]                  = eigenvalues[old];
          top_anchored[anchor_index(alpha, m, l, e)] = static_cast<unsigned char>(e < NHalfState);
          for (Index s = 0; s < NState; ++s) G_collect[alpha, m, l, s, e] = eigenvectors[s, old];
        }

        if (has_beam_source) {
          const rtepack::stokvec& beam = beam_stokes;
          rhs                          = 0.0;
          const Numeric epsilon        = m == 0 ? 1.0 : 2.0;
          for (Index i = 0; i < NQuad; ++i) {
            const rtepack::stokvec scattered = beam_phase_matrix[alpha, m, l, i] * beam;
            for (Index so = 0; so < stokes_dimension; ++so) {
              rhs[state_index(i, so)] = inv_mu_arr[i] * epsilon * omega_arr[l] * scattered[so] / (4.0 * Constant::pi);
            }
          }
          diagonal(A_real) += 1.0 / mu0;
          solve_inplace(rhs, A_real);
          for (Index i = 0; i < NQuad; ++i)
            for (Index s = 0; s < stokes_dimension; ++s) B_collect[alpha, m, l, i][s] = rhs[state_index(i, s)];
        }
      }
    }
  }
  // VDISORT CHANGE END
}

void main_data::set_scales() {
  ARTS_TIME_REPORT

  scaled_source_poly_coeffs = rtepack::stokvec{};
  for (Index l = 0; l < NLayers; ++l)
    for (Index p = 0; p < Nscoeffs; ++p)
      scaled_source_poly_coeffs[l, p] = (1.0 - omega_arr[l]) * source_poly_coeffs[l, p];
}

void main_data::source_function() {
  ARTS_TIME_REPORT
  source_collect = 0.0;
  if (not has_source_poly) return;

  // VDISORT CHANGE BEGIN: polynomial source-function coefficients are
  // four-vectors B.  The transfer-equation emission is (1 - omega) B, as in
  // scalar DISORT.  The m=0 combined systems split [I,Q] and [U,V] as in
  // paper Eq. (82).
  for (Index alpha = 0; alpha < 2; ++alpha) {
    const Index m = 0;
    for (Index l = 0; l < NLayers; ++l) {
      Matrix A(NState, NState, 0.0);
      for (Index i = 0; i < NQuad; ++i) {
        for (Index so = 0; so < stokes_dimension; ++so) {
          const Index row = state_index(i, so);
          for (Index j = 0; j < NQuad; ++j) {
            const Numeric factor = -0.5 * omega_arr[l] * W[j % N] * inv_mu_arr[i];
            for (Index si = 0; si < stokes_dimension; ++si)
              A[row, state_index(j, si)] += factor * phase_matrix[alpha, m, l, i, j][so, si];
          }
          A[row, row] += inv_mu_arr[i];
        }
      }

      Vector next(NState, 0.0);
      for (Index p = Nscoeffs - 1; p >= 0; --p) {
        Vector rhs(NState, 0.0);
        for (Index i = 0; i < NQuad; ++i) {
          for (Index s = 0; s < stokes_dimension; ++s) {
            const bool    active = alpha == cosine_mode ? s < 2 : s >= 2;
            const Numeric q      = active ? scaled_source_poly_coeffs[l, p][s] : 0.0;
            rhs[state_index(i, s)] =
                inv_mu_arr[i] * q + source_coordinate_scale[l] * static_cast<Numeric>(p + 1) * next[state_index(i, s)];
          }
        }
        Matrix work{A};
        solve_inplace(rhs, work);
        for (Index i = 0; i < NQuad; ++i)
          for (Index s = 0; s < stokes_dimension; ++s) {
            const Index state                    = state_index(i, s);
            source_collect[alpha, m, l, i, p][s] = rhs[state];
            next[state]                          = rhs[state];
          }
      }
    }
  }
  // VDISORT CHANGE END
}

void main_data::transmission() {
  // VDISORT CHANGE BEGIN: complex modes are anchored at the nearest layer
  // boundary and evaluated directly by homogeneous().  This is the stable
  // complex counterpart of DISORT's real expK_collect table.
  // VDISORT CHANGE END
}

void main_data::solve_for_coefs() {
  ARTS_TIME_REPORT

  const Index                  equation_count = NLayers * NState;
  const Index                  bandwidth      = 3 * NHalfState - 1;
  matpack::complex_band_matrix lhs(bandwidth, bandwidth, equation_count, equation_count);
  ComplexVector                rhs(equation_count);
  Matrix                       reflection(NHalfState, NHalfState);
  Vector                       direct_reflection(NHalfState);

  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      lhs.zero();
      rhs = 0.0;

      // VDISORT CHANGE BEGIN: reflection is a matrix coupling Stokes components.
      reflection        = 0.0;
      direct_reflection = 0.0;
      if (m < NBDRF) {
        rtepack::muelmat_matrix raw(N, N, rtepack::muelmat{0.0});
        brdf_fourier_modes[m](alpha, raw, mu_arr[Range{0, N}], mu_arr[Range{0, N}]);
        for (Index i = 0; i < N; ++i)
          for (Index so = 0; so < stokes_dimension; ++so)
            for (Index j = 0; j < N; ++j)
              for (Index si = 0; si < stokes_dimension; ++si)
                reflection[state_index(i, so), state_index(j, si)] =
                    Constant::pi * (m == 0 ? 1.0 : 0.5) * W[j] * mu_arr[j] * raw[i, j][so, si];

        if (has_beam_source) {
          rtepack::muelmat_matrix beam_raw(N, 1, rtepack::muelmat{0.0});
          const Vector            beam_mu{mu0};
          brdf_fourier_modes[m](alpha, beam_raw, mu_arr[Range{0, N}], beam_mu);
          const Numeric           attenuation = std::exp(-tau_arr.back() / mu0);
          const rtepack::stokvec& beam        = beam_stokes;
          for (Index i = 0; i < N; ++i) {
            const rtepack::stokvec reflected = 0.5 * mu0 * attenuation * (beam_raw[i, 0] * beam);
            for (Index s = 0; s < stokes_dimension; ++s) direct_reflection[state_index(i, s)] = reflected[s];
          }
        }
      }
      // VDISORT CHANGE END

      // Top boundary: prescribed downward diffuse Stokes field.
      for (Index i = 0; i < N; ++i) {
        for (Index s = 0; s < stokes_dimension; ++s) {
          const Index row   = state_index(i, s);
          const Index state = state_index(N + i, s);
          rhs[row]          = boundary_down[alpha, m, i][s] - particular(alpha, m, 0, state, 0.0);
          for (Index e = 0; e < NState; ++e) lhs[row, e] = homogeneous(alpha, m, 0, state, e, 0.0);
        }
      }

      // Layer-interface continuity, paper Eqs. (118)-(121).
      for (Index l = 0; l < NLayers - 1; ++l) {
        const Numeric tau  = tau_arr[l];
        const Index   row0 = NHalfState + l * NState;
        for (Index state = 0; state < NState; ++state) {
          const Index row = row0 + state;
          rhs[row]        = particular(alpha, m, l + 1, state, tau) - particular(alpha, m, l, state, tau);
          for (Index e = 0; e < NState; ++e) {
            lhs[row, l * NState + e]       = homogeneous(alpha, m, l, state, e, tau);
            lhs[row, (l + 1) * NState + e] = -homogeneous(alpha, m, l + 1, state, e, tau);
          }
        }
      }

      // Lower boundary, including the polarized BRDF of paper Eq. (124).
      const Index   last = NLayers - 1;
      const Numeric tau  = tau_arr.back();
      for (Index i = 0; i < N; ++i) {
        for (Index s = 0; s < stokes_dimension; ++s) {
          const Index hrow      = state_index(i, s);
          const Index row       = equation_count - NHalfState + hrow;
          const Index pos_state = state_index(i, s);
          Numeric     rhs_value =
              boundary_up[alpha, m, i][s] + direct_reflection[hrow] - particular(alpha, m, last, pos_state, tau);
          for (Index j = 0; j < N; ++j)
            for (Index si = 0; si < stokes_dimension; ++si) {
              const Index hin        = state_index(j, si);
              const Index neg_state  = state_index(N + j, si);
              rhs_value             += reflection[hrow, hin] * particular(alpha, m, last, neg_state, tau);
            }
          rhs[row] = rhs_value;

          for (Index e = 0; e < NState; ++e) {
            Complex value = homogeneous(alpha, m, last, pos_state, e, tau);
            for (Index j = 0; j < N; ++j)
              for (Index si = 0; si < stokes_dimension; ++si)
                value -=
                    reflection[hrow, state_index(j, si)] * homogeneous(alpha, m, last, state_index(N + j, si), e, tau);
            lhs[row, last * NState + e] = value;
          }
        }
      }

      const int info = lhs.solve(rhs);
      ARTS_USER_ERROR_IF(
          info != 0, "VDISORT boundary solve failed for alpha={}, Fourier mode={} (LAPACK info={})", alpha, m, info);
      for (Index l = 0; l < NLayers; ++l)
        for (Index e = 0; e < NState; ++e) GC_collect[alpha, m, l, e] = rhs[l * NState + e];
    }
  }
}

void main_data::rad_field() {
  // VDISORT CHANGE BEGIN: retain both combined modes at every layer bottom.
  // The physical [I,Q,U,V] field is reconstructed by u() using paper Eq. (78).
  Matrix field(2 * NFourier, NState);
  for (Index l = 0; l < NLayers; ++l) {
    combined_field(field, tau_arr[l]);
    for (Index alpha = 0; alpha < 2; ++alpha)
      for (Index m = 0; m < NFourier; ++m)
        for (Index i = 0; i < NQuad; ++i)
          for (Index s = 0; s < stokes_dimension; ++s)
            um[l, alpha, m, i][s] = field[alpha * NFourier + m, state_index(i, s)];
  }
  // VDISORT CHANGE END
}

void main_data::update_all() {
  ARTS_TIME_REPORT
  has_source_poly = Nscoeffs > 0;
  has_beam_source = polarized_source(beam_stokes);
  check_input_size();
  check_input_value();
  set_scales();
  diagonalize();
  transmission();
  source_function();
  solve_for_coefs();
  rad_field();
}

void main_data::combined_field(MatrixView out, const Numeric tau) const {
  const Index layer = tau_index(tau);
  ARTS_USER_ERROR_IF((out.shape() != std::array<Index, 2>{2 * NFourier, NState}),
                     "Combined field has shape {:B,}, expected [{}, {}]",
                     out.shape(),
                     2 * NFourier,
                     NState);

  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      ComplexVector modal_amplitude(NState);
      for (Index e = 0; e < NState; ++e) {
        const Numeric anchor = top_anchored[anchor_index(alpha, m, layer, e)] ? layer_top(layer) : tau_arr[layer];
        modal_amplitude[e] = GC_collect[alpha, m, layer, e] * std::exp(K_collect[alpha, m, layer, e] * (tau - anchor));
      }
      for (Index state = 0; state < NState; ++state) {
        Complex value = particular(alpha, m, layer, state, tau);
        for (Index e = 0; e < NState; ++e) value += G_collect[alpha, m, layer, state, e] * modal_amplitude[e];
        const Numeric scale = 1.0 + std::abs(value.real());
        // VDISORT CHANGE: highly conservative scalar-limit cases can be very
        // ill-conditioned; conjugate eigenpairs may leave roundoff at a few
        // parts in 1e9 while the reconstructed physical field remains real.
        ARTS_USER_ERROR_IF(std::abs(value.imag()) > 2e-8 * scale,
                           "VDISORT left an uncancelled imaginary radiance {} at alpha={}, m={}, state={}, tau={}",
                           value,
                           alpha,
                           m,
                           state,
                           tau);
        out[alpha * NFourier + m, state] = value.real();
      }
    }
  }
}

void main_data::u(u_data& data, const Numeric tau, const Numeric phi) const {
  ARTS_USER_ERROR_IF(
      !std::isfinite(phi) or phi < 0.0 or phi >= Constant::two_pi, "phi must be finite and in [0, 2*pi), got {}", phi);
  Matrix combined(2 * NFourier, NState);
  combined_field(combined, tau);

  data.intensities.resize(NQuad);
  data.intensities = rtepack::stokvec{};
  for (Index m = 0; m < NFourier; ++m) {
    const Numeric c = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const Numeric s = std::sin(static_cast<Numeric>(m) * (phi0 - phi));
    for (Index i = 0; i < NQuad; ++i) {
      for (Index stokes = 0; stokes < 2; ++stokes)
        data.intensities[i][stokes] +=
            combined[m, state_index(i, stokes)] * c + combined[NFourier + m, state_index(i, stokes)] * s;
      // VDISORT CHANGE: U and V swap combined cosine/sine roles, paper Eq. (78).
      for (Index stokes = 2; stokes < 4; ++stokes)
        data.intensities[i][stokes] +=
            combined[NFourier + m, state_index(i, stokes)] * c + combined[m, state_index(i, stokes)] * s;
    }
  }
}

void main_data::u_user(user_u_data&                  data,
                       const Numeric                 tau,
                       const Numeric                 phi,
                       const ConstVectorView&        user_mu,
                       const phase_matrix_data&      user_phase_matrix,
                       const beam_phase_matrix_data& user_beam_phase_matrix) const {
  rtepack::stokvec_tensor3 result(1, 1, user_mu.size());
  ungridded_u_user(result, AscendingGrid{tau}, Vector{phi}, user_mu, user_phase_matrix, user_beam_phase_matrix);
  data.intensities.resize(user_mu.size());
  for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user) data.intensities[user] = result[0, 0, user];
}

void main_data::user_fourier_modes(ComplexTensor4&               modes,
                                   const AscendingGrid&          tau,
                                   const ConstVectorView&        user_mu,
                                   const phase_matrix_data&      user_phase_matrix,
                                   const beam_phase_matrix_data& user_beam_phase_matrix) const {
  ARTS_TIME_REPORT

  const Index ntau  = static_cast<Index>(tau.size());
  const Index nuser = static_cast<Index>(user_mu.size());
  ARTS_USER_ERROR_IF(tau.empty(), "At least one optical depth is required for user-angle radiances");
  ARTS_USER_ERROR_IF(tau.front() < 0.0 or tau.back() > tau_arr.back(),
                     "User-angle optical depths must be in [0, {}], got {:B,}",
                     tau_arr.back(),
                     tau);
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(user_mu,
                          [](const Numeric mu) { return !std::isfinite(mu) or mu == 0.0 or std::abs(mu) > 1.0; }),
      "User polar-angle cosines must be finite, nonzero, and in [-1, 1], got {:B,}",
      user_mu);
  const std::array<Index, 5> expected_phase_shape{2, NFourier, NLayers, nuser, NQuad};
  const std::array<Index, 4> expected_beam_shape{2, NFourier, NLayers, nuser};
  ARTS_USER_ERROR_IF(user_phase_matrix.shape() != expected_phase_shape,
                     "User phase matrices have shape {:B,}, expected [2, {}, {}, {}, {}] Mueller blocks",
                     user_phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     nuser,
                     NQuad);
  ARTS_USER_ERROR_IF(has_beam_source and user_beam_phase_matrix.shape() != expected_beam_shape,
                     "User beam phase matrices have shape {:B,}, expected [2, {}, {}, {}] Mueller blocks",
                     user_beam_phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     nuser);

  Vector     interpolation_values(N);
  const auto interpolate_boundary =
      [&](const rtepack::stokvec_tensor3& boundary, const Index alpha, const Index m, const Numeric abs_mu) {
        rtepack::stokvec result{};
        for (Index s = 0; s < stokes_dimension; ++s) {
          for (Index i = 0; i < N; ++i) interpolation_values[i] = boundary[alpha, m, i][s];
          result[s] = barycentric_interpolate(
              mu_arr[Range{0, N}], half_range_barycentric_weights, interpolation_values, abs_mu);
        }
        return result;
      };

  modes.resize(ntau, nuser, 2 * NFourier, stokes_dimension);
  modes = Complex{0.0};

  rtepack::muelmat_matrix raw(1, N);
  rtepack::muelmat_matrix beam_raw(1, 1);
  Vector                  outgoing(1);
  const Vector            beam_direction{mu0};
  const Numeric           atmosphere_bottom = tau_arr.back();

  for (Index iu = 0; iu < nuser; ++iu) {
    const Numeric mu       = user_mu[iu];
    const Numeric abs_mu   = std::abs(mu);
    const bool    downward = mu < 0.0;
    outgoing[0]            = abs_mu;

    for (Index alpha = 0; alpha < 2; ++alpha) {
      for (Index m = 0; m < NFourier; ++m) {
        rtepack::stokvec mode = interpolate_boundary(downward ? boundary_down : boundary_up, alpha, m, abs_mu);

        const Numeric boundary_tau = downward ? 0.0 : tau_arr.back();
        if (not downward and m < NBDRF) {
          brdf_fourier_modes[m](alpha, raw, outgoing, mu_arr[Range{0, N}]);
          for (Index j = 0; j < N; ++j) {
            mode +=
                Constant::pi * (m == 0 ? 1.0 : 0.5) * W[j] * mu_arr[j] * (raw[0, j] * um[NLayers - 1, alpha, m, N + j]);
          }
          if (has_beam_source) {
            brdf_fourier_modes[m](alpha, beam_raw, outgoing, beam_direction);
            mode += 0.5 * mu0 * std::exp(-atmosphere_bottom / mu0) * (beam_raw[0, 0] * beam_stokes);
          }
        }
        for (Index t = 0; t < ntau; ++t) {
          const Numeric attenuation = std::exp(-std::abs(tau[t] - boundary_tau) / abs_mu);
          for (Index s = 0; s < stokes_dimension; ++s) modes[t, iu, alpha * NFourier + m, s] += attenuation * mode[s];
        }
      }
    }
  }

  const auto limits = [&](const Index layer, const Numeric observation, const bool downward) {
    const Numeric top    = layer_top(layer);
    const Numeric bottom = tau_arr[layer];
    return std::pair{downward ? top : std::max(top, observation), downward ? std::min(bottom, observation) : bottom};
  };

  rtepack::stokvec_vector polynomial_coefficients(Nscoeffs);
  for (Index layer = 0; layer < NLayers; ++layer) {
    for (Index alpha = 0; alpha < 2; ++alpha) {
      for (Index m = 0; m < NFourier; ++m) {
        const Index mode = alpha * NFourier + m;
        for (Index iu = 0; iu < nuser; ++iu) {
          const Numeric mu       = user_mu[iu];
          const Numeric abs_mu   = std::abs(mu);
          const bool    downward = mu < 0.0;

          for (Index eigen = 0; eigen < NState; ++eigen) {
            std::array<Complex, stokes_dimension> source{};
            for (Index j = 0; j < NQuad; ++j) {
              const Numeric factor = 0.5 * omega_arr[layer] * W[j % N];
              const auto&   phase  = user_phase_matrix[alpha, m, layer, iu, j];
              for (Index so = 0; so < stokes_dimension; ++so)
                for (Index si = 0; si < stokes_dimension; ++si)
                  source[so] += factor * phase[so, si] * G_collect[alpha, m, layer, state_index(j, si), eigen];
            }
            const Complex coefficient = GC_collect[alpha, m, layer, eigen];
            if (coefficient == Complex{0.0}) continue;
            const Numeric anchor =
                top_anchored[anchor_index(alpha, m, layer, eigen)] ? layer_top(layer) : tau_arr[layer];
            for (Index t = 0; t < ntau; ++t) {
              const auto [lower, upper] = limits(layer, tau[t], downward);
              if (upper <= lower) continue;
              const Complex integral = user_angle_exponential_integral(
                  K_collect[alpha, m, layer, eigen], anchor, lower, upper, tau[t], abs_mu, downward);
              for (Index s = 0; s < stokes_dimension; ++s) modes[t, iu, mode, s] += coefficient * source[s] * integral;
            }
          }

          if (has_beam_source) {
            rtepack::stokvec source{};
            for (Index j = 0; j < NQuad; ++j)
              source += 0.5 * omega_arr[layer] * W[j % N] *
                        (user_phase_matrix[alpha, m, layer, iu, j] * B_collect[alpha, m, layer, j]);
            const Numeric epsilon  = m == 0 ? 1.0 : 2.0;
            source                += epsilon * omega_arr[layer] / (4.0 * Constant::pi) *
                                     (user_beam_phase_matrix[alpha, m, layer, iu] * beam_stokes);
            for (Index t = 0; t < ntau; ++t) {
              const auto [lower, upper] = limits(layer, tau[t], downward);
              if (upper <= lower) continue;
              const Complex integral =
                  user_angle_exponential_integral(-1.0 / mu0, 0.0, lower, upper, tau[t], abs_mu, downward);
              for (Index s = 0; s < stokes_dimension; ++s) modes[t, iu, mode, s] += source[s] * integral;
            }
          }

          if (has_source_poly) {
            polynomial_coefficients = rtepack::stokvec{};
            for (Index p = 0; p < Nscoeffs; ++p) {
              for (Index j = 0; j < NQuad; ++j)
                polynomial_coefficients[p] +=
                    0.5 * omega_arr[layer] * W[j % N] *
                    (user_phase_matrix[alpha, m, layer, iu, j] * source_collect[alpha, m, layer, j, p]);
              if (m == 0)
                for (Index s = 0; s < stokes_dimension; ++s) {
                  const bool active = alpha == cosine_mode ? s < 2 : s >= 2;
                  if (active) polynomial_coefficients[p][s] += scaled_source_poly_coeffs[layer, p][s];
                }
            }
            for (Index t = 0; t < ntau; ++t) {
              const auto [lower, upper] = limits(layer, tau[t], downward);
              if (upper <= lower) continue;
              const rtepack::stokvec integral = user_angle_polynomial_integral(polynomial_coefficients,
                                                                               source_coordinate_scale[layer],
                                                                               source_coordinate_offset[layer],
                                                                               lower,
                                                                               upper,
                                                                               tau[t],
                                                                               abs_mu,
                                                                               downward);
              for (Index s = 0; s < stokes_dimension; ++s) modes[t, iu, mode, s] += integral[s];
            }
          }
        }
      }
    }
  }
}

void main_data::ungridded_u_user(rtepack::stokvec_tensor3_view out,
                                 const AscendingGrid&          tau,
                                 const Vector&                 phi,
                                 const ConstVectorView&        user_mu,
                                 const phase_matrix_data&      user_phase_matrix,
                                 const beam_phase_matrix_data& user_beam_phase_matrix) const {
  ARTS_TIME_REPORT

  const std::array<Index, 3> expected_shape{
      static_cast<Index>(tau.size()), static_cast<Index>(phi.size()), static_cast<Index>(user_mu.size())};
  ARTS_USER_ERROR_IF(out.shape() != expected_shape,
                     "User-angle output has shape {:B,}, expected [{}, {}, {}] Stokes blocks",
                     out.shape(),
                     tau.size(),
                     phi.size(),
                     user_mu.size());
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          phi, [](const Numeric value) { return !std::isfinite(value) or value < 0.0 or value >= Constant::two_pi; }),
      "User azimuths must be finite and in [0, 2*pi), got {:B,}",
      phi);

  ComplexTensor4 modes;
  user_fourier_modes(modes, tau, user_mu, user_phase_matrix, user_beam_phase_matrix);
  out = rtepack::stokvec{};
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t) {
    for (Index iu = 0; iu < static_cast<Index>(user_mu.size()); ++iu) {
      for (Index m = 0; m < NFourier; ++m) {
        for (Index alpha = 0; alpha < 2; ++alpha)
          for (Index stokes = 0; stokes < stokes_dimension; ++stokes) {
            const Complex value = modes[t, iu, alpha * NFourier + m, stokes];
            const Numeric scale = 1.0 + std::abs(value.real());
            ARTS_USER_ERROR_IF(std::abs(value.imag()) > 2e-8 * scale,
                               "VDISORT left an uncancelled imaginary user radiance {} at tau={}, user={}, "
                               "alpha={}, m={}, stokes={}",
                               value,
                               tau[t],
                               iu,
                               alpha,
                               m,
                               stokes);
          }

        for (Index p = 0; p < static_cast<Index>(phi.size()); ++p) {
          const Numeric c = std::cos(static_cast<Numeric>(m) * (phi0 - phi[p]));
          const Numeric s = std::sin(static_cast<Numeric>(m) * (phi0 - phi[p]));
          for (Index stokes = 0; stokes < 2; ++stokes)
            out[t, p, iu][stokes] += modes[t, iu, m, stokes].real() * c + modes[t, iu, NFourier + m, stokes].real() * s;
          for (Index stokes = 2; stokes < stokes_dimension; ++stokes)
            out[t, p, iu][stokes] += modes[t, iu, NFourier + m, stokes].real() * c + modes[t, iu, m, stokes].real() * s;
        }
      }
    }
  }
}

void main_data::u0(u0_data& data, const Numeric tau) const {
  Matrix combined(2 * NFourier, NState);
  combined_field(combined, tau);
  data.u0.resize(NQuad);
  for (Index i = 0; i < NQuad; ++i) {
    data.u0[i][0] = combined[0, state_index(i, 0)];
    data.u0[i][1] = combined[0, state_index(i, 1)];
    data.u0[i][2] = combined[NFourier, state_index(i, 2)];
    data.u0[i][3] = combined[NFourier, state_index(i, 3)];
  }
}

Numeric main_data::flux_up(flux_data& data, const Numeric tau) const {
  u0_data field;
  u0(field, tau);
  data.u0        = std::move(field.u0);
  Numeric result = 0.0;
  for (Index i = 0; i < N; ++i) result += W[i] * mu_arr[i] * data.u0[i].I();
  return Constant::two_pi * result;
}

std::pair<Numeric, Numeric> main_data::flux_down(flux_data& data, const Numeric tau) const {
  u0_data field;
  u0(field, tau);
  data.u0         = std::move(field.u0);
  Numeric diffuse = 0.0;
  for (Index i = 0; i < N; ++i) diffuse += W[i] * mu_arr[i] * data.u0[N + i].I();
  const Numeric direct = has_beam_source ? mu0 * beam_stokes.I() * std::exp(-tau / mu0) : 0.0;
  return {Constant::two_pi * diffuse, direct};
}

flux_values main_data::flux(flux_data& data, const Numeric tau) const {
  u0_data field;
  u0(field, tau);
  data.u0 = std::move(field.u0);

  Numeric up = 0.0, down = 0.0, mean_intensity = 0.0;
  for (Index i = 0; i < N; ++i) {
    const Numeric upward    = data.u0[i].I();
    const Numeric downward  = data.u0[N + i].I();
    up                     += W[i] * mu_arr[i] * upward;
    down                   += W[i] * mu_arr[i] * downward;
    mean_intensity         += 0.5 * W[i] * (upward + downward);
  }

  const Numeric direct = has_beam_source ? mu0 * beam_stokes.I() * std::exp(-tau / mu0) : 0.0;
  if (has_beam_source) mean_intensity += beam_stokes.I() * std::exp(-tau / mu0) / (4.0 * Constant::pi);

  const Index   layer      = tau_index(tau);
  const Numeric source_tau = std::fma(source_coordinate_scale[layer], tau, source_coordinate_offset[layer]);
  Numeric       source     = 0.0;
  for (Index coefficient = Nscoeffs - 1; coefficient >= 0; --coefficient)
    source = std::fma(source, source_tau, source_poly_coeffs[layer, coefficient].I());

  return {.up           = Constant::two_pi * up,
          .down_diffuse = Constant::two_pi * down,
          .down_direct  = direct,
          .dfdt         = 4.0 * Constant::pi * (1.0 - omega_arr[layer]) * (mean_intensity - source)};
}

void main_data::gridded_u(Tensor4View out, const Vector& phi) const {
  ARTS_USER_ERROR_IF(
      (out.shape() != std::array<Index, 4>{NLayers, static_cast<Index>(phi.size()), NQuad, stokes_dimension}),
      "gridded_u output has shape {:B,}, expected [{}, {}, {}, 4]",
      out.shape(),
      NLayers,
      phi.size(),
      NQuad);

  const Index nphi = static_cast<Index>(phi.size());
  Matrix      cosine(nphi, NFourier);
  Matrix      sine(nphi, NFourier);
  for (Index p = 0; p < nphi; ++p) {
    ARTS_USER_ERROR_IF(!std::isfinite(phi[p]) or phi[p] < 0.0 or phi[p] >= Constant::two_pi,
                       "phi must be finite and in [0, 2*pi), got {}",
                       phi[p]);
    for (Index m = 0; m < NFourier; ++m) {
      const Numeric angle = static_cast<Numeric>(m) * (phi0 - phi[p]);
      cosine[p, m]        = std::cos(angle);
      sine[p, m]          = std::sin(angle);
    }
  }

  out = 0.0;
  for (Index l = 0; l < NLayers; ++l)
    for (Index p = 0; p < nphi; ++p)
      for (Index m = 0; m < NFourier; ++m) {
        const Numeric c = cosine[p, m];
        const Numeric s = sine[p, m];
        for (Index i = 0; i < NQuad; ++i) {
          for (Index stokes = 0; stokes < 2; ++stokes)
            out[l, p, i, stokes] += um[l, cosine_mode, m, i][stokes] * c + um[l, sine_mode, m, i][stokes] * s;
          // VDISORT paper Eq. (78): U and V swap combined cosine/sine roles.
          for (Index stokes = 2; stokes < stokes_dimension; ++stokes)
            out[l, p, i, stokes] += um[l, sine_mode, m, i][stokes] * c + um[l, cosine_mode, m, i][stokes] * s;
        }
      }
}

void main_data::gridded_flux(VectorView up, VectorView down, VectorView down_direct) const {
  ARTS_USER_ERROR_IF(up.size() != static_cast<Size>(NLayers) or down.size() != static_cast<Size>(NLayers) or
                         down_direct.size() != static_cast<Size>(NLayers),
                     "All gridded flux outputs must have size {}",
                     NLayers);
  for (Index l = 0; l < NLayers; ++l) {
    Numeric upward   = 0.0;
    Numeric downward = 0.0;
    for (Index i = 0; i < N; ++i) {
      upward   += W[i] * mu_arr[i] * um[l, cosine_mode, 0, i].I();
      downward += W[i] * mu_arr[i] * um[l, cosine_mode, 0, N + i].I();
    }
    up[l]          = Constant::two_pi * upward;
    down[l]        = Constant::two_pi * downward;
    down_direct[l] = has_beam_source ? mu0 * beam_stokes.I() * std::exp(-tau_arr[l] / mu0) : 0.0;
  }
}

void main_data::ungridded_u(Tensor4View out, const AscendingGrid& tau, const Vector& phi) const {
  ARTS_USER_ERROR_IF(
      (out.shape() !=
       std::array<Index, 4>{static_cast<Index>(tau.size()), static_cast<Index>(phi.size()), NQuad, stokes_dimension}),
      "ungridded_u output has shape {:B,}, expected [{}, {}, {}, 4]",
      out.shape(),
      tau.size(),
      phi.size(),
      NQuad);

  const Index nphi = static_cast<Index>(phi.size());
  if (tau.empty() or nphi == 0) return;

  Matrix cosine(nphi, NFourier);
  Matrix sine(nphi, NFourier);
  for (Index p = 0; p < nphi; ++p) {
    ARTS_USER_ERROR_IF(!std::isfinite(phi[p]) or phi[p] < 0.0 or phi[p] >= Constant::two_pi,
                       "phi must be finite and in [0, 2*pi), got {}",
                       phi[p]);
    for (Index m = 0; m < NFourier; ++m) {
      const Numeric angle = static_cast<Numeric>(m) * (phi0 - phi[p]);
      cosine[p, m]        = std::cos(angle);
      sine[p, m]          = std::sin(angle);
    }
  }

  out = 0.0;
  Matrix combined(2 * NFourier, NState);
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t) {
    combined_field(combined, tau[t]);
    for (Index p = 0; p < nphi; ++p)
      for (Index m = 0; m < NFourier; ++m) {
        const Numeric c = cosine[p, m];
        const Numeric s = sine[p, m];
        for (Index i = 0; i < NQuad; ++i) {
          for (Index stokes = 0; stokes < 2; ++stokes)
            out[t, p, i, stokes] +=
                combined[m, state_index(i, stokes)] * c + combined[NFourier + m, state_index(i, stokes)] * s;
          for (Index stokes = 2; stokes < stokes_dimension; ++stokes)
            out[t, p, i, stokes] +=
                combined[NFourier + m, state_index(i, stokes)] * c + combined[m, state_index(i, stokes)] * s;
        }
      }
  }
}

void main_data::ungridded_flux(VectorView up, VectorView down, VectorView down_direct, const AscendingGrid& tau) const {
  ARTS_USER_ERROR_IF(up.size() != tau.size() or down.size() != tau.size() or down_direct.size() != tau.size(),
                     "All ungridded flux outputs must have the same size as tau ({})",
                     tau.size());
  Matrix combined(2 * NFourier, NState);
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t) {
    combined_field(combined, tau[t]);
    Numeric upward   = 0.0;
    Numeric downward = 0.0;
    for (Index i = 0; i < N; ++i) {
      upward   += W[i] * mu_arr[i] * combined[0, state_index(i, 0)];
      downward += W[i] * mu_arr[i] * combined[0, state_index(N + i, 0)];
    }
    up[t]          = Constant::two_pi * upward;
    down[t]        = Constant::two_pi * downward;
    down_direct[t] = has_beam_source ? mu0 * beam_stokes.I() * std::exp(-tau[t] / mu0) : 0.0;
  }
}

rtepack::stokvec_tensor3_const_view main_data::layer_um(const Size layer) const {
  ARTS_USER_ERROR_IF(layer >= static_cast<Size>(NLayers), "Layer {} is out of range [0, {})", layer, NLayers);
  return um[layer];
}

void main_data::set_beam_source(rtepack::stokvec beam) {
  beam_stokes     = std::move(beam);
  has_beam_source = polarized_source(beam_stokes);
}

bool main_data::has_complex_eigensolutions(const Numeric tolerance) const {
  for (Index alpha = 0; alpha < 2; ++alpha)
    for (Index m = 0; m < NFourier; ++m)
      for (Index l = 0; l < NLayers; ++l)
        for (Index e = 0; e < NState; ++e) {
          const Complex value = K_collect[alpha, m, l, e];
          if (std::abs(value.imag()) > tolerance * (1.0 + std::abs(value.real()))) return true;
        }
  return false;
}
}  // namespace vdisort

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

namespace dc = disort_common;

namespace {
/** Return whether a beam Stokes vector carries positive incident intensity. */
[[nodiscard]] bool polarized_source(const rtepack::stokvec& stokes) { return stokes.I() > 0.0; }

/** Compute integral moments of a finite exponential segment by series or recurrence. */
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

/** Integrate a Stokes-valued source polynomial along one finite user ray. */
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

/** Form combined cosine/sine Mueller blocks from ordinary Fourier coefficients. */
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

/** Verify that decaying and growing eigenmodes can be assigned to opposite boundaries.
 *
 * Eigenvalues whose real part is negligible relative to the spectral scale
 * are neutral and may occupy either half.  Every non-neutral mode must have
 * the sign implied by its sorted position.
 */
void check_eigenvalue_direction_split(const ComplexVector&      eigenvalues,
                                      const std::vector<Index>& order,
                                      const Index               half_state,
                                      const Index               alpha,
                                      const Index               mode,
                                      const Index               layer) {
  Numeric scale = 1.0;
  for (const Complex value : eigenvalues) {
    if (not std::isfinite(value.real()) or not std::isfinite(value.imag()))
      arts::user_error(
          "VDISORT eigensystem ({}, {}, {}) contains the non-finite eigenvalue {}", alpha, mode, layer, value);
    scale = std::max(scale, std::abs(value));
  }
  constexpr Numeric relative_tolerance = 1.0e-10;
  const Numeric     tolerance          = relative_tolerance * scale;

  Index negative = 0, neutral = 0, positive = 0;
  for (const Complex value : eigenvalues) {
    if (value.real() < -tolerance)
      ++negative;
    else if (value.real() > tolerance)
      ++positive;
    else
      ++neutral;
  }
  if (negative > half_state or positive > half_state)
    arts::user_error(
        "VDISORT eigensystem ({}, {}, {}) cannot be split between the boundaries: "
        "{} negative, {} neutral, and {} positive modes for half-state size {} "
        "(neutral tolerance {})",
        alpha,
        mode,
        layer,
        negative,
        neutral,
        positive,
        half_state,
        tolerance);

  for (Index sorted = 0; sorted < 2 * half_state; ++sorted) {
    const Numeric real = eigenvalues[order[static_cast<std::size_t>(sorted)]].real();
    if ((sorted < half_state and real > tolerance) or (sorted >= half_state and real < -tolerance))
      arts::user_error(
          "VDISORT eigenvalue {} at sorted position {} has the wrong propagation direction "
          "for eigensystem ({}, {}, {}); neutral tolerance is {}",
          eigenvalues[order[static_cast<std::size_t>(sorted)]],
          sorted,
          alpha,
          mode,
          layer,
          tolerance);
  }
}

}  // namespace

namespace vdisort {
delta_m_correction_cache::delta_m_correction_cache(AscendingGrid                       physical_tau,
                                                   Vector                              omega,
                                                   Vector                              fraction,
                                                   const Numeric                       mu0,
                                                   const Numeric                       phi0,
                                                   rtepack::stokvec                    beam,
                                                   Vector                              user_mu,
                                                   Vector                              phi,
                                                   lab_phase_function                  original_phase,
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
  if (omega_.size() != physical_tau_.size() or fraction_.size() != physical_tau_.size())
    arts::user_error("Polarized delta-M arrays must all contain {} layers, got omega={} and fraction={}",
                     nlayers,
                     omega_.size(),
                     fraction_.size());
  if (mu0_ <= 0.0 or mu0_ > 1.0) arts::user_error("The delta-M beam cosine must be in (0, 1], got {}", mu0_);
  if (intermediate_mu <= 0 or intermediate_phi <= 0)
    arts::user_error("IMS quadrature sizes must be positive, got {} x {}", intermediate_mu, intermediate_phi);
  if (stdr::any_of(user_mu_, [](const Numeric mu) { return mu == 0.0 or std::abs(mu) > 1.0; }))
    arts::user_error("Delta-M user cosines must be nonzero and in [-1, 1], got {:B,}", user_mu_);

  scale_.resize(nlayers);
  Vector  scaled_tau(nlayers);
  Numeric physical_top = 0.0, scaled_top = 0.0;
  for (Index layer = 0; layer < nlayers; ++layer) {
    if (omega_[layer] < 0.0 or omega_[layer] > 1.0 or fraction_[layer] < 0.0 or fraction_[layer] >= 1.0)
      arts::user_error("Invalid delta-M omega/fraction in layer {}: {}, {}", layer, omega_[layer], fraction_[layer]);
    scale_[layer]      = 1.0 - omega_[layer] * fraction_[layer];
    scaled_top        += scale_[layer] * (physical_tau_[layer] - physical_top);
    scaled_tau[layer]  = scaled_top;
    physical_top       = physical_tau_[layer];
  }
  scaled_tau_ = AscendingGrid{std::move(scaled_tau)};

  // The original part of TMS; the transport part is the solver's (tms_source)
  tms_operator_.resize(nlayers, nphi, nuser);
  for (Index layer = 0; layer < nlayers; ++layer) {
    const Numeric transport_omega = omega_[layer] * (1.0 - fraction_[layer]) / scale_[layer];
    for (Index p = 0; p < nphi; ++p)
      for (Index u = 0; u < nuser; ++u) {
        const Numeric beam_factor  = user_mu_[u] > 0.0 ? mu0_ / (mu0_ + user_mu_[u]) : 1.0;
        tms_operator_[layer, p, u] = transport_omega * Constant::inv_pi * 0.25 * beam_factor *
                                     original_phase(layer, user_mu_[u], phi_[p], -mu0_, phi0_) /
                                     (1.0 - fraction_[layer]);
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

Numeric delta_m_correction_cache::scaled_tau(const Numeric physical_tau) const {
  if (physical_tau_.empty()) arts::user_error("Cannot evaluate an empty delta-M correction cache");
  if (physical_tau < 0.0 or physical_tau > physical_tau_.back())
    arts::user_error("Physical correction depth must be in [0, {}], got {}", physical_tau_.back(), physical_tau);
  const Index layer = std::min<Index>(
      std::distance(physical_tau_.begin(), stdr::lower_bound(physical_tau_, physical_tau)), physical_tau_.size() - 1);
  const Numeric physical_layer_top = layer == 0 ? 0.0 : physical_tau_[layer - 1];
  const Numeric scaled_layer_top   = layer == 0 ? 0.0 : scaled_tau_[layer - 1];
  return scaled_layer_top + scale_[layer] * (physical_tau - physical_layer_top);
}

void delta_m_correction_cache::check_transport_beam_phase(const beam_phase_matrix_data& transport_beam_phase) const {
  if (transport_beam_phase.empty()) return;
  if (transport_beam_phase.extent(0) != 2 or
      transport_beam_phase.extent(2) != static_cast<Index>(physical_tau_.size()) or
      transport_beam_phase.extent(3) != static_cast<Index>(user_mu_.size()))
    arts::user_error(
        "The transport beam phase matrices of TMS must be [2, NFourier, {}, {}] (layers, correction "
        "directions), got {:B,}",
        physical_tau_.size(),
        user_mu_.size(),
        transport_beam_phase.shape());
}

/* The beam's single scattering as the solver carries it is the transport
   phase matrix from the beam synthesised from its NFourier combined modes
   at the azimuth phi, applied to the beam,
     Z(mu, phi; -mu0, phi0) S_b = sum_m eps_m [C^m cos m(phi0 - phi) + S^m sin m(phi0 - phi)] S_b,
   with C^m in the rows I, Q of the cosine system and U, V of the sine
   system and S^m in the others (combine_beam_phase_matrices), times the
   factor of the original part. */
rtepack::stokvec delta_m_correction_cache::tms_source(const beam_phase_matrix_data& transport_beam_phase,
                                                      const Index                   layer,
                                                      const Index                   phi_index,
                                                      const Index                   user) const {
  rtepack::stokvec source = tms_operator_[layer, phi_index, user] * beam_;
  if (transport_beam_phase.empty()) return source;

  rtepack::muelmat Z{0.0};
  for (Index m = 0; m < transport_beam_phase.extent(1); ++m) {
    const Numeric eps = m == 0 ? 1.0 : 2.0;
    const Numeric c   = eps * std::cos(static_cast<Numeric>(m) * (phi0_ - phi_[phi_index]));
    const Numeric s   = eps * std::sin(static_cast<Numeric>(m) * (phi0_ - phi_[phi_index]));
    const auto&   Bc  = transport_beam_phase[cosine_mode, m, layer, user];
    const auto&   Bs  = transport_beam_phase[sine_mode, m, layer, user];
    for (Index a = 0; a < stokes_dimension; ++a)
      for (Index b = 0; b < stokes_dimension; ++b)
        Z[a, b] += a < 2 ? Bc[a, b] * c + Bs[a, b] * s : Bs[a, b] * c + Bc[a, b] * s;
  }
  const Numeric transport_omega = omega_[layer] * (1.0 - fraction_[layer]) / scale_[layer];
  const Numeric beam_factor     = user_mu_[user] > 0.0 ? mu0_ / (mu0_ + user_mu_[user]) : 1.0;
  return source - (transport_omega * Constant::inv_pi * 0.25 * beam_factor) * (Z * beam_);
}

void delta_m_correction_cache::TMS(rtepack::stokvec_vector&      tms,
                                   const Numeric                 tau,
                                   const Index                   phi_index,
                                   const beam_phase_matrix_data& transport_beam_phase) const {
  ARTS_TIME_REPORT

  check_transport_beam_phase(transport_beam_phase);
  const Numeric transport_tau = scaled_tau(tau);
  if (phi_index < 0 or phi_index >= static_cast<Index>(phi_.size()))
    arts::user_error("Correction azimuth index {} is outside [0, {})", phi_index, phi_.size());
  const Index   layer = std::min<Index>(std::distance(physical_tau_.begin(), stdr::lower_bound(physical_tau_, tau)),
                                        physical_tau_.size() - 1);
  const Numeric scaled_layer_top = layer == 0 ? 0.0 : scaled_tau_[layer - 1];

  tms.resize(user_mu_.size());
  tms = rtepack::stokvec{};
  for (Index u = 0; u < static_cast<Index>(user_mu_.size()); ++u) {
    const Numeric mu = user_mu_[u], abs_mu = std::abs(mu);
    if (mu > 0.0) {
      const Numeric at_observation = std::exp(-transport_tau / mu0_);
      const Numeric at_bottom  = std::exp((transport_tau - scaled_tau_[layer]) / abs_mu - scaled_tau_[layer] / mu0_);
      tms[u]                  += (at_observation - at_bottom) * tms_source(transport_beam_phase, layer, phi_index, u);
      for (Index l = layer + 1; l < static_cast<Index>(physical_tau_.size()); ++l) {
        const Numeric top = scaled_tau_[l - 1], bottom = scaled_tau_[l];
        const Numeric a  = std::exp((transport_tau - top) / abs_mu - top / mu0_);
        const Numeric b  = std::exp((transport_tau - bottom) / abs_mu - bottom / mu0_);
        tms[u]          += (a - b) * tms_source(transport_beam_phase, l, phi_index, u);
      }
    } else {
      const Numeric at_observation = std::exp(-transport_tau / mu0_);
      const Numeric at_top         = std::exp((scaled_layer_top - transport_tau) / abs_mu - scaled_layer_top / mu0_);
      tms[u] += dc::downward_tms_kernel(abs_mu, mu0_, transport_tau - scaled_layer_top, at_observation, at_top) *
                tms_source(transport_beam_phase, layer, phi_index, u);
      for (Index l = 0; l < layer; ++l) {
        const Numeric top = l == 0 ? 0.0 : scaled_tau_[l - 1], bottom = scaled_tau_[l];
        const Numeric a  = std::exp((bottom - transport_tau) / abs_mu - bottom / mu0_);
        const Numeric b  = std::exp((top - transport_tau) / abs_mu - top / mu0_);
        tms[u]          += dc::downward_tms_kernel(abs_mu, mu0_, bottom - top, a, b) *
                           tms_source(transport_beam_phase, l, phi_index, u);
      }
    }
  }
}

void delta_m_correction_cache::IMS(rtepack::stokvec_vector& ims,
                                   const Numeric            tau,
                                   const Index              phi_index,
                                   const ims_convention     convention) const {
  ARTS_TIME_REPORT

  static_cast<void>(scaled_tau(tau));
  if (phi_index < 0 or phi_index >= static_cast<Index>(phi_.size()))
    arts::user_error("Correction azimuth index {} is outside [0, {})", phi_index, phi_.size());
  const Index   layer = std::min<Index>(std::distance(physical_tau_.begin(), stdr::lower_bound(physical_tau_, tau)),
                                        physical_tau_.size() - 1);
  const Numeric physical_layer_top = layer == 0 ? 0.0 : physical_tau_[layer - 1];
  const Numeric layer_bottom       = physical_tau_[layer];
  const Numeric weight             = (tau - physical_layer_top) / (layer_bottom - physical_layer_top);
  const Numeric beam_theta         = std::acos(-mu0_);
  const Numeric sign               = convention == ims_convention::disort ? -1.0 : 1.0;

  ims.resize(user_mu_.size());
  ims = rtepack::stokvec{};
  for (Index u = 0; u < static_cast<Index>(user_mu_.size()); ++u) {
    const Numeric mu = user_mu_[u];
    if (mu > 0.0) continue;
    if (convention == ims_convention::disort and std::abs(beam_theta - std::acos(mu)) > Constant::pi / 18.0) continue;

    const Numeric abs_mu              = -mu;
    const auto    boundary_correction = [&](const Index boundary) {
      return sign * ims_scalar_[boundary] * dc::ims_chi(tau, abs_mu, ims_mu0_[boundary]) *
             (ims_operator_[boundary, phi_index, u] * beam_);
    };
    ims[u] = (1.0 - weight) * boundary_correction(layer) + weight * boundary_correction(layer + 1);
  }
}

void delta_m_correction_cache::gridded_TMS(rtepack::stokvec_tensor3_view tms,
                                           const beam_phase_matrix_data& transport_beam_phase) const {
  ARTS_TIME_REPORT

  check_transport_beam_phase(transport_beam_phase);
  const std::array<Index, 3> expected{
      static_cast<Index>(physical_tau_.size()), static_cast<Index>(phi_.size()), static_cast<Index>(user_mu_.size())};
  if (tms.shape() != expected)
    arts::user_error("gridded_TMS output has shape {:B,}, expected {:B,}", tms.shape(), expected);

  tms = rtepack::stokvec{};
  Vector beam_at_boundary(expected[0] + 1);
  beam_at_boundary[0] = 1.0;
  for (Index layer = 0; layer < expected[0]; ++layer)
    beam_at_boundary[layer + 1] = std::exp(-scaled_tau_[layer] / mu0_);

  Matrix ray_transport(expected[2], expected[0]);
  for (Index user = 0; user < expected[2]; ++user)
    for (Index layer = 0; layer < expected[0]; ++layer) {
      const Numeric top          = layer == 0 ? 0.0 : scaled_tau_[layer - 1];
      const Numeric thickness    = scaled_tau_[layer] - top;
      ray_transport[user, layer] = std::exp(-thickness / std::abs(user_mu_[user]));
    }

  for (Index p = 0; p < expected[1]; ++p)
    for (Index user = 0; user < expected[2]; ++user) {
      const Numeric abs_mu = std::abs(user_mu_[user]);
      if (user_mu_[user] < 0.0) {
        rtepack::stokvec downward{};
        for (Index layer = 0; layer < expected[0]; ++layer) {
          const Numeric top       = layer == 0 ? 0.0 : scaled_tau_[layer - 1];
          const Numeric bottom    = scaled_tau_[layer];
          const Numeric thickness = bottom - top;
          const Numeric transport = ray_transport[user, layer];
          const Numeric at_bottom = beam_at_boundary[layer + 1];
          const Numeric at_top    = transport * beam_at_boundary[layer];
          downward = transport * downward + dc::downward_tms_kernel(abs_mu, mu0_, thickness, at_bottom, at_top) *
                                                tms_source(transport_beam_phase, layer, p, user);
          tms[layer, p, user] = downward;
        }
      } else {
        rtepack::stokvec upward{};
        tms[expected[0] - 1, p, user] = upward;
        for (Index layer = expected[0] - 1; layer > 0; --layer) {
          const Numeric transport = ray_transport[user, layer];
          const Numeric at_top    = beam_at_boundary[layer];
          const Numeric at_bottom = transport * beam_at_boundary[layer + 1];
          upward = transport * upward + (at_top - at_bottom) * tms_source(transport_beam_phase, layer, p, user);
          tms[layer - 1, p, user] = upward;
        }
      }
    }
}

void delta_m_correction_cache::gridded_IMS(rtepack::stokvec_tensor3_view ims, const ims_convention convention) const {
  ARTS_TIME_REPORT

  const std::array<Index, 3> expected{
      static_cast<Index>(physical_tau_.size()), static_cast<Index>(phi_.size()), static_cast<Index>(user_mu_.size())};
  if (ims.shape() != expected)
    arts::user_error("gridded_IMS output has shape {:B,}, expected {:B,}", ims.shape(), expected);

  ims                      = rtepack::stokvec{};
  const Numeric beam_theta = std::acos(-mu0_);
  const Numeric sign       = convention == ims_convention::disort ? -1.0 : 1.0;
  for (Index user = 0; user < expected[2]; ++user) {
    const Numeric mu = user_mu_[user];
    if (mu > 0.0 or
        (convention == ims_convention::disort and std::abs(beam_theta - std::acos(mu)) > Constant::pi / 18.0))
      continue;
    const Numeric abs_mu = -mu;
    for (Index layer = 0; layer < expected[0]; ++layer) {
      const Index   boundary = layer + 1;
      const Numeric depth_factor =
          sign * ims_scalar_[boundary] * dc::ims_chi(physical_tau_[layer], abs_mu, ims_mu0_[boundary]);
      for (Index p = 0; p < expected[1]; ++p)
        ims[layer, p, user] = depth_factor * (ims_operator_[boundary, p, user] * beam_);
    }
  }
}

rtepack::stokvec_vector delta_m_correction_cache::evaluate(const Numeric                 tau,
                                                           const Index                   phi_index,
                                                           const beam_phase_matrix_data& transport_beam_phase,
                                                           const ims_convention          convention) const {
  rtepack::stokvec_vector tms, ims;
  TMS(tms, tau, phi_index, transport_beam_phase);
  IMS(ims, tau, phi_index, convention);
  for (Index user = 0; user < static_cast<Index>(tms.size()); ++user) tms[user] += ims[user];
  return tms;
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

namespace {
//! Throw unless a beam cosine fits the specular part: a mirror cannot reflect a beam
void check_specular_beam(const bool specular, const Numeric beam_mu) {
  if (specular and not std::isnan(beam_mu))
    arts::user_error(
        "A specular (mirror) surface part reflects the direct beam into the single mirror direction, an "
        "upward beam that VDISORT cannot represent, so a surface with one has no beam reflection (its "
        "beam cosine must be NaN, got {}); use a finite-width surface model such as Cox-Munk",
        beam_mu);
}
}  // namespace

surface_reflection::surface_reflection(
    const Index NBRDF, const Index n_out, const Index N, const Numeric beam_mu, const bool specular)
    : kernel(2, NBRDF, n_out, N, rtepack::muelmat{0.0}),
      beam_kernel(2, NBRDF, n_out, rtepack::muelmat{0.0}),
      mirror(specular ? n_out : 0, rtepack::muelmat{0.0}),
      beam_cosine(beam_mu) {
  dc::check_beam_cosine(beam_mu);
  check_specular_beam(specular, beam_mu);
}

surface_reflection::surface_reflection(rtepack::muelmat_tensor4 brdf,
                                       rtepack::muelmat_tensor3 beam_brdf,
                                       const Numeric            beam_mu,
                                       rtepack::muelmat_vector  specular)
    : kernel(std::move(brdf)), beam_kernel(std::move(beam_brdf)), mirror(std::move(specular)), beam_cosine(beam_mu) {
  if ((kernel.extent(0) != 2 or beam_kernel.shape() != std::array{Index{2}, modes(), directions()} or
       (has_specular() and static_cast<Index>(mirror.size()) != directions())))
    arts::user_error(
        "The surface reflection has the shapes {:B,}, {:B,} and {:B,}; expected [2, NBRDF, n_out, N], [2, NBRDF, n_out] "
        "and [n_out] or [0]",
        kernel.shape(),
        beam_kernel.shape(),
        mirror.shape());
  dc::check_beam_cosine(beam_mu);
  check_specular_beam(has_specular(), beam_mu);
}

surface_reflection surface_reflection::like(const Index n_out) const {
  return {modes(), n_out, streams(), beam_cosine, has_specular()};
}

void surface_reflection::set(const rtepack::muelmat_tensor4_const_view& brdf,
                             const rtepack::muelmat_tensor3_const_view& beam_brdf,
                             const Numeric                              beam_mu,
                             const rtepack::muelmat_vector_const_view&  specular) {
  // Without modes there is no kernel to set, whatever the directions
  const bool no_kernel = brdf.extent(1) == 0 and beam_brdf.extent(1) == 0 and modes() == 0;
  if ((not no_kernel and (brdf.shape() != kernel.shape() or beam_brdf.shape() != beam_kernel.shape())) or
      specular.size() != mirror.size())
    arts::user_error("The surface reflection has the shapes {:B,}, {:B,} and {:B,}, expected {:B,}, {:B,} and {:B,}",
                     brdf.shape(),
                     beam_brdf.shape(),
                     specular.shape(),
                     kernel.shape(),
                     beam_kernel.shape(),
                     mirror.shape());
  dc::check_beam_cosine(beam_mu);
  check_specular_beam(has_specular(), beam_mu);
  if (not no_kernel) {
    rtepack::muelmat_tensor4_view{kernel}      = brdf;
    rtepack::muelmat_tensor3_view{beam_kernel} = beam_brdf;
  }
  rtepack::muelmat_vector_view{mirror} = specular;
  beam_cosine                          = beam_mu;
}

void surface_reflection::set(const surface_reflection& other) {
  set(other.brdf(), other.beam_brdf(), other.beam_mu(), other.specular());
}

namespace {
//! Throw unless a surface reflects a beam of cosine mu0, or there is no beam
void check_beam_reflection(const surface_reflection& surface, const bool beam, const Numeric mu0) {
  if (not beam) return;
  if (surface.has_specular())
    arts::user_error(
        "A specular (mirror) surface reflects the direct beam into the single mirror direction, an "
        "upward beam that VDISORT cannot represent; use a finite-width surface model such as Cox-Munk");
  if (surface.modes() > 0 and not(surface.beam_mu() == mu0))
    arts::user_error(
        "The surface's reflection of the direct beam is for the beam cosine {} (NaN: computed without "
        "a beam), but the beam has the cosine {}.  Compute the surface reflection with the beam: in "
        "the settings, set the surface after the sun.",
        surface.beam_mu(),
        mu0);
}
}  // namespace

phase_matrix_data combine_phase_matrices(const rtepack::muelmat_tensor4& cosine, const rtepack::muelmat_tensor4& sine) {
  if (cosine.shape() != sine.shape())
    arts::user_error(
        "Cosine and sine phase matrices have different shapes: {:B,} and {:B,}", cosine.shape(), sine.shape());
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

beam_phase_matrix_data combine_beam_phase_matrices(const rtepack::muelmat_tensor3& cosine,
                                                   const rtepack::muelmat_tensor3& sine) {
  if (cosine.shape() != sine.shape())
    arts::user_error(
        "Cosine and sine beam phase matrices have different shapes: {:B,} and {:B,}", cosine.shape(), sine.shape());
  const auto [nfourier, nlayers, nout] = cosine.shape();
  beam_phase_matrix_data out(2, nfourier, nlayers, nout, rtepack::muelmat{0.0});

  // Unlike the diffuse operator of Eq. (81), the beam is a delta in azimuth
  // at phi0 with only cosine terms, for every Stokes component of the beam.
  // Its scattered source sum_m [C^m cos m(phi0 - phi) + S^m sin m(phi0 - phi)] b
  // feeds the cosine system (I^c, Q^c, U^s, V^s) with rows I, Q of C^m and
  // rows U, V of S^m, and the sine system (I^s, Q^s, U^c, V^c) with rows I, Q
  // of S^m and rows U, V of C^m, in all four columns.
  for (Index m = 0; m < nfourier; ++m) {
    for (Index l = 0; l < nlayers; ++l) {
      for (Index i = 0; i < nout; ++i) {
        for (Index r = 0; r < stokes_dimension; ++r) {
          for (Index c = 0; c < stokes_dimension; ++c) {
            out[cosine_mode, m, l, i][r, c] = r < 2 ? cosine[m, l, i][r, c] : sine[m, l, i][r, c];
            out[sine_mode, m, l, i][r, c]   = r < 2 ? sine[m, l, i][r, c] : cosine[m, l, i][r, c];
          }
        }
      }
    }
  }
  return out;
}

delta_m_transport_data delta_m_preprocess(const AscendingGrid&          physical_tau,
                                          const Vector&                 physical_omega,
                                          const Vector&                 fraction,
                                          const phase_matrix_data&      original_phase_matrix,
                                          const phase_matrix_data&      removed_phase_matrix,
                                          const beam_phase_matrix_data& original_beam_phase_matrix,
                                          const beam_phase_matrix_data& removed_beam_phase_matrix) {
  const Index nlayers = static_cast<Index>(physical_tau.size());
  if (nlayers == 0) arts::user_error("Polarized delta-M preprocessing requires at least one layer");
  if (static_cast<Index>(physical_omega.size()) != nlayers or static_cast<Index>(fraction.size()) != nlayers)
    arts::user_error("Polarized delta-M arrays must contain {} layers, got omega={} and fraction={}",
                     nlayers,
                     physical_omega.size(),
                     fraction.size());

  const auto phase_shape = original_phase_matrix.shape();
  if (phase_shape[0] != 2 or phase_shape[2] != nlayers or phase_shape[3] != phase_shape[4])
    arts::user_error(
        "Original polarized phase matrix has incompatible shape {:B,}; expected "
        "[2, NFourier, {}, NQuad, NQuad] Mueller blocks",
        phase_shape,
        nlayers);
  if (removed_phase_matrix.shape() != phase_shape)
    arts::user_error("Original and removed polarized phase matrices have different shapes: {:B,} and {:B,}",
                     phase_shape,
                     removed_phase_matrix.shape());

  const bool has_beam_phase = original_beam_phase_matrix.size() != 0 or removed_beam_phase_matrix.size() != 0;
  if (has_beam_phase) {
    const std::array expected_beam_shape{phase_shape[0], phase_shape[1], phase_shape[2], phase_shape[3]};
    if (original_beam_phase_matrix.shape() != expected_beam_shape or
        removed_beam_phase_matrix.shape() != expected_beam_shape)
      arts::user_error("Original and removed beam phase matrices must both have shape {:B,}; got {:B,} and {:B,}",
                       expected_beam_shape,
                       original_beam_phase_matrix.shape(),
                       removed_beam_phase_matrix.shape());
  }

  Vector  scaled_tau(nlayers), transport_omega(nlayers), coordinate_scale(nlayers), coordinate_offset(nlayers);
  Numeric physical_top = 0.0;
  Numeric scaled_top   = 0.0;
  for (Index layer = 0; layer < nlayers; ++layer) {
    const Numeric omega = physical_omega[layer];
    const Numeric f     = fraction[layer];
    if (not std::isfinite(omega) or omega < 0.0 or omega > 1.0 or not std::isfinite(f) or f < 0.0 or f >= 1.0)
      arts::user_error("Invalid polarized delta-M omega/fraction in layer {}: {}, {}", layer, omega, f);
    const Numeric scale        = 1.0 - omega * f;
    const Numeric scaled_base  = scaled_top;
    scaled_top                += scale * (physical_tau[layer] - physical_top);
    scaled_tau[layer]          = scaled_top;
    transport_omega[layer]     = omega * (1.0 - f) / scale;
    coordinate_scale[layer]    = 1.0 / scale;
    coordinate_offset[layer]   = physical_top - scaled_base / scale;
    physical_top               = physical_tau[layer];
  }

  phase_matrix_data transport_phase(
      phase_shape[0], phase_shape[1], phase_shape[2], phase_shape[3], phase_shape[4], rtepack::muelmat{0.0});
  for (Index alpha = 0; alpha < phase_shape[0]; ++alpha)
    for (Index mode = 0; mode < phase_shape[1]; ++mode)
      for (Index layer = 0; layer < nlayers; ++layer)
        for (Index out = 0; out < phase_shape[3]; ++out)
          for (Index in = 0; in < phase_shape[4]; ++in)
            for (Index row = 0; row < stokes_dimension; ++row)
              for (Index column = 0; column < stokes_dimension; ++column) {
                const Numeric original = original_phase_matrix[alpha, mode, layer, out, in][row, column];
                const Numeric removed  = removed_phase_matrix[alpha, mode, layer, out, in][row, column];
                if (not std::isfinite(original) or not std::isfinite(removed))
                  arts::user_error("Non-finite polarized delta-M phase value at [{}, {}, {}, {}, {}][{}, {}]",
                                   alpha,
                                   mode,
                                   layer,
                                   out,
                                   in,
                                   row,
                                   column);
                transport_phase[alpha, mode, layer, out, in][row, column] =
                    (original - fraction[layer] * removed) / (1.0 - fraction[layer]);
              }

  beam_phase_matrix_data transport_beam_phase;
  if (has_beam_phase) {
    transport_beam_phase.resize(phase_shape[0], phase_shape[1], phase_shape[2], phase_shape[3]);
    for (Index alpha = 0; alpha < phase_shape[0]; ++alpha)
      for (Index mode = 0; mode < phase_shape[1]; ++mode)
        for (Index layer = 0; layer < nlayers; ++layer)
          for (Index out = 0; out < phase_shape[3]; ++out)
            for (Index row = 0; row < stokes_dimension; ++row)
              for (Index column = 0; column < stokes_dimension; ++column) {
                const Numeric original = original_beam_phase_matrix[alpha, mode, layer, out][row, column];
                const Numeric removed  = removed_beam_phase_matrix[alpha, mode, layer, out][row, column];
                if (not std::isfinite(original) or not std::isfinite(removed))
                  arts::user_error("Non-finite polarized delta-M beam phase value at [{}, {}, {}, {}][{}, {}]",
                                   alpha,
                                   mode,
                                   layer,
                                   out,
                                   row,
                                   column);
                transport_beam_phase[alpha, mode, layer, out][row, column] =
                    (original - fraction[layer] * removed) / (1.0 - fraction[layer]);
              }
  }

  return {.tau                      = AscendingGrid{std::move(scaled_tau)},
          .omega                    = std::move(transport_omega),
          .phase_matrix             = std::move(transport_phase),
          .beam_phase_matrix        = std::move(transport_beam_phase),
          .source_coordinate_scale  = std::move(coordinate_scale),
          .source_coordinate_offset = std::move(coordinate_offset)};
}

main_data::main_data(const Index NLayers_,
                     const Index NQuad_,
                     const Index NFourier_,
                     const Index Nscoeffs_,
                     const Index NBDRF_,
                     const bool  specular)
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
      surface_refl(NBDRF, N, N, std::numeric_limits<Numeric>::quiet_NaN(), specular),
      beam_stokes{},
      beam_phase_matrix(2, NFourier, NLayers, NQuad, rtepack::muelmat{0.0}),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      half_range_barycentric_weights(N),
      scaled_source_poly_coeffs(NLayers, Nscoeffs),
      transport(NState, NState),
      transport_work(NState, NState),
      real_eigenvectors(NState, NState),
      eigenvalues_real(NState),
      eigenvalues_imag(NState),
      eigenvectors(NState, NState),
      eigenvalues(NState),
      eigen_order(static_cast<std::size_t>(NState)),
      eigen_partner(static_cast<std::size_t>(NState)),
      eigen_position(static_cast<std::size_t>(NState)),
      half_sum(NHalfState, NHalfState),
      half_difference(NHalfState, NHalfState),
      half_product(NHalfState, NHalfState),
      half_eigenvectors(NHalfState, NHalfState),
      half_images(NHalfState, NHalfState),
      half_eigenvalues_real(NHalfState),
      half_eigenvalues_imag(NHalfState),
      half_eigen_work(static_cast<std::size_t>(NHalfState)),
      eigen_work(static_cast<std::size_t>(NState)),
      solve_work(static_cast<std::size_t>(NState)),
      state_rhs(NState),
      state_next(NState),
      boundary_lhs(3 * NHalfState - 1, 3 * NHalfState - 1, NLayers * NState, NLayers * NState),
      boundary_rhs(NLayers * NState),
      mode_exponential(NLayers, 2, NState),
      surface_operator(NHalfState, NHalfState),
      direct_reflection(NHalfState),
      layer_field(2 * NFourier, NState),
      G_collect(2, NFourier, NLayers, NState, NState),
      K_collect(2, NFourier, NLayers, NState),
      GC_collect(2, NFourier, NLayers, NState),
      B_collect(2, NFourier, NLayers, NQuad),
      source_collect(2, NFourier, NLayers, NQuad, Nscoeffs),
      um(NLayers, 2, NFourier, NQuad),
      top_anchored(static_cast<std::size_t>(2 * NFourier * NLayers * NState), 0),
      conjugate_partner(static_cast<std::size_t>(2 * NFourier * NLayers * NState), -1),
      conservative_pair_index(static_cast<std::size_t>(NLayers), std::array<Index, 2>{-1, -1}),
      conservative_pair_kappa(NLayers) {
  if (NQuad <= 0 or NQuad % 2 != 0) arts::user_error("NQuad must be a positive even number, got {}", NQuad);
  dc::initialize_streams(mu_arr, inv_mu_arr, W);
  dc::barycentric_weights(half_range_barycentric_weights, mu_arr[Range{0, N}]);
}

main_data::main_data(Index                     NQuad_,
                     Index                     NFourier_,
                     AscendingGrid             tau_arr_,
                     Vector                    omega_arr_,
                     phase_matrix_data         phase_matrix_,
                     rtepack::stokvec_tensor3  boundary_up_,
                     rtepack::stokvec_tensor3  boundary_down_,
                     rtepack::stokvec_matrix   source_poly_coeffs_,
                     const surface_reflection& surface_,
                     Numeric                   mu0_,
                     rtepack::stokvec          beam_stokes_,
                     Numeric                   phi0_,
                     beam_phase_matrix_data    beam_phase_matrix_,
                     Vector                    source_coordinate_scale_,
                     Vector                    source_coordinate_offset_)
    : main_data(static_cast<Index>(tau_arr_.size()),
                NQuad_,
                NFourier_,
                source_poly_coeffs_.size() == 0 ? 0 : source_poly_coeffs_.shape()[1],
                surface_.modes(),
                surface_.has_specular()) {
  tau_arr            = std::move(tau_arr_);
  omega_arr          = std::move(omega_arr_);
  phase_matrix       = std::move(phase_matrix_);
  boundary_up        = std::move(boundary_up_);
  boundary_down      = std::move(boundary_down_);
  source_poly_coeffs = std::move(source_poly_coeffs_);
  if (source_coordinate_scale_.size() != 0) source_coordinate_scale = std::move(source_coordinate_scale_);
  if (source_coordinate_offset_.size() != 0) source_coordinate_offset = std::move(source_coordinate_offset_);
  surface_refl.set(surface_);
  mu0             = mu0_;
  beam_stokes     = std::move(beam_stokes_);
  phi0            = phi0_;
  has_beam_source = polarized_source(beam_stokes);
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
  const auto [negative, positive] = conservative_pair_index[static_cast<std::size_t>(layer)];
  if (alpha == cosine_mode and m == 0 and negative >= 0 and (eigen == negative or eigen == positive)) {
    const Numeric center  = 0.5 * (layer_top(layer) + tau_arr[layer]);
    const Numeric s       = tau - center;
    const Numeric kappa   = conservative_pair_kappa[layer];
    const auto    basis   = dc::centered_pair(kappa, s);
    const auto    columns = dc::centered_pair_columns(
        G_collect[alpha, m, layer, state, negative], G_collect[alpha, m, layer, state, positive], basis);
    return columns[eigen == negative ? 0 : 1];
  }

  const Numeric anchor = top_anchored[anchor_index(alpha, m, layer, eigen)] ? layer_top(layer) : tau_arr[layer];
  return G_collect[alpha, m, layer, state, eigen] * std::exp(K_collect[alpha, m, layer, eigen] * (tau - anchor));
}

Numeric main_data::particular(
    const Index alpha, const Index m, const Index layer, const Index state, const Numeric tau) const {
  Numeric value = 0.0;
  // Each layer source polynomial uses x = offset + scale * tau.
  const Numeric source_tau = std::fma(source_coordinate_scale[layer], tau, source_coordinate_offset[layer]);
  const Index   stream     = state / stokes_dimension;
  const Index   stokes     = state % stokes_dimension;
  value                    = dc::horner_polynomial(
      Nscoeffs, source_tau, [&](const Index p) { return source_collect[alpha, m, layer, stream, p][stokes]; });
  if (has_beam_source) value += B_collect[alpha, m, layer, stream][stokes] * std::exp(-tau / mu0);
  return value;
}

void main_data::check_input_size() const {
  dc::check_layer_input_sizes(NLayers, tau_arr, omega_arr);
  if ((source_poly_coeffs.shape() != std::array<Index, 2>{NLayers, Nscoeffs}))
    arts::user_error("source_poly_coeffs has shape {:B,}, expected [{}, {}] Stokes blocks",
                     source_poly_coeffs.shape(),
                     NLayers,
                     Nscoeffs);
  if (static_cast<Index>(source_coordinate_scale.size()) != NLayers or
      static_cast<Index>(source_coordinate_offset.size()) != NLayers)
    arts::user_error("Source-coordinate arrays have sizes {} and {}, expected {}",
                     source_coordinate_scale.size(),
                     source_coordinate_offset.size(),
                     NLayers);
  if ((phase_matrix.shape() != std::array<Index, 5>{2, NFourier, NLayers, NQuad, NQuad}))
    arts::user_error("phase_matrix has shape {:B,}, expected [2, {}, {}, {}, {}] Mueller blocks",
                     phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     NQuad,
                     NQuad);
  if ((boundary_up.shape() != std::array<Index, 3>{2, NFourier, N}))
    arts::user_error(
        "boundary_up has shape {:B,}, expected [2, {}, {}] Stokes blocks", boundary_up.shape(), NFourier, N);
  if ((boundary_down.shape() != std::array<Index, 3>{2, NFourier, N}))
    arts::user_error(
        "boundary_down has shape {:B,}, expected [2, {}, {}] Stokes blocks", boundary_down.shape(), NFourier, N);
  if ((beam_phase_matrix.shape() != std::array<Index, 4>{2, NFourier, NLayers, NQuad}))
    arts::user_error("beam_phase_matrix has shape {:B,}, expected [2, {}, {}, {}] Mueller blocks",
                     beam_phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     NQuad);
}

void main_data::check_input_value() const {
  dc::check_layer_input_values(tau_arr, omega_arr, mu_arr, mu0, phi0, has_beam_source);
  if (NFourier <= 0) arts::user_error("VDISORT requires at least one Fourier mode");
  if (beam_stokes[0] < 0.0) arts::user_error("The beam I component must be non-negative, got {}", beam_stokes[0]);
  const Numeric polarized_norm = std::hypot(beam_stokes[1], beam_stokes[2], beam_stokes[3]);
  if (polarized_norm > beam_stokes[0] * (1.0 + 1e-12))
    arts::user_error(
        "The beam Stokes vector is non-physical: sqrt(Q^2+U^2+V^2)={} > I={}", polarized_norm, beam_stokes.I());
  if (has_beam_source and mu0 == 0.0) arts::user_error("A direct beam requires mu0 > 0");
  if (NBDRF > NFourier) arts::user_error("There are {} BRDF modes but only {} Fourier modes", NBDRF, NFourier);
}

Index main_data::tau_index(const Numeric tau) const { return dc::layer_index(tau_arr, tau); }

bool main_data::mirror_symmetric(const Index alpha, const Index m, const Index layer) const {
  Numeric scale = 0.0, deviation = 0.0;
  for (Index i = 0; i < N; ++i) {
    for (Index j = 0; j < N; ++j) {
      const rtepack::muelmat same     = rtepack::mirror(phase_matrix[alpha, m, layer, i, j]);
      const rtepack::muelmat opposite = rtepack::mirror(phase_matrix[alpha, m, layer, i, N + j]);
      for (Size n = 0; n < 16; ++n) {
        deviation = std::max({deviation,
                              std::abs(same.data[n] - phase_matrix[alpha, m, layer, N + i, N + j].data[n]),
                              std::abs(opposite.data[n] - phase_matrix[alpha, m, layer, N + i, j].data[n])});
        scale     = std::max({scale, std::abs(same.data[n]), std::abs(opposite.data[n])});
      }
    }
  }
  return deviation <= 1e-12 * scale;
}

/* With S = diag(1, Delta), Delta flipping U and V of each downward stream,
   the mirror symmetry makes S A S = [alpha, beta~; -beta~, -alpha], where
   alpha and beta are the upward rows of A and beta~ = beta Delta.  Its
   eigenvalues come in pairs +-k, and with u = x + y and w = x - y for an
   eigenvector [x; y] of k,
     (alpha - beta~)(alpha + beta~) u = k^2 u,   w = (alpha + beta~) u / k.
   The eigenvectors of A are [x; Delta y] for k and [y; Delta x] for -k. */
bool main_data::half_size_eigensolution() {
  for (Index r = 0; r < NHalfState; ++r) {
    for (Index c = 0; c < NHalfState; ++c) {
      const Numeric a       = transport[r, c];
      const Numeric b       = c % stokes_dimension < 2 ? transport[r, NHalfState + c] : -transport[r, NHalfState + c];
      half_sum[r, c]        = a + b;
      half_difference[r, c] = a - b;
    }
  }
  mult(half_product, half_difference, half_sum);
  diagonalize_inplace(half_eigenvectors, half_eigenvalues_real, half_eigenvalues_imag, half_product, half_eigen_work);
  for (Index j = 0; j < NHalfState; ++j)
    if (half_eigenvalues_real[j] == 0.0 and half_eigenvalues_imag[j] == 0.0) return false;
  mult(half_images, half_sum, half_eigenvectors);

  const auto set_pair = [&](const Index e, const Complex k, const auto& u, const auto& image) {
    // The eigenpair of k at e and that of -k at e + 1
    eigenvalues[e]     = k;
    eigenvalues[e + 1] = -k;
    for (Index r = 0; r < NHalfState; ++r) {
      const Complex w                     = image(r) / k;
      const Complex x                     = 0.5 * (u(r) + w);
      const Complex y                     = 0.5 * (u(r) - w);
      const Numeric delta                 = r % stokes_dimension < 2 ? 1.0 : -1.0;
      eigenvectors[r, e]                  = x;
      eigenvectors[NHalfState + r, e]     = delta * y;
      eigenvectors[r, e + 1]              = y;
      eigenvectors[NHalfState + r, e + 1] = delta * x;
    }
  };
  const auto set_conjugate = [&](const Index e, const Index from) {
    eigenvalues[e] = std::conj(eigenvalues[from]);
    for (Index r = 0; r < NState; ++r) eigenvectors[r, e] = std::conj(eigenvectors[r, from]);
    eigen_partner[static_cast<std::size_t>(e)]    = from;
    eigen_partner[static_cast<std::size_t>(from)] = e;
  };

  stdr::fill(eigen_partner, Index{-1});
  for (Index j = 0; j < NHalfState; ++j) {
    const Index e = 2 * j;
    if (half_eigenvalues_imag[j] == 0.0) {
      const Numeric lambda = half_eigenvalues_real[j];
      const auto    u      = [&](const Index r) { return Complex{half_eigenvectors[r, j]}; };
      const auto    image  = [&](const Index r) { return Complex{half_images[r, j]}; };
      if (lambda > 0.0) {
        set_pair(e, Complex{std::sqrt(lambda)}, u, image);
      } else {
        // k = i sqrt(-lambda), and -k is its conjugate
        set_pair(e, Complex{0.0, std::sqrt(-lambda)}, u, image);
        set_conjugate(e + 1, e);
      }
    } else {
      // A conjugate pair of squared eigenvalues, the first with the positive imaginary part
      const auto u     = [&](const Index r) { return Complex{half_eigenvectors[r, j], half_eigenvectors[r, j + 1]}; };
      const auto image = [&](const Index r) { return Complex{half_images[r, j], half_images[r, j + 1]}; };
      set_pair(e, std::sqrt(Complex{half_eigenvalues_real[j], half_eigenvalues_imag[j]}), u, image);
      set_conjugate(e + 2, e);
      set_conjugate(e + 3, e + 1);
      ++j;
    }
  }
  return true;
}

void main_data::diagonalize() {
  ARTS_TIME_REPORT

  // The real, generally nonsymmetric 8N x 8N operator may have complex eigenpairs.  The real
  // eigensolver returns a conjugate pair as the real and imaginary parts of its first vector.
  stdr::fill(conservative_pair_index, std::array<Index, 2>{-1, -1});
  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      for (Index l = 0; l < NLayers; ++l) {
        transport     = 0.0;
        bool scatters = false;
        for (Index i = 0; i < NQuad; ++i) {
          for (Index so = 0; so < stokes_dimension; ++so) {
            const Index row = state_index(i, so);
            for (Index j = 0; j < NQuad; ++j) {
              const Numeric factor = -0.5 * omega_arr[l] * W[j % N] * inv_mu_arr[i];
              for (Index si = 0; si < stokes_dimension; ++si) {
                const Numeric phase                 = phase_matrix[alpha, m, l, i, j][so, si];
                transport[row, state_index(j, si)] += factor * phase;
                scatters                           |= phase != 0.0;
              }
            }
            transport[row, row] += inv_mu_arr[i];
          }
        }
        scatters &= omega_arr[l] != 0.0;

        if (not scatters) {
          // Without scattering the operator is diag(1 / mu): its eigenvectors are the states themselves
          stdr::fill(eigen_partner, Index{-1});
          eigenvectors = 0.0;
          for (Index e = 0; e < NState; ++e) {
            eigenvalues[e]     = transport[e, e];
            eigenvectors[e, e] = 1.0;
          }
        } else if (not(alpha == cosine_mode and m == 0) and mirror_symmetric(alpha, m, l) and
                   half_size_eigensolution()) {
          // The half-size problem gave the eigenpairs (the conservation pair, in the cosine system of order zero,
          // keeps the full problem, as k^2 loses the relative accuracy of a small k)
        } else {
          stdr::fill(eigen_partner, Index{-1});
          transport_work = transport;
          diagonalize_inplace(real_eigenvectors, eigenvalues_real, eigenvalues_imag, transport_work, eigen_work);
          for (Index e = 0; e < NState; ++e) {
            if (eigenvalues_imag[e] == 0.0) {
              eigenvalues[e] = eigenvalues_real[e];
              for (Index s = 0; s < NState; ++s) eigenvectors[s, e] = real_eigenvectors[s, e];
            } else {
              eigenvalues[e]     = Complex{eigenvalues_real[e], eigenvalues_imag[e]};
              eigenvalues[e + 1] = Complex{eigenvalues_real[e + 1], eigenvalues_imag[e + 1]};
              for (Index s = 0; s < NState; ++s) {
                eigenvectors[s, e]     = Complex{real_eigenvectors[s, e], real_eigenvectors[s, e + 1]};
                eigenvectors[s, e + 1] = Complex{real_eigenvectors[s, e], -real_eigenvectors[s, e + 1]};
              }
              eigen_partner[static_cast<std::size_t>(e)]     = e + 1;
              eigen_partner[static_cast<std::size_t>(e + 1)] = e;
              ++e;
            }
          }
        }

        std::iota(eigen_order.begin(), eigen_order.end(), Index{0});
        stdr::sort(eigen_order, [&](const Index a, const Index b) {
          if (eigenvalues[a].real() != eigenvalues[b].real()) return eigenvalues[a].real() < eigenvalues[b].real();
          return eigenvalues[a].imag() < eigenvalues[b].imag();
        });
        check_eigenvalue_direction_split(eigenvalues, eigen_order, NHalfState, alpha, m, l);

        for (Index e = 0; e < NState; ++e)
          eigen_position[static_cast<std::size_t>(eigen_order[static_cast<std::size_t>(e)])] = e;
        for (Index e = 0; e < NState; ++e) {
          const Index old                            = eigen_order[static_cast<std::size_t>(e)];
          const Index partner                        = eigen_partner[static_cast<std::size_t>(old)];
          K_collect[alpha, m, l, e]                  = eigenvalues[old];
          top_anchored[anchor_index(alpha, m, l, e)] = static_cast<unsigned char>(e < NHalfState);
          conjugate_partner[static_cast<std::size_t>(anchor_index(alpha, m, l, e))] =
              partner < 0 ? -1 : eigen_position[static_cast<std::size_t>(partner)];
          for (Index s = 0; s < NState; ++s) G_collect[alpha, m, l, s, e] = eigenvectors[s, old];
        }

        if (alpha == cosine_mode and m == 0) {
          Index negative = 0;
          for (Index e = 1; e < NHalfState; ++e)
            if (std::abs(K_collect[alpha, m, l, e]) < std::abs(K_collect[alpha, m, l, negative])) negative = e;
          Index positive = NHalfState;
          for (Index e = NHalfState + 1; e < NState; ++e)
            if (std::abs(K_collect[alpha, m, l, e]) < std::abs(K_collect[alpha, m, l, positive])) positive = e;

          const Complex lambda_negative = K_collect[alpha, m, l, negative];
          const Complex lambda_positive = K_collect[alpha, m, l, positive];
          const Numeric measured_kappa  = std::sqrt(std::abs(-lambda_negative * lambda_positive));

          // A conservative phase operator has both isotropic right invariance
          // and integrated-Stokes-I left invariance.
          Numeric conservation_residual = 0.0;
          Numeric conservation_scale    = 1.0;
          for (Index i = 0; i < NQuad; ++i) {
            for (Index so = 0; so < stokes_dimension; ++so) {
              Numeric scattered = 0.0;
              Numeric magnitude = so == 0 ? 1.0 : 0.0;
              for (Index j = 0; j < NQuad; ++j) {
                const Numeric term  = 0.5 * W[j % N] * phase_matrix[alpha, m, l, i, j][so, 0];
                scattered          += term;
                magnitude          += std::abs(term);
              }
              conservation_residual =
                  std::max(conservation_residual, std::abs(inv_mu_arr[i] * ((so == 0 ? 1.0 : 0.0) - scattered)));
              conservation_scale = std::max(conservation_scale, std::abs(inv_mu_arr[i]) * magnitude);
            }
          }
          Numeric energy_residual = 0.0;
          Numeric energy_scale    = 1.0;
          for (Index j = 0; j < NQuad; ++j) {
            for (Index si = 0; si < stokes_dimension; ++si) {
              Numeric integrated = 0.0;
              Numeric magnitude  = si == 0 ? 1.0 : 0.0;
              for (Index i = 0; i < NQuad; ++i) {
                const Numeric term  = 0.5 * W[i % N] * phase_matrix[alpha, m, l, i, j][0, si];
                integrated         += term;
                magnitude          += std::abs(term);
              }
              energy_residual = std::max(energy_residual, std::abs((si == 0 ? 1.0 : 0.0) - integrated));
              energy_scale    = std::max(energy_scale, magnitude);
            }
          }
          const bool normalized_conservative =
              conservation_residual <= 1.0e-10 * conservation_scale and energy_residual <= 1.0e-10 * energy_scale;
          const bool requested_stable_pair = dc::use_centered_pair(omega_arr[l], measured_kappa);
          if (normalized_conservative and requested_stable_pair) {
            Numeric kappa = omega_arr[l] == 1.0 ? 0.0 : measured_kappa;

            const auto normalization = [&](const Index eigen) {
              Complex value = 0.0;
              for (Index stream = 0; stream < NQuad; ++stream)
                value += 0.5 * W[stream % N] * G_collect[alpha, m, l, state_index(stream, 0), eigen];
              return value;
            };

            if (omega_arr[l] != 1.0) {
              // Near conservation, eigensolution roundoff can perturb two
              // real roots into a slightly complex conjugate pair.
              bool use_real_pair = measured_kappa <= std::numeric_limits<Numeric>::epsilon() * conservation_scale or
                                   std::max(std::abs(lambda_negative.imag()), std::abs(lambda_positive.imag())) >
                                       1.0e-8 * measured_kappa;
              if (not use_real_pair) {
                const Complex negative_norm = normalization(negative);
                const Complex positive_norm = normalization(positive);
                if (std::abs(negative_norm) <= std::numeric_limits<Numeric>::epsilon() or
                    std::abs(positive_norm) <= std::numeric_limits<Numeric>::epsilon())
                  arts::user_error("VDISORT could not normalize the conservative eigenpair in layer {}: {} and {}",
                                   l,
                                   negative_norm,
                                   positive_norm);
                Numeric pair_imag  = 0.0;
                Numeric pair_scale = 1.0;
                for (Index state = 0; state < NState; ++state) {
                  const Complex v_negative                = G_collect[alpha, m, l, state, negative] / negative_norm;
                  const Complex v_positive                = G_collect[alpha, m, l, state, positive] / positive_norm;
                  G_collect[alpha, m, l, state, negative] = 0.5 * (v_negative + v_positive);
                  G_collect[alpha, m, l, state, positive] = (v_positive - v_negative) / (2.0 * kappa);
                  pair_imag  = std::max({pair_imag,
                                         std::abs(G_collect[alpha, m, l, state, negative].imag()),
                                         std::abs(G_collect[alpha, m, l, state, positive].imag())});
                  pair_scale = std::max({pair_scale,
                                         std::abs(G_collect[alpha, m, l, state, negative].real()),
                                         std::abs(G_collect[alpha, m, l, state, positive].real())});
                }
                use_real_pair = pair_imag > 1.0e-8 * pair_scale;
              }

              if (use_real_pair) {
                Index real_negative = -1;
                Index real_positive = -1;
                for (Index eigen = 0; eigen < NState; ++eigen) {
                  if (eigenvalues_imag[eigen] != 0.0) continue;
                  if (eigenvalues_real[eigen] < 0.0 and
                      (real_negative < 0 or
                       std::abs(eigenvalues_real[eigen]) < std::abs(eigenvalues_real[real_negative])))
                    real_negative = eigen;
                  if (eigenvalues_real[eigen] > 0.0 and
                      (real_positive < 0 or
                       std::abs(eigenvalues_real[eigen]) < std::abs(eigenvalues_real[real_positive])))
                    real_positive = eigen;
                }
                if (real_negative < 0 or real_positive < 0)
                  arts::user_error("VDISORT could not resolve the real conservative eigenpair in layer {}", l);

                kappa = std::sqrt(std::abs(-eigenvalues_real[real_negative] * eigenvalues_real[real_positive]));
                const auto real_normalization = [&](const Index eigen) {
                  Numeric value = 0.0;
                  for (Index stream = 0; stream < NQuad; ++stream)
                    value += 0.5 * W[stream % N] * real_eigenvectors[state_index(stream, 0), eigen];
                  return value;
                };
                const Numeric negative_norm = real_normalization(real_negative);
                const Numeric positive_norm = real_normalization(real_positive);
                if (std::abs(negative_norm) <= std::numeric_limits<Numeric>::epsilon() or
                    std::abs(positive_norm) <= std::numeric_limits<Numeric>::epsilon())
                  arts::user_error("VDISORT could not normalize the real conservative eigenpair in layer {}: {} and {}",
                                   l,
                                   negative_norm,
                                   positive_norm);
                for (Index state = 0; state < NState; ++state) {
                  const Numeric v_negative                = real_eigenvectors[state, real_negative] / negative_norm;
                  const Numeric v_positive                = real_eigenvectors[state, real_positive] / positive_norm;
                  G_collect[alpha, m, l, state, negative] = 0.5 * (v_negative + v_positive);
                  G_collect[alpha, m, l, state, positive] = (v_positive - v_negative) / (2.0 * kappa);
                }
              }

              Numeric pair_residual = 0.0;
              Numeric pair_scale    = 1.0;
              for (Index row = 0; row < NState; ++row) {
                Complex ax        = 0.0;
                Complex ar        = 0.0;
                Numeric row_scale = 0.0;
                for (Index column = 0; column < NState; ++column) {
                  const Complex term_x  = transport[row, column] * G_collect[alpha, m, l, column, negative];
                  const Complex term_r  = transport[row, column] * G_collect[alpha, m, l, column, positive];
                  ax                   += term_x;
                  ar                   += term_r;
                  row_scale            += std::abs(term_x) + std::abs(term_r);
                }
                const Complex expected_ax = kappa * kappa * G_collect[alpha, m, l, row, positive];
                const Complex expected_ar = G_collect[alpha, m, l, row, negative];
                pair_residual = std::max({pair_residual, std::abs(ax - expected_ax), std::abs(ar - expected_ar)});
                pair_scale    = std::max(pair_scale, row_scale + std::abs(expected_ax) + std::abs(expected_ar));
              }
              if (pair_residual > 1.0e-8 * pair_scale)
                arts::user_error("VDISORT conservative eigenpair in layer {} has relative residual {}",
                                 l,
                                 pair_residual / pair_scale);
            } else {
              // Exact conservation has a size-two Jordan block with A X = 0,
              // A R = X, and the gauge <R_I> = 0.
              transport_work = transport;
              state_rhs      = 0.0;
              for (Index stream = 0; stream < NQuad; ++stream) state_rhs[state_index(stream, 0)] = 1.0;

              Index constraint_stream = 0;
              for (Index stream = 1; stream < NQuad; ++stream)
                if (std::abs(W[stream % N] * mu_arr[stream]) >
                    std::abs(W[constraint_stream % N] * mu_arr[constraint_stream]))
                  constraint_stream = stream;
              const Index constraint_row = state_index(constraint_stream, 0);
              for (Index column = 0; column < NState; ++column) transport_work[constraint_row, column] = 0.0;
              for (Index stream = 0; stream < NQuad; ++stream)
                transport_work[constraint_row, state_index(stream, 0)] = 0.5 * W[stream % N];
              state_rhs[constraint_row] = 0.0;
              solve_inplace(state_rhs, transport_work, solve_work);

              Numeric residual = 0.0;
              Numeric scale    = 1.0;
              for (Index row = 0; row < NState; ++row) {
                Numeric value     = -(row % stokes_dimension == 0 ? 1.0 : 0.0);
                Numeric row_scale = std::abs(value);
                for (Index column = 0; column < NState; ++column) {
                  value     += transport[row, column] * state_rhs[column];
                  row_scale += std::abs(transport[row, column] * state_rhs[column]);
                }
                residual = std::max(residual, std::abs(value));
                scale    = std::max(scale, row_scale);
              }
              if (residual > 1.0e-8 * scale)
                arts::user_error(
                    "VDISORT could not construct the exact conservative mode in layer {}: "
                    "relative residual {}",
                    l,
                    residual / scale);

              for (Index state = 0; state < NState; ++state) {
                G_collect[alpha, m, l, state, negative] = state % stokes_dimension == 0 ? 1.0 : 0.0;
                G_collect[alpha, m, l, state, positive] = state_rhs[state];
              }
            }
            // The pair is real now; a complex one can only have been its own conjugate
            auto& negative_partner = conjugate_partner[static_cast<std::size_t>(anchor_index(alpha, m, l, negative))];
            auto& positive_partner = conjugate_partner[static_cast<std::size_t>(anchor_index(alpha, m, l, positive))];
            if ((negative_partner >= 0 and negative_partner != positive) or
                (positive_partner >= 0 and positive_partner != negative))
              arts::user_error("VDISORT conservative eigenpair in layer {} splits a conjugate pair", l);
            negative_partner = -1;
            positive_partner = -1;

            conservative_pair_index[static_cast<std::size_t>(l)] = {negative, positive};
            conservative_pair_kappa[l]                           = kappa;
            K_collect[alpha, m, l, negative]                     = -kappa;
            K_collect[alpha, m, l, positive]                     = kappa;
          }
        }

        if (has_beam_source) {
          const rtepack::stokvec& beam                      = beam_stokes;
          bool                    has_scattered_beam_source = false;
          const Numeric           epsilon                   = m == 0 ? 1.0 : 2.0;
          for (Index i = 0; i < NQuad; ++i) {
            const rtepack::stokvec scattered = beam_phase_matrix[alpha, m, l, i] * beam;
            for (Index so = 0; so < stokes_dimension; ++so) {
              state_rhs[state_index(i, so)] =
                  inv_mu_arr[i] * epsilon * omega_arr[l] * scattered[so] / (4.0 * Constant::pi);
              has_scattered_beam_source |= state_rhs[state_index(i, so)] != 0.0;
            }
          }
          if (not has_scattered_beam_source) {
            B_collect[alpha, m, l] = rtepack::stokvec{};
            continue;
          }
          diagonal(transport) += 1.0 / mu0;
          solve_inplace(state_rhs, transport, solve_work);
          for (Index i = 0; i < NQuad; ++i)
            for (Index s = 0; s < stokes_dimension; ++s) B_collect[alpha, m, l, i][s] = state_rhs[state_index(i, s)];
        }
      }
    }
  }
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
  if (not has_source_poly) return;

  // Transfer-equation emission is (1 - omega) B; the m=0 systems split
  // [I,Q] and [U,V] according to paper Eq. (82).  Only m = 0 has a source,
  // so the other modes keep the zeros of the constructor.
  for (Index alpha = 0; alpha < 2; ++alpha) {
    const Index m = 0;
    for (Index l = 0; l < NLayers; ++l) {
      bool zero_emission = true;
      for (Index p = 0; p < Nscoeffs; ++p) zero_emission &= scaled_source_poly_coeffs[l, p].is_zero();
      // Exact conservative scattering has zero absorption emission and a
      // singular transport matrix.
      if (zero_emission) {
        source_collect[alpha, m, l] = rtepack::stokvec{};
        continue;
      }

      transport = 0.0;
      for (Index i = 0; i < NQuad; ++i) {
        for (Index so = 0; so < stokes_dimension; ++so) {
          const Index row = state_index(i, so);
          for (Index j = 0; j < NQuad; ++j) {
            const Numeric factor = -0.5 * omega_arr[l] * W[j % N] * inv_mu_arr[i];
            for (Index si = 0; si < stokes_dimension; ++si)
              transport[row, state_index(j, si)] += factor * phase_matrix[alpha, m, l, i, j][so, si];
          }
          transport[row, row] += inv_mu_arr[i];
        }
      }

      // One LU decomposition serves every polynomial degree
      transport_work = transport;
      lu_decompose_inplace(transport_work, solve_work);
      state_next = 0.0;
      for (Index p = Nscoeffs - 1; p >= 0; --p) {
        for (Index i = 0; i < NQuad; ++i) {
          for (Index s = 0; s < stokes_dimension; ++s) {
            const bool    active         = alpha == cosine_mode ? s < 2 : s >= 2;
            const Numeric q              = active ? scaled_source_poly_coeffs[l, p][s] : 0.0;
            state_rhs[state_index(i, s)] = inv_mu_arr[i] * q + source_coordinate_scale[l] *
                                                                   static_cast<Numeric>(p + 1) *
                                                                   state_next[state_index(i, s)];
          }
        }
        lu_solve_inplace(state_rhs, transport_work, solve_work);
        for (Index i = 0; i < NQuad; ++i)
          for (Index s = 0; s < stokes_dimension; ++s) {
            const Index state                    = state_index(i, s);
            source_collect[alpha, m, l, i, p][s] = state_rhs[state];
            state_next[state]                    = state_rhs[state];
          }
      }
    }
  }
}

void main_data::transmission() {
  // Complex modes are anchored at their nearest layer boundary.
}

void main_data::solve_for_coefs() {
  ARTS_TIME_REPORT

  check_beam_reflection(surface_refl, has_beam_source, mu0);

  const Index equation_count = NLayers * NState;

  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      // Reflection couples all Stokes components.
      if (m < NBDRF) {
        for (Index i = 0; i < N; ++i)
          for (Index so = 0; so < stokes_dimension; ++so)
            for (Index j = 0; j < N; ++j)
              for (Index si = 0; si < stokes_dimension; ++si)
                surface_operator[state_index(i, so), state_index(j, si)] = Constant::pi * (m == 0 ? 1.0 : 0.5) * W[j] *
                                                                           mu_arr[j] *
                                                                           surface_refl.brdf()[alpha, m, i, j][so, si];

        if (has_beam_source) {
          const Numeric           attenuation = std::exp(-tau_arr.back() / mu0);
          const rtepack::stokvec& beam        = beam_stokes;
          for (Index i = 0; i < N; ++i) {
            const rtepack::stokvec reflected = 0.5 * mu0 * attenuation * (surface_refl.beam_brdf()[alpha, m, i] * beam);
            for (Index s = 0; s < stokes_dimension; ++s) direct_reflection[state_index(i, s)] = reflected[s];
          }
        } else {
          direct_reflection = 0.0;
        }
      } else {
        surface_operator  = 0.0;
        direct_reflection = 0.0;
      }
      // The specular part reflects each downward stream into the upward stream of the same cosine, in every mode
      for (Index i = 0; i < static_cast<Index>(surface_refl.specular().size()); ++i) {
        const rtepack::muelmat& R = surface_refl.specular(i);
        for (Index so = 0; so < stokes_dimension; ++so)
          for (Index si = 0; si < stokes_dimension; ++si)
            surface_operator[state_index(i, so), state_index(i, si)] += R[so, si];
      }
      const Index   last   = NLayers - 1;
      const Numeric bottom = tau_arr.back();

      // The right-hand side.  Top boundary: prescribed downward diffuse Stokes field.
      for (Index i = 0; i < N; ++i) {
        for (Index s = 0; s < stokes_dimension; ++s) {
          const Index row   = state_index(i, s);
          const Index state = state_index(N + i, s);
          boundary_rhs[row] = boundary_down[alpha, m, i][s] - particular(alpha, m, 0, state, 0.0);
        }
      }

      // Layer-interface continuity, paper Eqs. (118)-(121).
      for (Index l = 0; l < NLayers - 1; ++l) {
        const Numeric tau  = tau_arr[l];
        const Index   row0 = NHalfState + l * NState;
        for (Index state = 0; state < NState; ++state)
          boundary_rhs[row0 + state] = particular(alpha, m, l + 1, state, tau) - particular(alpha, m, l, state, tau);
      }

      /* Lower boundary, including the polarized BRDF of paper Eq. (124), and the radiance beneath the surface
         transmitted with 1 - R: for a Planck radiance B e_I, Kirchhoff's emission B (e_I - R e_I). */
      for (Index i = 0; i < N; ++i) {
        for (Index s = 0; s < stokes_dimension; ++s) {
          const Index hrow      = state_index(i, s);
          const Index pos_state = state_index(i, s);
          Numeric     rhs_value =
              boundary_up[alpha, m, i][s] + direct_reflection[hrow] - particular(alpha, m, last, pos_state, bottom);
          for (Index j = 0; j < N; ++j)
            for (Index si = 0; si < stokes_dimension; ++si) {
              const Index hin        = state_index(j, si);
              const Index neg_state  = state_index(N + j, si);
              rhs_value             += surface_operator[hrow, hin] *
                                       (particular(alpha, m, last, neg_state, bottom) - boundary_up[alpha, m, j][si]);
            }
          boundary_rhs[equation_count - NHalfState + hrow] = rhs_value;
        }
      }

      // Nothing drives this mode, so its solution is zero and there is no system to solve
      if (stdr::all_of(boundary_rhs, [](const Numeric x) { return x == 0.0; })) {
        GC_collect[alpha, m] = Complex{};
        continue;
      }

      // Each mode's exponential at each layer's top (edge 0) and bottom (edge 1), as homogeneous() has it
      for (Index l = 0; l < NLayers; ++l) {
        for (Index e = 0; e < NState; ++e) {
          const Numeric anchor      = top_anchored[anchor_index(alpha, m, l, e)] ? layer_top(l) : tau_arr[l];
          mode_exponential[l, 0, e] = std::exp(K_collect[alpha, m, l, e] * (layer_top(l) - anchor));
          mode_exponential[l, 1, e] = std::exp(K_collect[alpha, m, l, e] * (tau_arr[l] - anchor));
        }
      }
      const auto mode = [&](const Index l, const Index state, const Index e, const Index edge) -> Complex {
        const auto [negative, positive] = conservative_pair_index[static_cast<std::size_t>(l)];
        if (alpha == cosine_mode and m == 0 and negative >= 0 and (e == negative or e == positive))
          return homogeneous(alpha, m, l, state, e, edge == 0 ? layer_top(l) : tau_arr[l]);
        return G_collect[alpha, m, l, state, e] * mode_exponential[l, edge, e];
      };

      /* The field is real, so the coefficients of a conjugate pair of modes are conjugate:
         c H + c* H* = 2 Re(c) Re H - 2 Im(c) Im H.  The member with the positive imaginary
         part represents the pair, its real part in its own column and its imaginary part in
         the partner's, which is skipped. */
      const auto partner = [&](const Index l, const Index e) {
        return conjugate_partner[static_cast<std::size_t>(anchor_index(alpha, m, l, e))];
      };
      const auto skipped = [&](const Index l, const Index e) {
        return partner(l, e) >= 0 and K_collect[alpha, m, l, e].imag() < 0.0;
      };
      const auto set = [&](const Index row, const Index l, const Index e, const Complex value) {
        const Index p = partner(l, e);
        if (p < 0) {
          boundary_lhs[row, l * NState + e] = value.real();
        } else {
          boundary_lhs[row, l * NState + e] = 2.0 * value.real();
          boundary_lhs[row, l * NState + p] = -2.0 * value.imag();
        }
      };

      // The band storage holds the previous LU factors
      boundary_lhs.zero();
      for (Index i = 0; i < N; ++i) {
        for (Index s = 0; s < stokes_dimension; ++s) {
          const Index row   = state_index(i, s);
          const Index state = state_index(N + i, s);
          for (Index e = 0; e < NState; ++e)
            if (not skipped(0, e)) set(row, 0, e, mode(0, state, e, 0));
        }
      }
      for (Index l = 0; l < NLayers - 1; ++l) {
        const Index row0 = NHalfState + l * NState;
        for (Index state = 0; state < NState; ++state) {
          const Index row = row0 + state;
          for (Index e = 0; e < NState; ++e) {
            if (not skipped(l, e)) set(row, l, e, mode(l, state, e, 1));
            if (not skipped(l + 1, e)) set(row, l + 1, e, -mode(l + 1, state, e, 0));
          }
        }
      }
      for (Index i = 0; i < N; ++i) {
        for (Index s = 0; s < stokes_dimension; ++s) {
          const Index hrow      = state_index(i, s);
          const Index row       = equation_count - NHalfState + hrow;
          const Index pos_state = state_index(i, s);
          for (Index e = 0; e < NState; ++e) {
            if (skipped(last, e)) continue;
            Complex value = mode(last, pos_state, e, 1);
            for (Index j = 0; j < N; ++j)
              for (Index si = 0; si < stokes_dimension; ++si)
                value -= surface_operator[hrow, state_index(j, si)] * mode(last, state_index(N + j, si), e, 1);
            set(row, last, e, value);
          }
        }
      }

      const int info = boundary_lhs.solve(boundary_rhs);
      if (info != 0)
        arts::user_error(
            "VDISORT boundary solve failed for alpha={}, Fourier mode={} (LAPACK info={})", alpha, m, info);

      // Back to the complex coefficients of every mode.  A partner anchored at the other layer edge (possible
      // only for Re k = 0) gets the conjugate times its exponential between the two anchors.
      for (Index l = 0; l < NLayers; ++l) {
        for (Index e = 0; e < NState; ++e) {
          const Index p = partner(l, e);
          if (p < 0) {
            GC_collect[alpha, m, l, e] = boundary_rhs[l * NState + e];
          } else if (not skipped(l, e)) {
            const Complex c            = {boundary_rhs[l * NState + e], boundary_rhs[l * NState + p]};
            GC_collect[alpha, m, l, e] = c;
            const Numeric anchor_e     = top_anchored[anchor_index(alpha, m, l, e)] ? layer_top(l) : tau_arr[l];
            const Numeric anchor_p     = top_anchored[anchor_index(alpha, m, l, p)] ? layer_top(l) : tau_arr[l];
            GC_collect[alpha, m, l, p] =
                anchor_e == anchor_p
                    ? std::conj(c)
                    : std::conj(c) * std::exp(std::conj(K_collect[alpha, m, l, e]) * (anchor_p - anchor_e));
          }
        }
      }
    }
  }
}

void main_data::rad_field() {
  // Both combined modes are retained; paper Eq. (78) reconstructs [I,Q,U,V].
  for (Index l = 0; l < NLayers; ++l) {
    combined_field(layer_field, tau_arr[l]);
    for (Index alpha = 0; alpha < 2; ++alpha)
      for (Index m = 0; m < NFourier; ++m)
        for (Index i = 0; i < NQuad; ++i)
          for (Index s = 0; s < stokes_dimension; ++s)
            um[l, alpha, m, i][s] = layer_field[alpha * NFourier + m, state_index(i, s)];
  }
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
  if ((out.shape() != std::array<Index, 2>{2 * NFourier, NState}))
    arts::user_error("Combined field has shape {:B,}, expected [{}, {}]", out.shape(), 2 * NFourier, NState);

  // G_collect is contiguous in the innermost eigenmode index; centered pairs
  // occupy stored X and R columns.
  ComplexVector modal_amplitude(NState);
  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      for (Index eigen = 0; eigen < NState; ++eigen) {
        const Numeric anchor = top_anchored[anchor_index(alpha, m, layer, eigen)] ? layer_top(layer) : tau_arr[layer];
        modal_amplitude[eigen] =
            GC_collect[alpha, m, layer, eigen] * std::exp(K_collect[alpha, m, layer, eigen] * (tau - anchor));
      }

      const auto [negative, positive] = conservative_pair_index[static_cast<std::size_t>(layer)];
      if (alpha == cosine_mode and m == 0 and negative >= 0) {
        const Numeric center     = 0.5 * (layer_top(layer) + tau_arr[layer]);
        const Numeric s          = tau - center;
        const Numeric kappa      = conservative_pair_kappa[layer];
        const auto    basis      = dc::centered_pair(kappa, s);
        const auto    amplitudes = dc::centered_pair_amplitudes(
            GC_collect[alpha, m, layer, negative], GC_collect[alpha, m, layer, positive], basis);
        modal_amplitude[negative] = amplitudes[0];
        modal_amplitude[positive] = amplitudes[1];
      }

      for (Index state = 0; state < NState; ++state) {
        Complex value = particular(alpha, m, layer, state, tau);
        for (Index eigen = 0; eigen < NState; ++eigen)
          value += G_collect[alpha, m, layer, state, eigen] * modal_amplitude[eigen];
        const Numeric scale = 1.0 + std::abs(value.real());
        // Conservative scalar-limit eigenpairs may leave small imaginary
        // roundoff although the reconstructed physical field is real.
        if (std::abs(value.imag()) > 2e-8 * scale)
          arts::user_error("VDISORT left an uncancelled imaginary radiance {} at alpha={}, m={}, state={}, tau={}",
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
  if (!std::isfinite(phi) or phi < 0.0 or phi >= Constant::two_pi)
    arts::user_error("phi must be finite and in [0, 2*pi), got {}", phi);
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
      // U and V swap combined cosine/sine roles in paper Eq. (78).
      for (Index stokes = 2; stokes < 4; ++stokes)
        data.intensities[i][stokes] +=
            combined[NFourier + m, state_index(i, stokes)] * c + combined[m, state_index(i, stokes)] * s;
    }
  }
}

void main_data::u_user(user_u_data&                    data,
                       const Numeric                   tau,
                       const Numeric                   phi,
                       const ConstVectorView&          user_mu,
                       const phase_matrix_data&        user_phase_matrix,
                       const beam_phase_matrix_data&   user_beam_phase_matrix,
                       const rtepack::stokvec_tensor3& user_boundary,
                       const surface_reflection&       user_surface) const {
  rtepack::stokvec_tensor3 result(1, 1, user_mu.size());
  ungridded_u_user(result,
                   AscendingGrid{tau},
                   Vector{phi},
                   user_mu,
                   user_phase_matrix,
                   user_beam_phase_matrix,
                   user_boundary,
                   user_surface);
  data.intensities.resize(user_mu.size());
  for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user) data.intensities[user] = result[0, 0, user];
}

void main_data::check_correction_compatibility(const delta_m_correction_cache& correction,
                                               const bool                      quadrature_angles) const {
  const auto& correction_tau = correction.scaled_tau();
  if (correction_tau.size() != tau_arr.size())
    arts::user_error(
        "The delta-M correction has {} scaled layers, but VDISORT has {}", correction_tau.size(), tau_arr.size());
  for (Index layer = 0; layer < NLayers; ++layer) {
    const Numeric scale = std::max({1.0, std::abs(correction_tau[layer]), std::abs(tau_arr[layer])});
    if (std::abs(correction_tau[layer] - tau_arr[layer]) > 32.0 * std::numeric_limits<Numeric>::epsilon() * scale)
      arts::user_error("The delta-M correction scaled depth {} in layer {} differs from VDISORT depth {}",
                       correction_tau[layer],
                       layer,
                       tau_arr[layer]);
  }

  if (not quadrature_angles) return;
  const auto& correction_mu = correction.user_mu();
  if (correction_mu.size() != mu_arr.size())
    arts::user_error("The delta-M correction has {} directions, but VDISORT has {} quadrature streams",
                     correction_mu.size(),
                     mu_arr.size());
  for (Index stream = 0; stream < NQuad; ++stream) {
    const Numeric scale = std::max({1.0, std::abs(correction_mu[stream]), std::abs(mu_arr[stream])});
    if (std::abs(correction_mu[stream] - mu_arr[stream]) > 32.0 * std::numeric_limits<Numeric>::epsilon() * scale)
      arts::user_error("The delta-M correction cosine {} at stream {} differs from VDISORT cosine {}",
                       correction_mu[stream],
                       stream,
                       mu_arr[stream]);
  }
}

void main_data::TMS(rtepack::stokvec_vector&        tms,
                    const Numeric                   physical_tau,
                    const Index                     phi_index,
                    const delta_m_correction_cache& correction) const {
  check_correction_compatibility(correction, true);
  correction.TMS(tms, physical_tau, phi_index, beam_phase_matrix);
}

void main_data::IMS(rtepack::stokvec_vector&        ims,
                    const Numeric                   physical_tau,
                    const Index                     phi_index,
                    const delta_m_correction_cache& correction,
                    const ims_convention            convention) const {
  check_correction_compatibility(correction, false);
  correction.IMS(ims, physical_tau, phi_index, convention);
}

void main_data::gridded_TMS(rtepack::stokvec_tensor3_view tms, const delta_m_correction_cache& correction) const {
  check_correction_compatibility(correction, true);
  correction.gridded_TMS(tms, beam_phase_matrix);
}

void main_data::gridded_IMS(rtepack::stokvec_tensor3_view   ims,
                            const delta_m_correction_cache& correction,
                            const ims_convention            convention) const {
  check_correction_compatibility(correction, false);
  correction.gridded_IMS(ims, convention);
}

void main_data::u_user_corr(user_u_data&                    data,
                            rtepack::stokvec_vector&        ims,
                            rtepack::stokvec_vector&        tms,
                            const Numeric                   physical_tau,
                            const Index                     phi_index,
                            const delta_m_correction_cache& correction,
                            const phase_matrix_data&        user_phase_matrix,
                            const beam_phase_matrix_data&   user_beam_phase_matrix,
                            const surface_reflection&       user_surface,
                            const ims_convention            convention) const {
  ARTS_TIME_REPORT

  check_correction_compatibility(correction, false);
  correction.TMS(tms, physical_tau, phi_index, user_beam_phase_matrix);
  correction.IMS(ims, physical_tau, phi_index, convention);
  u_user(data,
         correction.scaled_tau(physical_tau),
         correction.phi()[phi_index],
         correction.user_mu(),
         user_phase_matrix,
         user_beam_phase_matrix,
         {},
         user_surface);
  for (Index user = 0; user < static_cast<Index>(data.intensities.size()); ++user)
    data.intensities[user] += tms[user] + ims[user];
}

void main_data::user_fourier_modes(ComplexTensor4&                 modes,
                                   const AscendingGrid&            tau,
                                   const ConstVectorView&          user_mu,
                                   const phase_matrix_data&        user_phase_matrix,
                                   const beam_phase_matrix_data&   user_beam_phase_matrix,
                                   const rtepack::stokvec_tensor3& user_boundary,
                                   const surface_reflection&       user_surface) const {
  ARTS_TIME_REPORT

  const Index ntau  = static_cast<Index>(tau.size());
  const Index nuser = static_cast<Index>(user_mu.size());
  if (tau.empty()) arts::user_error("At least one optical depth is required for user-angle radiances");
  if (tau.front() < 0.0 or tau.back() > tau_arr.back())
    arts::user_error("User-angle optical depths must be in [0, {}], got {:B,}", tau_arr.back(), tau);
  if (stdr::any_of(user_mu, [](const Numeric mu) { return !std::isfinite(mu) or mu == 0.0 or std::abs(mu) > 1.0; }))
    arts::user_error("User polar-angle cosines must be finite, nonzero, and in [-1, 1], got {:B,}", user_mu);
  if (not user_boundary.empty() and user_boundary.shape() != (std::array<Index, 3>{2, NFourier, nuser}))
    arts::user_error("The user-direction boundary radiances have shape {:B,}, expected [2, {}, {}] Stokes vectors",
                     user_boundary.shape(),
                     NFourier,
                     nuser);
  const std::array<Index, 5> expected_phase_shape{2, NFourier, NLayers, nuser, NQuad};
  const std::array<Index, 4> expected_beam_shape{2, NFourier, NLayers, nuser};
  if (user_phase_matrix.shape() != expected_phase_shape)
    arts::user_error("User phase matrices have shape {:B,}, expected [2, {}, {}, {}, {}] Mueller blocks",
                     user_phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     nuser,
                     NQuad);
  if (has_beam_source and user_beam_phase_matrix.shape() != expected_beam_shape)
    arts::user_error("User beam phase matrices have shape {:B,}, expected [2, {}, {}, {}] Mueller blocks",
                     user_beam_phase_matrix.shape(),
                     NFourier,
                     NLayers,
                     nuser);

  // The surface reflects into the upward user directions, so it is needed there
  const bool specular = surface_refl.has_specular();
  if ((NBDRF > 0 or specular) and stdr::any_of(user_mu, Cmp::gt(0.0))) {
    if ((user_surface.brdf().shape() != std::array<Index, 4>{2, NBDRF, nuser, N}) or
        (user_surface.beam_brdf().shape() != std::array<Index, 3>{2, NBDRF, nuser}) or
        (specular and static_cast<Index>(user_surface.specular().size()) != nuser))
      arts::user_error(
          "The surface reflects into the upward user directions, so their reflection must be given: shapes "
          "[2, {}, {}, {}], [2, {}, {}] and [{}] (specular part) are needed, got {:B,}, {:B,} and {:B,}",
          NBDRF,
          nuser,
          N,
          NBDRF,
          nuser,
          specular ? nuser : 0,
          user_surface.brdf().shape(),
          user_surface.beam_brdf().shape(),
          user_surface.specular().shape());
    check_beam_reflection(user_surface, has_beam_source, mu0);
  }

  /* A specular surface part reflects the downward radiance at the surface into the upward direction of the same
     cosine.  So the upward boundary value of an upward user direction mu includes R(mu) I_down(-mu) at the
     surface, from the user direction -mu, whose phase matrices only the caller has: every upward user direction
     must come with its downward partner. */
  std::vector<Index> partner(nuser, -1);  // the downward partner of an upward user direction, as an index of below
  ComplexTensor4     partner_modes;       // [1, partners, 2 NFourier, Stokes] at the surface
  if (specular) {
    std::vector<Index> partners;
    for (Index iu = 0; iu < nuser; ++iu) {
      if (user_mu[iu] < 0.0) continue;
      Index found = -1;
      for (Index id = 0; id < nuser; ++id)
        if (user_mu[id] < 0.0 and std::abs(user_mu[id] + user_mu[iu]) <= 64.0 * std::numeric_limits<Numeric>::epsilon())
          found = id;
      if (found < 0)
        arts::user_error(
            "The surface reflects specularly, so the upward radiance at the user cosine {} reflects the "
            "downward radiance at {} at the surface: include {} in the user cosines, with its phase "
            "matrices",
            user_mu[iu],
            -user_mu[iu],
            -user_mu[iu]);
      partner[iu] = static_cast<Index>(partners.size());
      partners.push_back(found);
    }
    if (not partners.empty()) {
      const Index            np = static_cast<Index>(partners.size());
      Vector                 partner_mu(np);
      phase_matrix_data      partner_phase(2, NFourier, NLayers, np, NQuad);
      beam_phase_matrix_data partner_beam;
      if (has_beam_source) partner_beam.resize(2, NFourier, NLayers, np);
      rtepack::stokvec_tensor3 partner_boundary;
      if (not user_boundary.empty()) partner_boundary.resize(2, NFourier, np);
      for (Index p = 0; p < np; ++p) {
        partner_mu[p] = user_mu[partners[p]];
        for (Index alpha = 0; alpha < 2; ++alpha)
          for (Index m = 0; m < NFourier; ++m)
            for (Index layer = 0; layer < NLayers; ++layer) {
              for (Index j = 0; j < NQuad; ++j)
                partner_phase[alpha, m, layer, p, j] = user_phase_matrix[alpha, m, layer, partners[p], j];
              if (has_beam_source)
                partner_beam[alpha, m, layer, p] = user_beam_phase_matrix[alpha, m, layer, partners[p]];
            }
        if (not user_boundary.empty())
          for (Index alpha = 0; alpha < 2; ++alpha)
            for (Index m = 0; m < NFourier; ++m) partner_boundary[alpha, m, p] = user_boundary[alpha, m, partners[p]];
      }
      user_fourier_modes(partner_modes,
                         AscendingGrid{Vector{tau_arr.back()}},
                         partner_mu,
                         partner_phase,
                         partner_beam,
                         partner_boundary,
                         surface_reflection{});
    }
  }

  const auto interpolate_boundary =
      [&](const rtepack::stokvec_tensor3& boundary, const Index alpha, const Index m, const Numeric abs_mu) {
        return dc::barycentric_interpolate(
            mu_arr[Range{0, N}], half_range_barycentric_weights, boundary[alpha, m], abs_mu);
      };

  modes.resize(ntau, nuser, 2 * NFourier, stokes_dimension);
  modes = Complex{0.0};

  const Numeric atmosphere_bottom = tau_arr.back();

  for (Index iu = 0; iu < nuser; ++iu) {
    const Numeric mu       = user_mu[iu];
    const Numeric abs_mu   = std::abs(mu);
    const bool    downward = mu < 0.0;

    for (Index alpha = 0; alpha < 2; ++alpha) {
      for (Index m = 0; m < NFourier; ++m) {
        // The boundary radiance where the ray starts: given, or interpolated from the streams
        rtepack::stokvec mode = user_boundary.empty()
                                    ? interpolate_boundary(downward ? boundary_down : boundary_up, alpha, m, abs_mu)
                                    : user_boundary[alpha, m, iu];

        const Numeric boundary_tau = downward ? 0.0 : tau_arr.back();
        if (not downward and (m < NBDRF or specular)) {
          // The reflection of the downward field, and the radiance beneath the surface transmitted with 1 - R
          const rtepack::stokvec beneath = mode;
          if (m < NBDRF) {
            for (Index j = 0; j < N; ++j) {
              mode += Constant::pi * (m == 0 ? 1.0 : 0.5) * W[j] * mu_arr[j] *
                      (user_surface.brdf()[alpha, m, iu, j] *
                       (um[NLayers - 1, alpha, m, N + j] - boundary_up[alpha, m, j]));
            }
            if (has_beam_source)
              mode += 0.5 * mu0 * std::exp(-atmosphere_bottom / mu0) *
                      (user_surface.beam_brdf()[alpha, m, iu] * beam_stokes);
          }
          if (specular) {
            rtepack::stokvec down;
            for (Index s = 0; s < stokes_dimension; ++s)
              down[s] = partner_modes[0, partner[iu], alpha * NFourier + m, s].real();
            mode += user_surface.specular(iu) * (down - beneath);
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

          const auto [pair_negative, pair_positive] = conservative_pair_index[static_cast<std::size_t>(layer)];
          const bool has_pair                       = alpha == cosine_mode and m == 0 and pair_negative >= 0;
          for (Index eigen = 0; eigen < NState; ++eigen) {
            if (has_pair and (eigen == pair_negative or eigen == pair_positive)) continue;
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
              const Complex integral = dc::user_angle_exponential_integral(
                  K_collect[alpha, m, layer, eigen], anchor, lower, upper, tau[t], abs_mu, downward);
              for (Index s = 0; s < stokes_dimension; ++s) modes[t, iu, mode, s] += coefficient * source[s] * integral;
            }
          }

          if (has_pair) {
            std::array<Complex, stokes_dimension> source_x{};
            std::array<Complex, stokes_dimension> source_r{};
            for (Index j = 0; j < NQuad; ++j) {
              const Numeric factor = 0.5 * omega_arr[layer] * W[j % N];
              const auto&   phase  = user_phase_matrix[alpha, m, layer, iu, j];
              for (Index so = 0; so < stokes_dimension; ++so) {
                for (Index si = 0; si < stokes_dimension; ++si) {
                  source_x[so] +=
                      factor * phase[so, si] * G_collect[alpha, m, layer, state_index(j, si), pair_negative];
                  source_r[so] +=
                      factor * phase[so, si] * G_collect[alpha, m, layer, state_index(j, si), pair_positive];
                }
              }
            }

            const Complex c0     = GC_collect[alpha, m, layer, pair_negative];
            const Complex c1     = GC_collect[alpha, m, layer, pair_positive];
            const Numeric center = 0.5 * (layer_top(layer) + tau_arr[layer]);
            const Numeric kappa  = conservative_pair_kappa[layer];
            for (Index t = 0; t < ntau; ++t) {
              const auto [lower, upper] = limits(layer, tau[t], downward);
              if (upper <= lower) continue;

              const auto [integral_c, integral_s] =
                  dc::centered_pair_integrals(kappa, center, lower, upper, tau[t], abs_mu, downward);
              for (Index s = 0; s < stokes_dimension; ++s) {
                const Complex ac       = source_x[s] * c0 + source_r[s] * c1;
                const Complex as       = source_x[s] * c1 + kappa * kappa * source_r[s] * c0;
                modes[t, iu, mode, s] += ac * integral_c + as * integral_s;
              }
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
              const Numeric integral =
                  dc::user_angle_exponential_integral(-1.0 / mu0, 0.0, lower, upper, tau[t], abs_mu, downward);
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

void main_data::ungridded_u_user(rtepack::stokvec_tensor3_view   out,
                                 const AscendingGrid&            tau,
                                 const Vector&                   phi,
                                 const ConstVectorView&          user_mu,
                                 const phase_matrix_data&        user_phase_matrix,
                                 const beam_phase_matrix_data&   user_beam_phase_matrix,
                                 const rtepack::stokvec_tensor3& user_boundary,
                                 const surface_reflection&       user_surface) const {
  ARTS_TIME_REPORT

  const std::array<Index, 3> expected_shape{
      static_cast<Index>(tau.size()), static_cast<Index>(phi.size()), static_cast<Index>(user_mu.size())};
  if (out.shape() != expected_shape)
    arts::user_error("User-angle output has shape {:B,}, expected [{}, {}, {}] Stokes blocks",
                     out.shape(),
                     tau.size(),
                     phi.size(),
                     user_mu.size());
  if (stdr::any_of(
          phi, [](const Numeric value) { return !std::isfinite(value) or value < 0.0 or value >= Constant::two_pi; }))
    arts::user_error("User azimuths must be finite and in [0, 2*pi), got {:B,}", phi);

  ComplexTensor4 modes;
  user_fourier_modes(modes, tau, user_mu, user_phase_matrix, user_beam_phase_matrix, user_boundary, user_surface);
  out = rtepack::stokvec{};
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t) {
    for (Index iu = 0; iu < static_cast<Index>(user_mu.size()); ++iu) {
      for (Index m = 0; m < NFourier; ++m) {
        for (Index alpha = 0; alpha < 2; ++alpha)
          for (Index stokes = 0; stokes < stokes_dimension; ++stokes) {
            const Complex value = modes[t, iu, alpha * NFourier + m, stokes];
            const Numeric scale = 1.0 + std::abs(value.real());
            if (std::abs(value.imag()) > 2e-8 * scale)
              arts::user_error(
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

void main_data::u_corr(u_data&                         data,
                       rtepack::stokvec_vector&        ims,
                       rtepack::stokvec_vector&        tms,
                       const Numeric                   physical_tau,
                       const Index                     phi_index,
                       const delta_m_correction_cache& correction,
                       const ims_convention            convention) const {
  ARTS_TIME_REPORT

  check_correction_compatibility(correction, true);
  correction.TMS(tms, physical_tau, phi_index, beam_phase_matrix);
  correction.IMS(ims, physical_tau, phi_index, convention);
  u(data, correction.scaled_tau(physical_tau), correction.phi()[phi_index]);
  for (Index stream = 0; stream < NQuad; ++stream) data.intensities[stream] += tms[stream] + ims[stream];
}

flux_values main_data::flux(flux_data& data, const Numeric tau) const {
  u0_data field;
  u0(field, tau);
  data.u0 = std::move(field.u0);

  const auto diffuse = dc::integrate_diffuse(mu_arr[Range{0, N}], W, [&](const Index i) { return data.u0[i].I(); });

  const Numeric direct         = has_beam_source ? dc::direct_beam_flux(beam_stokes.I(), mu0, tau) : 0.0;
  Numeric       mean_intensity = diffuse.mean_intensity;
  if (has_beam_source) mean_intensity += dc::direct_beam_radiance(beam_stokes.I(), mu0, tau) / (4.0 * Constant::pi);

  const Index   layer      = tau_index(tau);
  const Numeric source_tau = std::fma(source_coordinate_scale[layer], tau, source_coordinate_offset[layer]);
  const Numeric source     = dc::horner_polynomial(
      Nscoeffs, source_tau, [&](const Index coefficient) { return source_poly_coeffs[layer, coefficient].I(); });

  return {.up           = Constant::two_pi * diffuse.upward,
          .down_diffuse = Constant::two_pi * diffuse.downward,
          .down_direct  = direct,
          .dfdt         = 4.0 * Constant::pi * (1.0 - omega_arr[layer]) * (mean_intensity - source)};
}

void main_data::gridded_u(Tensor4View out, const Vector& phi) const {
  if ((out.shape() != std::array<Index, 4>{NLayers, static_cast<Index>(phi.size()), NQuad, stokes_dimension}))
    arts::user_error(
        "gridded_u output has shape {:B,}, expected [{}, {}, {}, 4]", out.shape(), NLayers, phi.size(), NQuad);

  const Index nphi = static_cast<Index>(phi.size());
  Matrix      cosine(nphi, NFourier);
  Matrix      sine(nphi, NFourier);
  for (Index p = 0; p < nphi; ++p) {
    if (!std::isfinite(phi[p]) or phi[p] < 0.0 or phi[p] >= Constant::two_pi)
      arts::user_error("phi must be finite and in [0, 2*pi), got {}", phi[p]);
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
          // U and V swap combined cosine/sine roles in paper Eq. (78).
          for (Index stokes = 2; stokes < stokes_dimension; ++stokes)
            out[l, p, i, stokes] += um[l, sine_mode, m, i][stokes] * c + um[l, cosine_mode, m, i][stokes] * s;
        }
      }
}

void main_data::gridded_u_corr(Tensor4View                     out,
                               rtepack::stokvec_tensor3_view   tms,
                               rtepack::stokvec_tensor3_view   ims,
                               const delta_m_correction_cache& correction,
                               const ims_convention            convention) const {
  ARTS_TIME_REPORT

  check_correction_compatibility(correction, true);
  correction.gridded_TMS(tms, beam_phase_matrix);
  correction.gridded_IMS(ims, convention);
  gridded_u(out, correction.phi());
  for (Index layer = 0; layer < NLayers; ++layer)
    for (Index p = 0; p < static_cast<Index>(correction.phi().size()); ++p)
      for (Index stream = 0; stream < NQuad; ++stream)
        for (Index stokes = 0; stokes < stokes_dimension; ++stokes)
          out[layer, p, stream, stokes] += tms[layer, p, stream][stokes] + ims[layer, p, stream][stokes];
}

void main_data::gridded_flux(VectorView up, VectorView down_diffuse, VectorView down_direct, VectorView dfdt) const {
  if (up.size() != static_cast<Size>(NLayers) or down_diffuse.size() != static_cast<Size>(NLayers) or
      down_direct.size() != static_cast<Size>(NLayers) or dfdt.size() != static_cast<Size>(NLayers))
    arts::user_error("All gridded flux outputs must have size {}", NLayers);
  for (Index l = 0; l < NLayers; ++l) {
    const auto diffuse =
        dc::integrate_diffuse(mu_arr[Range{0, N}], W, [&](const Index i) { return um[l, cosine_mode, 0, i].I(); });
    up[l]           = Constant::two_pi * diffuse.upward;
    down_diffuse[l] = Constant::two_pi * diffuse.downward;
    down_direct[l]  = has_beam_source ? dc::direct_beam_flux(beam_stokes.I(), mu0, tau_arr[l]) : 0.0;

    Numeric mean_intensity = diffuse.mean_intensity;
    if (has_beam_source)
      mean_intensity += dc::direct_beam_radiance(beam_stokes.I(), mu0, tau_arr[l]) / (4.0 * Constant::pi);
    const Numeric source_tau = std::fma(source_coordinate_scale[l], tau_arr[l], source_coordinate_offset[l]);
    const Numeric source     = dc::horner_polynomial(
        Nscoeffs, source_tau, [&](const Index coefficient) { return source_poly_coeffs[l, coefficient].I(); });
    dfdt[l] = 4.0 * Constant::pi * (1.0 - omega_arr[l]) * (mean_intensity - source);
  }
}

void main_data::ungridded_u(Tensor4View out, const AscendingGrid& tau, const Vector& phi) const {
  if ((out.shape() !=
       std::array<Index, 4>{static_cast<Index>(tau.size()), static_cast<Index>(phi.size()), NQuad, stokes_dimension}))
    arts::user_error(
        "ungridded_u output has shape {:B,}, expected [{}, {}, {}, 4]", out.shape(), tau.size(), phi.size(), NQuad);

  const Index nphi = static_cast<Index>(phi.size());
  if (tau.empty() or nphi == 0) return;

  Matrix cosine(nphi, NFourier);
  Matrix sine(nphi, NFourier);
  for (Index p = 0; p < nphi; ++p) {
    if (!std::isfinite(phi[p]) or phi[p] < 0.0 or phi[p] >= Constant::two_pi)
      arts::user_error("phi must be finite and in [0, 2*pi), got {}", phi[p]);
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

void main_data::ungridded_flux(
    VectorView up, VectorView down_diffuse, VectorView down_direct, VectorView dfdt, const AscendingGrid& tau) const {
  if (up.size() != tau.size() or down_diffuse.size() != tau.size() or down_direct.size() != tau.size() or
      dfdt.size() != tau.size())
    arts::user_error("All ungridded flux outputs must have the same size as tau ({})", tau.size());
  Matrix combined(2 * NFourier, NState);
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t) {
    combined_field(combined, tau[t]);
    const auto diffuse =
        dc::integrate_diffuse(mu_arr[Range{0, N}], W, [&](const Index i) { return combined[0, state_index(i, 0)]; });
    up[t]           = Constant::two_pi * diffuse.upward;
    down_diffuse[t] = Constant::two_pi * diffuse.downward;
    down_direct[t]  = has_beam_source ? dc::direct_beam_flux(beam_stokes.I(), mu0, tau[t]) : 0.0;

    Numeric mean_intensity = diffuse.mean_intensity;
    if (has_beam_source)
      mean_intensity += dc::direct_beam_radiance(beam_stokes.I(), mu0, tau[t]) / (4.0 * Constant::pi);
    const Index   layer      = tau_index(tau[t]);
    const Numeric source_tau = std::fma(source_coordinate_scale[layer], tau[t], source_coordinate_offset[layer]);
    const Numeric source     = dc::horner_polynomial(
        Nscoeffs, source_tau, [&](const Index coefficient) { return source_poly_coeffs[layer, coefficient].I(); });
    dfdt[t] = 4.0 * Constant::pi * (1.0 - omega_arr[layer]) * (mean_intensity - source);
  }
}

rtepack::stokvec_tensor3_const_view main_data::layer_um(const Size layer) const {
  if (layer >= static_cast<Size>(NLayers)) arts::user_error("Layer {} is out of range [0, {})", layer, NLayers);
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

#pragma once

#include <arts_constants.h>
#include <legendre.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

#include "disort.h"
#include "vdisort.h"

/**
 * VDISORT SCALAR TEST PORT BEGIN
 *
 * Compatibility surface used only by vdisort-test.cpp.  It deliberately has
 * the small part of disort::main_data's API exercised by the scalar DISORT
 * suite, while every transport calculation is performed by
 * vdisort::main_data with only the Stokes-I component populated.
 *
 * Keeping this adapter separate lets vdisort-test.cpp include the original
 * scalar case bodies and reference arrays without editing or copying them.
 */
namespace vdisort_scalar_test {
using BDRF = ::disort::BDRF;

struct u_data {
  Vector intensities;
};

struct u0_data {
  Vector u0;
};

struct tms_data {};
struct flux_data {};

class main_data {
  Index         nquad_{};
  Index         nleg_{};
  Index         nfourier_{};
  Index         nlayers_{};
  Index         n_{};
  AscendingGrid tau_{};
  Vector        omega_{};
  Vector        f_{};
  Matrix        legendre_{};
  Vector        scale_tau_{};
  Vector        scaled_tau_with_zero_{};
  Vector        scaled_omega_{};
  Matrix        weighted_legendre_{};
  Matrix        weighted_scaled_legendre_{};
  Numeric       mu0_{};
  Numeric       i0_{};
  Numeric       i0_original_{};
  Numeric       phi0_{};
  Numeric       omega_average_{};
  Numeric       f_average_{};
  Numeric       scaled_mu0_{};
  Vector        legendre_residue_average_{};

  std::unique_ptr<vdisort::main_data> model_{};

  [[nodiscard]] Index layer_index(const Numeric tau) const {
    ARTS_USER_ERROR_IF(tau < 0.0 or tau > tau_.back(), "tau ({}) must be in [0, {}]", tau, tau_.back());
    const Index layer = std::distance(tau_.begin(), std::ranges::lower_bound(tau_, tau));
    return std::min(layer, nlayers_ - 1);
  }

  [[nodiscard]] Numeric scaled_tau(const Numeric tau) const {
    const Index layer = layer_index(tau);
    return scaled_tau_with_zero_[layer + 1] - (tau_[layer] - tau) * scale_tau_[layer];
  }

  [[nodiscard]] static Numeric scattering_angle_cosine(const Numeric mu,
                                                       const Numeric phi,
                                                       const Numeric incident_mu,
                                                       const Numeric incident_phi) {
    return mu * incident_mu +
           std::sqrt(1.0 - incident_mu * incident_mu) * std::cos(incident_phi - phi) * std::sqrt(1.0 - mu * mu);
  }

  void scalar_u(Vector& out, const Numeric tau, const Numeric phi) const {
    vdisort::u_data field;
    model_->u(field, scaled_tau(tau), phi);
    out.resize(nquad_);
    for (Index stream = 0; stream < nquad_; ++stream) out[stream] = i0_original_ * field.intensities[stream, 0];
  }

  void truncation_correction(Vector& out, const Numeric tau, const Numeric phi) const {
    out.resize(nquad_);
    out = 0.0;
    if (i0_ <= 0.0) return;

    const Index   layer         = layer_index(tau);
    const Numeric tau_scaled    = scaled_tau(tau);
    const Numeric layer_top     = scaled_tau_with_zero_[layer];
    const Numeric layer_bottom  = scaled_tau_with_zero_[layer + 1];
    const Vector& quadrature_mu = model_->mu();

    for (Index stream = 0; stream < nquad_; ++stream) {
      const Numeric nu          = scattering_angle_cosine(quadrature_mu[stream], phi, -mu0_, phi0_);
      const Numeric p_true      = Legendre::legendre_sum(weighted_legendre_[layer], nu);
      const Numeric p_truncated = Legendre::legendre_sum(weighted_scaled_legendre_[layer], nu);
      const Numeric b = scaled_omega_[layer] * i0_ / (4.0 * Constant::pi) * (mu0_ / (mu0_ + quadrature_mu[stream])) *
                        (p_true / (1.0 - f_[layer]) - p_truncated);

      Numeric path_difference;
      if (stream < n_) {
        path_difference = std::exp(-tau_scaled / mu0_) -
                          std::exp((tau_scaled - layer_bottom) / quadrature_mu[stream] - layer_bottom / mu0_);
      } else {
        const Numeric abs_mu = -quadrature_mu[stream];
        path_difference = std::exp(-tau_scaled / mu0_) - std::exp((layer_top - tau_scaled) / abs_mu - layer_top / mu0_);
      }
      out[stream] = b * path_difference;
    }
  }

  void improved_multiple_scattering_correction(Vector& out, const Numeric tau, const Numeric phi) const {
    out.resize(n_);
    out = 0.0;
    if (i0_ <= 0.0 or f_average_ == 0.0) return;

    const Vector& quadrature_mu = model_->mu();
    for (Index i = 0; i < n_; ++i) {
      const Numeric nu  = scattering_angle_cosine(quadrature_mu[i + n_], phi, -mu0_, phi0_);
      const Numeric x   = 1.0 / quadrature_mu[i] - 1.0 / scaled_mu0_;
      const Numeric chi = 1.0 / (quadrature_mu[i] * scaled_mu0_ * x) *
                          ((tau - 1.0 / x) * std::exp(-tau / scaled_mu0_) + std::exp(-tau / quadrature_mu[i]) / x);
      out[i] = i0_ / (4.0 * Constant::pi) * Math::pow2(omega_average_ * f_average_) /
               (1.0 - omega_average_ * f_average_) * Legendre::legendre_sum(legendre_residue_average_, nu) * chi;
    }
  }

 public:
  main_data(const Index       nquad,
            const Index       nleg,
            const Index       nfourier,
            AscendingGrid     tau,
            Vector            omega,
            Matrix            legendre,
            Matrix            boundary_up,
            Matrix            boundary_down,
            Vector            f,
            Matrix            source_poly,
            std::vector<BDRF> brdf,
            const Numeric     mu0,
            const Numeric     beam,
            const Numeric     phi0)
      : nquad_(nquad),
        nleg_(nleg),
        nfourier_(nfourier),
        nlayers_(static_cast<Index>(tau.size())),
        n_(nquad / 2),
        tau_(std::move(tau)),
        omega_(std::move(omega)),
        f_(std::move(f)),
        legendre_(std::move(legendre)),
        scale_tau_(nlayers_),
        scaled_tau_with_zero_(nlayers_ + 1, 0.0),
        scaled_omega_(nlayers_),
        weighted_legendre_(nlayers_, legendre_.ncols(), 0.0),
        weighted_scaled_legendre_(nlayers_, nleg_, 0.0),
        mu0_(mu0),
        i0_(beam),
        i0_original_(1.0),
        phi0_(phi0),
        legendre_residue_average_(legendre_.ncols(), 0.0) {
    ARTS_USER_ERROR_IF(nquad_ <= 0 or nquad_ % 2 != 0, "NQuad must be a positive even number, got {}", nquad_);
    ARTS_USER_ERROR_IF(
        nleg_ > legendre_.ncols(), "NLeg ({}) exceeds the {} supplied Legendre coefficients", nleg_, legendre_.ncols());

    if (f_.empty()) {
      f_.resize(nlayers_);
      f_ = 0.0;
    }

    bool has_source   = source_poly.ncols() != 0;
    bool has_boundary = false;
    for (auto value = boundary_up.elem_begin(); value != boundary_up.elem_end(); ++value)
      if (*value != 0.0) has_boundary = true;
    for (auto value = boundary_down.elem_begin(); value != boundary_down.elem_end(); ++value)
      if (*value != 0.0) has_boundary = true;
    if (beam > 0.0 and not has_source and not has_boundary) {
      // Match disort::main_data's beam-only normalization exactly.
      i0_          = 1.0;
      i0_original_ = beam;
    }

    for (Index layer = 0; layer < nlayers_; ++layer) {
      scale_tau_[layer]                = 1.0 - omega_[layer] * f_[layer];
      scaled_tau_with_zero_[layer + 1] = tau_[layer] * scale_tau_[layer];
      scaled_omega_[layer]             = omega_[layer] * (1.0 - f_[layer]) / scale_tau_[layer];
      for (Index ell = 0; ell < legendre_.ncols(); ++ell)
        weighted_legendre_[layer, ell] = static_cast<Numeric>(2 * ell + 1) * legendre_[layer, ell];
      for (Index ell = 0; ell < nleg_; ++ell)
        weighted_scaled_legendre_[layer, ell] =
            static_cast<Numeric>(2 * ell + 1) * (legendre_[layer, ell] - f_[layer]) / (1.0 - f_[layer]);
    }

    Numeric omega_tau_sum   = 0.0;
    Numeric tau_sum         = 0.0;
    Numeric f_omega_tau_sum = 0.0;
    for (Index layer = 0; layer < nlayers_; ++layer) {
      omega_tau_sum   += omega_[layer] * tau_[layer];
      tau_sum         += tau_[layer];
      f_omega_tau_sum += f_[layer] * omega_[layer] * tau_[layer];
    }
    omega_average_ = omega_tau_sum / tau_sum;
    if (std::isnormal(f_omega_tau_sum)) {
      f_average_ = f_omega_tau_sum / omega_tau_sum;
      for (Index ell = 0; ell < legendre_.ncols(); ++ell) {
        Numeric numerator = 0.0;
        for (Index layer = 0; layer < nlayers_; ++layer) {
          const Numeric residue  = ell < nleg_ ? f_[layer] : legendre_[layer, ell];
          numerator             += residue * omega_[layer] * tau_[layer];
        }
        const Numeric x                = numerator / f_omega_tau_sum;
        legendre_residue_average_[ell] = static_cast<Numeric>(2 * ell + 1) * (2.0 * x - x * x);
      }
      scaled_mu0_ = mu0_ / (1.0 - omega_average_ * f_average_);
    } else {
      f_average_  = 0.0;
      scaled_mu0_ = mu0_;
    }

    Vector mu(nquad_), weights(n_);
    Legendre::PositiveDoubleGaussLegendre(mu[Range{0, n_}], weights);
    for (Index i = 0; i < n_; ++i) mu[n_ + i] = -mu[i];

    Tensor7 phase(2, nfourier_, nlayers_, nquad_, nquad_, 4, 4, 0.0);
    Tensor6 beam_phase(2, nfourier_, nlayers_, nquad_, 4, 4, 0.0);
    for (Index mode = 0; mode < std::min(nfourier_, nleg_); ++mode) {
      for (Index layer = 0; layer < nlayers_; ++layer) {
        for (Index out = 0; out < nquad_; ++out) {
          Numeric beam_value = 0.0;
          for (Index ell = mode; ell < nleg_; ++ell) {
            const Numeric normalization =
                Legendre::tgamma_ratio(static_cast<Numeric>(ell - mode + 1), static_cast<Numeric>(ell + mode + 1));
            const Numeric coefficient = weighted_scaled_legendre_[layer, ell] * normalization;
            beam_value +=
                coefficient * Legendre::assoc_legendre(ell, mode, mu[out]) * Legendre::assoc_legendre(ell, mode, -mu0_);
          }
          beam_phase[vdisort::cosine_mode, mode, layer, out, 0, 0] = beam_value;

          for (Index in = 0; in < nquad_; ++in) {
            Numeric value = 0.0;
            for (Index ell = mode; ell < nleg_; ++ell) {
              const Numeric normalization =
                  Legendre::tgamma_ratio(static_cast<Numeric>(ell - mode + 1), static_cast<Numeric>(ell + mode + 1));
              value += weighted_scaled_legendre_[layer, ell] * normalization *
                       Legendre::assoc_legendre(ell, mode, mu[out]) * Legendre::assoc_legendre(ell, mode, mu[in]);
            }
            phase[vdisort::cosine_mode, mode, layer, out, in, 0, 0] = value;
          }
        }
      }
    }

    Tensor4 vector_up(2, nfourier_, n_, 4, 0.0), vector_down(2, nfourier_, n_, 4, 0.0);
    for (Index mode = 0; mode < nfourier_; ++mode)
      for (Index i = 0; i < n_; ++i) {
        vector_up[vdisort::cosine_mode, mode, i, 0]   = boundary_up[mode, i];
        vector_down[vdisort::cosine_mode, mode, i, 0] = boundary_down[mode, i];
      }

    Tensor3 vector_source(nlayers_, source_poly.ncols(), 4, 0.0);
    Vector  source_coordinate_scale(nlayers_, 1.0);
    Vector  source_coordinate_offset(nlayers_, 0.0);
    for (Index layer = 0; layer < nlayers_; ++layer) {
      source_coordinate_scale[layer]  = 1.0 / scale_tau_[layer];
      source_coordinate_offset[layer] = tau_[layer] - scaled_tau_with_zero_[layer + 1] / scale_tau_[layer];
      for (Index power = 0; power < source_poly.ncols(); ++power)
        vector_source[layer, power, 0] = (1.0 - omega_[layer]) * source_poly[layer, power];
    }

    std::vector<vdisort::BDRF> vector_brdf;
    vector_brdf.reserve(brdf.size());
    for (Index mode = 0; mode < static_cast<Index>(brdf.size()); ++mode) {
      const BDRF scalar_mode = brdf[mode];
      const auto cosine = [scalar_mode, mode](
                              rtepack::muelmat_matrix_view out,
                              const ConstVectorView&       mu_out,
                              const ConstVectorView&       positive_mu_in) {
        Vector negative_mu_in(positive_mu_in.size());
        for (Index i = 0; i < static_cast<Index>(positive_mu_in.size()); ++i) negative_mu_in[i] = -positive_mu_in[i];
        Matrix scalar(out.nrows(), out.ncols(), 0.0);
        scalar_mode(scalar, mu_out, negative_mu_in);
        out                       = rtepack::muelmat{0.0};
        const bool    direct_beam = positive_mu_in.data_handle() != mu_out.data_handle();
        const Numeric scale       = direct_beam ? positive_mu_in[0] * Constant::inv_pi
                                                : static_cast<Numeric>(mode == 0 ? 2 : 1) * Constant::inv_pi;
        for (Index i = 0; i < scalar.nrows(); ++i)
          for (Index j = 0; j < scalar.ncols(); ++j) out[i, j][0, 0] = scale * scalar[i, j];
      };
      const auto sine = [](rtepack::muelmat_matrix_view out, const ConstVectorView&, const ConstVectorView&) {
        out = rtepack::muelmat{0.0};
      };
      vector_brdf.push_back(vdisort::BDRF{.cosine = {cosine}, .sine = {sine}});
    }

    Vector scaled_tau_values(nlayers_);
    for (Index layer = 0; layer < nlayers_; ++layer) scaled_tau_values[layer] = scaled_tau_with_zero_[layer + 1];
    Vector beam_stokes{std::max(0.0, i0_), 0.0, 0.0, 0.0};
    model_ = std::make_unique<vdisort::main_data>(nquad_,
                                                  nfourier_,
                                                  AscendingGrid{std::move(scaled_tau_values)},
                                                  scaled_omega_,
                                                  std::move(phase),
                                                  std::move(vector_up),
                                                  std::move(vector_down),
                                                  std::move(vector_source),
                                                  std::move(vector_brdf),
                                                  mu0_,
                                                  std::move(beam_stokes),
                                                  phi0_,
                                                  std::move(beam_phase),
                                                  std::move(source_coordinate_scale),
                                                  std::move(source_coordinate_offset));
  }

  main_data(main_data&&) noexcept            = default;
  main_data& operator=(main_data&&) noexcept = default;
  main_data(const main_data&)                = delete;
  main_data& operator=(const main_data&)     = delete;

  [[nodiscard]] const Vector& mu() const { return model_->mu(); }

  void u(u_data& data, const Numeric tau, const Numeric phi) const { scalar_u(data.intensities, tau, phi); }

  void u0(u0_data& data, const Numeric tau) const {
    vdisort::u0_data field;
    model_->u0(field, scaled_tau(tau));
    data.u0.resize(nquad_);
    for (Index stream = 0; stream < nquad_; ++stream) data.u0[stream] = i0_original_ * field.u0[stream, 0];
  }

  void u_corr(u_data& data, Vector& ims, tms_data&, const Numeric tau, const Numeric phi) const {
    scalar_u(data.intensities, tau, phi);
    Vector tms;
    truncation_correction(tms, tau, phi);
    improved_multiple_scattering_correction(ims, tau, phi);
    for (Index stream = 0; stream < n_; ++stream) data.intensities[stream] += i0_original_ * tms[stream];
    for (Index stream = n_; stream < nquad_; ++stream)
      data.intensities[stream] += i0_original_ * (tms[stream] + ims[stream - n_]);
  }

  [[nodiscard]] Numeric flux_up(flux_data&, const Numeric tau) const {
    vdisort::flux_data scratch;
    return i0_original_ * model_->flux_up(scratch, scaled_tau(tau));
  }

  [[nodiscard]] std::pair<Numeric, Numeric> flux_down(flux_data&, const Numeric tau) const {
    vdisort::flux_data scratch;
    const auto [diffuse, direct_scaled] = model_->flux_down(scratch, scaled_tau(tau));
    const Numeric direct_original       = i0_ > 0.0 ? i0_ * mu0_ * std::exp(-tau / mu0_) : 0.0;
    return {i0_original_ * (diffuse - direct_original + direct_scaled), i0_original_ * i0_ * direct_original};
  }

  void ungridded_u(Tensor3View out, const AscendingGrid& tau, const Vector& phi) const {
    ARTS_USER_ERROR_IF(
        (out.shape() != std::array<Index, 3>{static_cast<Index>(tau.size()), static_cast<Index>(phi.size()), nquad_}),
        "Scalar-port output has shape {:B,}, expected [{}, {}, {}]",
        out.shape(),
        tau.size(),
        phi.size(),
        nquad_);
    u_data data;
    for (Index t = 0; t < static_cast<Index>(tau.size()); ++t)
      for (Index p = 0; p < static_cast<Index>(phi.size()); ++p) {
        u(data, tau[t], phi[p]);
        out[t, p] = data.intensities;
      }
  }
};
}  // namespace vdisort_scalar_test

// VDISORT SCALAR TEST PORT END

#include "vdisort.h"

#include <arts_constants.h>
#include <debug.h>
#include <legendre.h>
#include <rtepack_multitype.h>
#include <time_report.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <numeric>
#include <ranges>

namespace {
[[nodiscard]] bool polarized_source(const Vector& stokes) { return stokes.size() == 4 and stokes[0] > 0.0; }

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
void BDRF::operator()(const Index            alpha,
                      rtepack::muelmat_matrix_view out,
                      const ConstVectorView& mu_out,
                      const ConstVectorView& mu_in) const {
  ARTS_USER_ERROR_IF(alpha != cosine_mode and alpha != sine_mode,
                     "The VDISORT BRDF mode must be 0 (cosine) or 1 (sine), got {}",
                     alpha);
  (alpha == cosine_mode ? cosine : sine)(out, mu_out, mu_in);
}

phase_matrix_data combine_phase_matrices(const rtepack::muelmat_tensor4& cosine,
                                         const rtepack::muelmat_tensor4& sine) {
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
                     "The last two phase-matrix dimensions must both be 4, got {:B,}", cosine.shape());
  rtepack::muelmat_tensor4 c(nfourier, nlayers, nout, nin), s(nfourier, nlayers, nout, nin);
  for (Index m = 0; m < nfourier; ++m)
    for (Index l = 0; l < nlayers; ++l)
      for (Index i = 0; i < nout; ++i)
        for (Index j = 0; j < nin; ++j) {
          c[m, l, i, j] = to_muelmat(cosine[m, l, i, j]);
          s[m, l, i, j] = to_muelmat(sine[m, l, i, j]);
        }
  const auto combined = combine_phase_matrices(c, s);
  Tensor7 out(2, nfourier, nlayers, nout, nin, stokes_dimension, stokes_dimension, 0.0);
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
                     "The last two beam phase-matrix dimensions must both be 4, got {:B,}", cosine.shape());
  rtepack::muelmat_tensor3 c(nfourier, nlayers, nout), s(nfourier, nlayers, nout);
  for (Index m = 0; m < nfourier; ++m)
    for (Index l = 0; l < nlayers; ++l)
      for (Index i = 0; i < nout; ++i) {
        c[m, l, i] = to_muelmat(cosine[m, l, i]);
        s[m, l, i] = to_muelmat(sine[m, l, i]);
      }
  const auto combined = combine_beam_phase_matrices(c, s);
  Tensor6 out(2, nfourier, nlayers, nout, stokes_dimension, stokes_dimension, 0.0);
  for (Index a = 0; a < 2; ++a)
    for (Index m = 0; m < nfourier; ++m)
      for (Index l = 0; l < nlayers; ++l)
        for (Index i = 0; i < nout; ++i)
          for (Index so = 0; so < stokes_dimension; ++so)
            for (Index si = 0; si < stokes_dimension; ++si)
              out[a, m, l, i, so, si] = combined[a, m, l, i][so, si];
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
      source_poly_coeffs(NLayers, Nscoeffs, stokes_dimension),
      source_coordinate_scale(NLayers, 1.0),
      source_coordinate_offset(NLayers, 0.0),
      phase_matrix(2, NFourier, NLayers, NQuad, NQuad, rtepack::muelmat{0.0}),
      boundary_up(2, NFourier, N, stokes_dimension),
      boundary_down(2, NFourier, N, stokes_dimension),
      brdf_fourier_modes(NBDRF),
      beam_stokes(stokes_dimension, 0.0),
      beam_phase_matrix(2, NFourier, NLayers, NQuad, rtepack::muelmat{0.0}),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      G_collect(2, NFourier, NLayers, NState, NState),
      K_collect(2, NFourier, NLayers, NState),
      GC_collect(2, NFourier, NLayers, NState),
      B_collect(2, NFourier, NLayers, NState),
      source_collect(2, NFourier, NLayers, NState, Nscoeffs),
      um(NLayers, 2, NFourier, NState),
      top_anchored(static_cast<std::size_t>(2 * NFourier * NLayers * NState), 0) {
  ARTS_USER_ERROR_IF(NQuad <= 0 or NQuad % 2 != 0, "NQuad must be a positive even number, got {}", NQuad);
  Legendre::PositiveDoubleGaussLegendre(mu_arr[Range{0, N}], W);
  std::transform(mu_arr.begin(), mu_arr.begin() + N, mu_arr.begin() + N, [](const Numeric x) { return -x; });
  std::transform(mu_arr.begin(), mu_arr.end(), inv_mu_arr.begin(), [](const Numeric x) { return 1.0 / x; });
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
                std::move(boundary_up_),
                std::move(boundary_down_),
                std::move(source_poly_coeffs_),
                std::move(brdf_fourier_modes_),
                mu0_,
                std::move(beam_stokes_),
                phi0_,
                to_beam_phase_matrix_data(beam_phase_matrix_),
                std::move(source_coordinate_scale_),
                std::move(source_coordinate_offset_)) {}

main_data::main_data(Index                  NQuad_,
                     Index                  NFourier_,
                     AscendingGrid          tau_arr_,
                     Vector                 omega_arr_,
                     phase_matrix_data      phase_matrix_,
                     Tensor4                boundary_up_,
                     Tensor4                boundary_down_,
                     Tensor3                source_poly_coeffs_,
                     std::vector<BDRF>      brdf_fourier_modes_,
                     Numeric                mu0_,
                     Vector                 beam_stokes_,
                     Numeric                phi0_,
                     beam_phase_matrix_data beam_phase_matrix_,
                     Vector                 source_coordinate_scale_,
                     Vector                 source_coordinate_offset_)
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
  for (Index p = Nscoeffs - 1; p >= 0; --p)
    value = std::fma(value, source_tau, source_collect[alpha, m, layer, state, p]);
  if (has_beam_source) value += B_collect[alpha, m, layer, state] * std::exp(-tau / mu0);
  return value;
}

void main_data::check_input_size() const {
  ARTS_USER_ERROR_IF(
      static_cast<Index>(tau_arr.size()) != NLayers, "tau_arr has size {}, expected {}", tau_arr.size(), NLayers);
  ARTS_USER_ERROR_IF(
      static_cast<Index>(omega_arr.size()) != NLayers, "omega_arr has size {}, expected {}", omega_arr.size(), NLayers);
  ARTS_USER_ERROR_IF((source_poly_coeffs.shape() != std::array<Index, 3>{NLayers, Nscoeffs, stokes_dimension}),
                     "source_poly_coeffs has shape {:B,}, expected [{}, {}, 4]",
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
  ARTS_USER_ERROR_IF((boundary_up.shape() != std::array<Index, 4>{2, NFourier, N, stokes_dimension}),
                     "boundary_up has shape {:B,}, expected [2, {}, {}, 4]",
                     boundary_up.shape(),
                     NFourier,
                     N);
  ARTS_USER_ERROR_IF((boundary_down.shape() != std::array<Index, 4>{2, NFourier, N, stokes_dimension}),
                     "boundary_down has shape {:B,}, expected [2, {}, {}, 4]",
                     boundary_down.shape(),
                     NFourier,
                     N);
  ARTS_USER_ERROR_IF(beam_stokes.size() != stokes_dimension, "beam_stokes has size {}, expected 4", beam_stokes.size());
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
  ARTS_USER_ERROR_IF(std::ranges::any_of(omega_arr, [](const Numeric x) { return x < 0.0 or x >= 1.0; }),
                     "omega_arr must be in [0, 1), got {:B,}",
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
  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      for (Index l = 0; l < NLayers; ++l) {
        Matrix A_real(NState, NState, 0.0);
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

        ComplexMatrix A(NState, NState);
        for (Index i = 0; i < NState; ++i)
          for (Index j = 0; j < NState; ++j) A[i, j] = A_real[i, j];

        ComplexMatrix eigenvectors(NState, NState);
        ComplexVector eigenvalues(NState);
        ::diagonalize(eigenvectors, eigenvalues, A);

        std::vector<Index> order(static_cast<std::size_t>(NState));
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
          const rtepack::stokvec beam = rtepack::to_stokvec(beam_stokes);
          Vector        rhs(NState, 0.0);
          const Numeric epsilon = m == 0 ? 1.0 : 2.0;
          for (Index i = 0; i < NQuad; ++i) {
            const rtepack::stokvec scattered = beam_phase_matrix[alpha, m, l, i] * beam;
            for (Index so = 0; so < stokes_dimension; ++so) {
              rhs[state_index(i, so)] =
                  inv_mu_arr[i] * epsilon * omega_arr[l] * scattered[so] / (4.0 * Constant::pi);
            }
          }
          diagonal(A_real) += 1.0 / mu0;
          solve_inplace(rhs, A_real);
          for (Index s = 0; s < NState; ++s) B_collect[alpha, m, l, s] = rhs[s];
        }
      }
    }
  }
  // VDISORT CHANGE END
}

void main_data::source_function() {
  ARTS_TIME_REPORT
  source_collect = 0.0;
  if (not has_source_poly) return;

  // VDISORT CHANGE BEGIN: polynomial source coefficients are four-vectors.
  // The m=0 combined systems split [I,Q] and [U,V] as in paper Eq. (82).
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
            const bool    active   = alpha == cosine_mode ? s < 2 : s >= 2;
            const Numeric q        = active ? source_poly_coeffs[l, p, s] : 0.0;
            rhs[state_index(i, s)] = inv_mu_arr[i] * q + static_cast<Numeric>(p + 1) * next[state_index(i, s)];
          }
        }
        Matrix work{A};
        solve_inplace(rhs, work);
        for (Index s = 0; s < NState; ++s) {
          source_collect[alpha, m, l, s, p] = rhs[s];
          next[s]                           = rhs[s];
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

  const Index equation_count = NLayers * NState;
  const Index bandwidth      = 3 * NHalfState - 1;

  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index m = 0; m < NFourier; ++m) {
      matpack::complex_band_matrix lhs(bandwidth, bandwidth, equation_count, equation_count);
      ComplexVector       rhs(equation_count, 0.0);

      // VDISORT CHANGE BEGIN: reflection is a matrix coupling Stokes components.
      Matrix reflection(NHalfState, NHalfState, 0.0);
      Vector direct_reflection(NHalfState, 0.0);
      if (m < NBDRF) {
        rtepack::muelmat_matrix raw(N, N, rtepack::muelmat{0.0});
        brdf_fourier_modes[m](alpha, raw, mu_arr[Range{0, N}], mu_arr[Range{0, N}]);
        for (Index i = 0; i < N; ++i)
          for (Index so = 0; so < stokes_dimension; ++so)
            for (Index j = 0; j < N; ++j)
              for (Index si = 0; si < stokes_dimension; ++si)
                reflection[state_index(i, so), state_index(j, si)] =
                    Constant::pi * W[j] * mu_arr[j] * raw[i, j][so, si];

        if (has_beam_source) {
          rtepack::muelmat_matrix beam_raw(N, 1, rtepack::muelmat{0.0});
          const Vector beam_mu{mu0};
          brdf_fourier_modes[m](alpha, beam_raw, mu_arr[Range{0, N}], beam_mu);
          const Numeric attenuation = std::exp(-tau_arr.back() / mu0);
          const rtepack::stokvec beam = rtepack::to_stokvec(beam_stokes);
          for (Index i = 0; i < N; ++i) {
            const rtepack::stokvec reflected = attenuation * (beam_raw[i, 0] * beam);
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
          rhs[row]          = boundary_down[alpha, m, i, s] - particular(alpha, m, 0, state, 0.0);
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
              boundary_up[alpha, m, i, s] + direct_reflection[hrow] - particular(alpha, m, last, pos_state, tau);
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
        for (Index state = 0; state < NState; ++state) um[l, alpha, m, state] = field[alpha * NFourier + m, state];
  }
  // VDISORT CHANGE END
}

void main_data::update_all() {
  ARTS_TIME_REPORT
  has_source_poly = Nscoeffs > 0;
  has_beam_source = polarized_source(beam_stokes);
  check_input_size();
  check_input_value();
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
      for (Index state = 0; state < NState; ++state) {
        Complex value = particular(alpha, m, layer, state, tau);
        for (Index e = 0; e < NState; ++e)
          value += GC_collect[alpha, m, layer, e] * homogeneous(alpha, m, layer, state, e, tau);
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
  ARTS_USER_ERROR_IF(phi < 0.0 or phi >= Constant::two_pi, "phi must be in [0, 2*pi), got {}", phi);
  Matrix combined(2 * NFourier, NState);
  combined_field(combined, tau);

  data.intensities.resize(NQuad, stokes_dimension);
  data.intensities = 0.0;
  for (Index m = 0; m < NFourier; ++m) {
    const Numeric c = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const Numeric s = std::sin(static_cast<Numeric>(m) * (phi0 - phi));
    for (Index i = 0; i < NQuad; ++i) {
      for (Index stokes = 0; stokes < 2; ++stokes)
        data.intensities[i, stokes] +=
            combined[m, state_index(i, stokes)] * c + combined[NFourier + m, state_index(i, stokes)] * s;
      // VDISORT CHANGE: U and V swap combined cosine/sine roles, paper Eq. (78).
      for (Index stokes = 2; stokes < 4; ++stokes)
        data.intensities[i, stokes] +=
            combined[NFourier + m, state_index(i, stokes)] * c + combined[m, state_index(i, stokes)] * s;
    }
  }
}

void main_data::u0(u0_data& data, const Numeric tau) const {
  Matrix combined(2 * NFourier, NState);
  combined_field(combined, tau);
  data.u0.resize(NQuad, stokes_dimension);
  for (Index i = 0; i < NQuad; ++i) {
    data.u0[i, 0] = combined[0, state_index(i, 0)];
    data.u0[i, 1] = combined[0, state_index(i, 1)];
    data.u0[i, 2] = combined[NFourier, state_index(i, 2)];
    data.u0[i, 3] = combined[NFourier, state_index(i, 3)];
  }
}

Numeric main_data::flux_up(flux_data& data, const Numeric tau) const {
  u0_data field;
  u0(field, tau);
  data.u0        = std::move(field.u0);
  Numeric result = 0.0;
  for (Index i = 0; i < N; ++i) result += W[i] * mu_arr[i] * data.u0[i, 0];
  return Constant::two_pi * result;
}

std::pair<Numeric, Numeric> main_data::flux_down(flux_data& data, const Numeric tau) const {
  u0_data field;
  u0(field, tau);
  data.u0         = std::move(field.u0);
  Numeric diffuse = 0.0;
  for (Index i = 0; i < N; ++i) diffuse += W[i] * mu_arr[i] * data.u0[N + i, 0];
  const Numeric direct = has_beam_source ? mu0 * beam_stokes[0] * std::exp(-tau / mu0) : 0.0;
  return {Constant::two_pi * diffuse, direct};
}

void main_data::gridded_u(Tensor4View out, const Vector& phi) const {
  ARTS_USER_ERROR_IF(
      (out.shape() != std::array<Index, 4>{NLayers, static_cast<Index>(phi.size()), NQuad, stokes_dimension}),
      "gridded_u output has shape {:B,}, expected [{}, {}, {}, 4]",
      out.shape(),
      NLayers,
      phi.size(),
      NQuad);
  u_data data;
  for (Index l = 0; l < NLayers; ++l)
    for (Index p = 0; p < static_cast<Index>(phi.size()); ++p) {
      u(data, tau_arr[l], phi[p]);
      out[l, p] = data.intensities;
    }
}

void main_data::gridded_flux(VectorView up, VectorView down, VectorView down_direct) const {
  ARTS_USER_ERROR_IF(up.size() != static_cast<Size>(NLayers) or down.size() != static_cast<Size>(NLayers) or
                         down_direct.size() != static_cast<Size>(NLayers),
                     "All gridded flux outputs must have size {}",
                     NLayers);
  flux_data data;
  for (Index l = 0; l < NLayers; ++l) {
    up[l]                             = flux_up(data, tau_arr[l]);
    std::tie(down[l], down_direct[l]) = flux_down(data, tau_arr[l]);
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
  u_data data;
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t)
    for (Index p = 0; p < static_cast<Index>(phi.size()); ++p) {
      u(data, tau[t], phi[p]);
      out[t, p] = data.intensities;
    }
}

void main_data::ungridded_flux(VectorView up, VectorView down, VectorView down_direct, const AscendingGrid& tau) const {
  ARTS_USER_ERROR_IF(up.size() != tau.size() or down.size() != tau.size() or down_direct.size() != tau.size(),
                     "All ungridded flux outputs must have the same size as tau ({})",
                     tau.size());
  flux_data data;
  for (Index t = 0; t < static_cast<Index>(tau.size()); ++t) {
    up[t]                             = flux_up(data, tau[t]);
    std::tie(down[t], down_direct[t]) = flux_down(data, tau[t]);
  }
}

ConstTensor3View main_data::layer_um(const Size layer) const {
  ARTS_USER_ERROR_IF(layer >= static_cast<Size>(NLayers), "Layer {} is out of range [0, {})", layer, NLayers);
  return um[layer];
}

void main_data::set_beam_source(Vector beam) {
  ARTS_USER_ERROR_IF(beam.size() != stokes_dimension, "A VDISORT beam source must have four Stokes components");
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

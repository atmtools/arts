#include "disort.h"

#include <arts_constants.h>
#include <arts_conversions.h>
#include <compare.h>
#include <debug.h>
#include <legendre.h>
#include <matpack.h>
#include <time_report.h>
#include <xml.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

namespace disort {
void BDRF::operator()(MatrixView x,
                      const ConstVectorView& a,
                      const ConstVectorView& b) const {
  ARTS_TIME_REPORT

  assert(static_cast<Size>(x.nrows()) == a.size());
  assert(static_cast<Size>(x.ncols()) == b.size());
  f(x, a, b);
}

Matrix BDRF::operator()(const Vector& a, const Vector& b) const {
  ARTS_TIME_REPORT

  Matrix x(a.size(), b.size());
  f(x, a, b);
  return x;
}

namespace {
constexpr Range rf(Size N) { return {0, N}; }
constexpr Range rb(Size N) { return {N, N}; }

void source_set_k1(mathscr_v_data& data,
                   const ConstMatrixView& G,
                   const ConstVectorView& inv_mu) {
  stdr::copy(inv_mu, data.k1.elem_begin());
  stdr::copy(elemwise_range(G), data.G.elem_begin());
  solve_inplace(data.k1, data.G, data.solve_work);
}

void source_set_k2(mathscr_v_data& data,
                   const Numeric tau,
                   const Numeric omega,
                   const ConstVectorView& source_poly_coeffs,
                   const ConstVectorView& K) {
  const Index Nk = K.size();
  const Index Nc = source_poly_coeffs.size();
  const Index n  = Nc - 1;

  Numeric x = 1 - omega;
  for (Index i = n; i >= 0; x *= tau, i--) data.cvec[i] = x;

  const auto leg = [n, &cvec = data.cvec, &source_poly_coeffs](Index i,
                                                               Index j) {
    return cvec[i] * Legendre::factorial(n - j) * source_poly_coeffs[n - j];
  };

  for (Index k = 0; k < Nk; k++) {
    data.k2[k] = 0.0;
    for (Index i = 0; i < Nc; i++) {
      Numeric fac = 1.0 / (std::pow(K[k], i + 1) * Legendre::factorial(n - i));

      for (Index j = 0; j < i; fac *= K[k], j++) data.k2[k] += leg(i, j) * fac;
      data.k2[k] += leg(i, i) * fac;
    }
  }
}

void source_scale_k2(mathscr_v_data& data) { data.k2 *= data.k1; }

void source_update_um(mathscr_v_data& data,
                      VectorView um,
                      const ConstMatrixView& G,
                      const Index Ni0   = 0,
                      const Numeric scl = 1.0,
                      const Numeric add = 0.0) {
  const Index Ni = um.size();
  mult(um, G[Range{Ni0, Ni}], data.k2, scl, add);
}

void mathscr_v(VectorView um,
               mathscr_v_data& data,
               const Numeric tau,
               const Numeric omega,
               const ConstVectorView& source_poly_coeffs,
               const ConstMatrixView& G,
               const ConstVectorView& K,
               const ConstVectorView& inv_mu,
               const Index Ni0   = 0,
               const Numeric scl = 1.0,
               const Numeric add = 0.0) {
  source_set_k1(data, G, inv_mu);
  source_set_k2(data, tau, omega, source_poly_coeffs, K);
  source_scale_k2(data);
  source_update_um(data, um, G, Ni0, scl, add);
}

void source_terms(MatrixView SRC0,
                  MatrixView SRC1,
                  VectorView SRCB,
                  mathscr_v_data& data,
                  const Vector& tau,
                  const Vector& omega,
                  const Matrix& Sc,
                  const ConstTensor3View& Gm,
                  const ConstMatrixView& Km,
                  const Vector& inv_mu_arr,
                  const Index N,
                  const Index NLayers) {
  source_set_k1(data, Gm[0], inv_mu_arr);
  source_set_k2(data, 0.0, omega[0], Sc[0], Km[0]);
  source_scale_k2(data);
  source_update_um(data, SRCB[rf(N)], Gm[0], N);

  for (Index l = 0; l < NLayers; l++) {
    source_set_k2(data, tau[l], omega[l], Sc[l], Km[l]);
    source_scale_k2(data);
    source_update_um(data, SRC0[l], Gm[l]);

    if (l < NLayers - 1) {
      source_set_k1(data, Gm[l + 1], inv_mu_arr);

      source_set_k2(data, tau[l], omega[l + 1], Sc[l + 1], Km[l + 1]);
      source_scale_k2(data);
      source_update_um(data, SRC1[l], Gm[l + 1]);
    }
  }

  const Index ln = NLayers - 1;
  source_set_k2(data, tau.back(), omega[ln], Sc[ln], Km[ln]);
  source_scale_k2(data);
  source_update_um(data, SRCB[rb(N)], Gm[ln]);
}
}  // namespace

void main_data::solve_for_coefs() {
  ARTS_TIME_REPORT

  const Index ln = NLayers - 1;

  //! FIXME: Original code is transposed, but I suspect it is a bug
  auto RHS_middle = RHS[Range{N, n - NQuad}].view_as(NLayers - 1, NQuad);

  for (Index m = 0; m < NFourier; m++) {
    const bool m_equals_0_bool = m == 0;
    const bool BDRF_bool       = m < NBDRF;
    const auto G_collect_m     = G_collect[m];
    const auto B_collect_m     = B_collect[m];

    if (BDRF_bool) {
      brdf_fourier_modes[m](mathscr_D_neg, mu_arr[rf(N)], mu_arr[rb(N)]),
          einsum<"ij", "", "ij", "j", "j">(
              R, 1 + m_equals_0_bool, mathscr_D_neg, mu_arr[rf(N)], W);
      if (has_beam_source) {
        brdf_fourier_modes[m](
            mathscr_X_pos.view_as(N, 1), mu_arr[rf(N)], ConstVectorView{-mu0});
        mathscr_X_pos *= mu0 * I0 / Constant::pi;
      }
    }

    const auto b_pos_m = b_pos[m];
    const auto b_neg_m = b_neg[m];

    // Fill RHS
    {
      ARTS_NAMED_TIME_REPORT("disort::rhs"s);

      if (has_source_poly and m_equals_0_bool) {
        for (Index i = 0; i < N; i++) RHS[i] = -SRCB[i];
        for (Index i = 0; i < N; i++) RHS[n - N + i] = -SRCB[i + N];
        std::transform(SRC1.elem_begin(),
                       SRC1.elem_end() - NQuad,
                       SRC0.elem_begin(),
                       RHS_middle.elem_begin(),
                       std::minus{});

        if (NBDRF > 0) {
          source_update_um(comp_data, jvec[rf(N)], G_collect_m[ln], N);
          mult(RHS[Range{n - N, N}], R, jvec[rf(N)], 1.0, 1.0);
        }
      } else {
        RHS = 0.0;
      }

      if (has_beam_source) {
        if (BDRF_bool) {
          stdr::copy(mathscr_X_pos, BDRF_RHS_contribution.begin());
          mult(BDRF_RHS_contribution, R, B_collect_m[ln, rf(N)], 1.0, 1.0);
        } else {
          BDRF_RHS_contribution = 0.0;
        }

        for (Index l = 0; l < ln; l++) {
          const Numeric scl = std::exp(-mu0 * scaled_tau_arr_with_0[l + 1]);
          for (Index j = 0; j < NQuad; j++) {
            RHS_middle[l, j] +=
                (B_collect_m[l + 1, j] - B_collect_m[l, j]) * scl;
          }
        }

        for (Index i = 0; i < N; i++) {
          RHS[i] += b_neg_m[i] - B_collect_m[0, N + i];
          RHS[n - N + i] +=
              b_pos_m[i] + (BDRF_RHS_contribution[i] - B_collect_m[ln, i]) *
                               std::exp(-scaled_tau_arr_with_0.back() / mu0);
        }
      } else {
        RHS[rf(N)]           += b_neg_m;
        RHS[Range{n - N, N}] += b_pos_m;
      }
    }

    // Fill LHS
    {
      ARTS_NAMED_TIME_REPORT("disort::lhs"s);

      if (BDRF_bool) {
        mult(BDRF_LHS, R, G_collect_m[ln, rb(N)]);
      } else if (m == NBDRF) {  // only once
        BDRF_LHS = 0;
      }

      for (Index j = 0; j < N; j++) {
        for (Index i = 0; i < N; i++) {
          LHSB[i, j]     = G_collect_m[0, i + N, j];
          LHSB[i, N + j] = G_collect_m[0, i + N, j + N] * expK_collect[m, 0, j];
          LHSB[n - N + i, n - 2 * N + j] =
              (G_collect_m[ln, i, j] - BDRF_LHS[i, j]) * expK_collect[m, ln, j];
          LHSB[n - N + i, n - N + j] =
              G_collect_m[ln, i, j + N] - BDRF_LHS[i, j + N];
        }
      }

      for (Index l = 0; l < ln; l++) {
        for (Index j = 0; j < N; j++) {
          const Numeric e1 = 1.0 / expK_collect[m, l, j + N];
          const Numeric e2 = 1.0 / expK_collect[m, l + 1, j + N];
          for (Index i = 0; i < N; i++) {
            LHSB[N + l * NQuad + i, l * NQuad + j] = G_collect_m[l, i, j] * e1;
            LHSB[2 * N + l * NQuad + i, l * NQuad + j] =
                G_collect_m[l, N + i, j] * e1;
            LHSB[N + l * NQuad + i, l * NQuad + 2 * NQuad - N + j] =
                -G_collect_m[l + 1, i, N + j] * e2;
            LHSB[2 * N + l * NQuad + i, l * NQuad + 2 * NQuad - N + j] =
                -G_collect_m[l + 1, N + i, N + j] * e2;
          }
        }

        for (Index i = 0; i < NQuad; i++) {
          for (Index j = 0; j < N; j++) {
            LHSB[N + l * NQuad + i, l * NQuad + N + j] =
                G_collect_m[l, i, N + j];
            LHSB[N + l * NQuad + i, l * NQuad + 2 * N + j] =
                -G_collect_m[l + 1, i, j];
          }
        }
      }
    }

    {
      ARTS_NAMED_TIME_REPORT("disort::solve-band"s);

      LHSB.solve(RHS);

      einsum<"ijm", "ijm", "im">(
          GC_collect[m], G_collect_m, RHS.view_as(NLayers, NQuad));
    }
  }
}

namespace {
Numeric poch(Index x, Index n) {
  return Legendre::tgamma_ratio(static_cast<Numeric>(x + n),
                                static_cast<Numeric>(x));
}
}  // namespace

void main_data::diagonalize() {
  ARTS_TIME_REPORT

  for (Index m = 0; m < NFourier; m++) {
    auto Km = K_collect[m];
    auto Gm = G_collect[m];
    auto Bm = B_collect[m];

    D_temp.resize(N, NLeg - m);
    X_temp.resize(NLeg - m);
    asso_leg_term_pos.resize(NLeg - m, N);
    asso_leg_term_neg.resize(NLeg - m, N);
    asso_leg_term_mu0.resize(NLeg - m);
    weighted_asso_Leg_coeffs_l.resize(NLeg - m);

    auto xtemp = X_temp.view_as(1, NLeg - m);

    const bool m_equals_0_bool = (m == 0);

    fac.resize(NLeg - m);
    for (Index i = m; i < NLeg; i++) fac[i - m] = poch(i + m + 1, -2 * m);

    for (Index i = m; i < NLeg; i++) {
      for (Index j = 0; j < N; j++) {
        asso_leg_term_pos[i - m, j] = Legendre::assoc_legendre(i, m, mu_arr[j]);
        asso_leg_term_neg[i - m, j] =
            asso_leg_term_pos[i - m, j] * ((i - m) % 2 ? -1.0 : 1.0);
      }
      asso_leg_term_mu0[i - m] = Legendre::assoc_legendre(i, m, -mu0);
    }

    const bool all_asso_leg_term_pos_finite =
        stdr::all_of(elemwise_range(asso_leg_term_pos),
                     [](auto& x) { return std::isfinite(x); });

    for (Index l = 0; l < NLayers; l++) {
      VectorView K = Km[l];
      MatrixView G = Gm[l];

      for (Index i = 0; i < NLeg - m; i++) {
        weighted_asso_Leg_coeffs_l[i] =
            weighted_scaled_Leg_coeffs[l, i + m] * fac[i];
      }

      const Numeric scaled_omega_l = scaled_omega_arr[l];

      if (scaled_omega_l != 0.0 or
          (all_asso_leg_term_pos_finite and
           stdr::any_of(elemwise_range(weighted_asso_Leg_coeffs_l),
                        Cmp::gt(0)))) {
        einsum<"ij", "j", "ji">(
            D_temp, weighted_asso_Leg_coeffs_l, asso_leg_term_pos);
        mult(D_pos, D_temp, asso_leg_term_pos, 0.5 * scaled_omega_l);
        mult(D_neg, D_temp, asso_leg_term_neg, 0.5 * scaled_omega_l);

        einsum<"ij", "i", "ij", "j">(sqr, inv_mu_arr[rf(N)], D_neg, W);
        einsum<"ij", "i", "ij", "j">(apb, inv_mu_arr[rf(N)], D_pos, W);
        diagonal(apb) -= inv_mu_arr[rf(N)];

        amb  = apb;  // still just alpha
        apb += sqr;  // sqr is beta
        amb -= sqr;

        VectorView eval = K[rf(N)];
        MatrixView evec = amb;
        MatrixView AB   = apb;

        mult(sqr, evec, AB);

        diagonalize_inplace(evec, eval, sqr, diag_work);

        for (Index i = 0; i < N; i++) {
          const Numeric sqrt_x = std::sqrt(std::abs(eval[i]));
          K[i]                 = -sqrt_x;
          K[i + N]             = sqrt_x;
        }

        mult(sqr, AB, evec);

        for (Index i = 0; i < N; i++) {
          for (Index j = 0; j < N; j++) {
            const Numeric a = evec[i, j];
            const Numeric b = sqr[i, j] / K[j];
            G[i, j]         = 0.5 * (a - b);
            G[i, j + N]     = 0.5 * (a + b);
            G[i + N, j]     = G[i, j + N];
            G[i + N, j + N] = G[i, j];
          }
        }

        if (has_beam_source) {
          einsum<"i", "i", "i", "">(
              X_temp,
              weighted_asso_Leg_coeffs_l,
              asso_leg_term_mu0,
              (scaled_omega_l * I0 * (2 - m_equals_0_bool) /
               (4 * Constant::pi)));

          mult(jvec[rf(N)].view_as(1, N), xtemp, asso_leg_term_pos, -1);
          jvec[rf(N)] *= inv_mu_arr[rf(N)];

          mult(jvec[rb(N)].view_as(1, N), xtemp, asso_leg_term_neg);
          jvec[rb(N)] *= inv_mu_arr[rf(N)];

          stdr::copy(elemwise_range(G), Gml.elem_begin());
          solve_inplace(jvec, Gml, solve_work);

          for (Index j = 0; j < NQuad; j++) jvec[j] *= mu0 / (1.0 + K[j] * mu0);

          mult(Bm[l], G, jvec, -1);
        }
      } else {
        for (Index i = 0; i < N; i++) {
          G[i + N, i] = 1;
          G[i, i + N] = 1;
          K[i]        = -inv_mu_arr[i];
          K[i + N]    = inv_mu_arr[i];
        }
      }
    }
  }
}

/** Computes the IMS factors
 * 
 * Dependent on:
 * - omega_arr
 * - tau_arr
 * - f_arr
 * - mu0
 * - Leg_coeffs_all
 *
 * Modifies:
 * - scaled_mu0
 * - Leg_coeffs_residue_avg
 * - omega_avg
 * - f_avg
 */
void main_data::set_ims_factors() {
  ARTS_TIME_REPORT

  const Numeric sum1 = dot(omega_arr, tau_arr.vec());
  omega_avg          = sum1 / sum(tau_arr.vec());

  const Numeric sum2 =
      einsum<Numeric, "", "i", "i", "i">({}, f_arr, omega_arr, tau_arr);
  if (std::isnormal(sum2)) {
    f_avg = sum2 / sum1;

    for (Index i = 0; i < NLeg_all; i++) {
      Numeric sum3 = 0.0;
      if (not f_arr.empty()) {
        if (i < NLeg) {
          for (Index j = 0; j < NLayers; j++) {
            sum3 += f_arr[j] * omega_arr[j] * tau_arr[j];
          }
        } else {
          for (Index j = 0; j < NLayers; j++) {
            sum3 += Leg_coeffs_all[j, i] * omega_arr[j] * tau_arr[j];
          }
        }
      }

      const Numeric x = sum3 / sum2;
      Leg_coeffs_residue_avg[i] =
          static_cast<Numeric>(2 * i + 1) * (2.0 * x - Math::pow2(x));
    }

    scaled_mu0 = mu0 / (1.0 - omega_avg * f_avg);
  } else {
    f_avg                  = 0.0;
    Leg_coeffs_residue_avg = 0.0;
    scaled_mu0             = mu0;
  }
}

void main_data::set_scales() {
  ARTS_TIME_REPORT

  std::transform(omega_arr.begin(),
                 omega_arr.end(),
                 f_arr.begin(),
                 scaled_omega_arr.begin(),
                 [](auto&& omega, auto&& f) { return 1.0 - omega * f; });

  for (Index i = 0; i < NLayers; i++) {
    scale_tau[i] = 1.0 - omega_arr[i] * f_arr[i];
  }

  scaled_tau_arr_with_0[0] = 0;
  einsum<"i", "i", "i">(
      scaled_tau_arr_with_0[Range{1, NLayers}], tau_arr, scale_tau);

  for (Index i = 0; i < NLayers; i++) {
    for (Index j = 0; j < NLeg; j++) {
      weighted_scaled_Leg_coeffs[i, j] = static_cast<Numeric>(2 * j + 1) *
                                         (Leg_coeffs_all[i, j] - f_arr[i]) /
                                         (1 - f_arr[i]);
    }
  }

  for (Index i = 0; i < NLayers; i++) {
    scaled_omega_arr[i] = omega_arr[i] * (1.0 - f_arr[i]) / scale_tau[i];
  }
}

void main_data::set_weighted_Leg_coeffs_all() {
  ARTS_TIME_REPORT

  for (Index j = 0; j < NLayers; j++) {
    for (Index i = 0; i < NLeg_all; i++) {
      weighted_Leg_coeffs_all[joker, i] =
          static_cast<Numeric>(2 * i + 1) * Leg_coeffs_all[j, i];
    }
  }
}

void main_data::set_beam_source(const Numeric I0_) {
  ARTS_TIME_REPORT

  has_beam_source = I0_ > 0;

  if (stdr::all_of(elemwise_range(b_pos), Cmp::eq(0)) and
      not has_source_poly and has_beam_source) {
    I0_orig = I0_;
    I0      = 1;
  } else {
    I0_orig = 1;
    I0      = I0_;
  }
}

void main_data::check_input_size() const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(static_cast<Index>(tau_arr.size()) != NLayers,
                     "{} vs {}",
                     tau_arr.size(),
                     NLayers);

  ARTS_USER_ERROR_IF(static_cast<Index>(omega_arr.size()) != NLayers,
                     "{} vs {}",
                     omega_arr.size(),
                     NLayers);

  ARTS_USER_ERROR_IF(
      (source_poly_coeffs.shape() != std::array{NLayers, Nscoeffs}),
      "{:B,} vs [{}, {}]",
      source_poly_coeffs.shape(),
      NLayers,
      Nscoeffs);

  ARTS_USER_ERROR_IF(static_cast<Index>(f_arr.size()) != NLayers,
                     "{} vs {}",
                     f_arr.size(),
                     NLayers);

  ARTS_USER_ERROR_IF((Leg_coeffs_all.shape() != std::array{NLayers, NLeg_all}),
                     "{:B,} vs [{}, {}]",
                     Leg_coeffs_all.shape(),
                     NLayers,
                     NLeg_all);

  ARTS_USER_ERROR_IF((b_pos.shape() != std::array{NFourier, N}),
                     "{:B,} vs [{}, {}]",
                     b_pos.shape(),
                     NFourier,
                     N)

  ARTS_USER_ERROR_IF((b_neg.shape() != std::array{NFourier, N}),
                     "{:B,} vs [{}, {}]",
                     b_neg.shape(),
                     NFourier,
                     N)

  ARTS_USER_ERROR_IF(
      brdf_fourier_modes.size() != static_cast<std::size_t>(NBDRF),
      "{} vs {}",
      brdf_fourier_modes.size(),
      NBDRF);
}

void main_data::check_input_value() const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau_arr.front() <= 0.0,
                     "tau_arr must be strictly positive, got {:B,}",
                     tau_arr);

  ARTS_USER_ERROR_IF(
      stdr::any_of(omega_arr,
                   [](auto&& omega) { return omega >= 1 or omega < 0; }),
      "omega_arr must be [0, 1), but got {:B,}",
      omega_arr);

  ARTS_USER_ERROR_IF(
      stdr::any_of(Leg_coeffs_all,
                   [](auto&& x) {
                     return x[0] != 1 or stdr::any_of(x, [](auto&& u) {
                              return std::abs<Numeric>(u) > 1;
                            });
                   }),
      "Leg_coeffs_all must have 1 in the first column and be [-1, 1] elsewhere, got {:B,}",
      Leg_coeffs_all);

  ARTS_USER_ERROR_IF(I0 < 0, "I0 must be non-negative, got {}", I0);

  ARTS_USER_ERROR_IF(phi0 < 0 or phi0 >= Constant::two_pi,
                     "phi0 must be [0, 2*pi), got {}",
                     phi0);

  ARTS_USER_ERROR_IF(
      stdr::any_of(f_arr, [](auto&& x) { return x > 1 or x < 0; }),
      "f_arr must be [0, 1], got {:B,}",
      f_arr);

  ARTS_USER_ERROR_IF(mu0 < 0 or mu0 > 1, "mu0 must be [0, 1], got {}", mu0);

  ARTS_USER_ERROR_IF(
      stdr::any_of(mu_arr[rf(N)],
                   [mu = mu0](auto&& x) { return std::abs(x - mu) < 1e-8; }),
      "mu0 in mu_arr, this creates a singularity.  Change NQuad or mu0. Got mu_arr {:B,} for mu0 {}",
      mu_arr,
      mu0);
}

void main_data::transmission() {
  ARTS_TIME_REPORT

  for (Index m = 0; m < NFourier; m++) {
    for (Index l = 0; l < NLayers; l++) {
      for (Index i = 0; i < NQuad; i++) {
        expK_collect[m, l, i] =
            std::exp(K_collect[m, l, i] *
                     (scaled_tau_arr_with_0[l + 1] - scaled_tau_arr_with_0[l]));
      }
    }
  }
}

void main_data::source_function() {
  ARTS_TIME_REPORT

  source_terms(SRC0,
               SRC1,
               SRCB,
               comp_data,
               tau_arr,
               omega_arr,
               source_poly_coeffs,
               G_collect[0],
               K_collect[0],
               inv_mu_arr,
               N,
               NLayers);
}

void main_data::update_all(const Numeric I0_) try {
  ARTS_TIME_REPORT

  check_input_value();

  set_weighted_Leg_coeffs_all();
  if (I0_ >= 0) set_beam_source(I0_);
  set_scales();
  set_ims_factors();
  diagonalize();
  transmission();
  source_function();
  solve_for_coefs();
}
ARTS_METHOD_ERROR_CATCH

main_data::main_data(const Index NLayers_,
                     const Index NQuad_,
                     const Index NLeg_,
                     const Index NFourier_,
                     const Index Nscoeffs_,
                     const Index NLeg_all_,
                     const Index NBDRF_)
    : NLayers(NLayers_),
      NQuad(NQuad_),
      NLeg(NLeg_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(Nscoeffs_),
      NLeg_all(NLeg_all_),
      NBDRF(NBDRF_),
      has_source_poly(Nscoeffs > 0),
      // User data
      tau_arr(NLayers),
      omega_arr(NLayers),
      f_arr(NLayers),
      source_poly_coeffs(NLayers, Nscoeffs),
      Leg_coeffs_all(NLayers, NLeg_all),
      b_pos(NFourier, N),
      b_neg(NFourier, N),
      brdf_fourier_modes(NBDRF),
      // Derived data
      scale_tau(NLayers),
      scaled_omega_arr(NLayers),
      scaled_tau_arr_with_0(NLayers + 1),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      Leg_coeffs_residue_avg(NLeg_all),
      weighted_scaled_Leg_coeffs(NLayers, NLeg),
      weighted_Leg_coeffs_all(NLayers, NLeg_all),
      GC_collect(NFourier, NLayers, NQuad, NQuad),
      G_collect(NFourier, NLayers, NQuad, NQuad),
      K_collect(NFourier, NLayers, NQuad),
      expK_collect(NFourier, NLayers, NQuad),
      B_collect(NFourier, NLayers, NQuad),
      // Pure compute allocations
      n(NQuad * NLayers),
      RHS(n),
      jvec(NQuad),
      fac(NLeg),
      weighted_asso_Leg_coeffs_l(NLeg),
      asso_leg_term_mu0(NLeg),
      X_temp(NLeg),
      mathscr_X_pos(N),
      E_Lm1L(N),
      E_lm1l(N),
      E_llp1(N),
      BDRF_RHS_contribution(N),
      SRCB(NQuad),
      SRC0(NLayers, NQuad),
      SRC1(NLayers, NQuad),
      Gml(NQuad, NQuad),
      BDRF_LHS(N, NQuad),
      R(N, N),
      mathscr_D_neg(N, N),
      D_pos(N, N),
      D_neg(N, N),
      apb(N, N),
      amb(N, N),
      sqr(N, N),
      asso_leg_term_pos(N, NLeg),
      asso_leg_term_neg(N, NLeg),
      D_temp(N, NLeg),
      solve_work(NQuad),
      diag_work(N),
      LHSB(3 * N - 1, 3 * N - 1, n, n),
      comp_data(NQuad, Nscoeffs) {
  ARTS_TIME_REPORT

  Legendre::PositiveDoubleGaussLegendre(mu_arr[rf(N)], W);

  std::transform(
      mu_arr.begin(), mu_arr.begin() + N, mu_arr.begin() + N, [](auto&& x) {
        return -x;
      });
  std::transform(
      mu_arr.begin(), mu_arr.end(), inv_mu_arr.begin(), [](auto&& x) {
        return 1.0 / x;
      });
}

main_data::main_data(const Index NQuad_,
                     const Index NLeg_,
                     const Index NFourier_,
                     AscendingGrid tau_arr_,
                     Vector omega_arr_,
                     Matrix Leg_coeffs_all_,
                     Matrix b_pos_,
                     Matrix b_neg_,
                     Vector f_arr_,
                     Matrix source_poly_coeffs_,
                     std::vector<BDRF> brdf_fourier_modes_,
                     Numeric mu0_,
                     Numeric I0_,
                     Numeric phi0_)
    : NLayers(tau_arr_.size()),
      NQuad(NQuad_),
      NLeg(NLeg_),
      NFourier(NFourier_),
      N(NQuad / 2),
      Nscoeffs(source_poly_coeffs_.ncols()),
      NLeg_all(Leg_coeffs_all_.ncols()),
      NBDRF(brdf_fourier_modes_.size()),
      has_source_poly(Nscoeffs > 0),
      has_beam_source(I0_ > 0),
      // User data
      tau_arr(std::move(tau_arr_)),
      omega_arr(std::move(omega_arr_)),
      f_arr(std::move(f_arr_)),
      source_poly_coeffs(std::move(source_poly_coeffs_)),
      Leg_coeffs_all(std::move(Leg_coeffs_all_)),
      b_pos(std::move(b_pos_)),
      b_neg(std::move(b_neg_)),
      brdf_fourier_modes(std::move(brdf_fourier_modes_)),
      mu0(mu0_),
      I0(I0_),
      phi0(phi0_),
      // Derived data
      scale_tau(NLayers),
      scaled_omega_arr(NLayers),
      scaled_tau_arr_with_0(NLayers + 1),
      mu_arr(NQuad),
      inv_mu_arr(NQuad),
      W(N),
      Leg_coeffs_residue_avg(NLeg_all),
      weighted_scaled_Leg_coeffs(NLayers, NLeg),
      weighted_Leg_coeffs_all(NLayers, NLeg_all),
      GC_collect(NFourier, NLayers, NQuad, NQuad),
      G_collect(NFourier, NLayers, NQuad, NQuad),
      K_collect(NFourier, NLayers, NQuad),
      expK_collect(NFourier, NLayers, NQuad),
      B_collect(NFourier, NLayers, NQuad),
      // Pure compute allocations
      n(NQuad * NLayers),
      RHS(n),
      jvec(NQuad),
      fac(NLeg),
      weighted_asso_Leg_coeffs_l(NLeg),
      asso_leg_term_mu0(NLeg),
      X_temp(NLeg),
      mathscr_X_pos(N),
      E_Lm1L(N),
      E_lm1l(N),
      E_llp1(N),
      BDRF_RHS_contribution(N),
      SRCB(NQuad),
      SRC0(NLayers, NQuad),
      SRC1(NLayers, NQuad),
      Gml(NQuad, NQuad),
      BDRF_LHS(N, NQuad),
      R(N, N),
      mathscr_D_neg(N, N),
      D_pos(N, N),
      D_neg(N, N),
      apb(N, N),
      amb(N, N),
      sqr(N, N),
      asso_leg_term_pos(N, NLeg),
      asso_leg_term_neg(N, NLeg),
      D_temp(N, NLeg),
      solve_work(NQuad),
      diag_work(N),
      LHSB(3 * N - 1, 3 * N - 1, n, n),
      comp_data(NQuad, Nscoeffs) {
  ARTS_TIME_REPORT

  Legendre::PositiveDoubleGaussLegendre(mu_arr[rf(N)], W);

  std::transform(
      mu_arr.begin(), mu_arr.begin() + N, mu_arr.begin() + N, [](auto&& x) {
        return -x;
      });
  std::transform(
      mu_arr.begin(), mu_arr.end(), inv_mu_arr.begin(), [](auto&& x) {
        return 1.0 / x;
      });

  check_input_size();
  update_all(I0_);
}

[[nodiscard]] Index main_data::tau_index(const Numeric tau) const {
  ARTS_TIME_REPORT

  const Index l =
      std::distance(tau_arr.begin(), stdr::lower_bound(tau_arr, tau));
  ARTS_USER_ERROR_IF(
      l == NLayers, "tau ({}) must be at most {}", tau, tau_arr.back());
  return l;
}

void main_data::u(u_data& data, const Numeric tau, const Numeric phi) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.exponent.resize(NFourier, NQuad);
  for (Index i = 0; i < NFourier; i++) {
    for (Index j = 0; j < N; j++) {
      data.exponent[i, j] =
          std::exp(K_collect[i, l, j] * (scaled_tau - scaled_tau_arr_lm1));
      data.exponent[i, j + N] =
          std::exp(K_collect[i, l, j + N] * (scaled_tau - scaled_tau_arr_l));
    }
  }

  data.um.resize(NFourier, NQuad);
  static_assert(
      matpack::einsum_optpath<"mi", "mij", "mj">(),
      "On Failure, the einsum has been changed to not use optimal path");
  einsum<"mi", "mij", "mj">(
      data.um, GC_collect[joker, l, joker, joker], data.exponent);

  if (has_beam_source) {
    for (Index m = 0; m < NFourier; m++) {
      for (Index i = 0; i < NQuad; i++) {
        data.um[m, i] += std::exp(-scaled_tau / mu0) * B_collect[m, l, i];
      }
    }
  }

  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.um[0],
              data.src,
              tau,
              omega_arr[l],
              source_poly_coeffs[l],
              G_collect[0, l],
              K_collect[0, l],
              inv_mu_arr,
              0,
              1.0,
              1.0);
  }

  data.intensities.resize(NQuad);
  data.intensities = 0.0;
  for (Index m = 0; m < NFourier; m++) {
    const Numeric cp = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const auto umm   = data.um[m];
    for (Index i = 0; i < NQuad; i++) {
      data.intensities[i] += umm[i] * cp;
    }
  }

  data.intensities *= I0_orig;
}

void main_data::u0(u0_data& data, const Numeric tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.exponent.resize(NQuad);
  for (Index j = 0; j < N; j++) {
    data.exponent[j] =
        std::exp(K_collect[0, l, j] * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[j + N] =
        std::exp(K_collect[0, l, j + N] * (scaled_tau - scaled_tau_arr_l));
  }

  data.u0.resize(NQuad);
  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.u0,
              data.src,
              tau,
              omega_arr[l],
              source_poly_coeffs[l],
              G_collect[0, l],
              K_collect[0, l],
              inv_mu_arr);
  } else {
    data.u0 = 0.0;
  }

  mult(data.u0, GC_collect[0, l, joker, joker], data.exponent, 1.0, 1.0);

  if (has_beam_source) {
    const auto tmp = B_collect[0, l, joker];
    std::transform(tmp.elem_begin(),
                   tmp.elem_end(),
                   data.u0.elem_begin(),
                   data.u0.elem_begin(),
                   [scl = std::exp(-scaled_tau / mu0)](auto&& x, auto&& y) {
                     return scl * x + y;
                   });
  }

  data.u0 *= I0_orig;
}

namespace {
Numeric calculate_nu(const Numeric mu,
                     const Numeric phi,
                     const Numeric mu_p,
                     const Numeric phi_p) {
  const Numeric scl = std::sqrt(1.0 - mu_p * mu_p) * std::cos(phi_p - phi);
  return mu * mu_p + scl * std::sqrt(1.0 - mu * mu);
}

void calculate_nu(Vector& nu,
                  const ConstVectorView& mu,
                  const Numeric phi,
                  const Numeric mu_p,
                  const Numeric phi_p) {
  nu.resize(mu.size());

  std::transform(
      mu.begin(),
      mu.end(),
      nu.begin(),
      [mu_p, scl = std::sqrt(1.0 - mu_p * mu_p) * std::cos(phi_p - phi)](
          auto&& x) { return x * mu_p + scl * std::sqrt(1.0 - x * x); });
}
}  // namespace

void main_data::TMS(tms_data& data,
                    const Numeric tau,
                    const Numeric phi) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  // mathscr_B
  calculate_nu(data.nu, mu_arr, phi, -mu0, phi0);

  data.mathscr_B.resize(NLayers, NQuad);
  for (Index j = 0; j < NLayers; j++) {
    for (Index i = 0; i < NQuad; i++) {
      const Numeric p_true =
          Legendre::legendre_sum(weighted_Leg_coeffs_all[j], data.nu[i]);
      const Numeric p_trun =
          Legendre::legendre_sum(weighted_scaled_Leg_coeffs[j], data.nu[i]);
      data.mathscr_B[j, i] = (scaled_omega_arr[j] * I0) / (4 * Constant::pi) *
                             (mu0 / (mu0 + mu_arr[i])) *
                             (p_true / (1.0 - f_arr[j]) - p_trun);
    }
  }

  data.TMS.resize(NQuad);
  const Numeric exptau = std::exp(-scaled_tau / mu0);
  for (Index i = 0; i < N; i++) {
    data.TMS[i] =
        exptau - std::exp((scaled_tau - scaled_tau_arr_l) / mu_arr[i] -
                          scaled_tau_arr_l / mu0);
    data.TMS[i + N] =
        exptau - std::exp((scaled_tau_arr_lm1 - scaled_tau) / mu_arr[i] -
                          scaled_tau_arr_lm1 / mu0);
  }
  data.TMS *= data.mathscr_B[l];

  if (tau_arr.size() > 1) {
    data.contribution_from_other_layers_pos.resize(N, NLayers);
    data.contribution_from_other_layers_pos = 0;
    data.contribution_from_other_layers_neg.resize(N, NLayers);
    data.contribution_from_other_layers_neg = 0;
    for (Index i = 0; i < NLayers; i++) {
      if (l > i) {
        // neg
        for (Index j = 0; j < N; j++) {
          data.contribution_from_other_layers_neg[j, i] =
              (data.mathscr_B[i, j + N] *
               (std::exp((scaled_tau_arr_with_0[i] - scaled_tau) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0) -
                std::exp((scaled_tau_arr_with_0[i] - scaled_tau) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0)));
        }
      } else if (l < i) {
        // pos
        for (Index j = 0; j < N; j++) {
          data.contribution_from_other_layers_pos[j, i] =
              (data.mathscr_B[i, j] *
               (std::exp((scaled_tau - scaled_tau_arr_with_0[i]) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0) -
                std::exp((scaled_tau - scaled_tau_arr_with_0[i]) / mu_arr[j] -
                         scaled_tau_arr_with_0[i] / mu0)));
        }
      } else {
        continue;
      }
    }

    for (Index i = 0; i < N; i++) {
      data.TMS[i + 0] += sum(data.contribution_from_other_layers_pos[i]);
      data.TMS[i + N] += sum(data.contribution_from_other_layers_neg[i]);
    }
  }
}

void main_data::IMS(Vector& ims, const Numeric tau, const Numeric phi) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);

  ims.resize(N);
  for (Index i = 0; i < N; i++) {
    const Numeric nu  = calculate_nu(mu_arr[i + N], phi, -mu0, phi0);
    const Numeric x   = 1.0 / mu_arr[i] - 1.0 / scaled_mu0;
    const Numeric chi = (1 / (mu_arr[i] * scaled_mu0 * x)) *
                        ((tau - 1 / x) * std::exp(-tau / scaled_mu0) +
                         std::exp(-tau / mu_arr[i]) / x);
    ims[i] = (I0 / (4 * Constant::pi) * Math::pow2(omega_avg * f_avg) /
              (1 - omega_avg * f_avg) *
              Legendre::legendre_sum(Leg_coeffs_residue_avg, nu)) *
             chi;
  }
}

void main_data::u_corr(u_data& u_data,
                       Vector& ims,
                       tms_data& tms_data,
                       const Numeric tau,
                       const Numeric phi) const {
  ARTS_TIME_REPORT

  TMS(tms_data, tau, phi);
  IMS(ims, tau, phi);
  u(u_data, tau, phi);

  for (Index i = 0; i < N; i++) {
    u_data.intensities[i] += I0_orig * tms_data.TMS[i];
  }

  for (Index i = N; i < NQuad; i++) {
    u_data.intensities[i] += I0_orig * (tms_data.TMS[i] + ims[i - N]);
  }
}

//! FIXME: This implementation should be improved
void main_data::gridded_TMS(Tensor3View tms, const Vector& phi) const {
  ARTS_TIME_REPORT

  const Index M = phi.size();
  tms_data t{};

  for (Index l = 0; l < NLayers; l++) {
    for (Index j = 0; j < M; j++) {
      TMS(t, tau_arr[l], phi[j]);
      tms[l, j] = t.TMS;
    }
  }
}

void main_data::gridded_IMS(Tensor3View ims, const Vector& phi) const {
  ARTS_TIME_REPORT

  const Index M = phi.size();

  for (Index l = 0; l < NLayers; l++) {
    for (Index j = 0; j < M; j++) {
      for (Index i = 0; i < N; i++) {
        const Numeric nu = calculate_nu(mu_arr[i + N], phi[j], -mu0, phi0);
        const Numeric x  = 1.0 / mu_arr[i] - 1.0 / scaled_mu0;
        const Numeric chi =
            (1 / (mu_arr[i] * scaled_mu0 * x)) *
            ((tau_arr[l] - 1.0 / x) * std::exp(-tau_arr[l] / scaled_mu0) +
             std::exp(-tau_arr[l] / mu_arr[i]) / x);
        ims[l, j, i] =
            (I0 / (4 * Constant::pi) * Math::pow2(omega_avg * f_avg) /
             (1 - omega_avg * f_avg) *
             Legendre::legendre_sum(Leg_coeffs_residue_avg, nu)) *
            chi;
      }
    }
  }
}

void main_data::gridded_u_corr(Tensor3View u_data,
                               Tensor3View tms,
                               Tensor3View ims,
                               const Vector& phi) const {
  ARTS_TIME_REPORT

  gridded_u(u_data, phi);

  if (has_beam_source) {
    gridded_TMS(tms, phi);
    gridded_IMS(ims, phi);

    tms                         *= I0_orig;
    u_data                      += tms;
    ims                         *= I0_orig;
    u_data[joker, joker, rb(N)] += ims;
  }
}

Numeric main_data::flux_up(flux_data& data, const Numeric tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);
  ARTS_USER_ERROR_IF(tau < 0,
                     "tau ({}) must be less than the last layer ({})",
                     tau,
                     tau_arr.back());

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.u0_pos.resize(N);
  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.u0_pos,
              data.src,
              tau,
              omega_arr[l],
              source_poly_coeffs[l],
              G_collect[0, l],
              K_collect[0, l],
              inv_mu_arr);
  } else {
    data.u0_pos = 0.0;
  }

  if (has_beam_source) {
    for (Index i = 0; i < N; i++) {
      data.u0_pos[i] += B_collect[0, l, i] * std::exp(-scaled_tau / mu0);
    }
  }

  data.exponent.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.exponent[i] =
        std::exp(K_collect[0, l, i] * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[i + N] =
        std::exp(K_collect[0, l, i + N] * (scaled_tau - scaled_tau_arr_l));
  }

  mult(data.u0_pos, GC_collect[0, l, rf(N), joker], data.exponent, 1.0, 1.0);

  return Constant::two_pi * I0_orig *
         einsum<Numeric, "", "i", "i", "i">({}, mu_arr[rf(N)], W, data.u0_pos);
}

std::pair<Numeric, Numeric> main_data::flux_down(flux_data& data,
                                                 const Numeric tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(tau < 0, "tau ({}) must be positive", tau);
  ARTS_USER_ERROR_IF(tau < 0,
                     "tau ({}) must be less than the last layer ({})",
                     tau,
                     tau_arr.back());

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.u0_neg.resize(N);
  if (has_source_poly) {
    data.src.resize(NQuad, Nscoeffs);
    mathscr_v(data.u0_neg,
              data.src,
              tau,
              omega_arr[l],
              source_poly_coeffs[l],
              G_collect[0, l],
              K_collect[0, l],
              inv_mu_arr,
              N);
  } else {
    data.u0_neg = 0.0;
  }

  const Numeric direct_beam =
      has_beam_source ? I0 * mu0 * std::exp(-tau / mu0) : 0;
  const Numeric direct_beam_scaled =
      has_beam_source ? I0 * mu0 * std::exp(-scaled_tau / mu0) : 0;
  if (has_beam_source) {
    for (Index i = 0; i < N; i++) {
      data.u0_neg[i] += B_collect[0, l, i + N] * std::exp(-scaled_tau / mu0);
    }
  }

  data.exponent.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.exponent[i] =
        std::exp(K_collect[0, l, i] * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[i + N] =
        std::exp(K_collect[0, l, i + N] * (scaled_tau - scaled_tau_arr_l));
  }

  mult(data.u0_neg, GC_collect[0, l, rb(N), joker], data.exponent, 1.0, 1.0);

  return {I0_orig * (Constant::two_pi * einsum<Numeric, "", "i", "i", "i">(
                                            {}, mu_arr[rf(N)], W, data.u0_neg) -
                     direct_beam + direct_beam_scaled),
          I0_orig * I0 * direct_beam};
}

void main_data::gridded_flux(VectorView flux_up,
                             VectorView flux_do,
                             VectorView flux_dd) const try {
  ARTS_TIME_REPORT

  Vector u0(NQuad);
  Vector exponent(NQuad, 1);

  for (Index l = 0; l < NLayers; l++) {
    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];

    if (has_source_poly) {
      u0 = SRC0[l];
    } else {
      u0 = 0.0;
    }

    const Numeric direct_beam =
        has_beam_source ? I0 * mu0 * std::exp(-tau_arr[l] / mu0) : 0;
    const Numeric direct_beam_scaled =
        has_beam_source ? I0 * mu0 * std::exp(-scaled_tau_arr_l / mu0) : 0;
    if (has_beam_source) {
      for (Index i = 0; i < NQuad; i++) {
        u0[i] += B_collect[0, l, i] * std::exp(-scaled_tau_arr_l / mu0);
      }
    }

    exponent[rf(N)] = expK_collect[0, l, rf(N)];
    mult(u0, GC_collect[0, l], exponent, 1.0, 1.0);

    flux_up[l] =
        Constant::two_pi * I0_orig *
        einsum<Numeric, "", "i", "i", "i">({}, mu_arr[rf(N)], W, u0[rf(N)]);
    flux_do[l] =
        I0_orig * (Constant::two_pi * einsum<Numeric, "", "i", "i", "i">(
                                          {}, mu_arr[rf(N)], W, u0[rb(N)]) -
                   direct_beam + direct_beam_scaled);
    flux_dd[l] = I0_orig * I0 * direct_beam;
  }
}
ARTS_METHOD_ERROR_CATCH

void main_data::gridded_u(Tensor3View out, const Vector& phi) const {
  ARTS_TIME_REPORT

  Matrix exponent(NFourier, NQuad, 1);
  Matrix um(NFourier, NQuad);

  const Index Nphi = phi.size();
  Matrix cp(Nphi, NFourier);
  for (Size p = 0; p < phi.size(); p++) {
    for (Index m = 0; m < NFourier; m++) {
      cp[p, m] = I0_orig * std::cos(static_cast<Numeric>(m) * (phi0 - phi[p]));
    }
  }

  for (Index l = 0; l < NLayers; l++) {
    const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];

    static_assert(
        matpack::einsum_optpath<"mi", "mij", "mj">(),
        "On Failure, the einsum has been changed to not use optimal path");
    exponent[joker, rf(N)] = expK_collect[joker, l, rf(N)];
    einsum<"mi", "mij", "mj">(um, GC_collect[joker, l], exponent);

    if (has_beam_source) {
      for (Index m = 0; m < NFourier; m++) {
        for (Index i = 0; i < NQuad; i++) {
          um[m, i] += std::exp(-scaled_tau_arr_l / mu0) * B_collect[m, l, i];
        }
      }
    }

    if (has_source_poly) um[0] += SRC0[l];

    static_assert(
        matpack::einsum_optpath<"pi", "im", "pm">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"pi", "im", "pm">(out[l], transpose(um), cp);
  }
}

void main_data::ungridded_flux(VectorView flux_up,
                               VectorView flux_do,
                               VectorView flux_dd,
                               const AscendingGrid& tau) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      tau.front() < 0, "the first tau ({}) must be positive", tau.front());
  ARTS_USER_ERROR_IF(tau.back() < tau_arr.back(),
                     "the last tau ({}) must be less than the last layer ({})",
                     tau.back(),
                     tau_arr.back());

  Vector u0(NQuad);
  Vector exponent(NQuad, 1);
  mathscr_v_data src(NQuad, Nscoeffs);

  Index l = tau_index(tau.front());
  for (Size il = 0; il < tau.size(); il++) {
    while (tau[il] > tau_arr[l]) l++;

    const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
    const Numeric scaled_tau =
        scaled_tau_arr_l - (tau_arr[l] - tau[il]) * scale_tau[l];

    if (has_source_poly) {
      mathscr_v(u0,
                src,
                tau[il],
                omega_arr[l],
                source_poly_coeffs[l],
                G_collect[0, l],
                K_collect[0, l],
                inv_mu_arr);
    } else {
      u0 = 0.0;
    }

    const Numeric direct_beam =
        has_beam_source ? I0 * mu0 * std::exp(-tau[il] / mu0) : 0;
    const Numeric direct_beam_scaled =
        has_beam_source ? I0 * mu0 * std::exp(-scaled_tau / mu0) : 0;
    if (has_beam_source) {
      for (Index i = 0; i < NQuad; i++) {
        u0[i] += B_collect[0, l, i] * std::exp(-scaled_tau / mu0);
      }
    }

    for (Index i = 0; i < N; i++) {
      exponent[i] =
          std::exp(K_collect[0, l, i] * (scaled_tau - scaled_tau_arr_lm1));
      exponent[i + N] =
          std::exp(K_collect[0, l, i + N] * (scaled_tau - scaled_tau_arr_l));
    }

    mult(u0, GC_collect[0, l, joker, joker], exponent, 1.0, 1.0);

    flux_up[il] =
        Constant::two_pi * I0_orig *
        einsum<Numeric, "", "i", "i", "i">({}, mu_arr[rf(N)], W, u0[rf(N)]);
    flux_do[il] =
        I0_orig * (Constant::two_pi * einsum<Numeric, "", "i", "i", "i">(
                                          {}, mu_arr[rf(N)], W, u0[rb(N)]) -
                   direct_beam + direct_beam_scaled);
    flux_dd[il] = I0_orig * I0 * direct_beam;
  }
}

void main_data::ungridded_u(Tensor3View out,
                            const AscendingGrid& tau,
                            const Vector& phi) const {
  ARTS_TIME_REPORT

  ARTS_USER_ERROR_IF(
      tau.front() < 0, "the first tau ({}) must be positive", tau.front());
  ARTS_USER_ERROR_IF(tau.back() < tau_arr.back(),
                     "the last tau ({}) must be less than the last layer ({})",
                     tau.back(),
                     tau_arr.back());

  Matrix exponent(NFourier, NQuad, 1);
  Matrix um(NFourier, NQuad);
  mathscr_v_data src(NQuad, Nscoeffs);

  const Index Nphi = phi.size();
  Matrix cp(Nphi, NFourier);
  for (Size p = 0; p < phi.size(); p++) {
    for (Index m = 0; m < NFourier; m++) {
      cp[p, m] = I0_orig * std::cos(static_cast<Numeric>(m) * (phi0 - phi[p]));
    }
  }

  Index l = tau_index(tau.front());
  for (Size il = 0; il < tau.size(); il++) {
    while (tau[il] > tau_arr[l]) l++;

    const Numeric scaled_tau_arr_l   = scaled_tau_arr_with_0[l + 1];
    const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
    const Numeric scaled_tau =
        scaled_tau_arr_l - (tau_arr[l] - tau[il]) * scale_tau[l];

    for (Index i = 0; i < NFourier; i++) {
      for (Index j = 0; j < N; j++) {
        exponent[i, j] =
            std::exp(K_collect[i, l, j] * (scaled_tau - scaled_tau_arr_lm1));
        exponent[i, j + N] =
            std::exp(K_collect[i, l, j + N] * (scaled_tau - scaled_tau_arr_l));
      }
    }

    static_assert(
        matpack::einsum_optpath<"mi", "mij", "mj">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"mi", "mij", "mj">(um, GC_collect[joker, l, joker, joker], exponent);

    if (has_beam_source) {
      for (Index m = 0; m < NFourier; m++) {
        for (Index i = 0; i < NQuad; i++) {
          um[m, i] += std::exp(-scaled_tau / mu0) * B_collect[m, l, i];
        }
      }
    }

    if (has_source_poly) {
      mathscr_v(um[0],
                src,
                tau[il],
                omega_arr[l],
                source_poly_coeffs[l],
                G_collect[0, l],
                K_collect[0, l],
                inv_mu_arr,
                0,
                1.0,
                1.0);
    }

    static_assert(
        matpack::einsum_optpath<"pi", "im", "pm">(),
        "On Failure, the einsum has been changed to not use optimal path");
    einsum<"pi", "im", "pm">(out[il], transpose(um), cp);
  }
}
}  // namespace disort

void DisortSettings::resize(Index quadrature_dimension_,
                            Index legendre_polynomial_dimension_,
                            Index fourier_mode_dimension_,
                            Index nfreq_,
                            Index nlay_) {
  quadrature_dimension          = quadrature_dimension_;
  legendre_polynomial_dimension = legendre_polynomial_dimension_;
  fourier_mode_dimension        = fourier_mode_dimension_;
  nfreq                         = nfreq_;
  nlay                          = nlay_;

  solar_source.resize(nfreq);
  solar_zenith_angle.resize(nfreq);
  solar_azimuth_angle.resize(nfreq);
  bidirectional_reflectance_distribution_functions.resize(nfreq, 0);
  optical_thicknesses.resize(nfreq, nlay);
  single_scattering_albedo.resize(nfreq, nlay);
  fractional_scattering.resize(nfreq, nlay);
  source_polynomial.resize(nfreq, nlay, 0);
  legendre_coefficients.resize(nfreq, nlay, legendre_polynomial_dimension);
  positive_boundary_condition.resize(
      nfreq, fourier_mode_dimension, quadrature_dimension / 2);
  negative_boundary_condition.resize(
      nfreq, fourier_mode_dimension, quadrature_dimension / 2);
}

void DisortSettings::check() const {
  ARTS_USER_ERROR_IF(
      solar_source.shape() != std::array{nfreq} or
          solar_zenith_angle.shape() != std::array{nfreq} or
          solar_azimuth_angle.shape() != std::array{nfreq} or
          (bidirectional_reflectance_distribution_functions.shape() !=
           std::array{
               nfreq,
               bidirectional_reflectance_distribution_functions.ncols()}) or
          (optical_thicknesses.shape() != std::array{nfreq, nlay}) or
          (single_scattering_albedo.shape() != std::array{nfreq, nlay}) or
          (fractional_scattering.shape() != std::array{nfreq, nlay}) or
          (source_polynomial.shape() !=
           std::array{nfreq, nlay, source_polynomial.ncols()}) or
          (legendre_coefficients.shape() !=
           std::array{nfreq, nlay, legendre_coefficients.ncols()}) or
          (positive_boundary_condition.shape() !=
           std::array{
               nfreq, fourier_mode_dimension, quadrature_dimension / 2}) or
          (negative_boundary_condition.shape() !=
           std::array{
               nfreq, fourier_mode_dimension, quadrature_dimension / 2}) or
          legendre_polynomial_dimension > legendre_coefficients.ncols(),
      R"-x-(Input is incorrect.

{:s}

Also note that the reduced Legendre polynomial dimension is {}.  It must be at most {}.
)-x-",
      *this,
      legendre_polynomial_dimension,
      legendre_coefficients.ncols());
}

disort::main_data DisortSettings::init() const try {
  check();
  return disort::main_data(
      nlay,
      quadrature_dimension,
      legendre_coefficients.ncols(),
      fourier_mode_dimension,
      source_polynomial.ncols(),
      legendre_polynomial_dimension,
      bidirectional_reflectance_distribution_functions.ncols());
}
ARTS_METHOD_ERROR_CATCH

disort::main_data& DisortSettings::set(disort::main_data& dis, Index iv) const
    try {
  using Conversion::cosd;
  using Conversion::deg2rad;

  for (Index i = 0;
       i < bidirectional_reflectance_distribution_functions.ncols();
       i++) {
    dis.brdf_modes()[i] =
        bidirectional_reflectance_distribution_functions[iv, i];
  }

  dis.tau(optical_thicknesses[iv]);
  dis.solar_zenith()        = cosd(solar_zenith_angle[iv]);
  dis.beam_azimuth()        = deg2rad(solar_azimuth_angle[iv]);
  dis.omega()               = single_scattering_albedo[iv];
  dis.f()                   = fractional_scattering[iv];
  dis.all_legendre_coeffs() = legendre_coefficients[iv];
  dis.positive_boundary()   = positive_boundary_condition[iv];
  dis.negative_boundary()   = negative_boundary_condition[iv];
  dis.source_poly()         = source_polynomial[iv];

  dis.update_all(solar_source[iv]);

  return dis;
}
ARTS_METHOD_ERROR_CATCH

disort::main_data& DisortSettings::set_cdisort(disort::main_data& dis, Index iv) const
    try {
  using Conversion::cosd;
  using Conversion::deg2rad;

  for (Index i = 0;
       i < bidirectional_reflectance_distribution_functions.ncols();
       i++) {
    dis.brdf_modes()[i] =
        bidirectional_reflectance_distribution_functions[iv, i];
  }

  dis.tau(optical_thicknesses[iv]);
  dis.solar_zenith()        = cosd(solar_zenith_angle[iv]);
  dis.beam_azimuth()        = deg2rad(solar_azimuth_angle[iv]);
  dis.omega()               = single_scattering_albedo[iv];
  dis.f()                   = fractional_scattering[iv];
  dis.all_legendre_coeffs() = legendre_coefficients[iv];
  dis.positive_boundary()   = positive_boundary_condition[iv];
  dis.negative_boundary()   = negative_boundary_condition[iv];

  return dis;
}
ARTS_METHOD_ERROR_CATCH

void xml_io_stream<DisortBDRF>::read(std::istream& is,
                                     DisortBDRF& x,
                                     bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.f, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<DisortBDRF>::write(std::ostream& os,
                                      const DisortBDRF& x,
                                      bofstream* pbofs,
                                      std::string_view name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.f, pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<DisortSettings>::read(std::istream& is_xml,
                                         DisortSettings& v,
                                         bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is_xml);
  tag.check_name(type_name);

  xml_read_from_stream(is_xml, v.quadrature_dimension, pbifs);
  xml_read_from_stream(is_xml, v.legendre_polynomial_dimension, pbifs);
  xml_read_from_stream(is_xml, v.fourier_mode_dimension, pbifs);
  xml_read_from_stream(is_xml, v.nfreq, pbifs);
  xml_read_from_stream(is_xml, v.nlay, pbifs);
  xml_read_from_stream(is_xml, v.solar_azimuth_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_zenith_angle, pbifs);
  xml_read_from_stream(is_xml, v.solar_source, pbifs);
  xml_read_from_stream(
      is_xml, v.bidirectional_reflectance_distribution_functions, pbifs);
  xml_read_from_stream(is_xml, v.optical_thicknesses, pbifs);
  xml_read_from_stream(is_xml, v.single_scattering_albedo, pbifs);
  xml_read_from_stream(is_xml, v.fractional_scattering, pbifs);
  xml_read_from_stream(is_xml, v.source_polynomial, pbifs);
  xml_read_from_stream(is_xml, v.legendre_coefficients, pbifs);
  xml_read_from_stream(is_xml, v.positive_boundary_condition, pbifs);
  xml_read_from_stream(is_xml, v.negative_boundary_condition, pbifs);

  tag.read_from_stream(is_xml);
  tag.check_end_name(type_name);
}

void xml_io_stream<DisortSettings>::write(std::ostream& os_xml,
                                          const DisortSettings& v,
                                          bofstream* pbofs,
                                          std::string_view) {
  XMLTag open_tag;
  XMLTag close_tag;

  open_tag.name = type_name;
  open_tag.write_to_stream(os_xml);

  xml_write_to_stream(
      os_xml, v.quadrature_dimension, pbofs, "quadrature_dimension");
  xml_write_to_stream(os_xml,
                      v.legendre_polynomial_dimension,
                      pbofs,
                      "legendre_polynomial_dimension");
  xml_write_to_stream(
      os_xml, v.fourier_mode_dimension, pbofs, "fourier_mode_dimension");
  xml_write_to_stream(os_xml, v.nfreq, pbofs, "nfreq");
  xml_write_to_stream(os_xml, v.nlay, pbofs, "nlay");
  xml_write_to_stream(os_xml, v.solar_azimuth_angle, pbofs, "solaz");
  xml_write_to_stream(os_xml, v.solar_zenith_angle, pbofs, "solza");
  xml_write_to_stream(os_xml, v.solar_source, pbofs, "solsrc");
  xml_write_to_stream(os_xml,
                      v.bidirectional_reflectance_distribution_functions,
                      pbofs,
                      "BRDF");
  xml_write_to_stream(os_xml, v.optical_thicknesses, pbofs, "Tau");
  xml_write_to_stream(os_xml, v.single_scattering_albedo, pbofs, "albedo");
  xml_write_to_stream(
      os_xml, v.fractional_scattering, pbofs, "fractional_scattering");
  xml_write_to_stream(os_xml, v.source_polynomial, pbofs, "source_polynomial");
  xml_write_to_stream(
      os_xml, v.legendre_coefficients, pbofs, "legendre_coefficients");
  xml_write_to_stream(os_xml,
                      v.positive_boundary_condition,
                      pbofs,
                      "positive_boundary_condition");
  xml_write_to_stream(os_xml,
                      v.negative_boundary_condition,
                      pbofs,
                      "negative_boundary_condition");

  close_tag.name = type_name;
  close_tag.write_to_end_stream(os_xml);
}

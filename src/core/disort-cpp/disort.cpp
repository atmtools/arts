#include "disort.h"

#include <gsl_gauss_legendre.h>

#include <algorithm>
#include <functional>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "arts_constants.h"
#include "configtypes.h"
#include "debug.h"
#include "legendre.h"
#include "matpack_data.h"
#include "matpack_einsum.h"
#include "matpack_iter.h"
#include "matpack_view.h"
#include "rational.h"

void print(auto&& x) {
  auto ptr = x.elem_begin();
  while (ptr != x.elem_end()) {
    std::cout << *ptr << ",";
    ptr++;
  }
  std::cout << '\n';
}

namespace disort {
void mathscr_v_add(ExhaustiveVectorView um,
                   Matrix& mathscr_v_coeffs,
                   const Numeric tau,
                   const ExhaustiveConstVectorView& s_poly_coeffs,
                   const ExhaustiveConstMatrixView& G,
                   const ExhaustiveConstVectorView& K,
                   const ExhaustiveConstMatrixView& G_inv,
                   const ExhaustiveConstVectorView& mu_arr) {
  const Index Ni = G.nrows();
  const Index Nk = G.ncols();
  const Index Nc = s_poly_coeffs.size();
  const Index Nj = mu_arr.size();

  mathscr_v_coeffs.resize(Nc, Nk);
  mathscr_v_coeffs = 0.0;
  for (Index i : std::ranges::iota_view(0, Nc)) {
    auto ci = mathscr_v_coeffs[i];
    for (Index k = 0; k < 16; k++) {
      auto& cik = ci[k];
      for (Index j = 0; j <= i; j++) {
        const Numeric fs = (Legendre::factorial(Nc - 1 - j) /
                            Legendre::factorial(Nc - 1 - i)) *
                           s_poly_coeffs[1 - j];
        cik += fs * std::pow(K[k], -(i - j + 1));
      }
    }
  }

  for (Index i = 0; i < Ni; i++) {
    for (Index k = 0; k < Nk; k++) {
      for (Index c = 0; c < Nc; c++) {
        for (Index j = 0; j < Nj; j++) {
          um[i] += G(i, k) * std::pow(tau, Nc - 1 - c) *
                   mathscr_v_coeffs(c, k) * G_inv(k, j) / mu_arr[j];
        }
      }
    }
  }
}

void solve_for_coefs(Tensor4& GC_collect,
                     const Index NFourier,
                     const Tensor4& G_collect,
                     const Tensor3& K_collect,
                     const Tensor3& B_collect,
                     const Tensor3& G_inv_collect_0,
                     const Vector& tau_arr,
                     const Vector& scaled_tau_arr_with_0,
                     const Vector& mu_arr_pos,
                     const Vector& mu_arr_pos_times_W,
                     const Vector& mu_arr,
                     const Index N,
                     const Index NQuad,
                     const Index NLayers,
                     const Index NBDRF,
                     const bool multilayer_bool,
                     const std::vector<BDRF>& BDRF_Fourier_modes,
                     const Numeric mu0,
                     const Numeric I0,
                     const bool beam_source_bool,
                     const Matrix& b_pos,
                     const Matrix& b_neg,
                     const bool scalar_b_pos,
                     const bool scalar_b_neg,
                     const Matrix& s_poly_coeffs,
                     bool iso_source_bool) {
  GC_collect.resize(NFourier, NLayers, NQuad, NQuad);

  Matrix R(NBDRF ? N : 0, N);
  Matrix mathscr_D_neg(NBDRF ? N : 0, N);
  Vector mathscr_X_pos(NBDRF ? N : 0);
  const Vector mu_arr_neg = einsum<Vector, "i", "i", "">({N}, mu_arr_pos, -1.0);

  Vector b_pos_m(N);
  Vector b_neg_m(N);

  Matrix comp_matrix;
  Vector mathscr_v_temporary(NQuad);
  Vector mathscr_v_contribution(NQuad * NLayers);

  Vector BDRF_RHS_contribution(N);

  Matrix RHS_middle(NQuad, NLayers - 1);
  Vector RHS(NQuad * NLayers);

  Matrix C_m(NLayers, NQuad);
  Vector E_Lm1L(N);
  Vector E_lm1l(N);
  Vector E_llp1(N);

  Matrix BDRF_LHS_contribution_neg(N, N);
  Matrix BDRF_LHS_contribution_pos(N, N);
  Matrix LHS(NLayers * NQuad, NLayers * NQuad);

  for (Index m = 0; m < NFourier; m++) {
    const bool m_equals_0_bool = m == 0;
    const bool BDRF_bool = m < NBDRF;
    const auto G_collect_m = G_collect[m];
    const auto K_collect_m = K_collect[m];
    const auto B_collect_m =
        beam_source_bool ? B_collect[m] : ExhaustiveConstMatrixView{};

    if (BDRF_bool) {
      BDRF_Fourier_modes[m](mathscr_D_neg, mu_arr_pos, mu_arr_neg),
          einsum<"ij", "ij", "j", "">(
              R, mathscr_D_neg, mu_arr_pos_times_W, 1 + m_equals_0_bool);
      if (beam_source_bool) {
        BDRF_Fourier_modes[m](mathscr_X_pos.reshape_as(N, 1),
                              mu_arr_pos,
                              ExhaustiveConstVectorView{-mu0});
        mathscr_X_pos *= mu0 * I0 / Constant::pi;
      }
    }

    if (scalar_b_pos) {
      if (m_equals_0_bool) {
        b_pos_m = b_pos(0, 0);
      } else {
        b_pos_m = 0.0;
      }
    } else {
      b_pos_m = b_pos(joker, m);
    }

    if (scalar_b_neg) {
      if (m_equals_0_bool) {
        b_neg_m = b_neg(0, 0);
      } else {
        b_neg_m = 0.0;
      }
    } else {
      b_neg_m = b_neg(joker, m);
    }

    // Fill RHS
    RHS = 0.0;

    if (iso_source_bool and m_equals_0_bool) {
      mathscr_v_contribution.slice(0, N) = 0.0;
      mathscr_v_add(mathscr_v_contribution.slice(0, N),
                    comp_matrix,
                    0.0,
                    s_poly_coeffs[0],
                    G_collect_m[0].slice(N, N),
                    K_collect_m[0],
                    G_inv_collect_0[0],
                    mu_arr);
      mathscr_v_contribution.slice(0, N) *= -1;

      if (multilayer_bool) {
        for (Index l = 0; l < NLayers - 1; l++) {
          mathscr_v_temporary = 0;
          mathscr_v_add(mathscr_v_temporary,
                        comp_matrix,
                        tau_arr[l],
                        s_poly_coeffs[l + 1],
                        G_collect_m[l + 1],
                        K_collect_m[l + 1],
                        G_inv_collect_0[l + 1],
                        mu_arr);
          mathscr_v_contribution(Range(l + N, NQuad, NLayers - 1)) -=
              mathscr_v_temporary;
          mathscr_v_temporary = 0.0;
          mathscr_v_add(mathscr_v_temporary,
                        comp_matrix,
                        tau_arr[l],
                        s_poly_coeffs[l],
                        G_collect_m[l],
                        K_collect_m[l],
                        G_inv_collect_0[l],
                        mu_arr);
          mathscr_v_contribution(Range(l + N, NQuad, NLayers - 1)) +=
              mathscr_v_temporary;
        }
      }
      mathscr_v_contribution.slice(mathscr_v_contribution.size() - N, N) = 0;
      mathscr_v_add(
          mathscr_v_contribution.slice(mathscr_v_contribution.size() - N, N),
          comp_matrix,
          tau_arr.back(),
          s_poly_coeffs[NLayers - 1],
          G_collect_m[NLayers - 1].slice(N, N),
          K_collect_m[NLayers - 1],
          G_inv_collect_0[NLayers - 1],
          mu_arr);

      if (NBDRF > 0) {
        // Add R@m to m using lapack
        mathscr_v_temporary =
            mathscr_v_contribution.slice(mathscr_v_contribution.size() - N, N);
        mult(mathscr_v_contribution.slice(mathscr_v_contribution.size() - N, N),
             R,
             mathscr_v_temporary,
             1.0,
             1.0);
      }
    } else {
      mathscr_v_contribution = 0.0;
    }

    if (beam_source_bool) {
      if (BDRF_bool) {
        BDRF_RHS_contribution = mathscr_X_pos.reshape_as(N);
        mult(BDRF_RHS_contribution,
             R,
             B_collect_m[NLayers - 1].slice(0, N),
             1.0,
             1.0);
      } else {
        BDRF_RHS_contribution = 0.0;
      }

      if (multilayer_bool) {
        for (Index l = 0; l < NLayers - 1; l++) {
          for (Index j = 0; j < NQuad; j++) {
            RHS_middle(j, l) = (B_collect_m(l + 1, j) - B_collect_m(l, j)) *
                               std::exp(-mu0 * scaled_tau_arr_with_0[l + 1]);
          }
        }
      }

      RHS = mathscr_v_contribution;
      RHS.slice(N, RHS_middle.size()) += RHS_middle.flat_view();
      for (Index i = 0; i < N; i++) {
        RHS[i] += b_neg_m[i] - B_collect_m(0, N + i);
        RHS[RHS.size() - N + i] +=
            b_pos_m[i] +
            (BDRF_RHS_contribution[i] - B_collect_m(NLayers - 1, i)) *
                std::exp(-scaled_tau_arr_with_0.back() / mu0);
      }
    } else {
      RHS = mathscr_v_contribution;
      RHS.slice(0, N) += b_neg_m;
      RHS.slice(RHS.size() - N, N) += b_pos_m;
    }

    // Fill LHS (THIS SHOULD BE SPARSE BUT IS NOT BECAUSE ARTS HAS NO SPARSE SOLVE!)
    LHS = 0.0;

    const auto G_0_np = G_collect_m(0, Range(N, N), Range(N, N));
    const auto G_L_pn = G_collect_m(NLayers - 1, Range(0, N), Range(0, N));
    const auto G_L_nn = G_collect_m(NLayers - 1, Range(N, N), Range(0, N));
    const auto G_L_pp = G_collect_m(NLayers - 1, Range(0, N), Range(N, N));
    const auto G_L_np = G_collect_m(NLayers - 1, Range(N, N), Range(N, N));
    for (Index i = 0; i < N; i++) {
      E_Lm1L[i] =
          std::exp(K_collect_m(K_collect_m.nrows() - 1, i) *
                   (scaled_tau_arr_with_0[scaled_tau_arr_with_0.size() - 1] -
                    scaled_tau_arr_with_0[scaled_tau_arr_with_0.size() - 2]));
    }

    if (BDRF_bool) {
      mult(BDRF_LHS_contribution_neg, R, G_L_nn);
      mult(BDRF_LHS_contribution_pos, R, G_L_np);
    } else {
      BDRF_LHS_contribution_neg = 0;
      BDRF_LHS_contribution_pos = 0;
    }

    LHS(Range(0, N), Range(0, N)) = G_collect_m(0, Range(N, N), Range(0, N));

    for (Index i = 0; i < N; i++) {
      for (Index j = 0; j < N; j++) {
        LHS(i, N + j) = G_0_np(i, j) *
                        std::exp(K_collect_m(0, j) * scaled_tau_arr_with_0[1]);
        LHS(LHS.nrows() - N + i, LHS.ncols() - 2 * N + j) =
            (G_L_pn(i, j) - BDRF_LHS_contribution_neg(i, j)) * E_Lm1L[j];
        LHS(LHS.nrows() - N + i, LHS.ncols() - N + j) =
            G_L_pp(i, j) - BDRF_LHS_contribution_pos(i, j);
      }
    }

    for (Index l = 0; l < NLayers - 1; l++) {
      const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
      const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
      const Numeric scaled_tau_arr_lp1 = scaled_tau_arr_with_0[l + 2];
      // Postive eigenvalues
      const auto K_l_pos = K_collect_m[l].slice(N, N);
      const auto K_lp1_pos = K_collect_m[l + 1].slice(N, N);

      for (Index i = 0; i < N; i++) {
        E_lm1l[i] =
            std::exp(K_l_pos[i] * (scaled_tau_arr_lm1 - scaled_tau_arr_l));
        E_llp1[i] =
            std::exp(K_lp1_pos[i] * (scaled_tau_arr_l - scaled_tau_arr_lp1));
      }

      einsum<"ij", "ij", "j">(LHS(Range(N + l * NQuad, N), Range(l * NQuad, N)),
                              G_collect_m(l, Range(0, N), Range(0, N)),
                              E_lm1l);
      einsum<"ij", "ij", "j">(
          LHS(Range(2 * N + l * NQuad, N), Range(l * NQuad, N)),
          G_collect_m(l, Range(N, N), Range(0, N)),
          E_lm1l);
      einsum<"ij", "", "ij", "j">(
          LHS(Range(N + l * NQuad, N), Range(l * NQuad + 2 * NQuad - N, N)),
          -1,
          G_collect_m(l + 1, Range(0, N), Range(N, N)),
          E_llp1);
      einsum<"ij", "", "ij", "j">(
          LHS(Range(2 * N + l * NQuad, N), Range(l * NQuad + 2 * NQuad - N, N)),
          -1,
          G_collect_m(l + 1, Range(N, N), Range(N, N)),
          E_llp1);
      LHS(Range(N + l * NQuad, NQuad), Range(l * NQuad + N, N)) =
          G_collect_m(l, joker, Range(N, N));
      einsum<"ij", "", "ij">(
          LHS(Range(N + l * NQuad, NQuad), Range(l * NQuad + 2 * N, N)),
          -1,
          G_collect_m(l + 1, joker, Range(0, N)));
    }

    //print(RHS);

    inplace_solve(RHS, LHS);
    einsum<"ijm", "ijm", "im">(
        GC_collect[m], G_collect_m, RHS.reshape_as(NLayers, NQuad));
  }
}

Numeric poch(Numeric x, Numeric n) { return Legendre::tgamma_ratio(x + n, x); }

void diagonalize(Tensor4& G_collect_,
                 Tensor3& K_collect_,
                 Tensor3& B_collect_,
                 Tensor3& G_inv_collect_0,
                 const Index NFourier,
                 const Vector& scaled_omega_arr,
                 const Vector& mu_arr_pos,
                 const Vector& mu_arr,
                 const Vector& M_inv,
                 const Vector& W,
                 const Index N,
                 const Index NQuad,
                 const Index NLeg,
                 const Index NLayers,
                 const Matrix& weighted_scaled_Leg_coeffs,
                 const Numeric mu0,
                 const Numeric I0,
                 const bool beam_source_bool,
                 const bool iso_source_bool) {
  G_collect_.resize(NFourier, NLayers, NQuad, NQuad);
  K_collect_.resize(NFourier, NLayers, NQuad);
  B_collect_.resize(NFourier, NLayers, NQuad);
  G_inv_collect_0.resize(NLayers, NQuad, NQuad);

  auto G_collect = G_collect_.reshape_as(NFourier * NLayers, NQuad, NQuad);
  auto K_collect = K_collect_.reshape_as(NFourier * NLayers, NQuad);
  auto B_collect = B_collect_.reshape_as(NFourier * NLayers, NQuad);

  Tensor3 alpha_arr(NFourier * NLayers, N, N);
  Tensor3 beta_arr(NFourier * NLayers, N, N);
  Matrix X_tilde_arr(NFourier * NLayers, NQuad);
  std::vector<Index> no_shortcut_indices;
  no_shortcut_indices.reserve(NFourier * NLayers);
  std::vector<std::pair<Index, Index>> no_shortcut_indices_0;
  no_shortcut_indices_0.reserve(NLayers);

  Vector fac(NLeg);
  Vector signs(NLeg);
  Vector weighted_asso_Leg_coeffs_l(NLeg);

  Vector asso_leg_term_mu0(NLeg);
  Matrix asso_leg_term_pos(N, NLeg);
  Matrix asso_leg_term_neg(N, NLeg);

  Matrix D_temp(N, NLeg);
  Matrix D_pos(N, N);
  Matrix D_neg(N, N);

  Vector X_temp(NLeg);
  Vector X_pos(N);
  Vector X_neg(N);
  auto xpos = X_pos.reshape_as(1, N);
  auto xneg = X_neg.reshape_as(1, N);

  Matrix apb(N, N);
  Matrix amb(N, N);
  Matrix sqr(N, N);

  Matrix GmG(N, NQuad);
  Matrix P(N, N);
  Vector eigr(N);
  Vector eigi(N);
  Matrix G_inv(NQuad, NQuad);

  matpack::matpack_data<Index, 1> ells_all(NLeg);
  for (Index i = 0; i < NLeg; i++) {
    ells_all[i] = i;
  }

  Index ind = 0;
  for (Index m = 0; m < NFourier; m++) {
    // All resizeing:
    const auto ells = ells_all.slice(m, NLeg - m);
    D_temp.resize(N, NLeg - m);
    X_temp.resize(NLeg - m);
    auto xtemp = X_temp.reshape_as(1, NLeg - m);
    fac.resize(ells.size());
    asso_leg_term_pos.resize(NLeg - m, N);
    asso_leg_term_neg.resize(NLeg - m, N);
    asso_leg_term_mu0.resize(NLeg - m);
    weighted_asso_Leg_coeffs_l.resize(ells.size());
    signs.resize(NLeg - m);

    const bool m_equals_0_bool = (m == 0);

    for (auto i : ells) {
      fac[i - m] =
          poch(static_cast<Numeric>(i + m + 1), static_cast<Numeric>(-2 * m));
    }

    for (Index i = 0; i < NLeg - m; i++) {
      signs[i] = i % 2 ? -1 : 1;
    }

    for (Index i = m; i < NLeg; i++) {
      for (Index j = 0; j < N; j++) {
        asso_leg_term_pos(i - m, j) =
            Legendre::assoc_legendre(i, m, mu_arr_pos[j]);
        asso_leg_term_neg(i - m, j) =
            asso_leg_term_pos(i - m, j) * signs[i - m];
      }
      asso_leg_term_mu0[i - m] = Legendre::assoc_legendre(ells[i - m], m, -mu0);
    }

    const bool all_asso_leg_term_pos_finite =
        std::all_of(asso_leg_term_pos.elem_begin(),
                    asso_leg_term_pos.elem_end(),
                    [](auto& x) { return std::isfinite(x); });

    for (Index l = 0; l < NLayers; l++) {
      for (Index i = 0; i < NLeg - m; i++) {
        weighted_asso_Leg_coeffs_l[i] =
            weighted_scaled_Leg_coeffs(l, ells[i]) * fac[i];
      }

      const Numeric scaled_omega_l = scaled_omega_arr[l];

      if (all_asso_leg_term_pos_finite and
          std::any_of(weighted_asso_Leg_coeffs_l.elem_begin(),
                      weighted_asso_Leg_coeffs_l.elem_end(),
                      Cmp::gt(0))) {
        einsum<"ij", "j", "ji">(
            D_temp, weighted_asso_Leg_coeffs_l, asso_leg_term_pos);
        mult(D_pos, D_temp, asso_leg_term_pos);
        mult(D_neg, D_temp, asso_leg_term_neg);
        D_pos *= 0.5 * scaled_omega_l;
        D_neg *= 0.5 * scaled_omega_l;

        if (beam_source_bool) {
          einsum<"i", "i", "i", "">(
              X_temp,
              weighted_asso_Leg_coeffs_l,
              asso_leg_term_mu0,
              (scaled_omega_l * I0 * (2 - m_equals_0_bool) /
               (4 * Constant::pi)));

          mult(xpos, xtemp, asso_leg_term_pos);
          mult(xneg, xtemp, asso_leg_term_neg);
        }

        auto alpha = alpha_arr[ind];
        auto beta = beta_arr[ind];
        no_shortcut_indices.push_back(ind);

        einsum<"ij", "i", "ij", "j">(alpha, M_inv, D_pos, W);
        alpha.diagonal() -= M_inv;

        einsum<"ij", "i", "ij", "j">(beta, M_inv, D_neg, W);

        einsum<"i", "", "i", "i">(
            X_tilde_arr[ind].slice(0, N), -1, M_inv, X_pos);

        einsum<"i", "i", "i">(X_tilde_arr[ind].slice(N, N), M_inv, X_neg);

        if (iso_source_bool and m_equals_0_bool) {
          no_shortcut_indices_0.emplace_back(ind , l);
        }
      } else {
        auto G = G_collect[ind];
        for (Index i = 0; i < N; i++) {
          G(i + N, i) = 1;
          G(i, i + N) = 1;
        }
        for (Index i = 0; i < NQuad; i++) {
          K_collect(ind, i) = -1 / mu_arr[i];
        }
        if (iso_source_bool and m_equals_0_bool) {
          G_inv_collect_0[l] = G;
        }
      }
      ++ind;
    }
  }

  for (auto i : no_shortcut_indices) {
    auto K = K_collect[i];
    auto G = G_collect[i];

    apb = alpha_arr[i];
    amb = alpha_arr[i];
    apb += beta_arr[i];
    amb -= beta_arr[i];
    mult(sqr, amb, apb);

    auto GpG = G(Range(0, N), Range(0, N));
    ::diagonalize(GpG, eigr, eigi, sqr);
    G(Range(0, N), Range(N, N)) = GpG;

    for (Index j = 0; j < N; j++) {
      const Numeric sqrt_x = std::sqrt(eigr[j]);
      K[j] = -sqrt_x;
      K[j + N] = sqrt_x;
    }

    mult(GmG, apb, G.slice(0, N));
    for (Index j = 0; j < NQuad; j++) {
      GmG(joker, j) /= K[j];
    }

    G.slice(N, N) = G.slice(0, N);
    G.slice(0, N) -= GmG;
    G.slice(N, N) += GmG;
    G *= 0.5;

    // Get inverse by the layer
    inv(G_inv, G);
    if (not no_shortcut_indices_0.empty() and no_shortcut_indices_0.front().first == i) {
      G_inv_collect_0[no_shortcut_indices_0.front().second] = G_inv;
      no_shortcut_indices_0.erase(no_shortcut_indices_0.begin());
    }
    if (iso_source_bool) {
      auto B = B_collect[i];
      const auto X_tilde = X_tilde_arr[i];

      for (Index pos = 0; pos < NQuad; pos++) {
        Numeric sum = 0.0;
        for (Index j = 0; j < NQuad; j++) {
          for (Index k = 0; k < NQuad; k++) {
            sum -= G(pos, j) / (1.0 / mu0 + K[j]) * G_inv(j, k) * X_tilde[k];
          }
        }
        B[pos] = sum;
      }
    }
  }
}

void DoubleGaussLegendre(Vector& x, Vector& w, Index n) {
  const Index N = n / 2;

  Vector part_x(N);
  Vector part_w(N);
  GSL::Integration::GaussLegendre(part_x, part_w, n);

  x.resize(n);
  w.resize(n);
  for (Index i = 0; i < N; i++) {
    x[i] = -part_x[N - 1 - i];
    x[N + i] = part_x[i];
    w[i] = part_w[N - 1 - i];
    w[N + i] = part_w[i];
  }
  x += 1.0;
  x *= 0.5;
  w *= 0.5;
}

main_data::main_data(Index NQuad_,
                     Index NLeg_,
                     Index NFourier_,
                     Vector tau_arr_,
                     Vector omega_arr_,
                     Matrix Leg_coeffs_all_,
                     Matrix b_pos_,
                     Matrix b_neg_,
                     Vector f_arr_,
                     Matrix s_poly_coeffs_,
                     std::vector<BDRF> fourier_modes_,
                     Numeric mu0_,
                     Numeric I0_,
                     Numeric phi0_)
    : NLayers(tau_arr_.size()),
      NLeg_all(Leg_coeffs_all_.ncols()),
      NQuad(NQuad_),
      N(NQuad / 2),
      NFourier(NFourier_),
      Nscoeffs(s_poly_coeffs_.nrows()),
      NBDRF(fourier_modes_.size()),
      NLeg(NLeg_),
      tau_arr(std::move(tau_arr_)),
      omega_arr(std::move(omega_arr_)),
      f_arr(std::move(f_arr_)),
      Leg_coeffs_all(std::move(Leg_coeffs_all_)),
      s_poly_coeffs(std::move(s_poly_coeffs_)),
      b_pos(std::move(b_pos_)),
      b_neg(std::move(b_neg_)),
      fourier_modes(std::move(fourier_modes_)),
      mu0(mu0_),
      I0(I0_),
      phi0(phi0_),
      beam_source_bool(I0 > 0),
      iso_source_bool(Nscoeffs > 0),
      multilayer_bool(NLayers > 1) {
  // Index correctness checks
  ARTS_USER_ERROR_IF(NLeg <= 0, "NLeg must be positive");
  ARTS_USER_ERROR_IF(NLeg > NLeg_all,
                     "NLeg must be less than or equal to NLeg_all");
  ARTS_USER_ERROR_IF(NQuad < 2, "NQuad must be at least 2");
  ARTS_USER_ERROR_IF(NQuad % 2 != 0, "NQuad must be even");
  ARTS_USER_ERROR_IF(NFourier <= 0, "NFourier must be positive");
  ARTS_USER_ERROR_IF(NFourier > NLeg,
                     "NFourier must be less than or equal to NLeg");
  ARTS_USER_ERROR_IF(NQuad < NLeg,
                     "NQuad must be at least NLeg (for accuracy reasons)");

  // Input size checks
  ARTS_USER_ERROR_IF(tau_arr.size() != NLayers, "tau_arr layers mismatch");
  ARTS_USER_ERROR_IF(omega_arr.size() != NLayers, "omega_arr layers mismatch");
  ARTS_USER_ERROR_IF(not(f_arr.size() == NLayers or f_arr.size() == 0),
                     "f_arr layers mismatch");
  ARTS_USER_ERROR_IF(Leg_coeffs_all.nrows() != NLayers,
                     "Leg_coeffs_all layers mismatch");
  ARTS_USER_ERROR_IF(
      s_poly_coeffs.nrows() != NLayers and s_poly_coeffs.ncols() != 0,
      "s_poly_coeffs shape mismatch");

  // Memory allocation
  scale_tau.resize(NLayers);
  scaled_omega_arr.resize(NLayers);
  thickness_arr.resize(NLayers);
  scaled_tau_arr_with_0.resize(NLayers + 1);
  W.resize(N);
  mu_arr.resize(NQuad);
  M_inv.resize(N);
  mu_arr_pos.resize(N);
  weighted_Leg_coeffs_all.resize(NLayers, NLeg_all);
  Leg_coeffs.resize(NLayers, NLeg);
  weighted_scaled_Leg_coeffs.resize(NLayers, NLeg);

  // Thickness is the difference between tau values, the first difference is towards 0
  thickness_arr.resize(NLayers);
  std::adjacent_difference(
      tau_arr.begin(), tau_arr.end(), thickness_arr.begin());

  // Legendre coefficients selection
  weighted_Leg_coeffs_all = Leg_coeffs_all;
  for (Index i = 0; i < NLeg_all; i++) {
    weighted_Leg_coeffs_all(joker, i) *= 2 * i + 1;
  }
  Leg_coeffs = Leg_coeffs_all(joker, Range(0, NLeg));

  // Quadrature points
  DoubleGaussLegendre(mu_arr_pos, W, N);
  std::transform(
      mu_arr_pos.begin(), mu_arr_pos.end(), M_inv.begin(), [](auto&& x) {
        return 1.0 / x;
      });
  std::transform(
      mu_arr_pos.begin(), mu_arr_pos.end(), mu_arr.begin(), [](auto&& x) {
        return x;
      });
  std::transform(
      mu_arr_pos.begin(), mu_arr_pos.end(), mu_arr.begin() + N, [](auto&& x) {
        return -x;
      });

  scalar_b_pos = b_pos.size() == 1;
  if (not scalar_b_pos) {
    ARTS_USER_ERROR_IF((b_pos.shape() != std::array{N, NFourier}),
                       "Bad shape for b_pos, should be (",
                       N,
                       ", ",
                       NFourier,
                       ") is ",
                       matpack::shape_help<2>{b_pos.shape()})
  }

  scalar_b_neg = b_neg.size() == 1;
  if (not scalar_b_neg) {
    ARTS_USER_ERROR_IF((b_neg.shape() != std::array{N, NFourier}),
                       "Bad shape for b_neg, should be (",
                       N,
                       ", ",
                       NFourier,
                       ") is ",
                       matpack::shape_help<2>{b_neg.shape()})
  }

  // Origin of I0
  if ((scalar_b_pos and b_pos(0, 0) == 0) and not iso_source_bool and
      beam_source_bool) {
    I0_orig = I0;
    I0 = 1;
  } else {
    I0_orig = 1;
  }

  // Scaling
  if (not f_arr.empty()) {
    std::transform(omega_arr.begin(),
                   omega_arr.end(),
                   f_arr.begin(),
                   scaled_omega_arr.begin(),
                   [](auto&& omega, auto&& f) { return 1.0 - omega * f; });

    einsum<"i", "i", "i", "">(scale_tau, omega_arr, f_arr, -1);

    scale_tau += 1;

    scaled_tau_arr_with_0[0] = 0;

    Numeric part = 0.0;
    for (Index i = 0; i < NLayers; i++) {
      part += scale_tau[i] * thickness_arr[i];
      scaled_tau_arr_with_0[i + 1] = part;
    }

    for (Index i = 0; i < NLayers; i++) {
      for (Index j = 0; j < NLeg; j++) {
        weighted_scaled_Leg_coeffs(i, j) = static_cast<Numeric>(2 * j + 1) *
                                           (Leg_coeffs(i, j) - f_arr[i]) /
                                           (1 - f_arr[i]);
      }
    }

    for (Index i = 0; i < NLayers; i++) {
      scaled_omega_arr[i] = omega_arr[i] * (1.0 - f_arr[i]) / scale_tau[i];
    }
  } else {
    scale_tau = 1.0;
    scaled_tau_arr_with_0[0] = 0;
    scaled_tau_arr_with_0(Range(1, NLayers)) = tau_arr;
    weighted_scaled_Leg_coeffs = weighted_Leg_coeffs_all(joker, Range(0, NLeg));
    scaled_omega_arr = omega_arr;
  }

  // Value checks
  ARTS_USER_ERROR_IF(std::ranges::any_of(tau_arr, Cmp::le(0.0)),
                     "tau_arr must be positive");
  ARTS_USER_ERROR_IF(std::ranges::any_of(thickness_arr, Cmp::le(0.0)),
                     "thickness_arr must be positive");
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(omega_arr,
                          [](auto&& omega) { return omega >= 1 or omega < 0; }),
      "omega_arr must be [0, 1)");
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(Leg_coeffs_all,
                          [](auto&& x) {
                            return x[0] != 1 or
                                   std::ranges::any_of(x, [](auto&& u) {
                                     return std::abs<Numeric>(u) > 1;
                                   });
                          }),
      "Leg_coeffs_all must have 1 in the first column and be [-1, 1] elsewhere");
  ARTS_USER_ERROR_IF(I0 < 0, "I0 must be non-negative");
  ARTS_USER_ERROR_IF(mu0 < 0 or mu0 > 1, "mu0 must be [0, 1]");
  ARTS_USER_ERROR_IF(phi0 < 0 and phi0 < Constant::two_pi,
                     "phi0 must be [0, 2*pi)");
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(f_arr, [](auto&& x) { return x > 1 or x < 0; }),
      "f_arr must be [0, 1]");
  ARTS_USER_ERROR_IF(
      std::ranges::any_of(
          mu_arr_pos, [mu = mu0](auto&& x) { return std::abs(x - mu) < 1e-8; }),
      "mu0 in mu_arr_pos, this creates a singularity.  Change NQuad or mu0.");

  // IMS
  const Numeric sum1 = omega_arr * tau_arr;
  omega_avg = sum1 / sum(tau_arr);
  const Numeric sum2 = f_arr.empty() ? 0.0
                                     : einsum<Numeric, "", "i", "i", "i">(
                                           {}, f_arr, omega_arr, tau_arr);
  f_avg = sum2 / sum1;
  Leg_coeffs_residue = Leg_coeffs_all;
  if (f_arr.empty()) {
    Leg_coeffs_residue = 0;
  } else {
    for (Index i = 0; i < NLayers; i++) {
      for (Index j = 0; j < NLeg; j++) {
        Leg_coeffs_residue(i, j) = f_arr[i];
      }
    }
  }
  Leg_coeffs_residue_avg.resize(NLeg_all);
  for (Index i = 0; i < NLeg_all; i++) {
    Numeric sum3 = 0.0;
    for (Index j = 0; j < NLayers; j++) {
      sum3 += Leg_coeffs_residue(j, i) * omega_arr[j] * tau_arr[j];
    }
    Leg_coeffs_residue_avg[i] = sum3 / sum2;
    Leg_coeffs_residue_avg[i] =
        static_cast<Numeric>(2 * i + 1) *
        (2 * Leg_coeffs_residue_avg[i] - Math::pow2(Leg_coeffs_residue_avg[i]));
  }
  scaled_mu0 = mu0 / (1 - omega_avg * f_avg);

  // Diagonalize
  diagonalize(G_collect,
              K_collect,
              B_collect,
              G_inv_collect_0,
              NFourier,
              scaled_omega_arr,
              mu_arr_pos,
              mu_arr,
              M_inv,
              W,
              N,
              NQuad,
              NLeg,
              NLayers,
              weighted_scaled_Leg_coeffs,
              mu0,
              I0,
              beam_source_bool,
              iso_source_bool);

  solve_for_coefs(GC_collect,
                  NFourier,
                  G_collect,
                  K_collect,
                  B_collect,
                  G_inv_collect_0,
                  tau_arr,
                  scaled_tau_arr_with_0,
                  mu_arr_pos,
                  einsum<Vector, "i", "i", "i">({N}, mu_arr_pos, W),
                  mu_arr,
                  N,
                  NQuad,
                  NLayers,
                  NBDRF,
                  multilayer_bool,
                  fourier_modes,
                  mu0,
                  I0,
                  beam_source_bool,
                  b_pos,
                  b_neg,
                  scalar_b_pos,
                  scalar_b_neg,
                  s_poly_coeffs,
                  iso_source_bool);
}

[[nodiscard]] Index main_data::tau_index(const Numeric tau) const {
  const Index l =
      std::distance(tau_arr.begin(), std::ranges::lower_bound(tau_arr, tau));
  ARTS_USER_ERROR_IF(
      l == NLayers, "tau (", tau, ") must be at most ", tau_arr.back());
  return l;
}

void main_data::u(u_data& data,
                  const Numeric tau,
                  const Numeric phi,
                  const bool return_fourier_error) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.exponent.resize(NFourier, NQuad);
  for (Index i = 0; i < NFourier; i++) {
    for (Index j = 0; j < N; j++) {
      data.exponent(i, j) =
          std::exp(K_collect(i, l, j) * (scaled_tau - scaled_tau_arr_lm1));
      data.exponent(i, j + N) =
          std::exp(K_collect(i, l, j + N) * (scaled_tau - scaled_tau_arr_l));
    }
  }

  data.um.resize(NFourier, NQuad);
  einsum<"mi", "mij", "mj">(
      data.um, GC_collect(joker, l, joker, joker), data.exponent);
      print(GC_collect);

  if (beam_source_bool) {
    const auto tmp = B_collect(joker, l, joker);
    std::transform(tmp.elem_begin(),
                   tmp.elem_end(),
                   data.um.elem_begin(),
                   data.um.elem_begin(),
                   [scl = std::exp(-scaled_tau / mu0)](auto&& x, auto&& y) {
                     return scl * x + y;
                   });
  }

  if (iso_source_bool) {
    mathscr_v_add(data.um[0],
                  data.mathscr_v_coeffs,
                  tau,
                  s_poly_coeffs[l],
                  G_collect[0][l],
                  K_collect[0][l],
                  G_inv_collect_0[l],
                  mu_arr);
  }

  data.intensities.resize(NQuad);
  data.intensities = 0.0;
  for (Index m = 0; m < NFourier; m++) {
    const Numeric cp = std::cos(static_cast<Numeric>(m) * (phi0 - phi));
    const auto umm = data.um[m];
    for (Index i = 0; i < NQuad; i++) {
      data.intensities[i] += umm[i] * cp;
    }
  }

  if (return_fourier_error) {
    data.ulast.resize(NQuad);
    einsum<"i", "ij", "j">(data.ulast,
                           GC_collect(NFourier - 1, l, joker, joker),
                           data.exponent[NFourier - 1]);
    if (beam_source_bool) {
      const auto tmp = B_collect(NFourier - 1, l, joker);
      std::transform(tmp.elem_begin(),
                     tmp.elem_end(),
                     data.ulast.elem_begin(),
                     data.ulast.elem_begin(),
                     [scl = std::exp(-scaled_tau / mu0)](auto&& x, auto&& y) {
                       return scl * x + y;
                     });
    }

    throw std::runtime_error(
        "Not implemented yet, cannot figure out max-statement");
  }

  data.intensities *= I0_orig;
}

void main_data::u0(u0_data& data, const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.exponent.resize(NQuad);
  for (Index j = 0; j < N; j++) {
    data.exponent[j] =
        std::exp(K_collect(0, l, j) * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[j + N] =
        std::exp(K_collect(0, l, j + N) * (scaled_tau - scaled_tau_arr_l));
  }

  data.u0.resize(NQuad);
  einsum<"i", "ij", "j">(
      data.u0, GC_collect(0, l, joker, joker), data.exponent);
  if (beam_source_bool) {
    const auto tmp = B_collect(0, l, joker);
    std::transform(tmp.elem_begin(),
                   tmp.elem_end(),
                   data.u0.elem_begin(),
                   data.u0.elem_begin(),
                   [scl = std::exp(-scaled_tau / mu0)](auto&& x, auto&& y) {
                     return scl * x + y;
                   });
  }

  if (iso_source_bool) {
    mathscr_v_add(data.u0,
                  data.mathscr_v_coeffs,
                  tau,
                  s_poly_coeffs[l],
                  G_collect[0][l],
                  K_collect[0][l],
                  G_inv_collect_0[l],
                  mu_arr);
  }

  data.u0 *= I0_orig;
}

void calculate_nu(Vector& nu,
                  const ExhaustiveConstVectorView& mu,
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

void main_data::TMS(tms_data& data,
                    const Numeric tau,
                    const Numeric phi) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  // mathscr_B
  calculate_nu(data.nu, mu_arr, phi, -mu0, phi0);

  data.p_true.resize(NLayers, NQuad);
  data.p_trun.resize(NLayers, NQuad);

  for (Index j = 0; j < NLayers; j++) {
    for (Index i = 0; i < NQuad; i++) {
      data.p_true(j, i) =
          Legendre::legendre_sum(weighted_Leg_coeffs_all[j], data.nu[i]);
      data.p_trun(j, i) =
          Legendre::legendre_sum(weighted_scaled_Leg_coeffs[j], data.nu[i]);
    }
  }

  data.mathscr_B.resize(NLayers, NQuad);
  if (f_arr.empty()) {
    for (Index j = 0; j < NLayers; j++) {
      for (Index i = 0; i < NQuad; i++) {
        data.mathscr_B(j, i) = (scaled_omega_arr[j] * I0) / (4 * Constant::pi) *
                               (mu0 / (mu0 + mu_arr[i])) *
                               (data.p_true(j, i) - data.p_trun(j, i));
      }
    }
  } else {
    for (Index j = 0; j < NLayers; j++) {
      for (Index i = 0; i < NQuad; i++) {
        data.mathscr_B(j, i) =
            (scaled_omega_arr[j] * I0) / (4 * Constant::pi) *
            (mu0 / (mu0 + mu_arr[i])) *
            (data.p_true(j, i) / (1 - f_arr[j]) - data.p_trun(j, i));
      }
    }
  }

  data.TMS_correction_pos.resize(N);
  data.TMS_correction_pos = std::exp(-scaled_tau / mu0);
  for (Index i = 0; i < N; i++) {
    data.TMS_correction_pos[i] -=
        std::exp((scaled_tau - scaled_tau_arr_l) / mu_arr_pos[i] -
                 scaled_tau_arr_l / mu0);
  }

  data.TMS_correction_neg.resize(N);
  data.TMS_correction_neg = std::exp(-scaled_tau / mu0);
  for (Index i = 0; i < N; i++) {
    data.TMS_correction_neg[i] -=
        std::exp((scaled_tau_arr_lm1 - scaled_tau) / mu_arr_pos[i] -
                 scaled_tau_arr_lm1 / mu0);
  }

  data.TMS.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.TMS[i + 0] = data.mathscr_B(l, i + 0) * data.TMS_correction_pos[i];
    data.TMS[i + N] = data.mathscr_B(l, i + N) * data.TMS_correction_neg[i];
  }

  if (multilayer_bool) {
    data.contribution_from_other_layers_pos.resize(N, NLayers);
    data.contribution_from_other_layers_pos = 0;
    data.contribution_from_other_layers_neg.resize(N, NLayers);
    data.contribution_from_other_layers_neg = 0;
    for (Index i = 0; i < NLayers; i++) {
      if (l > i) {
        // neg
        for (Index j = 0; j < N; j++) {
          data.contribution_from_other_layers_neg(j, i) =
              (data.mathscr_B(i, j + N) *
               (std::exp((scaled_tau_arr_with_0[i] - scaled_tau) /
                             mu_arr_pos[j] -
                         scaled_tau_arr_with_0[i] / mu0) -
                std::exp((scaled_tau_arr_with_0[i] - scaled_tau) /
                             mu_arr_pos[j] -
                         scaled_tau_arr_with_0[i] / mu0)));
        }
      } else if (l < i) {
        // pos
        for (Index j = 0; j < N; j++) {
          data.contribution_from_other_layers_pos(j, i) =
              (data.mathscr_B(i, j) *
               (std::exp((scaled_tau - scaled_tau_arr_with_0[i]) /
                             mu_arr_pos[j] -
                         scaled_tau_arr_with_0[i] / mu0) -
                std::exp((scaled_tau - scaled_tau_arr_with_0[i]) /
                             mu_arr_pos[j] -
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

void main_data::IMS(ims_data& data,
                    const Numeric tau,
                    const Numeric phi) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");

  // IMS nu is just for mu_arr_neg, which is just the second half of mu_arr
  calculate_nu(data.nu, mu_arr.slice(N, N), phi, -mu0, phi0);

  data.x.resize(N);
  for (Index i = 0; i < N; i++) {
    data.x[i] = 1.0 / mu_arr_pos[i] - 1.0 / scaled_mu0;
  }

  data.chi.resize(N);
  for (Index i = 0; i < N; i++) {
    data.chi[i] = (1 / (mu_arr_pos[i] * scaled_mu0 * data.x[i])) *
                  ((tau - 1 / data.x[i]) * std::exp(-tau / scaled_mu0) +
                   std::exp(-tau / mu_arr_pos[i]) / data.x[i]);
  }

  data.IMS.resize(N);
  for (Index i = 0; i < N; i++) {
    data.IMS[i] = (I0 / (4 * Constant::pi) * Math::pow2(omega_avg * f_avg) /
                   (1 - omega_avg * f_avg) *
                   Legendre::legendre_sum(Leg_coeffs_residue_avg, data.nu[i])) *
                  data.chi[i];
  }
}

void main_data::u_corr(u_data& u_data,
                       ims_data& ims_data,
                       tms_data& tms_data,
                       const Numeric tau,
                       const Numeric phi,
                       const bool return_fourier_error) const {
  TMS(tms_data, tau, phi);
  IMS(ims_data, tau, phi);
  u(u_data, tau, phi, return_fourier_error);

  for (Index i = 0; i < N; i++) {
    u_data.intensities[i] += I0_orig * tms_data.TMS[i];
  }

  for (Index i = N; i < NQuad; i++) {
    u_data.intensities[i] += I0_orig * (tms_data.TMS[i] + ims_data.IMS[i - N]);
  }

  if (return_fourier_error) {
    throw std::runtime_error("Not implemented yet");
  }
}

Numeric main_data::flux_up(flux_data& data, const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");
  ARTS_USER_ERROR_IF(tau > tau_arr.back(),
                     "tau (",
                     tau,
                     ") must be less than the last layer (",
                     tau_arr.back(),
                     ")");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  data.mathscr_v.resize(NQuad);
  data.mathscr_v = 0.0;
  if (iso_source_bool) {
    mathscr_v_add(data.mathscr_v,
                  data.mathscr_v_coeffs,
                  tau,
                  s_poly_coeffs[l],
                  G_collect[0][l],
                  K_collect[0][l],
                  G_inv_collect_0[l],
                  mu_arr);
  }

  data.direct_beam_contribution.resize(N);
  if (beam_source_bool) {
    einsum<"i", "i", "">(data.direct_beam_contribution,
                         B_collect(0, l, joker).slice(0, N),
                         std::exp(-scaled_tau / mu0));
  } else {
    data.direct_beam_contribution = 0.0;
  }

  data.exponent.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.exponent[i] =
        std::exp(K_collect(0, l, i) * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[i + N] =
        std::exp(K_collect(0, l, i + N) * (scaled_tau - scaled_tau_arr_l));
  }

  data.u0_pos.resize(N);
  einsum<"i", "ij", "j">(
      data.u0_pos, GC_collect(0, l, Range(0, N), joker), data.exponent);
  data.u0_pos += data.mathscr_v.slice(0, N);
  data.u0_pos += data.direct_beam_contribution;

  return einsum<Numeric, "", "", "i", "i", "i">(
      {}, Constant::two_pi * I0_orig, mu_arr_pos, W, data.u0_pos);
}

std::pair<Numeric, Numeric> main_data::flux_down(flux_data& data,
                                                 const Numeric tau) const {
  ARTS_USER_ERROR_IF(tau < 0, "tau (", tau, ") must be positive");
  ARTS_USER_ERROR_IF(tau > tau_arr.back(),
                     "tau (",
                     tau,
                     ") must be less than the last layer (",
                     tau_arr.back(),
                     ")");

  const Index l = tau_index(tau);

  const Numeric scaled_tau_arr_l = scaled_tau_arr_with_0[l + 1];
  const Numeric scaled_tau_arr_lm1 = scaled_tau_arr_with_0[l];
  const Numeric scaled_tau =
      scaled_tau_arr_l - (tau_arr[l] - tau) * scale_tau[l];

  if (iso_source_bool) {
    data.mathscr_v.resize(NQuad);
    data.mathscr_v = 0.0;
    mathscr_v_add(data.mathscr_v,
                  data.mathscr_v_coeffs,
                  tau,
                  s_poly_coeffs[l],
                  G_collect[0][l],
                  K_collect[0][l],
                  G_inv_collect_0[l],
                  mu_arr);
    data.mathscr_v.slice(0, N) = data.mathscr_v.slice(N, N);
    data.mathscr_v.resize(N);
  } else {
    data.mathscr_v.resize(N);
    data.mathscr_v = 0.0;
  }

  const Numeric direct_beam =
      beam_source_bool ? I0 * mu0 * std::exp(-tau / mu0) : 0;
  const Numeric direct_beam_scaled =
      beam_source_bool ? I0 * mu0 * std::exp(-scaled_tau / mu0) : 0;
  data.direct_beam_contribution.resize(N);
  if (beam_source_bool) {
    einsum<"i", "i", "">(data.direct_beam_contribution,
                         B_collect(0, l, joker).slice(N, N),
                         std::exp(-scaled_tau / mu0));
  } else {
    data.direct_beam_contribution = 0.0;
  }

  data.exponent.resize(NQuad);
  for (Index i = 0; i < N; i++) {
    data.exponent[i] =
        std::exp(K_collect(0, l, i) * (scaled_tau - scaled_tau_arr_lm1));
    data.exponent[i + N] =
        std::exp(K_collect(0, l, i + N) * (scaled_tau - scaled_tau_arr_l));
  }

  data.u0_neg.resize(N);
  einsum<"i", "ij", "j">(
      data.u0_neg, GC_collect(0, l, Range(N, N), joker), data.exponent);
  data.u0_neg += data.mathscr_v.slice(0, N);
  data.u0_neg += data.direct_beam_contribution;

  return {I0_orig * (einsum<Numeric, "", "", "i", "i", "i">(
                         {}, Constant::two_pi, mu_arr_pos, W, data.u0_neg) -
                     direct_beam + direct_beam_scaled),
          I0_orig * I0 * direct_beam};
}
}  // namespace disort
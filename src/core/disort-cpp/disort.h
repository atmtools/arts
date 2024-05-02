#pragma once

#include <matpack.h>

#include <functional>

#include "lin_alg.h"
#include "matpack_view.h"

namespace disort {
struct BDRF {
  std::function<void(ExhaustiveMatrixView,
                     const ExhaustiveConstVectorView&,
                     const ExhaustiveConstVectorView&)>
      f;
  void operator()(ExhaustiveMatrixView x,
                  const ExhaustiveConstVectorView& a,
                  const ExhaustiveConstVectorView& b) const {
    ARTS_ASSERT(x.nrows() == a.size());
    ARTS_ASSERT(x.ncols() == b.size());
    f(x, a, b);
  }
};

struct mathscr_v_data {
  Matrix src;
  Matrix G;
  solve_workdata solve_work;
  Vector k1, k2;
  Vector cvec;
  mathscr_v_data(const Index Nk = 0, const Index Nc = 0)
      : src(Nk, Nc), G(Nk, Nk), solve_work(Nk), k1(Nk), k2(Nk), cvec(Nc) {}
  void resize(const Index Nk, const Index Nc) {
    src.resize(Nk, Nc);
    G.resize(Nk, Nk);
    solve_work.resize(Nk);
    k1.resize(Nk);
    k2.resize(Nk);
    cvec.resize(Nc);
  }
};

struct u_data {
  Matrix exponent;
  Matrix um;
  mathscr_v_data src;
  Vector intensities;
  Vector ulast;
};

struct u0_data {
  mathscr_v_data src;
  Vector exponent;
  Vector u0;
};

struct tms_data {
  Vector nu;
  Vector TMS;
  Matrix mathscr_B;
  Matrix contribution_from_other_layers_pos;
  Matrix contribution_from_other_layers_neg;
};

struct flux_data {
  Vector mathscr_v;
  Vector exponent;
  Vector direct_beam_contribution;
  mathscr_v_data src;
  Vector u0_pos;
  Vector u0_neg;
};

class main_data {
  const Index NLayers;
  const Index NQuad;
  const Index NLeg;
  const Index NFourier;
  const Index N;
  const Index Nscoeffs;
  const Index NLeg_all;
  const Index NBDRF;

  const bool has_beam_source{};
  const bool has_source_poly{};
  const bool is_multilayer{};

  //! User inputs
  AscendingGrid tau_arr{};                 // [NLayers]
  Vector omega_arr{};                      // [NLayers]
  Vector f_arr{};                          // [NLayers] or [0]
  Matrix source_poly_coeffs{};             // [NLayers, Nscoeffs] or [0, 0]
  Matrix Leg_coeffs_all{};                 // [NLayers, NLeg_all]
  Matrix b_pos{};                          // [1, 1] or [N, NFourier]
  Matrix b_neg{};                          // [1, 1] or [N, NFourier]
  std::vector<BDRF> brdf_fourier_modes{};  // [NBDRF]
  Numeric mu0{};
  Numeric I0{};
  Numeric phi0{};

  //! Derived values
  Vector scale_tau{};                   // [NLayers]
  Vector scaled_omega_arr{};            // [NLayers]
  Vector scaled_tau_arr_with_0{};       // [NLayers + 1]
  Vector mu_arr{};                      // [NQuad]
  Vector inv_mu_arr{};                  // [NQuad]
  Vector W{};                           // [N]
  Vector Leg_coeffs_residue_avg{};      // [NLeg_all]
  Matrix weighted_scaled_Leg_coeffs{};  // [NLayers, NLeg]
  Matrix weighted_Leg_coeffs_all{};     // [NLayers, NLeg_all]
  Tensor4 GC_collect{};                 // [NFourier, NLayers, NQuad, NQuad]
  Tensor4 G_collect{};                  // [NFourier, NLayers, NQuad, NQuad]
  Tensor3 K_collect{};                  // [NFourier, NLayers, NQuad]
  Tensor3 B_collect{};                  // [NFourier, NLayers, NQuad]
  Numeric I0_orig{};
  Numeric f_avg{};
  Numeric omega_avg{};
  Numeric scaled_mu0{};

  //! Internal compute data
  const Index n;                        // NQuad * NLayers;
  Vector RHS{};                         // [n]
  Vector jvec{};                        // [NQuad]
  Vector fac{};                         // [NLeg]
  Vector weighted_asso_Leg_coeffs_l{};  // [NLeg]
  Vector asso_leg_term_mu0{};           // [NLeg]
  Vector X_temp{};                      // [NLeg]
  solve_workdata solve_work{};          // [NQuad]
  diagonalize_workdata diag_work{};     // [N]
  Vector mathscr_X_pos{};               // [N]
  Vector b_pos_m{};                     // [N]
  Vector b_neg_m{};                     // [N]
  Vector E_Lm1L{};                      // [N]
  Vector E_lm1l{};                      // [N]
  Vector E_llp1{};                      // [N]
  Vector BDRF_RHS_contribution{};       // [N]
  matpack::band_matrix LHSB;            // [n, n]
  mathscr_v_data
      comp_data;  // [NQuad, Nscoeffs] + [NQuad, NQuad] + 3 * [Nquad] + [Nscoeffs]
  Matrix Gml{};                        // [NQuad, NQuad]
  Matrix BDRF_LHS_contribution_neg{};  // [N, N]
  Matrix BDRF_LHS_contribution_pos{};  // [N, N]
  Matrix R{};                          // [N, N]
  Matrix mathscr_D_neg{};              // [N, N]
  Matrix D_pos{};                      // [N, N]
  Matrix D_neg{};                      // [N, N]
  Matrix apb{};                        // [N, N]
  Matrix amb{};                        // [N, N]
  Matrix sqr{};                        // [N, N]
  Matrix asso_leg_term_pos{};          // [N, NLeg]
  Matrix asso_leg_term_neg{};          // [N, NLeg]
  Matrix D_temp{};                     // [N, NLeg]

 public:
  main_data(const Index NQuad,
            const Index NLeg,
            const Index NFourier,
            AscendingGrid tau_arr,
            Vector omega_arr,
            Matrix Leg_coeffs_all,
            Matrix b_pos,
            Matrix b_neg,
            Vector f_arr,
            Matrix source_poly_coeffs,
            std::vector<BDRF> brdf_fourier_modes,
            Numeric mu0,
            Numeric I0,
            Numeric phi0);

  [[nodiscard]] Index quads() const { return mu_arr.size(); }
  [[nodiscard]] Index fouriers() const { return GC_collect.nbooks(); }
  [[nodiscard]] Index layers() const { return tau_arr.size(); }
  [[nodiscard]] Index legalls() const {
    return weighted_Leg_coeffs_all.ncols();
  }
  [[nodiscard]] Index scoeffs() const { return source_poly_coeffs.nrows(); }

  [[nodiscard]] Index tau_index(const Numeric tau) const;

  void TMS(tms_data& data, const Numeric tau, const Numeric phi) const;

  void IMS(Vector& ims, const Numeric tau, const Numeric phi) const;

  void u(u_data& data,
         const Numeric tau,
         const Numeric phi,
         const bool return_fourier_error = false) const;

  void u0(u0_data& data, const Numeric tau) const;

  void u_corr(u_data& u_data,
              Vector& ims,
              tms_data& tms_data,
              const Numeric tau,
              const Numeric phi,
              const bool return_fourier_error = false) const;

  [[nodiscard]] Numeric flux_up(flux_data&, const Numeric tau) const;
  [[nodiscard]] std::pair<Numeric, Numeric> flux_down(flux_data&,
                                                      const Numeric tau) const;

  void set_ims_factors();
  void set_scales();
  void diagonalize();
  void solve_for_coefs();
  void set_weighted_Leg_coeffs_all();
  void set_beam(const Numeric I0);
  void update_all(const Numeric I0);

  void check_input_size() const;
  void check_input_value() const;

  ExhaustiveConstVectorView mu() const { return mu_arr; }
  ExhaustiveConstVectorView tau() const { return tau_arr; }
};
}  // namespace disort

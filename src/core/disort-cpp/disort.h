#pragma once

#include <matpack.h>

#include <functional>

#include "matpack_view.h"

namespace disort {
void mathscr_v_add(ExhaustiveVectorView um,
                   Matrix& mathscr_v_coeffs,
                   const Numeric tau,
                   const ExhaustiveConstVectorView& s_poly_coeffs,
                   const ExhaustiveConstMatrixView& G,
                   const ExhaustiveConstVectorView& K,
                   const ExhaustiveConstMatrixView& G_inv,
                   const ExhaustiveConstVectorView& mu_arr);

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

struct u_data {
  Matrix exponent;
  Matrix um;
  Matrix mathscr_v_coeffs;
  Vector intensities;
  Vector ulast;
};

struct u0_data {
  Matrix mathscr_v_compdata;
  Vector exponent;
  Vector u0;
};

struct tms_data {
  Vector nu;
  Vector TMS_correction_pos;
  Vector TMS_correction_neg;
  Vector TMS;
  Matrix mathscr_B;
  Matrix contribution_from_other_layers_pos;
  Matrix contribution_from_other_layers_neg;
};

struct flux_data {
  Vector mathscr_v;
  Vector exponent;
  Vector direct_beam_contribution;
  Matrix mathscr_v_compdata;
  Vector u0_pos;
  Vector u0_neg;
};

class main_data {
  AscendingGrid tau_arr{};      // [NLayers]
  Vector omega_arr{};           // [NLayers]
  Vector f_arr{};               // [NLayers] or [0]
  Matrix Leg_coeffs_all{};      // [NLayers, NLeg_all]
  Matrix source_poly_coeffs{};  // [Nscoeffs, NLayers] or [0, 0]
  Matrix b_pos{};               // [NQuad/2, NFourier] or [1, 1]
  Matrix b_neg{};               // [NQuad/2, NFourier] or [1, 1]
  Numeric mu0{};
  Numeric I0{};
  Numeric phi0{};

  bool has_beam_source{};

  // When source function exists
  Tensor3 G_inv_collect_0{};  // [NLayers, NQuad, NQuad]

  // For all solvers
  Tensor4 GC_collect{};            // [NFourier, NLayers, NQuad, NQuad]
  Tensor4 G_collect{};             // [NFourier, NLayers, NQuad, NQuad]
  Tensor3 K_collect{};             // [NFourier, NLayers, NQuad]
  Tensor3 B_collect{};             // [NFourier, NLayers, NQuad]
  Vector scale_tau{};              // [NLayers]
  Vector scaled_tau_arr_with_0{};  // [NLayers + 1]
  Vector mu_arr{};                 // [NQuad]
  Numeric I0_orig{};

  // For flux
  Vector W{};  // [NQuad/2]

  // For TMS
  Vector scaled_omega_arr{};            // [NLayers]
  Matrix weighted_scaled_Leg_coeffs{};  // [NLayers, NLeg]
  Matrix weighted_Leg_coeffs_all{};     // [NLayers, NLeg_all]

  // For IMS
  Vector Leg_coeffs_residue_avg{};  // [NLeg_all]
  Numeric f_avg{};
  Numeric omega_avg{};
  Numeric scaled_mu0{};

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
            const std::vector<BDRF>& brdf_fourier_modes,
            Numeric mu0,
            Numeric I0,
            Numeric phi0);

  [[nodiscard]] Index quads() const { return mu_arr.size(); }
  [[nodiscard]] Index fouriers() const { return GC_collect.nbooks(); }
  [[nodiscard]] Index layers() const { return tau_arr.size(); }
  [[nodiscard]] Index legalls() const { return Leg_coeffs_all.ncols(); }
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
};
}  // namespace disort

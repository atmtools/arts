#pragma once

#include <matpack.h>

#include "matpack_view.h"

namespace disort {
struct BDRF {
  void operator()(ExhaustiveMatrixView,
                                   const ExhaustiveConstVectorView&,
                                   const ExhaustiveConstVectorView&) const {

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
  Matrix mathscr_v_coeffs;
  Vector exponent;
  Vector u0;
};

struct tms_data {
  Vector nu;
  Vector TMS_correction_pos;
  Vector TMS_correction_neg;
  Vector TMS;
  Matrix p_true;
  Matrix p_trun;
  Matrix mathscr_B;
  Matrix contribution_from_other_layers_pos;
  Matrix contribution_from_other_layers_neg;
};

struct ims_data {
  Vector nu;
  Vector x;
  Vector chi;
  Vector IMS;
};

struct flux_data {
  Vector mathscr_v;
  Vector exponent;
  Vector direct_beam_contribution;
  Matrix mathscr_v_coeffs;
  Vector u0_pos;
  Vector u0_neg;
};

class main_data {
  Index NLayers;
  Index NLeg_all;  // Propbably not needed
  Index NQuad;
  Index N;
  Index NFourier;
  Index Nscoeffs;
  Index NBDRF;
  Index NLeg;

  Vector tau_arr;                   // [NLayers]
  Vector omega_arr;                 // [NLayers]
  Vector f_arr;                     // [NLayers] or [0]
  Matrix Leg_coeffs_all;            // [NLayers, NLeg_all]
  Matrix s_poly_coeffs;             // [Nscoeffs, NLayers] or [0, 0]
  Matrix b_pos;                     // [NQuad/2, NFourier] or [1, 1]
  Matrix b_neg;                     // [NQuad/2, NFourier] or [1, 1]
  std::vector<BDRF> fourier_modes;  // [NBDRF]
  Numeric mu0;
  Numeric I0;
  Numeric phi0;

  bool beam_source_bool;
  bool iso_source_bool;
  bool multilayer_bool;
  bool scalar_b_pos;
  bool scalar_b_neg;

  Vector scale_tau;                   // [NLayers]
  Vector scaled_omega_arr;            // [NLayers]
  Vector thickness_arr;               // [NLayers]
  Vector scaled_tau_arr_with_0;       // [NLayers + 1]
  Vector W;                           // [NQuad/2]
  Vector mu_arr;                      // [NQuad]
  Vector M_inv;                       // [NQuad/2]
  Vector mu_arr_pos;                  // [NQuad/2]
  Matrix weighted_Leg_coeffs_all;     // [NLayers, NLeg_all]
  Matrix Leg_coeffs;                  // [NLayers, NLeg]
  Matrix weighted_scaled_Leg_coeffs;  // [NLayers, NLeg]

  Tensor3 K_collect;        // [NFourier, NLayers, NQuad]
  Tensor4 G_collect;        // [NFourier, NLayers, NQuad, NQuad]
  Tensor3 G_inv_collect_0;  // [NLayers, NQuad, NQuad]
  Tensor3 B_collect;        // [NFourier, NLayers, NQuad]
  Tensor4 GC_collect;       // [NFourier, NLayers, NQuad, NQuad]

  Numeric I0_orig;
  Numeric scaled_mu0;
  Numeric omega_avg;
  Numeric f_avg;
  Matrix Leg_coeffs_residue;      // [NLayers, NLeg_all]
  Vector Leg_coeffs_residue_avg;  // [NLeg_all]

 public:
  main_data(Index NQuad,
            Vector tau_arr,
            Vector omega_arr,
            Matrix Leg_coeffs_all,
            Matrix b_pos,
            Matrix b_neg,
            Vector f_arr,
            Matrix s_poly_coeffs,
            std::vector<BDRF> fourier_modes,
            Numeric mu0,
            Numeric I0,
            Numeric phi0);

  [[nodiscard]] Size mem() const {
    return sizeof(Numeric) *
               (NLayers * 7 + 7 + NQuad * 3 + NQuad / 2 +
                NLayers * NLeg_all * 2 + NQuad / 2 * NFourier * 2 +
                Nscoeffs * NLayers + NLayers * NLeg * 2) +
           sizeof(BDRF) * NBDRF * fourier_modes.size() + sizeof(main_data);
  }

  [[nodiscard]] Index tau_index(const Numeric tau) const;

  void TMS(tms_data& data, const Numeric tau, const Numeric phi) const;

  void IMS(ims_data& data, const Numeric tau, const Numeric phi) const;

  void u(u_data& data,
         const Numeric tau,
         const Numeric phi,
         const bool return_fourier_error) const;

  void u0(u0_data& data, const Numeric tau) const;

  void u_corr(u_data& u_data,
              ims_data& ims_data,
              tms_data& tms_data,
              const Numeric tau,
              const Numeric phi,
              const bool return_fourier_error) const;

  [[nodiscard]] Numeric flux_up(flux_data&, const Numeric tau) const;
  [[nodiscard]] std::pair<Numeric, Numeric> flux_down(flux_data&, const Numeric tau) const;
};
}  // namespace disort

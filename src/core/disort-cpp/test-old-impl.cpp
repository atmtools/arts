#include <disort.h>

#include <iostream>
#include <ranges>

#include "artstime.h"
#include "cdisort/cdisort.h"
#include "configtypes.h"
#include "sorted_grid.h"

constexpr Index NQuad = 40;
constexpr Numeric mu0 = 1.0;
constexpr Numeric phi0 = 0.0;
constexpr Numeric I0 = 1.0;

static const Vector dtau{
    9.23190158e-03, 7.99612323e-03, 6.88202828e-03, 5.95522020e-03,
    5.15497874e-03, 4.48289869e-03, 3.87281039e-03, 3.33206516e-03,
    2.91754647e-03, 2.50311436e-03, 2.18882921e-03, 1.85688485e-03,
    1.68289362e-03, 1.38491055e-03, 1.19477815e-03, 1.06759017e-03,
    9.31668444e-04, 7.96194092e-04, 6.87129066e-04, 6.00884704e-04,
    5.15784517e-04, 4.46068334e-04, 3.89421522e-04, 3.38275212e-04,
    2.90595816e-04, 2.52112344e-04, 2.19027249e-04, 1.87839604e-04,
    1.64492888e-04, 1.41520881e-04, 1.22438643e-04, 1.06831714e-04,
    9.16487997e-05, 8.01885000e-05, 6.90564204e-05, 5.98078820e-05,
    5.20330908e-05, 4.46637040e-05, 3.92476383e-05, 3.36617904e-05,
    2.93364306e-05, 2.53316710e-05, 2.17883030e-05, 1.87557628e-05,
    1.61806567e-05, 1.42263200e-05, 1.22168463e-05, 1.08360138e-05,
    9.37360948e-06, 7.96154891e-06, 6.90798933e-06, 6.01614287e-06,
    5.23961682e-06, 4.49903445e-06, 3.87102132e-06, 3.38136787e-06,
    2.95388431e-06, 2.57694312e-06, 2.14562856e-06, 1.88488693e-06,
    1.65618497e-06, 1.45555265e-06, 1.25028031e-06, 1.04045012e-06,
    9.13966746e-07, 8.03029789e-07, 7.05712474e-07, 6.09834839e-07,
    4.99424639e-07, 4.38625362e-07, 3.85307215e-07, 3.38542192e-07,
    2.97518153e-07, 2.50322767e-07, 2.14207336e-07, 1.88162824e-07,
    1.65319907e-07, 1.45281750e-07, 1.27701031e-07, 1.07430593e-07};
static const Index NLayers = dtau.size();

const Vector omega(NLayers, 1.0 - 1e-6);  // not 1.0

void oldimpl(bool print_results) {
  DebugTime oldimpl{"oldimpl"};

  disort_state ds;
  disort_output out;

  ds.accur = 0.005;
  ds.flag.prnt[0] = FALSE;
  ds.flag.prnt[1] = FALSE;
  ds.flag.prnt[2] = FALSE;
  ds.flag.prnt[3] = FALSE;
  ds.flag.prnt[4] = TRUE;
  ds.flag.usrtau = FALSE;
  ds.flag.usrang = TRUE;
  ds.flag.spher = FALSE;
  ds.flag.general_source = FALSE;
  ds.flag.output_uum = FALSE;
  ds.flag.brdf_type = BRDF_NONE;
  ds.flag.ibcnd = GENERAL_BC;
  ds.flag.usrang = TRUE;
  ds.flag.planck = FALSE;
  ds.flag.onlyfl = FALSE;
  ds.flag.lamber = TRUE;
  ds.flag.quiet = FALSE;
  ds.flag.intensity_correction = FALSE;
  ds.flag.old_intensity_correction = FALSE;

  ds.nlyr = static_cast<int>(NLayers);
  ds.nstr = static_cast<int>(NQuad);
  ds.nphase = ds.nstr;
  ds.nmom = ds.nstr;
  ds.numu = static_cast<int>(1);
  ds.nphi = static_cast<int>(1);

  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds, &out);

  for (Index i = 0; i < NLayers; i++) {
    ds.dtauc[i] = dtau[i];
    ds.ssalb[i] = omega[i];
    ds.pmom[0 + i * (ds.nmom_nstr + 1)] = 1.0;
    ds.pmom[2 + i * (ds.nmom_nstr + 1)] = 0.1;
  }

  ds.umu[0] = 0.5;
  ds.phi[0] = 0.0;
  ds.bc.umu0 = mu0;
  ds.bc.phi0 = phi0;
  ds.bc.fbeam = I0;

  c_disort(&ds, &out);

  if (!print_results) {
    return;
  }

  std::cout << "cdisort out:\n";
  for (Index i = 0; i < NLayers; i++) {
    std::cout << out.rad[i].rfldn << ',' << out.rad[i].rfldir << ','
              << out.rad[i].flup << ',' << '\n';
  }
}

void newimpl(bool print_results) {
  DebugTime newimpl{"newimpl"};

  const auto tau = []() {
    Vector out(NLayers);
    out[0] = dtau[0];
    for (Index i = 1; i < out.size(); i++) {
      out[i] = out[i - 1] + dtau[i];
    }
    return AscendingGrid{out};
  }();

  const Matrix Leg_coeffs_all = []() {
    Matrix out(80, 3);
    for (Index i = 0; i < NLayers; i++) {
      out[i] = {1.0, 0.0, 0.1};
    }
    return out;
  }();

  const Matrix b_pos(1, 1, 0);
  const Matrix b_neg(1, 1, 0);
  const Vector f_arr{};
  const std::vector<disort::BDRF> BDRF_Fourier_modes{};
  const Matrix s_poly_coeffs(0, 0);

  // DebugTime setup{"setup"};
  disort::main_data dis(NQuad,
                        3,
                        3,
                        tau,
                        omega,
                        Leg_coeffs_all,
                        b_pos,
                        b_neg,
                        f_arr,
                        s_poly_coeffs,
                        BDRF_Fourier_modes,
                        mu0,
                        I0,
                        phi0);

  if (print_results) {
    std::cout << "cpp-disort out:\n";
    // DebugTime comput{"comput"};
    disort::flux_data fd;
    disort::u_data ud;
    const auto [a0, b0] = dis.flux_down(fd, 0);
    const auto c0 = dis.flux_up(fd, 0);
    std::cout << a0 << ',' << b0 << ',' << c0 << ',' << '\n';
    for (auto t : tau | std::ranges::views::take(NLayers - 1)) {
      const auto [a, b] = dis.flux_down(fd, t);
      const auto c = dis.flux_up(fd, t);
      std::cout << a << ',' << b << ',' << c << ',' << '\n';
    }
  }
}

int main() {
  bool print = false;
  oldimpl(print);
  newimpl(print);
}

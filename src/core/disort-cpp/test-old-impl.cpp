#include <disort.h>

#include <iostream>
#include <ranges>

#include "artstime.h"
#include "cdisort/cdisort.h"
#include "configtypes.h"
#include "rng.h"
#include "sorted_grid.h"

struct Timing {
  std::string_view name;
  Timing(const char* c) : name(c) {}
  TimeStep dt{};
  template <typename Function>
  void operator()(Function&& f) {
    Time start{};
    f();
    Time end{};
    dt = end - start;
  }
};

inline std::ostream& operator<<(std::ostream& os, const Array<Timing>& vt) {
  for (auto& t : vt) {
    if (t.name.contains('\n') or t.name.contains(' ') or t.name.empty())
      throw std::runtime_error(var_string("bad name: ", '"', t.name, '"'));
    if (t.name not_eq "dummy") {
      os << std::setprecision(15) << t.name << " " << t.dt.count() << '\n';
    }
  }
  return os;
}

constexpr Index NQuad    = 40;
constexpr Index NFourier = 3;
constexpr Numeric mu0    = 1.0;
constexpr Numeric phi0   = 0.0;
constexpr Numeric I0     = 1.0;

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

void oldimpl(bool print_results = false) {
  disort_state ds;
  disort_output out;

  ds.accur                         = 0.005;
  ds.flag.prnt[0]                  = FALSE;
  ds.flag.prnt[1]                  = FALSE;
  ds.flag.prnt[2]                  = FALSE;
  ds.flag.prnt[3]                  = FALSE;
  ds.flag.prnt[4]                  = TRUE;
  ds.flag.usrtau                   = FALSE;
  ds.flag.usrang                   = TRUE;
  ds.flag.spher                    = FALSE;
  ds.flag.general_source           = FALSE;
  ds.flag.output_uum               = FALSE;
  ds.flag.brdf_type                = BRDF_NONE;
  ds.flag.ibcnd                    = GENERAL_BC;
  ds.flag.usrang                   = TRUE;
  ds.flag.planck                   = FALSE;
  ds.flag.onlyfl                   = FALSE;
  ds.flag.lamber                   = TRUE;
  ds.flag.quiet                    = FALSE;
  ds.flag.intensity_correction     = FALSE;
  ds.flag.old_intensity_correction = FALSE;

  ds.nlyr   = static_cast<int>(NLayers);
  ds.nstr   = static_cast<int>(NQuad);
  ds.nphase = ds.nstr;
  ds.nmom   = ds.nstr;
  ds.numu   = static_cast<int>(1);
  ds.nphi   = static_cast<int>(1);

  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds, &out);

  for (Index i = 0; i < NLayers; i++) {
    ds.dtauc[i]                         = dtau[i];
    ds.ssalb[i]                         = omega[i];
    ds.pmom[0 + i * (ds.nmom_nstr + 1)] = 1.0;
    ds.pmom[2 + i * (ds.nmom_nstr + 1)] = 0.1;
  }

  ds.umu[0]   = 0.5;
  ds.phi[0]   = 0.0;
  ds.bc.umu0  = mu0;
  ds.bc.phi0  = phi0;
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

void newimpl(bool print_results = false) {
  const auto tau = []() {
    Vector out(NLayers);
    out[0] = dtau[0];
    for (Index i = 1; i < out.size(); i++) {
      out[i] = out[i - 1] + dtau[i];
    }
    return AscendingGrid{out};
  }();

  const Matrix Leg_coeffs_all = []() {
    Matrix out(NLayers, NFourier);
    for (Index i = 0; i < NLayers; i++) {
      out[i] = {1.0, 0.0, 0.1};
    }
    return out;
  }();

  const Matrix b_pos(NFourier, NQuad / 2, 0);
  const Matrix b_neg(NFourier, NQuad / 2, 0);
  const Vector f_arr(NLayers, 0);
  const std::vector<disort::BDRF> BDRF_Fourier_modes{};
  const Matrix s_poly_coeffs(NLayers, 0);

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
    const auto c0       = dis.flux_up(fd, 0);
    std::cout << a0 << ',' << b0 << ',' << c0 << ',' << '\n';
    for (auto t : tau | std::ranges::views::take(NLayers - 1)) {
      const auto [a, b] = dis.flux_down(fd, t);
      const auto c      = dis.flux_up(fd, t);
      std::cout << a << ',' << b << ',' << c << ',' << '\n';
    }
  }
}

std::pair<Numeric, Numeric> absrel(ExhaustiveVectorView v1,
                                   const ExhaustiveConstVectorView& v2) {
  for (double i : v1) {
    if (i == 0) std::cerr << "ERROR\n";
  }

  v1               -= v2;
  const Numeric a_  = std::abs(std::ranges::max(
      v1, [](auto a, auto b) { return std::abs(a) < std::abs(b); }));
  v1               /= v2;
  const Numeric b_  = std::abs(std::ranges::max(
      v1, [](auto a, auto b) { return std::abs(a) < std::abs(b); }));
  return {a_, b_};
}

void test_flat() try {
  const Index NLayers_ = 100;
  const Index NQuad_   = 28;

  RandomNumberGenerator<> rng;
  auto draw = rng.get(0.00001, 0.99999);
  std::vector<Numeric> taus{0.5};
  for (Index i = 0; i < NLayers_ - 1; i++) taus.push_back(taus.back() + draw());
  for (Index i = 0; i < NLayers_; i++) {
    taus[i] *= 20.0 / taus.back();
  }

  const AscendingGrid tau_arr{taus};
  Vector omega_arr(NLayers_);
  for (auto& o : omega_arr) o = draw();
  Matrix Leg_coeffs_all(tau_arr.size(), 32);
  for (auto&& v : Leg_coeffs_all)
    v = {1.00000000e+00, 7.50000000e-01, 5.62500000e-01, 4.21875000e-01,
         3.16406250e-01, 2.37304688e-01, 1.77978516e-01, 1.33483887e-01,
         1.00112915e-01, 7.50846863e-02, 5.63135147e-02, 4.22351360e-02,
         3.16763520e-02, 2.37572640e-02, 1.78179480e-02, 1.33634610e-02,
         1.00225958e-02, 7.51694682e-03, 5.63771011e-03, 4.22828259e-03,
         3.17121194e-03, 2.37840895e-03, 1.78380672e-03, 1.33785504e-03,
         1.00339128e-03, 7.52543458e-04, 5.64407594e-04, 4.23305695e-04,
         3.17479271e-04, 2.38109454e-04, 1.78582090e-04, 1.33936568e-04};

  //const Numeric mu0 = 0.6;
  //const Numeric I0 = Constant::pi / mu0;
  //const Numeric phi0 = 0.9 * Constant::pi;
  Matrix b_neg(NQuad_, NQuad_ / 2, 0);
  b_neg[0] = 1;
  Matrix b_pos(NQuad_, NQuad_ / 2, 0);
  b_pos[0] = 1;
  const std::vector<disort::BDRF> BDRF_Fourier_modes{
      disort::BDRF{[](auto c, auto&, auto&) { c = 1; }}};
  Matrix s_poly_coeffs(tau_arr.size(), 2);
  for (auto&& v : s_poly_coeffs) v = {172311.79936609, -102511.4417051};
  const Vector f_arr{Leg_coeffs_all(joker, NQuad_)};

  // Optional (unused)
  const Index NLeg      = NQuad_;
  const Index NFourier_ = NQuad_;

  const disort::main_data dis(NQuad_,
                              NLeg,
                              NFourier_,
                              tau_arr,
                              omega_arr,
                              Leg_coeffs_all,
                              b_pos,
                              b_neg,
                              f_arr,
                              s_poly_coeffs,
                              BDRF_Fourier_modes,
                              mu0,
                              I0,
                              phi0);

  const Index NP   = 500;
  const Vector phi = uniform_grid(0, NP, Constant::pi / NP);

  Tensor3 u1(NLayers_, NP, NQuad_), u2(NLayers_, NP, NQuad_),
      u3(NLayers_, NP, NQuad_);
  {
    DebugTime norm{"1-by-1 U"};
    disort::u_data data_u;
    for (Index i = 0; i < NLayers_; i++) {
      for (Index j = 0; j < NP; j++) {
        dis.u(data_u, tau_arr[i], phi[j]);
        u1(i, j, joker) = data_u.intensities;
      }
    }
  }

  {
    DebugTime norm{"Gridded U"};
    dis.gridded_u(u2, phi);
  }

  {
    DebugTime norm{"Ungridded U"};
    dis.ungridded_u(u3, tau_arr, phi);
  }

  const auto [u2_abs, u2_rel] = absrel(u2.flat_view(), u1.flat_view());
  std::cout << "u abs-max gridded:   " << u2_abs << '\n';
  std::cout << "u abs-rel gridded:   " << u2_rel << '\n';

  const auto [u3_abs, u3_rel] = absrel(u3.flat_view(), u1.flat_view());
  std::cout << "u abs-max ungridded: " << u3_abs << '\n';
  std::cout << "u abs-rel ungridded: " << u3_rel << '\n';

  Vector fu1(NLayers_), fu2(NLayers_), fu3(NLayers_);
  Vector fd1(NLayers_), fd2(NLayers_), fd3(NLayers_);
  Vector fb1(NLayers_), fb2(NLayers_), fb3(NLayers_);
  {
    DebugTime norm{"1-by-1 Flux"};
    disort::flux_data data_flux;
    for (Index i = 0; i < NLayers_; i++) {
      fu1[i]      = dis.flux_up(data_flux, tau_arr[i]);
      auto [d, b] = dis.flux_down(data_flux, tau_arr[i]);
      fd1[i]      = d;
      fb1[i]      = b;
    }
  }

  {
    DebugTime norm{"Gridded Flux"};
    dis.gridded_flux(fu2, fd2, fb2);
  }

  {
    DebugTime norm{"Ungridded Flux"};
    dis.ungridded_flux(fu3, fd3, fb3, tau_arr);
  }

  const auto [fu2_abs, fu2_rel] = absrel(fu2.flat_view(), fu1.flat_view());
  const auto [fd2_abs, fd2_rel] = absrel(fd2.flat_view(), fd1.flat_view());
  const auto [fb2_abs, fb2_rel] = absrel(fb2.flat_view(), fb1.flat_view());
  std::cout << "fu abs-max gridded:   " << fu2_abs << '\n';
  std::cout << "fu abs-rel gridded:   " << fu2_rel << '\n';
  std::cout << "fd abs-max gridded:   " << fd2_abs << '\n';
  std::cout << "fd abs-rel gridded:   " << fd2_rel << '\n';
  std::cout << "fb abs-max gridded:   " << fb2_abs << '\n';
  std::cout << "fb abs-rel gridded:   " << fb2_rel << '\n';

  const auto [fu3_abs, fu3_rel] = absrel(fu3.flat_view(), fu1.flat_view());
  const auto [fd3_abs, fd3_rel] = absrel(fd3.flat_view(), fd1.flat_view());
  const auto [fb3_abs, fb3_rel] = absrel(fb3.flat_view(), fb1.flat_view());
  std::cout << "fu abs-max ungridded: " << fu3_abs << '\n';
  std::cout << "fu abs-rel ungridded: " << fu3_rel << '\n';
  std::cout << "fd abs-max ungridded: " << fd3_abs << '\n';
  std::cout << "fd abs-rel ungridded: " << fd3_rel << '\n';
  std::cout << "fb abs-max ungridded: " << fb3_abs << '\n';
  std::cout << "fb abs-rel ungridded: " << fb3_rel << '\n';
} catch (std::exception& e) {
  throw std::runtime_error(var_string("Error in test-test:\n", e.what()));
}

void handle_opt(const char* c, bool& print, Index& N) {
  std::string s = c;
  std::cout << s << '\n';

  if (s == "print")
    print = true;
  else
    N = std::stoi(s);
}

int main(int argc, char** argv) {
  bool print = false;
  Index N    = 100;
  Index argi = 1;
  while (argc > argi) {
    handle_opt(argv[argi], print, N);
    argi++;
  }

  if (print) {
    oldimpl(print);
    newimpl(print);
  } else {
    Array<Timing> ts1;
    ts1.reserve(N);
    Array<Timing> ts2;
    ts2.reserve(N);

    for (Index i = 0; i < N; i++) {
      ts1.emplace_back("old")([]() { oldimpl(); });
    }
    std::cout << "old best: "
              << std::ranges::min(
                     ts1,
                     [](auto a, auto b) { return a.count() < b.count(); },
                     &Timing::dt)
                     .dt
              << '\n';

    for (Index i = 0; i < N; i++) {
      ts2.emplace_back("new")([]() { newimpl(); });
    }
    std::cout << "new best: "
              << std::ranges::min(
                     ts2,
                     [](auto a, auto b) { return a.count() < b.count(); },
                     &Timing::dt)
                     .dt
              << '\n';
  }

  test_flat();
}

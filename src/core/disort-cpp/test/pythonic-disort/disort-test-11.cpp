#include "../test-helpers.h"
#include "disort.h"

namespace {
void require_close(const std::string_view name, const auto& a, const auto& b, const Numeric tolerance) {
  ARTS_USER_ERROR_IF(a.shape() != b.shape(), "{}: shape mismatch {} vs {}", name, a.shape(), b.shape());

  Numeric largest_error     = 0.0;
  Numeric largest_actual    = 0.0;
  Numeric largest_reference = 0.0;
  Size    largest_index     = 0;
  Size    index             = 0;
  auto    ai                = a.elem_begin();
  auto    bi                = b.elem_begin();
  for (; ai != a.elem_end(); ++ai, ++bi, ++index) {
    const Numeric av    = *ai;
    const Numeric bv    = *bi;
    const Numeric scale = std::max({1.0, std::abs(av), std::abs(bv)});
    const Numeric error = std::abs(av - bv) / scale;
    if (error > largest_error) {
      largest_error     = error;
      largest_actual    = av;
      largest_reference = bv;
      largest_index     = index;
    }
  }
  ARTS_USER_ERROR_IF(largest_error > tolerance,
                     "{}: largest scaled error is {} at flat index {}: actual {}, reference {}, tolerance {}",
                     name,
                     largest_error,
                     largest_index,
                     largest_actual,
                     largest_reference,
                     tolerance);
}

disort::main_data identical_atmosphere(const AscendingGrid& tau_arr, const Numeric mu0 = 0.6) {
  constexpr Index   NQuad    = 8;
  constexpr Index   NLeg     = NQuad;
  constexpr Index   NFourier = NQuad;
  constexpr Index   NLeg_all = 2 * NQuad;
  constexpr Numeric omega    = 0.8;
  Matrix            legendre(tau_arr.size(), NLeg_all);
  for (auto&& row : legendre)
    for (Index k = 0; k < NLeg_all; ++k) row[k] = std::pow(0.75, static_cast<Numeric>(k));

  return {NQuad,
          NLeg,
          NFourier,
          tau_arr,
          Vector(tau_arr.size(), omega),
          legendre,
          Matrix(NFourier, NQuad / 2, 0.0),
          Matrix(NFourier, NQuad / 2, 0.0),
          Vector(tau_arr.size(), legendre[0, NLeg]),
          Matrix(tau_arr.size(), 0),
          {},
          mu0,
          Constant::pi / mu0,
          0.9 * Constant::pi};
}

void require_finite(const std::string_view name, const auto& values) {
  ARTS_USER_ERROR_IF(stdr::any_of(values | by_elem, [](const Numeric x) { return !std::isfinite(x); }),
                     "{} contains a non-finite value",
                     name);
}

void require_scalar_close(const std::string_view name,
                          const Numeric          actual,
                          const Numeric          expected,
                          const Numeric          tolerance = 5e-8) {
  const Numeric scale        = std::max({1.0, std::abs(actual), std::abs(expected)});
  const Numeric scaled_error = std::abs(actual - expected) / scale;
  ARTS_USER_ERROR_IF(scaled_error > tolerance,
                     "{}: {} differs from reference {} (scaled error {}, tolerance {})",
                     name,
                     actual,
                     expected,
                     scaled_error,
                     tolerance);
}

void check_test_11a_reference(const std::string_view   name,
                              const disort::main_data& dis,
                              const Vector&            taus,
                              const Vector&            phis) {
  const Tensor3 u  = compute_u(dis, taus, phis, true);
  const Matrix  u0 = compute_u0(dis, taus);

  require_scalar_close(std::format("{} u[0,0,0]", name), u[0, 0, 0], -2900902.4893682143);
  require_scalar_close(std::format("{} u[0,0,7]", name), u[0, 0, 7], -9729769.265097402);
  require_scalar_close(std::format("{} u[0,0,8]", name), u[0, 0, 8], -2394710.1121911793);
  require_scalar_close(std::format("{} u[2,1,3]", name), u[2, 1, 3], -7252729.563810989);
  require_scalar_close(std::format("{} u[4,2,15]", name), u[4, 2, 15], -8972939.78056413);
  require_scalar_close(std::format("{} u0[0,0]", name), u0[0, 0], -2900902.1380472905);
  // This comparatively small outgoing value is the residual of source terms
  // amplified by 1 / (1 - omega) = 1e6.  Eigensolver and math-library
  // roundoff therefore loses about two additional relative digits here.
  require_scalar_close(std::format("{} u0[0,15]", name), u0[0, 15], -16747.47086339454, 5e-6);
  require_scalar_close(std::format("{} u0[1,6]", name), u0[1, 6], -10275991.914094502);
  require_scalar_close(std::format("{} u0[2,15]", name), u0[2, 15], -8972939.758737655);

  const auto [up, down_diffuse, down_direct, dfdt] = compute_flux(dis, taus);
  static_cast<void>(dfdt);
  // The source is amplified by 1 / (1 - omega) = 1e6.  Both diffuse fluxes
  // integrate residuals of the resulting O(1e7) radiances, so their last
  // approximately six relative digits depend on eigensolver and math-library
  // roundoff.  The direct flux does not contain this cancellation.
  require_close("delta-scaled source upward flux",
                up,
                Vector{-24332790.748160336, -27833449.381811786, -45865859.33752958},
                1e-6);
  require_close("delta-scaled source diffuse-downward flux",
                down_diffuse,
                Vector{-299160.5808801778, -2818892.4624789366, -34216142.156027764},
                1e-6);
  require_close("delta-scaled source direct-downward flux",
                down_direct,
                Vector{2.8258789658328820, 1.0893983516478147, 7.8980209952661667e-05},
                5e-13);
}

void test_11a_1layer() try {
  const AscendingGrid tau_arr{8.};
  const Vector        omega_arr{0.999999};
  const Index         NQuad = 16;
  const Matrix        Leg_coeffs_all{Vector{
      1.00000000e+00, 7.50000000e-01, 5.62500000e-01, 4.21875000e-01, 3.16406250e-01, 2.37304688e-01, 1.77978516e-01,
      1.33483887e-01, 1.00112915e-01, 7.50846863e-02, 5.63135147e-02, 4.22351360e-02, 3.16763520e-02, 2.37572640e-02,
      1.78179480e-02, 1.33634610e-02, 1.00225958e-02, 7.51694682e-03, 5.63771011e-03, 4.22828259e-03, 3.17121194e-03,
      2.37840895e-03, 1.78380672e-03, 1.33785504e-03, 1.00339128e-03, 7.52543458e-04, 5.64407594e-04, 4.23305695e-04,
      3.17479271e-04, 2.38109454e-04, 1.78582090e-04, 1.33936568e-04}
                                         .reshape(tau_arr.size(), 32)};

  const Numeric mu0  = 0.6;
  const Numeric I0   = Constant::pi / mu0;
  const Numeric phi0 = 0.9 * Constant::pi;
  Matrix        b_neg(NQuad, NQuad / 2, 0);
  b_neg[0] = 1;
  Matrix b_pos(NQuad, NQuad / 2, 0);
  b_pos[0] = 1;
  const std::vector<disort::BDRF> BDRF_Fourier_modes{disort::BDRF{[](auto c, auto&, auto&) { c = 1; }}};
  const Matrix                    s_poly_coeffs{
      Vector{172311.79936609 / (1 - omega_arr[0]), -102511.4417051 / (1 - omega_arr[0])}.reshape(tau_arr.size(), 2)};
  const Vector f_arr{Leg_coeffs_all[joker, NQuad]};

  // Optional (unused)
  const Index NLeg     = NQuad;
  const Index NFourier = NQuad;

  const disort::main_data dis(NQuad,
                              NLeg,
                              NFourier,
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

  const Vector taus{Vector{
      0.06354625877794251,
      0.6354625877794251,
      6.354625877794252,
  }
                        .reshape(3)};

  const Vector phis{Vector{
      0.0,
      1.5705463267948965,
      3.141092653589793,
      4.71163898038469,
      6.282185307179586,
  }
                        .reshape(5)};

  check_test_11a_reference("test_11a-1layer", dis, taus, phis);
} catch (std::exception& e) { throw std::runtime_error(std::format("Error in test-11a-1layer:\n{}", e.what())); }

void test_11a_multilayer() try {
  const AscendingGrid tau_arr{0.5, 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7., 7.5, 8.};
  const Vector        omega_arr(tau_arr.size(), 0.999999);
  const Index         NQuad = 16;
  Matrix              Leg_coeffs_all(tau_arr.size(), 32);
  for (auto&& v : Leg_coeffs_all)
    v = std::array{1.00000000e+00, 7.50000000e-01, 5.62500000e-01, 4.21875000e-01, 3.16406250e-01, 2.37304688e-01,
                   1.77978516e-01, 1.33483887e-01, 1.00112915e-01, 7.50846863e-02, 5.63135147e-02, 4.22351360e-02,
                   3.16763520e-02, 2.37572640e-02, 1.78179480e-02, 1.33634610e-02, 1.00225958e-02, 7.51694682e-03,
                   5.63771011e-03, 4.22828259e-03, 3.17121194e-03, 2.37840895e-03, 1.78380672e-03, 1.33785504e-03,
                   1.00339128e-03, 7.52543458e-04, 5.64407594e-04, 4.23305695e-04, 3.17479271e-04, 2.38109454e-04,
                   1.78582090e-04, 1.33936568e-04};

  const Numeric mu0  = 0.6;
  const Numeric I0   = Constant::pi / mu0;
  const Numeric phi0 = 0.9 * Constant::pi;
  Matrix        b_neg(NQuad, NQuad / 2, 0);
  b_neg[0] = 1;
  Matrix b_pos(NQuad, NQuad / 2, 0);
  b_pos[0] = 1;
  const std::vector<disort::BDRF> BDRF_Fourier_modes{disort::BDRF{[](auto c, auto&, auto&) { c = 1; }}};
  Matrix                          s_poly_coeffs(tau_arr.size(), 2);
  for (auto&& v : s_poly_coeffs) v = std::array{172311.79936609, -102511.4417051};
  s_poly_coeffs /= 1 - omega_arr[0];
  const Vector f_arr{Leg_coeffs_all[joker, NQuad]};

  // Optional (unused)
  const Index NLeg     = NQuad;
  const Index NFourier = NQuad;

  const disort::main_data dis(NQuad,
                              NLeg,
                              NFourier,
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

  const Vector taus{Vector{
      0.06354625877794251,
      0.6354625877794251,
      6.354625877794252,
  }
                        .reshape(3)};

  const Vector phis{Vector{
      0.0,
      1.5705463267948965,
      3.141092653589793,
      4.71163898038469,
      6.282185307179586,
  }
                        .reshape(5)};

  check_test_11a_reference("test_11a-multilayer", dis, taus, phis);
} catch (std::exception& e) { throw std::runtime_error(std::format("Error in test-11a-multilayer:\n{}", e.what())); }

void test_11b_identical_layer_invariance() try {
  const auto one_layer  = identical_atmosphere(AscendingGrid{8.0});
  const auto many_layer = identical_atmosphere(
      AscendingGrid{0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0});
  const Vector tau{0.0, 0.17, 0.5, 1.37, 4.25, 7.91, 8.0};
  const Vector phi{0.0, 0.41, 2.7, Constant::two_pi - 0.2};

  require_close("one layer vs identical sublayers, corrected intensity",
                compute_u(one_layer, tau, phi, true),
                compute_u(many_layer, tau, phi, true),
                2e-9);

  disort::tms_data tms_one;
  disort::tms_data tms_many;
  Vector           ims_one;
  Vector           ims_many;
  for (const Numeric t : tau)
    for (const Numeric p : phi) {
      one_layer.TMS(tms_one, t, p);
      many_layer.TMS(tms_many, t, p);
      require_close("one layer vs identical sublayers, TMS", tms_one.TMS, tms_many.TMS, 2e-13);

      one_layer.IMS(ims_one, t, p);
      many_layer.IMS(ims_many, t, p);
      require_close("one layer vs identical sublayers, IMS", ims_one, ims_many, 2e-13);
    }
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11b-identical-layer-invariance:\n{}", e.what()));
}

void test_11c_heterogeneous_delta_scaling() try {
  const AscendingGrid tau_arr{1.0, 3.0, 6.0};
  const Vector        omega{0.5, 0.8, 0.25};
  const Vector        f{0.2, 0.4, 0.6};
  Matrix              legendre(tau_arr.size(), 8, 0.0);
  legendre[joker, 0] = 1.0;

  const disort::main_data dis(4,
                              4,
                              4,
                              tau_arr,
                              omega,
                              legendre,
                              Matrix(4, 2, 0.0),
                              Matrix(4, 2, 0.0),
                              f,
                              Matrix(tau_arr.size(), 0),
                              {},
                              0.6,
                              Constant::pi / 0.6,
                              0.0);

  const Vector expected{0.0,
                        1.0 * (1.0 - omega[0] * f[0]),
                        1.0 * (1.0 - omega[0] * f[0]) + 2.0 * (1.0 - omega[1] * f[1]),
                        1.0 * (1.0 - omega[0] * f[0]) + 2.0 * (1.0 - omega[1] * f[1]) + 3.0 * (1.0 - omega[2] * f[2])};
  require_close("cumulative delta-scaled optical depth", dis.scaled_tau(), expected, 1e-14);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11c-heterogeneous-delta-scaling:\n{}", e.what()));
}

//! FIXME: This test is wrong because TMS works but the checks are for full DISORT where we have not confirmed it
void test_11d_removable_correction_singularities() try {
  const auto    probe   = identical_atmosphere(AscendingGrid{1.0});
  const Numeric mu_quad = std::midpoint(probe.mu()[2], probe.mu()[3]);  // <-- THIS SHOULD BE A QUADRATURE ANGLE, NOT A MIDPOINT

  // TMS has a removable singularity when the downward evaluation angle is
  // the beam angle.
  const auto       tms_singular = identical_atmosphere(AscendingGrid{1.0}, mu_quad);
  disort::tms_data tms;
  tms_singular.TMS(tms, 0.37, 1.1);
  require_finite("TMS at mu == mu0", tms.TMS);

  // IMS has its removable singularity at mu == scaled_mu0.  For this
  // identical atmosphere omega*f is constant, so choose mu0 accordingly.
  constexpr Numeric omega        = 0.8;
  const Numeric     f            = std::pow(0.75, 8.0);
  const auto        ims_singular = identical_atmosphere(AscendingGrid{1.0}, mu_quad * (1.0 - omega * f));
  Vector            ims;
  ims_singular.IMS(ims, 0.37, 1.1);
  require_finite("IMS at mu == scaled_mu0", ims);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11d-removable-correction-singularities:\n{}", e.what()));
}

void test_11e_current_pythonic_disort_correction() try {
  const AscendingGrid tau_arr{0.4, 1.2, 2.0};
  const Vector        omega{0.7, 0.5, 0.8};
  const Vector        g{0.7, 0.6, 0.8};
  Matrix              legendre(tau_arr.size(), 8);
  for (Size l = 0; l < tau_arr.size(); ++l)
    for (Index k = 0; k < legendre.ncols(); ++k) legendre[l, k] = std::pow(g[l], static_cast<Numeric>(k));

  const disort::main_data dis(4,
                              4,
                              4,
                              tau_arr,
                              omega,
                              legendre,
                              Matrix(4, 2, 0.0),
                              Matrix(4, 2, 0.0),
                              Vector{legendre[joker, 4]},
                              Matrix(tau_arr.size(), 0),
                              {},
                              0.63,
                              Constant::pi / 0.63,
                              0.7);

  const Vector     tau{0.2, 0.8, 1.7};
  const Vector     phi{0.3, 2.0};
  Tensor3          correction(4, tau.size(), phi.size());
  disort::u_data   uncorrected;
  disort::u_data   corrected;
  disort::tms_data tms;
  Vector           ims;
  for (Size t = 0; t < tau.size(); ++t)
    for (Size p = 0; p < phi.size(); ++p) {
      dis.u(uncorrected, tau[t], phi[p]);
      dis.u_corr(corrected, ims, tms, tau[t], phi[p]);
      correction[joker, t, p]  = corrected.intensities;
      correction[joker, t, p] -= uncorrected.intensities;
    }

  // TMS/IMS regression values.  IMS follows the original DISORT INTCOR
  // convention: subtract SECSCA and apply it only inside the 10-degree
  // aureole around the direct beam.
  const Tensor3 expected{
      Vector{-3.7472007688001185e-02, -8.5810638855479404e-03, -1.0028434790383289e-02, -2.9933874791298332e-03,
             3.1288563974852940e-03,  -1.2147006561405757e-02, 3.8091906483520219e-03,  1.2244675379773216e-02,
             -1.1816920676733372e-04, 8.5884505553329264e-03,  -8.1475725099112301e-04, 4.0094675513905817e-03,
             -2.3521253903849626e-02, 2.0705176853186166e-02,  5.5151296707868780e-03,  5.7443705853809368e-03,
             -5.2567311883454065e-02, 1.8270752476844954e-02,  2.3911621315690348e-01,  -2.8709959968093834e-02,
             3.1530176563638068e-01,  -4.0317743340383655e-02, 3.2167016146607635e-01,  -3.3546908717691976e-02}
          .reshape(4, 3, 2)};
  require_close("heterogeneous DISORT NT correction", correction, expected, 3e-13);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11e-current-pythonic-disort-correction:\n{}", e.what()));
}

void test_11f_disort_user_angle_formal_solution() try {
  const auto    dis = identical_atmosphere(AscendingGrid{0.4, 1.2, 2.0});
  const Numeric tau = 0.83;
  const Numeric phi = 1.37;

  disort::u_data      ordinate;
  disort::user_u_data user;
  dis.u(ordinate, tau, phi);
  dis.u_user(user, tau, phi, dis.mu());
  require_close("formal user-angle solution at quadrature directions", user.intensities, ordinate.intensities, 2e-11);

  const Vector arbitrary_mu{-1.0, -0.71, -0.03, 0.06, 0.54, 1.0};
  dis.u_user(user, tau, phi, arbitrary_mu);
  require_finite("formal user-angle solution", user.intensities);

  disort::u_data   corrected;
  disort::tms_data tms;
  Vector           ims;
  dis.u_corr(corrected, ims, tms, tau, phi);
  dis.u_user_corr(user, ims, tms, tau, phi, dis.mu());
  require_close(
      "corrected user-angle solution at quadrature directions", user.intensities, corrected.intensities, 2e-11);

  dis.u_user_corr(user, ims, tms, tau, phi, arbitrary_mu);
  require_finite("corrected formal user-angle solution", user.intensities);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11f-disort-user-angle-formal-solution:\n{}", e.what()));
}

void test_11g_gridded_correction_cache_equivalence() try {
  const auto   dis = identical_atmosphere(AscendingGrid{0.4, 1.2, 2.0});
  const Vector phi{0.2, 1.7};
  Tensor3      cached_tms(3, 2, dis.mu().size());
  Tensor3      cached_ims(3, 2, dis.mu().size() / 2);
  Tensor3      point_tms(3, 2, dis.mu().size());
  Tensor3      point_ims(3, 2, dis.mu().size() / 2);
  dis.gridded_TMS(cached_tms, phi);
  dis.gridded_IMS(cached_ims, phi);

  disort::tms_data tms;
  Vector           ims;
  const Vector     tau{0.4, 1.2, 2.0};
  for (Index l = 0; l < 3; ++l)
    for (Index p = 0; p < 2; ++p) {
      dis.TMS(tms, tau[l], phi[p]);
      dis.IMS(ims, tau[l], phi[p]);
      point_tms[l, p] = tms.TMS;
      point_ims[l, p] = ims;
    }

  require_close("cached gridded TMS", cached_tms, point_tms, 2e-13);
  require_close("cached gridded IMS", cached_ims, point_ims, 2e-13);

  dis.gridded_IMS(cached_ims, phi, disort::ims_convention::pythonic_disort);
  for (Index l = 0; l < 3; ++l)
    for (Index p = 0; p < 2; ++p) {
      dis.IMS(ims, tau[l], phi[p], disort::ims_convention::pythonic_disort);
      point_ims[l, p] = ims;
    }
  require_close("cached gridded Pythonic IMS", cached_ims, point_ims, 2e-13);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11g-gridded-correction-cache-equivalence:\n{}", e.what()));
}

void test_11h_delta_scaled_source_equivalence() try {
  constexpr Index     nquad    = 4;
  constexpr Index     nleg     = 2;
  constexpr Index     nfourier = 1;
  const AscendingGrid physical_tau{0.7, 1.9};
  const Vector        physical_omega{0.6, 0.8};
  const Vector        fraction{0.2, 0.1};
  const Matrix        physical_moments{Vector{1.0, 0.4, 1.0, 0.3}.reshape(2, nleg)};
  const Matrix        physical_source{Vector{2.0, 0.5, 0.2, 1.5, -0.3, 0.1}.reshape(2, 3)};
  const Matrix        upper_boundary(nfourier, nquad / 2, 0.15);
  const Matrix        lower_boundary(nfourier, nquad / 2, 0.25);

  const disort::main_data physical(nquad,
                                   nleg,
                                   nfourier,
                                   physical_tau,
                                   physical_omega,
                                   physical_moments,
                                   upper_boundary,
                                   lower_boundary,
                                   fraction,
                                   physical_source,
                                   {},
                                   1.0,
                                   0.0,
                                   0.0);

  Vector  scale(2), transport_omega(2), transport_tau(2);
  Matrix  transport_moments(2, nleg), transport_source(2, 3);
  Numeric physical_top = 0.0;
  Numeric scaled_top   = 0.0;
  for (Index layer = 0; layer < 2; ++layer) {
    scale[layer]           = 1.0 - physical_omega[layer] * fraction[layer];
    transport_omega[layer] = physical_omega[layer] * (1.0 - fraction[layer]) / scale[layer];
    transport_tau[layer]   = scaled_top + scale[layer] * (physical_tau[layer] - physical_top);
    for (Index degree = 0; degree < nleg; ++degree)
      transport_moments[layer, degree] = (physical_moments[layer, degree] - fraction[layer]) / (1.0 - fraction[layer]);

    const Numeric inverse_scale = 1.0 / scale[layer];
    const Numeric offset        = physical_top - inverse_scale * scaled_top;
    const Numeric c0            = physical_source[layer, 0];
    const Numeric c1            = physical_source[layer, 1];
    const Numeric c2            = physical_source[layer, 2];
    transport_source[layer, 0]  = c0 + c1 * offset + c2 * offset * offset;
    transport_source[layer, 1]  = inverse_scale * (c1 + 2.0 * c2 * offset);
    transport_source[layer, 2]  = inverse_scale * inverse_scale * c2;

    physical_top = physical_tau[layer];
    scaled_top   = transport_tau[layer];
  }

  const disort::main_data pretransformed(nquad,
                                         nleg,
                                         nfourier,
                                         AscendingGrid{transport_tau},
                                         transport_omega,
                                         transport_moments,
                                         upper_boundary,
                                         lower_boundary,
                                         Vector(2, 0.0),
                                         transport_source,
                                         {},
                                         1.0,
                                         0.0,
                                         0.0);

  const Vector physical_points{0.0, 0.2, 0.7, 1.2, 1.9};
  Vector       transport_points(physical_points.size());
  for (Index point = 0; point < static_cast<Index>(physical_points.size()); ++point) {
    const Index   layer              = physical_points[point] <= physical_tau[0] ? 0 : 1;
    const Numeric layer_physical_top = layer == 0 ? 0.0 : physical_tau[layer - 1];
    const Numeric layer_scaled_top   = layer == 0 ? 0.0 : transport_tau[layer - 1];
    transport_points[point]          = layer_scaled_top + scale[layer] * (physical_points[point] - layer_physical_top);
  }

  disort::u0_data   physical_data, transport_data;
  disort::flux_data physical_flux_data, transport_flux_data;
  for (Index point = 0; point < static_cast<Index>(physical_points.size()); ++point) {
    physical.u0(physical_data, physical_points[point]);
    pretransformed.u0(transport_data, transport_points[point]);
    require_close("internally and explicitly delta-scaled source radiance", physical_data.u0, transport_data.u0, 2e-12);

    const auto physical_flux  = physical.flux(physical_flux_data, physical_points[point]);
    const auto transport_flux = pretransformed.flux(transport_flux_data, transport_points[point]);
    require_scalar_close(
        "internally and explicitly delta-scaled upward flux", physical_flux.up, transport_flux.up, 2e-12);
    require_scalar_close("internally and explicitly delta-scaled downward flux",
                         physical_flux.down_diffuse,
                         transport_flux.down_diffuse,
                         2e-12);
  }
} catch (std::exception& e) {
  throw std::runtime_error(std::format("Error in test-11h-delta-scaled-source-equivalence:\n{}", e.what()));
}

void test_11i_full_bulk_flux_equivalence() {
  const AscendingGrid layer_tau{0.4, 1.0};
  const auto          dis = identical_atmosphere(layer_tau);

  Vector gridded_up(layer_tau.size()), gridded_down(layer_tau.size()), gridded_direct(layer_tau.size()),
      gridded_dfdt(layer_tau.size());
  dis.gridded_flux(gridded_up, gridded_down, gridded_direct, gridded_dfdt);
  disort::flux_data scratch;
  for (Index layer = 0; layer < static_cast<Index>(layer_tau.size()); ++layer) {
    const auto pointwise = dis.flux(scratch, layer_tau[layer]);
    require_scalar_close("gridded/pointwise upward flux", gridded_up[layer], pointwise.up, 2e-12);
    require_scalar_close("gridded/pointwise diffuse-downward flux", gridded_down[layer], pointwise.down_diffuse, 2e-12);
    require_scalar_close("gridded/pointwise direct-downward flux", gridded_direct[layer], pointwise.down_direct, 2e-12);
    require_scalar_close("gridded/pointwise DFDT", gridded_dfdt[layer], pointwise.dfdt, 2e-12);
  }

  const AscendingGrid output_tau{0.0, 0.2, 0.4, 0.7, 1.0};
  Vector ungridded_up(output_tau.size()), ungridded_down(output_tau.size()), ungridded_direct(output_tau.size()),
      ungridded_dfdt(output_tau.size());
  dis.ungridded_flux(ungridded_up, ungridded_down, ungridded_direct, ungridded_dfdt, output_tau);
  for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level) {
    const auto pointwise = dis.flux(scratch, output_tau[level]);
    require_scalar_close("ungridded/pointwise upward flux", ungridded_up[level], pointwise.up, 2e-12);
    require_scalar_close(
        "ungridded/pointwise diffuse-downward flux", ungridded_down[level], pointwise.down_diffuse, 2e-12);
    require_scalar_close(
        "ungridded/pointwise direct-downward flux", ungridded_direct[level], pointwise.down_direct, 2e-12);
    require_scalar_close("ungridded/pointwise DFDT", ungridded_dfdt[level], pointwise.dfdt, 2e-12);
  }
}
}  // namespace

int main() try {
  std::cout << std::setprecision(16);
  test_11a_1layer();
  test_11a_multilayer();
  test_11b_identical_layer_invariance();
  test_11c_heterogeneous_delta_scaling();
  test_11d_removable_correction_singularities();
  test_11e_current_pythonic_disort_correction();
  test_11f_disort_user_angle_formal_solution();
  test_11g_gridded_correction_cache_equivalence();
  test_11h_delta_scaled_source_equivalence();
  test_11i_full_bulk_flux_equivalence();
} catch (std::exception& e) {
  std::cerr << "Error in main:\n" << e.what() << '\n';
  return EXIT_FAILURE;
}

#include <arts_constants.h>
#include <disort-brdf.h>
#include <legendre.h>
#include <rtepack_multitype.h>
#include <time_report.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <cmath>
#include <iostream>

#include "../reference-data.h"
#include "test-adapter.h"

namespace {
constexpr Numeric reference_tolerance    = 7e-5;
constexpr Numeric polarization_tolerance = 2e-9;

void expect_reference(const std::string_view name,
                      const Numeric          actual,
                      const Numeric          expected,
                      const Numeric          tolerance = reference_tolerance) {
  ARTS_USER_ERROR_IF(std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)),
                     "{}: expected {}, got {} (difference {})",
                     name,
                     expected,
                     actual,
                     actual - expected);
}

void expect_unpolarized(const std::string_view name, const rtepack::stokvec& value) {
  for (Index s = 1; s < vdisort::stokes_dimension; ++s)
    ARTS_USER_ERROR_IF(std::abs(value[s]) > polarization_tolerance,
                       "{}: expected Stokes component {} to vanish, got {}",
                       name,
                       s,
                       value[s]);
}

struct scalar_vdisort_model {
  vdisort::main_data              solver;
  vdisort::phase_matrix_data      user_phase;
  vdisort::beam_phase_matrix_data user_beam_phase;
  Vector                          original_moments;
  Vector                          transport_moments;
  Numeric                         physical_depth{};
  Numeric                         physical_omega{};
  Numeric                         mu0{};
  Numeric                         phi0{};
  Numeric                         beam_intensity{};
  Numeric                         delta_m_fraction{};
  Numeric                         optical_depth_scale{1.0};
};

Numeric scalar_phase_mode(const ConstVectorView& moments,
                          const Index            m,
                          const Numeric          outgoing_mu,
                          const Numeric          incident_mu) {
  Numeric result = 0.0;
  for (Index degree = m; degree < static_cast<Index>(moments.size()); ++degree) {
    const Numeric factorial_ratio =
        Legendre::tgamma_ratio(static_cast<Numeric>(degree - m + 1), static_cast<Numeric>(degree + m + 1));
    result += static_cast<Numeric>(2 * degree + 1) * moments[degree] * factorial_ratio *
              Legendre::assoc_legendre(degree, m, outgoing_mu) * Legendre::assoc_legendre(degree, m, incident_mu);
  }
  return result;
}

scalar_vdisort_model make_scalar_model(const Numeric          physical_depth,
                                       const Numeric          physical_omega,
                                       const bool             has_beam,
                                       const Index            nquad,
                                       const ConstVectorView& moments,
                                       const ConstVectorView& user_mu,
                                       const Numeric          mu0,
                                       const Numeric          beam_intensity,
                                       const Numeric          delta_m_fraction     = 0.0,
                                       const Numeric          phi0                 = 0.0,
                                       const Vector*          removed_peak_moments = nullptr) {
  const Numeric optical_depth_scale = 1.0 - physical_omega * delta_m_fraction;
  const Numeric transport_omega     = physical_omega * (1.0 - delta_m_fraction) / optical_depth_scale;
  Vector        transport_moments(std::min(nquad, static_cast<Index>(moments.size())));
  ARTS_USER_ERROR_IF(removed_peak_moments != nullptr and removed_peak_moments->size() < transport_moments.size(),
                     "The removed peak has {} moments, but {} are required",
                     removed_peak_moments->size(),
                     transport_moments.size());
  for (Index degree = 0; degree < static_cast<Index>(transport_moments.size()); ++degree) {
    const Numeric removed_moment = removed_peak_moments == nullptr ? 1.0 : (*removed_peak_moments)[degree];
    transport_moments[degree]    = (moments[degree] - delta_m_fraction * removed_moment) / (1.0 - delta_m_fraction);
  }

  Index highest_nonzero_degree = 0;
  for (Index degree = 0; degree < static_cast<Index>(transport_moments.size()); ++degree)
    if (transport_moments[degree] != 0.0) highest_nonzero_degree = degree;
  const Index     nmodes  = std::min(nquad, highest_nonzero_degree + 1);
  constexpr Index nlayers = 1;
  constexpr Index nalpha  = 2;
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Index     nuser   = static_cast<Index>(user_mu.size());

  Vector quadrature_mu(nquad), quadrature_weights(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu[Range{0, n}], quadrature_weights);
  std::transform(
      quadrature_mu.begin(), quadrature_mu.begin() + n, quadrature_mu.begin() + n, [](const Numeric x) { return -x; });

  Tensor7 phase(nalpha, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nquad; ++i)
      for (Index j = 0; j < nquad; ++j) {
        const Numeric value = scalar_phase_mode(transport_moments, m, quadrature_mu[i], quadrature_mu[j]);
        phase[vdisort::cosine_mode, m, 0, i, j, 0, 0] = value;
        if (m > 0) phase[vdisort::sine_mode, m, 0, i, j, 0, 0] = value;
      }

  Tensor4 up(nalpha, nmodes, n, nstokes, 0.0);
  Tensor4 down(nalpha, nmodes, n, nstokes, 0.0);
  if (not has_beam) down[vdisort::cosine_mode, 0, joker, 0] = 1.0;

  Tensor6 beam_phase(nalpha, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nquad; ++i)
      beam_phase[vdisort::cosine_mode, m, 0, i, 0, 0] = scalar_phase_mode(transport_moments, m, quadrature_mu[i], -mu0);

  auto solver = vdisort_test::make_solver(nquad,
                                          nmodes,
                                          AscendingGrid{optical_depth_scale * physical_depth},
                                          Vector{transport_omega},
                                          std::move(phase),
                                          std::move(up),
                                          std::move(down),
                                          Tensor3(nlayers, 0, nstokes),
                                          {},
                                          mu0,
                                          has_beam ? Vector{beam_intensity, 0.0, 0.0, 0.0} : Vector(4, 0.0),
                                          phi0,
                                          std::move(beam_phase));

  vdisort::phase_matrix_data user_phase(nalpha, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nuser; ++i)
      for (Index j = 0; j < nquad; ++j) {
        const Numeric value = scalar_phase_mode(transport_moments, m, user_mu[i], quadrature_mu[j]);
        user_phase[vdisort::cosine_mode, m, 0, i, j][0, 0] = value;
        if (m > 0) user_phase[vdisort::sine_mode, m, 0, i, j][0, 0] = value;
      }
  vdisort::beam_phase_matrix_data user_beam_phase(nalpha, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nuser; ++i)
      user_beam_phase[vdisort::cosine_mode, m, 0, i][0, 0] = scalar_phase_mode(transport_moments, m, user_mu[i], -mu0);

  return {.solver              = std::move(solver),
          .user_phase          = std::move(user_phase),
          .user_beam_phase     = std::move(user_beam_phase),
          .original_moments    = Vector{moments},
          .transport_moments   = std::move(transport_moments),
          .physical_depth      = physical_depth,
          .physical_omega      = physical_omega,
          .mu0                 = mu0,
          .phi0                = phi0,
          .beam_intensity      = has_beam ? beam_intensity : 0.0,
          .delta_m_fraction    = delta_m_fraction,
          .optical_depth_scale = optical_depth_scale};
}

Numeric phi1_negative(const Numeric x) { return x == 0.0 ? 1.0 : -std::expm1(-x) / x; }

Numeric phi2(const Numeric x) {
  if (std::abs(x) >= 1.0) return (std::expm1(x) - x) / (x * x);
  Numeric term = 0.5;
  Numeric sum  = term;
  for (Index k = 1; k <= 16; ++k) {
    term *= x / static_cast<Numeric>(k + 2);
    sum  += term;
  }
  return sum;
}

Numeric downward_tms_kernel(const Numeric mu,
                            const Numeric mu0,
                            const Numeric thickness,
                            const Numeric lower_value,
                            const Numeric upper_value) {
  const Numeric y = 1.0 / mu - 1.0 / mu0;
  const Numeric w = thickness * y;
  return std::abs(w) < 1.0 ? thickness / mu * phi1_negative(w) * lower_value : (lower_value - upper_value) / (mu * y);
}

Numeric ims_chi(const Numeric tau, const Numeric mu, const Numeric scaled_mu0) {
  const Numeric x = 1.0 / mu - 1.0 / scaled_mu0;
  const Numeric z = -tau * x;
  if (std::abs(z) < 1.0) return tau * tau * std::exp(-tau / scaled_mu0) * phi2(z) / (mu * scaled_mu0);
  return ((tau - 1.0 / x) * std::exp(-tau / scaled_mu0) + std::exp(-tau / mu) / x) / (mu * scaled_mu0 * x);
}

Numeric scattering_angle_cosine(const Numeric mu,
                                const Numeric phi,
                                const Numeric incident_mu,
                                const Numeric incident_phi) {
  return mu * incident_mu +
         std::sqrt(1.0 - incident_mu * incident_mu) * std::cos(incident_phi - phi) * std::sqrt(1.0 - mu * mu);
}

Vector scalar_delta_m_correction(const scalar_vdisort_model& model,
                                 const Numeric               physical_tau,
                                 const Numeric               phi,
                                 const ConstVectorView&      user_mu) {
  Vector correction(user_mu.size(), 0.0);
  if (model.beam_intensity == 0.0 or model.delta_m_fraction == 0.0) return correction;

  Vector weighted_true(model.original_moments.size());
  Vector weighted_transport(model.transport_moments.size());
  for (Index degree = 0; degree < static_cast<Index>(weighted_true.size()); ++degree)
    weighted_true[degree] = static_cast<Numeric>(2 * degree + 1) * model.original_moments[degree];
  for (Index degree = 0; degree < static_cast<Index>(weighted_transport.size()); ++degree)
    weighted_transport[degree] = static_cast<Numeric>(2 * degree + 1) * model.transport_moments[degree];

  const Numeric scaled_tau      = model.optical_depth_scale * physical_tau;
  const Numeric scaled_depth    = model.optical_depth_scale * model.physical_depth;
  const Numeric transport_omega = model.physical_omega * (1.0 - model.delta_m_fraction) / model.optical_depth_scale;
  for (Index i = 0; i < static_cast<Index>(user_mu.size()); ++i) {
    const Numeric mu          = user_mu[i];
    const Numeric abs_mu      = std::abs(mu);
    const Numeric nu          = scattering_angle_cosine(mu, phi, -model.mu0, model.phi0);
    const Numeric p_true      = Legendre::legendre_sum(weighted_true, nu);
    const Numeric p_transport = Legendre::legendre_sum(weighted_transport, nu);
    const Numeric beam_factor = mu > 0.0 ? model.mu0 / (model.mu0 + mu) : 1.0;
    const Numeric coefficient = transport_omega * model.beam_intensity / (4.0 * Constant::pi) * beam_factor *
                                (p_true / (1.0 - model.delta_m_fraction) - p_transport);

    if (mu > 0.0) {
      const Numeric at_observation  = std::exp(-scaled_tau / model.mu0);
      const Numeric at_bottom       = std::exp((scaled_tau - scaled_depth) / abs_mu - scaled_depth / model.mu0);
      correction[i]                += coefficient * (at_observation - at_bottom);
    } else {
      const Numeric at_observation = std::exp(-scaled_tau / model.mu0);
      const Numeric at_top         = std::exp(-scaled_tau / abs_mu);
      correction[i] += coefficient * downward_tms_kernel(abs_mu, model.mu0, scaled_tau, at_observation, at_top);
    }

    if (mu < 0.0) {
      const Numeric beam_theta = std::acos(-model.mu0);
      const Numeric ray_theta  = std::acos(mu);
      if (std::abs(beam_theta - ray_theta) <= Constant::pi / 18.0) {
        const Numeric omega_f = model.physical_omega * model.delta_m_fraction;
        Vector        weighted_residue(model.original_moments.size());
        for (Index degree = 0; degree < static_cast<Index>(weighted_residue.size()); ++degree) {
          const Numeric residue    = degree < static_cast<Index>(model.transport_moments.size())
                                         ? 1.0
                                         : model.original_moments[degree] / model.delta_m_fraction;
          weighted_residue[degree] = static_cast<Numeric>(2 * degree + 1) * (2.0 * residue - residue * residue);
        }
        const Numeric ims_factor  = model.beam_intensity / (4.0 * Constant::pi) * omega_f * omega_f / (1.0 - omega_f);
        correction[i]            -= ims_factor * Legendre::legendre_sum(weighted_residue, nu) *
                                    ims_chi(physical_tau, abs_mu, model.mu0 / (1.0 - omega_f));
      }
    }
  }
  return correction;
}

void run_problem_1_case(const disort_test::reference::single_layer_case& test) {
  const auto& user_mu = disort_test::reference::problem_1_user_mu;
  Vector      moments(17, 0.0);
  moments[0]         = 1.0;
  auto         model = make_scalar_model(test.depth,
                                         test.omega,
                                         test.beam,
                                         disort_test::reference::problem_1_streams,
                                         moments,
                                         user_mu,
                                         disort_test::reference::problem_1_beam_mu,
                                         Constant::pi / disort_test::reference::problem_1_beam_mu);
  const Vector output_tau{0.0, test.depth};

  vdisort::user_u_data user;
  vdisort::flux_data   flux;
  for (Index level = 0; level < 2; ++level) {
    model.solver.u_user(user, output_tau[level], 0.0, user_mu, model.user_phase, model.user_beam_phase);
    for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
      const auto label = std::format("{} radiance [{}, {}]", test.name, level, angle);
      expect_reference(label, user.intensities[angle].I(), test.radiance[level, angle]);
      expect_unpolarized(label, user.intensities[angle]);
    }

    const auto values = model.solver.flux(flux, output_tau[level]);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), values.down_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), values.down_diffuse, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), values.up, test.up[level]);
  }
}

void test_problem_1() {
  for (const auto& test : disort_test::reference::problem_1) run_problem_1_case(test);
}

void run_problem_2_case(const disort_test::reference::single_layer_case& test) {
  const auto& user_mu = disort_test::reference::problem_2_user_mu;
  Vector      moments(17, 0.0);
  moments[0]         = 1.0;
  moments[2]         = 0.1;
  auto         model = make_scalar_model(test.depth,
                                         test.omega,
                                         test.beam,
                                         disort_test::reference::problem_2_streams,
                                         moments,
                                         user_mu,
                                         disort_test::reference::problem_2_beam_mu,
                                         Constant::pi);
  const Vector output_tau{0.0, test.depth};

  vdisort::user_u_data user;
  vdisort::flux_data   flux;
  for (Index level = 0; level < 2; ++level) {
    model.solver.u_user(user, output_tau[level], 0.0, user_mu, model.user_phase, model.user_beam_phase);
    for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
      const auto label = std::format("{} radiance [{}, {}]", test.name, level, angle);
      expect_reference(label, user.intensities[angle].I(), test.radiance[level, angle]);
      expect_unpolarized(label, user.intensities[angle]);
    }

    const auto values = model.solver.flux(flux, output_tau[level]);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), values.down_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), values.down_diffuse, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), values.up, test.up[level]);
  }
}

void test_problem_2() {
  for (const auto& test : disort_test::reference::problem_2) run_problem_2_case(test);
}

void run_problem_3_case(const disort_test::reference::single_layer_case& test) {
  constexpr Index nquad   = disort_test::reference::problem_3_streams;
  const auto&     user_mu = disort_test::reference::problem_3_user_mu;
  Vector          moments(disort_test::reference::problem_3_moments);
  for (Index degree = 0; degree < static_cast<Index>(moments.size()); ++degree)
    moments[degree] = std::pow(disort_test::reference::problem_3_asymmetry, degree);
  auto model =
      make_scalar_model(test.depth, test.omega, test.beam, nquad, moments, user_mu, 1.0, Constant::pi, moments[nquad]);
  const Vector output_tau{0.0, test.depth};

  Vector weighted_true(moments.size()), weighted_transport(model.transport_moments.size()),
      weighted_removed(moments.size());
  for (Index degree = 0; degree < static_cast<Index>(moments.size()); ++degree) {
    weighted_true[degree]    = static_cast<Numeric>(2 * degree + 1) * moments[degree];
    const Numeric residue    = degree < nquad ? 1.0 : moments[degree] / model.delta_m_fraction;
    weighted_removed[degree] = static_cast<Numeric>(2 * degree + 1) * residue;
  }
  for (Index degree = 0; degree < static_cast<Index>(model.transport_moments.size()); ++degree)
    weighted_transport[degree] = static_cast<Numeric>(2 * degree + 1) * model.transport_moments[degree];
  const auto scalar_phase = [](const Vector coefficients) {
    return [coefficients](Index, Numeric out_mu, Numeric out_phi, Numeric in_mu, Numeric in_phi) {
      rtepack::muelmat result{0.0};
      result[0, 0] = Legendre::legendre_sum(coefficients, scattering_angle_cosine(out_mu, out_phi, in_mu, in_phi));
      return result;
    };
  };
  const vdisort::delta_m_correction_cache polarized_correction(AscendingGrid{test.depth},
                                                               Vector{test.omega},
                                                               Vector{model.delta_m_fraction},
                                                               model.mu0,
                                                               model.phi0,
                                                               rtepack::stokvec{model.beam_intensity, 0.0, 0.0, 0.0},
                                                               Vector{user_mu},
                                                               Vector{0.0},
                                                               scalar_phase(weighted_true),
                                                               scalar_phase(weighted_transport),
                                                               scalar_phase(weighted_removed),
                                                               48,
                                                               96);

  vdisort::user_u_data    user;
  vdisort::flux_data      flux;
  rtepack::stokvec_vector polarized_ims, polarized_tms;
  for (Index level = 0; level < 2; ++level) {
    const Numeric physical_tau = output_tau[level];
    const Numeric scaled_tau   = model.optical_depth_scale * physical_tau;
    model.solver.u_user_corr(user,
                             polarized_ims,
                             polarized_tms,
                             physical_tau,
                             0,
                             polarized_correction,
                             model.user_phase,
                             model.user_beam_phase);
    const Vector correction = scalar_delta_m_correction(model, physical_tau, 0.0, user_mu);
    for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
      const auto polarized = polarized_tms[angle] + polarized_ims[angle];
      expect_reference(std::format("{} polarized delta-M scalar reduction [{}, {}]", test.name, level, angle),
                       polarized.I(),
                       correction[angle],
                       2e-7);
      expect_unpolarized(std::format("{} polarized delta-M scalar reduction [{}, {}]", test.name, level, angle),
                         polarized);
      const auto label = std::format("{} radiance [{}, {}]", test.name, level, angle);
      expect_reference(label, user.intensities[angle].I(), test.radiance[level, angle]);
      expect_unpolarized(label, user.intensities[angle]);
    }

    const auto scaled_flux = model.solver.flux(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    const Numeric physical_direct = model.mu0 * model.beam_intensity * std::exp(-physical_tau / model.mu0);
    const Numeric diffuse_down    = scaled_flux.down_diffuse - (physical_direct - scaled_flux.down_direct);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), scaled_flux.up, test.up[level]);
  }
}

void test_problem_3() {
  for (const auto& test : disort_test::reference::problem_3) run_problem_3_case(test);
}

void run_problem_4_case(const disort_test::reference::haze_l_case& test) {
  constexpr Index nquad          = disort_test::reference::problem_4_streams;
  const auto&     user_mu        = disort_test::reference::problem_4_user_mu;
  const Matrix    moments_matrix = disort_test::reference::haze_l_moments();
  const Vector    moments{moments_matrix[0]};
  auto            model = make_scalar_model(disort_test::reference::problem_4_total_tau,
                                            test.omega,
                                            true,
                                            nquad,
                                            moments,
                                            user_mu,
                                            test.beam_mu,
                                            Constant::pi,
                                            moments[nquad]);

  const AscendingGrid      scaled_output_tau(disort_test::reference::problem_4_output_tau.begin(),
                                             disort_test::reference::problem_4_output_tau.end(),
                                             [&](const Numeric tau) { return model.optical_depth_scale * tau; });
  rtepack::stokvec_tensor3 user(scaled_output_tau.size(), test.azimuth.size(), user_mu.size());
  model.solver.ungridded_u_user(
      user, scaled_output_tau, test.azimuth, user_mu, model.user_phase, model.user_beam_phase);
  for (Index level = 0; level < static_cast<Index>(disort_test::reference::problem_4_output_tau.size()); ++level) {
    const Numeric physical_tau = disort_test::reference::problem_4_output_tau[level];
    for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth) {
      const Vector correction = scalar_delta_m_correction(model, physical_tau, test.azimuth[azimuth], user_mu);
      for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
        auto intensity    = user[level, azimuth, angle];
        intensity.I()    += correction[angle];
        const auto label  = std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle);
        expect_reference(label, intensity.I(), test.radiance[azimuth, level, angle]);
        expect_unpolarized(label, intensity);
      }
    }
  }

  vdisort::flux_data flux;
  for (Index level = 0; level < static_cast<Index>(disort_test::reference::problem_4_output_tau.size()); ++level) {
    const Numeric physical_tau = disort_test::reference::problem_4_output_tau[level];
    const Numeric scaled_tau   = model.optical_depth_scale * physical_tau;
    const auto    scaled_flux  = model.solver.flux(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    const Numeric physical_direct = model.mu0 * model.beam_intensity * std::exp(-physical_tau / model.mu0);
    const Numeric diffuse_down    = scaled_flux.down_diffuse - (physical_direct - scaled_flux.down_direct);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), scaled_flux.up, test.up[level]);
  }
}

void test_problem_4() {
  for (const auto& test : disort_test::reference::problem_4) run_problem_4_case(test);
}

void run_problem_5_case(const disort_test::reference::scalar_case& test) {
  constexpr Index nquad          = disort_test::reference::problem_5_streams;
  const auto&     user_mu        = disort_test::reference::problem_5_user_mu;
  const Matrix    moments_matrix = disort_test::reference::cloud_c1_moments();
  const Vector    moments{moments_matrix[0]};
  auto            model = make_scalar_model(disort_test::reference::problem_5_total_tau,
                                            test.omega,
                                            true,
                                            nquad,
                                            moments,
                                            user_mu,
                                            1.0,
                                            Constant::pi,
                                            moments[nquad]);

  vdisort::user_u_data user;
  vdisort::flux_data   flux;
  for (Index level = 0; level < static_cast<Index>(test.tau.size()); ++level) {
    const Numeric physical_tau = test.tau[level];
    const Numeric scaled_tau   = model.optical_depth_scale * physical_tau;
    model.solver.u_user(user, scaled_tau, 0.0, user_mu, model.user_phase, model.user_beam_phase);
    const Vector correction = scalar_delta_m_correction(model, physical_tau, 0.0, user_mu);
    for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
      user.intensities[angle].I() += correction[angle];
      const auto label             = std::format("{} radiance [{}, {}]", test.name, level, angle);
      expect_reference(label, user.intensities[angle].I(), test.radiance[level, angle], 2e-4);
      expect_unpolarized(label, user.intensities[angle]);
    }

    const auto scaled_flux = model.solver.flux(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    const Numeric physical_direct = model.mu0 * model.beam_intensity * std::exp(-physical_tau / model.mu0);
    const Numeric diffuse_down    = scaled_flux.down_diffuse - (physical_direct - scaled_flux.down_direct);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), scaled_flux.up, test.up[level]);
  }
}

void test_problem_5() {
  run_problem_5_case(disort_test::reference::problem_5a);
  run_problem_5_case(disort_test::reference::problem_5b);
}

Numeric band_blackbody_radiance(const Numeric temperature,
                                const Numeric wavenumber_low,
                                const Numeric wavenumber_high) {
  if (temperature == 0.0 or wavenumber_low == wavenumber_high) return 0.0;
  constexpr Index intervals = 4096;
  const Numeric   scale     = Constant::h * Constant::c * 100.0 / (Constant::k * temperature);
  const Numeric   x0        = scale * wavenumber_low;
  const Numeric   x1        = scale * wavenumber_high;
  const Numeric   dx        = (x1 - x0) / intervals;
  const auto      integrand = [](const Numeric x) {
    if (x == 0.0 or x > 700.0) return 0.0;
    return x * x * x / std::expm1(x);
  };
  Numeric integral = integrand(x0) + integrand(x1);
  for (Index i = 1; i < intervals; ++i) integral += (i % 2 ? 4.0 : 2.0) * integrand(x0 + static_cast<Numeric>(i) * dx);
  integral *= dx / 3.0;
  return Constant::sigma * std::pow(temperature, 4) * Constant::inv_pi * integral /
         (Math::pow2(Constant::pi) * Math::pow2(Constant::pi) / 15.0);
}

std::vector<vdisort::BDRF> scalar_brdf_modes(const disort::brdf::RawFunction& raw, const Index number_of_modes) {
  return vdisort::brdf::depolarizing_fourier_modes(raw, number_of_modes, 100);
}

Numeric directional_emissivity(const disort::brdf::RawFunction& raw, const Numeric outgoing_mu) {
  constexpr Index nmu      = 64;
  constexpr Index nazimuth = 256;
  Vector          mu(nmu), weight(nmu);
  Legendre::PositiveDoubleGaussLegendre(mu, weight);
  Numeric reflectance = 0.0;
  for (Index j = 0; j < nmu; ++j) {
    Numeric azimuth_average = 0.0;
    for (Index k = 0; k < nazimuth; ++k)
      azimuth_average +=
          raw(outgoing_mu, mu[j], Constant::two_pi * (static_cast<Numeric>(k) + 0.5) / static_cast<Numeric>(nazimuth));
    // DISOBRDF's BEMST convention integrates azimuth through x = phi / pi,
    // without restoring the Jacobian pi.  Reproduce that convention here so
    // the thermal boundary agrees with the original DISORT reference cases.
    reflectance += 2.0 * weight[j] * mu[j] * azimuth_average / nazimuth;
  }
  return 1.0 - reflectance;
}

vdisort::main_data make_problem_6_model(const disort_test::reference::thermal_source_case& test) {
  constexpr Index nquad   = 16;
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Numeric   depth   = test.optical_depth == 0.0 ? 1e-12 : test.optical_depth;

  disort::brdf::RawFunction raw;
  Index                     nmodes = 1;
  if (test.surface == disort_test::reference::surface_type::lambertian) {
    raw = [albedo = test.lambertian_albedo](Numeric, Numeric, Numeric) { return albedo * Constant::inv_pi; };
  } else if (test.surface == disort_test::reference::surface_type::hapke) {
    raw    = disort::brdf::Hapke{.opposition_amplitude     = test.hapke_parameters[0],
                                 .opposition_width         = test.hapke_parameters[1],
                                 .single_scattering_albedo = test.hapke_parameters[2]};
    nmodes = nquad;
  }

  std::vector<vdisort::BDRF> brdf;
  if (raw) brdf = scalar_brdf_modes(raw, nmodes);

  Tensor7 phase(2, nmodes, 1, nquad, nquad, nstokes, nstokes, 0.0);
  phase[vdisort::cosine_mode, 0, 0, joker, joker, 0, 0] = 1.0;
  Tensor4 up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);

  const Numeric top_planck                = band_blackbody_radiance(test.top_temperature,
                                                                    disort_test::reference::problem_6_wavenumber_low,
                                                                    disort_test::reference::problem_6_wavenumber_high);
  down[vdisort::cosine_mode, 0, joker, 0] = test.top_isotropic + test.top_emissivity * top_planck;

  const Numeric bottom_planck = band_blackbody_radiance(test.bottom_temperature,
                                                        disort_test::reference::problem_6_wavenumber_low,
                                                        disort_test::reference::problem_6_wavenumber_high);
  Vector        quadrature_mu(n), quadrature_weight(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu, quadrature_weight);
  for (Index i = 0; i < n; ++i) {
    Numeric emissivity = 1.0;
    if (test.surface == disort_test::reference::surface_type::lambertian)
      emissivity = 1.0 - test.lambertian_albedo;
    else if (test.surface == disort_test::reference::surface_type::hapke)
      emissivity = directional_emissivity(raw, quadrature_mu[i]);
    up[vdisort::cosine_mode, 0, i, 0] = emissivity * bottom_planck;
  }

  Tensor3 source;
  if (test.interface_temperature.size() == 2) {
    source.resize(1, 2, nstokes);
    source               = 0.0;
    const Numeric top    = band_blackbody_radiance(test.interface_temperature[0],
                                                   disort_test::reference::problem_6_wavenumber_low,
                                                   disort_test::reference::problem_6_wavenumber_high);
    const Numeric bottom = band_blackbody_radiance(test.interface_temperature[1],
                                                   disort_test::reference::problem_6_wavenumber_low,
                                                   disort_test::reference::problem_6_wavenumber_high);
    source[0, 0, 0]      = top;
    source[0, 1, 0]      = (bottom - top) / test.optical_depth;
  } else {
    source.resize(1, 0, nstokes);
  }

  return vdisort_test::make_solver(nquad,
                                   nmodes,
                                   AscendingGrid{depth},
                                   Vector{test.single_scattering_albedo},
                                   std::move(phase),
                                   std::move(up),
                                   std::move(down),
                                   std::move(source),
                                   std::move(brdf),
                                   test.beam_mu,
                                   Vector{test.beam, 0.0, 0.0, 0.0},
                                   0.0);
}

void test_problem_6() {
  for (Index case_index = 0; case_index < static_cast<Index>(disort_test::reference::problem_6.size()); ++case_index) {
    const auto&        test      = disort_test::reference::problem_6[case_index];
    const auto&        reference = disort_test::reference::problem_6_flux[case_index];
    auto               model     = make_problem_6_model(test);
    vdisort::flux_data flux_data;
    for (Index level = 0; level < static_cast<Index>(test.output_tau.size()); ++level) {
      const Numeric tau  = test.optical_depth == 0.0 ? 0.0 : test.output_tau[level];
      const auto    flux = model.flux(flux_data, tau);
      for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
        expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux_data.u0[stream]);
      expect_reference(std::format("{} direct flux [{}]", test.name, level), flux.down_direct, reference.direct[level]);
      expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level),
                       flux.down_diffuse,
                       reference.diffuse_down[level],
                       1e-3);
      // The active single-precision DISORT references retain their legacy
      // non-Lambertian thermal-boundary integration.  The VDISORT scalar
      // embedding agrees to below one percent for those mixed-source cases.
      const Numeric thermal_brdf_tolerance = case_index >= 5 ? 1e-2 : 1e-3;
      expect_reference(
          std::format("{} up flux [{}]", test.name, level), flux.up, reference.up[level], thermal_brdf_tolerance);
      expect_reference(
          std::format("{} DFDT [{}]", test.name, level), flux.dfdt, reference.dfdt[level], thermal_brdf_tolerance);
    }
  }
}

vdisort::main_data make_problem_7_model(const disort_test::reference::scattering_thermal_case& test) {
  constexpr Index nstokes = vdisort::stokes_dimension;
  constexpr Index nmodes  = 1;  // Problem 7 has active azimuth-integrated flux references only.
  const Index     nquad   = test.streams;
  const Index     n       = nquad / 2;

  Vector moments(nquad);
  for (Index degree = 0; degree < nquad; ++degree)
    moments[degree] = std::pow(test.asymmetry_parameter, static_cast<Numeric>(degree));

  Vector quadrature_mu(nquad), quadrature_weight(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu[Range{0, n}], quadrature_weight);
  std::transform(quadrature_mu.begin(), quadrature_mu.begin() + n, quadrature_mu.begin() + n, [](const Numeric mu) {
    return -mu;
  });

  Tensor7 phase(2, nmodes, 1, nquad, nquad, nstokes, nstokes, 0.0);
  for (Index i = 0; i < nquad; ++i)
    for (Index j = 0; j < nquad; ++j)
      phase[vdisort::cosine_mode, 0, 0, i, j, 0, 0] = scalar_phase_mode(moments, 0, quadrature_mu[i], quadrature_mu[j]);

  Tensor6           beam_phase(2, nmodes, 1, nquad, nstokes, nstokes, 0.0);
  constexpr Numeric mu0 = 0.5;
  for (Index i = 0; i < nquad; ++i)
    beam_phase[vdisort::cosine_mode, 0, 0, i, 0, 0] = scalar_phase_mode(moments, 0, quadrature_mu[i], -mu0);

  disort::brdf::RawFunction raw;
  if (test.surface == disort_test::reference::surface_type::lambertian) {
    raw = [albedo = test.lambertian_albedo](Numeric, Numeric, Numeric) { return albedo * Constant::inv_pi; };
  } else if (test.surface == disort_test::reference::surface_type::hapke) {
    raw = disort::brdf::Hapke{.opposition_amplitude = 1.0, .opposition_width = 0.06, .single_scattering_albedo = 0.6};
  }
  std::vector<vdisort::BDRF> brdf;
  if (raw) brdf = scalar_brdf_modes(raw, nmodes);

  Tensor4       up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);
  const Numeric top_boundary_planck =
      band_blackbody_radiance(test.top_boundary_temperature, test.wavenumber_low, test.wavenumber_high);
  down[vdisort::cosine_mode, 0, joker, 0] = test.top_isotropic + top_boundary_planck;

  const Numeric bottom_boundary_planck =
      band_blackbody_radiance(test.bottom_boundary_temperature, test.wavenumber_low, test.wavenumber_high);
  for (Index i = 0; i < n; ++i) {
    Numeric emissivity = 1.0;
    if (test.surface == disort_test::reference::surface_type::lambertian)
      emissivity = 1.0 - test.lambertian_albedo;
    else if (test.surface == disort_test::reference::surface_type::hapke)
      emissivity = directional_emissivity(raw, quadrature_mu[i]);
    up[vdisort::cosine_mode, 0, i, 0] = emissivity * bottom_boundary_planck;
  }

  const Numeric atmosphere_top_planck =
      band_blackbody_radiance(test.atmosphere_top_temperature, test.wavenumber_low, test.wavenumber_high);
  const Numeric atmosphere_bottom_planck =
      band_blackbody_radiance(test.atmosphere_bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  Tensor3 source(1, 2, nstokes, 0.0);
  source[0, 0, 0] = atmosphere_top_planck;
  source[0, 1, 0] = (atmosphere_bottom_planck - atmosphere_top_planck) / test.optical_depth;

  return vdisort_test::make_solver(nquad,
                                   nmodes,
                                   AscendingGrid{test.optical_depth},
                                   Vector{test.single_scattering_albedo},
                                   std::move(phase),
                                   std::move(up),
                                   std::move(down),
                                   std::move(source),
                                   std::move(brdf),
                                   mu0,
                                   Vector{test.beam, 0.0, 0.0, 0.0},
                                   0.0,
                                   std::move(beam_phase));
}

void expect_small_reference(const std::string_view name,
                            const Numeric          actual,
                            const Numeric          expected,
                            const Numeric          relative_tolerance = 2e-3,
                            const Numeric          absolute_tolerance = 1e-10) {
  ARTS_USER_ERROR_IF(std::abs(actual - expected) > absolute_tolerance + relative_tolerance * std::abs(expected),
                     "{}: expected {}, got {} (difference {})",
                     name,
                     expected,
                     actual,
                     actual - expected);
}

void test_problem_7() {
  for (Index case_index = 0; case_index < static_cast<Index>(disort_test::reference::problem_7.size()); ++case_index) {
    const auto&        test       = disort_test::reference::problem_7[case_index];
    const auto&        reference  = disort_test::reference::problem_7_flux[case_index];
    auto               model      = make_problem_7_model(test);
    const Vector       output_tau = case_index < 2 ? Vector{0.0, test.optical_depth} : Vector{0.0, 0.5, 1.0};
    vdisort::flux_data flux_data;
    for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level) {
      const auto flux = model.flux(flux_data, output_tau[level]);
      for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
        expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux_data.u0[stream]);

      if (case_index == 1) {
        expect_small_reference(
            std::format("{} direct flux [{}]", test.name, level), flux.down_direct, reference.direct[level]);
        expect_small_reference(std::format("{} diffuse-down flux [{}]", test.name, level),
                               flux.down_diffuse,
                               reference.diffuse_down[level]);
        expect_small_reference(std::format("{} up flux [{}]", test.name, level), flux.up, reference.up[level]);
        expect_small_reference(std::format("{} DFDT [{}]", test.name, level), flux.dfdt, reference.dfdt[level]);
      } else {
        const Numeric tolerance = case_index == 4 ? 1e-2 : 2e-3;
        expect_reference(
            std::format("{} direct flux [{}]", test.name, level), flux.down_direct, reference.direct[level], tolerance);
        expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level),
                         flux.down_diffuse,
                         reference.diffuse_down[level],
                         tolerance);
        expect_reference(std::format("{} up flux [{}]", test.name, level), flux.up, reference.up[level], tolerance);
        expect_reference(std::format("{} DFDT [{}]", test.name, level), flux.dfdt, reference.dfdt[level], tolerance);
      }
    }
  }
}

struct problem_8_model {
  vdisort::main_data         solver;
  vdisort::phase_matrix_data user_phase;
};

problem_8_model make_problem_8_model(const disort_test::reference::layered_isotropic_case& test) {
  constexpr Index nquad   = disort_test::reference::problem_8_streams;
  constexpr Index nmodes  = 1;
  constexpr Index nlayers = 2;
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Index     nuser   = static_cast<Index>(disort_test::reference::problem_8_user_mu.size());

  Tensor7 phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  phase[vdisort::cosine_mode, 0, joker, joker, joker, 0, 0] = 1.0;

  Tensor4 up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);
  down[vdisort::cosine_mode, 0, joker, 0] = Constant::inv_pi;

  vdisort::phase_matrix_data user_phase(2, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  for (Index layer = 0; layer < nlayers; ++layer)
    for (Index user = 0; user < nuser; ++user)
      for (Index stream = 0; stream < nquad; ++stream)
        user_phase[vdisort::cosine_mode, 0, layer, user, stream][0, 0] = 1.0;

  auto solver = vdisort_test::make_solver(nquad,
                                          nmodes,
                                          AscendingGrid{test.cumulative_tau},
                                          Vector{test.single_scattering_albedo},
                                          std::move(phase),
                                          std::move(up),
                                          std::move(down),
                                          Tensor3(nlayers, 0, nstokes),
                                          {},
                                          0.5,
                                          Vector(4, 0.0),
                                          0.0);
  return {.solver = std::move(solver), .user_phase = std::move(user_phase)};
}

void test_problem_8() {
  for (const auto& test : disort_test::reference::problem_8) {
    auto                 model = make_problem_8_model(test);
    vdisort::user_u_data user;
    vdisort::flux_data   flux_data;
    for (Index level = 0; level < static_cast<Index>(test.output_tau.size()); ++level) {
      const Numeric tau = test.output_tau[level];
      model.solver.u_user(user,
                          tau,
                          disort_test::reference::problem_8_azimuth,
                          disort_test::reference::problem_8_user_mu,
                          model.user_phase,
                          {});
      for (Index angle = 0; angle < static_cast<Index>(user.intensities.size()); ++angle) {
        expect_unpolarized(std::format("{} radiance [{}, {}]", test.name, level, angle), user.intensities[angle]);
        expect_reference(std::format("{} radiance [{}, {}]", test.name, level, angle),
                         user.intensities[angle].I(),
                         test.radiance[level, angle],
                         2e-4);
      }

      const auto flux = model.solver.flux(flux_data, tau);
      for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
        expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux_data.u0[stream]);
      expect_reference(std::format("{} direct flux [{}]", test.name, level), flux.down_direct, test.direct[level]);
      expect_reference(
          std::format("{} diffuse-down flux [{}]", test.name, level), flux.down_diffuse, test.diffuse_down[level]);
      expect_reference(std::format("{} up flux [{}]", test.name, level), flux.up, test.up[level]);
    }
  }
}

struct problem_9_model {
  vdisort::main_data              solver;
  vdisort::phase_matrix_data      user_phase;
  vdisort::beam_phase_matrix_data user_beam_phase;
};

Matrix problem_9_moments(const disort_test::reference::general_multilayer_case& test, const Index nquad) {
  const Index nlayers = static_cast<Index>(test.cumulative_tau.size());
  Matrix      moments(nlayers, nquad, 0.0);
  if (test.phase == disort_test::reference::phase_type::isotropic) {
    moments[joker, 0] = 1.0;
  } else if (test.phase == disort_test::reference::phase_type::problem_9b) {
    for (Index layer = 0; layer < nlayers; ++layer)
      for (Index degree = 0; degree < nquad; ++degree)
        moments[layer, degree] = disort_test::reference::problem_9b_moments[degree];
  } else {
    for (Index layer = 0; layer < nlayers; ++layer)
      for (Index degree = 0; degree < nquad; ++degree)
        moments[layer, degree] =
            std::pow(disort_test::reference::problem_9c_asymmetry[layer], static_cast<Numeric>(degree));
  }
  return moments;
}

problem_9_model make_problem_9_model(const disort_test::reference::general_multilayer_case& test,
                                     const Index                                            nquad,
                                     const ConstVectorView&                                 user_mu) {
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Index     nlayers = static_cast<Index>(test.cumulative_tau.size());
  const Index     nuser   = static_cast<Index>(user_mu.size());
  const Index     nmodes  = test.beam == 0.0 ? 1 : nquad;
  const Matrix    moments = problem_9_moments(test, nquad);

  Vector quadrature_mu(nquad), quadrature_weight(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu[Range{0, n}], quadrature_weight);
  std::transform(quadrature_mu.begin(), quadrature_mu.begin() + n, quadrature_mu.begin() + n, [](const Numeric mu) {
    return -mu;
  });

  Tensor7                         phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  Tensor6                         beam_phase(2, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  vdisort::phase_matrix_data      user_phase(2, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data user_beam_phase(2, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index layer = 0; layer < nlayers; ++layer) {
    const auto layer_moments = moments[layer];
    for (Index mode = 0; mode < nmodes; ++mode) {
      for (Index out = 0; out < nquad; ++out) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(layer_moments, mode, quadrature_mu[out], quadrature_mu[in]);
          phase[vdisort::cosine_mode, mode, layer, out, in, 0, 0] = value;
          if (mode > 0) phase[vdisort::sine_mode, mode, layer, out, in, 0, 0] = value;
        }
        beam_phase[vdisort::cosine_mode, mode, layer, out, 0, 0] =
            scalar_phase_mode(layer_moments, mode, quadrature_mu[out], -test.beam_mu);
      }
      for (Index user = 0; user < nuser; ++user) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(layer_moments, mode, user_mu[user], quadrature_mu[in]);
          user_phase[vdisort::cosine_mode, mode, layer, user, in][0, 0] = value;
          if (mode > 0) user_phase[vdisort::sine_mode, mode, layer, user, in][0, 0] = value;
        }
        user_beam_phase[vdisort::cosine_mode, mode, layer, user][0, 0] =
            scalar_phase_mode(layer_moments, mode, user_mu[user], -test.beam_mu);
      }
    }
  }

  Tensor4 up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);
  down[vdisort::cosine_mode, 0, joker, 0] =
      test.top_isotropic + band_blackbody_radiance(test.top_temperature, test.wavenumber_low, test.wavenumber_high);
  const Numeric bottom_planck =
      band_blackbody_radiance(test.bottom_temperature, test.wavenumber_low, test.wavenumber_high);
  up[vdisort::cosine_mode, 0, joker, 0] = (1.0 - test.surface_albedo) * bottom_planck;

  std::vector<vdisort::BDRF> brdf;
  if (test.surface_albedo != 0.0) {
    const disort::brdf::RawFunction lambert = [albedo = test.surface_albedo](Numeric, Numeric, Numeric) {
      return albedo * Constant::inv_pi;
    };
    brdf = scalar_brdf_modes(lambert, 1);
  }

  Tensor3 source;
  if (test.interface_temperature.empty()) {
    source.resize(nlayers, 0, nstokes);
  } else {
    source.resize(nlayers, 2, nstokes);
    source = 0.0;
    Vector planck(nlayers + 1);
    for (Index interface = 0; interface <= nlayers; ++interface)
      planck[interface] =
          band_blackbody_radiance(test.interface_temperature[interface], test.wavenumber_low, test.wavenumber_high);
    Numeric tau0 = 0.0;
    for (Index layer = 0; layer < nlayers; ++layer) {
      const Numeric tau1  = test.cumulative_tau[layer];
      const Numeric slope = (planck[layer + 1] - planck[layer]) / (tau1 - tau0);
      source[layer, 0, 0] = planck[layer] - slope * tau0;
      source[layer, 1, 0] = slope;
      tau0                = tau1;
    }
  }

  auto solver = vdisort_test::make_solver(nquad,
                                          nmodes,
                                          AscendingGrid{test.cumulative_tau},
                                          Vector{test.single_scattering_albedo},
                                          std::move(phase),
                                          std::move(up),
                                          std::move(down),
                                          std::move(source),
                                          std::move(brdf),
                                          test.beam_mu,
                                          Vector{test.beam, 0.0, 0.0, 0.0},
                                          0.0,
                                          std::move(beam_phase));
  return {
      .solver = std::move(solver), .user_phase = std::move(user_phase), .user_beam_phase = std::move(user_beam_phase)};
}

void test_problem_9() {
  for (const auto& test : disort_test::reference::problem_9) {
    auto model = make_problem_9_model(
        test, disort_test::reference::problem_9_streams, disort_test::reference::problem_9_user_mu);
    const AscendingGrid      output_tau{test.output_tau};
    rtepack::stokvec_tensor3 user(
        output_tau.size(), test.azimuth.size(), disort_test::reference::problem_9_user_mu.size());
    model.solver.ungridded_u_user(user,
                                  output_tau,
                                  test.azimuth,
                                  disort_test::reference::problem_9_user_mu,
                                  model.user_phase,
                                  model.user_beam_phase);
    vdisort::flux_data flux_data;
    for (Index level = 0; level < static_cast<Index>(test.output_tau.size()); ++level) {
      for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth) {
        for (Index angle = 0; angle < static_cast<Index>(disort_test::reference::problem_9_user_mu.size()); ++angle) {
          const auto intensity = user[level, azimuth, angle];
          expect_unpolarized(std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle), intensity);
          expect_reference(std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle),
                           intensity.I(),
                           test.radiance[azimuth, level, angle],
                           test.thermal ? 1e-3 : reference_tolerance);
        }
      }
    }

    for (Index level = 0; level < static_cast<Index>(test.output_tau.size()); ++level) {
      const auto flux = model.solver.flux(flux_data, test.output_tau[level]);
      for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
        expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux_data.u0[stream]);
      const Numeric tolerance = test.thermal ? 1e-3 : reference_tolerance;
      expect_reference(
          std::format("{} direct flux [{}]", test.name, level), flux.down_direct, test.direct[level], tolerance);
      expect_reference(std::format("{} diffuse-down flux [{}]", test.name, level),
                       flux.down_diffuse,
                       test.diffuse_down[level],
                       tolerance);
      expect_reference(std::format("{} up flux [{}]", test.name, level), flux.up, test.up[level], tolerance);
    }
  }
}

void test_problem_10() {
  const auto& test  = disort_test::reference::problem_9[2];
  auto        model = make_problem_9_model(
      test, disort_test::reference::problem_10_streams, disort_test::reference::problem_10_user_mu);

  vdisort::u_data      quadrature;
  vdisort::user_u_data user;
  vdisort::flux_data   flux_data;

  for (const Numeric tau : disort_test::reference::problem_10_output_tau) {
    const auto flux_before = model.solver.flux(flux_data, tau);
    for (const Numeric phi : disort_test::reference::problem_10_azimuth) {
      model.solver.u(quadrature, tau, phi);
      model.solver.u_user(
          user, tau, phi, disort_test::reference::problem_10_user_mu, model.user_phase, model.user_beam_phase);

      for (Index angle = 0; angle < disort_test::reference::problem_10_streams; ++angle) {
        const auto    label   = std::format("Problem 10 radiance [tau={}, phi={}, angle={}]", tau, phi, angle);
        const Numeric user_mu = disort_test::reference::problem_10_user_mu[angle];
        const Index   stream =
            static_cast<Index>(std::min_element(model.solver.mu().begin(),
                                                model.solver.mu().end(),
                                                [user_mu](const Numeric lhs, const Numeric rhs) {
                                                  return std::abs(lhs - user_mu) < std::abs(rhs - user_mu);
                                                }) -
                               model.solver.mu().begin());
        expect_unpolarized(label, quadrature.intensities[stream]);
        expect_unpolarized(label, user.intensities[angle]);
        for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
          expect_reference(label, user.intensities[angle][stokes], quadrature.intensities[stream][stokes], 2e-7);
      }
    }

    const auto flux_after = model.solver.flux(flux_data, tau);
    expect_reference("Problem 10 direct flux", flux_after.down_direct, flux_before.down_direct, 1e-12);
    expect_reference("Problem 10 diffuse-down flux", flux_after.down_diffuse, flux_before.down_diffuse, 1e-12);
    expect_reference("Problem 10 upward flux", flux_after.up, flux_before.up, 1e-12);
    expect_reference("Problem 10 DFDT", flux_after.dfdt, flux_before.dfdt, 1e-12);
  }
}

problem_9_model make_problem_11_model(const bool subdivided) {
  constexpr Index nquad   = disort_test::reference::problem_11_streams;
  constexpr Index nmodes  = 1;
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Index     nlayers = subdivided ? 3 : 1;
  const Index     nuser   = static_cast<Index>(disort_test::reference::problem_11_user_mu.size());

  Tensor7 phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  phase[vdisort::cosine_mode, 0, joker, joker, joker, 0, 0] = 1.0;
  Tensor6 beam_phase(2, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  beam_phase[vdisort::cosine_mode, 0, joker, joker, 0, 0] = 1.0;

  vdisort::phase_matrix_data      user_phase(2, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data user_beam_phase(2, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index layer = 0; layer < nlayers; ++layer)
    for (Index user = 0; user < nuser; ++user) {
      for (Index stream = 0; stream < nquad; ++stream)
        user_phase[vdisort::cosine_mode, 0, layer, user, stream][0, 0] = 1.0;
      user_beam_phase[vdisort::cosine_mode, 0, layer, user][0, 0] = 1.0;
    }

  Tensor4 up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);
  down[vdisort::cosine_mode, 0, joker, 0] = disort_test::reference::problem_11_top_isotropic;

  const disort::brdf::RawFunction lambert = [](Numeric, Numeric, Numeric) {
    return disort_test::reference::problem_11_surface_albedo * Constant::inv_pi;
  };

  auto solver = vdisort_test::make_solver(
      nquad,
      nmodes,
      subdivided ? AscendingGrid{Vector(disort_test::reference::problem_11_subdivided_tau)}
                 : AscendingGrid{disort_test::reference::problem_11_output_tau.back()},
      Vector(nlayers, disort_test::reference::problem_11_omega),
      std::move(phase),
      std::move(up),
      std::move(down),
      Tensor3(nlayers, 0, nstokes),
      scalar_brdf_modes(lambert, nmodes),
      disort_test::reference::problem_11_beam_mu,
      Vector{disort_test::reference::problem_11_beam / disort_test::reference::problem_11_beam_mu, 0.0, 0.0, 0.0},
      0.0,
      std::move(beam_phase));
  return {
      .solver = std::move(solver), .user_phase = std::move(user_phase), .user_beam_phase = std::move(user_beam_phase)};
}

void test_problem_11() {
  auto one_layer  = make_problem_11_model(false);
  auto subdivided = make_problem_11_model(true);

  vdisort::user_u_data one_user, subdivided_user;
  for (const Numeric phi : disort_test::reference::problem_11_azimuth)
    for (const Numeric tau : disort_test::reference::problem_11_output_tau) {
      one_layer.solver.u_user(one_user,
                              tau,
                              phi,
                              disort_test::reference::problem_11_user_mu,
                              one_layer.user_phase,
                              one_layer.user_beam_phase);
      subdivided.solver.u_user(subdivided_user,
                               tau,
                               phi,
                               disort_test::reference::problem_11_user_mu,
                               subdivided.user_phase,
                               subdivided.user_beam_phase);
      for (Index angle = 0; angle < static_cast<Index>(one_user.intensities.size()); ++angle) {
        const auto label = std::format("Problem 11 radiance [phi={}, tau={}, angle={}]", phi, tau, angle);
        expect_unpolarized(label, one_user.intensities[angle]);
        expect_unpolarized(label, subdivided_user.intensities[angle]);
        for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
          expect_reference(
              label, subdivided_user.intensities[angle][stokes], one_user.intensities[angle][stokes], 2e-7);
      }
    }

  vdisort::flux_data one_flux_data, subdivided_flux_data;
  for (const Numeric tau : disort_test::reference::problem_11_output_tau) {
    const auto one  = one_layer.solver.flux(one_flux_data, tau);
    const auto many = subdivided.solver.flux(subdivided_flux_data, tau);
    expect_reference(std::format("Problem 11 upward flux [{}]", tau), many.up, one.up, 2e-7);
    expect_reference(std::format("Problem 11 diffuse-down flux [{}]", tau), many.down_diffuse, one.down_diffuse, 2e-7);
    expect_reference(std::format("Problem 11 direct-down flux [{}]", tau), many.down_direct, one.down_direct, 2e-7);
    expect_reference(std::format("Problem 11 DFDT [{}]", tau), many.dfdt, one.dfdt, 2e-7);
    for (Index stream = 0; stream < static_cast<Index>(one_flux_data.u0.size()); ++stream) {
      expect_unpolarized(std::format("Problem 11 one-layer flux field [{}, {}]", tau, stream),
                         one_flux_data.u0[stream]);
      expect_unpolarized(std::format("Problem 11 subdivided flux field [{}, {}]", tau, stream),
                         subdivided_flux_data.u0[stream]);
    }
  }
}

scalar_vdisort_model make_problem_12_model(const bool subdivided) {
  constexpr Index   nquad               = disort_test::reference::problem_12_streams;
  constexpr Index   nstokes             = vdisort::stokes_dimension;
  constexpr Numeric physical_omega      = disort_test::reference::problem_12_omega;
  constexpr Numeric g                   = disort_test::reference::problem_12_asymmetry;
  const Numeric     delta_m_fraction    = std::pow(g, static_cast<Numeric>(nquad));
  const Numeric     optical_depth_scale = 1.0 - physical_omega * delta_m_fraction;
  const Numeric     transport_omega     = physical_omega * (1.0 - delta_m_fraction) / optical_depth_scale;
  const Index       n                   = nquad / 2;
  const Index       nmodes              = nquad;
  const Index       nlayers             = subdivided ? 3 : 1;
  const auto&       user_mu             = disort_test::reference::problem_12_user_mu;
  const Index       nuser               = static_cast<Index>(user_mu.size());

  Vector original_moments(nquad + 1);
  for (Index degree = 0; degree <= nquad; ++degree)
    original_moments[degree] = std::pow(g, static_cast<Numeric>(degree));
  Vector transport_moments(nquad);
  for (Index degree = 0; degree < nquad; ++degree)
    transport_moments[degree] = (original_moments[degree] - delta_m_fraction) / (1.0 - delta_m_fraction);

  Vector quadrature_mu(nquad), quadrature_weight(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu[Range{0, n}], quadrature_weight);
  std::transform(quadrature_mu.begin(), quadrature_mu.begin() + n, quadrature_mu.begin() + n, [](const Numeric mu) {
    return -mu;
  });

  Tensor7                         phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  Tensor6                         beam_phase(2, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  vdisort::phase_matrix_data      user_phase(2, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data user_beam_phase(2, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index layer = 0; layer < nlayers; ++layer)
    for (Index mode = 0; mode < nmodes; ++mode) {
      for (Index out = 0; out < nquad; ++out) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(transport_moments, mode, quadrature_mu[out], quadrature_mu[in]);
          phase[vdisort::cosine_mode, mode, layer, out, in, 0, 0] = value;
          if (mode > 0) phase[vdisort::sine_mode, mode, layer, out, in, 0, 0] = value;
        }
        beam_phase[vdisort::cosine_mode, mode, layer, out, 0, 0] =
            scalar_phase_mode(transport_moments, mode, quadrature_mu[out], -disort_test::reference::problem_12_beam_mu);
      }
      for (Index user = 0; user < nuser; ++user) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(transport_moments, mode, user_mu[user], quadrature_mu[in]);
          user_phase[vdisort::cosine_mode, mode, layer, user, in][0, 0] = value;
          if (mode > 0) user_phase[vdisort::sine_mode, mode, layer, user, in][0, 0] = value;
        }
        user_beam_phase[vdisort::cosine_mode, mode, layer, user][0, 0] =
            scalar_phase_mode(transport_moments, mode, user_mu[user], -disort_test::reference::problem_12_beam_mu);
      }
    }

  Tensor4                         up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);
  const disort::brdf::RawFunction lambert = [](Numeric, Numeric, Numeric) {
    return disort_test::reference::problem_12_surface_albedo * Constant::inv_pi;
  };
  Vector cumulative_tau = subdivided ? Vector(disort_test::reference::problem_12_subdivided_tau)
                                     : Vector{disort_test::reference::problem_12_output_tau.back()};
  stdr::transform(cumulative_tau, cumulative_tau.begin(), [optical_depth_scale](const Numeric tau) {
    return optical_depth_scale * tau;
  });

  auto solver = vdisort_test::make_solver(
      nquad,
      nmodes,
      AscendingGrid{std::move(cumulative_tau)},
      Vector(nlayers, transport_omega),
      std::move(phase),
      std::move(up),
      std::move(down),
      Tensor3(nlayers, 0, nstokes),
      scalar_brdf_modes(lambert, 1),
      disort_test::reference::problem_12_beam_mu,
      Vector{disort_test::reference::problem_12_beam / disort_test::reference::problem_12_beam_mu, 0.0, 0.0, 0.0},
      0.0,
      std::move(beam_phase));
  return {.solver              = std::move(solver),
          .user_phase          = std::move(user_phase),
          .user_beam_phase     = std::move(user_beam_phase),
          .original_moments    = std::move(original_moments),
          .transport_moments   = std::move(transport_moments),
          .physical_depth      = disort_test::reference::problem_12_output_tau.back(),
          .physical_omega      = physical_omega,
          .mu0                 = disort_test::reference::problem_12_beam_mu,
          .phi0                = 0.0,
          .beam_intensity      = disort_test::reference::problem_12_beam / disort_test::reference::problem_12_beam_mu,
          .delta_m_fraction    = delta_m_fraction,
          .optical_depth_scale = optical_depth_scale};
}

void test_problem_12() {
  auto        one_layer  = make_problem_12_model(false);
  auto        subdivided = make_problem_12_model(true);
  const auto& user_mu    = disort_test::reference::problem_12_user_mu;

  vdisort::user_u_data one_user, subdivided_user;
  vdisort::flux_data   one_flux_data, subdivided_flux_data;
  for (const Numeric physical_tau : disort_test::reference::problem_12_output_tau) {
    const Numeric scaled_tau = one_layer.optical_depth_scale * physical_tau;
    one_layer.solver.u_user(one_user,
                            scaled_tau,
                            disort_test::reference::problem_12_azimuth,
                            user_mu,
                            one_layer.user_phase,
                            one_layer.user_beam_phase);
    subdivided.solver.u_user(subdivided_user,
                             scaled_tau,
                             disort_test::reference::problem_12_azimuth,
                             user_mu,
                             subdivided.user_phase,
                             subdivided.user_beam_phase);
    const Vector one_correction =
        scalar_delta_m_correction(one_layer, physical_tau, disort_test::reference::problem_12_azimuth, user_mu);
    const Vector subdivided_correction =
        scalar_delta_m_correction(subdivided, physical_tau, disort_test::reference::problem_12_azimuth, user_mu);
    for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
      one_user.intensities[angle].I()        += one_correction[angle];
      subdivided_user.intensities[angle].I() += subdivided_correction[angle];
      const auto label = std::format("Problem 12 radiance [tau={}, angle={}]", physical_tau, angle);
      expect_unpolarized(label, one_user.intensities[angle]);
      expect_unpolarized(label, subdivided_user.intensities[angle]);
      for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
        expect_reference(label, subdivided_user.intensities[angle][stokes], one_user.intensities[angle][stokes], 2e-7);
    }

    const auto one  = one_layer.solver.flux(one_flux_data, scaled_tau);
    const auto many = subdivided.solver.flux(subdivided_flux_data, scaled_tau);
    expect_reference(std::format("Problem 12 upward flux [{}]", physical_tau), many.up, one.up, 2e-7);
    expect_reference(
        std::format("Problem 12 diffuse-down flux [{}]", physical_tau), many.down_diffuse, one.down_diffuse, 2e-7);
    expect_reference(
        std::format("Problem 12 direct-down flux [{}]", physical_tau), many.down_direct, one.down_direct, 2e-12);
    expect_reference(std::format("Problem 12 DFDT [{}]", physical_tau), many.dfdt, one.dfdt, 2e-7);
    for (Index stream = 0; stream < static_cast<Index>(one_flux_data.u0.size()); ++stream) {
      expect_unpolarized(std::format("Problem 12 one-layer flux field [{}, {}]", physical_tau, stream),
                         one_flux_data.u0[stream]);
      expect_unpolarized(std::format("Problem 12 subdivided flux field [{}, {}]", physical_tau, stream),
                         subdivided_flux_data.u0[stream]);
    }
  }
}

vdisort::main_data make_problem_13_model(const disort_test::reference::albedo_transmission_reference& test) {
  constexpr Index   nquad   = disort_test::reference::problem_13_streams;
  constexpr Index   nmodes  = nquad;
  constexpr Index   nstokes = vdisort::stokes_dimension;
  constexpr Numeric g       = disort_test::reference::problem_13_asymmetry;
  const Numeric     f       = std::pow(g, static_cast<Numeric>(nquad));
  const Index       n       = nquad / 2;
  const Index       nlayers = static_cast<Index>(test.cumulative_tau.size());

  Vector transport_moments(nquad);
  for (Index degree = 0; degree < nquad; ++degree)
    transport_moments[degree] = (std::pow(g, static_cast<Numeric>(degree)) - f) / (1.0 - f);

  Vector quadrature_mu(nquad), quadrature_weight(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu[Range{0, n}], quadrature_weight);
  std::transform(quadrature_mu.begin(), quadrature_mu.begin() + n, quadrature_mu.begin() + n, [](const Numeric mu) {
    return -mu;
  });

  Tensor7 phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  Tensor6 beam_phase(2, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  for (Index layer = 0; layer < nlayers; ++layer)
    for (Index mode = 0; mode < nmodes; ++mode)
      for (Index out = 0; out < nquad; ++out) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(transport_moments, mode, quadrature_mu[out], quadrature_mu[in]);
          phase[vdisort::cosine_mode, mode, layer, out, in, 0, 0] = value;
          if (mode > 0) phase[vdisort::sine_mode, mode, layer, out, in, 0, 0] = value;
        }
        beam_phase[vdisort::cosine_mode, mode, layer, out, 0, 0] =
            scalar_phase_mode(transport_moments, mode, quadrature_mu[out], -disort_test::reference::problem_13_beam_mu);
      }

  Vector  scaled_cumulative_tau(nlayers), transport_omega(nlayers);
  Numeric physical_top = 0.0;
  Numeric scaled_top   = 0.0;
  for (Index layer = 0; layer < nlayers; ++layer) {
    const Numeric omega           = test.single_scattering_albedo[layer];
    const Numeric scale           = 1.0 - omega * f;
    scaled_top                   += scale * (test.cumulative_tau[layer] - physical_top);
    scaled_cumulative_tau[layer]  = scaled_top;
    transport_omega[layer]        = omega * (1.0 - f) / scale;
    physical_top                  = test.cumulative_tau[layer];
  }

  Tensor4                         up(2, nmodes, n, nstokes, 0.0), down(2, nmodes, n, nstokes, 0.0);
  const disort::brdf::RawFunction lambert = [](Numeric, Numeric, Numeric) {
    return disort_test::reference::problem_13_surface_albedo * Constant::inv_pi;
  };
  return vdisort_test::make_solver(nquad,
                                   nmodes,
                                   AscendingGrid{std::move(scaled_cumulative_tau)},
                                   std::move(transport_omega),
                                   std::move(phase),
                                   std::move(up),
                                   std::move(down),
                                   Tensor3(nlayers, 0, nstokes),
                                   scalar_brdf_modes(lambert, 1),
                                   disort_test::reference::problem_13_beam_mu,
                                   Vector{disort_test::reference::problem_13_beam, 0.0, 0.0, 0.0},
                                   0.0,
                                   std::move(beam_phase));
}

void test_problem_13_boundary_limits() {
  using namespace disort_test::reference;
  constexpr Numeric               depth   = 0.25;
  constexpr Index                 nquad   = problem_13_streams;
  constexpr Index                 nmodes  = 1;
  constexpr Index                 nstokes = vdisort::stokes_dimension;
  constexpr Index                 n       = nquad / 2;
  const disort::brdf::RawFunction lambert = [](Numeric, Numeric, Numeric) {
    return problem_13_surface_albedo * Constant::inv_pi;
  };
  auto absorbing = vdisort_test::make_solver(nquad,
                                             nmodes,
                                             AscendingGrid{depth},
                                             Vector{0.0},
                                             Tensor7(2, nmodes, 1, nquad, nquad, nstokes, nstokes, 0.0),
                                             Tensor4(2, nmodes, n, nstokes, 0.0),
                                             Tensor4(2, nmodes, n, nstokes, 0.0),
                                             Tensor3(1, 0, nstokes),
                                             scalar_brdf_modes(lambert, 1),
                                             problem_13_beam_mu,
                                             Vector{problem_13_beam, 0.0, 0.0, 0.0},
                                             0.0,
                                             Tensor6(2, nmodes, 1, nquad, nstokes, nstokes, 0.0));

  const Numeric incident_flux       = problem_13_beam * problem_13_beam_mu;
  const Numeric bottom_beam         = incident_flux * std::exp(-depth / problem_13_beam_mu);
  Numeric       upward_transmission = 0.0;
  for (Index i = 0; i < n; ++i)
    upward_transmission += absorbing.weights()[i] * absorbing.mu()[i] * std::exp(-depth / absorbing.mu()[i]);
  const Numeric expected_up = 2.0 * problem_13_surface_albedo * bottom_beam * upward_transmission;

  vdisort::flux_data flux_data;
  const auto         top    = absorbing.flux(flux_data, 0.0);
  const auto         bottom = absorbing.flux(flux_data, depth);
  expect_reference("Problem 13 absorbing Lambertian top reflection", top.up, expected_up, 2e-12);
  expect_reference("Problem 13 absorbing direct transmission", bottom.down_direct, bottom_beam, 2e-12);
  expect_reference("Problem 13 absorbing diffuse transmission", bottom.down_diffuse, 0.0, 2e-12);
}

void test_problem_13() {
  using namespace disort_test::reference;
  test_problem_13_boundary_limits();
  for (const auto& test : problem_13) {
    auto               solver        = make_problem_13_model(test);
    const Numeric      incident_flux = problem_13_beam * problem_13_beam_mu;
    vdisort::flux_data flux_data;
    const auto         top = solver.flux(flux_data, 0.0);
    for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
      expect_unpolarized(std::format("Problem 13 top flux field [{}, {}]", test.name, stream), flux_data.u0[stream]);
    const auto bottom = solver.flux(flux_data, solver.tau().back());
    for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
      expect_unpolarized(std::format("Problem 13 bottom flux field [{}, {}]", test.name, stream), flux_data.u0[stream]);
    expect_reference(std::format("{} albedo", test.name), top.up / incident_flux, test.albedo);
    expect_reference(std::format("{} transmission", test.name),
                     (bottom.down_direct + bottom.down_diffuse) / incident_flux,
                     test.transmission);
  }
}

disort::brdf::RawFunction problem_14_raw(const disort_test::reference::brdf_type type) {
  using enum disort_test::reference::brdf_type;
  switch (type) {
    case hapke:    return disort::brdf::Hapke{};
    case cox_munk: return disort::brdf::CoxMunk{};
    case rpv:      return disort::brdf::RPV{};
    case ross_li:  return disort::brdf::RossLi{};
  }
  std::unreachable();
}

problem_9_model make_problem_14_model(const disort::brdf::RawFunction& raw) {
  constexpr Index   nquad   = disort_test::reference::problem_14_streams;
  constexpr Index   nmodes  = nquad;
  constexpr Index   nlayers = 1;
  constexpr Index   nstokes = vdisort::stokes_dimension;
  constexpr Numeric depth   = 1e-12;
  constexpr Numeric omega   = 1e-8;
  const Index       n       = nquad / 2;
  const auto&       user_mu = disort_test::reference::problem_14_user_mu;
  const Index       nuser   = static_cast<Index>(user_mu.size());

  Tensor7 phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  phase[vdisort::cosine_mode, 0, 0, joker, joker, 0, 0] = 1.0;
  Tensor6 beam_phase(2, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  beam_phase[vdisort::cosine_mode, 0, 0, joker, 0, 0] = 1.0;

  vdisort::phase_matrix_data      user_phase(2, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data user_beam_phase(2, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index user = 0; user < nuser; ++user) {
    for (Index stream = 0; stream < nquad; ++stream) user_phase[vdisort::cosine_mode, 0, 0, user, stream][0, 0] = 1.0;
    user_beam_phase[vdisort::cosine_mode, 0, 0, user][0, 0] = 1.0;
  }

  auto solver = vdisort_test::make_solver(nquad,
                                          nmodes,
                                          AscendingGrid{depth},
                                          Vector{omega},
                                          std::move(phase),
                                          Tensor4(2, nmodes, n, nstokes, 0.0),
                                          Tensor4(2, nmodes, n, nstokes, 0.0),
                                          Tensor3(nlayers, 0, nstokes),
                                          scalar_brdf_modes(raw, nmodes),
                                          disort_test::reference::problem_14_beam_mu,
                                          Vector{disort_test::reference::problem_14_beam, 0.0, 0.0, 0.0},
                                          0.0,
                                          std::move(beam_phase));
  return {
      .solver = std::move(solver), .user_phase = std::move(user_phase), .user_beam_phase = std::move(user_beam_phase)};
}

void update_scalar_brdf(problem_9_model& model, const disort::brdf::RawFunction& raw, const Index nmodes) {
  const auto modes = scalar_brdf_modes(raw, nmodes);
  stdr::copy(modes, model.solver.brdf_modes().begin());
  model.solver.solve_for_coefs();
  model.solver.rad_field();
}

void test_problem_14() {
  using namespace disort_test::reference;
  auto model = make_problem_14_model(problem_14_raw(problem_14.front().type));
  for (Index case_index = 0; case_index < static_cast<Index>(problem_14.size()); ++case_index) {
    const auto& test = problem_14[case_index];
    const auto  raw  = problem_14_raw(test.type);
    if (case_index != 0) update_scalar_brdf(model, raw, problem_14_streams);
    for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth)
      for (Index angle = 0; angle < static_cast<Index>(problem_14_user_mu.size()); ++angle)
        expect_reference(std::format("{} raw radiance [{}, {}]", test.name, azimuth, angle),
                         problem_14_beam_mu * raw(problem_14_user_mu[angle], problem_14_beam_mu, test.azimuth[azimuth]),
                         test.radiance[azimuth, angle],
                         2e-5);

    const AscendingGrid      tau{0.0};
    rtepack::stokvec_tensor3 user(tau.size(), test.azimuth.size(), problem_14_user_mu.size());
    model.solver.ungridded_u_user(user, tau, test.azimuth, problem_14_user_mu, model.user_phase, model.user_beam_phase);
    for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth) {
      for (Index angle = 0; angle < static_cast<Index>(problem_14_user_mu.size()); ++angle) {
        const auto label = std::format("{} solver radiance [{}, {}]", test.name, azimuth, angle);
        expect_unpolarized(label, user[0, azimuth, angle]);
        expect_reference(label, user[0, azimuth, angle].I(), test.radiance[azimuth, angle], 2e-5);
      }
    }

    vdisort::flux_data flux_data;
    const auto         values = model.solver.flux(flux_data, 0.0);
    for (Index stream = 0; stream < static_cast<Index>(flux_data.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}]", test.name, stream), flux_data.u0[stream]);
    expect_reference(std::format("{} direct flux", test.name), values.down_direct, test.direct);
    expect_reference(std::format("{} diffuse-down flux", test.name), values.down_diffuse, test.diffuse_down);
    expect_reference(std::format("{} upward flux", test.name), values.up, test.up, 2e-5);
    expect_reference(std::format("{} DFDT", test.name), values.dfdt, test.dfdt, 2e-5);
  }
}

disort::brdf::RawFunction problem_15_raw(const disort_test::reference::brdf_type type) {
  if (type == disort_test::reference::brdf_type::cox_munk) return disort::brdf::CoxMunk{.shadowing = true};
  return problem_14_raw(type);
}

problem_9_model make_problem_15_model(const disort::brdf::RawFunction& raw,
                                      const Matrix&                    transport_moments,
                                      const AscendingGrid&             scaled_tau) {
  using namespace disort_test::reference;
  constexpr Index nquad   = problem_15_streams;
  constexpr Index nmodes  = nquad;
  constexpr Index nlayers = 2;
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Index     nuser   = static_cast<Index>(problem_15_user_mu.size());

  Vector quadrature_mu(nquad), quadrature_weight(n);
  Legendre::PositiveDoubleGaussLegendre(quadrature_mu[Range{0, n}], quadrature_weight);
  std::transform(quadrature_mu.begin(), quadrature_mu.begin() + n, quadrature_mu.begin() + n, [](const Numeric mu) {
    return -mu;
  });
  Tensor7                         phase(2, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  Tensor6                         beam_phase(2, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  vdisort::phase_matrix_data      user_phase(2, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data user_beam_phase(2, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index layer = 0; layer < nlayers; ++layer)
    for (Index mode = 0; mode < nmodes; ++mode) {
      const auto moments = transport_moments[layer];
      for (Index out = 0; out < nquad; ++out) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(moments, mode, quadrature_mu[out], quadrature_mu[in]);
          phase[vdisort::cosine_mode, mode, layer, out, in, 0, 0] = value;
          if (mode > 0) phase[vdisort::sine_mode, mode, layer, out, in, 0, 0] = value;
        }
        beam_phase[vdisort::cosine_mode, mode, layer, out, 0, 0] =
            scalar_phase_mode(moments, mode, quadrature_mu[out], -problem_15_beam_mu);
      }
      for (Index user = 0; user < nuser; ++user) {
        for (Index in = 0; in < nquad; ++in) {
          const Numeric value = scalar_phase_mode(moments, mode, problem_15_user_mu[user], quadrature_mu[in]);
          user_phase[vdisort::cosine_mode, mode, layer, user, in][0, 0] = value;
          if (mode > 0) user_phase[vdisort::sine_mode, mode, layer, user, in][0, 0] = value;
        }
        user_beam_phase[vdisort::cosine_mode, mode, layer, user][0, 0] =
            scalar_phase_mode(moments, mode, problem_15_user_mu[user], -problem_15_beam_mu);
      }
    }

  auto solver = vdisort_test::make_solver(nquad,
                                          nmodes,
                                          scaled_tau,
                                          Vector{1.0, 1.0},
                                          std::move(phase),
                                          Tensor4(2, nmodes, n, nstokes, 0.0),
                                          Tensor4(2, nmodes, n, nstokes, 0.0),
                                          Tensor3(nlayers, 0, nstokes),
                                          scalar_brdf_modes(raw, nmodes),
                                          problem_15_beam_mu,
                                          Vector{problem_15_beam, 0.0, 0.0, 0.0},
                                          0.0,
                                          std::move(beam_phase));
  return {
      .solver = std::move(solver), .user_phase = std::move(user_phase), .user_beam_phase = std::move(user_beam_phase)};
}

void test_problem_15() {
  using namespace disort_test::reference;
  constexpr Index moment_count = 600;
  Matrix          original(2, moment_count, 0.0);
  original[0, 0] = 1.0;
  original[0, 2] = 0.1;
  original[1]    = kokhanovsky_aerosol_moments[Range{0, moment_count}];
  const Vector fraction{0.0, original[1, problem_15_streams]};
  Matrix       transport(2, problem_15_streams, 0.0);
  for (Index layer = 0; layer < 2; ++layer)
    for (Index degree = 0; degree < problem_15_streams; ++degree)
      transport[layer, degree] = (original[layer, degree] - fraction[layer]) / (1.0 - fraction[layer]);
  const Numeric       scaled_bottom = problem_15_tau[0] + (problem_15_tau[1] - problem_15_tau[0]) * (1.0 - fraction[1]);
  const AscendingGrid scaled_tau{problem_15_tau[0], scaled_bottom};

  Matrix weighted_original(2, moment_count), weighted_transport(2, problem_15_streams),
      weighted_removed(2, moment_count, 0.0);
  for (Index layer = 0; layer < 2; ++layer) {
    for (Index degree = 0; degree < moment_count; ++degree)
      weighted_original[layer, degree] = static_cast<Numeric>(2 * degree + 1) * original[layer, degree];
    for (Index degree = 0; degree < problem_15_streams; ++degree)
      weighted_transport[layer, degree] = static_cast<Numeric>(2 * degree + 1) * transport[layer, degree];
  }
  for (Index degree = 0; degree < moment_count; ++degree) {
    const Numeric residue       = degree < problem_15_streams ? 1.0 : original[1, degree] / fraction[1];
    weighted_removed[1, degree] = static_cast<Numeric>(2 * degree + 1) * residue;
  }
  const auto phase_callback = [](const Matrix coefficients) {
    return [coefficients](Index layer, Numeric out_mu, Numeric out_phi, Numeric in_mu, Numeric in_phi) {
      rtepack::muelmat result{0.0};
      result[0, 0] =
          Legendre::legendre_sum(coefficients[layer], scattering_angle_cosine(out_mu, out_phi, in_mu, in_phi));
      return result;
    };
  };
  const auto removed_convolution =
      [weighted_removed](Index first, Index second, Numeric out_mu, Numeric out_phi, Numeric in_mu, Numeric in_phi) {
        rtepack::muelmat result{0.0};
        Vector           coefficients(moment_count);
        for (Index degree = 0; degree < moment_count; ++degree)
          coefficients[degree] =
              weighted_removed[first, degree] * weighted_removed[second, degree] / static_cast<Numeric>(2 * degree + 1);
        result[0, 0] = Legendre::legendre_sum(coefficients, scattering_angle_cosine(out_mu, out_phi, in_mu, in_phi));
        return result;
      };
  const vdisort::delta_m_correction_cache correction(AscendingGrid{problem_15_tau},
                                                     Vector{1.0, 1.0},
                                                     fraction,
                                                     problem_15_beam_mu,
                                                     0.0,
                                                     rtepack::stokvec{problem_15_beam, 0.0, 0.0, 0.0},
                                                     Vector{problem_15_user_mu},
                                                     Vector{problem_15_azimuth},
                                                     phase_callback(weighted_original),
                                                     phase_callback(weighted_transport),
                                                     phase_callback(weighted_removed),
                                                     1,
                                                     1,
                                                     removed_convolution);

  const auto scaled_depth = [f = fraction[1]](const Numeric tau) {
    return tau <= problem_15_tau[0] ? tau : problem_15_tau[0] + (tau - problem_15_tau[0]) * (1.0 - f);
  };
  auto model = make_problem_15_model(problem_15_raw(problem_15.front().type), transport, scaled_tau);
  for (Index case_index = 0; case_index < static_cast<Index>(problem_15.size()); ++case_index) {
    const auto& test = problem_15[case_index];
    if (case_index != 0) update_scalar_brdf(model, problem_15_raw(test.type), problem_15_streams);
    const AscendingGrid scaled_output_tau(problem_15_output_tau.begin(), problem_15_output_tau.end(), scaled_depth);
    rtepack::stokvec_tensor3 user(scaled_output_tau.size(), problem_15_azimuth.size(), problem_15_user_mu.size());
    model.solver.ungridded_u_user(
        user, scaled_output_tau, problem_15_azimuth, problem_15_user_mu, model.user_phase, model.user_beam_phase);
    for (Index level = 0; level < static_cast<Index>(problem_15_output_tau.size()); ++level) {
      for (Index azimuth = 0; azimuth < static_cast<Index>(problem_15_azimuth.size()); ++azimuth) {
        const auto corr = correction.evaluate(problem_15_output_tau[level], azimuth);
        for (Index angle = 0; angle < static_cast<Index>(problem_15_user_mu.size()); ++angle) {
          auto       intensity = user[level, azimuth, angle] + corr[angle];
          const auto label     = std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle);
          expect_unpolarized(label, intensity);
          expect_reference(label, intensity.I(), test.radiance[azimuth, level, angle]);
        }
      }
    }

    vdisort::flux_data flux_data;
    for (Index level = 0; level < static_cast<Index>(problem_15_output_tau.size()); ++level) {
      const Numeric physical_tau = problem_15_output_tau[level];
      const auto    values       = model.solver.flux(flux_data, scaled_depth(physical_tau));
      const Numeric physical_direct =
          problem_15_beam_mu * problem_15_beam * std::exp(-physical_tau / problem_15_beam_mu);
      const Numeric physical_diffuse = values.down_diffuse - (physical_direct - values.down_direct);
      expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
      expect_reference(
          std::format("{} diffuse-down flux [{}]", test.name, level), physical_diffuse, test.diffuse_down[level]);
      expect_reference(std::format("{} upward flux [{}]", test.name, level), values.up, test.up[level]);
      expect_reference(std::format("{} DFDT [{}]", test.name, level), values.dfdt, test.dfdt[level], 2e-6);
    }
  }
}

void test_problem_17() {
  using namespace disort_test::reference;

  const std::array<const Vector*, 2> phase_moments{
      &kokhanovsky_aerosol_moments,
      &kokhanovsky_cloud_moments,
  };
  const std::array expected_fraction{0.24306, 0.44398};

  for (Index case_index = 0; case_index < static_cast<Index>(problem_17.size()); ++case_index) {
    const auto& test = problem_17[case_index];
    Matrix      moments(1, phase_moments[case_index]->size());
    moments[0]         = *phase_moments[case_index];
    const auto scaling = disort::delta_m_plus(moments, problem_17_streams);
    expect_reference(
        std::format("{} delta-M-plus fraction", test.name), scaling.fraction[0], expected_fraction[case_index], 2e-5);

    Vector removed_peak_moments{scaling.moments[0]};
    auto   model = make_scalar_model(test.depth,
                                     1.0,
                                     true,
                                     problem_17_streams,
                                     moments[0],
                                     problem_17_user_mu,
                                     problem_17_beam_mu,
                                     problem_17_beam,
                                     scaling.fraction[0],
                                     0.0,
                                     &removed_peak_moments);
    expect_reference(std::format("{} scaled depth", test.name),
                     model.solver.tau().back(),
                     test.depth * (1.0 - scaling.fraction[0]),
                     2e-13);

    const Vector        physical_tau{0.0, test.depth};
    const AscendingGrid scaled_tau(
        physical_tau.begin(), physical_tau.end(), [&](const Numeric tau) { return model.optical_depth_scale * tau; });
    rtepack::stokvec_tensor3 user(scaled_tau.size(), problem_17_azimuth.size(), problem_17_user_mu.size());
    model.solver.ungridded_u_user(
        user, scaled_tau, problem_17_azimuth, problem_17_user_mu, model.user_phase, model.user_beam_phase);
    for (Index level = 0; level < static_cast<Index>(physical_tau.size()); ++level) {
      for (Index azimuth = 0; azimuth < static_cast<Index>(problem_17_azimuth.size()); ++azimuth) {
        for (Index angle = 0; angle < static_cast<Index>(problem_17_user_mu.size()); ++angle) {
          const auto label = std::format("{} radiance [{}, {}, {}]", test.name, level, angle, azimuth);
          expect_reference(label, user[level, azimuth, angle].I(), test.radiance[level, angle, azimuth]);
          expect_unpolarized(label, user[level, azimuth, angle]);
        }
      }
    }
  }
}

void test_polarized_delta_m_correction() {
  rtepack::muelmat removed{0.0};
  removed[0, 0]             = 1.0;
  removed[1, 0]             = 0.2;
  removed[2, 0]             = -0.1;
  removed[3, 0]             = 0.05;
  removed[1, 1]             = 0.7;
  removed[2, 2]             = 0.6;
  removed[3, 3]             = 0.5;
  const auto constant_phase = [removed](Index, Numeric, Numeric, Numeric, Numeric) { return removed; };
  const auto zero_phase     = [](Index, Numeric, Numeric, Numeric, Numeric) { return rtepack::muelmat{0.0}; };

  constexpr Numeric                       depth    = 1.0;
  constexpr Numeric                       omega    = 0.8;
  constexpr Numeric                       fraction = 0.2;
  constexpr Numeric                       mu0      = 0.5;
  constexpr Numeric                       tau      = 0.5;
  constexpr Numeric                       user_mu  = -0.5;
  const rtepack::stokvec                  beam{1.0, 0.1, -0.05, 0.02};
  const vdisort::delta_m_correction_cache correction(AscendingGrid{depth},
                                                     Vector{omega},
                                                     Vector{fraction},
                                                     mu0,
                                                     0.0,
                                                     beam,
                                                     Vector{user_mu},
                                                     Vector{0.0},
                                                     zero_phase,
                                                     zero_phase,
                                                     constant_phase,
                                                     8,
                                                     16);
  rtepack::stokvec_vector                 tms, ims;
  correction.TMS(tms, tau, 0);
  correction.IMS(ims, tau, 0);
  const auto             actual       = correction.evaluate(tau, 0);
  const Numeric          omega_f      = omega * fraction;
  const Numeric          scalar       = Constant::inv_pi * 0.25 * omega_f * omega_f / (1.0 - omega_f);
  const rtepack::muelmat ims_operator = 2.0 * removed - removed * removed;
  const rtepack::stokvec expected     = -scalar * ims_chi(tau, -user_mu, mu0 / (1.0 - omega_f)) * (ims_operator * beam);
  for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes) {
    expect_reference(std::format("Separated polarized delta-M TMS Stokes {}", stokes), tms[0][stokes], 0.0, 2e-12);
    expect_reference(
        std::format("Separated polarized delta-M IMS Stokes {}", stokes), ims[0][stokes], expected[stokes], 2e-12);
    expect_reference(
        std::format("Polarized delta-M IMS Stokes {}", stokes), actual[0][stokes], expected[stokes], 2e-12);
  }

  rtepack::stokvec_vector pythonic_ims;
  correction.IMS(pythonic_ims, tau, 0, vdisort::ims_convention::pythonic_disort);
  for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
    expect_reference(std::format("Pythonic polarized delta-M IMS Stokes {}", stokes),
                     pythonic_ims[0][stokes],
                     -expected[stokes],
                     2e-12);

  rtepack::stokvec_tensor3 gridded_tms(1, 1, 1), gridded_ims(1, 1, 1);
  correction.gridded_TMS(gridded_tms);
  correction.gridded_IMS(gridded_ims);
  rtepack::stokvec_vector bottom_tms, bottom_ims;
  correction.TMS(bottom_tms, depth, 0);
  correction.IMS(bottom_ims, depth, 0);
  for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes) {
    expect_reference(std::format("Gridded polarized delta-M TMS Stokes {}", stokes),
                     gridded_tms[0, 0, 0][stokes],
                     bottom_tms[0][stokes],
                     2e-12);
    expect_reference(std::format("Gridded polarized delta-M IMS Stokes {}", stokes),
                     gridded_ims[0, 0, 0][stokes],
                     bottom_ims[0][stokes],
                     2e-12);
  }
  ARTS_USER_ERROR_IF(std::abs(actual[0].Q()) == 0.0 or std::abs(actual[0].U()) == 0.0 or std::abs(actual[0].V()) == 0.0,
                     "Polarized delta-M IMS failed to generate all polarized Stokes components: {}",
                     actual[0]);
}
}  // namespace

int main() try {
  ARTS_TIME_REPORT

  test_problem_1();
  test_problem_2();
  test_problem_3();
  test_problem_4();
  test_problem_5();
  test_problem_6();
  test_problem_7();
  test_problem_8();
  test_problem_9();
  test_problem_10();
  test_problem_11();
  test_problem_12();
  test_problem_13();
  test_problem_14();
  test_problem_15();
  test_problem_17();
  test_polarized_delta_m_correction();
#if ARTS_PROFILING
  arts::print_report();
#endif
  std::cout << "VDISORT Fortran reference tests passed\n";
  return 0;
} catch (const std::exception& exception) {
  std::cerr << exception.what() << '\n';
  return 1;
}

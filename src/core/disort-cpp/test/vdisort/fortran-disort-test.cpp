#include <arts_constants.h>
#include <disort-brdf.h>
#include <legendre.h>
#include <vdisort.h>

#include <cmath>
#include <iostream>

#include "../reference-data.h"

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
                                       const Numeric          delta_m_fraction = 0.0,
                                       const Numeric          phi0             = 0.0) {
  const Numeric optical_depth_scale = 1.0 - physical_omega * delta_m_fraction;
  const Numeric transport_omega     = physical_omega * (1.0 - delta_m_fraction) / optical_depth_scale;
  Vector        transport_moments(std::min(nquad, static_cast<Index>(moments.size())));
  for (Index degree = 0; degree < static_cast<Index>(transport_moments.size()); ++degree)
    transport_moments[degree] = (moments[degree] - delta_m_fraction) / (1.0 - delta_m_fraction);

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

  vdisort::main_data solver(nquad,
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

    const auto [diffuse_down, direct] = model.solver.flux_down(flux, output_tau[level]);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);

    const Numeric up = model.solver.flux_up(flux, output_tau[level]);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), up, test.up[level]);
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

    const auto [diffuse_down, direct] = model.solver.flux_down(flux, output_tau[level]);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);

    const Numeric up = model.solver.flux_up(flux, output_tau[level]);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), up, test.up[level]);
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

  vdisort::user_u_data user;
  vdisort::flux_data   flux;
  for (Index level = 0; level < 2; ++level) {
    const Numeric physical_tau = output_tau[level];
    const Numeric scaled_tau   = model.optical_depth_scale * physical_tau;
    model.solver.u_user(user, scaled_tau, 0.0, user_mu, model.user_phase, model.user_beam_phase);
    const Vector correction = scalar_delta_m_correction(model, physical_tau, 0.0, user_mu);
    for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
      user.intensities[angle].I() += correction[angle];
      const auto label             = std::format("{} radiance [{}, {}]", test.name, level, angle);
      expect_reference(label, user.intensities[angle].I(), test.radiance[level, angle]);
      expect_unpolarized(label, user.intensities[angle]);
    }

    const auto [scaled_diffuse_down, scaled_direct] = model.solver.flux_down(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    const Numeric physical_direct = model.mu0 * model.beam_intensity * std::exp(-physical_tau / model.mu0);
    const Numeric diffuse_down    = scaled_diffuse_down - (physical_direct - scaled_direct);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    const Numeric up = model.solver.flux_up(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), up, test.up[level]);
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

  vdisort::user_u_data user;
  for (Index azimuth = 0; azimuth < static_cast<Index>(test.azimuth.size()); ++azimuth) {
    for (Index level = 0; level < static_cast<Index>(disort_test::reference::problem_4_output_tau.size()); ++level) {
      const Numeric physical_tau = disort_test::reference::problem_4_output_tau[level];
      const Numeric scaled_tau   = model.optical_depth_scale * physical_tau;
      model.solver.u_user(user, scaled_tau, test.azimuth[azimuth], user_mu, model.user_phase, model.user_beam_phase);
      const Vector correction = scalar_delta_m_correction(model, physical_tau, test.azimuth[azimuth], user_mu);
      for (Index angle = 0; angle < static_cast<Index>(user_mu.size()); ++angle) {
        user.intensities[angle].I() += correction[angle];
        const auto label             = std::format("{} radiance [{}, {}, {}]", test.name, azimuth, level, angle);
        expect_reference(label, user.intensities[angle].I(), test.radiance[azimuth, level, angle]);
        expect_unpolarized(label, user.intensities[angle]);
      }
    }
  }

  vdisort::flux_data flux;
  for (Index level = 0; level < static_cast<Index>(disort_test::reference::problem_4_output_tau.size()); ++level) {
    const Numeric physical_tau                      = disort_test::reference::problem_4_output_tau[level];
    const Numeric scaled_tau                        = model.optical_depth_scale * physical_tau;
    const auto [scaled_diffuse_down, scaled_direct] = model.solver.flux_down(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    const Numeric physical_direct = model.mu0 * model.beam_intensity * std::exp(-physical_tau / model.mu0);
    const Numeric diffuse_down    = scaled_diffuse_down - (physical_direct - scaled_direct);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    const Numeric up = model.solver.flux_up(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), up, test.up[level]);
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

    const auto [scaled_diffuse_down, scaled_direct] = model.solver.flux_down(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    const Numeric physical_direct = model.mu0 * model.beam_intensity * std::exp(-physical_tau / model.mu0);
    const Numeric diffuse_down    = scaled_diffuse_down - (physical_direct - scaled_direct);
    expect_reference(std::format("{} direct flux [{}]", test.name, level), physical_direct, test.direct[level]);
    expect_reference(
        std::format("{} diffuse-down flux [{}]", test.name, level), diffuse_down, test.diffuse_down[level]);
    const Numeric up = model.solver.flux_up(flux, scaled_tau);
    for (Index stream = 0; stream < static_cast<Index>(flux.u0.size()); ++stream)
      expect_unpolarized(std::format("{} flux field [{}, {}]", test.name, level, stream), flux.u0[stream]);
    expect_reference(std::format("{} up flux [{}]", test.name, level), up, test.up[level]);
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
  constexpr Index nazimuth = 100;
  Vector          azimuth_node(nazimuth), azimuth_weight(nazimuth);
  Legendre::PositiveDoubleGaussLegendre(azimuth_node, azimuth_weight);
  std::vector<vdisort::BDRF> result;
  result.reserve(static_cast<std::size_t>(number_of_modes));
  for (Index m = 0; m < number_of_modes; ++m) {
    const auto coefficient = [raw, m, azimuth_node, azimuth_weight](rtepack::muelmat_matrix_view out,
                                                                    const ConstVectorView&       outgoing,
                                                                    const ConstVectorView&       incoming) {
      out = rtepack::muelmat{0.0};
      for (Index i = 0; i < static_cast<Index>(outgoing.size()); ++i)
        for (Index j = 0; j < static_cast<Index>(incoming.size()); ++j) {
          Numeric coefficient = 0.0;
          for (Index k = 0; k < nazimuth; ++k) {
            const Numeric azimuth  = Constant::pi * azimuth_node[k];
            coefficient           += azimuth_weight[k] * raw(outgoing[i], std::abs(incoming[j]), azimuth) *
                                     std::cos(static_cast<Numeric>(m) * azimuth);
          }
          out[i, j][0, 0] = 2.0 * (m == 0 ? 1.0 : 2.0) * coefficient;
        }
    };
    const vdisort::BDRF::func_t sine{[](rtepack::muelmat_matrix_view out,
                                        const ConstVectorView&,
                                        const ConstVectorView&) { out = rtepack::muelmat{0.0}; }};
    result.push_back(vdisort::BDRF{.cosine = {coefficient}, .sine = sine});
  }
  return result;
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

  return vdisort::main_data(nquad,
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
  const Numeric absorption = 1.0 - test.single_scattering_albedo;
  Tensor3       source(1, 2, nstokes, 0.0);
  source[0, 0, 0] = absorption * atmosphere_top_planck;
  source[0, 1, 0] = absorption * (atmosphere_bottom_planck - atmosphere_top_planck) / test.optical_depth;

  return vdisort::main_data(nquad,
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

  vdisort::main_data solver(nquad,
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
}  // namespace

int main() try {
  test_problem_1();
  test_problem_2();
  test_problem_3();
  test_problem_4();
  test_problem_5();
  test_problem_6();
  test_problem_7();
  test_problem_8();
  std::cout << "VDISORT Fortran reference tests passed\n";
  return 0;
} catch (const std::exception& exception) {
  std::cerr << exception.what() << '\n';
  return 1;
}

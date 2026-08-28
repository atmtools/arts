#include <arts_constants.h>
#include <legendre.h>
#include <vdisort.h>

#include <cmath>
#include <iostream>

#include "../reference-data.h"

namespace {
constexpr Numeric reference_tolerance    = 7e-5;
constexpr Numeric polarization_tolerance = 2e-9;

void expect_reference(const std::string_view name, const Numeric actual, const Numeric expected) {
  ARTS_USER_ERROR_IF(std::abs(actual - expected) > reference_tolerance * std::max(1.0, std::abs(expected)),
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

scalar_vdisort_model make_scalar_model(const disort_test::reference::single_layer_case& test,
                                       const Index                                      nquad,
                                       const ConstVectorView&                           moments,
                                       const ConstVectorView&                           user_mu,
                                       const Numeric                                    mu0,
                                       const Numeric                                    beam_intensity) {
  Index highest_nonzero_degree = 0;
  for (Index degree = 0; degree < static_cast<Index>(moments.size()); ++degree)
    if (moments[degree] != 0.0) highest_nonzero_degree = degree;
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
        const Numeric value = scalar_phase_mode(moments, m, quadrature_mu[i], quadrature_mu[j]);
        phase[vdisort::cosine_mode, m, 0, i, j, 0, 0] = value;
        if (m > 0) phase[vdisort::sine_mode, m, 0, i, j, 0, 0] = value;
      }

  Tensor4 up(nalpha, nmodes, n, nstokes, 0.0);
  Tensor4 down(nalpha, nmodes, n, nstokes, 0.0);
  if (not test.beam) down[vdisort::cosine_mode, 0, joker, 0] = 1.0;

  Tensor6 beam_phase(nalpha, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nquad; ++i)
      beam_phase[vdisort::cosine_mode, m, 0, i, 0, 0] = scalar_phase_mode(moments, m, quadrature_mu[i], -mu0);

  vdisort::main_data solver(nquad,
                            nmodes,
                            AscendingGrid{test.depth},
                            Vector{test.omega},
                            std::move(phase),
                            std::move(up),
                            std::move(down),
                            Tensor3(nlayers, 0, nstokes),
                            {},
                            mu0,
                            test.beam ? Vector{beam_intensity, 0.0, 0.0, 0.0} : Vector(4, 0.0),
                            0.0,
                            std::move(beam_phase));

  vdisort::phase_matrix_data user_phase(nalpha, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nuser; ++i)
      for (Index j = 0; j < nquad; ++j) {
        const Numeric value = scalar_phase_mode(moments, m, user_mu[i], quadrature_mu[j]);
        user_phase[vdisort::cosine_mode, m, 0, i, j][0, 0] = value;
        if (m > 0) user_phase[vdisort::sine_mode, m, 0, i, j][0, 0] = value;
      }
  vdisort::beam_phase_matrix_data user_beam_phase(nalpha, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index m = 0; m < nmodes; ++m)
    for (Index i = 0; i < nuser; ++i)
      user_beam_phase[vdisort::cosine_mode, m, 0, i][0, 0] = scalar_phase_mode(moments, m, user_mu[i], -mu0);

  return {std::move(solver), std::move(user_phase), std::move(user_beam_phase)};
}

void run_problem_1_case(const disort_test::reference::single_layer_case& test) {
  const auto& user_mu = disort_test::reference::problem_1_user_mu;
  Vector      moments(17, 0.0);
  moments[0]         = 1.0;
  auto         model = make_scalar_model(test,
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
  auto         model = make_scalar_model(test,
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
}  // namespace

int main() try {
  test_problem_1();
  test_problem_2();
  std::cout << "VDISORT Fortran reference tests passed\n";
  return 0;
} catch (const std::exception& exception) {
  std::cerr << exception.what() << '\n';
  return 1;
}

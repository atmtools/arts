#include <arts_constants.h>
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

struct scalar_isotropic_model {
  vdisort::main_data              solver;
  vdisort::phase_matrix_data      user_phase;
  vdisort::beam_phase_matrix_data user_beam_phase;
};

scalar_isotropic_model make_scalar_isotropic_model(const disort_test::reference::single_layer_case& test,
                                                   const ConstVectorView&                           user_mu) {
  constexpr Index nquad   = disort_test::reference::problem_1_streams;
  constexpr Index nmodes  = 1;
  constexpr Index nlayers = 1;
  constexpr Index nalpha  = 2;
  constexpr Index nstokes = vdisort::stokes_dimension;
  const Index     n       = nquad / 2;
  const Index     nuser   = static_cast<Index>(user_mu.size());

  Tensor7 phase(nalpha, nmodes, nlayers, nquad, nquad, nstokes, nstokes, 0.0);
  phase[vdisort::cosine_mode, 0, 0, joker, joker, 0, 0] = 1.0;

  Tensor4 up(nalpha, nmodes, n, nstokes, 0.0);
  Tensor4 down(nalpha, nmodes, n, nstokes, 0.0);
  if (not test.beam) down[vdisort::cosine_mode, 0, joker, 0] = 1.0;

  Tensor6 beam_phase(nalpha, nmodes, nlayers, nquad, nstokes, nstokes, 0.0);
  beam_phase[vdisort::cosine_mode, 0, 0, joker, 0, 0] = 1.0;

  const Numeric      mu0 = disort_test::reference::problem_1_beam_mu;
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
                            test.beam ? Vector{Constant::pi / mu0, 0.0, 0.0, 0.0} : Vector(4, 0.0),
                            0.0,
                            std::move(beam_phase));

  vdisort::phase_matrix_data user_phase(nalpha, nmodes, nlayers, nuser, nquad, rtepack::muelmat{0.0});
  for (Index i = 0; i < nuser; ++i)
    for (Index j = 0; j < nquad; ++j) user_phase[vdisort::cosine_mode, 0, 0, i, j][0, 0] = 1.0;
  vdisort::beam_phase_matrix_data user_beam_phase(nalpha, nmodes, nlayers, nuser, rtepack::muelmat{0.0});
  for (Index i = 0; i < nuser; ++i) user_beam_phase[vdisort::cosine_mode, 0, 0, i][0, 0] = 1.0;

  return {std::move(solver), std::move(user_phase), std::move(user_beam_phase)};
}

void run_problem_1_case(const disort_test::reference::single_layer_case& test) {
  const auto&  user_mu = disort_test::reference::problem_1_user_mu;
  auto         model   = make_scalar_isotropic_model(test, user_mu);
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
}  // namespace

int main() try {
  test_problem_1();
  std::cout << "VDISORT Fortran reference tests passed\n";
  return 0;
} catch (const std::exception& exception) {
  std::cerr << exception.what() << '\n';
  return 1;
}

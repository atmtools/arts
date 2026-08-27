#include "vdisort.h"

#include <arts_constants.h>

#include <cmath>
#include <iostream>

#include "disort.h"
#include "vdisort-scalar-test-adapter.h"

namespace {
constexpr Numeric tolerance = 2e-9;

void expect_close(const Numeric actual, const Numeric expected, const std::string_view what) {
  const Numeric scale = 1.0 + std::abs(expected);
  ARTS_USER_ERROR_IF(std::abs(actual - expected) > tolerance * scale,
                     "{}: expected {}, got {} (difference {})",
                     what,
                     expected,
                     actual,
                     actual - expected);
}

vdisort::main_data make_vdisort(const Index                nquad,
                                AscendingGrid              tau,
                                Vector                     omega,
                                Tensor7                    phase,
                                Tensor4                    up,
                                Tensor4                    down,
                                Tensor3                    source     = {},
                                Vector                     beam       = Vector(4, 0.0),
                                Tensor6                    beam_phase = {},
                                std::vector<vdisort::BDRF> brdf       = {}) {
  if (source.size() == 0) source.resize(tau.size(), 0, 4);
  const Index nfourier = phase.shape()[1];
  return vdisort::main_data(nquad,
                            nfourier,
                            std::move(tau),
                            std::move(omega),
                            std::move(phase),
                            std::move(up),
                            std::move(down),
                            std::move(source),
                            std::move(brdf),
                            0.5,
                            std::move(beam),
                            0.0,
                            std::move(beam_phase));
}

void test_polarized_absorption() {
  constexpr Index nquad  = 2;
  constexpr Index n      = nquad / 2;
  constexpr Index modes  = 1;
  constexpr Index layers = 2;
  Tensor7         phase(2, modes, layers, nquad, nquad, 4, 4, 0.0);
  Tensor4         up(2, modes, n, 4, 0.0), down(2, modes, n, 4, 0.0);

  const Vector top{2.0, -0.3, 0.2, -0.1};
  const Vector bottom{1.0, 0.25, -0.15, 0.05};
  for (Index s = 0; s < 2; ++s) {
    down[vdisort::cosine_mode, 0, 0, s] = top[s];
    up[vdisort::cosine_mode, 0, 0, s]   = bottom[s];
  }
  for (Index s = 2; s < 4; ++s) {
    down[vdisort::sine_mode, 0, 0, s] = top[s];
    up[vdisort::sine_mode, 0, 0, s]   = bottom[s];
  }

  const Numeric depth = 0.7;
  auto          model = make_vdisort(
      nquad, AscendingGrid{0.25, depth}, Vector{0.0, 0.0}, std::move(phase), std::move(up), std::move(down));
  vdisort::u_data field;
  const Numeric   tau = 0.3;
  model.u(field, tau, 1.3);
  const Numeric mu = model.mu()[0];
  for (Index s = 0; s < 4; ++s) {
    expect_close(field.intensities[0, s], bottom[s] * std::exp((tau - depth) / mu), "upward absorbing Stokes");
    expect_close(field.intensities[1, s], top[s] * std::exp(-tau / mu), "downward absorbing Stokes");
  }
}

void test_scalar_limit() {
  constexpr Index     nquad = 4;
  constexpr Index     n     = nquad / 2;
  constexpr Index     modes = 1;
  const AscendingGrid tau{0.8};
  const Vector        omega{0.4};

  Tensor7 phase(2, modes, 1, nquad, nquad, 4, 4, 0.0);
  for (Index i = 0; i < nquad; ++i)
    for (Index j = 0; j < nquad; ++j) phase[vdisort::cosine_mode, 0, 0, i, j, 0, 0] = 1.0;

  Tensor4 vup(2, modes, n, 4, 0.0), vdown(2, modes, n, 4, 0.0);
  vup[vdisort::cosine_mode, 0, joker, 0]   = 0.5;
  vdown[vdisort::cosine_mode, 0, joker, 0] = 2.0;
  Tensor6 beam_phase(2, modes, 1, nquad, 4, 4, 0.0);
  beam_phase[vdisort::cosine_mode, 0, 0, joker, 0, 0] = 1.0;
  auto vector_model                                   = make_vdisort(nquad,
                                                                     tau,
                                                                     omega,
                                                                     std::move(phase),
                                                                     std::move(vup),
                                                                     std::move(vdown),
                                                                     {},
                                                                     Vector{1.0, 0.0, 0.0, 0.0},
                                                                     std::move(beam_phase));

  Matrix legendre(1, nquad, 0.0);
  legendre[0, 0] = 1.0;
  Matrix            sup(modes, n, 0.5), sdown(modes, n, 2.0), source(1, 0);
  disort::main_data scalar_model(nquad,
                                 nquad,
                                 modes,
                                 tau,
                                 omega,
                                 std::move(legendre),
                                 std::move(sup),
                                 std::move(sdown),
                                 Vector(1, 0.0),
                                 std::move(source),
                                 {},
                                 0.5,
                                 1.0,
                                 0.0);

  for (const Numeric optical_depth : {0.0, 0.2, 0.8}) {
    vdisort::u0_data vu;
    disort::u0_data  su;
    vector_model.u0(vu, optical_depth);
    scalar_model.u0(su, optical_depth);
    for (Index i = 0; i < nquad; ++i) {
      expect_close(vu.u0[i, 0], su.u0[i], "scalar-limit I");
      expect_close(vu.u0[i, 1], 0.0, "scalar-limit Q");
      expect_close(vu.u0[i, 2], 0.0, "scalar-limit U");
      expect_close(vu.u0[i, 3], 0.0, "scalar-limit V");
    }
  }
}

void test_scalar_linear_source_limit() {
  constexpr Index     nquad = 4;
  constexpr Index     modes = 1;
  const AscendingGrid tau{0.8};
  const Vector        omega{0.4};
  Matrix              legendre(1, nquad, 0.0);
  legendre[0, 0] = 1.0;
  Matrix up(modes, nquad / 2, 0.0), down(modes, nquad / 2, 0.0);
  Matrix source{Vector{2.0, -0.7}.reshape(1, 2)};

  disort::main_data scalar(
      nquad, nquad, modes, tau, omega, legendre, up, down, Vector(1, 0.0), source, {}, 0.5, 0.0, 0.0);
  vdisort_scalar_test::main_data vector(nquad,
                                        nquad,
                                        modes,
                                        tau,
                                        omega,
                                        std::move(legendre),
                                        std::move(up),
                                        std::move(down),
                                        Vector(1, 0.0),
                                        std::move(source),
                                        {},
                                        0.5,
                                        0.0,
                                        0.0);
  for (const Numeric optical_depth : {0.0, 0.2, 0.8}) {
    disort::u_data              scalar_field;
    vdisort_scalar_test::u_data vector_field;
    scalar.u(scalar_field, optical_depth, 0.0);
    vector.u(vector_field, optical_depth, 0.0);
    for (Index stream = 0; stream < nquad; ++stream)
      expect_close(vector_field.intensities[stream], scalar_field.intensities[stream], "linear-source scalar limit");
  }
}

void test_conservative_reflecting_source_limit() {
  constexpr Index     nquad = 4;
  constexpr Index     modes = 1;
  const AscendingGrid tau{8.0};
  const Vector        omega{1.0 - 1e-6};
  Matrix              legendre(1, nquad + 1, 0.0);
  legendre[0, 0] = 1.0;
  Matrix                    up(modes, nquad / 2, 0.0), down(modes, nquad / 2, 0.0);
  Matrix                    source{Vector{1.0 / (1.0 - omega[0]), 0.0}.reshape(1, 2)};
  std::vector<disort::BDRF> brdf{disort::BDRF{[](auto value, auto&, auto&) { value = 1.0; }}};

  disort::main_data scalar(
      nquad, nquad, modes, tau, omega, legendre, up, down, Vector(1, 0.0), source, brdf, 0.5, 0.0, 0.0);
  vdisort_scalar_test::main_data vector(nquad,
                                        nquad,
                                        modes,
                                        tau,
                                        omega,
                                        std::move(legendre),
                                        std::move(up),
                                        std::move(down),
                                        Vector(1, 0.0),
                                        std::move(source),
                                        std::move(brdf),
                                        0.5,
                                        0.0,
                                        0.0);
  for (const Numeric optical_depth : {0.2, 4.0, 8.0}) {
    disort::u_data              scalar_field;
    vdisort_scalar_test::u_data vector_field;
    scalar.u(scalar_field, optical_depth, 0.0);
    vector.u(vector_field, optical_depth, 0.0);
    for (Index stream = 0; stream < nquad; ++stream)
      expect_close(vector_field.intensities[stream],
                   scalar_field.intensities[stream],
                   "conservative reflecting source scalar limit");
  }
}

void test_vector_source() {
  constexpr Index nquad = 2;
  constexpr Index n     = nquad / 2;
  Tensor7         phase(2, 1, 1, nquad, nquad, 4, 4, 0.0);
  Tensor4         up(2, 1, n, 4, 0.0), down(2, 1, n, 4, 0.0);
  Tensor3         source(1, 1, 4, 0.0);
  const Vector    q{1.0, -0.2, 0.1, 0.05};
  source[0, 0] = q;

  const Numeric depth = 0.6;
  auto          model = make_vdisort(
      nquad, AscendingGrid{depth}, Vector{0.0}, std::move(phase), std::move(up), std::move(down), std::move(source));
  const Numeric   tau = 0.25;
  const Numeric   mu  = model.mu()[0];
  vdisort::u_data field;
  model.u(field, tau, 0.0);
  for (Index s = 0; s < 4; ++s) {
    expect_close(field.intensities[0, s], q[s] * (1.0 - std::exp(-(depth - tau) / mu)), "upward vector source");
    expect_close(field.intensities[1, s], q[s] * (1.0 - std::exp(-tau / mu)), "downward vector source");
  }
}

void test_polarized_brdf() {
  constexpr Index nquad = 2;
  constexpr Index n     = nquad / 2;
  Tensor7         phase(2, 1, 1, nquad, nquad, 4, 4, 0.0);
  Tensor4         up(2, 1, n, 4, 0.0), down(2, 1, n, 4, 0.0);
  down[vdisort::cosine_mode, 0, 0, 0] = 1.0;
  down[vdisort::cosine_mode, 0, 0, 1] = 0.2;

  const Numeric reflectance = 0.3;
  const auto    callback    = [reflectance](rtepack::muelmat_matrix_view value, const auto&, const auto&) {
    value              = rtepack::muelmat{0.0};
    const Index blocks = std::min(value.nrows(), value.ncols());
    for (Index i = 0; i < blocks; ++i)
      for (Index s = 0; s < 4; ++s) value[i, i][s, s] = reflectance;
  };
  std::vector<vdisort::BDRF> brdf{{.cosine = {callback}, .sine = {callback}}};

  const Numeric    depth = 0.2;
  auto             model = make_vdisort(nquad,
                                        AscendingGrid{depth},
                                        Vector{0.0},
                                        std::move(phase),
                                        std::move(up),
                                        std::move(down),
                                        {},
                                        Vector(4, 0.0),
                                        {},
                                        std::move(brdf));
  vdisort::u0_data field;
  model.u0(field, depth);
  const Numeric mu               = model.mu()[0];
  const Numeric reflected_factor = Constant::pi * mu * reflectance;
  expect_close(field.u0[0, 0], reflected_factor * std::exp(-depth / mu), "polarized BRDF I");
  expect_close(field.u0[0, 1], 0.2 * reflected_factor * std::exp(-depth / mu), "polarized BRDF Q");
}

void test_complex_uv_eigenmodes() {
  constexpr Index nquad = 2;
  constexpr Index n     = nquad / 2;
  constexpr Index modes = 1;
  Tensor7         phase(2, modes, 1, nquad, nquad, 4, 4, 0.0);
  for (Index alpha = 0; alpha < 2; ++alpha) {
    for (Index i = 0; i < nquad; ++i) {
      phase[alpha, 0, 0, i, i, 2, 3] = 4.0;
      phase[alpha, 0, 0, i, i, 3, 2] = -4.0;
    }
  }

  Tensor4 up(2, modes, n, 4, 0.0), down(2, modes, n, 4, 0.0);
  down[vdisort::sine_mode, 0, 0, 2] = 0.2;
  up[vdisort::sine_mode, 0, 0, 3]   = 0.1;
  auto model = make_vdisort(nquad, AscendingGrid{0.4}, Vector{0.5}, std::move(phase), std::move(up), std::move(down));
  ARTS_USER_ERROR_IF(not model.has_complex_eigensolutions(),
                     "The U-V coupling test did not produce the expected complex eigenpairs");

  vdisort::u_data field;
  model.u(field, 0.2, 0.7);
  for (Index i = 0; i < nquad; ++i)
    for (Index s = 0; s < 4; ++s) {
      const Numeric value = field.intensities[i, s];
      ARTS_USER_ERROR_IF(not std::isfinite(value), "Complex-eigenmode test produced {}", value);
    }
}

void test_combined_matrix_transform() {
  Tensor6 cosine(2, 1, 1, 1, 4, 4, 0.0), sine(2, 1, 1, 1, 4, 4, 0.0);
  cosine[1, 0, 0, 0, 0, 2] = 3.0;
  sine[1, 0, 0, 0, 0, 2]   = 5.0;
  const Tensor7 combined   = vdisort::combine_phase_matrices(cosine, sine);
  expect_close(combined[vdisort::cosine_mode, 1, 0, 0, 0, 0, 2], -5.0, "Eq. 81 cosine sign");
  expect_close(combined[vdisort::sine_mode, 1, 0, 0, 0, 0, 2], 5.0, "Eq. 81 sine sign");

  rtepack::muelmat_tensor4 native_cosine(2, 1, 1, 1, rtepack::muelmat{0.0});
  rtepack::muelmat_tensor4 native_sine(2, 1, 1, 1, rtepack::muelmat{0.0});
  native_cosine[1, 0, 0, 0][0, 2] = 3.0;
  native_sine[1, 0, 0, 0][0, 2]   = 5.0;
  const auto native_combined = vdisort::combine_phase_matrices(native_cosine, native_sine);
  expect_close(native_combined[vdisort::cosine_mode, 1, 0, 0, 0][0, 2], -5.0, "native Eq. 81 cosine sign");
  expect_close(native_combined[vdisort::sine_mode, 1, 0, 0, 0][0, 2], 5.0, "native Eq. 81 sine sign");
}
}  // namespace

// VDISORT SCALAR TEST PORT BEGIN: compile the complete, unchanged scalar
// driver against vdisort_scalar_test::main_data.  Renaming its main function
// makes it callable from the polarized VDISORT test executable.
#define disort vdisort_scalar_test
#define main   run_all_scalar_disort_tests_through_vdisort
#include "disotest.cpp"
#undef main
#undef disort
// VDISORT SCALAR TEST PORT END

int main() try {
  test_polarized_absorption();
  test_scalar_limit();
  test_scalar_linear_source_limit();
  test_conservative_reflecting_source_limit();
  test_vector_source();
  test_polarized_brdf();
  test_complex_uv_eigenmodes();
  test_combined_matrix_transform();
  ARTS_USER_ERROR_IF(run_all_scalar_disort_tests_through_vdisort() != EXIT_SUCCESS,
                     "The scalar DISORT suite failed through the VDISORT port");
  std::cout << "vdisort tests passed\n";
  return 0;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}

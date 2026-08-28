#include <arts_constants.h>
#include <rtepack_multitype.h>
#include <vdisort.h>

#include <cmath>
#include <iostream>

#include "disort.h"
namespace {
constexpr Numeric tolerance = 2e-9;

struct two_stream_solution {
  Complex up;
  Complex down;
};

two_stream_solution analytic_two_stream(const Complex phase,
                                        const Numeric omega,
                                        const Numeric mu,
                                        const Numeric depth,
                                        const Numeric tau,
                                        const Complex top_down,
                                        const Complex bottom_up) {
  const Complex c = 0.5 * omega * phase;
  const Complex a = (1.0 - c) / mu;
  const Complex b = c / mu;
  const Complex k = std::sqrt(a * a - b * b);

  const auto propagator = [&](const Numeric x) {
    const Complex ch = std::cosh(k * x);
    const Complex sh = std::sinh(k * x) / k;
    return std::array<Complex, 4>{ch + a * sh, -b * sh, b * sh, ch - a * sh};
  };

  const auto    full   = propagator(depth);
  const Complex top_up = (bottom_up - full[1] * top_down) / full[0];
  const auto    at_tau = propagator(tau);
  return {.up = at_tau[0] * top_up + at_tau[1] * top_down, .down = at_tau[2] * top_up + at_tau[3] * top_down};
}

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

vdisort::main_data make_native_two_stream(const Numeric              depth,
                                          const Numeric              omega,
                                          vdisort::phase_matrix_data phase,
                                          rtepack::stokvec           top,
                                          rtepack::stokvec           bottom) {
  rtepack::stokvec_tensor3 up(2, 1, 1), down(2, 1, 1);
  up[vdisort::cosine_mode, 0, 0]   = {bottom.I(), bottom.Q(), 0.0, 0.0};
  down[vdisort::cosine_mode, 0, 0] = {top.I(), top.Q(), 0.0, 0.0};
  up[vdisort::sine_mode, 0, 0]     = {0.0, 0.0, bottom.U(), bottom.V()};
  down[vdisort::sine_mode, 0, 0]   = {0.0, 0.0, top.U(), top.V()};
  rtepack::stokvec_matrix source(1, 0);
  return vdisort::main_data(2,
                            1,
                            AscendingGrid{depth},
                            Vector{omega},
                            std::move(phase),
                            std::move(up),
                            std::move(down),
                            std::move(source),
                            {},
                            0.5,
                            rtepack::stokvec{},
                            0.0);
}

void test_analytic_iq_two_stream() {
  constexpr Numeric      depth = 0.8;
  constexpr Numeric      omega = 0.6;
  constexpr Numeric      mu    = 0.5;
  constexpr Numeric      p     = 0.7;
  constexpr Numeric      q     = 0.2;
  const rtepack::stokvec top{1.2, -0.25, 0.0, 0.0};
  const rtepack::stokvec bottom{0.4, 0.15, 0.0, 0.0};

  vdisort::phase_matrix_data phase(2, 1, 1, 2, 2, rtepack::muelmat{0.0});
  for (Index out = 0; out < 2; ++out) {
    for (Index in = 0; in < 2; ++in) {
      phase[vdisort::cosine_mode, 0, 0, out, in][0, 0] = p;
      phase[vdisort::cosine_mode, 0, 0, out, in][0, 1] = q;
      phase[vdisort::cosine_mode, 0, 0, out, in][1, 0] = q;
      phase[vdisort::cosine_mode, 0, 0, out, in][1, 1] = p;
    }
  }
  auto model = make_native_two_stream(depth, omega, std::move(phase), top, bottom);

  for (const Numeric tau : {0.0, 0.31, depth}) {
    const auto plus  = analytic_two_stream(p + q, omega, mu, depth, tau, top.I() + top.Q(), bottom.I() + bottom.Q());
    const auto minus = analytic_two_stream(p - q, omega, mu, depth, tau, top.I() - top.Q(), bottom.I() - bottom.Q());
    vdisort::u0_data field;
    model.u0(field, tau);
    expect_close(field.u0[0].I(), 0.5 * (plus.up.real() + minus.up.real()), "analytic I/Q upward I");
    expect_close(field.u0[0].Q(), 0.5 * (plus.up.real() - minus.up.real()), "analytic I/Q upward Q");
    expect_close(field.u0[1].I(), 0.5 * (plus.down.real() + minus.down.real()), "analytic I/Q downward I");
    expect_close(field.u0[1].Q(), 0.5 * (plus.down.real() - minus.down.real()), "analytic I/Q downward Q");
  }
}

void test_analytic_uv_two_stream() {
  constexpr Numeric      depth = 0.65;
  constexpr Numeric      omega = 0.55;
  constexpr Numeric      mu    = 0.5;
  constexpr Numeric      p     = 0.6;
  constexpr Numeric      q     = 0.23;
  const rtepack::stokvec top{0.0, 0.0, 0.3, -0.12};
  const rtepack::stokvec bottom{0.0, 0.0, -0.18, 0.21};

  vdisort::phase_matrix_data phase(2, 1, 1, 2, 2, rtepack::muelmat{0.0});
  for (Index out = 0; out < 2; ++out) {
    for (Index in = 0; in < 2; ++in) {
      phase[vdisort::sine_mode, 0, 0, out, in][2, 2] = p;
      phase[vdisort::sine_mode, 0, 0, out, in][2, 3] = q;
      phase[vdisort::sine_mode, 0, 0, out, in][3, 2] = -q;
      phase[vdisort::sine_mode, 0, 0, out, in][3, 3] = p;
    }
  }
  auto model = make_native_two_stream(depth, omega, std::move(phase), top, bottom);
  ARTS_USER_ERROR_IF(not model.has_complex_eigensolutions(), "Analytic U/V case did not exercise complex eigenvalues");

  for (const Numeric tau : {0.0, 0.27, depth}) {
    const auto expected = analytic_two_stream(
        Complex{p, -q}, omega, mu, depth, tau, Complex{top.U(), top.V()}, Complex{bottom.U(), bottom.V()});
    vdisort::u0_data field;
    model.u0(field, tau);
    expect_close(field.u0[0].U(), expected.up.real(), "analytic U/V upward U");
    expect_close(field.u0[0].V(), expected.up.imag(), "analytic U/V upward V");
    expect_close(field.u0[1].U(), expected.down.real(), "analytic U/V downward U");
    expect_close(field.u0[1].V(), expected.down.imag(), "analytic U/V downward V");
  }
}

void test_analytic_polarized_beam_two_stream() {
  constexpr Numeric      depth = 0.7;
  constexpr Numeric      omega = 0.45;
  constexpr Numeric      mu    = 0.5;
  constexpr Numeric      mu0   = 0.8;
  const rtepack::stokvec beam{1.0, 0.2, -0.15, 0.1};
  const rtepack::muelmat scatter{0.8, 0.1, 0.0, 0.0, -0.05, 0.6, 0.0, 0.0, 0.0, 0.0, 0.7, 0.12, 0.0, 0.0, -0.08, 0.5};

  vdisort::phase_matrix_data      phase(2, 1, 1, 2, 2, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data beam_phase(2, 1, 1, 2, rtepack::muelmat{0.0});
  for (Index i = 0; i < 2; ++i) {
    for (Index so = 0; so < 2; ++so)
      for (Index si = 0; si < 2; ++si) beam_phase[vdisort::cosine_mode, 0, 0, i][so, si] = scatter[so, si];
    for (Index so = 2; so < 4; ++so)
      for (Index si = 2; si < 4; ++si) beam_phase[vdisort::sine_mode, 0, 0, i][so, si] = scatter[so, si];
  }

  rtepack::stokvec_tensor3 boundary(2, 1, 1);
  rtepack::stokvec_matrix  source(1, 0);
  vdisort::main_data       model(2,
                                 1,
                                 AscendingGrid{depth},
                                 Vector{omega},
                                 std::move(phase),
                                 boundary,
                                 std::move(boundary),
                                 std::move(source),
                                 {},
                                 mu0,
                                 beam,
                                 0.0,
                                 std::move(beam_phase));
  const rtepack::stokvec   scattered = scatter * beam;

  const auto expected = [&](const Index stream, const Index stokes, const Numeric tau) {
    const Numeric inverse_mu = stream == 0 ? 1.0 / mu : -1.0 / mu;
    const Numeric rhs        = inverse_mu * omega * scattered[stokes] / (4.0 * Constant::pi);
    const Numeric b          = rhs / (inverse_mu + 1.0 / mu0);
    if (stream == 0) { return b * (std::exp(-tau / mu0) - std::exp(-depth / mu0 + (tau - depth) / mu)); }
    return b * (std::exp(-tau / mu0) - std::exp(-tau / mu));
  };

  for (const Numeric tau : {0.0, 0.29, depth}) {
    vdisort::u0_data field;
    model.u0(field, tau);
    for (Index i = 0; i < 2; ++i)
      for (Index s = 0; s < 4; ++s) expect_close(field.u0[i][s], expected(i, s, tau), "analytic polarized beam source");
  }
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
    expect_close(field.intensities[0][s], bottom[s] * std::exp((tau - depth) / mu), "upward absorbing Stokes");
    expect_close(field.intensities[1][s], top[s] * std::exp(-tau / mu), "downward absorbing Stokes");
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
      expect_close(vu.u0[i][0], su.u0[i], "scalar-limit I");
      expect_close(vu.u0[i][1], 0.0, "scalar-limit Q");
      expect_close(vu.u0[i][2], 0.0, "scalar-limit U");
      expect_close(vu.u0[i][3], 0.0, "scalar-limit V");
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
  Tensor7 phase(2, modes, 1, nquad, nquad, 4, 4, 0.0);
  phase[vdisort::cosine_mode, 0, 0, joker, joker, 0, 0] = 1.0;
  Tensor4 vector_up(2, modes, nquad / 2, 4, 0.0), vector_down(2, modes, nquad / 2, 4, 0.0);
  Tensor3 vector_source(1, 2, 4, 0.0);
  vector_source[0, 0, 0] = (1.0 - omega[0]) * source[0, 0];
  vector_source[0, 1, 0] = (1.0 - omega[0]) * source[0, 1];
  auto vector            = make_vdisort(
      nquad, tau, omega, std::move(phase), std::move(vector_up), std::move(vector_down), std::move(vector_source));
  for (const Numeric optical_depth : {0.0, 0.2, 0.8}) {
    disort::u_data  scalar_field;
    vdisort::u_data vector_field;
    scalar.u(scalar_field, optical_depth, 0.0);
    vector.u(vector_field, optical_depth, 0.0);
    for (Index stream = 0; stream < nquad; ++stream)
      expect_close(
          vector_field.intensities[stream].I(), scalar_field.intensities[stream], "linear-source scalar limit");
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
  Tensor7 phase(2, modes, 1, nquad, nquad, 4, 4, 0.0);
  phase[vdisort::cosine_mode, 0, 0, joker, joker, 0, 0] = 1.0;
  Tensor4 vector_up(2, modes, nquad / 2, 4, 0.0), vector_down(2, modes, nquad / 2, 4, 0.0);
  Tensor3 vector_source(1, 2, 4, 0.0);
  vector_source[0, 0, 0] = (1.0 - omega[0]) * source[0, 0];
  vector_source[0, 1, 0] = (1.0 - omega[0]) * source[0, 1];
  const auto reflection  = [](rtepack::muelmat_matrix_view value, const auto&, const auto&) {
    value = rtepack::muelmat{0.0};
    for (Index i = 0; i < value.nrows(); ++i)
      for (Index j = 0; j < value.ncols(); ++j) value[i, j][0, 0] = 2.0 * Constant::inv_pi;
  };
  std::vector<vdisort::BDRF> vector_brdf{{.cosine = {reflection}, .sine = {reflection}}};
  auto                       vector = make_vdisort(nquad,
                                                   tau,
                                                   omega,
                                                   std::move(phase),
                                                   std::move(vector_up),
                                                   std::move(vector_down),
                                                   std::move(vector_source),
                                                   Vector(4, 0.0),
                                                   {},
                                                   std::move(vector_brdf));
  for (const Numeric optical_depth : {0.2, 4.0, 8.0}) {
    disort::u_data  scalar_field;
    vdisort::u_data vector_field;
    scalar.u(scalar_field, optical_depth, 0.0);
    vector.u(vector_field, optical_depth, 0.0);
    for (Index stream = 0; stream < nquad; ++stream)
      expect_close(vector_field.intensities[stream].I(),
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
    expect_close(field.intensities[0][s], q[s] * (1.0 - std::exp(-(depth - tau) / mu)), "upward vector source");
    expect_close(field.intensities[1][s], q[s] * (1.0 - std::exp(-tau / mu)), "downward vector source");
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
  expect_close(field.u0[0][0], reflected_factor * std::exp(-depth / mu), "polarized BRDF I");
  expect_close(field.u0[0][1], 0.2 * reflected_factor * std::exp(-depth / mu), "polarized BRDF Q");
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
      const Numeric value = field.intensities[i][s];
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
  const auto native_combined      = vdisort::combine_phase_matrices(native_cosine, native_sine);
  expect_close(native_combined[vdisort::cosine_mode, 1, 0, 0, 0][0, 2], -5.0, "native Eq. 81 cosine sign");
  expect_close(native_combined[vdisort::sine_mode, 1, 0, 0, 0][0, 2], 5.0, "native Eq. 81 sine sign");
}

void test_spectral_phase_matrix_split() {
  rtepack::specmat_matrix spectral_result(2, 3, rtepack::specmat{Complex{0.0, 0.0}});
  spectral_result[1, 2][3, 1] = Complex{4.5, -0.75};
  const auto split            = vdisort::phase_matrix_fourier_split(spectral_result);
  ARTS_USER_ERROR_IF((split.cosine.shape() != spectral_result.shape() or split.sine.shape() != spectral_result.shape()),
                     "Splitting the phase matrix changed its shape");
  expect_close(split.cosine[1, 2][3, 1], 4.5, "spectral phase cosine coefficient");
  expect_close(split.sine[1, 2][3, 1], 0.75, "spectral phase sine coefficient");
}
}  // namespace

int main() try {
  test_analytic_iq_two_stream();
  test_analytic_uv_two_stream();
  test_analytic_polarized_beam_two_stream();
  test_polarized_absorption();
  test_scalar_limit();
  test_scalar_linear_source_limit();
  test_conservative_reflecting_source_limit();
  test_vector_source();
  test_polarized_brdf();
  test_complex_uv_eigenmodes();
  test_combined_matrix_transform();
  test_spectral_phase_matrix_split();
  std::cout << "vdisort tests passed\n";
  return 0;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}

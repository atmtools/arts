#include <arts_constants.h>
#include <rtepack_multitype.h>
#include <vdisort-brdf.h>
#include <vdisort.h>

#include <cmath>
#include <iostream>

#include "disort-brdf.h"
#include "disort.h"
#include "test-adapter.h"

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
                                Tensor3                    source        = {},
                                Vector                     beam          = Vector(4, 0.0),
                                Tensor6                    beam_phase    = {},
                                std::vector<vdisort::BDRF> brdf          = {},
                                Vector                     source_scale  = {},
                                Vector                     source_offset = {}) {
  if (source.size() == 0) source.resize(tau.size(), 0, 4);
  const Index nfourier = phase.shape()[1];
  return vdisort_test::make_solver(nquad,
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
                                   std::move(beam_phase),
                                   std::move(source_scale),
                                   std::move(source_offset));
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

  rtepack::stokvec_tensor3 boundary_up(2, 1, 1), boundary_down(2, 1, 1);
  rtepack::stokvec_matrix  source(1, 0);
  vdisort::main_data       model(2,
                                 1,
                                 AscendingGrid{depth},
                                 Vector{omega},
                                 std::move(phase),
                                 std::move(boundary_up),
                                 std::move(boundary_down),
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

  const Vector               user_mu{-1.0, -0.5, -0.1, 0.1, 0.5, 1.0};
  vdisort::phase_matrix_data user_phase(2, modes, 1, user_mu.size(), nquad, rtepack::muelmat{0.0});
  for (Index iu = 0; iu < static_cast<Index>(user_mu.size()); ++iu)
    for (Index j = 0; j < nquad; ++j) user_phase[vdisort::cosine_mode, 0, 0, iu, j][0, 0] = 1.0;
  vdisort::beam_phase_matrix_data user_beam_phase(2, modes, 1, user_mu.size(), rtepack::muelmat{0.0});
  for (Index iu = 0; iu < static_cast<Index>(user_mu.size()); ++iu)
    user_beam_phase[vdisort::cosine_mode, 0, 0, iu][0, 0] = 1.0;

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

    vdisort::user_u_data vector_user;
    disort::user_u_data  scalar_user;
    vector_model.u_user(vector_user, optical_depth, 0.0, user_mu, user_phase, user_beam_phase);
    scalar_model.u_user(scalar_user, optical_depth, 0.0, user_mu);
    for (Index i = 0; i < static_cast<Index>(user_mu.size()); ++i) {
      expect_close(vector_user.intensities[i].I(), scalar_user.intensities[i], "user-angle scalar-limit I");
      expect_close(vector_user.intensities[i].Q(), 0.0, "user-angle scalar-limit Q");
      expect_close(vector_user.intensities[i].U(), 0.0, "user-angle scalar-limit U");
      expect_close(vector_user.intensities[i].V(), 0.0, "user-angle scalar-limit V");
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
  vector_source[0, 0, 0] = source[0, 0];
  vector_source[0, 1, 0] = source[0, 1];
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

    disort::flux_data  scalar_flux_data;
    vdisort::flux_data vector_flux_data;
    const auto         scalar_flux = scalar.flux(scalar_flux_data, optical_depth);
    const auto         vector_flux = vector.flux(vector_flux_data, optical_depth);
    expect_close(vector_flux.up, scalar_flux.up, "linear-source upward flux");
    expect_close(vector_flux.down_diffuse, scalar_flux.down_diffuse, "linear-source downward flux");
    expect_close(vector_flux.down_direct, scalar_flux.down_direct, "linear-source direct flux");
    expect_close(vector_flux.dfdt, scalar_flux.dfdt, "linear-source DFDT");
  }
}

void test_affine_source_coordinate() {
  constexpr Index   nquad             = 4;
  constexpr Index   modes             = 1;
  constexpr Numeric depth             = 0.9;
  constexpr Numeric omega             = 0.35;
  constexpr Numeric intercept         = 0.8;
  constexpr Numeric slope             = -0.3;
  constexpr Numeric coordinate_scale  = 2.5;
  constexpr Numeric coordinate_offset = -0.4;

  const auto make_model = [=](const Numeric scale,
                              const Numeric offset,
                              const Numeric constant_coefficient,
                              const Numeric linear_coefficient) {
    Tensor7 phase(2, modes, 1, nquad, nquad, 4, 4, 0.0);
    phase[vdisort::cosine_mode, 0, 0, joker, joker, 0, 0] = 1.0;
    Tensor4 up(2, modes, nquad / 2, 4, 0.0), down(2, modes, nquad / 2, 4, 0.0);
    Tensor3 source(1, 2, 4, 0.0);
    source[0, 0, 0] = constant_coefficient;
    source[0, 1, 0] = linear_coefficient;
    return make_vdisort(nquad,
                        AscendingGrid{depth},
                        Vector{omega},
                        std::move(phase),
                        std::move(up),
                        std::move(down),
                        std::move(source),
                        Vector(4, 0.0),
                        {},
                        {},
                        Vector{scale},
                        Vector{offset});
  };

  auto physical = make_model(1.0, 0.0, intercept, slope);
  auto affine   = make_model(coordinate_scale,
                             coordinate_offset,
                             intercept - slope * coordinate_offset / coordinate_scale,
                             slope / coordinate_scale);

  for (const Numeric optical_depth : {0.0, 0.23, depth}) {
    vdisort::u_data physical_field;
    vdisort::u_data affine_field;
    physical.u(physical_field, optical_depth, 0.7);
    affine.u(affine_field, optical_depth, 0.7);
    for (Index stream = 0; stream < nquad; ++stream)
      for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
        expect_close(affine_field.intensities[stream][stokes],
                     physical_field.intensities[stream][stokes],
                     "affine source-coordinate quadrature field");
  }

  const Vector               user_mu{-0.75, 0.35};
  vdisort::phase_matrix_data user_phase(2, modes, 1, user_mu.size(), nquad, rtepack::muelmat{0.0});
  for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user)
    for (Index stream = 0; stream < nquad; ++stream) user_phase[vdisort::cosine_mode, 0, 0, user, stream][0, 0] = 1.0;
  const AscendingGrid      output_tau{0.0, 0.23, depth};
  const Vector             phi{0.0, 0.7};
  rtepack::stokvec_tensor3 physical_user(output_tau.size(), phi.size(), user_mu.size());
  rtepack::stokvec_tensor3 affine_user(output_tau.size(), phi.size(), user_mu.size());
  physical.ungridded_u_user(physical_user, output_tau, phi, user_mu, user_phase);
  affine.ungridded_u_user(affine_user, output_tau, phi, user_mu, user_phase);
  for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level)
    for (Index azimuth = 0; azimuth < static_cast<Index>(phi.size()); ++azimuth)
      for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user)
        for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
          expect_close(affine_user[level, azimuth, user][stokes],
                       physical_user[level, azimuth, user][stokes],
                       "affine source-coordinate user radiance");
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
  vector_source[0, 0, 0] = source[0, 0];
  vector_source[0, 1, 0] = source[0, 1];
  const auto reflection  = [](rtepack::muelmat_matrix_view value, const auto&, const auto&) {
    value = rtepack::muelmat{0.0};
    for (Index i = 0; i < value.nrows(); ++i)
      for (Index j = 0; j < value.ncols(); ++j) value[i, j][0, 0] = 2.0 * Constant::inv_pi;
  };
  std::vector<vdisort::BDRF> vector_brdf{
      {.cosine = {reflection}, .sine = {reflection}, .beam_cosine = {}, .beam_sine = {}}};
  auto vector = make_vdisort(nquad,
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
  const Vector    source_function{1.0, -0.2, 0.1, 0.05};
  source[0, 0] = source_function;

  const Numeric depth = 0.6;
  const Numeric omega = 0.3;
  auto          model = make_vdisort(
      nquad, AscendingGrid{depth}, Vector{omega}, std::move(phase), std::move(up), std::move(down), std::move(source));
  const Numeric   tau = 0.25;
  const Numeric   mu  = model.mu()[0];
  vdisort::u_data field;
  model.u(field, tau, 0.0);
  for (Index s = 0; s < 4; ++s) {
    expect_close(field.intensities[0][s],
                 (1.0 - omega) * source_function[s] * (1.0 - std::exp(-(depth - tau) / mu)),
                 "upward vector source function");
    expect_close(field.intensities[1][s],
                 (1.0 - omega) * source_function[s] * (1.0 - std::exp(-tau / mu)),
                 "downward vector source function");
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
  std::vector<vdisort::BDRF> brdf{{.cosine = {callback}, .sine = {callback}, .beam_cosine = {}, .beam_sine = {}}};

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

void test_bulk_quadrature_equivalence() {
  constexpr Index nquad  = 2;
  constexpr Index n      = nquad / 2;
  constexpr Index modes  = 2;
  constexpr Index layers = 2;

  Tensor7 phase(2, modes, layers, nquad, nquad, 4, 4, 0.0);
  Tensor4 up(2, modes, n, 4, 0.0), down(2, modes, n, 4, 0.0);

  up[vdisort::cosine_mode, 0, 0]   = Vector{1.1, -0.2, 0.0, 0.0};
  up[vdisort::sine_mode, 0, 0]     = Vector{0.0, 0.0, 0.3, -0.1};
  down[vdisort::cosine_mode, 0, 0] = Vector{0.7, 0.15, 0.0, 0.0};
  down[vdisort::sine_mode, 0, 0]   = Vector{0.0, 0.0, -0.25, 0.08};

  up[vdisort::cosine_mode, 1, 0]   = Vector{0.12, -0.04, 0.07, 0.02};
  up[vdisort::sine_mode, 1, 0]     = Vector{-0.03, 0.05, -0.06, 0.09};
  down[vdisort::cosine_mode, 1, 0] = Vector{-0.08, 0.02, 0.11, -0.04};
  down[vdisort::sine_mode, 1, 0]   = Vector{0.06, -0.01, 0.03, 0.05};

  Tensor6             beam_phase(2, modes, layers, nquad, 4, 4, 0.0);
  const AscendingGrid layer_tau{0.3, 0.9};
  auto                model = make_vdisort(nquad,
                                           layer_tau,
                                           Vector{0.0, 0.0},
                                           std::move(phase),
                                           std::move(up),
                                           std::move(down),
                                           {},
                                           Vector{0.7, 0.1, -0.05, 0.02},
                                           std::move(beam_phase));

  const Vector phi{0.0, 0.4, 1.7};
  Tensor4      gridded(layer_tau.size(), phi.size(), nquad, 4);
  model.gridded_u(gridded, phi);
  for (Index layer = 0; layer < static_cast<Index>(layer_tau.size()); ++layer)
    for (Index azimuth = 0; azimuth < static_cast<Index>(phi.size()); ++azimuth) {
      vdisort::u_data pointwise;
      model.u(pointwise, layer_tau[layer], phi[azimuth]);
      for (Index stream = 0; stream < nquad; ++stream)
        for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
          expect_close(gridded[layer, azimuth, stream, stokes],
                       pointwise.intensities[stream][stokes],
                       "gridded/pointwise polarized radiance");
    }

  const AscendingGrid output_tau{0.0, 0.17, 0.3, 0.56, 0.9};
  Tensor4             ungridded(output_tau.size(), phi.size(), nquad, 4);
  model.ungridded_u(ungridded, output_tau, phi);
  for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level)
    for (Index azimuth = 0; azimuth < static_cast<Index>(phi.size()); ++azimuth) {
      vdisort::u_data pointwise;
      model.u(pointwise, output_tau[level], phi[azimuth]);
      for (Index stream = 0; stream < nquad; ++stream)
        for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
          expect_close(ungridded[level, azimuth, stream, stokes],
                       pointwise.intensities[stream][stokes],
                       "ungridded/pointwise polarized radiance");
    }

  ARTS_USER_ERROR_IF(std::abs(ungridded[1, 1, 0, 2]) < 1e-12 or std::abs(ungridded[1, 1, 0, 3]) < 1e-12 or
                         std::abs(ungridded[1, 0, 0, 0] - ungridded[1, 1, 0, 0]) < 1e-12,
                     "Bulk equivalence model did not exercise U, V, and m>0 radiance");

  Vector gridded_up(layers), gridded_down(layers), gridded_direct(layers), gridded_dfdt(layers);
  model.gridded_flux(gridded_up, gridded_down, gridded_direct, gridded_dfdt);
  for (Index layer = 0; layer < layers; ++layer) {
    vdisort::flux_data scratch;
    const auto         pointwise = model.flux(scratch, layer_tau[layer]);
    expect_close(gridded_up[layer], pointwise.up, "gridded/pointwise upward flux");
    expect_close(gridded_down[layer], pointwise.down_diffuse, "gridded/pointwise downward flux");
    expect_close(gridded_direct[layer], pointwise.down_direct, "gridded/pointwise direct flux");
    expect_close(gridded_dfdt[layer], pointwise.dfdt, "gridded/pointwise DFDT");
  }

  Vector ungridded_up(output_tau.size()), ungridded_down(output_tau.size()), ungridded_direct(output_tau.size()),
      ungridded_dfdt(output_tau.size());
  model.ungridded_flux(ungridded_up, ungridded_down, ungridded_direct, ungridded_dfdt, output_tau);
  for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level) {
    vdisort::flux_data scratch;
    const auto         pointwise = model.flux(scratch, output_tau[level]);
    expect_close(ungridded_up[level], pointwise.up, "ungridded/pointwise upward flux");
    expect_close(ungridded_down[level], pointwise.down_diffuse, "ungridded/pointwise downward flux");
    expect_close(ungridded_direct[level], pointwise.down_direct, "ungridded/pointwise direct flux");
    expect_close(ungridded_dfdt[level], pointwise.dfdt, "ungridded/pointwise DFDT");
  }
}

void test_delta_m_correction_api_overlap() {
  constexpr Index        nquad          = 4;
  constexpr Index        nfourier       = 1;
  constexpr Numeric      physical_depth = 1.0;
  constexpr Numeric      omega          = 0.5;
  constexpr Numeric      fraction       = 0.2;
  constexpr Numeric      scale          = 1.0 - omega * fraction;
  const rtepack::stokvec beam{0.7, 0.1, -0.05, 0.02};

  Tensor7 phase(2, nfourier, 1, nquad, nquad, 4, 4, 0.0);
  Tensor4 up(2, nfourier, nquad / 2, 4, 0.0), down(2, nfourier, nquad / 2, 4, 0.0);
  Tensor6 beam_phase(2, nfourier, 1, nquad, 4, 4, 0.0);
  auto    model = make_vdisort(nquad,
                               AscendingGrid{scale * physical_depth},
                               Vector{0.0},
                               std::move(phase),
                               std::move(up),
                               std::move(down),
                               {},
                               Vector{beam.I(), beam.Q(), beam.U(), beam.V()},
                               std::move(beam_phase));

  rtepack::muelmat removed{0.0};
  removed[0, 0]               = 1.0;
  removed[1, 0]               = 0.2;
  removed[2, 0]               = -0.1;
  removed[3, 0]               = 0.05;
  removed[1, 1]               = 0.7;
  removed[2, 2]               = 0.6;
  removed[3, 3]               = 0.5;
  const auto   original_phase = [removed](Index, Numeric, Numeric, Numeric, Numeric) { return removed; };
  const auto   zero_phase     = [](Index, Numeric, Numeric, Numeric, Numeric) { return rtepack::muelmat{0.0}; };
  const Vector phi{0.0, 0.7};
  const vdisort::delta_m_correction_cache correction(AscendingGrid{physical_depth},
                                                     Vector{omega},
                                                     Vector{fraction},
                                                     0.5,
                                                     0.0,
                                                     beam,
                                                     Vector{model.mu()},
                                                     Vector{phi},
                                                     original_phase,
                                                     zero_phase,
                                                     original_phase,
                                                     8,
                                                     16);

  constexpr Numeric       physical_tau = 0.4;
  vdisort::u_data         quadrature;
  rtepack::stokvec_vector ims, tms;
  model.TMS(tms, physical_tau, 1, correction);
  model.IMS(ims, physical_tau, 1, correction);
  const auto separated_tms = tms;
  const auto separated_ims = ims;
  model.u_corr(quadrature, ims, tms, physical_tau, 1, correction);
  bool polarized = false;
  for (Index stream = 0; stream < nquad; ++stream) {
    for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes) {
      expect_close((separated_tms[stream] + separated_ims[stream])[stokes],
                   (tms[stream] + ims[stream])[stokes],
                   "separated TMS/IMS API");
      expect_close(quadrature.intensities[stream][stokes],
                   (tms[stream] + ims[stream])[stokes],
                   "quadrature corrected-radiance API");
    }
    polarized = polarized or std::abs(quadrature.intensities[stream].Q()) > 1e-12 or
                std::abs(quadrature.intensities[stream].U()) > 1e-12 or
                std::abs(quadrature.intensities[stream].V()) > 1e-12;
  }
  ARTS_USER_ERROR_IF(
      not polarized, "The polarized correction API produced no Q, U, or V signal: {:B,}", quadrature.intensities);

  vdisort::phase_matrix_data      user_phase(2, nfourier, 1, nquad, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data user_beam_phase(2, nfourier, 1, nquad, rtepack::muelmat{0.0});
  vdisort::user_u_data            user;
  rtepack::stokvec_vector         user_ims, user_tms;
  model.u_user_corr(user, user_ims, user_tms, physical_tau, 1, correction, user_phase, user_beam_phase);
  for (Index stream = 0; stream < nquad; ++stream)
    for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
      expect_close(user.intensities[stream][stokes],
                   quadrature.intensities[stream][stokes],
                   "user-angle corrected-radiance API");

  Tensor4                  corrected(1, phi.size(), nquad, vdisort::stokes_dimension);
  rtepack::stokvec_tensor3 gridded_tms(1, phi.size(), nquad), gridded_ims(1, phi.size(), nquad);
  model.gridded_TMS(gridded_tms, correction);
  model.gridded_IMS(gridded_ims, correction);
  model.gridded_u_corr(corrected, gridded_tms, gridded_ims, correction);
  for (Index p = 0; p < static_cast<Index>(phi.size()); ++p) {
    model.u_corr(quadrature, ims, tms, physical_depth, p, correction);
    for (Index stream = 0; stream < nquad; ++stream)
      for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes) {
        expect_close(
            corrected[0, p, stream, stokes], quadrature.intensities[stream][stokes], "gridded corrected-radiance API");
        expect_close(gridded_tms[0, p, stream][stokes], tms[stream][stokes], "gridded TMS API");
        expect_close(gridded_ims[0, p, stream][stokes], ims[stream][stokes], "gridded IMS API");
      }
  }
}

void test_combined_matrix_transform() {
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

void test_combined_surface_models() {
  constexpr Numeric lambertian_albedo = 0.8;
  constexpr Numeric fresnel_fraction  = 0.02;
  const Vector      outgoing{0.2, 0.7};
  const Vector      incoming{0.2, 0.7};

  const auto lambert = vdisort::brdf::lambertian_fourier_modes(lambertian_albedo, 3);
  for (Index mode = 0; mode < 3; ++mode) {
    rtepack::muelmat_matrix value(2, 2, rtepack::muelmat{0.0});
    lambert[mode](vdisort::cosine_mode, value, outgoing, incoming);
    for (Index i = 0; i < value.nrows(); ++i)
      for (Index j = 0; j < value.ncols(); ++j)
        for (Index so = 0; so < vdisort::stokes_dimension; ++so)
          for (Index si = 0; si < vdisort::stokes_dimension; ++si)
            expect_close(value[i, j][so, si],
                         mode == 0 and so == 0 and si == 0 ? 2.0 * lambertian_albedo * Constant::inv_pi : 0.0,
                         "depolarizing Lambertian mode");
  }

  const auto fresnel  = vdisort::brdf::fresnel_fourier_modes(1.5, 3);
  const auto combined = vdisort::brdf::fresnel_lambertian_fourier_modes(fresnel_fraction, lambertian_albedo, 1.5, 3);
  for (Index mode = 0; mode < 3; ++mode) {
    for (Index alpha = 0; alpha < 2; ++alpha) {
      rtepack::muelmat_matrix actual(2, 2, rtepack::muelmat{0.0});
      rtepack::muelmat_matrix specular(2, 2, rtepack::muelmat{0.0});
      rtepack::muelmat_matrix diffuse(2, 2, rtepack::muelmat{0.0});
      combined[mode](alpha, actual, outgoing, incoming);
      fresnel[mode](alpha, specular, outgoing, incoming);
      lambert[mode](alpha, diffuse, outgoing, incoming);
      for (Index i = 0; i < actual.nrows(); ++i)
        for (Index j = 0; j < actual.ncols(); ++j)
          for (Index so = 0; so < vdisort::stokes_dimension; ++so)
            for (Index si = 0; si < vdisort::stokes_dimension; ++si)
              expect_close(actual[i, j][so, si],
                           fresnel_fraction * specular[i, j][so, si] + (1.0 - fresnel_fraction) * diffuse[i, j][so, si],
                           "Fresnel/Lambertian mixture");
    }
  }

  constexpr Numeric cox_fraction = 0.35;
  const Vector      beam_mu{0.4};
  const auto        cox = vdisort::brdf::cox_munk_fourier_modes(5.0, 1.34, true, 2, 24);
  const auto        rough =
      vdisort::brdf::cox_munk_lambertian_fourier_modes(cox_fraction, lambertian_albedo, 5.0, 1.34, true, 2, 24);
  for (Index mode = 0; mode < 2; ++mode) {
    for (Index alpha = 0; alpha < 2; ++alpha) {
      rtepack::muelmat_matrix actual(2, 1, rtepack::muelmat{0.0});
      rtepack::muelmat_matrix specular(2, 1, rtepack::muelmat{0.0});
      rtepack::muelmat_matrix diffuse(2, 1, rtepack::muelmat{0.0});
      rough[mode].beam(alpha, actual, outgoing, beam_mu);
      cox[mode].beam(alpha, specular, outgoing, beam_mu);
      lambert[mode].beam(alpha, diffuse, outgoing, beam_mu);
      for (Index i = 0; i < actual.nrows(); ++i)
        for (Index so = 0; so < vdisort::stokes_dimension; ++so)
          for (Index si = 0; si < vdisort::stokes_dimension; ++si)
            expect_close(actual[i, 0][so, si],
                         cox_fraction * specular[i, 0][so, si] + (1.0 - cox_fraction) * diffuse[i, 0][so, si],
                         "Cox-Munk/Lambertian beam mixture");
    }
  }
}

void test_scalar_and_polarized_cox_munk_overlap() {
  for (const Complex refractive_index : {Complex{1.34, 0.0}, Complex{0.75, 0.0}, Complex{1.34, 0.08}})
    for (const bool shadowing : {false, true})
      for (const Numeric outgoing_mu : {0.15, 0.5, 0.9})
        for (const Numeric incoming_mu : {0.2, 0.65})
          for (const Numeric azimuth : {0.0, 0.7, 2.4}) {
            const disort::brdf::CoxMunk scalar{
                .wind_speed = 7.0, .refractive_index = refractive_index, .shadowing = shadowing};
            const vdisort::brdf::CoxMunk polarized{
                .wind_speed = 7.0, .refractive_index = refractive_index, .shadowing = shadowing};
            expect_close(polarized(outgoing_mu, incoming_mu, azimuth)[0, 0],
                         scalar(outgoing_mu, incoming_mu, azimuth),
                         "scalar/polarized Cox-Munk M00 overlap");
          }

  const auto absorbing_fresnel = vdisort::brdf::Fresnel{Complex{1.5, 0.1}}(0.5);
  ARTS_USER_ERROR_IF(std::abs(absorbing_fresnel[2, 3]) <= 1.0e-12,
                     "A complex refractive index produced no Fresnel U/V phase coupling");
}

void test_delta_m_preprocessing() {
  constexpr Index     nfourier = 2;
  constexpr Index     nlayers  = 2;
  constexpr Index     nquad    = 2;
  const AscendingGrid physical_tau{1.0, 3.0};
  const Vector        physical_omega{0.5, 0.8};
  const Vector        fraction{0.2, 0.25};

  vdisort::phase_matrix_data      original(2, nfourier, nlayers, nquad, nquad, rtepack::muelmat{0.0});
  vdisort::phase_matrix_data      removed(2, nfourier, nlayers, nquad, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data original_beam(2, nfourier, nlayers, nquad, rtepack::muelmat{0.0});
  vdisort::beam_phase_matrix_data removed_beam(2, nfourier, nlayers, nquad, rtepack::muelmat{0.0});
  for (Index alpha = 0; alpha < 2; ++alpha)
    for (Index mode = 0; mode < nfourier; ++mode)
      for (Index layer = 0; layer < nlayers; ++layer)
        for (Index out = 0; out < nquad; ++out) {
          const Numeric layer_value                    = static_cast<Numeric>(layer);
          original_beam[alpha, mode, layer, out][0, 0] = 4.0 + layer_value;
          removed_beam[alpha, mode, layer, out][0, 0]  = 1.0 + layer_value;
          for (Index in = 0; in < nquad; ++in) {
            original[alpha, mode, layer, out, in][0, 0] = 3.0 + layer_value;
            removed[alpha, mode, layer, out, in][0, 0]  = 0.5 + layer_value;
          }
        }

  const auto transport = vdisort::delta_m_preprocess(
      physical_tau, physical_omega, fraction, original, removed, original_beam, removed_beam);
  expect_close(transport.tau[0], 0.9, "delta-M first scaled boundary");
  expect_close(transport.tau[1], 2.5, "delta-M second scaled boundary");
  expect_close(transport.omega[0], 0.5 * 0.8 / 0.9, "delta-M first transport albedo");
  expect_close(transport.omega[1], 0.8 * 0.75 / 0.8, "delta-M second transport albedo");
  expect_close(transport.source_coordinate_scale[0], 1.0 / 0.9, "delta-M first source scale");
  expect_close(transport.source_coordinate_offset[0], 0.0, "delta-M first source offset");
  expect_close(transport.source_coordinate_scale[1], 1.0 / 0.8, "delta-M second source scale");
  expect_close(transport.source_coordinate_offset[1], -0.125, "delta-M second source offset");
  for (Index layer = 0; layer < nlayers; ++layer) {
    const Numeric layer_value = static_cast<Numeric>(layer);
    const Numeric expected_phase =
        (3.0 + layer_value - fraction[layer] * (0.5 + layer_value)) / (1.0 - fraction[layer]);
    const Numeric expected_beam = (4.0 + layer_value - fraction[layer] * (1.0 + layer_value)) / (1.0 - fraction[layer]);
    expect_close(transport.phase_matrix[1, 1, layer, 1, 0][0, 0], expected_phase, "delta-M transport phase");
    expect_close(transport.beam_phase_matrix[1, 1, layer, 1][0, 0], expected_beam, "delta-M transport beam phase");
  }
}

void test_depolarizing_surface_catalogue() {
  constexpr Index modes    = 3;
  constexpr Index nazimuth = 64;
  const Vector    outgoing{0.2, 0.7};
  const Vector    incoming{0.3, 0.8};

  const auto compare = [&](const vdisort::brdf::ScalarRawFunction& raw,
                           const std::vector<vdisort::BDRF>&       named,
                           const std::string_view                  name) {
    const auto generic = vdisort::brdf::depolarizing_fourier_modes(raw, modes, nazimuth);
    for (Index mode = 0; mode < modes; ++mode) {
      for (Index alpha = 0; alpha < 2; ++alpha) {
        rtepack::muelmat_matrix actual(outgoing.size(), incoming.size(), rtepack::muelmat{0.0});
        rtepack::muelmat_matrix expected(outgoing.size(), incoming.size(), rtepack::muelmat{0.0});
        named[mode](alpha, actual, outgoing, incoming);
        generic[mode](alpha, expected, outgoing, incoming);
        for (Index out = 0; out < actual.nrows(); ++out)
          for (Index in = 0; in < actual.ncols(); ++in)
            for (Index row = 0; row < vdisort::stokes_dimension; ++row)
              for (Index column = 0; column < vdisort::stokes_dimension; ++column)
                expect_close(actual[out, in][row, column], expected[out, in][row, column], name);
      }
    }
  };

  const disort::brdf::Hapke hapke;
  expect_close(hapke(0.4, 0.7, 0.3), disort_common::hapke_brdf(0.4, 0.7, 0.3, 1.0, 0.06, 0.6), "shared Hapke kernel");
  compare(hapke, vdisort::brdf::hapke_fourier_modes(1.0, 0.06, 0.6, modes, nazimuth), "depolarizing Hapke embedding");
  const disort::brdf::RPV rpv;
  expect_close(
      rpv(0.4, 0.7, 0.3), disort_common::rpv_brdf(0.4, 0.7, 0.3, 0.027, 0.647, -0.169, 0.1), "shared RPV kernel");
  compare(
      rpv, vdisort::brdf::rpv_fourier_modes(0.027, 0.647, -0.169, 0.1, modes, nazimuth), "depolarizing RPV embedding");
  const disort::brdf::RossLi ross_li;
  expect_close(ross_li(0.4, 0.7, 0.3),
               disort_common::ross_li_brdf(0.4, 0.7, 0.3, 0.091, 0.02, 0.01, 1.5 * Constant::pi / 180.0),
               "shared Ross-Li kernel");
  compare(ross_li,
          vdisort::brdf::ross_li_fourier_modes(0.091, 0.02, 0.01, 1.5 * Constant::pi / 180.0, modes, nazimuth),
          "depolarizing Ross-Li embedding");
}

void test_eigenvalue_direction_check() {
  constexpr Index nquad = 2;
  Tensor7         phase(2, 1, 1, nquad, nquad, 4, 4, 0.0);
  for (Index alpha = 0; alpha < 2; ++alpha)
    for (Index stokes = 0; stokes < 4; ++stokes) phase[alpha, 0, 0, 1, 1, stokes, stokes] = 10.0;
  Tensor4 up(2, 1, 1, 4, 0.0), down(2, 1, 1, 4, 0.0);

  try {
    static_cast<void>(
        make_vdisort(nquad, AscendingGrid{1.0}, Vector{0.5}, std::move(phase), std::move(up), std::move(down)));
  } catch (const std::exception& error) {
    ARTS_USER_ERROR_IF(
        std::string_view{error.what()}.find("cannot be split between the boundaries") == std::string_view::npos,
        "Unexpected eigenvalue-direction error: {}",
        error.what());
    return;
  }
  ARTS_USER_ERROR("VDISORT accepted an eigenspectrum without a valid propagation-direction split");
}
}  // namespace

int main() try {
  test_analytic_iq_two_stream();
  test_analytic_uv_two_stream();
  test_analytic_polarized_beam_two_stream();
  test_polarized_absorption();
  test_scalar_limit();
  test_scalar_linear_source_limit();
  test_affine_source_coordinate();
  test_conservative_reflecting_source_limit();
  test_vector_source();
  test_polarized_brdf();
  test_complex_uv_eigenmodes();
  test_bulk_quadrature_equivalence();
  test_delta_m_correction_api_overlap();
  test_combined_matrix_transform();
  test_spectral_phase_matrix_split();
  test_combined_surface_models();
  test_scalar_and_polarized_cox_munk_overlap();
  test_delta_m_preprocessing();
  test_depolarizing_surface_catalogue();
  test_eigenvalue_direction_check();
  std::cout << "vdisort tests passed\n";
  return 0;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}

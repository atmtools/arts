#include <arts_constants.h>
#include <disort.h>
#include <rtepack_multitype.h>
#include <vdisort.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string_view>

namespace {
constexpr Numeric top_downward  = 0.2;
constexpr Numeric bottom_upward = 0.7;
constexpr Numeric total_depth   = 3.0;

void expect_close(const std::string_view label,
                  const Numeric          actual,
                  const Numeric          expected,
                  const Numeric          tolerance = 5.0e-9) {
  if (not std::isfinite(actual) or std::abs(actual - expected) > tolerance * std::max(1.0, std::abs(expected)))
    throw std::runtime_error(std::format(
        "{}: expected {:.17g}, got {:.17g} (difference {:.3g})", label, expected, actual, actual - expected));
}

void expect_finite(const std::string_view label, const Numeric value) {
  if (not std::isfinite(value)) throw std::runtime_error(std::format("{} is not finite: {}", label, value));
}

Numeric sinhc(const Numeric x) {
  if (std::abs(x) > 1.0e-5) return std::sinh(x) / x;
  const Numeric x2 = x * x;
  return 1.0 + x2 * (1.0 / 6.0 + x2 * (1.0 / 120.0 + x2 / 5040.0));
}

struct two_stream_solution {
  Numeric up;
  Numeric down;
};

two_stream_solution analytic_two_stream(const Numeric tau,
                                        const Numeric omega,
                                        const Numeric phase,
                                        const Numeric top_down  = top_downward,
                                        const Numeric bottom_up = bottom_upward) {
  // VDISORT uses one double-Gauss node per hemisphere for NQuad=2, mu=1/2.
  constexpr Numeric mu = 0.5;
  const Numeric     c  = 0.5 * omega * phase;
  const Numeric     a  = (1.0 - c) / mu;
  const Numeric     b  = c / mu;
  const Numeric     k  = std::sqrt(std::max(0.0, a * a - b * b));

  const auto propagator = [&](const Numeric x) {
    const Numeric ch = std::cosh(k * x);
    const Numeric sh = x * sinhc(k * x);
    return std::array<Numeric, 4>{ch + a * sh, -b * sh, b * sh, ch - a * sh};
  };

  const auto    full   = propagator(total_depth);
  const Numeric top_up = (bottom_up - full[1] * top_down) / full[0];
  const auto    at_tau = propagator(tau);
  return {.up = at_tau[0] * top_up + at_tau[1] * top_down, .down = at_tau[2] * top_up + at_tau[3] * top_down};
}

vdisort::main_data make_model(const Numeric omega, const Numeric phase = 1.0, rtepack::stokvec_matrix source = {}) {
  vdisort::phase_matrix_data phase_matrix(2, 1, 1, 2, 2, rtepack::muelmat{0.0});
  for (Index outgoing = 0; outgoing < 2; ++outgoing)
    for (Index incoming = 0; incoming < 2; ++incoming)
      phase_matrix[vdisort::cosine_mode, 0, 0, outgoing, incoming][0, 0] = phase;

  rtepack::stokvec_tensor3 boundary_up(2, 1, 1), boundary_down(2, 1, 1);
  boundary_up                                   = rtepack::stokvec{};
  boundary_down                                 = rtepack::stokvec{};
  boundary_up[vdisort::cosine_mode, 0, 0].I()   = bottom_upward;
  boundary_down[vdisort::cosine_mode, 0, 0].I() = top_downward;

  if (source.size() == 0) source.resize(1, 0);
  return vdisort::main_data(2,
                            1,
                            AscendingGrid{total_depth},
                            Vector{omega},
                            std::move(phase_matrix),
                            std::move(boundary_up),
                            std::move(boundary_down),
                            std::move(source),
                            {},
                            0.5,
                            rtepack::stokvec{},
                            0.0);
}

vdisort::main_data make_layered_model(const AscendingGrid& depths, const Numeric omega) {
  const Index                layers = depths.size();
  vdisort::phase_matrix_data phase_matrix(2, 1, layers, 2, 2, rtepack::muelmat{0.0});
  for (Index layer = 0; layer < layers; ++layer)
    for (Index outgoing = 0; outgoing < 2; ++outgoing)
      for (Index incoming = 0; incoming < 2; ++incoming)
        phase_matrix[vdisort::cosine_mode, 0, layer, outgoing, incoming][0, 0] = 1.0;

  rtepack::stokvec_tensor3 boundary_up(2, 1, 1), boundary_down(2, 1, 1);
  boundary_up                                   = rtepack::stokvec{};
  boundary_down                                 = rtepack::stokvec{};
  boundary_up[vdisort::cosine_mode, 0, 0].I()   = bottom_upward;
  boundary_down[vdisort::cosine_mode, 0, 0].I() = top_downward;
  rtepack::stokvec_matrix source(layers, 0);
  return vdisort::main_data(2,
                            1,
                            AscendingGrid{depths},
                            Vector(layers, omega),
                            std::move(phase_matrix),
                            std::move(boundary_up),
                            std::move(boundary_down),
                            std::move(source),
                            {},
                            0.5,
                            rtepack::stokvec{},
                            0.0);
}

std::pair<vdisort::main_data, disort::main_data> make_multistream_scalar_pair(const Index nquad, const Numeric omega) {
  const Index                n = nquad / 2;
  vdisort::phase_matrix_data vector_phase(2, 1, 1, nquad, nquad, rtepack::muelmat{0.0});
  for (Index outgoing = 0; outgoing < nquad; ++outgoing)
    for (Index incoming = 0; incoming < nquad; ++incoming)
      vector_phase[vdisort::cosine_mode, 0, 0, outgoing, incoming][0, 0] = 1.0;

  rtepack::stokvec_tensor3 vector_up(2, 1, n), vector_down(2, 1, n);
  vector_up   = rtepack::stokvec{};
  vector_down = rtepack::stokvec{};
  for (Index stream = 0; stream < n; ++stream) {
    vector_up[vdisort::cosine_mode, 0, stream].I()   = bottom_upward;
    vector_down[vdisort::cosine_mode, 0, stream].I() = top_downward;
  }
  rtepack::stokvec_matrix vector_source(1, 0);
  vdisort::main_data      vector_model(nquad,
                                       1,
                                       AscendingGrid{total_depth},
                                       Vector{omega},
                                       std::move(vector_phase),
                                       std::move(vector_up),
                                       std::move(vector_down),
                                       std::move(vector_source),
                                       {},
                                       0.5,
                                       rtepack::stokvec{},
                                       0.0);

  Matrix scalar_moments(1, nquad, 0.0);
  scalar_moments[0, 0] = 1.0;
  Matrix            scalar_up(1, n, bottom_upward), scalar_down(1, n, top_downward);
  disort::main_data scalar_model(nquad,
                                 nquad,
                                 1,
                                 AscendingGrid{total_depth},
                                 Vector{omega},
                                 std::move(scalar_moments),
                                 std::move(scalar_up),
                                 std::move(scalar_down),
                                 Vector{0.0},
                                 Matrix(1, 0),
                                 {},
                                 1.0,
                                 0.0,
                                 0.0);
  return {std::move(vector_model), std::move(scalar_model)};
}

void check_quadrature_and_flux_sweep() {
  constexpr std::array epsilon{0.0, 1.0e-14, 1.0e-12, 1.0e-10, 1.0e-8, 2.0e-8};
  constexpr std::array depths{0.0, 0.3, 1.2, total_depth};

  for (const Numeric eps : epsilon) {
    const auto model = make_model(1.0 - eps);
    for (const Numeric tau : depths) {
      const auto expected = analytic_two_stream(tau, 1.0 - eps, 1.0);

      vdisort::u0_data u0;
      model.u0(u0, tau);
      expect_close("u0 upward I", u0.u0[0].I(), expected.up);
      expect_close("u0 downward I", u0.u0[1].I(), expected.down);
      for (Index stream = 0; stream < 2; ++stream)
        for (Index stokes = 1; stokes < vdisort::stokes_dimension; ++stokes)
          expect_close("u0 inactive Stokes component", u0.u0[stream][stokes], 0.0);

      vdisort::u_data u;
      model.u(u, tau, 0.37);
      expect_close("u upward I", u.intensities[0].I(), expected.up);
      expect_close("u downward I", u.intensities[1].I(), expected.down);

      vdisort::flux_data flux_data;
      const auto         flux = model.flux(flux_data, tau);
      expect_close("upward flux", flux.up, Constant::pi * expected.up);
      expect_close("downward flux", flux.down_diffuse, Constant::pi * expected.down);
      expect_close("direct flux", flux.down_direct, 0.0);
      expect_finite("DFDT", flux.dfdt);
    }
  }
}

Numeric analytic_conservative_user(const Numeric tau, const Numeric mu) {
  // At exact conservation the two quadrature-stream sum is linear and their
  // difference is constant.  The arbitrary-angle formal solution therefore
  // has an elementary linear-source particular solution.
  const Numeric q0             = (bottom_upward - top_downward) / (1.0 + total_depth);
  const Numeric s0             = q0 + 2.0 * top_downward;
  const auto    mean_intensity = [&](const Numeric x) { return 0.5 * s0 + q0 * x; };
  const Numeric boundary_tau   = mu > 0.0 ? total_depth : 0.0;
  const Numeric boundary       = mu > 0.0 ? bottom_upward : top_downward;
  return mean_intensity(tau) + mu * q0 +
         (boundary - mean_intensity(boundary_tau) - mu * q0) * std::exp((tau - boundary_tau) / mu);
}

void check_exact_conservative_user_angles() {
  const auto   model = make_model(1.0);
  const Vector user_mu{-0.91, -0.37, 0.19, 0.73};

  vdisort::phase_matrix_data user_phase(2, 1, 1, user_mu.size(), 2, rtepack::muelmat{0.0});
  for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user)
    for (Index incoming = 0; incoming < 2; ++incoming)
      user_phase[vdisort::cosine_mode, 0, 0, user, incoming][0, 0] = 1.0;

  constexpr std::array depths{0.0, 0.41, 1.7, total_depth};
  for (const Numeric tau : depths) {
    vdisort::user_u_data user_field;
    model.u_user(user_field, tau, 0.63, user_mu, user_phase);
    for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user) {
      expect_close("exact-conservative user I",
                   user_field.intensities[user].I(),
                   analytic_conservative_user(tau, user_mu[user]));
      for (Index stokes = 1; stokes < vdisort::stokes_dimension; ++stokes)
        expect_close("exact-conservative user inactive Stokes", user_field.intensities[user][stokes], 0.0);
    }
  }

  Vector output_depths(depths.size());
  stdr::copy(depths, output_depths.begin());
  const AscendingGrid      output_tau{std::move(output_depths)};
  const Vector             phi{0.0, 0.8};
  rtepack::stokvec_tensor3 bulk(output_tau.size(), phi.size(), user_mu.size());
  model.ungridded_u_user(bulk, output_tau, phi, user_mu, user_phase);
  for (Index level = 0; level < static_cast<Index>(output_tau.size()); ++level)
    for (Index azimuth = 0; azimuth < static_cast<Index>(phi.size()); ++azimuth)
      for (Index user = 0; user < static_cast<Index>(user_mu.size()); ++user)
        expect_close("bulk exact-conservative user I",
                     bulk[level, azimuth, user].I(),
                     analytic_conservative_user(output_tau[level], user_mu[user]));
}

void check_no_false_conservative_pair() {
  // omega=1 alone does not imply a conservative transport operator.  A phase
  // kernel with P_II=0.7 has ordinary nonzero eigenvalues and must retain the
  // ordinary modal path.
  constexpr Numeric phase = 0.7;
  const auto        model = make_model(1.0, phase);
  for (const Numeric tau : {0.0, 0.3, 1.2, total_depth}) {
    const auto       expected = analytic_two_stream(tau, 1.0, phase);
    vdisort::u0_data field;
    model.u0(field, tau);
    expect_close("nonconservative phase upward", field.u0[0].I(), expected.up, 5.0e-10);
    expect_close("nonconservative phase downward", field.u0[1].I(), expected.down, 5.0e-10);
  }
}

void check_simple_zero_is_not_conservative_pair() {
  // Both rows map an isotropic field back to itself, but the unequal column
  // sums do not conserve angularly integrated intensity.  The transport
  // matrix has one ordinary zero eigenvalue, not a size-two Jordan pair.
  vdisort::phase_matrix_data phase(2, 1, 1, 2, 2, rtepack::muelmat{0.0});
  phase[vdisort::cosine_mode, 0, 0, 0, 0][0, 0] = 2.0;
  phase[vdisort::cosine_mode, 0, 0, 0, 1][0, 0] = 0.0;
  phase[vdisort::cosine_mode, 0, 0, 1, 0][0, 0] = 1.0;
  phase[vdisort::cosine_mode, 0, 0, 1, 1][0, 0] = 1.0;

  rtepack::stokvec_tensor3 boundary_up(2, 1, 1), boundary_down(2, 1, 1);
  boundary_up                                   = rtepack::stokvec{};
  boundary_down                                 = rtepack::stokvec{};
  boundary_up[vdisort::cosine_mode, 0, 0].I()   = bottom_upward;
  boundary_down[vdisort::cosine_mode, 0, 0].I() = top_downward;
  rtepack::stokvec_matrix source(1, 0);
  vdisort::main_data      model(2,
                                1,
                                AscendingGrid{total_depth},
                                Vector{1.0},
                                std::move(phase),
                                std::move(boundary_up),
                                std::move(boundary_down),
                                std::move(source),
                                {},
                                0.5,
                                rtepack::stokvec{},
                                0.0);

  for (const Numeric tau : {0.0, 0.3, 1.2, total_depth}) {
    vdisort::u0_data field;
    model.u0(field, tau);
    expect_close("simple-zero upward", field.u0[0].I(), bottom_upward, 5.0e-10);
    expect_close("simple-zero downward",
                 field.u0[1].I(),
                 bottom_upward + (top_downward - bottom_upward) * std::exp(-tau),
                 5.0e-10);
  }
}

void check_polarized_generalized_pair() {
  // The intensity kernel conserves energy.  Its QI entries cancel for an
  // isotropic unpolarized field, but drive Q from the conserved intensity
  // current.  The exact generalized eigenvector must therefore be solved in
  // the full Stokes space rather than copied from scalar DISORT.
  constexpr Numeric r        = 0.35;
  constexpr Numeric d        = 0.2;
  constexpr Numeric top_q    = -0.1;
  constexpr Numeric bottom_q = 0.15;

  vdisort::phase_matrix_data phase(2, 1, 1, 2, 2, rtepack::muelmat{0.0});
  for (Index outgoing = 0; outgoing < 2; ++outgoing) {
    for (Index incoming = 0; incoming < 2; ++incoming) {
      auto& block = phase[vdisort::cosine_mode, 0, 0, outgoing, incoming];
      block[0, 0] = 1.0;
      block[1, 0] = incoming == 0 ? d : -d;
      block[1, 1] = r;
    }
  }

  rtepack::stokvec_tensor3 boundary_up(2, 1, 1), boundary_down(2, 1, 1);
  boundary_up                               = rtepack::stokvec{};
  boundary_down                             = rtepack::stokvec{};
  boundary_up[vdisort::cosine_mode, 0, 0]   = {bottom_upward, bottom_q, 0.0, 0.0};
  boundary_down[vdisort::cosine_mode, 0, 0] = {top_downward, top_q, 0.0, 0.0};
  rtepack::stokvec_matrix source(1, 0);
  vdisort::main_data      model(2,
                                1,
                                AscendingGrid{total_depth},
                                Vector{1.0},
                                std::move(phase),
                                std::move(boundary_up),
                                std::move(boundary_down),
                                std::move(source),
                                {},
                                0.5,
                                rtepack::stokvec{},
                                0.0);

  const Numeric intensity_current = (bottom_upward - top_downward) / (1.0 + total_depth);
  const Numeric q_particular      = d * intensity_current / (2.0 * (1.0 - r));
  for (const Numeric tau : {0.0, 0.3, 1.2, total_depth}) {
    const auto       expected_i = analytic_two_stream(tau, 1.0, 1.0);
    const auto       expected_q = analytic_two_stream(tau, 1.0, r, top_q - q_particular, bottom_q - q_particular);
    vdisort::u0_data field;
    model.u0(field, tau);
    expect_close("polarized conservative upward I", field.u0[0].I(), expected_i.up);
    expect_close("polarized conservative downward I", field.u0[1].I(), expected_i.down);
    expect_close("polarized conservative upward Q", field.u0[0].Q(), q_particular + expected_q.up);
    expect_close("polarized conservative downward Q", field.u0[1].Q(), q_particular + expected_q.down);
    for (Index stream = 0; stream < 2; ++stream) {
      expect_close("polarized conservative U", field.u0[stream].U(), 0.0);
      expect_close("polarized conservative V", field.u0[stream].V(), 0.0);
    }
  }
}

void check_exact_zero_emission_and_repeat_update() {
  rtepack::stokvec_matrix source(1, 2);
  source[0, 0]     = {123.0, -7.0, 5.0, 2.0};
  source[0, 1]     = {-45.0, 3.0, -9.0, 1.0};
  auto with_source = make_model(1.0, 1.0, std::move(source));
  auto no_source   = make_model(1.0);

  for (const Numeric tau : {0.0, 0.9, total_depth}) {
    vdisort::u0_data actual, expected;
    with_source.u0(actual, tau);
    no_source.u0(expected, tau);
    for (Index stream = 0; stream < 2; ++stream)
      for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
        expect_close("zero exact-conservative emission", actual.u0[stream][stokes], expected.u0[stream][stokes]);
  }

  vdisort::u0_data before, after;
  no_source.u0(before, 1.3);
  no_source.update_all();
  expect_close("omega retained by repeat update", no_source.omega()[0], 1.0, 0.0);
  no_source.u0(after, 1.3);
  for (Index stream = 0; stream < 2; ++stream)
    for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
      expect_close("repeat exact-conservative update", after.u0[stream][stokes], before.u0[stream][stokes]);
}

void check_layer_splitting() {
  for (const Numeric epsilon : {0.0, 1.0e-12}) {
    const auto one   = make_layered_model(AscendingGrid{total_depth}, 1.0 - epsilon);
    const auto split = make_layered_model(AscendingGrid{0.4, 1.1, total_depth}, 1.0 - epsilon);
    for (const Numeric tau : {0.0, 0.3, 0.8, 1.7, total_depth}) {
      vdisort::u0_data one_field, split_field;
      one.u0(one_field, tau);
      split.u0(split_field, tau);
      for (Index stream = 0; stream < 2; ++stream)
        for (Index stokes = 0; stokes < vdisort::stokes_dimension; ++stokes)
          expect_close(
              "conservative layer splitting", split_field.u0[stream][stokes], one_field.u0[stream][stokes], 1.0e-8);

      vdisort::flux_data one_scratch, split_scratch;
      const auto         one_flux   = one.flux(one_scratch, tau);
      const auto         split_flux = split.flux(split_scratch, tau);
      expect_close("conservative split upward flux", split_flux.up, one_flux.up, 1.0e-8);
      expect_close("conservative split downward flux", split_flux.down_diffuse, one_flux.down_diffuse, 1.0e-8);
      expect_close("conservative split DFDT", split_flux.dfdt, one_flux.dfdt, 1.0e-8);
    }
  }
}

void check_multistream_scalar_limit() {
  for (const Index nquad : {4, 8}) {
    for (const Numeric epsilon : {0.0, 1.0e-14, 1.0e-10, 1.0e-8, 2.0e-8}) {
      const auto [vector_model, scalar_model] = make_multistream_scalar_pair(nquad, 1.0 - epsilon);
      for (const Numeric tau : {0.0, 0.3, 1.2, total_depth}) {
        vdisort::u0_data vector_field;
        disort::u0_data  scalar_field;
        vector_model.u0(vector_field, tau);
        scalar_model.u0(scalar_field, tau);
        for (Index stream = 0; stream < nquad; ++stream) {
          expect_close("multistream scalar-limit I", vector_field.u0[stream].I(), scalar_field.u0[stream], 2.0e-8);
          for (Index stokes = 1; stokes < vdisort::stokes_dimension; ++stokes)
            expect_close("multistream scalar-limit inactive Stokes", vector_field.u0[stream][stokes], 0.0);
        }

        vdisort::flux_data vector_scratch;
        disort::flux_data  scalar_scratch;
        const auto         vector_flux = vector_model.flux(vector_scratch, tau);
        const auto         scalar_flux = scalar_model.flux(scalar_scratch, tau);
        expect_close("multistream scalar-limit upward flux", vector_flux.up, scalar_flux.up, 2.0e-8);
        expect_close(
            "multistream scalar-limit downward flux", vector_flux.down_diffuse, scalar_flux.down_diffuse, 2.0e-8);
        expect_close("multistream scalar-limit DFDT", vector_flux.dfdt, scalar_flux.dfdt, 2.0e-8);
      }
    }
  }
}
}  // namespace

int main() try {
  check_quadrature_and_flux_sweep();
  check_exact_conservative_user_angles();
  check_no_false_conservative_pair();
  check_simple_zero_is_not_conservative_pair();
  check_polarized_generalized_pair();
  check_exact_zero_emission_and_repeat_update();
  check_layer_splitting();
  check_multistream_scalar_limit();
  return 0;
} catch (const std::exception& error) {
  std::cerr << "Conservative VDISORT test failed: " << error.what() << '\n';
  return 1;
}

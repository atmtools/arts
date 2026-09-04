#include <rtepack_rtestep.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
constexpr Numeric tolerance = 1e-15;

bool close(const rtepack::stokvec& lhs, const rtepack::stokvec& rhs) {
  return std::max({std::abs(lhs.I() - rhs.I()),
                   std::abs(lhs.Q() - rhs.Q()),
                   std::abs(lhs.U() - rhs.U()),
                   std::abs(lhs.V() - rhs.V())}) < tolerance;
}

bool test_current_rte_emission() {
  constexpr Numeric transmission = 0.25;
  constexpr Numeric evolution    = 0.60;

  TransmittanceMatrix tramat;
  tramat.option = TransmittanceOption::linsrc;
  tramat.T.resize(1, 2);
  tramat.L.resize(1, 2);
  tramat.P.resize(1, 2);
  tramat.dT.resize(2, 1, 2, 1);
  tramat.dL.resize(2, 1, 2, 1);
  tramat.T       = rtepack::muelmat::id();
  tramat.L       = rtepack::muelmat::id();
  tramat.P       = rtepack::muelmat::id();
  tramat.dT      = rtepack::muelmat::constant(0.0);
  tramat.dL      = rtepack::muelmat::constant(0.0);
  tramat.T[0, 1] = transmission;
  tramat.L[0, 1] = evolution;

  SourceVector source;
  source.J.resize(1, 2);
  source.dJ.resize(1, 2, 1);
  source.J = rtepack::stokvec{};
  const rtepack::stokvec observer_deriv{2.0, 3.0, 5.0, 7.0};
  const rtepack::stokvec background_deriv{11.0, 13.0, 17.0, 19.0};
  source.dJ[0, 0, 0] = observer_deriv;
  source.dJ[0, 1, 0] = background_deriv;

  rtepack::stokvec_vector  radiance(1, rtepack::stokvec{});
  rtepack::stokvec_tensor3 jacobian(1, 2, 1, rtepack::stokvec{});
  rtepack::rte_emission(radiance, jacobian, tramat, source);

  const rtepack::stokvec expected_observer   = (1.0 - evolution) * observer_deriv;
  const rtepack::stokvec expected_background = (evolution - transmission) * background_deriv;
  if (not close(jacobian[0, 0, 0], expected_observer) or not close(jacobian[0, 1, 0], expected_background)) {
    std::cerr << "rte_emission assigns linear-source derivatives to the wrong endpoint\n";
    return false;
  }
  return true;
}

bool test_legacy_rte_emission() {
  constexpr Numeric transmission = 0.25;
  constexpr Numeric evolution    = 0.60;

  std::vector<rtepack::muelmat_vector>  Ts(2), Ls(2), Ps(2);
  std::vector<rtepack::muelmat_tensor3> dTs(2), dLs(2);
  std::vector<rtepack::stokvec_vector>  Js(2);
  std::vector<rtepack::stokvec_matrix>  dJs(2);
  for (Size i = 0; i < 2; ++i) {
    Ts[i].resize(1);
    Ls[i].resize(1);
    Ps[i].resize(1);
    Ts[i] = Ls[i] = Ps[i] = rtepack::muelmat::id();
    dTs[i].resize(2, 1, 1);
    dLs[i].resize(2, 1, 1);
    dTs[i] = dLs[i] = rtepack::muelmat::constant(0.0);
    Js[i].resize(1);
    Js[i] = rtepack::stokvec{};
    dJs[i].resize(1, 1);
  }
  Ts[1][0] = transmission;
  Ls[1][0] = evolution;

  const rtepack::stokvec observer_deriv{2.0, 3.0, 5.0, 7.0};
  const rtepack::stokvec background_deriv{11.0, 13.0, 17.0, 19.0};
  dJs[0][0, 0] = observer_deriv;
  dJs[1][0, 0] = background_deriv;

  rtepack::stokvec_vector              radiance;
  std::vector<rtepack::stokvec_matrix> jacobian;
  const rtepack::stokvec_vector        background(1, rtepack::stokvec{});
  rtepack::two_level_linear_evolution_step_by_step_full(radiance, jacobian, Ts, Ls, Ps, dTs, dLs, Js, dJs, background);

  const rtepack::stokvec expected_observer   = (1.0 - evolution) * observer_deriv;
  const rtepack::stokvec expected_background = (evolution - transmission) * background_deriv;
  if (not close(jacobian[0][0, 0], expected_observer) or not close(jacobian[1][0, 0], expected_background)) {
    std::cerr << "legacy linear-evolution step assigns source derivatives to the wrong endpoint\n";
    return false;
  }
  return true;
}

bool test_thin_layer_remainders() {
  constexpr Numeric optical_depth = 1e-18;
  constexpr Numeric source_scale  = 1e20;

  std::vector<rtepack::propmat>        K(2, rtepack::propmat{1.0});
  std::vector<rtepack::propmat_vector> dK(2);
  for (auto& value : dK) value.resize(1);
  Vector  distance(1, optical_depth);
  Tensor3 dr(2, 1, 1, 0.0);

  SourceVector source;
  source.J.resize(1, 2);
  source.dJ.resize(1, 2, 1);
  source.J           = rtepack::stokvec{source_scale};
  source.dJ          = rtepack::stokvec{};
  source.dJ[0, 0, 0] = rtepack::stokvec{source_scale};
  source.dJ[0, 1, 0] = rtepack::stokvec{source_scale};

  const Numeric expected_radiance = -std::expm1(-optical_depth) * source_scale;
  // Evaluate phi_1(-tau)-1 without a subtraction from one; its first omitted
  // term is 1e-56 after scaling here.
  const Numeric stable_phi_m1                = -0.5 * optical_depth + optical_depth * optical_depth / 6.0;
  const Numeric expected_observer_jacobian   = -stable_phi_m1 * source_scale;
  const Numeric expected_background_jacobian = (stable_phi_m1 - std::expm1(-optical_depth)) * source_scale;

  constexpr std::array options{TransmittanceOption::constant,
                               TransmittanceOption::linsrc,
                               TransmittanceOption::linprop,
                               TransmittanceOption::magop,
                               TransmittanceOption::magop_linsrc};
  for (const auto option : options) {
    TransmittanceMatrix tramat;
    tramat.init(K, dK, distance, dr, option);

    rtepack::stokvec_vector  radiance(1, rtepack::stokvec{});
    rtepack::stokvec_tensor3 jacobian(1, 2, 1, rtepack::stokvec{});
    rtepack::rte_emission(radiance, jacobian, tramat, source);

    const Numeric observer_expected   = option == TransmittanceOption::constant or option == TransmittanceOption::magop
                                            ? 0.5 * expected_radiance
                                            : expected_observer_jacobian;
    const Numeric background_expected = option == TransmittanceOption::constant or option == TransmittanceOption::magop
                                            ? 0.5 * expected_radiance
                                            : expected_background_jacobian;
    const auto    near                = [](const Numeric lhs, const Numeric rhs) {
      return std::abs(lhs - rhs) <= 2e-14 * std::max(1.0, std::abs(rhs));
    };
    if (not near(radiance[0].I(), expected_radiance) or not near(jacobian[0, 0, 0].I(), observer_expected) or
        not near(jacobian[0, 1, 0].I(), background_expected)) {
      std::cerr << "Thin-layer source cancellation in rte_emission for option " << static_cast<int>(option) << '\n';
      return false;
    }

    rtepack::stokvec_matrix path(1, 2, rtepack::stokvec{});
    rtepack::rte_emission_path(path, tramat, source);
    if (not near(path[0, 0].I(), expected_radiance)) {
      std::cerr << "Thin-layer source cancellation in rte_emission_path for option " << static_cast<int>(option)
                << '\n';
      return false;
    }
  }

  return true;
}
}  // namespace

int main() {
  return test_current_rte_emission() and test_legacy_rte_emission() and test_thin_layer_remainders() ? 0 : 1;
}

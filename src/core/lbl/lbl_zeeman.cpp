/**
 * @file   zeemandata.cc
 * @author Richard Larsson <larsson (at) mps.mpg.de>
 * @date   2018-04-06
 * 
 * @brief Implementations of Zeeman modeling
 * 
 * This file serves to implement Zeeman splitting
 * using various up-to-speed methods
 */

#include "lbl_zeeman.h"

#include "debug.h"
#include "double_imanip.h"
#include "quantum_numbers.h"
#include "wigner_functions.h"

namespace lbl::zeeman {
constexpr Numeric get_lande_spin_constant(
    const Species::Species species) noexcept {
  switch (species) {
    case Species::fromShortName("O2"):
      return 2.002064;
    case Species::fromShortName("NO"):
      return 2.00071;
    case Species::fromShortName("OH"):
      return 2.00089;
    case Species::fromShortName("ClO"):
      return 2.00072;
    case Species::fromShortName("SO"):
      return 2.002106;
    default:
      break;
  }
  return 2.00231930436182;
}

constexpr Numeric get_lande_lambda_constant() noexcept { return 1.0; }

data SimpleG(const Quantum::Number::ValueList& qns,
             const Numeric& GS,
             const Numeric& GL) noexcept {
  if (Quantum::Number::vamdcCheck(qns, Quantum::Number::VAMDC::hunda) and
      qns.has(QuantumNumberType::Omega,
              QuantumNumberType::J,
              QuantumNumberType::Lambda,
              QuantumNumberType::S)) {
    auto& Omega = qns[QuantumNumberType::Omega];
    auto& J = qns[QuantumNumberType::J];
    auto& Lambda = qns[QuantumNumberType::Lambda];
    auto& S = qns[QuantumNumberType::S];
    return {
        .gu = SimpleGCaseA(Omega.upp(), J.upp(), Lambda.upp(), S.upp(), GS, GL),
        .gl =
            SimpleGCaseA(Omega.low(), J.low(), Lambda.low(), S.low(), GS, GL)};
  }

  if (Quantum::Number::vamdcCheck(qns, Quantum::Number::VAMDC::hundb) and
      qns.has(QuantumNumberType::N,
              QuantumNumberType::J,
              QuantumNumberType::Lambda,
              QuantumNumberType::S)) {
    auto& N = qns[QuantumNumberType::N];
    auto& J = qns[QuantumNumberType::J];
    auto& Lambda = qns[QuantumNumberType::Lambda];
    auto& S = qns[QuantumNumberType::S];
    return {
        .gu = SimpleGCaseB(N.upp(), J.upp(), Lambda.upp(), S.upp(), GS, GL),
        .gl = SimpleGCaseB(N.low(), J.low(), Lambda.low(), S.low(), GS, GL)};
  }

  return {};
}

data GetSimpleModel(const QuantumIdentifier& qid) ARTS_NOEXCEPT {
  const Numeric GS = get_lande_spin_constant(qid.Species());
  const Numeric GL = get_lande_lambda_constant();
  return SimpleG(qid.val, GS, GL);
}

Numeric case_b_g_coefficient_o2(Rational J,
                                Rational N,
                                Numeric GS,
                                Numeric GR,
                                Numeric GLE,
                                Numeric B,
                                Numeric D,
                                Numeric H,
                                Numeric gB,
                                Numeric gD,
                                Numeric gH,
                                Numeric lB,
                                Numeric lD,
                                Numeric lH) {
  using Math::pow2, Math::pow3;
  using std::atan2, std::cos, std::sin;

  if (J.isUndefined() or N.isUndefined()) return NAN;
  if (J == 0) return 0;

  auto nom = (lB + lD * (J * J + J + 1) + lH * pow2(J * J + J + 1)) *
             (2 * sqrt(J * J + J) / (2 * J + 1));

  auto denom =
      B * J * (J - 1) - D * pow2(J * (J - 1)) + H * pow3(J * (J - 1)) +
      (gB + gD * J * (J - 1) + gH * pow2(J * (J - 1))) * (J - 1) +
      (lB + lD * J * (J - 1) + lH * pow2(J * (J - 1))) *
          (2. / 3. - 2 * J / (2 * J + 1)) -
      (B * (J + 2) * (J + 1) - D * pow2((J + 2) * (J + 1)) +
       H * pow3((J + 2) * (J + 1)) -
       (gB + gD * (J + 2) * (J + 1) + gH * pow2((J + 2) * (J + 1))) * (J + 2) +
       (lB + lD * (J + 2) * (J + 1) + lH * pow2((J + 2) * (J + 1))) *
           (2. / 3. - 2 * (J + 1) / (2 * J + 1)));

  auto phi = atan2(2 * nom, denom) / 2;

  if (J == N) return (GS + GR) / (J * (J + 1)) - GR;
  if (J < N)
    return (GS + GR) * (pow2(cos(phi)) / J - pow2(sin(phi)) / (J + 1)) +
           2 * GLE * cos(2 * phi) / (2 * J + 1) - GR;
  return (GS + GR) * (pow2(sin(phi)) / J - pow2(cos(phi)) / (J + 1)) -
         2 * GLE * cos(2 * phi) / (2 * J + 1) - GR;
}

constexpr Numeric closed_shell_trilinear(Rational k,
                                         Rational j,
                                         Numeric gperp,
                                         Numeric gpara) {
  using Math::pow2;
  if (k.isUndefined() or j.isUndefined() or j == 0) return 0;
  return gperp + (gperp + gpara) * (pow2(k) / (j * (j + 1)));
}

data GetAdvancedModel(const QuantumIdentifier& qid) ARTS_NOEXCEPT {
  if (qid.Isotopologue() == "O2-66") {
    if (qid.val.has(QuantumNumberType::J,
                    QuantumNumberType::N,
                    QuantumNumberType::v1)) {
      if (qid.val[QuantumNumberType::v1].low() == 0 and
          qid.val[QuantumNumberType::v1].upp() == 0) {
        constexpr Numeric GS = 2.002084;
        constexpr Numeric GLE = 2.77e-3;
        constexpr Numeric GR = -1.16e-4;
        constexpr Numeric B = 43100.44276e6;
        constexpr Numeric D = 145.1271e3;
        constexpr Numeric H = 49e-3;
        constexpr Numeric lB = 59501.3438e6;
        constexpr Numeric lD = 58.3680e3;
        constexpr Numeric lH = 290.8e-3;
        constexpr Numeric gB = -252.58634e6;
        constexpr Numeric gD = -243.42;
        constexpr Numeric gH = -1.46e-3;

        const auto& J = qid.val[QuantumNumberType::J];
        const auto& N = qid.val[QuantumNumberType::N];
        auto JU = J.upp();
        auto NU = N.upp();
        auto JL = J.low();
        auto NL = N.low();

        Numeric gu = case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        Numeric gl = case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {.gu = gu, .gl = gl};
      }
    }
  } else if (qid.Isotopologue() == "O2-68") {
    if (qid.val.has(QuantumNumberType::J,
                    QuantumNumberType::N,
                    QuantumNumberType::v1)) {
      if (qid.val[QuantumNumberType::v1].low() == 0 and
          qid.val[QuantumNumberType::v1].upp() == 0) {
        constexpr Numeric GS = 2.002025;
        constexpr Numeric GLE = 2.813e-3;
        constexpr Numeric GR = -1.26e-4;
        constexpr Numeric B = 40707.38657e6;
        constexpr Numeric D = 129.4142e3;
        constexpr Numeric H = 0;
        constexpr Numeric lB = 59499.0375e6;
        constexpr Numeric lD = 54.9777e3;
        constexpr Numeric lH = 272.1e-3;
        constexpr Numeric gB = -238.51530e6;
        constexpr Numeric gD = -217.77;
        constexpr Numeric gH = -1.305e-3;

        const auto& J = qid.val[QuantumNumberType::J];
        const auto& N = qid.val[QuantumNumberType::N];
        auto JU = J.upp();
        auto NU = N.upp();
        auto JL = J.low();
        auto NL = N.low();

        Numeric gu = case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        Numeric gl = case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {.gu = gu, .gl = gl};
      }
    }
  } else if (qid.Isotopologue() == "CO-26") {
    constexpr Numeric gperp =
        -0.2689 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    return {.gu = gperp, .gl = gperp};
  } else if (qid.Isotopologue() == "OCS-622") {
    constexpr Numeric gperp =
        -.02889 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        0 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    if (qid.val.has(QuantumNumberType::J, QuantumNumberType::Ka)) {
      const auto& J = qid.val[QuantumNumberType::J];
      const auto& Ka = qid.val[QuantumNumberType::Ka];
      auto JU = J.upp();
      auto KU = Ka.upp();
      auto JL = J.low();
      auto KL = Ka.low();

      return {.gu = closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.Isotopologue() == "OCS-624") {
    constexpr Numeric gperp =
        -.0285 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        -.061 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    if (qid.val.has(QuantumNumberType::J, QuantumNumberType::Ka)) {
      const auto& J = qid.val[QuantumNumberType::J];
      const auto& Ka = qid.val[QuantumNumberType::Ka];
      auto JU = J.upp();
      auto KU = Ka.upp();
      auto JL = J.low();
      auto KL = Ka.low();

      return {.gu = closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.Isotopologue() == "CO2-626") {
    constexpr Numeric gperp =
        -.05508 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        0 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    if (qid.val.has(QuantumNumberType::J, QuantumNumberType::Ka)) {
      const auto& J = qid.val[QuantumNumberType::J];
      const auto& Ka = qid.val[QuantumNumberType::Ka];
      auto JU = J.upp();
      auto KU = Ka.upp();
      auto JL = J.low();
      auto KL = Ka.low();

      return {.gu = closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  }

  // Set to zero otherwise as we practically say "there's no Zeeman effect" then
  return {.gu = 0, .gl = 0};
}

model::model(const QuantumIdentifier& qid) noexcept {
  model m = GetAdvancedModel(qid);
  if (m.empty()) m = GetSimpleModel(qid);
  *this = m;
}

Numeric model::Strength(Rational Ju,
                        Rational Jl,
                        pol type,
                        Index n) const ARTS_NOEXCEPT {
  ARTS_ASSERT(type not_eq Polarization::None);
  using Math::pow2;

  auto ml = Ml(Ju, Jl, type, n);
  auto mu = Mu(Ju, Jl, type, n);
  auto dm = Rational(dM(type));
  return type == pol::no
             ? 0.0
             : polarization_factor(type) *
                   pow2(wigner3j(Jl, Rational(1), Ju, ml, -dm, -mu));
}

Numeric model::Strength(const QuantumNumberValueList& qn,
                        pol type,
                        Index n) const ARTS_NOEXCEPT {
  if (type == pol::no) return 1.0;

  const auto& J = qn[QuantumNumberType::J];

  return Strength(J.upp(), J.low(), type, n);
}

Numeric model::Splitting(const QuantumNumberValueList& qn,
                         pol type,
                         Index n) const noexcept {
  if (type == pol::no) return 0.0;

  const auto& J = qn[QuantumNumberType::J];

  return Splitting(J.upp(), J.low(), type, n);
}

Index model::size(const QuantumNumberValueList& qn, pol type) const noexcept {
  if (type == pol::no) return 1;

  const auto& J = qn[QuantumNumberType::J];

  return zeeman::size(J.upp(), J.low(), type);
}

std::ostream& operator<<(std::ostream& os, const model& m) {
  os << m.mdata.gu << ' ' << m.mdata.gl;
  return os;
}

std::istream& operator>>(std::istream& is, model& m) {
  is >> double_imanip() >> m.mdata.gu >> m.mdata.gl;
  return is;
}
}  // namespace lbl::zeeman

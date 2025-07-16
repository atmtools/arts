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

#include <arts_constexpr_math.h>
#include <double_imanip.h>
#include <wigner_functions.h>

#include <utility>

namespace {
constexpr Numeric get_lande_spin_constant(const SpeciesEnum species) noexcept {
  switch (species) {
    case to<SpeciesEnum>("O2"):  return 2.002064;
    case to<SpeciesEnum>("NO"):  return 2.00071;
    case to<SpeciesEnum>("OH"):  return 2.00089;
    case to<SpeciesEnum>("ClO"): return 2.00072;
    case to<SpeciesEnum>("SO"):  return 2.002106;
    default:                     break;
  }
  return 2.00231930436182;
}

constexpr Numeric get_lande_lambda_constant() noexcept { return 1.0; }

lbl::zeeman::data SimpleG(const QuantumState& qns,
                          const Numeric& GS,
                          const Numeric& GL) noexcept {
  if (Quantum::vamdcCheck(qns, Quantum::VAMDC::hunda) and
      qns.contains(QuantumNumberType::Omega) and
      qns.contains(QuantumNumberType::J) and
      qns.contains(QuantumNumberType::Lambda) and
      qns.contains(QuantumNumberType::S)) {
    auto& Omega  = qns.at(QuantumNumberType::Omega);
    auto& J      = qns.at(QuantumNumberType::J);
    auto& Lambda = qns.at(QuantumNumberType::Lambda);
    auto& S      = qns.at(QuantumNumberType::S);
    return {.gu = lbl::zeeman::SimpleGCaseA(Omega.upper.get<Rational>(),
                                            J.upper.get<Rational>(),
                                            Lambda.upper.get<Rational>(),
                                            S.upper.get<Rational>(),
                                            GS,
                                            GL),
            .gl = lbl::zeeman::SimpleGCaseA(Omega.lower.get<Rational>(),
                                            J.lower.get<Rational>(),
                                            Lambda.lower.get<Rational>(),
                                            S.lower.get<Rational>(),
                                            GS,
                                            GL)};
  }

  if (Quantum::vamdcCheck(qns, Quantum::VAMDC::hundb) and
      qns.contains(QuantumNumberType::N) and
      qns.contains(QuantumNumberType::J) and
      qns.contains(QuantumNumberType::Lambda) and
      qns.contains(QuantumNumberType::S)) {
    auto& N      = qns.at(QuantumNumberType::N);
    auto& J      = qns.at(QuantumNumberType::J);
    auto& Lambda = qns.at(QuantumNumberType::Lambda);
    auto& S      = qns.at(QuantumNumberType::S);
    return {.gu = lbl::zeeman::SimpleGCaseB(N.upper.get<Rational>(),
                                            J.upper.get<Rational>(),
                                            Lambda.upper.get<Rational>(),
                                            S.upper.get<Rational>(),
                                            GS,
                                            GL),
            .gl = lbl::zeeman::SimpleGCaseB(N.lower.get<Rational>(),
                                            J.lower.get<Rational>(),
                                            Lambda.lower.get<Rational>(),
                                            S.lower.get<Rational>(),
                                            GS,
                                            GL)};
  }

  return {};
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
}  // namespace

namespace lbl::zeeman {
data GetSimpleModel(const QuantumIdentifier& qid) {
  const Numeric GS = ::get_lande_spin_constant(qid.isot.spec);
  const Numeric GL = ::get_lande_lambda_constant();
  return ::SimpleG(qid.state, GS, GL);
}

data GetAdvancedModel(const QuantumIdentifier& qid) {
  if (qid.isot == "O2-66"_isot) {
    if (qid.state.contains(QuantumNumberType::J) and
        qid.state.contains(QuantumNumberType::N) and
        qid.state.contains(QuantumNumberType::v)) {
      if (qid.state.at(QuantumNumberType::v).lower.get<Rational>() == 0 and
          qid.state.at(QuantumNumberType::v).upper.get<Rational>() == 0) {
        constexpr Numeric GS  = 2.002084;
        constexpr Numeric GLE = 2.77e-3;
        constexpr Numeric GR  = -1.16e-4;
        constexpr Numeric B   = 43100.44276e6;
        constexpr Numeric D   = 145.1271e3;
        constexpr Numeric H   = 49e-3;
        constexpr Numeric lB  = 59501.3438e6;
        constexpr Numeric lD  = 58.3680e3;
        constexpr Numeric lH  = 290.8e-3;
        constexpr Numeric gB  = -252.58634e6;
        constexpr Numeric gD  = -243.42;
        constexpr Numeric gH  = -1.46e-3;

        const auto& J = qid.state.at(QuantumNumberType::J);
        const auto& N = qid.state.at(QuantumNumberType::N);
        auto JU       = J.upper.get<Rational>();
        auto NU       = N.upper.get<Rational>();
        auto JL       = J.lower.get<Rational>();
        auto NL       = N.lower.get<Rational>();

        Numeric gu = ::case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        Numeric gl = ::case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {.gu = gu, .gl = gl};
      }
    }
  } else if (qid.isot == "O2-68"_isot) {
    if (qid.state.contains(QuantumNumberType::J) and
        qid.state.contains(QuantumNumberType::N) and
        qid.state.contains(QuantumNumberType::v)) {
      if (qid.state.at(QuantumNumberType::v).lower.get<Rational>() == 0 and
          qid.state.at(QuantumNumberType::v).upper.get<Rational>() == 0) {
        constexpr Numeric GS  = 2.002025;
        constexpr Numeric GLE = 2.813e-3;
        constexpr Numeric GR  = -1.26e-4;
        constexpr Numeric B   = 40707.38657e6;
        constexpr Numeric D   = 129.4142e3;
        constexpr Numeric H   = 0;
        constexpr Numeric lB  = 59499.0375e6;
        constexpr Numeric lD  = 54.9777e3;
        constexpr Numeric lH  = 272.1e-3;
        constexpr Numeric gB  = -238.51530e6;
        constexpr Numeric gD  = -217.77;
        constexpr Numeric gH  = -1.305e-3;

        const auto& J = qid.state.at(QuantumNumberType::J);
        const auto& N = qid.state.at(QuantumNumberType::N);
        auto JU       = J.upper.get<Rational>();
        auto NU       = N.upper.get<Rational>();
        auto JL       = J.lower.get<Rational>();
        auto NL       = N.lower.get<Rational>();

        Numeric gu = ::case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        Numeric gl = ::case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {.gu = gu, .gl = gl};
      }
    }
  } else if (qid.isot == "CO-26"_isot) {
    constexpr Numeric gperp =
        -0.2689 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    return {.gu = gperp, .gl = gperp};
  } else if (qid.isot == "OCS-622"_isot) {
    constexpr Numeric gperp =
        -.02889 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        0 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    if (qid.state.contains(QuantumNumberType::J) and
        qid.state.contains(QuantumNumberType::Ka)) {
      const auto& J  = qid.state.at(QuantumNumberType::J);
      const auto& Ka = qid.state.at(QuantumNumberType::K);
      auto JU        = J.upper.get<Rational>();
      auto KU        = Ka.upper.get<Rational>();
      auto JL        = J.lower.get<Rational>();
      auto KL        = Ka.lower.get<Rational>();

      return {.gu = ::closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = ::closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.isot == "OCS-624"_isot) {
    constexpr Numeric gperp =
        -.0285 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        -.061 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    if (qid.state.contains(QuantumNumberType::J) and
        qid.state.contains(QuantumNumberType::Ka)) {
      const auto& J  = qid.state.at(QuantumNumberType::J);
      const auto& Ka = qid.state.at(QuantumNumberType::K);
      auto JU        = J.upper.get<Rational>();
      auto KU        = Ka.upper.get<Rational>();
      auto JL        = J.lower.get<Rational>();
      auto KL        = Ka.lower.get<Rational>();

      return {.gu = ::closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = ::closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.isot == "CO2-626"_isot) {
    constexpr Numeric gperp =
        -.05508 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        0 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    if (qid.state.contains(QuantumNumberType::J) and qid.state.contains( QuantumNumberType::Ka)) {
      const auto& J  = qid.state.at(QuantumNumberType::J);
      const auto& Ka = qid.state.at(QuantumNumberType::K);
      auto JU        = J.upper.get<Rational>();
      auto KU        = Ka.upper.get<Rational>();
      auto JL        = J.lower.get<Rational>();
      auto KL        = Ka.lower.get<Rational>();

      return {.gu = ::closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = ::closed_shell_trilinear(KL, JL, gperp, gpara)};
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

Numeric model::Strength(Rational Ju, Rational Jl, pol type, Index n) const {
  using Math::pow2;

  auto ml = Ml(Ju, Jl, type, n);
  auto mu = Mu(Ju, Jl, type, n);
  auto dm = Rational(dM(type));
  return type == pol::no
             ? 1.0
             : polarization_factor(type) *
                   pow2(wigner3j(Jl, Rational(1), Ju, ml, -dm, -mu));
}

Numeric model::Strength(const QuantumState& qn,
                        pol type,
                        Index n) const {
  if (type == pol::no) return 1.0;

  const auto& J = qn.at(QuantumNumberType::J);

  return Strength(J.upper.get<Rational>(), J.lower.get<Rational>(), type, n);
}

Numeric model::Splitting(const QuantumState& qn,
                         pol type,
                         Index n) const noexcept {
  if (type == pol::no) return 0.0;

  const auto& J = qn.at(QuantumNumberType::J);

  return Splitting(J.upper.get<Rational>(), J.lower.get<Rational>(), type, n);
}

Index model::size(const QuantumState& qn, pol type) const noexcept {
  if (on) {
    if (type == pol::no) return 0;

    const auto& J = qn.at(QuantumNumberType::J);

    return zeeman::size(J.upper.get<Rational>(), J.lower.get<Rational>(), type);
  }

  return static_cast<Index>(type == pol::no);
}

std::istream& operator>>(std::istream& is, model& m) try {
  Index i;
  is >> i >> double_imanip() >> m.mdata.gu >> m.mdata.gl;
  m.on = static_cast<bool>(i);
  return is;
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading Zeeman model data: {}", e.what()));
}

magnetic_angles::magnetic_angles(const Vector3 mag, const Vector2 los)
    : u(mag[0]),
      v(mag[1]),
      w(mag[2]),
      sa(std::sin(Conversion::deg2rad(los[1]))),
      ca(std::cos(Conversion::deg2rad(los[1]))),
      sz(std::sin(Conversion::deg2rad(los[0]))),
      cz(std::cos(Conversion::deg2rad(los[0]))),
      H(std::hypot(u, v, w)),
      uct(ca * sz * v + cz * w + sa * sz * u),
      duct(u * sa * cz + v * ca * cz - w * sz) {}

Numeric magnetic_angles::theta() const {
  return H == 0 ? 0 : std::acos(uct / H);
}

Numeric magnetic_angles::dtheta_du() const {
  using Math::pow2;
  using Math::pow3;
  using std::sqrt;
  const Numeric rat = pow2(uct / H);
  const Numeric nom = u * uct - sa * sz * pow2(H);
  return (H == 0.0 or rat == 1.0) ? 0 : nom / (sqrt(1.0 - rat) * pow3(H));
}

Numeric magnetic_angles::dtheta_dv() const {
  using Math::pow2;
  using Math::pow3;
  using std::sqrt;
  const Numeric rat = pow2(uct / H);
  const Numeric nom = v * uct - ca * sz * pow2(H);
  return (H == 0.0 or rat == 1.0) ? 0 : nom / (sqrt(1.0 - rat) * pow3(H));
}

Numeric magnetic_angles::dtheta_dw() const {
  using Math::pow2;
  using Math::pow3;
  using std::sqrt;
  const Numeric rat = pow2(uct / H);
  const Numeric nom = w * uct - cz * pow2(H);
  return (H == 0.0 or rat == 1.0) ? 0 : nom / (sqrt(1.0 - rat) * pow3(H));
}

Numeric magnetic_angles::eta() const {
  return std::atan2(u * ca - v * sa, duct);
}

Numeric magnetic_angles::deta_du() const {
  return H == 0 ? 0
                : (-ca * sz * w + cz * v) /
                      (Math::pow2(ca * u - sa * v) + Math::pow2(duct));
}

Numeric magnetic_angles::deta_dv() const {
  return H == 0 ? 0
                : (-cz * u + sa * sz * w) /
                      (Math::pow2(ca * u - sa * v) + Math::pow2(duct));
}

Numeric magnetic_angles::deta_dw() const {
  return H == 0 ? 0
                : sz * (ca * u - sa * v) /
                      (Math::pow2(ca * u - sa * v) + Math::pow2(duct));
}

Propmat norm_view(pol p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta = ma.theta();
  const Numeric eta   = ma.eta();

  const Numeric CT  = std::cos(theta);
  const Numeric ST  = std::sin(theta);
  const Numeric C2E = std::cos(2 * eta);
  const Numeric S2E = std::sin(2 * eta);

  switch (p) {
    case pol::pi:
      return {ST * ST,
              -ST * ST * C2E,
              -S2E * ST * ST,
              0,
              0,
              -2 * S2E * ST * ST,
              2 * ST * ST * C2E};
    case pol::sm:
      return {CT * CT + 1,
              ST * ST * C2E,
              S2E * ST * ST,
              2 * CT,
              4 * CT,
              2 * S2E * ST * ST,
              -2 * ST * ST * C2E};
    case pol::sp:
      return {CT * CT + 1,
              ST * ST * C2E,
              S2E * ST * ST,
              -2 * CT,
              -4 * CT,
              2 * S2E * ST * ST,
              -2 * ST * ST * C2E};
    case pol::no: return {1, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}

Propmat dnorm_view_du(pol p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta  = ma.theta();
  const Numeric eta    = ma.eta();
  const Numeric dtheta = ma.dtheta_du();
  const Numeric deta   = ma.deta_du();

  const Numeric CT  = std::cos(theta);
  const Numeric ST  = std::sin(theta);
  const Numeric S2T = 2 * ST * CT;
  const Numeric C2E = std::cos(2 * eta);
  const Numeric S2E = std::sin(2 * eta);

  switch (p) {
    case pol::pi:
      return {S2T * dtheta,
              2 * (S2E * ST * deta - C2E * CT * dtheta) * ST,
              2 * (-S2E * CT * dtheta - ST * C2E * deta) * ST,
              0,
              0,
              4 * (-S2E * CT * dtheta - ST * C2E * deta) * ST,
              4 * (-S2E * ST * deta + C2E * CT * dtheta) * ST};
    case pol::sm:
      return {-S2T * dtheta,
              2 * (-S2E * ST * deta + C2E * CT * dtheta) * ST,
              2 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              -2 * ST * dtheta,
              -4 * ST * dtheta,
              4 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              4 * (S2E * ST * deta - C2E * CT * dtheta) * ST};
    case pol::sp:
      return {-S2T * dtheta,
              2 * (-S2E * ST * deta + C2E * CT * dtheta) * ST,
              2 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              2 * ST * dtheta,
              4 * ST * dtheta,
              4 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              4 * (S2E * ST * deta - C2E * CT * dtheta) * ST};
    case pol::no: return {0, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}

Propmat dnorm_view_dv(pol p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta  = ma.theta();
  const Numeric eta    = ma.eta();
  const Numeric dtheta = ma.dtheta_dv();
  const Numeric deta   = ma.deta_dv();

  const Numeric CT  = std::cos(theta);
  const Numeric ST  = std::sin(theta);
  const Numeric S2T = 2 * ST * CT;
  const Numeric C2E = std::cos(2 * eta);
  const Numeric S2E = std::sin(2 * eta);

  switch (p) {
    case pol::pi:
      return {S2T * dtheta,
              2 * (S2E * ST * deta - C2E * CT * dtheta) * ST,
              2 * (-S2E * CT * dtheta - ST * C2E * deta) * ST,
              0,
              0,
              4 * (-S2E * CT * dtheta - ST * C2E * deta) * ST,
              4 * (-S2E * ST * deta + C2E * CT * dtheta) * ST};
    case pol::sm:
      return {-S2T * dtheta,
              2 * (-S2E * ST * deta + C2E * CT * dtheta) * ST,
              2 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              -2 * ST * dtheta,
              -4 * ST * dtheta,
              4 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              4 * (S2E * ST * deta - C2E * CT * dtheta) * ST};
    case pol::sp:
      return {-S2T * dtheta,
              2 * (-S2E * ST * deta + C2E * CT * dtheta) * ST,
              2 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              2 * ST * dtheta,
              4 * ST * dtheta,
              4 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              4 * (S2E * ST * deta - C2E * CT * dtheta) * ST};
    case pol::no: return {0, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}

Propmat dnorm_view_dw(pol p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta  = ma.theta();
  const Numeric eta    = ma.eta();
  const Numeric dtheta = ma.dtheta_dw();
  const Numeric deta   = ma.deta_dw();

  const Numeric CT  = std::cos(theta);
  const Numeric ST  = std::sin(theta);
  const Numeric S2T = 2 * ST * CT;
  const Numeric C2E = std::cos(2 * eta);
  const Numeric S2E = std::sin(2 * eta);

  switch (p) {
    case pol::pi:
      return {S2T * dtheta,
              2 * (S2E * ST * deta - C2E * CT * dtheta) * ST,
              2 * (-S2E * CT * dtheta - ST * C2E * deta) * ST,
              0,
              0,
              4 * (-S2E * CT * dtheta - ST * C2E * deta) * ST,
              4 * (-S2E * ST * deta + C2E * CT * dtheta) * ST};
    case pol::sm:
      return {-S2T * dtheta,
              2 * (-S2E * ST * deta + C2E * CT * dtheta) * ST,
              2 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              -2 * ST * dtheta,
              -4 * ST * dtheta,
              4 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              4 * (S2E * ST * deta - C2E * CT * dtheta) * ST};
    case pol::sp:
      return {-S2T * dtheta,
              2 * (-S2E * ST * deta + C2E * CT * dtheta) * ST,
              2 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              2 * ST * dtheta,
              4 * ST * dtheta,
              4 * (S2E * CT * dtheta + ST * C2E * deta) * ST,
              4 * (S2E * ST * deta - C2E * CT * dtheta) * ST};
    case pol::no: return {0, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}
}  // namespace lbl::zeeman

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
#include <arts_conversions.h>
#include <double_imanip.h>
#include <wigner_functions.h>

#include <utility>

using enum QuantumNumberType;

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
      qns.contains(Omega) and qns.contains(J) and qns.contains(Lambda) and
      qns.contains(S)) {
    auto& o = qns.at(Omega);
    auto& j = qns.at(J);
    auto& l = qns.at(Lambda);
    auto& s = qns.at(S);
    return {.gu = lbl::zeeman::SimpleGCaseA(
                o.upper, j.upper, l.upper, s.upper, GS, GL),
            .gl = lbl::zeeman::SimpleGCaseA(
                o.lower, j.lower, l.lower, s.lower, GS, GL)};
  }

  if (Quantum::vamdcCheck(qns, Quantum::VAMDC::hundb) and qns.contains(N) and
      qns.contains(J) and qns.contains(Lambda) and qns.contains(S)) {
    auto& n = qns.at(N);
    auto& j = qns.at(J);
    auto& l = qns.at(Lambda);
    auto& s = qns.at(S);
    return {.gu = lbl::zeeman::SimpleGCaseB(
                n.upper, j.upper, l.upper, s.upper, GS, GL),
            .gl = lbl::zeeman::SimpleGCaseB(
                n.lower, j.lower, l.lower, s.lower, GS, GL)};
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
    if (qid.state.contains(J) and qid.state.contains(N) and
        qid.state.contains(v)) {
      if (qid.state.at(v).lower.get<Rational>() == 0 and
          qid.state.at(v).upper.get<Rational>() == 0) {
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

        const auto& j     = qid.state.at(J);
        const auto& n     = qid.state.at(N);
        const Rational JU = j.upper;
        const Rational NU = n.upper;
        const Rational JL = j.lower;
        const Rational NL = n.lower;

        Numeric gu = ::case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        Numeric gl = ::case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {.gu = gu, .gl = gl};
      }
    }
  } else if (qid.isot == "O2-68"_isot) {
    if (qid.state.contains(J) and qid.state.contains(N) and
        qid.state.contains(v)) {
      if (qid.state.at(v).lower.get<Rational>() == 0 and
          qid.state.at(v).upper.get<Rational>() == 0) {
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

        const auto& j     = qid.state.at(J);
        const auto& n     = qid.state.at(N);
        const Rational JU = j.upper;
        const Rational NU = n.upper;
        const Rational JL = j.lower;
        const Rational NL = n.lower;

        Numeric gu = ::case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        Numeric gl = ::case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {.gu = gu, .gl = gl};
      }
    }
  } else if (qid.isot == "CO-26"_isot) {
    // Flygare and Benson 1971
    constexpr Numeric gperp =
        -0.2689 / Constant::mass_ratio_electrons_per_proton;

    return {.gu = gperp, .gl = gperp};
  } else if (qid.isot == "OCS-622"_isot) {
    // Flygare and Benson 1971
    constexpr Numeric gperp =
        -.02889 / Constant::mass_ratio_electrons_per_proton;
    constexpr Numeric gpara = 0 / Constant::mass_ratio_electrons_per_proton;
    if (qid.state.contains(J) and qid.state.contains(Ka)) {
      const auto& j     = qid.state.at(J);
      const auto& Ka    = qid.state.at(K);
      const Rational JU = j.upper;
      const Rational KU = Ka.upper;
      const Rational JL = j.lower;
      const Rational KL = Ka.lower;

      return {.gu = ::closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = ::closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.isot == "OCS-624"_isot) {
    // Flygare and Benson 1971
    constexpr Numeric gperp =
        -.0285 / Constant::mass_ratio_electrons_per_proton;
    constexpr Numeric gpara = -.061 / Constant::mass_ratio_electrons_per_proton;

    if (qid.state.contains(J) and qid.state.contains(Ka)) {
      const auto& j     = qid.state.at(J);
      const auto& Ka    = qid.state.at(K);
      const Rational JU = j.upper;
      const Rational KU = Ka.upper;
      const Rational JL = j.lower;
      const Rational KL = Ka.lower;

      return {.gu = ::closed_shell_trilinear(KU, JU, gperp, gpara),
              .gl = ::closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.isot == "CO2-626"_isot) {
    // Flygare and Benson 1971
    constexpr Numeric gperp =
        -.05508 / Constant::mass_ratio_electrons_per_proton;
    constexpr Numeric gpara = 0 / Constant::mass_ratio_electrons_per_proton;

    if (qid.state.contains(J) and qid.state.contains(Ka)) {
      const auto& j     = qid.state.at(J);
      const auto& Ka    = qid.state.at(K);
      const Rational JU = j.upper;
      const Rational KU = Ka.upper;
      const Rational JL = j.lower;
      const Rational KL = Ka.lower;

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

Numeric model::Strength(Rational Ju,
                        Rational Jl,
                        ZeemanPolarization type,
                        Index n) const {
  using Math::pow2;

  if (type == ZeemanPolarization::no) return 1.0;

  auto ml = Ml(Ju, Jl, type, n);
  auto mu = Mu(Ju, Jl, type, n);

  if (abs(ml) > Jl or abs(mu) > Ju) return 0.0;

  auto dm         = Rational(dM(type));
  const Numeric C = polarization_factor(type);
  return C * pow2(wigner3j(Jl, Rational(1), Ju, ml, dm, -mu));
}

Numeric model::Strength(const QuantumState& qn,
                        ZeemanPolarization type,
                        Index n) const {
  if (type == ZeemanPolarization::no) return 1.0;

  const auto& j = qn.at(J);

  return Strength(j.upper, j.lower, type, n);
}

Numeric model::Splitting(const QuantumState& qn,
                         ZeemanPolarization type,
                         Index n) const noexcept {
  if (type == ZeemanPolarization::no) return 0.0;

  const auto& j = qn.at(J);
  return Splitting(j.upper, j.lower, type, n);
}

Index model::size(const QuantumState& qn,
                  ZeemanPolarization type) const noexcept {
  if (on) {
    if (type == ZeemanPolarization::no) return 0;

    const auto& j = qn.at(J);

    return zeeman::size(j.upper, j.lower, type);
  }

  return static_cast<Index>(type == ZeemanPolarization::no);
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
      uct(sz * sa * u + sz * ca * v + cz * w),
      duct(u * sa * cz + v * ca * cz - w * sz) {
  /**
    * Defines the geometry for the Zeeman effect angles theta and eta.
    * The coordinate system is (x=East, y=North, z=Up).
    * LOS (k) is the line-of-sight vector, B is the magnetic field vector.
    * uct = B · k, which gives theta = acos((B·k)/|B|).
    *
    * k_x = sin(z) * sin(a)
    * k_y = sin(z) * cos(a)
    * k_z = cos(z)
    *
    * B_x = u
    * B_y = v
    * B_z = w
    *
    * For eta, we define a basis on the plane perpendicular to k:
    * e1 ('projected up') = (-cz*sa, -cz*ca, sz)
    * e2 ('right')        = (ca, -sa, 0)
    *
    * The components of the projected B-field on this basis are:
    * B1 = B · e1 = -u*cz*sa - v*cz*ca + w*sz
    * B2 = B · e2 = u*ca - v*sa
    *
    * Note: the code's `duct` variable is equal to -B1.
    * The clockwise angle eta is -atan2(B2, B1), which is implemented as
    * -atan2( (u*ca-v*sa), -duct ).
    */
}

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
  return -std::atan2(ca * u - sa * v, -duct);
}

Numeric magnetic_angles::deta_du() const {
  return H == 0 ? 0
                : (cz * v - ca * sz * w) /
                      (Math::pow2(ca * u - sa * v) + Math::pow2(duct));
}

Numeric magnetic_angles::deta_dv() const {
  return H == 0 ? 0
                : (sa * sz * w - cz * u) /
                      (Math::pow2(ca * u - sa * v) + Math::pow2(duct));
}

Numeric magnetic_angles::deta_dw() const {
  return H == 0 ? 0
                : sz * (ca * u - sa * v) /
                      (Math::pow2(ca * u - sa * v) + Math::pow2(duct));
}

Propmat norm_view(ZeemanPolarization p, Vector3 mag, Vector2 los) {
  /* The propagation matrix elemets for different polarizaitions.

  The 7-element returned vector represents the symmetric matrix
  [[A,B,C,D],[B,A,U,V],[C,-U,A,W],[D,-V,-W,A]], in order [A,B,C,D,U,V,W].
  PI is scaled by 1/2, SM and SP by 1/4 elsewhere.

  PI:
        sin^2(theta)              ; -sin^2(theta) cos(2 eta) ;  sin^2(theta) sin(2 eta) ;  0                       ;
       -sin^2(theta) cos(2 eta)   ;  sin^2(theta)            ;  0                       ;  sin^2(theta) sin(2 eta) ;
        sin^2(theta) sin(2 eta)   ;  0                       ;  sin^2(theta)            ;  sin^2(theta) cos(2 eta) ;
        0                         ; -sin^2(theta) sin(2 eta) ; -sin^2(theta) cos(2 eta) ;  sin^2(theta)            ;

  SM:
        1 + cos^2(theta)          ;  sin^2(theta) cos(2 eta) ; -sin^2(theta) sin(2 eta) ;  2 cos(theta)            ;
        sin^2(theta) cos(2 eta)   ;  1 + cos^2(theta)        ; -2 cos(theta)            ; -sin^2(theta) sin(2 eta) ;
       -sin^2(theta) sin(2 eta)   ;  2 cos(theta)            ;  1 + cos^2(theta)        ; -sin^2(theta) cos(2 eta) ;
        2 cos(theta)              ;  sin^2(theta) sin(2 eta) ;  sin^2(theta) cos(2 eta) ;  1 + cos^2(theta)        ;

  SP:
        1 + cos^2(theta)          ;  sin^2(theta) cos(2 eta) ; -sin^2(theta) sin(2 eta) ; -2 cos(theta)            ;
        sin^2(theta) cos(2 eta)   ;  1 + cos^2(theta)        ;  2 cos(theta)            ; -sin^2(theta) sin(2 eta) ;
       -sin^2(theta) sin(2 eta)   ; -2 cos(theta)            ;  1 + cos^2(theta)        ; -sin^2(theta) cos(2 eta) ;
       -2 cos(theta)              ;  sin^2(theta) sin(2 eta) ;  sin^2(theta) cos(2 eta) ;  1 + cos^2(theta)        ;
  */

  const magnetic_angles ma(mag, los);
  const Numeric theta = ma.theta();
  const Numeric eta   = ma.eta();
  const Numeric CT    = std::cos(theta);
  const Numeric ST2   = Math::pow2(std::sin(theta));
  const Numeric Q     = ST2 * std::cos(2 * eta);
  const Numeric U     = ST2 * std::sin(2 * eta);
  switch (p) {
    using enum ZeemanPolarization;
    case pi: return {ST2, -Q, U, 0, 0, U, Q};
    case sm: return {2 - ST2, Q, -U, 2 * CT, -2 * CT, -U, -Q};
    case sp: return {2 - ST2, Q, -U, -2 * CT, 2 * CT, -U, -Q};
    case no: return {1, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}

Propmat dnorm_view_du(ZeemanPolarization p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta  = ma.theta();
  const Numeric eta    = ma.eta();
  const Numeric dtheta = ma.dtheta_du();
  const Numeric deta   = ma.deta_du();
  const Numeric CT     = std::cos(theta);
  const Numeric ST     = std::sin(theta);
  const Numeric CE     = std::cos(2 * eta);
  const Numeric SE     = std::sin(2 * eta);
  const Numeric ST2    = Math::pow2(ST);
  const Numeric dST2   = 2 * dtheta * ST * CT;
  const Numeric dQ     = 2 * dtheta * ST * CE * CT - 2 * deta * SE * ST2;
  const Numeric dU     = 2 * deta * ST2 * CE + 2 * dtheta * SE * ST * CT;
  const Numeric dCT    = -dtheta * ST;

  switch (p) {
    using enum ZeemanPolarization;
    case pi: return {dST2, -dQ, dU, 0, 0, dU, dQ};
    case sm: return {-dST2, dQ, -dU, 2 * dCT, -2 * dCT, -dU, -dQ};
    case sp: return {-dST2, dQ, -dU, -2 * dCT, 2 * dCT, -dU, -dQ};
    case no: return {0, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}

Propmat dnorm_view_dv(ZeemanPolarization p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta  = ma.theta();
  const Numeric eta    = ma.eta();
  const Numeric dtheta = ma.dtheta_dv();
  const Numeric deta   = ma.deta_dv();
  const Numeric CT     = std::cos(theta);
  const Numeric ST     = std::sin(theta);
  const Numeric CE     = std::cos(2 * eta);
  const Numeric SE     = std::sin(2 * eta);
  const Numeric ST2    = Math::pow2(ST);
  const Numeric dST2   = 2 * dtheta * ST * CT;
  const Numeric dQ     = 2 * dtheta * ST * CE * CT - 2 * deta * SE * ST2;
  const Numeric dU     = 2 * deta * ST2 * CE + 2 * dtheta * SE * ST * CT;
  const Numeric dCT    = -dtheta * ST;

  switch (p) {
    using enum ZeemanPolarization;
    case pi: return {dST2, -dQ, dU, 0, 0, dU, dQ};
    case sm: return {-dST2, dQ, -dU, 2 * dCT, -2 * dCT, -dU, -dQ};
    case sp: return {-dST2, dQ, -dU, -2 * dCT, 2 * dCT, -dU, -dQ};
    case no: return {0, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}

Propmat dnorm_view_dw(ZeemanPolarization p, Vector3 mag, Vector2 los) {
  const magnetic_angles ma(mag, los);
  const Numeric theta  = ma.theta();
  const Numeric eta    = ma.eta();
  const Numeric dtheta = ma.dtheta_dw();
  const Numeric deta   = ma.deta_dw();
  const Numeric CT     = std::cos(theta);
  const Numeric ST     = std::sin(theta);
  const Numeric CE     = std::cos(2 * eta);
  const Numeric SE     = std::sin(2 * eta);
  const Numeric ST2    = Math::pow2(ST);
  const Numeric dST2   = 2 * dtheta * ST * CT;
  const Numeric dQ     = 2 * dtheta * ST * CE * CT - 2 * deta * SE * ST2;
  const Numeric dU     = 2 * deta * ST2 * CE + 2 * dtheta * SE * ST * CT;
  const Numeric dCT    = -dtheta * ST;

  switch (p) {
    using enum ZeemanPolarization;
    case pi: return {dST2, -dQ, dU, 0, 0, dU, dQ};
    case sm: return {-dST2, dQ, -dU, 2 * dCT, -2 * dCT, -dU, -dQ};
    case sp: return {-dST2, dQ, -dU, -2 * dCT, 2 * dCT, -dU, -dQ};
    case no: return {0, 0, 0, 0, 0, 0, 0};
  }

  std::unreachable();
}
}  // namespace lbl::zeeman

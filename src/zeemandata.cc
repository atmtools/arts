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

#include "zeemandata.h"

#include <cstdio>
#include <stdexcept>

#include "arts_constexpr_math.h"
#include "arts_conversions.h"
#include "debug.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "species_info.h"
#include "wigner_functions.h"

Zeeman::SplittingData SimpleG(const Quantum::Number::ValueList& qns,
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
    return {.gu = Zeeman::SimpleGCaseA(
                Omega.upp(), J.upp(), Lambda.upp(), S.upp(), GS, GL),
            .gl = Zeeman::SimpleGCaseA(
                Omega.low(), J.low(), Lambda.low(), S.low(), GS, GL)};
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
    return {.gu = Zeeman::SimpleGCaseB(
                N.upp(), J.upp(), Lambda.upp(), S.upp(), GS, GL),
            .gl = Zeeman::SimpleGCaseB(
                N.low(), J.low(), Lambda.low(), S.low(), GS, GL)};
  }

  return {NAN, NAN};
}

Zeeman::Model Zeeman::GetSimpleModel(const QuantumIdentifier& qid)
    ARTS_NOEXCEPT {
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

Zeeman::Model Zeeman::GetAdvancedModel(const QuantumIdentifier& qid)
    ARTS_NOEXCEPT {
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

        auto JU = qid.val[QuantumNumberType::J].upp();
        auto NU = qid.val[QuantumNumberType::N].upp();
        Numeric gu = case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        auto JL = qid.val[QuantumNumberType::J].low();
        auto NL = qid.val[QuantumNumberType::N].low();
        Numeric gl = case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {gu, gl};
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

        auto JU = qid.val[QuantumNumberType::J].upp();
        auto NU = qid.val[QuantumNumberType::N].upp();
        Numeric gu = case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        auto JL = qid.val[QuantumNumberType::J].low();
        auto NL = qid.val[QuantumNumberType::N].low();
        Numeric gl = case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return {gu, gl};
      }
    }
  } else if (qid.Isotopologue() == "CO-26") {
    constexpr Numeric gperp =
        -0.2689 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    return {gperp, gperp};
  } else if (qid.Isotopologue() == "OCS-622") {
    constexpr Numeric gperp =
        -.02889 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        0 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    if (qid.val.has(QuantumNumberType::J, QuantumNumberType::Ka)) {
      auto JU = qid.val[QuantumNumberType::J].upp();
      auto KU = qid.val[QuantumNumberType::Ka].upp();
      auto JL = qid.val[QuantumNumberType::J].low();
      auto KL = qid.val[QuantumNumberType::Ka].low();

      return {closed_shell_trilinear(KU, JU, gperp, gpara),
              closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.Isotopologue() == "OCS-624") {
    constexpr Numeric gperp =
        -.0285 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        -.061 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    if (qid.val.has(QuantumNumberType::J, QuantumNumberType::Ka)) {
      auto JU = qid.val[QuantumNumberType::J].upp();
      auto KU = qid.val[QuantumNumberType::Ka].upp();
      auto JL = qid.val[QuantumNumberType::J].low();
      auto KL = qid.val[QuantumNumberType::Ka].low();

      return {closed_shell_trilinear(KU, JU, gperp, gpara),
              closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  } else if (qid.Isotopologue() == "CO2-626") {
    constexpr Numeric gperp =
        -.05508 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    constexpr Numeric gpara =
        0 /
        Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971

    if (qid.val.has(QuantumNumberType::J, QuantumNumberType::Ka)) {
      auto JU = qid.val[QuantumNumberType::J].upp();
      auto KU = qid.val[QuantumNumberType::Ka].upp();
      auto JL = qid.val[QuantumNumberType::J].low();
      auto KL = qid.val[QuantumNumberType::Ka].low();

      return {closed_shell_trilinear(KU, JU, gperp, gpara),
              closed_shell_trilinear(KL, JL, gperp, gpara)};
    }
  }

  // Set to zero otherwise as we practically say "there's no Zeeman effect" then
  return {0, 0};
}

Zeeman::Model::Model(const QuantumIdentifier& qid) noexcept {
  Model m = GetAdvancedModel(qid);
  if (m.empty()) m = GetSimpleModel(qid);
  *this = m;
}

Zeeman::Derived Zeeman::FromGrids(
    Numeric u, Numeric v, Numeric w, Numeric z, Numeric a) noexcept {
  using Math::pow2, Math::pow3, Conversion::sind, Conversion::cosd;
  using std::acos, std::hypot, std::sqrt, std::atan2;

  Derived output;
  const Numeric sa = sind(a);
  const Numeric ca = cosd(a);
  const Numeric sz = sind(z);
  const Numeric cz = cosd(z);
  const Numeric H = hypot(u, v, w);
  const Numeric uct = ca * sz * v + cz * w + sa * sz * u;
  const Numeric duct = u * sa * cz + v * ca * cz - w * sz;

  output.H = H;
  if (H == 0) {
    output.dH_du = 0;
    output.dH_dv = 0;
    output.dH_dw = 0;
  } else {
    const Numeric inv = 1.0 / H;
    output.dH_du = u * inv;
    output.dH_dv = v * inv;
    output.dH_dw = w * inv;
  }

  output.theta = H == 0 ? 0 : acos(uct / H);
  const Numeric rat = pow2(uct / H);
  if (H == 0.0 or rat == 1.0) {
    output.dtheta_du = 0;
    output.dtheta_dv = 0;
    output.dtheta_dw = 0;
  } else {
    const Numeric inv = 1.0 / (sqrt(1.0 - rat) * pow3(H));
    const Numeric nomu = u * uct - sa * sz * pow2(H);
    const Numeric nomv = v * uct - ca * sz * pow2(H);
    const Numeric nomw = w * uct - cz * pow2(H);
    output.dtheta_du = nomu * inv;
    output.dtheta_dv = nomv * inv;
    output.dtheta_dw = nomw * inv;
  }

  output.eta = -atan2(ca * u - sa * v, -duct);
  if (H == 0.0) {
    output.deta_du = 0;
    output.deta_dv = 0;
    output.deta_dw = 0;
  } else {
    const Numeric inv = 1.0 / (pow2(ca * u - sa * v) + pow2(duct));
    output.deta_du = (cz * v - ca * sz * w) * inv;
    output.deta_dv = (sa * sz * w - cz * u) * inv;
    output.deta_dw = sz * (ca * u - sa * v) * inv;
  }

  return output;
}

namespace Zeeman {
Numeric Model::Strength(Rational Ju,
                        Rational Jl,
                        Zeeman::Polarization type,
                        Index n) const ARTS_NOEXCEPT {
  if (type == Zeeman::Polarization::None) return 1.0;

  using Math::pow2;

  auto ml = Ml(Ju, Jl, type, n);
  auto mu = Mu(Ju, Jl, type, n);

  if (abs(ml) > Jl or abs(mu) > Ju) return 0.0;

  auto dm = Rational(dM(type));
  const Numeric C = PolarizationFactor(type);
  return C * pow2(wigner3j(Jl, Rational(1), Ju, ml, dm, -mu));
}

std::ostream& operator<<(std::ostream& os, const Model& m) {
  os << m.mdata.gu << ' ' << m.mdata.gl;
  return os;
}

std::istream& operator>>(std::istream& is, Model& m) {
  is >> double_imanip() >> m.mdata.gu >> m.mdata.gl;
  return is;
}

std::ostream& operator<<(bofstream& bof, const Model& m) {
  bof << m.mdata.gu << m.mdata.gl;
  return bof;
}

std::istream& operator>>(bifstream& bif, Model& m) {
  bif >> m.mdata.gu >> m.mdata.gl;
  return bif;
}

AllPolarizationVectors AllPolarization(Numeric theta, Numeric eta) noexcept {
  const Numeric CT = std::cos(theta);
  const Numeric ST2 = Math::pow2(std::sin(theta));
  const Numeric BW = ST2 * std::cos(2 * eta);
  const Numeric CV = ST2 * std::sin(2 * eta);

  return {.sm = {2 - ST2, -BW, -CV, -2 * CT, -2 * CT, CV, -BW},
          .pi = {ST2, BW, CV, 0, 0, -CV, BW},
          .sp = {2 - ST2, -BW, -CV, 2 * CT, 2 * CT, CV, -BW}};
}

AllPolarizationVectors AllPolarization_dtheta(Numeric theta,
                                              const Numeric eta) noexcept {
  const Numeric ST = std::sin(theta);
  const Numeric dST2 = 2 * std::cos(theta) * ST;
  const Numeric dBW = std::cos(2 * eta) * dST2;
  const Numeric dCV = std::sin(2 * eta) * dST2;
  const Numeric dCT = -ST;

  return {.sm = {-dST2, -dBW, -dCV, -2 * dCT, -2 * dCT, dCV, -dBW},
          .pi = {dST2, dBW, dCV, 0, 0, -dCV, dBW},
          .sp = {-dST2, -dBW, -dCV, 2 * dCT, 2 * dCT, dCV, -dBW}};
}

AllPolarizationVectors AllPolarization_deta(Numeric theta,
                                            Numeric eta) noexcept {
  const Numeric ST2 = Math::pow2(std::sin(theta));
  const Numeric dBW = -2 * std::sin(2 * eta) * ST2;
  const Numeric dCV = 2 * std::cos(2 * eta) * ST2;

  return {.sm = {0, -dBW, -dCV, 0, 0, dCV, -dBW},
          .pi = {0, dBW, dCV, 0, 0, -dCV, dBW},
          .sp = {0, -dBW, -dCV, 0, 0, dCV, -dBW}};
}

const PolarizationVector& SelectPolarization(const AllPolarizationVectors& data,
                                             Polarization type) {
  switch (type) {
    case Polarization::SigmaMinus:
      return data.sm;
    case Polarization::Pi:
      return data.pi;
    case Polarization::SigmaPlus:
      return data.sp;
    case Polarization::None:;
  }
  throw std::logic_error("Unknown polarization type");
}

void sum(PropagationMatrix& pm,
         const ComplexVectorView& abs,
         const PolarizationVector& polvec,
         const bool do_phase) ARTS_NOEXCEPT {
  ARTS_ASSERT(pm.NumberOfZenithAngles() == 1)
  ARTS_ASSERT(pm.NumberOfAzimuthAngles() == 1)
  ARTS_ASSERT(pm.NumberOfFrequencies() == abs.nelem())
  ARTS_ASSERT(do_phase ? pm.NumberOfNeededVectors() == 7
                       : pm.NumberOfNeededVectors() == 4)

  const ExhaustiveConstVectorView pol_real(polvec.att);
  const ExhaustiveConstVectorView pol_imag(polvec.dis);

  MatrixView out = pm.Data()(0, 0, joker, joker);
  matpack::eigen::mat(out).leftCols<4>().noalias() +=
      matpack::eigen::row_vec(abs.real()) * matpack::eigen::col_vec(pol_real);
  if (do_phase)
    matpack::eigen::mat(out).rightCols<3>().noalias() +=
        matpack::eigen::row_vec(abs.imag()) * matpack::eigen::col_vec(pol_imag);
}

void dsum(PropagationMatrix& pm,
          const ComplexVectorView& abs,
          const ComplexVectorView& dabs,
          const PolarizationVector& polvec,
          const PolarizationVector& dpolvec_dtheta,
          const PolarizationVector& dpolvec_deta,
          const Numeric dH,
          const Numeric dt,
          const Numeric de,
          const bool do_phase) ARTS_NOEXCEPT {
  ARTS_ASSERT(pm.NumberOfZenithAngles() == 1)
  ARTS_ASSERT(pm.NumberOfAzimuthAngles() == 1)
  ARTS_ASSERT(pm.NumberOfFrequencies() == abs.nelem())
  ARTS_ASSERT(do_phase ? pm.NumberOfNeededVectors() == 7
                       : pm.NumberOfNeededVectors() == 4)

  const ExhaustiveConstVectorView pol_real(polvec.att);
  const ExhaustiveConstVectorView pol_imag(polvec.dis);
  const ExhaustiveConstVectorView dpolvec_dtheta_real(dpolvec_dtheta.att);
  const ExhaustiveConstVectorView dpolvec_dtheta_imag(dpolvec_dtheta.dis);
  const ExhaustiveConstVectorView dpolvec_deta_real(dpolvec_deta.att);
  const ExhaustiveConstVectorView dpolvec_deta_imag(dpolvec_deta.dis);

  auto da_r = (dt * matpack::eigen::col_vec(dpolvec_dtheta_real) +
               de * matpack::eigen::col_vec(dpolvec_deta_real));
  auto da_i = (dt * matpack::eigen::col_vec(dpolvec_dtheta_imag) +
               de * matpack::eigen::col_vec(dpolvec_deta_imag));

  MatrixView out = pm.Data()(0, 0, joker, joker);
  matpack::eigen::mat(out).leftCols<4>().noalias() +=
      dH * matpack::eigen::row_vec(dabs.real()) *
          matpack::eigen::col_vec(pol_real) +
      matpack::eigen::row_vec(abs.real()) * da_r;
  if (do_phase)
    matpack::eigen::mat(out).rightCols<3>().noalias() +=
        dH * matpack::eigen::row_vec(dabs.imag()) *
            matpack::eigen::col_vec(pol_imag) +
        matpack::eigen::row_vec(abs.imag()) * da_i;
}
}  // namespace Zeeman

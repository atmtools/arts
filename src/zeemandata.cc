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
    return {Zeeman::SimpleGCaseA(
                Omega.upp(), J.upp(), Lambda.upp(), S.upp(), GS, GL),
            Zeeman::SimpleGCaseA(
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
    return {
        Zeeman::SimpleGCaseB(N.upp(), J.upp(), Lambda.upp(), S.upp(), GS, GL),
        Zeeman::SimpleGCaseB(N.low(), J.low(), Lambda.low(), S.low(), GS, GL)};
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

Eigen::Vector3d los_xyz_by_uvw_local(Numeric u, Numeric v, Numeric w) {
  return Eigen::Vector3d(v, u, w).normalized();
}

Eigen::Vector3d los_xyz_by_za_local(Numeric z, Numeric a) {
  using std::cos;
  using std::sin;
  return Eigen::Vector3d(cos(a) * sin(z), sin(a) * sin(z), cos(z)).normalized();
}

Eigen::Vector3d ev_xyz_by_za_local(Numeric z, Numeric a) {
  using std::cos;
  using std::sin;
  return Eigen::Vector3d(cos(a) * cos(z), sin(a) * cos(z), -sin(z)).normalized();
}

Zeeman::Derived Zeeman::FromGrids(
    Numeric u, Numeric v, Numeric w, Numeric z, Numeric a) noexcept {
  Derived output;

  // If there is no magnetic field, bailout quickly
  output.H = std::hypot(u, v, w);
  if (output.H == 0) {
    output = FromPreDerived(0, 0, 0);
  } else {
    // XYZ vectors normalized
    const Eigen::Vector3d n = los_xyz_by_za_local(z, a);
    const Eigen::Vector3d ev = ev_xyz_by_za_local(z, a);
    const Eigen::Vector3d nH = los_xyz_by_uvw_local(u, v, w);

    // Normalized vector (which are also the magnetic field derivatives)
    output.dH_dv = nH[0];
    output.dH_du = nH[1];
    output.dH_dw = nH[2];

    // Compute theta (and its derivatives if possible)
    const Numeric cos_theta = n.dot(nH);
    const Numeric sin_theta = std::sqrt(1 - Math::pow2(cos_theta));
    output.theta = std::acos(cos_theta);
    if (sin_theta not_eq 0) {
      const Eigen::Vector3d dtheta =
          (nH * cos_theta - n) / (output.H * sin_theta);
      output.dtheta_dv = dtheta[0];
      output.dtheta_du = dtheta[1];
      output.dtheta_dw = dtheta[2];
    } else {
      output.dtheta_dv = 0;
      output.dtheta_du = 0;
      output.dtheta_dw = 0;
    }

    // Compute eta (and its derivatives if possible)
    const Eigen::Vector3d inplane = nH - nH.dot(n) * n;
    const Numeric y = ev.cross(inplane).dot(n);
    const Numeric x = ev.dot(inplane);
    output.eta = std::atan2(y, x);
    if (x not_eq 0 or y not_eq 0) {
      const Eigen::Vector3d deta =
          n.cross(nH) / (output.H * (Math::pow2(x) + Math::pow2(y)));
      output.deta_dv = deta[0];
      output.deta_du = deta[1];
      output.deta_dw = deta[2];
    } else {
      output.deta_dv = 0;
      output.deta_du = 0;
      output.deta_dw = 0;
    }
  }

  return output;
}

namespace Zeeman {
Numeric Model::Strength(Rational Ju,
                        Rational Jl,
                        Zeeman::Polarization type,
                        Index n) const ARTS_NOEXCEPT {
  ARTS_ASSERT(type not_eq Zeeman::Polarization::None);
  using Math::pow2;

  auto ml = Ml(Ju, Jl, type, n);
  auto mu = Mu(Ju, Jl, type, n);
  auto dm = Rational(dM(type));
  return PolarizationFactor(type) *
         pow2(wigner3j(Jl, Rational(1), Ju, ml, -dm, -mu));
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
  const Numeric ST = std::sin(theta), CT = std::cos(theta), ST2 = ST * ST,
                CT2 = CT * CT, ST2C2E = ST2 * std::cos(2 * eta),
                ST2S2E = ST2 * std::sin(2 * eta);

  AllPolarizationVectors pv;
  pv.sm = PolarizationVector(
      1 + CT2, ST2C2E, ST2S2E, 2 * CT, 4 * CT, 2 * ST2S2E, -2 * ST2C2E);
  pv.pi =
      PolarizationVector(ST2, -ST2C2E, -ST2S2E, 0, 0, -2 * ST2S2E, 2 * ST2C2E);
  pv.sp = PolarizationVector(
      1 + CT2, ST2C2E, ST2S2E, -2 * CT, -4 * CT, 2 * ST2S2E, -2 * ST2C2E);
  return pv;
}

AllPolarizationVectors AllPolarization_dtheta(Numeric theta,
                                              const Numeric eta) noexcept {
  const Numeric ST = std::sin(theta), CT = std::cos(theta),
                C2E = std::cos(2 * eta), S2E = std::sin(2 * eta), dST = CT,
                dST2 = 2 * ST * dST, dCT = -ST, dST2C2E = dST2 * C2E,
                dST2S2E = dST2 * S2E, dCT2 = 2 * CT * dCT;

  AllPolarizationVectors pv;
  pv.sm = PolarizationVector(
      dCT2, dST2C2E, dST2S2E, 2 * dCT, 4 * dCT, 2 * dST2S2E, -2 * dST2C2E);
  pv.pi = PolarizationVector(
      dST2, -dST2C2E, -dST2S2E, 0, 0, -2 * dST2S2E, 2 * dST2C2E);
  pv.sp = PolarizationVector(
      dCT2, dST2C2E, dST2S2E, -2 * dCT, -4 * dCT, 2 * dST2S2E, -2 * dST2C2E);
  return pv;
}

AllPolarizationVectors AllPolarization_deta(Numeric theta,
                                            Numeric eta) noexcept {
  const Numeric ST = std::sin(theta), ST2 = ST * ST, C2E = std::cos(2 * eta),
                S2E = std::sin(2 * eta), dST2C2E = -2 * ST2 * S2E,
                dST2S2E = 2 * ST2 * C2E;

  AllPolarizationVectors pv;
  pv.sm =
      PolarizationVector(0, dST2C2E, dST2S2E, 0, 0, 2 * dST2S2E, -2 * dST2C2E);
  pv.pi = PolarizationVector(
      0, -dST2C2E, -dST2S2E, 0, 0, -2 * dST2S2E, 2 * dST2C2E);
  pv.sp =
      PolarizationVector(0, dST2C2E, dST2S2E, 0, 0, 2 * dST2S2E, -2 * dST2C2E);
  return pv;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
const PolarizationVector& SelectPolarization(const AllPolarizationVectors& data,
                                             Polarization type) noexcept {
  switch (type) {
    case Polarization::SigmaMinus:
      return data.sm;
    case Polarization::Pi:
      return data.pi;
    case Polarization::SigmaPlus:
      return data.sp;
    case Polarization::None:
      return data.sm;  //! This should never be reached
  }
}
#pragma GCC diagnostic pop

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
  matpack::eigen::mat(out).leftCols<4>().noalias() += matpack::eigen::row_vec(abs.real()) * matpack::eigen::col_vec(pol_real);
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

  auto da_r =
      (dt * matpack::eigen::col_vec(dpolvec_dtheta_real) + de * matpack::eigen::col_vec(dpolvec_deta_real));
  auto da_i =
      (dt * matpack::eigen::col_vec(dpolvec_dtheta_imag) + de * matpack::eigen::col_vec(dpolvec_deta_imag));

  MatrixView out = pm.Data()(0, 0, joker, joker);
  matpack::eigen::mat(out).leftCols<4>().noalias() +=
      dH * matpack::eigen::row_vec(dabs.real()) * matpack::eigen::col_vec(pol_real) + matpack::eigen::row_vec(abs.real()) * da_r;
  if (do_phase)
    matpack::eigen::mat(out).rightCols<3>().noalias() +=
        dH * matpack::eigen::row_vec(dabs.imag()) * matpack::eigen::col_vec(pol_imag) + matpack::eigen::row_vec(abs.imag()) * da_i;
}
}  // namespace Zeeman

/* Copyright (C) 2018 Richard Larsson

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
   USA. */

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
#include "abs_species_tags.h"
#include "species_info.h"

Zeeman::Model Zeeman::GetSimpleModel(const QuantumIdentifier& qid) noexcept {
  const Numeric GS = get_lande_spin_constant(qid.Species());
  const Numeric GL = get_lande_lambda_constant();
  const Numeric gu = SimpleG(qid.UpperQuantumNumbers(), GS, GL);
  const Numeric gl = SimpleG(qid.LowerQuantumNumbers(), GS, GL);
  return Model(gu, gl);
}

Numeric case_b_g_coefficient_o2(Rational j,
                                Rational n,
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
  using Constant::pow2;
  using Constant::pow3;
  using std::atan2;
  using std::cos;
  using std::sin;
  using std::sqrt;

  if (j.isUndefined() or n.isUndefined())
    return NAN;
  else if (j == 0)
    return 0;

  auto J = j.toNumeric();

  auto nom = (lB + lD * (J * J + J + 1) + lH * pow2(J * J + J + 1)) * (2 * sqrt(J * J + J) / (2 * J + 1));

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

  if (j == n)
    return (GS + GR) / (J * (J + 1)) - GR;
  else if (j < n)
    return (GS + GR) * (pow2(cos(phi)) / J - pow2(sin(phi)) / (J + 1)) +
           2 * GLE * cos(2 * phi) / (2 * J + 1) - GR;
  else /*if(j > n)*/
    return (GS + GR) * (pow2(sin(phi)) / J - pow2(cos(phi)) / (J + 1)) -
           2 * GLE * cos(2 * phi) / (2 * J + 1) - GR;
}

constexpr Numeric closed_shell_trilinear(Rational k,
                                         Rational j,
                                         Numeric gperp,
                                         Numeric gpara)
{
  if (k.isUndefined() or j.isUndefined() or j == 0)
    return 0;
  else
    return gperp + (gperp + gpara) * Numeric(k*k / (j*(j+1)));
}

Zeeman::Model Zeeman::GetAdvancedModel(const QuantumIdentifier& qid) noexcept {
  const auto name = qid.SpeciesName();
  if (name == "O2-66") {
    if (qid.LowerQuantumNumber(QuantumNumberType::v1) == 0 and
        qid.UpperQuantumNumber(QuantumNumberType::v1) == 0) {
        Numeric GS = 2.002084;
        Numeric GLE = 2.77e-3;
        Numeric GR = -1.16e-4;
        Numeric B = 43100.44276e6;
        Numeric D = 145.1271e3;
        Numeric H = 49e-3;
        Numeric lB = 59501.3438e6;
        Numeric lD = 58.3680e3;
        Numeric lH = 290.8e-3;
        Numeric gB = -252.58634e6;
        Numeric gD = -243.42;
        Numeric gH = -1.46e-3;

        auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
        auto NU = qid.UpperQuantumNumber(QuantumNumberType::N);
        Numeric gu = case_b_g_coefficient_o2(
            JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
        auto NL = qid.LowerQuantumNumber(QuantumNumberType::N);
        Numeric gl = case_b_g_coefficient_o2(
            JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
        return Model(gu, gl);
      }
    } else if (name == "O2-68") {
    if (qid.LowerQuantumNumber(QuantumNumberType::v1) == 0 and
        qid.UpperQuantumNumber(QuantumNumberType::v1) == 0) {
      Numeric GS = 2.002025;
      Numeric GLE = 2.813e-3;
      Numeric GR = -1.26e-4;
      Numeric B = 40707.38657e6;
      Numeric D = 129.4142e3;
      Numeric H = 0;
      Numeric lB = 59499.0375e6;
      Numeric lD = 54.9777e3;
      Numeric lH = 272.1e-3;
      Numeric gB = -238.51530e6;
      Numeric gD = -217.77;
      Numeric gH = -1.305e-3;

      auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
      auto NU = qid.UpperQuantumNumber(QuantumNumberType::N);
      Numeric gu = case_b_g_coefficient_o2(
          JU, NU, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
      auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
      auto NL = qid.LowerQuantumNumber(QuantumNumberType::N);
      Numeric gl = case_b_g_coefficient_o2(
          JL, NL, GS, GR, GLE, B, D, H, gB, gD, gH, lB, lD, lH);
      return Model(gu, gl);
    }
  } else if (name == "CO-26") {
      Numeric gperp = -0.2689 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      
      return Model(gperp, gperp);
  } else if (name == "OCS-622") {
      Numeric gperp = -.02889 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      Numeric gpara = 0 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      
      auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
      auto KU = qid.UpperQuantumNumber(QuantumNumberType::Ka);
      auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
      auto KL = qid.LowerQuantumNumber(QuantumNumberType::Ka);
      
      return Model(closed_shell_trilinear(KU, JU, gperp, gpara),
                   closed_shell_trilinear(KL, JL, gperp, gpara));
  } else if (name == "OCS-624") {
    Numeric gperp = -.0285 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    Numeric gpara = -.061 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    
    auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
    auto KU = qid.UpperQuantumNumber(QuantumNumberType::Ka);
    auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
    auto KL = qid.LowerQuantumNumber(QuantumNumberType::Ka);
    
    return Model(closed_shell_trilinear(KU, JU, gperp, gpara),
                 closed_shell_trilinear(KL, JL, gperp, gpara));
  } else if (name == "CO2-626") {
    Numeric gperp = -.05508 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    Numeric gpara = 0 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
    
    auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
    auto KU = qid.UpperQuantumNumber(QuantumNumberType::Ka);
    auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
    auto KL = qid.LowerQuantumNumber(QuantumNumberType::Ka);
    
    return Model(closed_shell_trilinear(KU, JU, gperp, gpara),
                 closed_shell_trilinear(KL, JL, gperp, gpara));
  }
  
  // Take care of zeroes since they do not show up in replacement databases
  const bool upperzero = qid.UpperQuantumNumber(QuantumNumberType::J) == 0 or qid.UpperQuantumNumber(QuantumNumberType::F) == 0;
  const bool lowerzero = qid.LowerQuantumNumber(QuantumNumberType::J) == 0 or qid.LowerQuantumNumber(QuantumNumberType::F) == 0;
  return Model(upperzero ? 0 : NAN, lowerzero ? 0 : NAN);
}

Zeeman::Model::Model(const QuantumIdentifier& qid) noexcept {
  Model m = GetAdvancedModel(qid);
  if (m.empty()) m = GetSimpleModel(qid);
  *this = m;
}

Eigen::Vector3d los_xyz_by_uvw_local(Numeric u, Numeric v, Numeric w)
{
  return Eigen::Vector3d(v, u, w).normalized();
}

Eigen::Vector3d los_xyz_by_za_local(Numeric z, Numeric a)
{
  using std::cos;
  using std::sin;
  return Eigen::Vector3d(cos(a)*sin(z), sin(a)*sin(z), cos(z));
}

Eigen::Vector3d ev_xyz_by_za_local(Numeric z, Numeric a)
{
  using std::cos;
  using std::sin;
  return Eigen::Vector3d(cos(a)*cos(z), sin(a)*cos(z), -sin(z));
}

Zeeman::Derived Zeeman::FromGrids(Numeric u, Numeric v, Numeric w, Numeric z, Numeric a) noexcept
{
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
    const Numeric sin_theta = std::sqrt(1 - Constant::pow2(cos_theta));
    output.theta = std::acos(cos_theta);
    if (sin_theta not_eq 0) {
      const Eigen::Vector3d dtheta = (nH * cos_theta - n) / (output.H * sin_theta);
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
      const Eigen::Vector3d deta = n.cross(nH) / (output.H * (Constant::pow2(x) + Constant::pow2(y)));
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
Numeric Model::Strength(Rational Ju, Rational Jl, Zeeman::Polarization type, Index n) const {
  using Constant::pow2;
  
  auto ml = Ml(Ju, Jl, type, n);
  auto mu = Mu(Ju, Jl, type, n);
  auto dm = Rational(dM(type));
  return PolarizationFactor(type) * pow2(wigner3j(Jl, Rational(1), Ju, ml, -dm, -mu));
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

AllPolarizationVectors AllPolarization(Numeric theta,
                                       Numeric eta) noexcept {
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

AllPolarizationVectors AllPolarization_dtheta(
    Numeric theta, const Numeric eta) noexcept {
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
const PolarizationVector& SelectPolarization(
    const AllPolarizationVectors& data, Polarization type) noexcept {
  switch (type) {
    case Polarization::SigmaMinus:
      return data.sm;
    case Polarization::Pi:
      return data.pi;
    case Polarization::SigmaPlus:
      return data.sp;
  }
}
#pragma GCC diagnostic pop

void sum(PropagationMatrix& pm, const ComplexVector& abs, const PolarizationVector& polvec) {
  auto pol_real = polvec.attenuation();
  auto pol_imag = polvec.dispersion();
  
  for (Index iv=0; iv<pm.NumberOfFrequencies(); iv++) {
    pm.Kjj()[iv] += abs[iv].real() * pol_real[0];
    pm.K12()[iv] += abs[iv].real() * pol_real[1];
    pm.K13()[iv] += abs[iv].real() * pol_real[2];
    pm.K14()[iv] += abs[iv].real() * pol_real[3];
    pm.K23()[iv] += abs[iv].imag() * pol_imag[0];
    pm.K24()[iv] += abs[iv].imag() * pol_imag[1];
    pm.K34()[iv] += abs[iv].imag() * pol_imag[2];
  }
}

void dsum(PropagationMatrix& pm,
          const ComplexVector& abs,
          const ComplexVector& dabs,
          const PolarizationVector& polvec,
          const PolarizationVector& dpolvec_dtheta,
          const PolarizationVector& dpolvec_deta,
          const Numeric dH,
          const Numeric dt,
          const Numeric de) {
  auto pol_real = polvec.attenuation();
  auto pol_imag = polvec.dispersion();
  auto dp_dt_r = dpolvec_dtheta.attenuation();
  auto dp_dt_i = dpolvec_dtheta.dispersion();
  auto dp_de_r = dpolvec_deta.attenuation();
  auto dp_de_i = dpolvec_deta.dispersion();
  
  for (Index iv=0; iv<pm.NumberOfFrequencies(); iv++) {
    const auto da_r = (dt * dp_dt_r + de * dp_de_r);
    const auto da_i = (dt * dp_dt_i + de * dp_de_i);
    pm.Kjj()[iv] += dabs[iv].real() * dH * pol_real[0] + abs[iv].real() * da_r[0];
    pm.K12()[iv] += dabs[iv].real() * dH * pol_real[1] + abs[iv].real() * da_r[1];
    pm.K13()[iv] += dabs[iv].real() * dH * pol_real[2] + abs[iv].real() * da_r[2];
    pm.K14()[iv] += dabs[iv].real() * dH * pol_real[3] + abs[iv].real() * da_r[3];
    pm.K23()[iv] += dabs[iv].imag() * dH * pol_imag[0] + abs[iv].imag() * da_i[0];
    pm.K24()[iv] += dabs[iv].imag() * dH * pol_imag[1] + abs[iv].imag() * da_i[1];
    pm.K34()[iv] += dabs[iv].imag() * dH * pol_imag[2] + abs[iv].imag() * da_i[2];;
  }
}
}  // Zeeman

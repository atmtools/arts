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
  return Model({gu, gl});
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
  if (qid.SpeciesName() == "O2") {
    if (qid.Isotopologue() == SpeciesTag("O2-66").Isotopologue()) {
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
        return Model({gu, gl});
      }
    } else if (qid.Isotopologue() == SpeciesTag("O2-68").Isotopologue()) {
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
        return Model({gu, gl});
      }
    }
  } else if (qid.SpeciesName() == "CO") {
    if (qid.Isotopologue() == SpeciesTag("CO-26").Isotopologue()) {
      Numeric gperp = -0.2689 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      
      return Model({gperp, gperp});
    }
  } else if (qid.SpeciesName() == "OCS") {
    if (qid.Isotopologue() == SpeciesTag("OCS-622").Isotopologue()) {
      Numeric gperp = -.02889 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      Numeric gpara = 0 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      
      auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
      auto KU = qid.UpperQuantumNumber(QuantumNumberType::Ka);
      auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
      auto KL = qid.LowerQuantumNumber(QuantumNumberType::Ka);
      
      return Model({closed_shell_trilinear(KU, JU, gperp, gpara),
                    closed_shell_trilinear(KL, JL, gperp, gpara)});
    } else if (qid.Isotopologue() == SpeciesTag("OCS-624").Isotopologue()) {
      Numeric gperp = -.0285 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      Numeric gpara = -.061 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      
      auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
      auto KU = qid.UpperQuantumNumber(QuantumNumberType::Ka);
      auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
      auto KL = qid.LowerQuantumNumber(QuantumNumberType::Ka);
      
      return Model({closed_shell_trilinear(KU, JU, gperp, gpara),
                    closed_shell_trilinear(KL, JL, gperp, gpara)});
    }
  } else if (qid.SpeciesName() == "CO2") {
    if (qid.Isotopologue() == SpeciesTag("CO2-626").Isotopologue()) {
      Numeric gperp = -.05508 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      Numeric gpara = 0 / Constant::mass_ratio_electrons_per_proton;  // Flygare and Benson 1971
      
      auto JU = qid.UpperQuantumNumber(QuantumNumberType::J);
      auto KU = qid.UpperQuantumNumber(QuantumNumberType::Ka);
      auto JL = qid.LowerQuantumNumber(QuantumNumberType::J);
      auto KL = qid.LowerQuantumNumber(QuantumNumberType::Ka);
      
      return Model({closed_shell_trilinear(KU, JU, gperp, gpara),
                    closed_shell_trilinear(KL, JL, gperp, gpara)});
    }
  }
  
  // Take care of zeroes since they do not show up in replacement databases
  const bool upperzero = qid.UpperQuantumNumber(QuantumNumberType::J) == 0 or qid.UpperQuantumNumber(QuantumNumberType::F) == 0;
  const bool lowerzero = qid.LowerQuantumNumber(QuantumNumberType::J) == 0 or qid.LowerQuantumNumber(QuantumNumberType::F) == 0;
  return Model({upperzero ? 0 : NAN, lowerzero ? 0 : NAN});
}

Zeeman::Model::Model(const QuantumIdentifier& qid) noexcept {
  Model m = GetAdvancedModel(qid);
  if (m.empty()) m = GetSimpleModel(qid);
  *this = m;
}

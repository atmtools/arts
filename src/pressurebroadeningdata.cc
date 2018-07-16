/* Copyright (C) 2014 
 R *ichard Larsson <ric.larsson@gmail.com>
 
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

/** Contains additional functions for the pressure broadening data class
 * \file   pressurebroadeningdata.cc
 * 
 * \author Richard Larsson
 * \date   2014-11-06
 **/
#include "pressurebroadeningdata.h"
#include "linerecord.h"


inline Numeric test_pressure_shift(const Numeric& T, const Numeric& T0, const Numeric& A, const Numeric& d, const Numeric& m)
{
  return d * pow(T0/T, m) * (1 + A * log(T/T0));
}


inline Numeric test_pressure_broadening(const Numeric& T, const Numeric& T0, const Numeric& g, const Numeric& n)
{
  return g * pow(T0/T, n);
}


inline Numeric test_pressure_broadening_speed_term(const Numeric& T, const Numeric& T0, const Numeric& g, const Numeric& n)
{
  return g + n * (T - T0);
}

inline void voigt_test_params(Numeric& G0, Numeric& D0, 
                              const Numeric& aP, const Numeric& sP,  const Numeric& T, const Numeric& T0, 
                              const Numeric& sA, const Numeric& sg0, const Numeric& sn, 
                              const Numeric& sd0, const Numeric& sm,
                              const Numeric& aA, const Numeric& ag0, const Numeric& an, 
                              const Numeric& ad0, const Numeric& am)
{
  G0 = sP * test_pressure_broadening(T, T0, sg0, sn) + aP * test_pressure_broadening(T, T0, ag0, an);
  D0 = sP * test_pressure_shift(T, T0, sA, sd0, sm) + aP *test_pressure_shift(T, T0, aA, ad0, am);
}

inline void speed_dependent_test_params(Numeric& G0, Numeric& D0, Numeric& G2, Numeric& D2, 
                                        const Numeric& aP, const Numeric& sP, const Numeric& T, const Numeric& T0,
                                        const Numeric& sA, const Numeric& sg0, const Numeric& sn0,
                                        const Numeric& sg2, const Numeric& sn2, const Numeric& sd0,
                                        const Numeric& sm, const Numeric& sd2,
                                        const Numeric& aA, const Numeric& ag0, const Numeric& an0,
                                        const Numeric& ag2, const Numeric& an2, const Numeric& ad0,
                                        const Numeric& am, const Numeric& ad2)
{
  G0 = sP * test_pressure_broadening(T, T0, sg0, sn0) + aP * test_pressure_broadening(T, T0, ag0, an0);
  G2 = sP * test_pressure_broadening_speed_term(T, T0, sg2, sn2) + aP * test_pressure_broadening_speed_term(T, T0, ag2, an2);
  D0 = sP * test_pressure_shift(T, T0, sA, sd0, sm)   + aP * test_pressure_shift(T, T0, aA, ad0, am);
  D2 = sP * sd2 + aP * ad2;  // test_pressure_shift(T, T0, 0, d2, 0);
}

///////////////////////////////////////////
//  Broadening calculations below here
//////////////////////////////////////////

// Get broadening parameters
void PressureBroadeningData::GetPressureBroadeningParams(Numeric& gamma_0,
                                                         Numeric& gamma_2,
                                                         Numeric& eta,
                                                         Numeric& df_0,
                                                         Numeric& df_2,
                                                         Numeric& f_VC,
                                                         const Numeric& temperature,
                                                         const Numeric& ref_temperature,
                                                         const Numeric& pressure,
                                                         const Numeric& self_pressure,
                                                         const Index    this_species,
                                                         const Index    h2o_species,
                                                         const ArrayOfIndex& broad_spec_locations,
                                                         ConstVectorView vmrs) const
{
  gamma_0=gamma_2=eta=df_0=df_2=f_VC=0;
  switch(mtype) {
    case PB_NONE:
      // Note that this is oftentimes not wanted, but a valid case at low pressures
      break;
    case PB_AIR_BROADENING:
      GetAirBroadening(gamma_0, df_0, ref_temperature / temperature, pressure, self_pressure);
      break;
    case PB_AIR_AND_WATER_BROADENING:
      GetAirAndWaterBroadening(gamma_0, df_0, ref_temperature / temperature, pressure, self_pressure, 
                               this_species, h2o_species, vmrs);
      break;
    case PB_PLANETARY_BROADENING:
      GetPlanetaryBroadening(gamma_0, df_0, ref_temperature / temperature, pressure, self_pressure, broad_spec_locations, vmrs);
      break;
    case PB_SD_AIR_VOLUME:
      GetSDAirBroadening(gamma_0,gamma_2,df_0,df_2,ref_temperature / temperature,pressure);
      break;
    case PB_PURELY_FOR_TESTING:
      GetTestBroadening(gamma_0, gamma_2, df_0, vmrs, ref_temperature / temperature, pressure, h2o_species);
      break;
    case PB_HTP_AIR_VOLUME:
      throw std::runtime_error("Not implemented");
      break;
    case PB_VOIGT_TEST_WATER:
      voigt_test_params(gamma_0, df_0, 
                        pressure-self_pressure, self_pressure, temperature, ref_temperature, 
                        mdata[0][Index(TestParams::sA)], mdata[0][Index(TestParams::sg0)], 
                        mdata[0][Index(TestParams::sn0)], mdata[0][Index(TestParams::sd0)],
                        mdata[0][Index(TestParams::sm)],
                        mdata[0][Index(TestParams::aA)], mdata[0][Index(TestParams::ag0)], 
                        mdata[0][Index(TestParams::an0)], mdata[0][Index(TestParams::ad0)],
                        mdata[0][Index(TestParams::am)]);
      break;
    case PB_SD_TEST_WATER:
      speed_dependent_test_params(gamma_0, df_0, gamma_2, df_2, 
                                  pressure-self_pressure, self_pressure, temperature, ref_temperature,
                                  mdata[0][Index(TestParams::sA)], mdata[0][Index(TestParams::sg0)], 
                                  mdata[0][Index(TestParams::sn0)], mdata[0][Index(TestParams::sg2)], 
                                  mdata[0][Index(TestParams::sn2)], mdata[0][Index(TestParams::sd0)],
                                  mdata[0][Index(TestParams::sm)], mdata[0][Index(TestParams::sd2)],
                                  mdata[0][Index(TestParams::aA)], mdata[0][Index(TestParams::ag0)], 
                                  mdata[0][Index(TestParams::an0)], mdata[0][Index(TestParams::ag2)], 
                                  mdata[0][Index(TestParams::an2)], mdata[0][Index(TestParams::ad0)],
                                  mdata[0][Index(TestParams::am)], mdata[0][Index(TestParams::ad2)]);
      break;
  }
}

// Get temperature derivative of the broadening
void PressureBroadeningData::GetPressureBroadeningParams_dT(Numeric& dgamma_0_dT,
                                                            Numeric& dgamma_2_dT,
                                                            Numeric& deta_dT,
                                                            Numeric& ddf_0_dT,
                                                            Numeric& ddf_2_dT,
                                                            Numeric& df_VC_dT,
                                                            const Numeric& T,
                                                            const Numeric& T0,
                                                            const Numeric& pressure,
                                                            const Numeric& self_pressure,
                                                            const Index    this_species,
                                                            const Index    h2o_species,
                                                            const ArrayOfIndex& broad_spec_locations,
                                                            ConstVectorView vmrs) const
{
  dgamma_0_dT=dgamma_2_dT=deta_dT=ddf_0_dT=ddf_2_dT=df_VC_dT=0;
  switch(mtype)
  {
    case PB_NONE:
      // Note that this is oftentimes not wanted, but a valid case at low pressures
      break;
    case PB_AIR_BROADENING:
      GetAirBroadening_dT(dgamma_0_dT, ddf_0_dT, T, T0, pressure, self_pressure);
      break;
    case PB_AIR_AND_WATER_BROADENING:
      GetAirAndWaterBroadening_dT(dgamma_0_dT, ddf_0_dT, T, T0, pressure, self_pressure, this_species, h2o_species, vmrs);
      break;
    case PB_PLANETARY_BROADENING:
      GetPlanetaryBroadening_dT(dgamma_0_dT, ddf_0_dT, T, T0, pressure, self_pressure, broad_spec_locations, vmrs);
      break;
    case PB_SD_AIR_VOLUME:
      GetSDAirBroadening_dT(dgamma_0_dT, dgamma_2_dT, ddf_0_dT, ddf_2_dT, T, T0, pressure);
      break;
    case PB_HTP_AIR_VOLUME:
      GetHTPAirBroadening_dT(dgamma_0_dT, dgamma_2_dT, ddf_0_dT, ddf_2_dT, df_VC_dT, deta_dT, T, T0, pressure);
      break;
    default:
      throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
  }
}


void PressureBroadeningData::SetInternalDerivatives(ComplexVector& derivatives, 
                                                    const ArrayOfRetrievalQuantity& ppd, 
                                                    const QuantumIdentifier& QI,
                                                    const Numeric& theta,
                                                    const Numeric& pressure,
                                                    const Numeric& self_pressure,
                                                    const Index    this_species,
                                                    const Index    h2o_species,
                                                    ConstVectorView vmrs) const
{
  const Index nppd = ppd.nelem();
  
  ComplexVector res(11);
  Numeric results1, results2;
  Index ipd = 0;
  
  // nb that continue is here to not count wrongly the number of parameters
  for(Index iq = 0; iq < nppd; iq++)
  {
    if(ppd[iq] == JacPropMatType::VMR)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dSelfVMR(results1, results2, theta, pressure);
        res[ipd] = Complex(-results2, results1);
      }
      else 
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineGammaSelf)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dSelfGamma(results1, theta, self_pressure);
        res[ipd] = Complex(0.0, results1);
      }
      else 
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineGammaForeign)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dForeignGamma(results1, theta, pressure, self_pressure, 
                                                  this_species, h2o_species, vmrs);
        res[ipd] = Complex(0.0, results1);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineGammaWater)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dWaterGamma(results1, theta, pressure, 
                                                this_species, h2o_species, vmrs);
        res[ipd] = Complex(0.0, results1);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineShiftSelf)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dSelfPsf(results2, theta, self_pressure);
        res[ipd] = Complex(-results2, 0.0);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineShiftForeign)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dForeignPsf(results2, theta, pressure, self_pressure, 
                                                this_species, h2o_species, vmrs);
        res[ipd] = Complex(-results2, 0);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineShiftWater)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dWaterPsf(results2, theta, pressure, 
                                              this_species, h2o_species, vmrs);
        res[ipd] = Complex(-results2, 0);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineGammaSelfExp)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dSelfExponent(results1, results2, theta, self_pressure);
        res[ipd] = Complex(-results2, results1);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineGammaForeignExp)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dForeignExponent(results1, results2, theta, pressure, self_pressure, 
                                                     this_species, h2o_species, vmrs);
        res[ipd] = Complex(-results2, results1);
      }
      else
        continue;
    }
    else if(ppd[iq] == JacPropMatType::LineGammaWaterExp)
    {
      if(QI > ppd[iq].QuantumIdentity())
      {
        GetPressureBroadeningParams_dWaterExponent(results1, results2, theta, pressure, 
                                                   this_species, h2o_species, vmrs);
        res[ipd] = Complex(-results2, results1);
      }
      else
        continue;
    }
    else
      continue;
    
    // Only activate this when something hit the target
    ++ipd;
  }
  
  derivatives.resize(ipd);
  for(Index iq = 0; iq < ipd; iq++)
    derivatives[iq] = res[iq];
}



// Get catalog parameter derivatives:  self broadening gamma
void PressureBroadeningData::GetPressureBroadeningParams_dSelfGamma(Numeric& gamma_dSelf,
                                                                    const Numeric& theta,
                                                                    const Numeric& self_pressure) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            GetAirBroadening_dSelfGamma(gamma_dSelf, theta, self_pressure);
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dSelfGamma(gamma_dSelf,theta, self_pressure);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  foreign broadening gamma
void PressureBroadeningData::GetPressureBroadeningParams_dForeignGamma(Numeric& gamma_dForeign,
                                                                       const Numeric& theta,
                                                                       const Numeric& pressure,
                                                                       const Numeric& self_pressure,
                                                                       const Index    this_species,
                                                                       const Index    h2o_species,
                                                                       ConstVectorView vmrs) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            GetAirBroadening_dForeignGamma(gamma_dForeign, theta, pressure, self_pressure);
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dForeignGamma(gamma_dForeign, theta, pressure, self_pressure, 
                                                   this_species, h2o_species, vmrs);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  water broadening gamma
void PressureBroadeningData::GetPressureBroadeningParams_dWaterGamma(Numeric& gamma_dWater,
                                                                     const Numeric& theta,
                                                                     const Numeric& pressure,
                                                                     const Index    this_species,
                                                                     const Index    h2o_species,
                                                                     ConstVectorView vmrs) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            throw std::runtime_error("Air broadening calculation type does not support water broadening partial derivatives.\n"
            "Please check your catalog type and input lines to ensure that you are doing what you expect.\n");
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dWaterGamma(gamma_dWater, theta, pressure, 
                                                 this_species, h2o_species, vmrs);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  self broadening gamma
void PressureBroadeningData::GetPressureBroadeningParams_dSelfPsf(Numeric& psf_dSelf,
                                                                  const Numeric& theta,
                                                                  const Numeric& self_pressure) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            throw std::runtime_error("Air broadening calculations does not support self pressure shift partial derivatives.\n"
            "Please check your catalog type and input lines to ensure that you are doing what you expect.\n");
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dSelfPsf(psf_dSelf,theta, self_pressure);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  foreign broadening gamma
void PressureBroadeningData::GetPressureBroadeningParams_dForeignPsf(Numeric& psf_dForeign,
                                                                     const Numeric& theta,
                                                                     const Numeric& pressure,
                                                                     const Numeric& self_pressure,
                                                                     const Index    this_species,
                                                                     const Index    h2o_species,
                                                                     ConstVectorView vmrs) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            GetAirBroadening_dForeignPsf(psf_dForeign, theta, pressure);
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dForeignPsf(psf_dForeign, theta, pressure, self_pressure, 
                                                 this_species, h2o_species, vmrs);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  water broadening psf
void PressureBroadeningData::GetPressureBroadeningParams_dWaterPsf(Numeric& psf_dWater,
                                                                   const Numeric& theta,
                                                                   const Numeric& pressure,
                                                                   const Index    this_species,
                                                                   const Index    h2o_species,
                                                                   ConstVectorView vmrs) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            throw std::runtime_error("Air broadening calculation type does not support water broadening partial derivatives.\n"
            "Please check your catalog type and input lines to ensure that you are doing what you expect.\n");
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dWaterPsf(psf_dWater, theta, pressure, 
                                               this_species, h2o_species, vmrs);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  self broadening exponent
void PressureBroadeningData::GetPressureBroadeningParams_dSelfExponent(Numeric& gamma_dSelfExponent,
                                                                       Numeric& psf_dSelfExponent,
                                                                       const Numeric& theta,
                                                                       const Numeric& self_pressure) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            GetAirBroadening_dSelfExponent(gamma_dSelfExponent, psf_dSelfExponent, theta, self_pressure);
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dSelfExponent(gamma_dSelfExponent, psf_dSelfExponent, 
                                                   theta, self_pressure);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get catalog parameter derivatives:  foreign broadening exponent
void PressureBroadeningData::GetPressureBroadeningParams_dForeignExponent(Numeric& gamma_dForeignExponent,
                                                                          Numeric& psf_dForeignExponent,
                                                                          const Numeric& theta,
                                                                          const Numeric& pressure,
                                                                          const Numeric& self_pressure,
                                                                          const Index    this_species,
                                                                          const Index    h2o_species,
                                                                          ConstVectorView vmrs) const
{
    
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            GetAirBroadening_dForeignExponent(gamma_dForeignExponent, psf_dForeignExponent, 
                                              theta, pressure, self_pressure);
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dForeignExponent(gamma_dForeignExponent, psf_dForeignExponent,
                                                      theta, pressure, self_pressure, 
                                                      this_species, h2o_species, vmrs);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

void PressureBroadeningData::GetPressureBroadeningParams_dWaterExponent(Numeric& gamma_dWaterExponent,
                                                                        Numeric& psf_dWaterExponent,
                                                                        const Numeric& theta,
                                                                        const Numeric& pressure,
                                                                        const Index    this_species,
                                                                        const Index    h2o_species,
                                                                        ConstVectorView vmrs) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            throw std::runtime_error("Air broadening calculation type does not support water broadening partial derivatives.\n"
                "Please check your catalog type and input lines to ensure that you are doing what you expect.\n");
            break;
        case PB_AIR_AND_WATER_BROADENING:
            GetAirAndWaterBroadening_dWaterExponent(gamma_dWaterExponent, psf_dWaterExponent,
                                                    theta, pressure, 
                                                    this_species, h2o_species, vmrs);
            break;
        case PB_PLANETARY_BROADENING:
            throw std::runtime_error("Planetary broadening calculation type do not support "
                                     "pressure broadening partial derivatives.\n");
            break;
        default:
            throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

// Get VMR influence
void PressureBroadeningData::GetPressureBroadeningParams_dSelfVMR(Numeric& gamma_dvmr,
                                                                  Numeric& split_dvmr,
                                                                  const Numeric& theta,
                                                                  const Numeric& pressure) const
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
          gamma_dvmr = pressure * ( -mdata[2][0] * pow(theta,mdata[3][0]) + mdata[0][0] * pow(theta,mdata[1][0]));
          split_dvmr = 0.0;
          break;
        case PB_AIR_AND_WATER_BROADENING:
          gamma_dvmr = pressure * (- mdata[1][0] * pow(theta,mdata[1][1]) + mdata[0][0] * pow(theta,mdata[0][1]));
          split_dvmr = pressure * (- mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]));
          break;
        case PB_PLANETARY_BROADENING:
          gamma_dvmr = mdata[0][0] * pow(theta, mdata[1][0]) * pressure;
          split_dvmr = 0.0;
          break;
        default:
          throw std::runtime_error("You have defined an unknown broadening mechanism.\n");
    }
}

///////////////////////////////////////////////////////////////////////////////////////
// Get air broadening
///////////////////////////////////////////////////////////////////////////////////////

// This is the broadening used by ARTSCAT-3; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening(Numeric& gamma,
                                              Numeric& deltaf,
                                              const Numeric& theta,
                                              const Numeric& pressure,
                                              const Numeric& self_pressure) const
{
    gamma   =
      mdata[2][0] * pow(theta,mdata[3][0]) * (pressure-self_pressure)
      + mdata[0][0] * pow(theta,mdata[1][0]) *         self_pressure;
    deltaf  =
      mdata[4][0] * pressure
      * pow (theta,(Numeric)0.25+(Numeric)1.5*mdata[3][0]);
}

// This is the temperature derivative of the broadening used by ARTSCAT-3 type of broadening; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening_dT(Numeric& dgamma_dT,
                                                 Numeric& ddeltaf_dT,
                                                 const Numeric& T,
                                                 const Numeric& T0,
                                                 const Numeric& pressure,
                                                 const Numeric& self_pressure) const
{
    const Numeric theta = T0/T;
    
    dgamma_dT   = - (mdata[3][0] *
    mdata[2][0] * pow(theta,mdata[3][0]) * (pressure-self_pressure)
    + mdata[1][0] * 
    mdata[0][0] * pow(theta,mdata[1][0]) *           self_pressure) / T;
    
    ddeltaf_dT  = - ((Numeric)0.25+(Numeric)1.5*mdata[3][0])/T *
    mdata[4][0] * pressure
    * pow (theta,(Numeric)0.25+(Numeric)1.5*mdata[3][0]);
}

// This is the self broadening derivative of the broadening used by ARTSCAT-3; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening_dSelfGamma(Numeric& gamma_dSelf,
                                                         const Numeric& theta,
                                                         const Numeric& self_pressure) const
{
    gamma_dSelf =  pow(theta, mdata[1][0]) * self_pressure;
}

// This is the foreign broadening derivative of the broadening used by ARTSCAT-3; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening_dForeignGamma(Numeric& gamma_dForeign,
                                                            const Numeric& theta,
                                                            const Numeric& pressure,
                                                            const Numeric& self_pressure) const
{
    gamma_dForeign = pow(theta, mdata[3][0]) * (pressure-self_pressure);
}

// This is the foreign broadening derivative of the broadening used by ARTSCAT-3; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening_dForeignPsf(Numeric& psf_dForeign,
                                                          const Numeric& theta,
                                                          const Numeric& pressure) const
{
    psf_dForeign  = pressure * pow (theta,(Numeric)0.25+(Numeric)1.5*mdata[3][0]);
}

// This is the self broadening exponent derivative of the broadening used by ARTSCAT-3; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening_dSelfExponent(Numeric& gamma_dSelfExponent,
                                                            Numeric& psf_dSelfExponent,
                                                            const Numeric& theta,
                                                            const Numeric& self_pressure) const
{
    Numeric log_theta = log(theta);
    gamma_dSelfExponent = mdata[0][0] * pow(theta,mdata[1][0]) * self_pressure * log_theta;
    psf_dSelfExponent = 0.; //Note that we assume the broadening is from 'air' so no self regardless of species
}

// This is the foreign broadening exponent derivative of the broadening used by ARTSCAT-3; the "N2"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirBroadening_dForeignExponent(Numeric& gamma_dForeignExponent,
                                                               Numeric& psf_dForeignExponent,
                                                               const Numeric& theta,
                                                               const Numeric& pressure,
                                                               const Numeric& self_pressure) const
{
    Numeric log_theta = log(theta);
    gamma_dForeignExponent = mdata[2][0] * pow(theta,mdata[3][0]) * (pressure-self_pressure) * log_theta;
    psf_dForeignExponent  = 
    mdata[4][0] * pressure * pow (theta,(Numeric)0.25+(Numeric)1.5*mdata[3][0]) * log_theta * 1.5;
}

///////////////////////////////////////////////////////////////////////////////////////
// Get air and water broadening
///////////////////////////////////////////////////////////////////////////////////////

// This is the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening(Numeric& gamma,
                                                      Numeric& deltaf,
                                                      const Numeric& theta,
                                                      const Numeric& pressure,
                                                      const Numeric& self_pressure,
                                                      const Index    this_species,
                                                      const Index    h2o_species,
                                                      ConstVectorView vmrs) const
{
    if(this_species==h2o_species)    
    {
        gamma   =
        mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][0] * pow(theta,mdata[0][1]) *         self_pressure;
        
        deltaf  =
        mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) *         self_pressure;
    }
    else if(h2o_species==-1)
    {
        gamma   =
        mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][0] * pow(theta,mdata[0][1]) *         self_pressure;
        
        deltaf  =
        mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) *         self_pressure;
    }
    else
    {
        gamma   =
        mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure)
        + mdata[0][0] * pow(theta,mdata[0][1]) * self_pressure
        + mdata[2][0] * pow(theta,mdata[2][1]) * vmrs[h2o_species]*pressure;
        
        deltaf  =
        mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure)
        + mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) * self_pressure
        + mdata[2][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[2][1]) * vmrs[h2o_species]*pressure;
    }
}

// This is the temperature derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dT(Numeric& dgamma_dT,
                                                         Numeric& ddeltaf_dT,
                                                         const Numeric& T,
                                                         const Numeric& T0,
                                                         const Numeric& pressure,
                                                         const Numeric& self_pressure,
                                                         const Index    this_species,
                                                         const Index    h2o_species,
                                                         ConstVectorView vmrs) const
{
    const Numeric theta = T0/T;
    
    if(this_species==h2o_species)    
    {
        dgamma_dT   = - (
          mdata[1][1] * mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][1] * mdata[0][0] * pow(theta,mdata[0][1]) *           self_pressure) / T;
        
        ddeltaf_dT  = - (
          ((Numeric)0.25+(Numeric)1.5*mdata[1][1]) * mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure)
        + ((Numeric)0.25+(Numeric)1.5*mdata[0][1]) * mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) *           self_pressure) /T;
    }
    else if(h2o_species==-1)
    {
        dgamma_dT   = - (
          mdata[1][1] * mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure)
        + mdata[0][1] * mdata[0][0] * pow(theta,mdata[0][1]) *           self_pressure) /T;
        
        ddeltaf_dT  = - (
        ((Numeric)0.25+(Numeric)1.5*mdata[1][1]) * mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure)
        + ((Numeric)0.25+(Numeric)1.5*mdata[0][1]) * mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) *         self_pressure) / T;
    }
    else
    {
        dgamma_dT   = - (
        mdata[1][1] * mdata[1][0] * pow(theta,mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure)
        + mdata[0][1] * mdata[0][0] * pow(theta,mdata[0][1]) * self_pressure
        + mdata[2][1] * mdata[2][0] * pow(theta,mdata[2][1]) * vmrs[h2o_species]*pressure ) / T;
        
        ddeltaf_dT  = - (
        ((Numeric)0.25+(Numeric)1.5*mdata[1][1]) * mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure)
        + ((Numeric)0.25+(Numeric)1.5*mdata[0][1]) * mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) * self_pressure
        + ((Numeric)0.25+(Numeric)1.5*mdata[2][1]) * mdata[2][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[2][1]) * vmrs[h2o_species]*pressure) / T;
    }
}

// This is the self broadening derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dSelfGamma(Numeric& gamma_dSelf,
                                                                 const Numeric& theta,
                                                                 const Numeric& self_pressure) const
{
    gamma_dSelf   =  pow(theta,mdata[0][1]) * self_pressure;
}

// This is the foreign broadening derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dForeignGamma(Numeric& gamma_dForeign,
                                                                    const Numeric& theta,
                                                                    const Numeric& pressure,
                                                                    const Numeric& self_pressure,
                                                                    const Index    this_species,
                                                                    const Index    h2o_species,
                                                                    ConstVectorView vmrs) const
{
    if(this_species==h2o_species || h2o_species==-1)    
        gamma_dForeign = pow(theta,mdata[1][1]) * (pressure-self_pressure);
    else
        gamma_dForeign = pow(theta,mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure);
}

// This is the water broadening derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dWaterGamma(Numeric& gamma_dWater,
                                                                  const Numeric& theta,
                                                                  const Numeric& pressure,
                                                                  const Index    this_species,
                                                                  const Index    h2o_species,
                                                                  ConstVectorView vmrs) const
{
    if(this_species==h2o_species)
        throw std::runtime_error("Use \"Self broadening\" types of derivatives rather than water broadening for water lines.\n");
    else if(h2o_species==-1)
    {
        gamma_dWater = 0.0;
        
//         out2 << "You have no H2O in species but you want the water broadening derivative.  It is thus set to zero.\n";
    }
    else
        gamma_dWater = pow(theta,mdata[2][1]) * vmrs[h2o_species]*pressure;
}

// This is the self broadening derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dSelfPsf(Numeric& psf_dSelf,
                                                               const Numeric& theta,
                                                               const Numeric& self_pressure) const
{
    psf_dSelf = pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) * self_pressure;
}

// This is the foreign broadening derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dForeignPsf(Numeric& psf_dForeign,
                                                                  const Numeric& theta,
                                                                  const Numeric& pressure,
                                                                  const Numeric& self_pressure,
                                                                  const Index    this_species,
                                                                  const Index    h2o_species,
                                                                  ConstVectorView vmrs) const
{
    if(this_species==h2o_species || h2o_species==-1)    
        psf_dForeign = pow(theta,mdata[1][1]) * (pressure-self_pressure);
    else
        psf_dForeign = pow(theta,mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure);
}

// This is the water broadening derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dWaterPsf(Numeric& psf_dWater,
                                                                const Numeric& theta,
                                                                const Numeric& pressure,
                                                                const Index    this_species,
                                                                const Index    h2o_species,
                                                                ConstVectorView vmrs) const
{
    if(this_species==h2o_species)
        throw std::runtime_error("Use \"Self broadening\" types of derivatives rather than water broadening for water lines.\n");
    else if(h2o_species==-1)
    {
        psf_dWater = 0.0;
        
//         out2 << "You have no H2O in species but you want the water broadening derivative.  It is thus set to zero.\n";
    }
    else
        psf_dWater = pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[2][1]) * vmrs[h2o_species]*pressure;
}

// This is the self broadening exponent derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dSelfExponent(Numeric& gamma_dSelfExponent,
                                                                    Numeric& psf_dSelfExponent,
                                                                    const Numeric& theta,
                                                                    const Numeric& self_pressure) const
{
    Numeric log_theta = log(theta);
    gamma_dSelfExponent =  mdata[0][0] * pow(theta,mdata[0][1]) * self_pressure * log_theta;
    psf_dSelfExponent = mdata[0][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[0][1]) * self_pressure * 1.5 * log_theta;
}

// This is the foreign broadening exponent derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dForeignExponent(Numeric& gamma_dForeignExponent,
                                                                       Numeric& psf_dForeignExponent,
                                                                       const Numeric& theta,
                                                                       const Numeric& pressure,
                                                                       const Numeric& self_pressure,
                                                                       const Index    this_species,
                                                                       const Index    h2o_species,
                                                                       ConstVectorView vmrs) const
{
    Numeric log_theta = log(theta);
    if(this_species==h2o_species || h2o_species==-1) 
    {
        gamma_dForeignExponent = mdata[1][0]* pow(theta,mdata[1][1]) * (pressure-self_pressure) * log_theta;
        psf_dForeignExponent = mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure) * log_theta * 1.5;
    }
    else
    {
        gamma_dForeignExponent = mdata[1][0]* pow(theta,mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure) * log_theta;
        psf_dForeignExponent = mdata[1][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[1][1]) * (pressure-self_pressure-vmrs[h2o_species]*pressure) * log_theta * 1.5;
    }
}

// This is the water broadening exponent derivative of the broadening used by the "WA"-tag in ARTSCAT-5
void PressureBroadeningData::GetAirAndWaterBroadening_dWaterExponent(Numeric& gamma_dWaterExponent,
                                                                     Numeric& psf_dWaterExponent,
                                                                     const Numeric& theta,
                                                                     const Numeric& pressure,
                                                                     const Index    this_species,
                                                                     const Index    h2o_species,
                                                                     ConstVectorView vmrs) const
{
    if(this_species==h2o_species)
        throw std::runtime_error("Use \"Self broadening\" types of derivatives rather than water broadening for water lines.\n");
    else if(h2o_species==-1)
    {
        gamma_dWaterExponent = 0.0;
        psf_dWaterExponent = 0.;
        
//         out2 << "You have no H2O in species but you want the water broadening derivative.  It is thus set to zero.\n";
    }
    else
    {
        Numeric log_theta = log(theta);
        gamma_dWaterExponent = mdata[2][0] * pow(theta,mdata[2][1]) * vmrs[h2o_species]*pressure * log_theta;
        psf_dWaterExponent = mdata[2][2] * pow(theta,(Numeric)0.25+(Numeric)1.5*mdata[2][1]) * vmrs[h2o_species]*pressure * 1.5 * log_theta;
    }
}


///////////////////////////////////////////////////////////////////////////////////////
// Get all planetary broadening
///////////////////////////////////////////////////////////////////////////////////////

// This is the broadening used by ARTSCAT-4; the "AP"-tag in ARTSCAT-5
void PressureBroadeningData::GetPlanetaryBroadening(Numeric& gamma,
                                                    Numeric& deltaf,
                                                    const Numeric& theta,
                                                    const Numeric& pressure,
                                                    const Numeric& self_pressure,
                                                    const ArrayOfIndex& broad_spec_locations,
                                                    ConstVectorView vmrs) const
{
    // Number of broadening species:
    const Index nbs = LineRecord::NBroadSpec();
    assert(nbs==broad_spec_locations.nelem());
    
    // Split total pressure in self and foreign part:
    const Numeric foreign_pressure = pressure - self_pressure;
    
    // Calculate sum of VMRs of all available foreign broadening species (we need this
    // for normalization). The species "Self" will not be included in the sum!
    Numeric broad_spec_vmr_sum = 0;
    
    // Gamma is the line width. We first initialize gamma with the self width
    gamma = mdata[0][0] * pow(theta, mdata[1][0]) * self_pressure;
    
    // Set deltaf to 0
    deltaf = 0;
    
    // and treat foreign width separately:
    Numeric foreign_gamma = 0;
    
    // Add up foreign broadening species, where available:
    for (Index i=0; i<nbs; ++i) {
        if ( broad_spec_locations[i] < -1 ) {
            // -2 means that this broadening species is identical to Self.
            // Throw runtime errors if the parameters are not identical.
            if (mdata[2][i]!=mdata[0][0] ||
                mdata[3][i]!=mdata[1][0])
            {
                std::ostringstream os;
                os << "Inconsistency in LineRecord, self broadening and line "
                << "broadening for " << LineRecord::BroadSpecName(i) << "\n"
                << "should be identical.\n";
                throw std::runtime_error(os.str());
            }
        } else if ( broad_spec_locations[i] >= 0 ) {
            
            // Add to VMR sum:
            broad_spec_vmr_sum += vmrs[broad_spec_locations[i]];
            
            // foreign broadening:
            foreign_gamma +=  mdata[2][i] * pow(theta, mdata[3][i])
                * vmrs[broad_spec_locations[i]];
            
            // Delta f (not .25+1.5*foreign_broadening)
            deltaf += mdata[4][i]
            * pow( theta , (Numeric).25 + (Numeric)1.5*mdata[3][i] )
            * vmrs[broad_spec_locations[i]];
        }
    }
    
//     // Check that sum of self and all foreign VMRs is not too far from 1:
//     if ( abs(vmrs[this_species]+broad_spec_vmr_sum-1) > 0.1
//         && out2.sufficient_priority() )
//     {
//         std::ostringstream os;
//         os << "Warning: The total VMR of all your defined broadening\n"
//         << "species (including \"self\") is "
//         << vmrs[this_species]+broad_spec_vmr_sum
//         << ", more than 10% " << "different from 1.\n";
//         out2 << os.str();
//     }
    
    // Normalize foreign gamma and deltaf with the foreign VMR sum (but only if
    // we have any foreign broadening species):
    if (broad_spec_vmr_sum != 0.)
    {
        foreign_gamma /= broad_spec_vmr_sum;
        deltaf        /= broad_spec_vmr_sum;
    }
    else if (self_pressure > 0.)
        // If there are no foreign broadening species present, the best assumption
        // we can make is to use gamma_self in place of foreign_gamma. for deltaf
        // there is no equivalent solution, as we don't have a Delta_self and don't
        // know which other Delta we should apply (in this case delta_f gets 0,
        // which should be okayish):
    {
        foreign_gamma = gamma/self_pressure;
    }
    // It can happen that broad_spec_vmr_sum==0 AND p_self==0 (e.g., when p_grid
    // exceeds the given atmosphere and zero-padding is applied). In this case,
    // both gamma_foreign and deltaf are 0 and we leave it like that.
    
    // Multiply by pressure. For the width we take only the foreign pressure.
    // This is consistent with that we have scaled with the sum of all foreign
    // broadening VMRs. In this way we make sure that the total foreign broadening
    // scales with the total foreign pressure.
    foreign_gamma  *= foreign_pressure;
    
    // For the width, add foreign parts:
    gamma += foreign_gamma;
    
    // For the shift we simply take the total pressure, since there is no self part.
    deltaf *= pressure;
    
    // That's it, we're done.
}

// This is the temperature derivative of the broadening used by ARTSCAT-4; the "AP"-tag in ARTSCAT-5
void PressureBroadeningData::GetPlanetaryBroadening_dT(Numeric& dgamma_dT,
                                                    Numeric& ddeltaf_dT,
                                                    const Numeric& T,
                                                    const Numeric& T0,
                                                    const Numeric& pressure,
                                                    const Numeric& self_pressure,
                                                    const ArrayOfIndex& broad_spec_locations,
                                                    ConstVectorView vmrs) const
{
    const Numeric theta = T0/T;
    
    // Number of broadening species:
    const Index nbs = LineRecord::NBroadSpec();
    assert(nbs==broad_spec_locations.nelem());
    
    // Split total pressure in self and foreign part:
    const Numeric foreign_pressure = pressure - self_pressure;
    
    // Calculate sum of VMRs of all available foreign broadening species (we need this
    // for normalization). The species "Self" will not be included in the sum!
    Numeric broad_spec_vmr_sum = 0;
    
    // Gamma is the line width. We first initialize gamma with the self width
    dgamma_dT = - mdata[1][0] * mdata[0][0] * pow(theta, mdata[1][0]) * self_pressure / T;
    
    // Set deltaf to 0
    ddeltaf_dT = 0;
    
    // and treat foreign width separately:
    Numeric foreign_dgamma_dT = 0;
    
    // Add up foreign broadening species, where available:
    for (Index i=0; i<nbs; ++i) {
        if ( broad_spec_locations[i] < -1 ) {
            // -2 means that this broadening species is identical to Self.
            // Throw runtime errors if the parameters are not identical.
            if (mdata[2][i]!=mdata[0][0] ||
                mdata[3][i]!=mdata[1][0])
            {
                std::ostringstream os;
                os << "Inconsistency in LineRecord, self broadening and line "
                << "broadening for " << LineRecord::BroadSpecName(i) << "\n"
                << "should be identical.\n";
                throw std::runtime_error(os.str());
            }
        } else if ( broad_spec_locations[i] >= 0 ) {
            
            // Add to VMR sum:
            broad_spec_vmr_sum += vmrs[broad_spec_locations[i]];
            
            // foreign broadening:
            foreign_dgamma_dT += - mdata[3][i] * mdata[2][i] * pow(theta, mdata[3][i])
                * vmrs[broad_spec_locations[i]] / T;
            
            // Delta f (not .25+1.5*foreign_broadening)
            ddeltaf_dT += - ((Numeric).25 + (Numeric)1.5*mdata[3][i]) * mdata[4][i]
            * pow( theta , (Numeric).25 + (Numeric)1.5*mdata[3][i] )
            * vmrs[broad_spec_locations[i]] / T;
        }
    }
    
//     // Check that sum of self and all foreign VMRs is not too far from 1:
//     if ( abs(vmrs[this_species]+broad_spec_vmr_sum-1) > 0.1
//         && out2.sufficient_priority() )
//     {
//         std::ostringstream os;
//         os << "Warning: The total VMR of all your defined broadening\n"
//         << "species (including \"self\") is "
//         << vmrs[this_species]+broad_spec_vmr_sum
//         << ", more than 10% " << "different from 1.\n";
//         out2 << os.str();
//     }
    
    // Normalize foreign gamma and deltaf with the foreign VMR sum (but only if
    // we have any foreign broadening species):
    if (broad_spec_vmr_sum != 0.)
    {
        foreign_dgamma_dT /= broad_spec_vmr_sum;
        ddeltaf_dT        /= broad_spec_vmr_sum;
    }
    else if (self_pressure > 0.)
        // If there are no foreign broadening species present, the best assumption
        // we can make is to use gamma_self in place of foreign_gamma. for deltaf
        // there is no equivalent solution, as we don't have a Delta_self and don't
        // know which other Delta we should apply (in this case delta_f gets 0,
        // which should be okayish):
    {
        foreign_dgamma_dT = dgamma_dT/self_pressure;
    }
    // It can happen that broad_spec_vmr_sum==0 AND p_self==0 (e.g., when p_grid
    // exceeds the given atmosphere and zero-padding is applied). In this case,
    // both gamma_foreign and deltaf are 0 and we leave it like that.
    
    // Multiply by pressure. For the width we take only the foreign pressure.
    // This is consistent with that we have scaled with the sum of all foreign
    // broadening VMRs. In this way we make sure that the total foreign broadening
    // scales with the total foreign pressure.
    foreign_dgamma_dT  *= foreign_pressure;
    
    // For the width, add foreign parts:
    dgamma_dT += foreign_dgamma_dT;
    
    // For the shift we simply take the total pressure, since there is no self part.
    ddeltaf_dT *= pressure;
    
    // That's it, we're done.
}


///////////////////////////////////////////////////////////////////////////////////////
// Get SD broadening for air
///////////////////////////////////////////////////////////////////////////////////////

// This is the broadening used by the "SD-AIR"-tag in ARTSCAT-5
void PressureBroadeningData::GetSDAirBroadening(Numeric& gamma0,
                                                Numeric& gamma2,
                                                Numeric& delta0,
                                                Numeric& delta2,
                                                const Numeric& theta,
                                                const Numeric& pressure) const
{
    gamma0 = mdata[0][0] * pow(theta,mdata[1][0]) * pressure;
    
    gamma2 = mdata[0][1] * pow(theta,mdata[1][1]) * pressure;
    
    delta0 = mdata[0][2] * pow(theta,mdata[1][2]) * pressure;
    
    delta2 = mdata[0][3] * pow(theta,mdata[1][3]) * pressure;
}

// This is the broadening used by the "SD-AIR"-tag in ARTSCAT-5
void PressureBroadeningData::GetSDAirBroadening_dT(Numeric& dgamma0,
                                                   Numeric& dgamma2,
                                                   Numeric& ddelta0,
                                                   Numeric& ddelta2,
                                                   const Numeric& T,
                                                   const Numeric& T0,
                                                   const Numeric& pressure) const
{
  const Numeric theta = T0/T;
  
  dgamma0 = - (mdata[0][0] * pow(theta,mdata[1][0]) * pressure) / T * mdata[1][0];
  
  dgamma2 = - (mdata[0][1] * pow(theta,mdata[1][1]) * pressure) / T * mdata[1][1];
  
  ddelta0 = - (mdata[0][2] * pow(theta,mdata[1][2]) * pressure) / T * mdata[1][2];
  
  ddelta2 = - (mdata[0][3] * pow(theta,mdata[1][3]) * pressure) / T * mdata[1][3];
}



///////////////////////////////////////////////////////////////////////////////////////
// Get HTP broadening for air
///////////////////////////////////////////////////////////////////////////////////////

// This is the broadening used by the "HTP-AIR"-tag in ARTSCAT-5
void PressureBroadeningData::GetHTPAirBroadening(Numeric& gamma0,
                                                 Numeric& gamma2,
                                                 Numeric& delta0,
                                                 Numeric& delta2,
                                                 Numeric& fvc,
                                                 Numeric& eta,
                                                 const Numeric& theta,
                                                 const Numeric& pressure) const
{
  gamma0 = mdata[0][0] * pow(theta,mdata[1][0]) * pressure;
  
  gamma2 = mdata[0][1] * pow(theta,mdata[1][1]) * pressure;
  
  delta0 = mdata[0][2] * pow(theta,mdata[1][2]) * pressure;
  
  delta2 = mdata[0][3] * pow(theta,mdata[1][3]) * pressure;
  
  fvc = mdata[0][4] * pow(theta,mdata[1][4]) * pressure;
  
  eta = mdata[0][5] * pow(theta,mdata[1][5]) * pressure;
}


// This is the broadening used by the "HTP-AIR"-tag in ARTSCAT-5
void PressureBroadeningData::GetHTPAirBroadening_dT(Numeric& dgamma0,
                                                    Numeric& dgamma2,
                                                    Numeric& ddelta0,
                                                    Numeric& ddelta2,
                                                    Numeric& dfvc,
                                                    Numeric& deta,
                                                    const Numeric& T,
                                                    const Numeric& T0,
                                                    const Numeric& pressure) const
{
  const Numeric theta = T0/T;
  
  dgamma0 = - (mdata[0][0] * pow(theta,mdata[1][0]) * pressure) / T * mdata[1][0];
  
  dgamma2 = - (mdata[0][1] * pow(theta,mdata[1][1]) * pressure) / T * mdata[2][1];
  
  ddelta0 = - (mdata[0][2] * pow(theta,mdata[1][2]) * pressure) / T * mdata[1][2];
  
  ddelta2 = - (mdata[0][3] * pow(theta,mdata[1][3]) * pressure) / T * mdata[1][3];
  
  dfvc = - (mdata[0][4] * pow(theta,mdata[1][4]) * pressure) / T * mdata[1][4];
  
  deta = - (mdata[0][5] * pow(theta,mdata[1][5]) * pressure) / T * mdata[1][5];
}

///////////////////////////////////////////////////////////////////////////////////////
// Get SD broadening for air
///////////////////////////////////////////////////////////////////////////////////////

// This is the broadening used by the "Testing"-tag in ARTSCAT-5
void PressureBroadeningData::GetTestBroadening(Numeric& gamma0,
                                               Numeric& gamma2,
                                               Numeric& delta0,
                                               ConstVectorView vmrs,
                                               const Numeric& theta,
                                               const Numeric& pressure,
                                               const Index h2o_index) const
{
  if(h2o_index > -1)
  {
    gamma0 = (mdata[0][0] * pow(theta,mdata[1][0]) * (1 - vmrs[h2o_index]) + 
              mdata[2][0] * pow(theta,mdata[3][0]) *      vmrs[h2o_index]) * pressure ;
  
    delta0 = (mdata[4][0] * (1 - vmrs[h2o_index]) + 
              mdata[5][0] *      vmrs[h2o_index]) * pressure;
              
    gamma2 = (mdata[6][0] * (1 - vmrs[h2o_index]) + 
              mdata[7][0] *      vmrs[h2o_index]) * pressure;
  }
  else
  {
    gamma0 = mdata[0][0] * pow(theta,mdata[1][0]) * pressure;
    delta0 = mdata[4][0] * pressure;
    gamma2 = mdata[6][0] * pressure;
  }
  
}

///////////////////////////////////////////
//  Catalog interactions here
///////////////////////////////////////////

// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetAirBroadeningFromCatalog(const Numeric& sgam, 
                                                         const Numeric& nself,
                                                         const Numeric& agam,
                                                         const Numeric& nair,
                                                         const Numeric& air_pressure_DF,
                                                         const Numeric& dsgam, 
                                                         const Numeric& dnself,
                                                         const Numeric& dagam,
                                                         const Numeric& dnair,
                                                         const Numeric& dair_pressure_DF) 
{
    mtype = PB_AIR_BROADENING;
    mdata.resize(5);
    mdataerror.resize(5);
    for(Index ii=0;ii<5;ii++)
    {
        mdata[ii].resize(1);
        mdataerror[ii].resize(1);
    }
    mdata[0][0] = sgam;             // Self broadening gamma parameter
    mdata[1][0] = nself;            // Self broadening n parameter
    mdata[2][0] = agam;             // Air broadening gamma parameter
    mdata[3][0] = nair;             // Air broadening n parameter
    mdata[4][0] = air_pressure_DF;  // Pressure shift parameter
    
    mdataerror[0][0] = dsgam;             // Self broadening gamma parameter
    mdataerror[1][0] = dnself;            // Self broadening n parameter
    mdataerror[2][0] = dagam;             // Air broadening gamma parameter
    mdataerror[3][0] = dnair;             // Air broadening n parameter
    mdataerror[4][0] = dair_pressure_DF;  // Pressure shift parameter
}


// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetAirAndWaterBroadeningFromCatalog(const Numeric& sgam, 
                                                                 const Numeric& sn, 
                                                                 const Numeric& sdelta, 
                                                                 const Numeric& agam,
                                                                 const Numeric& an,
                                                                 const Numeric& adelta,
                                                                 const Numeric& wgam,
                                                                 const Numeric& wn,
                                                                 const Numeric& wdelta) 
{
    mtype = PB_AIR_AND_WATER_BROADENING;
    mdata.resize(3);
    mdataerror.resize(0);
    
    mdata[0].resize(3);
    mdata[0][0] = sgam;             // Self broadening gamma parameter
    mdata[0][1] = sn;               // Self broadening n parameter
    mdata[0][2] = sdelta;           // Self broadening shift parameter
    
    mdata[1].resize(3);
    mdata[1][0] = agam;             // Air broadening gamma parameter
    mdata[1][1] = an;               // Air broadening n parameter
    mdata[1][2] = adelta;           // Air broadening shift parameter
    
    mdata[2].resize(3);
    mdata[2][0] = wgam;             // Water broadening gamma parameter
    mdata[2][1] = wn;               // Water broadening n parameter
    mdata[2][2] = wdelta;           // Water broadening shift parameter
}


// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetPlanetaryBroadeningFromCatalog(const Numeric& sgam, 
                                                            const Numeric& nself,
                                                            const Vector&  foreign_gamma,
                                                            const Vector&  n_foreign,
                                                            const Vector&  foreign_pressure_DF) 
{
    // All vectors must have the same length
    assert(n_foreign.nelem()==foreign_gamma.nelem());
    assert(foreign_pressure_DF.nelem()==foreign_gamma.nelem());
    
    mtype = PB_PLANETARY_BROADENING;
    mdata.resize(5);
    mdataerror.resize(0);
    for(Index ii=0;ii<2;ii++) 
    {
        mdata[ii].resize(1);
    }
    mdata[0][0] = sgam;         // Self broadening gamma parameter
    mdata[1][0] = nself;        // Self broadening n parameter
    mdata[2] = foreign_gamma;   // Gas broadening gamma parameter per species
    mdata[3] = n_foreign;       // Gas broadening n parameter per species
    mdata[4] = foreign_pressure_DF;      // Pressure shift parameter per species
}

// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetSDAirFromCatalog(const Numeric& gamma0,
                                                 const Numeric& gamma0_exp,
                                                 const Numeric& gamma2,
                                                 const Numeric& gamma2_exp,
                                                 const Numeric& delta0,
                                                 const Numeric& delta0_exp,
                                                 const Numeric& delta2,
                                                 const Numeric& delta2_exp)
{
    mtype = PB_SD_AIR_VOLUME;
    mdata.resize(2);
    mdataerror.resize(0);
    for(Index ii=0;ii<2;ii++)
    {
        mdata[ii].resize(4);
    }
    
    mdata[0][0] = gamma0;
    mdata[1][0] = gamma0_exp;
    
    mdata[0][1] = gamma2;
    mdata[1][1] = gamma2_exp;
    
    mdata[0][2] = delta0;
    mdata[1][2] = delta0_exp;
    
    mdata[0][3] = delta2;
    mdata[1][3] = delta2_exp;
}

// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetHTPAirFromCatalog(const Numeric& gamma0,
                                                 const Numeric& gamma0_exp,
                                                 const Numeric& gamma2,
                                                 const Numeric& gamma2_exp,
                                                 const Numeric& delta0,
                                                 const Numeric& delta0_exp,
                                                 const Numeric& delta2,
                                                 const Numeric& delta2_exp,
                                                 const Numeric& fvc,
                                                 const Numeric& fvc_exp,
                                                 const Numeric& eta,
                                                 const Numeric& eta_exp)
{
  mtype = PB_HTP_AIR_VOLUME;
  mdata.resize(2);
  mdataerror.resize(0);
  for(Index ii=0;ii<2;ii++)
  {
    mdata[ii].resize(4);
  }
  
  mdata[0][0] = gamma0;
  mdata[1][0] = gamma0_exp;
  
  mdata[0][1] = gamma2;
  mdata[1][1] = gamma2_exp;
  
  mdata[0][2] = delta0;
  mdata[1][2] = delta0_exp;
  
  mdata[0][3] = delta2;
  mdata[1][3] = delta2_exp;
  
  mdata[0][4] = fvc;
  mdata[1][4] = fvc_exp;
  
  mdata[0][5] = eta;
  mdata[1][5] = eta_exp;
}

// Use these to insert the data in the required format from catalog readings
void PressureBroadeningData::SetTestFromCatalog(const Numeric& gamma0_air, 
                                                const Numeric& gamma0_air_exp,
                                                const Numeric& gamma0_water,
                                                const Numeric& gamma0_water_exp,
                                                const Numeric& gamma2_air,
                                                const Numeric& gamma2_water,
                                                const Numeric& delta0_air,
                                                const Numeric& delta0_water)
{
  mtype = PB_PURELY_FOR_TESTING;
  mdata.resize(8);
  mdataerror.resize(0);
  for(Index ii=0;ii<8;ii++)
  {
    mdata[ii].resize(1);
  }
  
  mdata[0][0] = gamma0_air;
  mdata[1][0] = gamma0_air_exp;
  
  mdata[2][0] = gamma0_water;
  mdata[3][0] = gamma0_water_exp;
  
  mdata[4][0] = delta0_air;
  mdata[5][0] = delta0_water;
  
  mdata[6][0] = gamma2_air;
  mdata[7][0] = gamma2_water;

}

///////////////////////////////////////////
// Change of internals formats
///////////////////////////////////////////

void PressureBroadeningData::ChangeSelf(const Numeric& change, 
                                        const Index this_species, 
                                        const Index h2o_species, 
                                        const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[0][0]+=change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[0][0]+=change;
            if(this_species==h2o_species)
                mdata[2][0]+=change;
            break;
        case PB_PLANETARY_BROADENING:
            mdata[0][0]+=change;
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]==-2)
                    mdata[2][ii]+=change;
            break;
        default:
          throw std::runtime_error("ChangeSelf: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeSelfExponent(const Numeric& change, 
                                                const Index this_species, 
                                                const Index h2o_species, 
                                                const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[1][0]+=change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[0][1]+=change;
            if(this_species==h2o_species)
                mdata[2][1]+=change;
            break;
        case PB_PLANETARY_BROADENING:
            mdata[1][0]+=change;
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]==-2)
                    mdata[3][ii]+=change;
            break;
        default:
          throw std::runtime_error("ChangeSelfExponent: Cannot recognize type");
    }
}
void PressureBroadeningData::SetSelf(const Numeric& new_value, 
                                     const Index this_species, 
                                     const Index h2o_species, 
                                     const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[0][0]=new_value;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[0][0]=new_value;
            if(this_species==h2o_species)
                mdata[2][0]=new_value;
            break;
        case PB_PLANETARY_BROADENING:
            mdata[0][0]=new_value;
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]==-2)
                    mdata[2][ii]=new_value;
            break;
        default:
          throw std::runtime_error("SetSelf: Cannot recognize type");
    }
}
void PressureBroadeningData::SetSelfExponent(const Numeric& new_value, 
                                             const Index this_species, 
                                             const Index h2o_species, 
                                             const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[1][0]=new_value;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[0][1]=new_value;
            if(this_species==h2o_species)
                mdata[2][1]=new_value;
            break;
        case PB_PLANETARY_BROADENING:
            mdata[1][0]=new_value;
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]==-2)
                    mdata[3][ii]=new_value;
            break;
        default:
          throw std::runtime_error("SetSelfExponent: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeSelfRelative(const Numeric& change, 
                                                const Index this_species, 
                                                const Index h2o_species, 
                                                const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[0][0]*=1.0e0+change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[0][0]*=1.0e0+change;
            if(this_species==h2o_species)
                mdata[2][0]*=1.0e0+change;
            break;
        case PB_PLANETARY_BROADENING:
            mdata[0][0]*=1.0e0+change;
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]==-2)
                    mdata[2][ii]*=1.0e0+change;
            break;
        default:
          throw std::runtime_error("ChangeSelfRelative: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeSelfExponentRelative(const Numeric& change, 
                                                        const Index this_species, 
                                                        const Index h2o_species, 
                                                        const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[1][0]*=1.0e0+change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[0][1]*=1.0e0+change;
            if(this_species==h2o_species)
                mdata[2][1]*=1.0e0+change;
            break;
        case PB_PLANETARY_BROADENING:
            mdata[1][0]*=1.0e0+change;
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]==-2)
                    mdata[3][ii]*=1.0e0+change;
            break;
        default:
          throw std::runtime_error("ChangeSelfExponentRelative: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeForeign(const Numeric& change, 
                                           const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[2][0]+=change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[1][0]+=change;
            break;
        case PB_PLANETARY_BROADENING:
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]!=-2)
                    mdata[2][ii]+=change;
            break;
        default:
          throw std::runtime_error("ChangeForeign: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeForeignExponent(const Numeric& change, 
                                                   const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[3][0]+=change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[1][1]+=change;
            break;
        case PB_PLANETARY_BROADENING:
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]!=-2)
                    mdata[3][ii]+=change;
            break;
        default:
          throw std::runtime_error("ChangeForeignExponent: Cannot recognize type");
    }
}
void PressureBroadeningData::SetForeign(const Numeric& new_value, 
                                        const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[2][0]=new_value;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[1][0]=new_value;
            break;
        case PB_PLANETARY_BROADENING:
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]!=-2)
                    mdata[2][ii]=new_value;
            break;
        default:
          throw std::runtime_error("SetForeign: Cannot recognize type");
    }
}
void PressureBroadeningData::SetForeignExponent(const Numeric& new_value, 
                                                const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[3][0]=new_value;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[1][1]=new_value;
            break;
        case PB_PLANETARY_BROADENING:
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]!=-2)
                    mdata[3][ii]=new_value;
            break;
        default:
          throw std::runtime_error("SetForeignExponent: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeForeignRelative(const Numeric& change, 
                                                   const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[2][0]*=1.0e0+change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[1][0]*=1.0e0+change;
            break;
        case PB_PLANETARY_BROADENING:
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]!=-2)
                    mdata[2][ii]*=1.0e0+change;
            break;
        default:
          throw std::runtime_error("ChangeForeignRelative: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeForeignExponentRelative(const Numeric& change, 
                                                           const ArrayOfIndex& broad_spec_locations)
{
    switch(mtype)
    {
        case PB_NONE:
            // Note that this is oftentimes not wanted, but a valid case at low pressures
            break;
        case PB_AIR_BROADENING:
            mdata[3][0]*=1.0e0+change;
            break;
        case PB_AIR_AND_WATER_BROADENING:
            mdata[1][1]*=1.0e0+change;
            break;
        case PB_PLANETARY_BROADENING:
            for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
                if(broad_spec_locations[ii]!=-2)

                    mdata[3][ii]*=1.0e0+change;
            break;
        default:
          throw std::runtime_error("ChangeForeignExponentRelative: Cannot recognize type");
    }
}
void PressureBroadeningData::ChangeForeignShiftRelative(const Numeric& change, 
                                                   const ArrayOfIndex& broad_spec_locations)
{
  switch(mtype)
  {
    case PB_NONE:
      // Note that this is oftentimes not wanted, but a valid case at low pressures
      break;
    case PB_AIR_BROADENING:
      mdata[4][0]*=1.0e0+change;
      break;
    case PB_AIR_AND_WATER_BROADENING:
      mdata[1][2]*=1.0e0+change;
      break;
    case PB_PLANETARY_BROADENING:
      for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
        if(broad_spec_locations[ii]!=-2)
          mdata[4][ii]*=1.0e0+change;
      break;
    default:
      throw std::runtime_error("ChangeForeignRelative: Cannot recognize type");
  }
}
void PressureBroadeningData::ChangeForeignShift(const Numeric& change, 
                                                const ArrayOfIndex& broad_spec_locations)
{
  switch(mtype)
  {
    case PB_NONE:
      // Note that this is oftentimes not wanted, but a valid case at low pressures
      break;
    case PB_AIR_BROADENING:
      mdata[4][0]+=change;
      break;
    case PB_AIR_AND_WATER_BROADENING:
      mdata[1][2]+=change;
      break;
    case PB_PLANETARY_BROADENING:
      for(Index ii=0;ii<broad_spec_locations.nelem();ii++)
        if(broad_spec_locations[ii]!=-2)
          mdata[4][ii]+=change;
      break;
    default:
      throw std::runtime_error("ChangeForeignRelative: Cannot recognize type");
  }
}

///////////////////////////////////////////
// Formating and readings here
///////////////////////////////////////////

Index PressureBroadeningData::ExpectedVectorLengthFromType() const
{
  switch(mtype) {
    case PB_NONE:
      return 0;
    case PB_AIR_BROADENING:
      return 10;
    case PB_AIR_AND_WATER_BROADENING:
      return 9;
    case PB_PLANETARY_BROADENING:
      return 20;
    case PB_PURELY_FOR_TESTING:
    case PB_SD_AIR_VOLUME:
      return 8;
    case PB_HTP_AIR_VOLUME:
      return 12;
    case PB_VOIGT_TEST_WATER:
    case PB_SD_TEST_WATER:
      return Index(TestParams::COUNT);
  }
  return -1;
}

void PressureBroadeningData::SetDataFromVectorWithKnownType(ConstVectorView input)
{
  if(input.nelem()!=ExpectedVectorLengthFromType())
      throw std::runtime_error("Input pressure broadening is of wrong length.\n");
  
  switch(mtype)
  {
    case PB_NONE:
      // The none case
      mdata.resize(0);
      mdataerror.resize(0);
      break;
    case PB_AIR_BROADENING:
      SetAirBroadeningFromCatalog(input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7], input[8], input[9]);
      break;
    case PB_AIR_AND_WATER_BROADENING:
      SetAirAndWaterBroadeningFromCatalog(input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7], input[8]);
      break;
    case PB_PLANETARY_BROADENING:
      SetPlanetaryBroadeningFromCatalog(input[0], input[7], input[Range(1,6)], input[Range(8,6)], input[Range(14,6)]);
      break;
    case PB_PURELY_FOR_TESTING:
    case PB_SD_AIR_VOLUME:
      SetSDAirFromCatalog(input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7]); 
      break;
    case PB_HTP_AIR_VOLUME:
      SetHTPAirFromCatalog(input[0], input[1], input[2], input[3], input[4], input[5], input[6], input[7], input[8], input[9], input[10], input[11]);
      break;
    case PB_VOIGT_TEST_WATER:
    case PB_SD_TEST_WATER:
      mdata.resize(1);
      mdata[0] = input;
      break;
  }
}

void PressureBroadeningData::StorageTag2SetType(const String & input)
{
    if(input == "NA") // The none case
        mtype=PB_NONE;
    else if(input == "N2") // Air Broadening is N2 broadening mostly...
        mtype=PB_AIR_BROADENING;
    else if(input == "WA") // Water and Air Broadening
        mtype=PB_AIR_AND_WATER_BROADENING;
    else if(input == "AP") // Planetary broadening
        mtype=PB_PLANETARY_BROADENING;
    else if(input == "SD-AIR") 
      mtype=PB_SD_AIR_VOLUME;
    else if(input == "HTP-AIR") 
      mtype=PB_HTP_AIR_VOLUME;
    else if(input == "TESTING") 
      mtype=PB_PURELY_FOR_TESTING;
    else if(input == "PB_SD_TEST_WATER") 
      mtype=PB_SD_TEST_WATER;
    else if(input == "PB_VOIGT_TEST_WATER") 
      mtype=PB_VOIGT_TEST_WATER;
    else
      throw std::runtime_error("StorageTag2SetType: Cannot recognize tag.\n");
}

void PressureBroadeningData::GetVectorFromData(Vector& output) const 
{
  const Index n = ExpectedVectorLengthFromType();
  output.resize(n);
  
  switch(mtype) {
    case PB_NONE:
      break;
    case PB_AIR_BROADENING:
      output[0]=mdata[0][0];
      output[1]=mdata[1][0];
      output[2]=mdata[2][0];
      output[3]=mdata[3][0];
      output[4]=mdata[4][0];
      output[5]=mdataerror[0][0];
      output[6]=mdataerror[1][0];
      output[7]=mdataerror[2][0];
      output[8]=mdataerror[3][0];
      output[9]=mdataerror[4][0];
      break;
    case PB_AIR_AND_WATER_BROADENING:
      output[0]=mdata[0][0];
      output[1]=mdata[0][1];
      output[2]=mdata[0][2];
      output[3]=mdata[1][0];
      output[4]=mdata[1][1];
      output[5]=mdata[1][2];
      output[6]=mdata[2][0];
      output[7]=mdata[2][1];
      output[8]=mdata[2][2];
      break;
    case PB_PLANETARY_BROADENING:
      output[0]=mdata[0][0];
      output[7]=mdata[1][0];
      output[Range(1,6)]=mdata[2];
      output[Range(8,6)]=mdata[3];
      output[Range(14,6)]=mdata[4];
      break;
    case PB_SD_AIR_VOLUME:
    case PB_HTP_AIR_VOLUME:
      for(Index i = 0; i < n; i++) {
        if(i%2)
          output[i] = mdata[1][i/2];
        else
          output[i] = mdata[0][i/2];
      }
      break;
    case PB_PURELY_FOR_TESTING:
      for(Index i = 0; i < n; i++)
        output[i] = mdata[i][0];
      break;
    case PB_VOIGT_TEST_WATER:
    case PB_SD_TEST_WATER:
      output = mdata[0];
      break;
  }
}


String PressureBroadeningData::Type2StorageTag() const
{
  String output;

  switch(mtype) {
    case PB_NONE:
      return "NA";
    case PB_AIR_BROADENING:
      return "N2";
    case PB_AIR_AND_WATER_BROADENING:
      return "WA";
    case  PB_SD_AIR_VOLUME:
      return "SD-AIR";
    case  PB_HTP_AIR_VOLUME:
      return "HTP-AIR";
    case PB_PLANETARY_BROADENING:
      return "AP";
    case PB_SD_TEST_WATER:
      return "PB_SD_TEST_WATER";
    case PB_VOIGT_TEST_WATER:
      return "PB_VOIGT_TEST_WATER";
    case PB_PURELY_FOR_TESTING:
      throw std::runtime_error("Cannot save pure testing version");
  }
  return "-1";
}

// T0 temperature type is just a constant
inline Numeric main_t0(const Numeric& x0) noexcept { return x0; }
inline Numeric dmain_dx0_t0() noexcept { return 1; }
inline Numeric dmain_dx1_t0() noexcept { return 0; }
inline Numeric dmain_dx2_t0() noexcept { return 0; }
inline Numeric dmain_dT0_t0() noexcept { return 0; }
inline Numeric dmain_dT__t0() noexcept { return 0; }

// T1 temperature type is standard HITRAN x0 (T0 / T) ^ x1
inline Numeric main_t1(const Numeric& TH, const Numeric& x0, const Numeric& x1) noexcept { return x0 * pow(TH, x1); }
inline Numeric dmain_dx0_t1(const Numeric& main, const Numeric& x0) noexcept { return main / x0; }
inline Numeric dmain_dx1_t1(const Numeric& main, const Numeric& TH) noexcept { return main * log(TH); }
inline Numeric dmain_dx2_t1() noexcept { return 0; }
inline Numeric dmain_dT0_t1(const Numeric& main, const Numeric& T0, const Numeric& x1) noexcept { return x1 * main / T0; }
inline Numeric dmain_dT__t1(const Numeric& main, const Numeric& T, const Numeric& x1) noexcept { return - dmain_dT0_t1(main, T, x1); }

// T2 temperature type is proposed for line shifts x0 (T0 / T) ^ x1 / (1 + x2 ln(T / T0))
inline Numeric main_t2(const Numeric& TH, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return x0 * pow(TH, x1) * (1 + x2 * log(1/TH)); }
inline Numeric dmain_dx0_t2(const Numeric& main, const Numeric& x0) noexcept { return dmain_dx0_t1(main, x0); }
inline Numeric dmain_dx1_t2(const Numeric& main, const Numeric& TH) noexcept { return dmain_dx1_t1(main, TH); }
inline Numeric dmain_dx2_t2(const Numeric& main, const Numeric& TH, const Numeric& x2) noexcept { return -log(1/TH)/(x2*log(1/TH) + 1) * main; }
inline Numeric dmain_dT0_t2(const Numeric& main, const Numeric& TH, const Numeric& T0, const Numeric& x1, const Numeric& x2) noexcept { return - (x1*(x2*log(1/TH) + 1) + x2)/(T0*(x1*log(1/TH) + 1)) * main; }
inline Numeric dmain_dT__t2(const Numeric& main, const Numeric& TH, const Numeric& T, const Numeric& x1, const Numeric& x2) noexcept { return (x1*(x2*log(1/TH) + 1) + x2)/(T*(x1*log(1/TH) + 1)) * main; }

// T3 temperature type is proposed for speed-dependent parameters x0 + x1 (T - T0)
inline Numeric main_t3(const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1) noexcept { return x0 + x1 * (T - T0); }
inline Numeric dmain_dx0_t3() noexcept { return 1; }
inline Numeric dmain_dx1_t3(const Numeric& T, const Numeric& T0) noexcept { return (T - T0); }
inline Numeric dmain_dx2_t3() noexcept { return 0; }
inline Numeric dmain_dT0_t3(const Numeric& x1) noexcept { return -x1; }
inline Numeric dmain_dT__t3(const Numeric& x1) noexcept { return x1; }

// T4 temperature type is proposed for line mixing of second order (x0 + x1 (T0 / T - 1)) (T0 / T) ^ x2
inline Numeric main_t4(const Numeric& TH, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return (x0 + x1 * (TH - 1)) * pow(TH, x2); }
inline Numeric dmain_dx0_t4(const Numeric& TH, const Numeric& x2) noexcept { return pow(TH, x2); }
inline Numeric dmain_dx1_t4(const Numeric& TH, const Numeric& x2) noexcept { return (TH - 1) * pow(TH, x2); }
inline Numeric dmain_dx2_t4(const Numeric& main, const Numeric& TH) noexcept { return main * log(TH); }
inline Numeric dmain_dT0_t4(const Numeric& main, const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return main * (x1*T0 + x2*(x0*T + x1*(T0 - T)))/(T*(x0*T + x1*(T0 - T))); }
inline Numeric dmain_dT__t4(const Numeric& main, const Numeric& T, const Numeric& T0, const Numeric& x0, const Numeric& x1, const Numeric& x2) noexcept { return main * (x1*T0 + x2*(x0*T + x1*(T0 - T)))/(T0*(x0*T + x1*(T0 - T))); }

//! Select line parameter based on parameter order
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param G0    The speed-independent pressure broadening term
 * \param D0    The speed-independent pressure shift term
 * \param G2    The speed-dependent pressure broadening term
 * \param D2    The speed-dependent pressure shift term
 * \param FVC   The velocity-changing collisional parameter
 * \param ETA   The correlation between velocity and rotational changes due to collisions
 * \param param Index enumerating the parameter
 * \param type  Enum class for knowing which parameter we are interested in
 * 
 * \return A reference pointing at G0, D0, G2, D2, FVC, or ETA.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
inline Numeric& select_line_shape_param(Numeric& G0, Numeric& D0, Numeric& G2, Numeric& D2,  Numeric& FVC, Numeric& ETA, const Index param, const LineFunctionData::LineShapeType type) {
  switch(type) {
    case LineFunctionData::LineShapeType::DP:
      return G0;
    case LineFunctionData::LineShapeType::LP:
      switch(param) {
        case Index(LineFunctionData::LorentzParam::G0): return G0;
        case Index(LineFunctionData::LorentzParam::D0): return D0;
      }
      break;
    case LineFunctionData::LineShapeType::VP:
      switch(param) {
        case Index(LineFunctionData::VoigtParam::G0): return G0;
        case Index(LineFunctionData::VoigtParam::D0): return D0;
      }
      break;
    case LineFunctionData::LineShapeType::SDVP:
      switch(param) {
        case Index(LineFunctionData::SpeedVoigtParam::G0): return G0;
        case Index(LineFunctionData::SpeedVoigtParam::D0): return D0;
        case Index(LineFunctionData::SpeedVoigtParam::G2): return G2;
        case Index(LineFunctionData::SpeedVoigtParam::D2): return D2;
      }
      break;
    case LineFunctionData::LineShapeType::HTP:
      switch(param) {
        case Index(LineFunctionData::HTPParam::G0):  return G0;
        case Index(LineFunctionData::HTPParam::D0):  return D0;
        case Index(LineFunctionData::HTPParam::G2):  return G2;
        case Index(LineFunctionData::HTPParam::D2):  return D2;
        case Index(LineFunctionData::HTPParam::FVC): return FVC;
        case Index(LineFunctionData::HTPParam::ETA): return ETA;
      }
      break;
  }
  throw std::runtime_error("Something very bad has happened, your type is either not defined or you are trying to access this function in a way that is not permitted");
}

//! Select line parameter based on parameter order
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param Y     The imaginary part of the line shape due to line mixing
 * \param G     The strength-altering due to line mixing
 * \param DV    The frequency shift due to line mixing
 * \param param Index enumerating the parameter
 * \param type  Enum class for knowing which parameter we are interested in
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
inline Numeric& select_line_mixing_param(Numeric& Y, Numeric& G, Numeric& DV, const Index param, const LineFunctionData::LineMixingType type) {
  switch(type) {
    case LineFunctionData::LineMixingType::None:
      return Y;
    case LineFunctionData::LineMixingType::Interp:
      return Y;
    case LineFunctionData::LineMixingType::LM1:
      switch(param) {
        case Index(LineFunctionData::FirstOrderParam::Y): return Y;
      }
      break;
    case LineFunctionData::LineMixingType::LM2:
      switch(param) {
        case Index(LineFunctionData::SecondOrderParam::Y): return Y;
        case Index(LineFunctionData::SecondOrderParam::G): return G;
        case Index(LineFunctionData::SecondOrderParam::DV): return DV;
      }
      break;
  }
  throw std::runtime_error("Something very bad has happened, your type is either not defined or you are trying to access this function in a way that is not permitted");
}


//! Special function for line mixing of LBLRTM type
/*! 
 * LBLRTM interpolates a list of parameters in their data to a given
 * range.  This method keeps track of the positions of the data to
 * complete the interpolation.
 * 
 * \param Y           The imaginary part of the line shape due to line mixing
 * \param G           The strength-altering due to line mixing
 * \param partial_vmr The VMR of the species in question
 * \param T           The atmospheric temperature
 * \param data        The associated data struction [T1, T2, T3, T4, Y1, Y2, Y3, Y4, G1, G2, G3, G4]
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
inline void special_line_mixing_aer(Numeric&Y, Numeric&G, const Numeric& partial_vmr, const Numeric& T, ConstVectorView data) {
  const Index n=data.nelem();
  if(n not_eq 12)
    throw std::runtime_error("Data for line mixing not matching AER LM style despite marked as such...\n must have this structure: [T1, T2, T3, T4, Y1, Y2, Y3, Y4, G1, G2, G3, G4]");
  
  // Perform standard lagrangian interpolation... FIXME: Test that this works
  Numeric Yx=0, Gx=0;
  for(Index i=0; i<4; i++) {
    Numeric xi=partial_vmr;
    #pragma omp simd
    for(Index j=0; j<4; j++)
      if(j not_eq i)
        xi *= (T - data[j]) / (data[i] - data[j]);
    Yx += xi * data[4 + i];
    Gx += xi * data[8 + i];
  }
  
  Y += Yx;
  G += Gx;
}


//! Compute the pressure broadening and line mixing parameters
/*! 
 * This helper code simply selects which of the input the
 * param-index is pointing at as a function of the type of
 * line shape under consideration
 * 
 * \param G0          The speed-independent pressure broadening term
 * \param D0          The speed-independent pressure shift term
 * \param G2          The speed-dependent pressure broadening term
 * \param D2          The speed-dependent pressure shift term
 * \param FVC         The velocity-changing collisional parameter
 * \param ETA         The correlation between velocity and rotational changes due to collisions
 * \param Y           The imaginary part of the line shape due to line mixing
 * \param G           The strength-altering due to line mixing
 * \param DV          The frequency shift due to line mixing
 * \param T0          The line reference temperature
 * \param T           The atmospheric temperature
 * \param P           The atmospheric pressure
 * \param self_vmr    The VMR of the species pressure
 * \param rtp_vmr     The VMR of all species in the atmosphere at a specific radiative transfer point
 * \param abs_species The list of all species in the atmosphere
 * 
 * \return A reference pointing at Y, G, or DV.
 * 
 * \author Richard Larsson
 * \date   2018-07-16
 */
void LineFunctionData::GetParams(Numeric& G0, Numeric& D0, 
                                 Numeric& G2, Numeric& D2, 
                                 Numeric& FVC, Numeric& ETA,
                                 Numeric& Y, Numeric& G, Numeric& DV,
                                 const Numeric& T0, const Numeric& T,
                                 const Numeric& P, const Numeric& self_vmr,
                                 const ConstVectorView rtp_vmr, const ArrayOfSpeciesTag& abs_species,
                                 const bool normalization) const
{
  // Set to zero
  G0 = D0 = G2 = D2 = FVC = ETA = Y = G = DV = 0;
  
  // Doppler broadening has no types...
  if(mp == LineShapeType::DP) return;
  
  // Theta is useful in many models so it gets precomputed
  const Numeric theta = T0/T;
  
  // Set up holders for partial and total VMR
  Numeric total_vmr=0, partial_vmr;
  
  // Add all species up
  for(Index i=0; i<mspecies.nelem(); i++) {
    if(i == 0 and mself)  // The first value might be self-broadening (use self_vmr)
      partial_vmr = self_vmr;
    else if(i == mspecies.nelem()-1 and mbath)  // The last value might be air-broadening (set total_vmr to 1)
      partial_vmr = 1 - total_vmr;
    else {  // Otherwise we have to find the species in the list of species tags
      Index this_species=-1;
      for(Index j=0; j<abs_species.nelem(); j++)
        if(abs_species[j].Species() == mspecies[i].Species())
          this_species = j;
        
      // If we cannot find the species, it does not exist and we are in a different atmosphere than covered by the line data
      if(this_species == -1) continue;  // Species does not exist
      partial_vmr = rtp_vmr[this_species];
    }
    // Sum up VMR
    total_vmr += partial_vmr;
    
    // Set up a current counter to keep track of the data
    Index current = 0;
    
    // Do the line shape parameters
    for(Index j=0; j<LineShapeTypeNelem(); j++) {
      // Selects one of G0, D0, G2, D2, FVC, ETA based on the order of the data for this line shape type
      Numeric& param = select_line_shape_param(G0, D0, G2, D2, FVC, ETA, j, mp);
      
      switch(mtypes[i][j]) {
        case TemperatureType::None:
          break;
        case TemperatureType::T0:
          param += partial_vmr * main_t0(mdata[i][current]);
          break;
        case TemperatureType::T1:
          param += partial_vmr * main_t1(theta, mdata[i][current], mdata[i][current+1]);
          break;
        case TemperatureType::T2:
          param += partial_vmr * main_t2(theta, mdata[i][current], mdata[i][current+1], mdata[i][current+2]);
          break;
        case TemperatureType::T3:
          param += partial_vmr * main_t3(T, T0, mdata[i][current], mdata[i][current+1]);
          break;
        case TemperatureType::T4:
          param += partial_vmr * main_t4(theta, mdata[i][current], mdata[i][current+1], mdata[i][current+2]);
          break;
        case TemperatureType::LM_AER:
          throw std::runtime_error("Not allowed for line shape parameters");
      }
      current += TemperatureTypeNelem(mtypes[i][j]);
    }
    
    // Do the line mixing parameters
    for(Index j=0; j<LineMixingTypeNelem(); j++) {
      Numeric& param = select_line_mixing_param(Y, G, DV, j, mlm);
      switch(mtypes[i][j]) {
        case TemperatureType::None:
          break;
        case TemperatureType::T0:
          param += partial_vmr * main_t0(mdata[i][current]);
          break;
        case TemperatureType::T1:
          param += partial_vmr * main_t1(theta, mdata[i][current], mdata[i][current+1]);
          break;
        case TemperatureType::T2:
          param += partial_vmr * main_t2(theta, mdata[i][current], mdata[i][current+1], mdata[i][current+2]);
          break;
        case TemperatureType::T3:
          param += partial_vmr * main_t3(T, T0, mdata[i][current], mdata[i][current+1]);
          break;
        case TemperatureType::T4:
          param += partial_vmr * main_t4(theta, mdata[i][current], mdata[i][current+1], mdata[i][current+2]);
          break;
        case TemperatureType::LM_AER:
          // Special case adding to both G and Y and nothing else!
          special_line_mixing_aer(Y, G, partial_vmr, T, mdata[i][Range(current, joker)]);
          break;
      }
      current += TemperatureTypeNelem(mtypes[i][j]);
    }
  }
  
  // If there was no species at all, then we cannot compute the pressure broadening
  if(total_vmr == 0) return;
  
  // Rescale to pressure with or without normalization
  const Numeric scale = normalization ? P / total_vmr : P;
  G0  *= scale;
  D0  *= scale;
  G2  *= scale;
  D2  *= scale;
  FVC *= scale;
  if(normalization) ETA /= total_vmr;
  Y   *= scale;
  G   *= scale * P;
  DV  *= scale * P;
}


//! Prints data so as operator>> can read it
std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd) {
  // Keep track of the size of the problem
  const Index nshapes=lfd.LineShapeTypeNelem(), nmixing=lfd.LineMixingTypeNelem();
  
  // Init the problem
  os << lfd.LineShapeType2String() << " " << lfd.LineMixingType2String() << " " << lfd.mspecies.nelem() << " ";
  
  // For all species we should now do mostly the same
  for(Index i=0; i<lfd.mspecies.nelem(); i++) {
    if(i==0 and lfd.mself)
      os << SELF << " ";  // SELF can be the first species
    else if(i==lfd.mspecies.nelem()-1 and lfd.mbath)
      os << BATH << " ";  // BATH can be the last species
    else
      os << lfd.mspecies[i].SpeciesNameMain() << " ";  // Otherwise we have a species defined in the assoc. SpeciesTag
    
    // Now we must count the data
    Index counter=0;
    
    // For everything that relates to shapes, do the same thing
    for(Index j=0; j<nshapes; j++) {
      os << lfd.TemperatureType2String(lfd.mtypes[i][j]) << " ";  // Name the temperature type
      for(Index k=0; k<lfd.TemperatureTypeNelem(lfd.mtypes[i][j]); k++)
        os << lfd.mdata[i][counter+k] << " ";  // Print the assoc. data
      counter += lfd.TemperatureTypeNelem(lfd.mtypes[i][j]);  // Count how much data has been printed
    }
    
    // For mixing we must take care with the interp. values but otherwise we act the same as prev. loop
    for(Index j=nshapes; j<nmixing+nshapes; j++) {
      os << lfd.TemperatureType2String(lfd.mtypes[i][j]) << " ";
      if(lfd.mlm == lfd.LineMixingType::Interp) {
        counter++;
        os << lfd.mdata[i].nelem() - counter << " ";
        for(; counter < lfd.mdata[i].nelem(); counter++)
          os << lfd.mdata[i][counter] << " ";
      }
      else {
        for(Index k=0; k<lfd.TemperatureTypeNelem(lfd.mtypes[i][j]); k++)
          os << lfd.mdata[i][counter+k] << " ";
        counter += lfd.TemperatureTypeNelem(lfd.mtypes[i][j]);
      }
    }
  }
  
  return os;
}

//! Reads data as created by operator<< can read it
std::istream& operator>>(std::istream& data, LineFunctionData& lfd) {
  lfd.mself = lfd.mbath = false;
  Index specs, c;
  Numeric n;
  String s;
  
  // The first tag should give the line shape scheme
  data >> s;
  lfd.StringSetLineShapeType(s);
  
  // The second tag should give the line mixing scheme
  data >> s;
  lfd.StringSetLineMixingType(s);
  
  // From the line shape and line mixing types, we know how many parameters are needed
  const Index count = lfd.LineShapeTypeNelem() + lfd.LineMixingTypeNelem();
  
  // The third tag should contain the number of species
  data >> specs;
  lfd.mspecies.resize(specs);
  lfd.mtypes = Array<Array<LineFunctionData::TemperatureType>>(specs, Array<LineFunctionData::TemperatureType>(count));
  lfd.mdata.resize(specs);
  
  if(lfd.mp not_eq LineFunctionData::LineShapeType::DP and not specs)
    throw std::runtime_error("Need at least one species for non-Doppler line shapes");
  
  // For all species, we need to set the methods to compute them
  for(Index i=0; i<specs; i++) {
    
    // This should be a species tag or one of the specials, SELF or BATH
    data >> s;
    if(s == SELF) {
      // If the species is self, then  we need to flag this
      lfd.mself = true;
      if(i not_eq 0)  // but self has to be first for consistent behavior
        throw std::runtime_error("Self broadening must be first, it is not\n");
    }
    else if(s == BATH) {
      // If the species is air, then we need to flag this
      lfd.mbath = true;
      if(i not_eq specs - 1)  // but air has to be last because it needs the rest's VMR
        throw std::runtime_error("Air/bath broadening must be last, it is not\n");
    }
    else {
      // Otherwise, we hope we find a species
      try {
        lfd.mspecies[i] = SpeciesTag(s);
      }
      catch(const std::runtime_error& e) {
        ostringstream os;
        os << "Encountered " << s << " in a position where a species should have been ";
        os << "defined.\nPlease check your pressure broadening data structure and ensure ";
        os << "that it follows the correct conventions.\n";
        os << "SpeciesTag error reads:  " << e.what();
        throw std::runtime_error(os.str());
      }
    }
    
    ArrayOfNumeric nums; nums.reserve(20); // buffers
    
    // For all parameters
    for(Index j=0; j<count; j++) {
      data >> s;  // Should contain a temperature tag
      lfd.StringSetTemperatureType(i, j, s);
      
      // Find the count of species, negative numbers means the next data value has this number
      c = lfd.TemperatureTypeNelem(lfd.mtypes[i][j]);
      if(c < 0)
        data >> c;
      
      // Add all new numbers for this line to the numbers
      for(Index k=0; k<c; k++) {
        data >> double_imanip() >> n;
        nums.push_back(n);
      }
    }
    
    // Set the data now that we know how many counts are required
    lfd.mdata[i].resize(nums.nelem());
    for(Index j=0; j<nums.nelem(); j++)
      lfd.mdata[i][j] = nums[j];
  }
  return data;
}

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

#include "zeemandata.h"
#include "species_info.h"


extern const Numeric PLANCK_CONST;
extern const Numeric BOHR_MAGNETON;
static const Numeric ZeemanSplittingConstant = BOHR_MAGNETON / PLANCK_CONST;


inline Numeric caseB(const Rational& N,
                     const Rational& J,
                     const Rational& Lambda,
                     const Rational& S,
                     const Numeric& GS,
                     const Numeric& GL) noexcept
{
  const Rational JJ = J*(J+1);
  const Rational NN = N*(N+1);
  const Rational SS = S*(S+1);
  const Rational LL = Lambda*Lambda;
  const Rational T1 = (JJ + SS - NN)           / JJ / 2;
  const Rational T2 = (JJ - SS + NN) * LL / NN / JJ / 2;
  
  Numeric g;
  if(JJ == 0)
    g = 0.0;
  else if(NN not_eq 0)
    g = GS * T1.toNumeric() + GL * T2.toNumeric();
  else
    g = GS * T1.toNumeric();
  return g;
}


inline Numeric caseA(const Rational& Omega,
                     const Rational& J,
                     const Rational& Lambda,
                     const Rational& Sigma,
                     const Numeric& GS,
                     const Numeric& GL) noexcept
{  
  const Rational JJ = J*(J+1);
  const Rational DIV = Omega / JJ;
  const Rational T1 = Sigma  * DIV;
  const Rational T2 = Lambda * DIV;
  
  Numeric g;
  if(JJ == 0)
    g = 0.0;
  else
    g = GS * T1.toNumeric() + GL * T2.toNumeric();
  return g; 
}


inline Numeric gHund(const QuantumNumbers& qns, const Numeric& GS, const Numeric& GL)
{
  switch(qns[QuantumNumberType::Hund].toIndex()) {
    case Index(Hund::CaseA):
      return caseA(qns[QuantumNumberType::Omega], 
                   qns[QuantumNumberType::J],
                   qns[QuantumNumberType::Lambda],
                   qns[QuantumNumberType::S], GS, GL);
    case Index(Hund::CaseB):
      return caseB(qns[QuantumNumberType::N], 
                   qns[QuantumNumberType::J],
                   qns[QuantumNumberType::Lambda],
                   qns[QuantumNumberType::S], GS, GL);
    default:
      throw std::runtime_error("cannot understand upper hund case");
  }
}


inline Numeric frequency_shift_per_teslaByHund(const QuantumNumbers& upper, const QuantumNumbers& lower, const Index species)
{
  // Find the constants
  const Numeric GS = get_lande_spin_constant(species);
  const Numeric GL = get_lande_lambda_constant();
  
  // Set the g*M factors
  const Numeric gMl = lower[QuantumNumberType::M].toNumeric() * gHund(lower, GS, GL);
  const Numeric gMu = upper[QuantumNumberType::M].toNumeric() * gHund(upper, GS, GL);
  
  // convert from energy state to frequency and be done with it
  return (gMl - gMu) * ZeemanSplittingConstant;
}


inline Numeric frequency_shift_per_teslaByGData(const QuantumNumbers& upper, 
                                                const QuantumNumbers& lower, 
                                                const Numeric& gu, 
                                                const Numeric& gl) noexcept
{
  // Set the g*M factors
  const Numeric gMu = gu * upper[QuantumNumberType::M].toNumeric();
  const Numeric gMl = gl * lower[QuantumNumberType::M].toNumeric();
  
  // convert from energy state to frequency and be done with it
  return (gMl - gMu) * ZeemanSplittingConstant;
}


inline bool hund_compatible(const QuantumNumbers& qns) noexcept
{
  switch(qns[QuantumNumberType::Hund].toIndex()) {
    case Index(Hund::CaseA):
      if(qns[QuantumNumberType::Omega].isUndefined() or
         qns[QuantumNumberType::J].isUndefined() or
         qns[QuantumNumberType::Lambda].isUndefined() or
         qns[QuantumNumberType::S].isUndefined())
        return false;
      break;
    case Index(Hund::CaseB):
      if(qns[QuantumNumberType::N].isUndefined() or
         qns[QuantumNumberType::J].isUndefined() or
         qns[QuantumNumberType::Lambda].isUndefined() or
         qns[QuantumNumberType::S].isUndefined())
        return false;
      break;
    default:
      return false;
  }
  
  return true;
}


inline bool J_compatible(const Rational& J1, const Rational& J2)
{
  if(abs((J1 - J2).toNumeric()) > 1.0)
    return false;
  return true;
}


inline ZeemanPolarizationType get_polarization(const Rational& Mu, const Rational& Ml)
{
  const Index dM = (Ml - Mu).toIndex();
  if(abs(dM) > 1)
    throw std::runtime_error("bad M-values. encountered too high change from upper to lower state");
  
  switch(dM) {
    case -1: return ZeemanPolarizationType::SigmaMinus;
    case  0: return ZeemanPolarizationType::Pi;
    case +1: return ZeemanPolarizationType::SigmaPlus;
    default: return ZeemanPolarizationType::None;
  } 
}


Numeric ZeemanEffectData::frequency_shift_per_tesla(const QuantumNumbers& upper, const QuantumNumbers& lower, const Index species) const
{
  switch(msplit) {
    case ZeemanSplittingType::None:
      return 0.0;
    case ZeemanSplittingType::ByHund:
      return frequency_shift_per_teslaByHund(upper, lower, species);
    case ZeemanSplittingType::ByGData:
      return frequency_shift_per_teslaByGData(upper, lower, mdata[Index(ByGDataPos::GU)], mdata[Index(ByGDataPos::GL)]);
    case ZeemanSplittingType::ByPrecalc:
      return mdata[Index(ByPrecalcPos::DF)];
    default:
      throw std::runtime_error("Unknown type for frequency shift");
  }
}


void ZeemanEffectData::convertNoneToHund(const QuantumNumbers& upper, const QuantumNumbers& lower)
{
  assert(msplit == ZeemanSplittingType::None);
  
  if(lower[QuantumNumberType::Hund].isUndefined() or upper[QuantumNumberType::Hund].isUndefined())
    throw std::runtime_error("Undefined Hund-cases encountered...");
  else if(not (hund_compatible(upper) and hund_compatible(lower)))
    throw std::runtime_error("line incompatible with Hund-style Zeeman calculations");
  else if(not J_compatible(lower[QuantumNumberType::J], upper[QuantumNumberType::J]))
    throw std::runtime_error("Too large difference in J-values of transition");
  
  msplit = ZeemanSplittingType::ByHund;
  mdata = Vector(Index(ByHundPos::LEN));
  
  assert(ok());
}


void ZeemanEffectData::setNumericalAndPolarization(const QuantumNumbers& upper, const QuantumNumbers& lower, const Index species)
{
  if(not (msplit == ZeemanSplittingType::ByPrecalc)) {
    mdata = Vector(Index(ByPrecalcPos::LEN), frequency_shift_per_tesla(upper, lower, species));
    msplit = ZeemanSplittingType::ByPrecalc;
  }
  
  mpolar = get_polarization(upper[QuantumNumberType::M], lower[QuantumNumberType::M]);
  
  assert(ok());
}

Index ZeemanEffectData::dM()
{
  switch(mpolar) {
    case ZeemanPolarizationType::SigmaMinus: return -1;
    case ZeemanPolarizationType::Pi:         return  0;
    case ZeemanPolarizationType::SigmaPlus:  return +1;
    case ZeemanPolarizationType::None: default: throw std::runtime_error("Unkown change in M");
  }
}


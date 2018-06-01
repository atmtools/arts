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
#include "wigner_functions.h"

Numeric ZeemanEffectData::SplittingConstant(const Index i) const
{
  extern const Numeric PLANCK_CONST;
  extern const Numeric BOHR_MAGNETON;
  static const Numeric C = BOHR_MAGNETON / PLANCK_CONST;
  
  return C * (mMl[i].toNumeric() * mgl - mMu[i].toNumeric() * mgu);
}


inline bool hund_compatible(const QuantumNumbers& qns) noexcept
{
  switch(Hund(qns[QuantumNumberType::Hund].toIndex())) {
    case Hund::CaseA:
      if(qns[QuantumNumberType::Omega].isUndefined() or
         qns[QuantumNumberType::J].isUndefined() or
         qns[QuantumNumberType::Lambda].isUndefined() or
         qns[QuantumNumberType::S].isUndefined())
        return false;
      break;
    case Hund::CaseB:
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
  if(JJ == 0)          g = 0.0;
  else if(NN not_eq 0) g = GS * T1.toNumeric() + GL * T2.toNumeric();
  else                 g = GS * T1.toNumeric();
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
  if(JJ == 0) g = 0.0;
  else        g = GS * T1.toNumeric() + GL * T2.toNumeric();
  return g; 
}


inline Numeric gHund(const QuantumNumbers& qns, const Numeric& GS, const Numeric& GL)
{
  if(not hund_compatible(qns))
    throw std::runtime_error("Bad quantum numbers for Zeeman via Hund approximation");
  
  switch(Hund(qns[QuantumNumberType::Hund].toIndex())) {
    case Hund::CaseA:
      return caseA(qns[QuantumNumberType::Omega], 
                   qns[QuantumNumberType::J],
                   qns[QuantumNumberType::Lambda],
                   qns[QuantumNumberType::S], GS, GL);
    case Hund::CaseB:
      return caseB(qns[QuantumNumberType::N], 
                   qns[QuantumNumberType::J],
                   qns[QuantumNumberType::Lambda],
                   qns[QuantumNumberType::S], GS, GL);
    default:
      throw std::runtime_error("cannot understand upper hund case");
  }
}

inline void are_Ju_and_Jl_compatible(const Rational& Ju, const Rational& Jl)
{
  Rational dJ = Jl - Ju;
  dJ.Simplify();
  if(dJ.Denom() not_eq 1)
    throw std::runtime_error("Delta J is not an index");
  else if(abs(dJ.toIndex()) > 1)
    throw std::runtime_error("Delta J is neither of -1, 0, nor 1");
  else if(Ju < 1)
    throw std::runtime_error("Upper J must be above 0");
}


ZeemanEffectData::ZeemanEffectData(const Numeric& gu, const Numeric& gl, const QuantumIdentifier& qi, const ZeemanPolarizationType polarization)
{
  mpolar = polarization;
  mgl = gl;
  mgu = gu;
  
  const Rational& Jl = qi.LowerQuantumNumbers()[QuantumNumberType::J];
  const Rational& Ju = qi.UpperQuantumNumbers()[QuantumNumberType::J];
  are_Ju_and_Jl_compatible(Ju, Jl);

  // Find M-vectors
  Rational end, start, dM;
  Numeric C=0;
  switch(mpolar) {
    case ZeemanPolarizationType::Pi:
      C = 1.5;
      dM = Rational(0, 1);
      end = min(Ju, Jl);
      start = -end;
      break;
    case ZeemanPolarizationType::SigmaPlus:
      C = .75;
      dM = Rational(1, 1);
      start = -Ju;
      if(Ju < Jl)       end = Ju + 1;
      else if(Ju == Jl) end = Ju;
      else              end = Jl;
      break;
    case ZeemanPolarizationType::SigmaMinus:
      C = .75;
      dM = Rational(-1, 1);
      end = Ju + 1;
      if(Ju < Jl)       start = -Ju;
      else if(Ju == Jl) start = -Ju + 1;
      else              start = -Ju + 2;
      break;
    case ZeemanPolarizationType::None:
      throw std::runtime_error("To developer: never initialize ZeemanEffectData without known polarization.");
      break;
  }
  
  // Set the computational data 
  mnelem = (end - start).toIndex() + 1;
  mS0.resize(mnelem);
  mMu.resize(mnelem);
  mMl.resize(mnelem);
  for(Index i=0; i<mnelem; i++) {
    mMu[i] = start + i;
    mMl[i] = mMu[i] + dM;
    mS0[i] = wigner3j(Jl, Rational(1), Ju, mMl[i], -dM, -mMu[i]);
    mS0[i] *= mS0[i] * C;
  }
}


ZeemanEffectData::ZeemanEffectData(const QuantumIdentifier& qi, const ZeemanPolarizationType polarization)
{
  const Numeric GS = get_lande_spin_constant(qi.Species());
  const Numeric GL = get_lande_lambda_constant();
  const Numeric gl = gHund(qi.LowerQuantumNumbers(), GS, GL);;
  const Numeric gu = gHund(qi.UpperQuantumNumbers(), GS, GL);;
  
//   ZeemanEffectData(gHund(qi.UpperQuantumNumbers(), GS, GL), gHund(qi.LowerQuantumNumbers(), GS, GL), qi, polarization);
  // All below here should be a perfect copy of the above initialization
  
  
  mpolar = polarization;
  mgl = gl;
  mgu = gu;
  
  const Rational& Jl = qi.LowerQuantumNumbers()[QuantumNumberType::J];
  const Rational& Ju = qi.UpperQuantumNumbers()[QuantumNumberType::J];
  are_Ju_and_Jl_compatible(Ju, Jl);
  
  // Find M-vectors
  Rational end, start, dM;
  Numeric C=0.75;
  switch(mpolar) {
    case ZeemanPolarizationType::Pi:
      C *= 2;
      dM = Rational(0, 1);
      end = min(Ju, Jl);
      start = -end;
      break;
    case ZeemanPolarizationType::SigmaPlus:
      dM = Rational(1, 1);
      start = -Ju;
      if(Ju < Jl)       end = Ju;
      else if(Ju == Jl) end = Ju-1;
      else              end = Jl-1;
      break;
    case ZeemanPolarizationType::SigmaMinus:
      dM = Rational(-1, 1);
      end = Ju;
      if(Ju < Jl)       start = -Ju;
      else if(Ju == Jl) start = -Ju + 1;
      else              start = -Ju + 2;
      break;
    case ZeemanPolarizationType::None:
      throw std::runtime_error("To developer: never initialize ZeemanEffectData without known polarization.  Throw warnings earlier than here!");
  }
  
  // Set the computational data 
  mnelem = (end - start).toIndex() + 1;
  mS0.resize(mnelem);
  mMu.resize(mnelem);
  mMl.resize(mnelem);
  for(Index i=0; i<mnelem; i++) {
    mMu[i] = start + i;
    mMl[i] = mMu[i] + dM;
    mS0[i] = wigner3j(Jl, Rational(1), Ju, mMl[i], -dM, -mMu[i]);
    mS0[i] *= mS0[i] * C;
  }
}

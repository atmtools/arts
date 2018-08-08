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

/** Contains Zeeman Effect data class
 * \file   zeemandata.h
 * 
 * \author Richard Larsson
 * \date   2018-04-06
 **/

#ifndef zeemandata_h
#define zeemandata_h

#include "matpackI.h"
#include "mystring.h"
#include "quantum.h"

enum class ZeemanPolarizationType : Index { SigmaMinus=-1, Pi=0, SigmaPlus=1, None };

class ZeemanEffectData 
{
public:
  
  ZeemanEffectData() : mpolar(ZeemanPolarizationType::None), mgu(0), mgl(0), mnelem(1), mMu(1, 0), mMl(1, 0), mS0(1, 1) {};
  ZeemanEffectData(const Numeric& gu, const Numeric& gl, const QuantumIdentifier& qi, const ZeemanPolarizationType polarization);
  ZeemanEffectData(const QuantumIdentifier& qi, const ZeemanPolarizationType polarization);
  ZeemanEffectData(const ZeemanEffectData& zed) : mpolar(zed.mpolar), mgu(zed.mgu), mgl(zed.mgl), mnelem(zed.mnelem), mMu(zed.mMu), mMl(zed.mMl), mS0(zed.mS0) {};
  ZeemanEffectData(ZeemanEffectData&& zed) : mpolar(std::move(zed.mpolar)), mgu(std::move(zed.mgu)), mgl(std::move(zed.mgl)), mnelem(std::move(zed.mnelem)), mMu(std::move(zed.mMu)), mMl(std::move(zed.mMl)), mS0(std::move(zed.mS0)) {};
  
  void init(const Numeric& gu, const Numeric& gl, const QuantumIdentifier& qi, const ZeemanPolarizationType polarization);
  
  ZeemanEffectData& operator=(ZeemanEffectData&& zed)
  {
    if (this not_eq &zed) {
      mpolar = std::move(zed.mpolar);
      mgu = std::move(zed.mgu);
      mgl = std::move(zed.mgl);
      mnelem = std::move(zed.mnelem);
      mMu = std::move(zed.mMu);
      mMl = std::move(zed.mMl);
      mS0 = std::move(zed.mS0);
    }
    return *this;
  }
  
  ZeemanEffectData& operator=(const ZeemanEffectData& zed)
  {
    mpolar = zed.mpolar;
    mgu = zed.mgu;
    mgl = zed.mgl;
    mnelem = zed.mnelem;
    mMu = zed.mMu;
    mMl = zed.mMl;
    mS0 = zed.mS0;
    return *this;
  }
  
  // Access to data for settings and iterations
  ZeemanPolarizationType PolarizationType() const {return mpolar;}
  const Numeric& UpperG() const { return mgu; }
  const Numeric& LowerG() const { return mgl; }
  Numeric& UpperG() { return mgu; }
  Numeric& LowerG() { return mgl; }
  Index nelem() const { return mnelem; }
  
  Numeric SplittingConstant(const Index i) const;
  Numeric StrengthScaling(const Index i) const {return mS0[i];}
  Numeric SumStrengthScale() const {return mS0.sum();}
  Rational Mu(const Index i) const {return mMu[i];}
  Rational Ml(const Index i) const {return mMl[i];}
  
  Vector Polarization(const Numeric& theta, const Numeric& eta) const;
  Vector dPolarization_dtheta(const Numeric& theta, const Numeric& eta) const;
  Vector dPolarization_deta(const Numeric& theta, const Numeric& eta) const;
  
  void SetPolarizationTypeFromString(const String& t)
  {
    if(t=="PI")      mpolar = ZeemanPolarizationType::Pi; 
    else if(t=="SM") mpolar = ZeemanPolarizationType::SigmaMinus; 
    else if(t=="SP") mpolar = ZeemanPolarizationType::SigmaPlus;
    else             mpolar = ZeemanPolarizationType::None;
  }
  
  String StringFromPolarizationType() const
  {
    switch(mpolar) { 
      case ZeemanPolarizationType::SigmaMinus: return "SM";
      case ZeemanPolarizationType::Pi:         return "PI";
      case ZeemanPolarizationType::SigmaPlus:  return "SP";
      case ZeemanPolarizationType::None: {}
    }
    return "-1";
  }
  
private:
  ZeemanPolarizationType mpolar;
  Numeric                mgu;
  Numeric                mgl;
  Index                  mnelem;
  std::vector<Rational>  mMu;
  std::vector<Rational>  mMl;
  Vector                 mS0;
};

#endif /* zeemandata_h */

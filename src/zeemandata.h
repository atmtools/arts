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

#include "wigner_functions.h"
#include "mystring.h"
#include "quantum.h"
#include "constants.h"

namespace Zeeman {
  enum class Polarization { SigmaMinus, Pi, SigmaPlus };
  
  inline String polarization2string(Polarization type) noexcept {
    switch(type) {
      case Polarization::SigmaMinus:  return "SM";
      case Polarization::Pi:          return "PI";
      case Polarization::SigmaPlus:   return "SP";
      default:                        return "#";
    };
  }
  
  inline Polarization string2polarization(const String& type) {
    if(type == "SM")      return Polarization::SigmaMinus;
    else if(type == "PI") return Polarization::Pi;
    else if(type == "SP") return Polarization::SigmaPlus;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  "
      << "See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }
  
  inline Index dM(Polarization type) {
    switch(type) {
      case Polarization::SigmaMinus:  return -1;
      case Polarization::Pi:          return  0;
      case Polarization::SigmaPlus:   return  1;
      default:                        return  0;
    }
  }
  
  inline Rational start(Rational Ju, Rational Jl, Polarization type) noexcept {
    switch(type) {
      case Polarization::SigmaMinus:  
        if(Ju < Jl)                   return -Ju;
        else if(Ju == Jl)             return -Ju + 1;
        else                          return -Ju + 2;
      case Polarization::Pi:          return -std::min(Ju, Jl);
      case Polarization::SigmaPlus:   return -Ju;
      default:                        return 0;
    }
  }
  
  inline Rational end(Rational Ju, Rational Jl, Polarization type) noexcept {
    switch(type) {
      case Polarization::SigmaMinus:  return Ju + 1;
      case Polarization::Pi:          return std::min(Ju, Jl);
      case Polarization::SigmaPlus:
        if(Ju < Jl)                   return Ju + 1;
        else if(Ju == Jl)             return Ju;
        else                          return Jl;
      default:                        return 0;
    }
  }
  
  inline Index nelem(Rational Ju, Rational Jl, Polarization type) noexcept {
    return (end(Ju, Jl, type) - start(Ju, Jl, type)).toIndex() + 1;
  }
  
  inline Rational Mu(Rational Ju, Rational Jl, Polarization type, Index n) noexcept {
    return start(Ju, Jl, type) + n;
  }
  
  inline Rational Ml(Rational Ju, Rational Jl, Polarization type, Index n) noexcept {
    return Mu(Ju, Jl, type, n) + dM(type);
  }
  
  inline Numeric PolarizationFactor(Polarization type) noexcept {
    switch(type) {
      case Polarization::SigmaMinus:  return .75;
      case Polarization::Pi:          return 1.5;
      case Polarization::SigmaPlus:   return .75;
      default:                        return 1.0;
    }
  }
  
  inline bool GoodHundData(const QuantumNumbers& qns) noexcept {
    if(qns[QuantumNumberType::Hund].isUndefined())
      return false;
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
  
  inline Numeric SimpleGCaseB(Rational N,
                              Rational J,
                              Rational Lambda,
                              Rational S,
                              Numeric GS,
                              Numeric GL) noexcept {
    auto JJ = J*(J+1);
    auto NN = N*(N+1);
    auto SS = S*(S+1);
    auto LL = Lambda*Lambda;
    
    if(JJ == 0)
      return 0.0;
    else if(NN not_eq 0) {
      auto T1 = ((JJ + SS - NN)           / JJ / 2).toNumeric();
      auto T2 = ((JJ - SS + NN) * LL / NN / JJ / 2).toNumeric();
      return GS * T1 + GL * T2;
    }
    else {
      auto T1 = ((JJ + SS - NN)           / JJ / 2).toNumeric();
      return GS * T1;
    }
  }

  inline Numeric SimpleGCaseA(Rational Omega,
                              Rational J,
                              Rational Lambda,
                              Rational Sigma,
                              Numeric GS,
                              Numeric GL) noexcept {  
    auto JJ = J*(J+1);
    
    if(JJ == 0)
      return 0.0;
    else {
      auto DIV = Omega / JJ;
      auto T1 = (Sigma  * DIV).toNumeric();
      auto T2 = (Lambda * DIV).toNumeric();
      return GS * T1 + GL * T2;
    }
  }
  
  inline Numeric SimpleG(const QuantumNumbers& qns, const Numeric& GS, const Numeric& GL) {
    if(not GoodHundData(qns)) {
      std::ostringstream os;
      os << "Bad quantum numbers for Zeeman via Hund approximation:\n" << qns;
      throw std::runtime_error(os.str());
    }
    
    switch(Hund(qns[QuantumNumberType::Hund].toIndex())) {
      case Hund::CaseA:
        return SimpleGCaseA(qns[QuantumNumberType::Omega], 
                            qns[QuantumNumberType::J],
                            qns[QuantumNumberType::Lambda],
                            qns[QuantumNumberType::S], GS, GL);
      case Hund::CaseB:
        return SimpleGCaseB(qns[QuantumNumberType::N], 
                            qns[QuantumNumberType::J],
                            qns[QuantumNumberType::Lambda],
                            qns[QuantumNumberType::S], GS, GL);
      default:
        throw std::runtime_error("cannot understand hund case");
    }
  }
  
  struct SplittingData {Numeric gu, gl;};
  
  class Model {
  private:
    SplittingData mdata;
    
  public:
    constexpr Model(SplittingData gs={0, 0}) noexcept : mdata(gs) {}
    Model(const QuantumIdentifier& qid);
    
    constexpr bool empty() const noexcept {return mdata.gu==mdata.gl and mdata.gu==0;}
    
    Numeric& gu() noexcept {return mdata.gu;}
    Numeric& gl() noexcept {return mdata.gl;}
    constexpr Numeric gu() const noexcept {return mdata.gu;}
    constexpr Numeric gl() const noexcept {return mdata.gl;}
    
    Numeric Strength(Rational Ju, Rational Jl, Polarization type, Index n) const {
      using Constant::pow2;
      
      auto ml = Ml(Ju, Jl, type, n);
      auto mu = Mu(Ju, Jl, type, n);
      auto dm = dM(type);
      return PolarizationFactor(type) * pow2(wigner3j(Jl,  1, Ju,
                                                      ml,-dm,-mu));
    }
    
    Numeric Splitting(Rational Ju, Rational Jl, Polarization type, Index n) const noexcept {
      using Constant::h;
      using Constant::bohr_magneton;
      constexpr Numeric C = bohr_magneton / h;
      
      return C * (Ml(Ju, Jl, type, n).toNumeric() * gl() - 
                  Mu(Ju, Jl, type, n).toNumeric() * gu());
    }
    
    friend inline std::ostream& operator<<(std::ostream& os, const Model& m);
    friend inline std::istream& operator>>(std::istream& is, Model& m);
  };
  
  Model GetSimpleModel(const QuantumIdentifier& qid);
  Model GetAdvancedModel(const QuantumIdentifier& qid);
  
  inline std::ostream& operator<<(std::ostream& os, const Model& m) {
    os << m.mdata.gu << ' ' << m.mdata.gl;
    return os;
  }
  
  inline std::istream& operator>>(std::istream& is, Model& m) {
    is >> m.mdata.gu >> m.mdata.gl;
    return is;
  }
  
  class PolarizationVector {
  private:
    Eigen::RowVector4d att;
    Eigen::RowVector3d dis;
    
  public:
    PolarizationVector(Numeric a=1, Numeric b=0, Numeric c=0, Numeric d=0, Numeric u=0, Numeric v=0, Numeric w=0) noexcept :
    att(a, b, c, d), dis(u, v, w) {};
    
    const Eigen::RowVector4d& attenuation() const noexcept {return att;}
    const Eigen::RowVector3d& dispersion() const noexcept {return dis;}
    Eigen::RowVector4d& attenuation() noexcept {return att;}
    Eigen::RowVector3d& dispersion() noexcept {return dis;}
    
    Eigen::Matrix4d matrix() const noexcept {
      return (Eigen::Matrix4d() << att[0], att[1], att[2], att[3],
                                   att[1], att[0], dis[0], dis[1],
                                   att[2],-dis[0], att[0], dis[2],
                                   att[3],-dis[1],-dis[2], att[0]).finished();
    }
  };
  
  struct AllPolarizationVectors {PolarizationVector sm, pi, sp;};
  
  inline AllPolarizationVectors AllPolarization(Numeric theta, Numeric eta) noexcept {
    const Numeric ST=std::sin(theta), CT=std::cos(theta), ST2=ST*ST, CT2=CT*CT, C2E=std::cos(2*eta), S2E=std::sin(2*eta), ST2C2E=ST2*C2E, ST2S2E=ST2*S2E;
    
    AllPolarizationVectors pv;
    pv.sm = PolarizationVector(1 + CT2,  ST2C2E,  ST2S2E,  2*CT,  4*CT,  2*ST2S2E, -2*ST2C2E);
    pv.pi = PolarizationVector(    ST2, -ST2C2E, -ST2S2E,     0,     0, -2*ST2S2E,  2*ST2C2E);
    pv.sp = PolarizationVector(1 + CT2,  ST2C2E,  ST2S2E, -2*CT, -4*CT,  2*ST2S2E, -2*ST2C2E);
    return pv;
  }
  
  inline AllPolarizationVectors AllPolarization_dtheta(Numeric theta, const Numeric eta) noexcept {
    const Numeric ST=std::sin(theta), CT=std::cos(theta), C2E=std::cos(2*eta), S2E=std::sin(2*eta),
    dST=CT, dST2=2*ST*dST, dCT=-ST, dST2C2E=dST2*C2E, dST2S2E=dST2*S2E, dCT2=2*CT*dCT;
    
    AllPolarizationVectors pv;
    pv.sm = PolarizationVector(dCT2,  dST2C2E,  dST2S2E,  2*dCT,  4*dCT,  2*dST2S2E, -2*dST2C2E);
    pv.pi = PolarizationVector(dST2, -dST2C2E, -dST2S2E,      0,      0, -2*dST2S2E,  2*dST2C2E);
    pv.sp = PolarizationVector(dCT2,  dST2C2E,  dST2S2E, -2*dCT, -4*dCT,  2*dST2S2E, -2*dST2C2E);
    return pv;
  }
  
  inline AllPolarizationVectors AllPolarization_deta(Numeric theta, Numeric eta) noexcept {
    const Numeric ST=std::sin(theta), ST2=ST*ST, C2E=std::cos(2*eta), S2E=std::sin(2*eta), dST2C2E=-2*ST2*S2E, dST2S2E=2*ST2*C2E;
    
    AllPolarizationVectors pv;
    pv.sm = PolarizationVector(0,  dST2C2E,  dST2S2E, 0, 0,  2*dST2S2E, -2*dST2C2E);
    pv.pi = PolarizationVector(0, -dST2C2E, -dST2S2E, 0, 0, -2*dST2S2E,  2*dST2C2E);
    pv.sp = PolarizationVector(0,  dST2C2E,  dST2S2E, 0, 0,  2*dST2S2E, -2*dST2C2E);
    return pv;
  }
  
  inline const PolarizationVector& SelectPolarization(const AllPolarizationVectors& data, Polarization mpolar) noexcept {
    switch(mpolar) {
      case Polarization::SigmaMinus:
        return data.sm;
      case Polarization::Pi:
        return data.pi;
      case Polarization::SigmaPlus:
        return data.sp;
    }
    std::terminate();
  }
};

#endif /* zeemandata_h */

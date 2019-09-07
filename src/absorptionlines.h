/* Copyright (C) 2019
 Richard Larsson <larsson@mps.mpg.de>

 
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

/** Contains the absorption namespace
 * @file   lineshapemodel.h
 * @author Richard Larsson
 * @date   2019-09-07
 * 
 * @brief  Contains the absorption lines implementation
 * 
 * This namespace contains classes to deal with absorption lines
 **/

#ifndef absorptionlines_h
#define absorptionlines_h

#include <vector>
#include "lineshapemodel.h"
#include "matpack.h"
#include "quantum.h"
#include "zeemandata.h"

/** Namespace to contain things required for absorption calculations */
namespace Absorption {
/** Describes the type of mirroring line effects
 * 
 * Each type but None has to have an implemented effect
 */
enum class MirroringType {
  None,             // No mirroring
  Lorentz,          // Mirror, but use Lorentz line shape
  SameAsLineShape,  // Mirror using the same line shape
  Manual,           // Mirror by having a line in the array of line record with negative F0
};  // MirroringType

/** Describes the type of normalization line effects
 *
 * Each type but None has to have an implemented effect
 */
enum class NormalizationType {
  None,  // Do not renormalize the line shape
  VVH,   // Renormalize with Van Vleck and Huber specifications
  VVW,   // Renormalize with Van Vleck and Weiskopf specifications
  RosenkranzQuadratic,  // Renormalize using Rosenkranz's quadratic specifications
};  // LineNormalizationType

/** Describes the type of population level counter
 *
 * The types here might require that different data is available at runtime absorption calculations
 */
enum class PopulationType {
  ByLTE,                      // Assume line is in LTE
  ByVibrationalTemperatures,  // Assume line is in NLTE described by vibrational temperatures
  ByPopulationDistribution,   // Assume line is in NLTE and the upper-to-lower ratio is known
};  // LinePopulationType

/** Describes the type of cutoff calculations */
enum class CutoffType {
  None,                // No cutoff frequency at all
  LineByLineOffset,    // The cutoff frequency is at SingleLine::F0 plus the cutoff frequency
  BandFixedFrequency,  // The curoff frequency is the cutoff frequency for all SingleLine(s)
};  // LineCutoffType

/** Computations and data for a single absorption line */
class SingleLine {
private:
  /** Central frequency */
  Numeric mF0;
  
  /** Reference intensity */
  Numeric mI0;
  
  /** Lower state energy level */
  Numeric mE0;
  
  /** Lower level statistical weight */
  Numeric mglow;
  
  /** Upper level statistical weight */
  Numeric mgupp;
  
  /** Einstein spontaneous emission coefficient */
  Numeric ma;
  
  /** Zeeman model */
  Zeeman::Model mzeeman;
  
  /** Line shape model */
  LineShape::Model2 mlineshape;
  
  /** Lower level quantum numbers */
  std::vector<Rational> mlowerquanta;
  
  /** Upper level quantum numbers */
  std::vector<Rational> mupperquanta;

public:
  /** Number of lineshape elements */
  Index LineShapeElems() const noexcept {return mlineshape.nelem();}
  
  /** Number of lower quantum numbers */
  Index LowerQuantumCount() const noexcept {return mlowerquanta.size();}
  
  /** Lower quantum number */
  Rational LowerQuantumNumber(size_t i) const noexcept {return mlowerquanta[i];}
  
  /** Number of upper quantum numbers */
  Index UpperQuantumCount() const noexcept {return mupperquanta.size();}
  
  /** Upper quantum number */
  Rational UpperQuantumNumber(size_t i) const noexcept {return mupperquanta[i];}
  
  /** Central frequency */
  Numeric F0() const noexcept {return mF0;}
  
  /** Lower level energy */
  Numeric E0() const noexcept {return mE0;}
  
  /** Reference line strength */
  Numeric I0() const noexcept {return mI0;}
  
  /** Einstein spontaneous emission */
  Numeric A() const noexcept {return ma;}
  
  /** Lower level statistical weight */
  Numeric g_low() const noexcept {return mglow;}
  
  /** Upper level statistical weight */
  Numeric g_upp() const noexcept {return mgupp;}
  
  /** Zeeman model */
  Zeeman::Model Zeeman() const noexcept {return mzeeman;}
  
  /** Line shape model */
  const LineShape::Model2& LineShape() const noexcept {return mlineshape;}
};  // SingleLine

class Lines {
private:
  /** Catalog ID */
  QuantumIdentifier mquantumidentity;
  
  /** Reference temperature for all parameters of the lines */
  Numeric mT0;
  
  /** cutoff frequency */
  Numeric mcutofffreq;
  
  /** linemixing limit */
  Numeric mlinemixinglimit;
  
  /** cutoff type, by band or by line */
  CutoffType mcutoff;
  
  /** Mirroring type */
  MirroringType mmirroring;
  
  /** Line population distribution */
  PopulationType mpopulation;
  
  /** Line normalization type */
  NormalizationType mnormalization;
  
  /** List of local quantum numbers, these must be defined */
  std::vector<QuantumNumberType> mlocalquanta;

  /** Type of line shape */
  LineShape::Type mlineshapetype;
  
  /** Does the line broadening have self broadening */
  bool mselfbroadening;
  
  /** Does the line broadening have bath broadening */
  bool mbathbroadening;
  
  /** A list of broadening species */
  std::vector<SpeciesTag> mbroadeningspecies;
  
  /** A list of individual lines */
  std::vector<SingleLine> mlines;
  
public:
  /** Species Index */
  Index Species() const noexcept {return mquantumidentity.Species();}
  
  /** Isotopologue Index */
  Index Isotopologue() const noexcept {return mquantumidentity.Isotopologue();}
  
  /** Quantum number lower level */
  Rational LowerQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept {
    for(size_t i=0; i<mlocalquanta.size(); i++)
      if(mlocalquanta[i] == qnt)
        return mlines[k].LowerQuantumNumber(i);
    return mquantumidentity.LowerQuantumNumber(qnt);
  }
  
  /** Quantum number upper level */
  Rational UpperQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept {
    for(size_t i=0; i<mlocalquanta.size(); i++)
      if(mlocalquanta[i] == qnt)
        return mlines[k].UpperQuantumNumber(i);
    return mquantumidentity.UpperQuantumNumber(qnt);
  }
  
  /** Checks if all defined quantum numbers in qid are equal to the lower levels of specified line */
  bool InLowerLevel(size_t k, const QuantumIdentifier& qid) const noexcept {
    if(mquantumidentity.Species() not_eq qid.Species())
      return false;
    if(mquantumidentity.Isotopologue() not_eq qid.Isotopologue())
      return false;
    if(qid.Type() == QuantumIdentifier::ALL)
      return true;
    else if(qid.Type() == QuantumIdentifier::NONE)
      return false;
    else if(qid.Type() == QuantumIdentifier::ENERGY_LEVEL) {
      for(size_t i=0; i<mlocalquanta.size(); i++) {
        const auto qn = qid.EnergyLevelQuantumNumbers()[mlocalquanta[i]];
        if(qn.isDefined() and qn not_eq mlines[k].LowerQuantumNumber(i)) {
          return false;
        }
      }
    }
    else if(qid.Type() == QuantumIdentifier::TRANSITION)
      return false;
    
    return qid.InLower(mquantumidentity.LowerQuantumId());
  }
  
  /** Checks if all defined quantum numbers in qid are equal to the upper levels of specified line */
  bool InUpperLevel(size_t k, const QuantumIdentifier& qid) const noexcept {
    if(mquantumidentity.Species() not_eq qid.Species())
      return false;
    if(mquantumidentity.Isotopologue() not_eq qid.Isotopologue())
      return false;
    
    if(qid.Type() == QuantumIdentifier::ALL)
      return true;
    else if(qid.Type() == QuantumIdentifier::NONE)
      return false;
    else if(qid.Type() == QuantumIdentifier::ENERGY_LEVEL) {
      for(size_t i=0; i<mlocalquanta.size(); i++) {
        const auto qn = qid.EnergyLevelQuantumNumbers()[mlocalquanta[i]];
        if(qn.isDefined() and qn not_eq mlines[k].UpperQuantumNumber(i)) {
          return false;
        }
      }
    }
    else if(qid.Type() == QuantumIdentifier::TRANSITION)
      return false;
    
    return qid.InUpper(mquantumidentity.UpperQuantumId());
  }
  
  /** Returns the number of Zeeman split lines */
  Index ZeemanCount(size_t k, Zeeman::Polarization type) const noexcept {
    return Zeeman::nelem(UpperQuantumNumber(k, QuantumNumberType::J),
                         LowerQuantumNumber(k, QuantumNumberType::J),
                         type);
  }
  
  /** Returns the strength of a Zeeman split line */
  Numeric ZeemanStrength(size_t k, Zeeman::Polarization type, Index i) const noexcept {
    return mlines[k].Zeeman().Strength(UpperQuantumNumber(k, QuantumNumberType::J),
                                       LowerQuantumNumber(k, QuantumNumberType::J),
                                       type, i);
  }
  
  /** Returns the splitting of a Zeeman split line */
  Numeric ZeemanSplitting(size_t k, Zeeman::Polarization type, Index i) const noexcept {
    return mlines[k].Zeeman().Splitting(UpperQuantumNumber(k, QuantumNumberType::J),
                                        LowerQuantumNumber(k, QuantumNumberType::J),
                                        type, i);
  }
  
  /** Central frequency */
  Numeric F0(size_t k) const noexcept {return mlines[k].F0();}
  
  /** Lower level energy */
  Numeric E0(size_t k) const noexcept {return mlines[k].E0();}
  
  /** Reference line strength */
  Numeric I0(size_t k) const noexcept {return mlines[k].I0();}
  
  /** Einstein spontaneous emission */
  Numeric A(size_t k) const noexcept {return mlines[k].A();}
  
  /** Lower level statistical weight */
  Numeric g_low(size_t k) const noexcept {return mlines[k].g_low();}
  
  /** Upper level statistical weight */
  Numeric g_upp(size_t k) const noexcept {return mlines[k].g_upp();}
  
  /** Returns mirroring style */
  MirroringType Mirroring() const noexcept {return mmirroring;}
  
  /** Returns normalization style */
  NormalizationType Normalization() const noexcept {return mnormalization;}
  
  /** Returns if the pressure should do line mixing */
  bool DoLineMixing(Numeric P) const noexcept {return mlinemixinglimit > P;}
  
  /** Returns line shape parameters */
  LineShape::Output ShapeParameters(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept {
    auto x = mlines[k].LineShape().GetParams(T, mT0, P, vmrs);
    if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
    return x;
  }
  
  /** Returns line shape parameters */
  LineShape::Output ShapeParameters_dT(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept {
    auto x = mlines[k].LineShape().GetTemperatureDerivs(T, mT0, P, vmrs);
    if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
    return x;
  }
  
  Index LineShapePos(const Index& spec) const noexcept {
    for(size_t i=size_t(mselfbroadening); i<mbroadeningspecies.size()-size_t(mbathbroadening); i++)
      if(spec == mbroadeningspecies[i].Species())
        return Index(i);
    return -1;
  }
  
  Index LineShapePos(const QuantumIdentifier& qi) const noexcept {
    LineShapePos(qi.Species());
  }
  
  /** Returns line shape parameters */
  LineShape::Output ShapeParameters_dVMR(size_t k, Numeric T, Numeric P, 
                                         const QuantumIdentifier& vmr_qi) const noexcept {
    const bool self = vmr_qi.Species() == mquantumidentity.Species();
    const auto& ls = mlines[k].LineShape();
    if (mselfbroadening and self) {
      auto x = ls.GetVMRDerivs(T, mT0, P, 0);
      
      if (mbathbroadening)
        x = LineShape::differenceOutput(x, ls.GetVMRDerivs(
            T, mT0, P, ls.nelem() - 1));
      
      if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
      return x;
    } else if (mbathbroadening and self)
      return {0, 0, 0, 0, 0, 0, 0, 0, 0};
    else {
      auto x = ls.GetVMRDerivs(T, mT0, P, LineShapePos(vmr_qi));
      
      if (mbathbroadening)
        x = LineShape::differenceOutput(x, ls.GetVMRDerivs(
            T, mT0, P, ls.nelem() - 1));
      
      if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
      return x;
    }
  }
  
  /***/
  Numeric ShapeParameterDerivative(size_t k, Numeric T, Numeric P, 
                                   const Vector& vmrs,
                                   const RetrievalQuantity& derivative) const noexcept {

    const auto self = derivative.Mode() == LineShape::self_broadening;
    const auto bath = derivative.Mode() == LineShape::bath_broadening;
    const auto& ls = mlines[k].LineShape();
    
    if(derivative.QuantumIdentity().Species() != Species() or
       derivative.QuantumIdentity().Isotopologue() != Isotopologue()))
      return 0;
    else if(self and mselfbroadening)
      return ls.GetInternalDeriv(
        T, mT0, P, 0, vmrs, derivative.PropMatType());
    else if(self)
      return ls.GetInternalDeriv(
        T, mT0, P, LineShapePos(SpeciesTag(derivative.Mode()).Species()), vmrs, derivative.PropMatType());
    else if(bath and mbathbroadening)
      return ls.GetInternalDeriv(
        T, mT0, P, ls.nelem() - 1, vmrs, derivative.PropMatType());
    else if(bath)
      return 0;
    else
      return ls.GetInternalDeriv(
        T, mT0, P, LineShapePos(SpeciesTag(derivative.Mode()).Species()), vmrs, derivative.PropMatType());
  }
  
  /** Returns cutoff frequency */
  Numeric Cutoff(size_t k) const noexcept {
    switch(mcutoff) {
      case CutoffType::LineByLineOffset:
        return F0(k) + mcutofffreq;
      case CutoffType::BandFixedFrequency:
        return mcutofffreq;
      case CutoffType::None:
        return -1;
    }
    return -1;
  }
};  // Lines
};  // Absorption

typedef Absorption::Lines AbsorptionLines;
typedef Array<AbsorptionLines> ArrayOfAbsorptionLines;
typedef Array<ArrayOfAbsorptionLines> ArrayOfArrayOfAbsorptionLines;

#endif  // absorptionlines_h

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
 * @file   lineshapemodel.cc
 * @author Richard Larsson
 * @date   2019-09-07
 * 
 * @brief  Contains the absorption lines implementation
 * 
 * This namespace contains classes to deal with absorption lines
 **/

#include "absorptionlines.h"

Rational Absorption::Lines::LowerQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept {
  for(size_t i=0; i<mlocalquanta.size(); i++)
    if(mlocalquanta[i] == qnt)
      return mlines[k].LowerQuantumNumber(i);
  
  return mquantumidentity.LowerQuantumNumber(qnt);
}

Rational Absorption::Lines::UpperQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept {
  for(size_t i=0; i<mlocalquanta.size(); i++)
    if(mlocalquanta[i] == qnt)
      return mlines[k].UpperQuantumNumber(i);
  
  return mquantumidentity.UpperQuantumNumber(qnt);
}

bool Absorption::Lines::InLowerLevel(size_t k,
                                     const QuantumIdentifier& qid) const noexcept {
  if(mquantumidentity.Species() not_eq qid.Species())
    return false;
  else if(mquantumidentity.Isotopologue() not_eq qid.Isotopologue())
    return false;
  else if(qid.Type() == QuantumIdentifier::ALL)
    return true;
  else if(qid.Type() == QuantumIdentifier::NONE)
    return false;
  else if(qid.Type() == QuantumIdentifier::TRANSITION)
    return false;
  else if(qid.Type() == QuantumIdentifier::ENERGY_LEVEL)
    for(size_t i=0; i<mlocalquanta.size(); i++)
      if(qid.EnergyLevelQuantumNumber(mlocalquanta[i]).isDefined() and
         qid.EnergyLevelQuantumNumber(mlocalquanta[i]) not_eq mlines[k].LowerQuantumNumber(i))
        return false;
  
  return qid.InLower(mquantumidentity.LowerQuantumId());
}

bool Absorption::Lines::InUpperLevel(size_t k,
                                     const QuantumIdentifier& qid) const noexcept {
  if(mquantumidentity.Species() not_eq qid.Species())
    return false;
  else if(mquantumidentity.Isotopologue() not_eq qid.Isotopologue())
    return false;
  else if(qid.Type() == QuantumIdentifier::ALL)
    return true;
  else if(qid.Type() == QuantumIdentifier::NONE)
    return false;
  else if(qid.Type() == QuantumIdentifier::TRANSITION)
    return false;
  else if(qid.Type() == QuantumIdentifier::ENERGY_LEVEL)
    for(size_t i=0; i<mlocalquanta.size(); i++)
      if(qid.EnergyLevelQuantumNumber(mlocalquanta[i]).isDefined() and
         qid.EnergyLevelQuantumNumber(mlocalquanta[i]) not_eq mlines[k].UpperQuantumNumber(i))
        return false;
  
  return qid.InUpper(mquantumidentity.UpperQuantumId());
}

LineShape::Output Absorption::Lines::ShapeParameters(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept {
  auto x = mlines[k].LineShape().GetParams(T, mT0, P, vmrs);
  
  if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
  
  return x;
}

LineShape::Output Absorption::Lines::ShapeParameters_dT(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept {
  auto x = mlines[k].LineShape().GetTemperatureDerivs(T, mT0, P, vmrs);
  
  if (not DoLineMixing(P)) x.Y = x.G = x.DV = 0;
  
  return x;
}

Index Absorption::Lines::LineShapePos(const Index& spec) const noexcept {
  if(mselfbroadening and spec == mquantumidentity.Species())
    return 0;
  
  for(Index i=Index(mselfbroadening); i<Index(mbroadeningspecies.size())-Index(mbathbroadening); i++)
    if(spec == mbroadeningspecies[i].Species())
      return Index(i);
  
  return mbathbroadening ? Index(mbroadeningspecies.size())-Index(mbathbroadening) : -1;
}

LineShape::Output Absorption::Lines::ShapeParameters_dVMR(size_t k, Numeric T, Numeric P, const QuantumIdentifier& vmr_qi) const noexcept {
  const auto self = vmr_qi.Species() == mquantumidentity.Species();
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

Numeric Absorption::Lines::ShapeParameter_dInternal(size_t k, Numeric T, Numeric P, const Vector& vmrs, const RetrievalQuantity& derivative) const noexcept {
  const auto self = derivative.Mode() == LineShape::self_broadening;
  const auto bath = derivative.Mode() == LineShape::bath_broadening;
  const auto& ls = mlines[k].LineShape();
  
  if(derivative.QuantumIdentity().Species() not_eq Species() or
      derivative.QuantumIdentity().Isotopologue() not_eq Isotopologue())
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

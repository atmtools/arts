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

#include "absorption.h"
#include "global_data.h"
#include "quantum_parser_hitran.h"

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

bool Absorption::Lines::InLowerGlobalLevel(const QuantumIdentifier& qid) const noexcept {
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
    return qid.InLower(mquantumidentity.LowerQuantumId());
  else
    return false;
}

bool Absorption::Lines::InUpperGlobalLevel(const QuantumIdentifier& qid) const noexcept {
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
    return qid.InLower(mquantumidentity.LowerQuantumId());
  else
    return false;
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
  // Is always first if this is self and self broadening exists
  if(mselfbroadening and spec == mquantumidentity.Species())
    return 0;
  
  // First and last might be artificial so they should not be checked
  for(Index i=Index(mselfbroadening); i<Index(mbroadeningspecies.size())-Index(mbathbroadening); i++) {
    if(spec == mbroadeningspecies[i].Species())
      return Index(i);
  }
  
  // At this point, the ID is not explicitly among the broadeners, but bath broadening means its VMR still might matter
  if(mbathbroadening)
    return Index(mbroadeningspecies.size())-Index(mbathbroadening);
  else
    return -1;
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

Absorption::SingleLineExternalRead Absorption::ReadFromArtscat3Stream(istream& is) {
  
  // Default data and values for this type
  SingleLineExternalRead data;
  data.selfbroadening = true;
  data.bathbroadening = true;
  data.lineshapetype = LineShape::Type::VP;
  
  // Global species lookup data:
  using global_data::species_data;

  // We need a species index sorted by Arts identifier. Keep this in a
  // static variable, so that we have to do this only once.  The ARTS
  // species index is ArtsMap[<Arts String>].
  static map<String, SpecIsoMap> ArtsMap;

  // Remember if this stuff has already been initialized:
  static bool hinit = false;

  if (!hinit) {
    for (Index i = 0; i < species_data.nelem(); ++i) {
      const SpeciesRecord& sr = species_data[i];

      for (Index j = 0; j < sr.Isotopologue().nelem(); ++j) {
        SpecIsoMap indicies(i, j);
        String buf = sr.Name() + "-" + sr.Isotopologue()[j].Name();
        ArtsMap[buf] = indicies;
      }
    }
    hinit = true;
  }

  // This always contains the rest of the line to parse. At the
  // beginning the entire line. Line gets shorter and shorter as we
  // continue to extract stuff from the beginning.
  String line;

  // Look for more comments?
  bool comment = true;

  while (comment) {
    // Return true if eof is reached:
    if (is.eof()) return data;

    // Throw runtime_error if stream is bad:
    if (!is) throw runtime_error("Stream bad.");

    // Read line from file into linebuffer:
    getline(is, line);

    // It is possible that we were exactly at the end of the file before
    // calling getline. In that case the previous eof() was still false
    // because eof() evaluates only to true if one tries to read after the
    // end of the file. The following check catches this.
    if (line.nelem() == 0 && is.eof()) return data;

    // @ as first character marks catalogue entry
    char c;
    extract(c, line, 1);

    // check for empty line
    if (c == '@') {
      comment = false;
    }
  }

  // read the arts identifier String
  istringstream icecream(line);

  String artsid;
  icecream >> artsid;

  if (artsid.length() != 0) {
    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(artsid);
    if (i == ArtsMap.end()) {
      ostringstream os;
      os << "ARTS Tag: " << artsid << " is unknown.";
      throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;

    // Set mspecies:
    data.quantumidentity.SetSpecies(id.Speciesindex());

    // Set misotopologue:
    data.quantumidentity.SetIsotopologue(id.Isotopologueindex());

    // Extract center frequency:
    icecream >> data.line.F0();

    Numeric psf;
    // Extract pressure shift:
    icecream >> psf;

    // Extract intensity:
    icecream >> data.line.I0();

    // Extract reference temperature for Intensity in K:
    icecream >> data.T0;

    // Extract lower state energy:
    icecream >> data.line.E0();

    // Extract air broadening parameters:
    Numeric agam, sgam;
    icecream >> agam;
    icecream >> sgam;

    // Extract temperature coefficient of broadening parameters:
    Numeric nair, nself;
    icecream >> nair;
    icecream >> nself;

    // Extract reference temperature for broadening parameter in K:
    Numeric tgam;
    icecream >> tgam;

    // Extract the aux parameters:
    Index naux;
    icecream >> naux;

    // resize the aux array and read it
    ArrayOfNumeric maux;
    maux.resize(naux);

    for (Index j = 0; j < naux; j++) {
      icecream >> maux[j];
      //cout << "maux" << j << " = " << maux[j] << "\n";
    }

    // Extract accuracies:
    Numeric dagam, dsgam, dnair, dnself, dpsf;
    try {
      Numeric mdf, mdi0;
      icecream >> mdf;
      icecream >> mdi0;
      icecream >> dagam;
      icecream >> dsgam;
      icecream >> dnair;
      icecream >> dnself;
      icecream >> dpsf;
    } catch (const std::runtime_error&) {
      // Nothing to do here, the accuracies are optional, so we
      // just set them to -1 and continue reading the next line of
      // the catalogue
      dagam = -1;
      dsgam = -1;
      dnair = -1;
      dnself = -1;
      dpsf = -1;
    }

    // Fix if tgam is different from ti0
    if (tgam != data.T0) {
      agam = agam * pow(tgam / data.T0, nair);
      sgam = sgam * pow(tgam / data.T0, nself);
      psf = psf * pow(tgam / data.T0, (Numeric).25 + (Numeric)1.5 * nair);
    }

    // Set line shape computer
    LineShape::Model lineshapemodel(sgam, nself, agam, nair, psf);
    
    data.line.LineShape() = LineShape::Model2(lineshapemodel.Data());
  }

  // That's it!
  data.bad = false;
  return data;
}

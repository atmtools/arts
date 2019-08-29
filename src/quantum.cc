/* Copyright (C) 2013
   Oliver Lemke <olemke@core-dump.info>

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

/** \file
    Classes to handle quantum numbers.

    \author Oliver Lemke
*/

#include <stdexcept>
#include "quantum.h"
#include "global_data.h"
#include "absorption.h"


bool QuantumNumbers::Compare(const QuantumNumbers& qn) const
{
    const QuantumContainer& qnumbers2 = qn.GetNumbers();

    bool match = true;

    Index qnri = 0;

    // Compare all quantum numbers in mqnumbers and qnumbers2
    while (match && qnri != Index(QuantumNumberType::FINAL_ENTRY))
    {
        // If one of the two numbers is undefined, it is considered as
        // a match.
        if (!mqnumbers[qnri].isUndefined()
            && !qnumbers2[qnri].isUndefined()
            && mqnumbers[qnri] != qnumbers2[qnri])
            match = false;

        qnri++;
    }

    return match;
}


// Tests if this is in other upper 
bool QuantumIdentifier::InUpper(const QuantumIdentifier& other) const
{
  if(mspecies not_eq other.mspecies)
    return false;
  if(miso not_eq other.miso)
    return false;
  
  if(mqtype == QuantumIdentifier::NONE or other.mqtype == QuantumIdentifier::NONE)
    return false;
  else if(mqtype == QuantumIdentifier::ALL or other.mqtype == QuantumIdentifier::ALL)
    return true;
  else if(mqtype not_eq QuantumIdentifier::ENERGY_LEVEL or other.mqtype not_eq QuantumIdentifier::TRANSITION)
    throw runtime_error("One of your inputs is bad.  You are using function comparing energy levels to the upper state of lines, but the types mismatch");
  
  Index qnri = 0;
  while (qnri not_eq Index(QuantumNumberType::FINAL_ENTRY)) {
    if(other.mqm[TRANSITION_UPPER_INDEX][qnri].isUndefined()) {
      if(not mqm[ENERGY_LEVEL_INDEX][qnri].isUndefined())
        return false;
    }
    else  {
      if(other.mqm[TRANSITION_UPPER_INDEX][qnri] not_eq mqm[ENERGY_LEVEL_INDEX][qnri])
        return false;
    }
    qnri++;
  }
  
  return true;
}



// Tests that if this is in other lower
bool QuantumIdentifier::InLower(const QuantumIdentifier& other) const
{
  if(mspecies not_eq other.mspecies)
    return false;
  if(miso not_eq other.miso)
    return false;
  
  if(mqtype == QuantumIdentifier::NONE or other.mqtype == QuantumIdentifier::NONE)
    return false;
  else if(mqtype == QuantumIdentifier::ALL or other.mqtype == QuantumIdentifier::ALL)
    return true;
  else if(mqtype not_eq QuantumIdentifier::ENERGY_LEVEL or other.mqtype not_eq QuantumIdentifier::TRANSITION)
    throw runtime_error("One of your inputs is bad.  You are using function comparing energy levels to the lower state of lines, but the types mismatch");
  
  Index qnri = 0;
  while (qnri not_eq Index(QuantumNumberType::FINAL_ENTRY)) {
    if(other.mqm[TRANSITION_LOWER_INDEX][qnri].isUndefined()) {
      if(not mqm[ENERGY_LEVEL_INDEX][qnri].isUndefined())
        return false;
    }
    else  {
      if(other.mqm[TRANSITION_LOWER_INDEX][qnri] not_eq mqm[ENERGY_LEVEL_INDEX][qnri])
        return false;
    }
    qnri++;
  }
  
  return true;
}


// Tests that all of this is in other 
bool QuantumIdentifier::In(const QuantumIdentifier& other) const
{
  if(mspecies not_eq other.mspecies)
    return false;
  if(miso not_eq other.miso)
    return false;
  
  if(mqtype == QuantumIdentifier::NONE or other.mqtype == QuantumIdentifier::NONE)
    return false;
  else if(mqtype == QuantumIdentifier::ALL or other.mqtype == QuantumIdentifier::ALL) {}
  else if(mqtype not_eq other.mqtype)
    throw std::runtime_error("Can never compare different types of identifiers with QID.In(QID), one of your inputs is of wrong QuantumIdentifier type");
  else if(mqtype == QuantumIdentifier::TRANSITION) {
    auto& other_low = other.mqm[TRANSITION_LOWER_INDEX];
    auto& other_upp = other.mqm[TRANSITION_UPPER_INDEX];
    auto& this_low = mqm[TRANSITION_LOWER_INDEX];
    auto& this_upp = mqm[TRANSITION_UPPER_INDEX];
    
    for(Index i=0; i<Index(QuantumNumberType::FINAL_ENTRY); i++) {
      if(other_low[i].isUndefined()) {}
      else if(this_low[i].isUndefined()) return false;
      else if(this_low[i] not_eq other_low[i]) return false;
      
      if(other_upp[i].isUndefined()) {}
      else if(this_upp[i].isUndefined()) return false;
      else if(this_upp[i] not_eq other_upp[i]) return false;
    }
  }
  else {
    auto& other_qn = other.mqm[ENERGY_LEVEL_INDEX];
    auto& this_qn = mqm[ENERGY_LEVEL_INDEX];
    
    for(Index i=0; i<Index(QuantumNumberType::FINAL_ENTRY); i++) {
      if(other_qn[i].isUndefined()) {}
      else if(this_qn[i].isUndefined()) return false;
      else if(this_qn[i] not_eq other_qn[i]) return false;
    }
  }
  
  return true;
}


bool IsValidQuantumNumberName(String name)
{
    bool valid = false;
    // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
if (name == #ID) valid = true

        INPUT_QUANTUM(J);
        else INPUT_QUANTUM(dJ);
        else INPUT_QUANTUM(M);
        else INPUT_QUANTUM(N);
        else INPUT_QUANTUM(dN);
        else INPUT_QUANTUM(S);
        else INPUT_QUANTUM(F);
        else INPUT_QUANTUM(K);
        else INPUT_QUANTUM(Ka);
        else INPUT_QUANTUM(Kc);
        else INPUT_QUANTUM(Omega);
        else INPUT_QUANTUM(i);
        else INPUT_QUANTUM(Lambda);
        else INPUT_QUANTUM(alpha);
        else INPUT_QUANTUM(Sym);
        else INPUT_QUANTUM(parity);
        else INPUT_QUANTUM(v1);
        else INPUT_QUANTUM(v2);
        else INPUT_QUANTUM(l2);
        else INPUT_QUANTUM(v3);
        else INPUT_QUANTUM(v4);
        else INPUT_QUANTUM(v5);
        else INPUT_QUANTUM(v6);
        else INPUT_QUANTUM(l);
        else INPUT_QUANTUM(pm);
        else INPUT_QUANTUM(r);
        else INPUT_QUANTUM(S_global);
        else INPUT_QUANTUM(X);
        else INPUT_QUANTUM(n_global);
        else INPUT_QUANTUM(C);
        else INPUT_QUANTUM(Hund);
#undef INPUT_QUANTUM
    return valid;
}

void ThrowIfQuantumNumberNameInvalid(String name)
{
    if (!IsValidQuantumNumberName(name))
    {
        ostringstream os;
        os << "Invalid quantum number: " << name;
        throw std::runtime_error(os.str());
    }
}

std::istream& operator>>(std::istream& is, QuantumNumbers& qn)
{
    String name;
    Rational r;

    is >> name >> r;

    // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
    if (name == #ID) qn.Set(QuantumNumberType::ID, r)

        INPUT_QUANTUM(J);
        else INPUT_QUANTUM(dJ);
        else INPUT_QUANTUM(M);
        else INPUT_QUANTUM(N);
        else INPUT_QUANTUM(dN);
        else INPUT_QUANTUM(S);
        else INPUT_QUANTUM(F);
        else INPUT_QUANTUM(K);
        else INPUT_QUANTUM(Ka);
        else INPUT_QUANTUM(Kc);
        else INPUT_QUANTUM(Omega);
        else INPUT_QUANTUM(i);
        else INPUT_QUANTUM(Lambda);
        else INPUT_QUANTUM(alpha);
        else INPUT_QUANTUM(Sym);
        else INPUT_QUANTUM(parity);
        else INPUT_QUANTUM(v1);
        else INPUT_QUANTUM(v2);
        else INPUT_QUANTUM(l2);
        else INPUT_QUANTUM(v3);
        else INPUT_QUANTUM(v4);
        else INPUT_QUANTUM(v5);
        else INPUT_QUANTUM(v6);
        else INPUT_QUANTUM(l);
        else INPUT_QUANTUM(pm);
        else INPUT_QUANTUM(r);
        else INPUT_QUANTUM(S_global);
        else INPUT_QUANTUM(X);
        else INPUT_QUANTUM(n_global);
        else INPUT_QUANTUM(C);
        else INPUT_QUANTUM(Hund);
    else
    {
        std::ostringstream os;
        os << "Unknown quantum number: " << name << " (" << r << ").";
        throw std::runtime_error(os.str());
    }

#undef INPUT_QUANTUM

    return is;
}


std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn)
{
    bool first = true;
    // Define a helper macro to save some typing.
#define OUTPUT_QUANTUM(ID) \
    if (!qn[Index(QuantumNumberType::ID)].isUndefined()) \
      { if (!first) os << " "; first = false; os << #ID << " " << qn[QuantumNumberType::ID]; }

        OUTPUT_QUANTUM(J);
        OUTPUT_QUANTUM(dJ);
        OUTPUT_QUANTUM(M);
        OUTPUT_QUANTUM(N);
        OUTPUT_QUANTUM(dN);
        OUTPUT_QUANTUM(S);
        OUTPUT_QUANTUM(F);
        OUTPUT_QUANTUM(K);
        OUTPUT_QUANTUM(Ka);
        OUTPUT_QUANTUM(Kc);
        OUTPUT_QUANTUM(Omega);
        OUTPUT_QUANTUM(i);
        OUTPUT_QUANTUM(Lambda);
        OUTPUT_QUANTUM(alpha);
        OUTPUT_QUANTUM(Sym);
        OUTPUT_QUANTUM(parity);
        OUTPUT_QUANTUM(v1);
        OUTPUT_QUANTUM(v2);
        OUTPUT_QUANTUM(l2);
        OUTPUT_QUANTUM(v3);
        OUTPUT_QUANTUM(v4);
        OUTPUT_QUANTUM(v5);
        OUTPUT_QUANTUM(v6);
        OUTPUT_QUANTUM(l);
        OUTPUT_QUANTUM(pm);
        OUTPUT_QUANTUM(r);
        OUTPUT_QUANTUM(S_global);
        OUTPUT_QUANTUM(X);
        OUTPUT_QUANTUM(n_global);
        OUTPUT_QUANTUM(C);
        OUTPUT_QUANTUM(Hund);

#undef OUTPUT_QUANTUM

    return os;
}


String QuantumIdentifier::TypeStr() const {
    String t;
    switch (mqtype) {
        case QuantumIdentifier::TRANSITION:
            t = "TR";
            break;
        case QuantumIdentifier::ENERGY_LEVEL:
            t = "EN";
            break;
        case QuantumIdentifier::ALL:
          t = "ALL";
          break;
        case QuantumIdentifier::NONE:
          t = "NONE";
          break;
        default:
            assert(0);
            break;
    }
    return t;
}

String QuantumIdentifier::SpeciesName() const { return species_name_from_species_index(mspecies); }


void QuantumIdentifier::SetTransition(const QuantumNumbers& upper, const QuantumNumbers& lower)
{
    mqtype = QuantumIdentifier::TRANSITION;
    mqm.resize(2);
    mqm[TRANSITION_UPPER_INDEX] = upper;
    mqm[TRANSITION_LOWER_INDEX] = lower;
}


void QuantumIdentifier::SetEnergyLevel(const QuantumNumbers& q)
{
    mqtype = QuantumIdentifier::ENERGY_LEVEL;
    mqm.resize(1);
    mqm[ENERGY_LEVEL_INDEX] = q;
}


void QuantumIdentifier::SetAll()
{
  mqtype = QuantumIdentifier::ALL;
  mqm.resize(0);
}


void QuantumIdentifier::SetFromString(String str)
{
    // Global species lookup data:
    using global_data::species_data;

    // We need a species index sorted by Arts identifier. Keep this in a
    // static variable, so that we have to do this only once.  The ARTS
    // species index is ArtsMap[<Arts String>].
    static map<String, SpecIsoMap> ArtsMap;

    // Remember if this stuff has already been initialized:
    static bool hinit = false;

    if ( !hinit )
    {
        for ( Index i=0; i<species_data.nelem(); ++i )
        {
            const SpeciesRecord& sr = species_data[i];

            for ( Index j=0; j<sr.Isotopologue().nelem(); ++j)
            {
                SpecIsoMap indicies(i,j);
                String buf = sr.Name()+"-"+sr.Isotopologue()[j].Name();

                ArtsMap[buf] = indicies;
            }
        }
        hinit = true;
    }

    std::istringstream is(str);
    String token;

    is >> token;

    // ok, now for the cool index map:
    // is this arts identifier valid?
    const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(token);
    if ( i == ArtsMap.end() )
    {
        ostringstream os;
        os << "ARTS Tag: " << token << " is unknown.";
        throw runtime_error(os.str());
    }

    SpecIsoMap id = i->second;
    SetSpecies(id.Speciesindex());
    SetIsotopologue(id.Isotopologueindex());

    is >> token;
    if (token == "TR")
    {
        SetType(QuantumIdentifier::TRANSITION);
        is >> token;
        if (token != "UP")
        {
            std::ostringstream os;
            os << "Expected 'UP', but got: " << token;
            throw std::runtime_error(os.str());
        }

        is >> token;
        Rational r;
        while (is)
        {
            ThrowIfQuantumNumberNameInvalid(token);
            is >> r;
            QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX].Set(token, r);
            is >> token;
            if (token == "LO") break;
        }

        if (!is)
        {
            std::ostringstream os;
            os << "Premature end of data, expected 'LO'.";
            throw std::runtime_error(os.str());
        }
        is >> token;
        while (is)
        {
            ThrowIfQuantumNumberNameInvalid(token);
            is >> r;
            QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX].Set(token, r);
            is >> token;
        }
    }
    else if (token == "EN")
    {
        SetType(QuantumIdentifier::ENERGY_LEVEL);

        is >> token;
        Rational r;
        while (is)
        {
            ThrowIfQuantumNumberNameInvalid(token);
            is >> r;
            QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX].Set(token, r);
            is >> token;
        }
    }
    else if (token == "ALL")
    {
      SetType(QuantumIdentifier::ALL);
    }
    else if (token == "NONE")
    {
      SetType(QuantumIdentifier::NONE);
    }
    else
    {
        std::ostringstream os;
        os << "Error parsing QuantumIdentifier. Expected TR or EN, but got: " << token << "\n"
        << "QI: " << str;
        throw std::runtime_error(os.str());
    }
}


void QuantumIdentifier::SetFromStringForCO2Band(String upper, String lower, String iso)
{
  
  assert(upper.nelem() == 5);
  assert(lower.nelem() == 5);
  assert(iso.nelem() == 3);
  
  std::ostringstream os;
  
  os << "CO2-" << iso << " TR " <<
  "UP " <<
  "v1 " << upper[0] << " v2 " << upper[1] << " l2 " << upper[2] << " v3 " << upper[3] << " r " << upper[4] << " " <<
  "LO " <<
  "v1 " << lower[0] << " v2 " << lower[1] << " l2 " << lower[2] << " v3 " << lower[3] << " r " << lower[4];
  
  SetFromString(os.str());
}


std::ostream& operator<<(std::ostream& os, const QuantumIdentifier& qi)
{
  using global_data::species_data;

  const  SpeciesRecord& spr = species_data[qi.Species()];
  
  os << spr.Name() << "-";
  if ( qi.Isotopologue() == spr.Isotopologue().nelem() )
    os << "*";
  else
    os << spr.Isotopologue()[qi.Isotopologue()].Name();
  os << " ";

  if (qi.Type() == QuantumIdentifier::TRANSITION)
      os << "TR UP " << qi.QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX]
         <<   " LO " << qi.QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX];
  else if (qi.Type() == QuantumIdentifier::ENERGY_LEVEL)
      os << "EN " << qi.QuantumMatch()[QuantumIdentifier::ENERGY_LEVEL_INDEX];
  else if (qi.Type() == QuantumIdentifier::ALL)
    os << "ALL";
  else if (qi.Type() == QuantumIdentifier::NONE)
    os << "NONE";
  else
      assert(0);

  return os;
}

bool QuantumIdentifier::any_quantumnumbers() const
{
  Index qni=0;
  switch(mqtype) {
    case QuantumIdentifier::ENERGY_LEVEL:
    case QuantumIdentifier::TRANSITION: 
      for(const auto& qns: mqm)
        do {
          if(not qns[qni].isUndefined())
            return true;
          qni++;
        } while(qni not_eq Index(QuantumNumberType::FINAL_ENTRY));
      break;
    case QuantumIdentifier::ALL:
    case QuantumIdentifier::NONE:
      break;
  }
  return false;
}

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

/**
 * @file quantum.cc
 * @author Oliver Lemke
 * @date 2013-04-12
 * 
 * @brief Classes to handle quantum numbers
 */

#include "quantum.h"
#include <stdexcept>
#include "absorption.h"
#include "global_data.h"
#include "special_interp.h"

void ThrowIfQuantumNumberNameInvalid(String name) {
  if (!IsValidQuantumNumberName(name)) {
    ostringstream os;
    os << "Invalid quantum number: " << name;
    throw std::runtime_error(os.str());
  }
}

std::istream& operator>>(std::istream& is, QuantumNumbers& qn) {
  String name;
  Rational r;

  is >> name >> r;

  qn.Set(name, r);
  
  return is;
}

std::ostream& operator<<(std::ostream& os, const QuantumNumbers& qn) {
  bool first = true;
  for (Index i=0; i<Index(QuantumNumberType::FINAL); i++) {
    if (qn[i].isDefined()) {
      if (first) {
        os << QuantumNumberType(i) << ' ' << qn[i];
        first = false;
      } else {
        os << ' ' << QuantumNumberType(i) << ' ' << qn[i];
      }
    }
  }

  return os;
}

String QuantumIdentifier::SpeciesName() const {
  // Species lookup data:
  using global_data::species_data;
  
  // A reference to the relevant record of the species data:
  const SpeciesRecord& spr = species_data[Species()];
  
  // First the species name:
  return spr.Name() + "-" +
  spr.Isotopologue()[Isotopologue()].Name();
}

Numeric QuantumIdentifier::SpeciesMass() const {
  // Species lookup data:
  using global_data::species_data;
  
  // A reference to the relevant record of the species data:
  const SpeciesRecord& spr = species_data[Species()];
  
  // First the species name:
  return spr.Isotopologue()[Isotopologue()].Mass();
}

void QuantumIdentifier::SetFromString(String str) {
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

  std::istringstream is(str);
  String token;

  is >> token;

  // ok, now for the cool index map:
  // is this arts identifier valid?
  const map<String, SpecIsoMap>::const_iterator i = ArtsMap.find(token);
  if (i == ArtsMap.end()) {
    ostringstream os;
    os << "ARTS Tag: " << token << " is unknown.";
    throw runtime_error(os.str());
  }

  SpecIsoMap id = i->second;
  Species(id.Speciesindex());
  Isotopologue(id.Isotopologueindex());

  is >> token;
  if (token == "TR") {
    mqtype = QuantumIdentifier::TRANSITION;
    is >> token;
    if (token != "UP") {
      std::ostringstream os;
      os << "Expected 'UP', but got: " << token;
      throw std::runtime_error(os.str());
    }

    is >> token;
    Rational r;
    while (is) {
      ThrowIfQuantumNumberNameInvalid(token);
      is >> r;
      QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX].Set(token, r);
      is >> token;
      if (token == "LO") break;
    }

    if (!is) {
      std::ostringstream os;
      os << "Premature end of data, expected 'LO'.";
      throw std::runtime_error(os.str());
    }
    is >> token;
    while (is) {
      ThrowIfQuantumNumberNameInvalid(token);
      is >> r;
      QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX].Set(token, r);
      is >> token;
    }
  } else if (token == "EN") {
    mqtype = QuantumIdentifier::ENERGY_LEVEL;

    is >> token;
    Rational r;
    while (is) {
      ThrowIfQuantumNumberNameInvalid(token);
      is >> r;
      QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX].Set(token, r);
      is >> token;
    }
  } else if (token == "ALL") {
    mqtype = QuantumIdentifier::ALL;
  } else if (token == "NONE") {
    mqtype = QuantumIdentifier::NONE;
  } else {
    std::ostringstream os;
    os << "Error parsing QuantumIdentifier. Expected TR or EN, but got: "
       << token << "\n"
       << "QI: " << str;
    throw std::runtime_error(os.str());
  }
}

void QuantumIdentifier::SetFromStringForCO2Band(String upper,
                                                String lower,
                                                String iso) {
  assert(upper.nelem() == 5);
  assert(lower.nelem() == 5);
  assert(iso.nelem() == 3);

  std::ostringstream os;

  os << "CO2-" << iso << " TR "
     << "UP "
     << "v1 " << upper[0] << " v2 " << upper[1] << " l2 " << upper[2] << " v3 "
     << upper[3] << " r " << upper[4] << " "
     << "LO "
     << "v1 " << lower[0] << " v2 " << lower[1] << " l2 " << lower[2] << " v3 "
     << lower[3] << " r " << lower[4];

  SetFromString(os.str());
}

std::ostream& operator<<(std::ostream& os, const QuantumIdentifier& qi) {
  using global_data::species_data;

  if (qi.Species() < 0 || qi.Isotopologue() < 0)
    return os;

  os << qi.SpeciesName() << ' ';

  if (qi.Type() == QuantumIdentifier::TRANSITION)
    os << "TR UP "
       << qi.QuantumMatch()[QuantumIdentifier::TRANSITION_UPPER_INDEX] << " LO "
       << qi.QuantumMatch()[QuantumIdentifier::TRANSITION_LOWER_INDEX];
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

Rational interpret_stringdata(const QuantumNumberType key, const String& val) {
  if (key == QuantumNumberType::parity) {
    if (val == "+")
      return 1_rat;
    else if (val == "-")
      return -1_rat;
  } else if (key == QuantumNumberType::ElectronState) {
    if (val == "X")
      return Rational(int('X'));
  } else if (key == QuantumNumberType::kronigParity) {
    if (val == "f")
      return Rational(int('f'));
    else if (val == "e")
      return Rational(int('e'));
  } else {
    return Rational(val);
  }
  
  return RATIONAL_UNDEFINED;
}

void update_id(QuantumIdentifier& qid, const std::vector<std::array<String, 2> >& upper_list, const std::vector<std::array<String, 2> >& lower_list)
{
  for (auto& keyval: upper_list) {
    auto key = string2quantumnumbertype(keyval[0]);
    if (key == QuantumNumberType::FINAL) {
      std::ostringstream os;
      os << "The key \"" << keyval[0] << "\" is an invalid input as a quantum number key";
      std::cout << "WARNING: " << os.str() << '\n';
    } else {
      auto val = interpret_stringdata(key, keyval[1]);
      if (val != RATIONAL_UNDEFINED) {
        qid.UpperQuantumNumber(key) = val;
      } else {
        std::ostringstream os;
        os << "The key \"" << keyval[0] << "\" and value \"" << keyval[1] << "\" are invalid input as a quantum number key and value pair";
        std::cout << "WARNING: " << os.str() << '\n';
      }
    }
  }
  
  for (auto& keyval: lower_list) {
    auto key = string2quantumnumbertype(keyval[0]);
    if (key == QuantumNumberType::FINAL) {
      std::ostringstream os;
      os << "The key \"" << keyval[0] << "\" is an invalid input as a quantum number key";
      std::cout << "WARNING: " << os.str() << '\n';
    } else {
      auto val = interpret_stringdata(key, keyval[1]);
      if (val != RATIONAL_UNDEFINED) {
        qid.LowerQuantumNumber(key) = val;
      } else {
        std::ostringstream os;
        os << "The key \"" << keyval[0] << "\" and value \"" << keyval[1] << "\" are invalid input as a quantum number key and value pair";
        std::cout << "WARNING: " << os.str() << '\n';
      }
    }
  }
}

String QuantumNumbers::toString() const
{
  std::ostringstream out;
  out << (*this) << ' ';
  String s=out.str();
  if(s.back() == ' ')
    s.pop_back();
  return s;
}

QuantumIdentifier::QuantumIdentifier(const Index spec,
                                     const Index isot,
                                     const std::vector<QuantumNumberType>& keys,
                                     const std::vector<Rational>& upper,
                                     const std::vector<Rational>& lower) :
mqtype(QuantumIdentifier::TRANSITION),
mspecies(spec),
miso(isot) {
  for(size_t i=0; i<keys.size(); i++) {
    mqm[TRANSITION_UPPER_INDEX][keys[i]] = upper[i];
    mqm[TRANSITION_LOWER_INDEX][keys[i]] = lower[i];
  }
}

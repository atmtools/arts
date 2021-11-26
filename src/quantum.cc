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
#include "absorption.h"
#include "debug.h"
#include "global_data.h"
#include "special_interp.h"
#include <stdexcept>

void ThrowIfQuantumNumberNameInvalid(String name) {
  ARTS_USER_ERROR_IF (!IsValidQuantumNumberName(name),
    "Invalid quantum number: ", name)
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

void Quantum::Identifier::SetFromString(String str) {
  std::istringstream is(str);
  String token;

  is >> token;

  spec_ind = Species::Tag(token).spec_ind;

  is >> token;
  if (token == "TR") {
    type = Quantum::IdentifierType::Transition;
    is >> token;
    ARTS_USER_ERROR_IF (token != "UP",
      "Expected 'UP', but got: ", token)

    is >> token;
    Rational r;
    while (is) {
      if (token == "LO") break;
      ThrowIfQuantumNumberNameInvalid(token);
      is >> r;
      Upper().Set(token, r);
      is >> token;
    }

    ARTS_USER_ERROR_IF (!is,
      "Premature end of data, expected 'LO'.")
    is >> token;
    while (is) {
      ThrowIfQuantumNumberNameInvalid(token);
      is >> r;
      Lower().Set(token, r);
      is >> token;
    }
  } else if (token == "EN") {
    type = Quantum::IdentifierType::EnergyLevel;

    is >> token;
    Rational r;
    while (is) {
      ThrowIfQuantumNumberNameInvalid(token);
      is >> r;
      Level().Set(token, r);
      is >> token;
    }
  } else if (token == "ALL") {
    type = Quantum::IdentifierType::All;
  } else if (token == "NONE") {
    type = Quantum::IdentifierType::None;
  } else {
    ARTS_USER_ERROR (
      "Error parsing QuantumIdentifier. Expected TR or EN, but got: ",
      token, "\n"
      "QI: ", str)
  }
}

String Quantum::Identifier::GetString() const
{
  std::ostringstream os;
  os << *this;
  return os.str();
}


Rational interpret_stringdata(const QuantumNumberType key, const String& val) {
  if (key == QuantumNumberType::parity) {
    if (val == "+")
      return 1;
    if (val == "-")
      return -1;
  }

  if (key == QuantumNumberType::ElectronState) {
    if (val == "X")
      return Rational(int('X'));
    if (val == "a")
      return Rational(int('a'));
    if (val == "b")
      return Rational(int('b'));
    if (val == "c")
      return Rational(int('c'));
    if (val == "A")
      return Rational(int('A'));
    if (val == "'")
      return Rational(int('\''));
    if (val == "B")
      return Rational(int('B'));
  }
  
  if (key == QuantumNumberType::kronigParity) {
    if (val == "f")
      return Rational(int('f'));
    if (val == "e")
      return Rational(int('e'));
  }

  try {
    return Rational(val);
  } catch (std::runtime_error& e) {
    ARTS_USER_ERROR("Key: ", key, "\nFailed\n", e.what());
  }
}

void update_id(QuantumIdentifier& qid, const std::vector<std::array<String, 2> >& upper_list, const std::vector<std::array<String, 2> >& lower_list)
{
  for (auto& keyval: upper_list) {
    auto key = string2quantumnumbertype(keyval[0]);
    if (key == QuantumNumberType::FINAL) {
    } else {
      auto val = interpret_stringdata(key, keyval[1]);
      if (val != RATIONAL_UNDEFINED) {
        qid.Upper()[key] = val;
      } else {
      }
    }
  }
  
  for (auto& keyval: lower_list) {
    auto key = string2quantumnumbertype(keyval[0]);
    if (key == QuantumNumberType::FINAL) {
    } else {
      auto val = interpret_stringdata(key, keyval[1]);
      if (val != RATIONAL_UNDEFINED) {
        qid.Lower()[key] = val;
      } else {
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

/* Copyright (C) 2012 
Richard Larsson <ric.larsson@gmail.com>

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
 * @file   rational.cc
 * @author Richard Larsson
 * @date   2012-10-31
 * 
 * @brief  Contains the rational class implmentations
 **/

#include "rational.h"

#include <ostream>
#include <stdexcept>
#include "mystring.h"

std::ostream& operator<<(std::ostream& os, const Rational& a) {
  Rational r = reduce_by_gcd(a);
  r.fixSign();

  if (r.Denom() == 1)
    os << r.Nom();
  else
    os << r.Nom() << "/" << r.Denom();
  return os;
}

std::istream& operator>>(std::istream& is, Rational& a) {
  String s;
  Index nom;
  Index denom;
  char* endptr;

  is >> s;

  try {
    ArrayOfString as;

    s.split(as, "/");

    if (as.nelem() == 1) {
      nom = strtol(s.c_str(), &endptr, 10);
      ARTS_USER_ERROR_IF (endptr != s.c_str() + s.nelem(),
                          "Error parsing rational number");
      a = Rational(nom, 1);
    } else if (as.nelem() == 2) {
      nom = strtol(as[0].c_str(), &endptr, 10);
      ARTS_USER_ERROR_IF (endptr != as[0].c_str() + as[0].nelem(),
                          "Error parsing rational number nominator");
      denom = strtol(as[1].c_str(), &endptr, 10);
      ARTS_USER_ERROR_IF (endptr != as[1].c_str() + as[1].nelem(),
                          "Error parsing rational number denominator");
      a = Rational(nom, denom);
    } else {
      ARTS_USER_ERROR_IF (true, "Error parsing rational number");
    }
  } catch (const std::runtime_error& e) {
    ARTS_USER_ERROR_IF (true,
                        "Error parsing rational number: ",
                        s, '\n', e.what())
  }

  return is;
}


Rational::Rational(const String& s)
{
  auto len = s.length();
  
  if (len) {
    auto dot_pos = s.find(".");
    auto slash_pos = s.find("/");
    if (len > dot_pos) {
      *this = numeric2rational(std::stod(s), len - dot_pos - 1);
    } else if (len > slash_pos) {
      *this = Rational(std::stoi(s.substr(0, slash_pos)),
                       std::stoi(s.substr(slash_pos + 1, len)));
    } else {
      *this = Rational(std::stoi(s));
    }
  } else {
    *this = RATIONAL_UNDEFINED;
  }
}


void Rational::simplify_in_place() noexcept
{
  Rational a = reduce_by_gcd(*this);
  mnom = a.Nom();
  mdenom = a.Denom();
  fixSign();
}

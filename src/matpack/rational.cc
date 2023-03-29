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

#include "debug.h"
#include "mystring.h"
#include <ostream>
#include <stdexcept>

std::ostream& operator<<(std::ostream& os, const Rational& a) {
  Rational r = reduce_by_gcd(a);
  r.fixSign();

  if (r.denom == 1)
    os << r.numer;
  else
    os << r.numer << "/" << r.denom;
  return os;
}

std::istream& operator>>(std::istream& is, Rational& a) {
  String s;

  is >> s;
  a = Rational(s);

  ARTS_USER_ERROR_IF(a.isUndefined(), "Cannot read ", s, " as a rational")

  return is;
}


Rational::Rational(const String& s) {
  auto len = s.length();
  
  if (len) {
    auto dot_pos = s.find(".");
    auto slash_pos = s.find("/");
    if (len > dot_pos) {
      *this = numeric2rational(std::stod(s), len - dot_pos - 1);
    } else if (len > slash_pos) {
      const String a{s.substr(0, slash_pos)};
      const String b{s.substr(slash_pos + 1, len)};
      try {
        *this = Rational(std::stoi(a), std::stoi(b));
      } catch (...) {
        ARTS_USER_ERROR("Cannot interpret either '", a, "' or '", b, "' as an integer (or both)");
      }
    } else {
      try {
        *this = Rational(std::stoi(s));
      } catch (...) {
        ARTS_USER_ERROR("Cannot interpret '", s, "' as an integer");
      }
    }
  } else {
    *this = RATIONAL_UNDEFINED;
  }
}


void Rational::simplify_in_place() noexcept {
  Rational a = reduce_by_gcd(*this);
  numer = a.numer;
  denom = a.denom;
  fixSign();
}

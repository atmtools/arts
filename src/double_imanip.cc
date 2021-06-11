/* Copyright (C) 2021 Oliver Lemke <oliver.lemke@uni-hamburg.de>

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

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/**
   \file  double_imanip.cc

   IO manipulation classes for parsing nan and inf.

   Not all compilers do support parsing of nan and inf.
   The code below is taken (and slightly modified) from the discussion at:
   http://stackoverflow.com/questions/11420263/is-it-possible-to-read-infinity-or-nan-values-using-input-streams

   \author Oliver Lemke
   \date 2021-06-03
*/

////////////////////////////////////////////////////////////////////////////
//   External declarations
////////////////////////////////////////////////////////////////////////////

#include "double_imanip.h"

#include <limits>


////////////////////////////////////////////////////////////////////////////
//   IO manipulation classes for parsing nan and inf
////////////////////////////////////////////////////////////////////////////

double_istream& double_istream::parse_on_fail(double& x, bool neg) {
  const char* exp[] = {"", "inf", "Inf", "nan", "NaN"};
  const char* e = exp[0];
  int l = 0;
  char inf[4] = "\0\0\0";
  char* c = inf;
  if (neg) *c++ = '-';
  in.clear();
  if (!(in >> *c).good()) return *this;
  switch (*c) {
    case 'i':
      e = exp[l = 1];
      break;
    case 'I':
      e = exp[l = 2];
      break;
    case 'n':
      e = exp[l = 3];
      break;
    case 'N':
      e = exp[l = 4];
      break;
  }
  while (*c == *e) {
    if ((e - exp[l]) == 2) break;
    ++e;
    if (!(in >> *++c).good()) break;
  }
  if (in.good() && *c == *e) {
    switch (l) {
      case 1:
      case 2:
        x = std::numeric_limits<double>::infinity();
        break;
      case 3:
      case 4:
        x = std::numeric_limits<double>::quiet_NaN();
        break;
    }
    if (neg) x = -x;
    return *this;
  } else if (!in.good()) {
    if (!in.fail()) return *this;
    in.clear();
    --c;
  }
  do {
    in.putback(*c);
  } while (c-- != inf);
  in.setstate(std::ios_base::failbit);
  return *this;
}

const double_imanip& operator>>(std::istream& in, const double_imanip& dm) {
  dm.in = &in;
  return dm;
}

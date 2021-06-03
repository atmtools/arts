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
   \file  double_imanip.h

   IO manipulation classes for parsing nan and inf.

   Not all compilers do support parsing of nan and inf.
   The code below is taken (and slightly modified) from the discussion at:
   http://stackoverflow.com/questions/11420263/is-it-possible-to-read-infinity-or-nan-values-using-input-streams

   \author Oliver Lemke
   \date 2021-06-03
*/

#ifndef double_imanip_h
#define double_imanip_h

#include <fstream>

/** Input stream class for doubles that correctly handles nan and inf. */
class double_istream {
 public:
  double_istream(std::istream& i) : in(i) {}

  double_istream& parse_on_fail(double& x, bool neg);

  double_istream& operator>>(double& x) {
    bool neg = false;
    char c;
    if (!in.good()) return *this;
    while (isspace(c = (char)in.peek())) in.get();
    if (c == '-') {
      neg = true;
    }
    in >> x;
    if (!in.fail()) return *this;
    return parse_on_fail(x, neg);
  }

 private:
  std::istream& in;
};

/** Input manipulator class for doubles to enable nan and inf parsing. */
class double_imanip {
 public:
  const double_imanip& operator>>(double& x) const {
    double_istream(*in) >> x;
    return *this;
  }
  std::istream& operator>>(const double_imanip&) const { return *in; }

  friend const double_imanip& operator>>(std::istream& in,
                                         const double_imanip& dm);

 private:
  mutable std::istream* in;
};

const double_imanip& operator>>(std::istream& in, const double_imanip& dm);

#endif

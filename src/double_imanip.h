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

   Fast double input stream with support for parsing nan and inf.

   \author Oliver Lemke, Richard Larsson
   \date 2022-05-23
*/

#ifndef double_imanip_h
#define double_imanip_h

#include <fstream>

#include "debug.h"

/** Input manipulator class for doubles to enable nan and inf parsing. */
class double_imanip {
 public:
  const double_imanip& operator>>(double& x) const {
    std::istream& is = *in;
    std::string buf;
    try {
      is >> buf;
      ARTS_USER_ERROR_IF(is.fail(), "Cannot read from stream")
      std::size_t n;
      x = std::stod(buf, &n);
      while (n++ not_eq buf.size()) is.unget();
    } catch (std::runtime_error& e) {
      ARTS_USER_ERROR(R"--(Failed conversion to double

)--",
                      e.what())
    } catch (std::invalid_argument&) {
      ARTS_USER_ERROR("The argument: \n\n'", buf, R"--('

is not convertible to a valid double.  At the very least it
cannot be converted to one using the standard string-to-double
routine
)--")
    } catch (std::out_of_range&) {
      ARTS_USER_ERROR("The argument: \n\n'",
                      buf,
                      R"--('

is not convertible to a normal double.  The value is
bad or simply subnormal (i.e., the absolute value is either
below ~)--",
                      std::numeric_limits<double>::min(),
                      " or above ~",
                      std::numeric_limits<double>::max(),
                      ")\n")
    }
    return *this;
  }

  std::istream& operator>>(const double_imanip&) const { return *in; }

  friend const double_imanip& operator>>(std::istream& in,
                                         const double_imanip& dm);

 private:
  mutable std::istream* in;
};

inline const double_imanip& operator>>(std::istream& in,
                                       const double_imanip& dm) {
  dm.in = &in;
  return dm;
}

#endif

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
#include <limits>
#include <stdexcept>

#include "debug.h"

#include <fast_float/fast_float.h>

/** Input manipulator class for doubles to enable nan and inf parsing. */
class double_imanip {
 public:
  const double_imanip& operator>>(double& x) const {
    std::istream& is = *in;
    std::string buf;

    // Read to the buffer
    is >> buf;
    ARTS_USER_ERROR_IF(is.fail(), "Cannot read from stream")

    // Actual conversion
    const auto res =
        fast_float::from_chars(buf.c_str(), buf.c_str() + buf.size(), x);

    // Error (only std::errc::invalid_argument possible)
    ARTS_USER_ERROR_IF(res.ec == std::errc::invalid_argument,
                       "The argument: \n\n'",
                       buf,
                       R"--('

is not convertible to a valid double.  At the very least it
cannot be converted to one using the standard string-to-double
routine
)--")

    // Put the stream to be where it is supposed to be
    std::size_t n = std::distance(buf.c_str(), res.ptr);
    while (n++ < buf.size()) is.unget();

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

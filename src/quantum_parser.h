/* Copyright (C) 2015, The ARTS Developers.

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

  Parser for quantum numbers from spectroscopic catalogs.

  \author Oliver Lemke
*/

#ifndef quantum_parser_h
#define quantum_parser_h

#include "quantum.h"

//! Function pointer type for quantum number parsing routines
typedef void (*QuantumParseFunction)(Rational& qn,
                                     String& s,
                                     const Index species);

//! Class mapping quantum numbers to parsing functions
class QuantumFieldDescription {
 public:
  QuantumFieldDescription() {}

  QuantumFieldDescription(QuantumNumberType quantum_id,
                          QuantumParseFunction qpfunc)
      : mquantum_id(quantum_id), mqpfunc(qpfunc) {}

  void Parse(QuantumNumbers& qnr, String& s, const Index species) const {
    Rational qn;
    mqpfunc(qn, s, species);
    if (!qn.isUndefined()) qnr.Set(mquantum_id, qn);
  }

 private:
  QuantumNumberType mquantum_id;
  QuantumParseFunction mqpfunc;
};

#endif /* quantum_parser_h */

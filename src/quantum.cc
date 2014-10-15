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

#include "quantum.h"

#include <stdexcept>
#include "mystring.h"


bool QuantumNumbers::Compare(const QuantumNumbers& qn) const
{
    const QuantumContainer& qnumbers2 = qn.GetNumbers();

    bool match = true;

    Index qnri = 0;

    // Compare all quantum numbers in mqnumbers and qnumbers2
    while (match && qnri != QN_FINAL_ENTRY)
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


std::istream& operator>>(std::istream& is, QuantumNumbers& qn)
{
    String name;
    Rational r;

    is >> name >> r;

    // Define a helper macro to save some typing.
#define INPUT_QUANTUM(ID) \
    if (name == #ID) qn.Set(QN_ ## ID, r)

    INPUT_QUANTUM(J);
    else INPUT_QUANTUM(N);
    else INPUT_QUANTUM(S);
    else INPUT_QUANTUM(F);
    else INPUT_QUANTUM(Omega);
    else INPUT_QUANTUM(K1);
    else INPUT_QUANTUM(K2);
    else INPUT_QUANTUM(v1);
    else INPUT_QUANTUM(v2);
    else INPUT_QUANTUM(v3);
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
    if (!qn[QN_ ## ID].isUndefined()) \
      { if (!first) os << " "; first = false; os << #ID << " " << qn[QN_ ## ID]; }

    OUTPUT_QUANTUM(J);
    OUTPUT_QUANTUM(N);
    OUTPUT_QUANTUM(S);
    OUTPUT_QUANTUM(F);
    OUTPUT_QUANTUM(Omega);
    OUTPUT_QUANTUM(K1);
    OUTPUT_QUANTUM(K2);
    OUTPUT_QUANTUM(v1);
    OUTPUT_QUANTUM(v2);
    OUTPUT_QUANTUM(v3);

#undef OUTPUT_QUANTUM

    return os;
}

std::ostream& operator<<(std::ostream& os, const QuantumNumberRecord& qr)
{
    os << "Upper: " << qr.Upper() << " ";
    os << "Lower: " << qr.Lower();

    return os;
}


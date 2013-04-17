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


bool QuantumNumbers::Compare(const QuantumNumbers& qn) const
{
    const ArrayOfRational& qnumbers2 = qn.GetNumbers();

    assert(mqnumbers.nelem() == qnumbers2.nelem());

    bool match = true;

    // Compare Upper Quantum Numbers
    ArrayOfRational::const_iterator qnr1it = mqnumbers.begin();
    ArrayOfRational::const_iterator qnr2it = qnumbers2.begin();

    while (match && qnr1it != mqnumbers.end() && qnr2it != qnumbers2.end())
    {
        if (*qnr1it != *qnr2it
            && !(*qnr1it).isUndefined()
            && !(*qnr2it).isUndefined())
            match = false;

        qnr1it++;
        qnr2it++;
    }

    return match;
}

ostream& operator<<(ostream& os, const QuantumNumbers& qn)
{
    os << qn.GetNumbers() << endl;

    return os;
}

ostream& operator<<(ostream& os, const QuantumNumberRecord& qr)
{
    os << "Upper: " << qr.Upper().GetNumbers() << " ";
    os << "Lower: " << qr.Lower().GetNumbers() << endl;

    return os;
}


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
    const QuantumContainer& qnumbers2 = qn.GetNumbers();

    bool match = true;

    QuantumContainer::const_iterator qnr1it = mqnumbers.begin();
    QuantumContainer::const_iterator qnr2it;

    // Loop over all quantum numbers in mqnumbers and use their keys
    // to find the corresponding numbers in qnumbers2
    while (match && qnr1it != mqnumbers.end())
    {
        qnr2it = qnumbers2.find(qnr1it->first);

        // If one of the two numbers is undefined, it is considered as
        // a match.
        if (qnr2it != qnumbers2.end()
            && !(qnr1it->second).isUndefined()
            && !(qnr2it->second).isUndefined()
            && qnr1it->second != qnr2it->second)
            match = false;

        qnr1it++;
    }

    return match;
}

ostream& operator<<(ostream& os, const QuantumNumbers& qn)
{
    for (Index i = 0; i < QN_FINAL_ENTRY; i++)
        os << qn[i] << " ";

    return os;
}

ostream& operator<<(ostream& os, const QuantumNumberRecord& qr)
{
    os << "Upper: " << qr.Upper() << " ";
    os << "Lower: " << qr.Lower();

    return os;
}


/* Copyright (C) 2003-2012 Oliver Lemke <olemke@core-dump.info>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

#include <iostream>
#include "arts.h"
#include "auto_md.h"
#include "matpackII.h"
#include "xml_io.h"
#include "exceptions.h"
#include "absorption.h"

extern Array<SpeciesRecord> species_data;

int
main (int /*argc*/, char * /*argv*/ [])
{
    QuantumNumbers q1;
    QuantumNumbers q2;

    q1.Set(QN_J, Rational(1,2));
    q1.Set(QN_S, Rational(1,2));

    q2.Set(QN_J, Rational(1,2));
    q2.Set(QN_N, Rational(1,3));
    q2.Set(QN_S, Rational(1,2));

    cout << "Compare q1==q2: " << q1.Compare(q2) << endl;
    cout << "Compare q2==q1: " << q2.Compare(q1) << endl;

    define_species_data();
    define_species_map();

    Verbosity v(2,2,2);

    ArrayOfLineRecord abs_lines;
    Timer timer;

    timerStart(timer, v);
    abs_linesReadFromHitran(abs_lines,
                            "/Users/olemke/Dropbox/Hacking/sat/catalogue/HITRAN2008/HITRAN08.par",
                            1, 3e12, v);
//    118e9, 119e9, v);
    timerStop(timer, v);

    Print(timer, 1, v);


    Index species = species_index_from_species_name("O2");
    Index iso = 0;
    QuantumNumberRecord qnr;
    qnr.SetUpper(QN_J, Rational(57, 1));
    qnr.SetLower(QN_J, Rational(58, 1));

    ArrayOfIndex matches;

    timerStart(timer, v);
    find_matching_lines(matches, abs_lines, species, iso, qnr);
    timerStop(timer, v);
    Print(timer, 1, v);

    cout << "Matches: " << matches.nelem() << endl;

    for (Index i = 0; i < matches.nelem(); i++)
    {
        cout << abs_lines[matches[i]] << endl;
        cout << abs_lines[matches[i]].QuantumNumbers() << endl;
    }

    cout << qnr << endl;

    return (0);
}

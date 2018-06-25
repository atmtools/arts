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
#include "m_xml.h"
#include "xml_io.h"
#include "exceptions.h"
#include "absorption.h"

int
main (int /*argc*/, char * /*argv*/ [])
{
    try {
        QuantumNumbers q1;
        QuantumNumbers q2;

        q1.Set(QuantumNumberType::J, Rational(1,2));
        q1.Set(QuantumNumberType::S, Rational(1,2));

        q2.Set(QuantumNumberType::J, Rational(1,2));
        q2.Set(QuantumNumberType::N, Rational(1,3));
        q2.Set(QuantumNumberType::S, Rational(1,2));

        cout << "Compare q1==q2: " << q1.Compare(q2) << endl;
        cout << "Compare q2==q1: " << q2.Compare(q1) << endl;

        ostringstream os;
        os << q2;
        cout << "q1: " << q1 << endl;
        cout << "q2: " << q2 << endl;

        WriteXML("ascii", q1, "quantum.xml", 0,
                 "q1", "", "", Verbosity(0, 2, 0));

        QuantumNumbers q3;
        ReadXML(q3, "q3", "quantum.xml", "", Verbosity(0,2,0));
        cout << "q3: " << q3 << endl;

        cout << endl << "========================================" << endl << endl;

        define_species_data();
        define_species_map();

        Verbosity v(2,2,2);

        ArrayOfLineRecord abs_lines;
        Timer timer;

        timerStart(timer, v);
        abs_linesReadFromHitran(abs_lines,
                                "/Users/olemke/Dropbox/Hacking/sat/catalogue/HITRAN2008/HITRAN08.par",
                                1, 1.1876e+11, v);
        //    1, 3e12, v);
        //    118e9, 119e9, v);
        timerStop(timer, v);

        Print(timer, 1, v);


        SpeciesTag stag("O2-66");
        QuantumNumberRecord qnr;
        qnr.SetUpper(QuantumNumberType::J, Rational(1, 1));
        qnr.SetLower(QuantumNumberType::J, Rational(1, 1));

        cout << "qnr:  " << qnr << endl;

        WriteXML("ascii", qnr, "quantumrecord.xml", 0,
                 "qnr", "", "", Verbosity(0, 2, 0));

        QuantumNumberRecord qnr2;
        ReadXML(qnr2, "qnr2", "quantumrecord.xml", "", Verbosity(0,2,0));
        cout << "qnr2: " << qnr2 << endl;

        
        ArrayOfIndex matches;

        timerStart(timer, v);
        find_matching_lines(matches, abs_lines, stag.Species(), stag.Isotopologue(), qnr);
        timerStop(timer, v);
        Print(timer, 1, v);

        cout << "========================================" << endl << endl;

        cout << "Matches: " << matches.nelem() << endl;

        for (Index i = 0; i < matches.nelem(); i++)
        {
            cout << abs_lines[matches[i]] << endl;
            cout << abs_lines[matches[i]].QuantumIdentity() << endl;
        }

        cout << "========================================" << endl << endl;
        
        cout << "Search key: " << endl;
        cout << "Species: " << stag.Name() << " species: " << stag.Species() << " iso: " << stag.Isotopologue() << endl;
        cout << qnr << endl;
        
    }
    catch (const std::runtime_error &e)
    {
        cout << e.what() << endl;
    }


    return (0);
}

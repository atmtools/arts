/* Copyright (C) 2013
   Oliver Lemke  <olemke@core-dump.info>

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

#include "arts.h"
#include "absorption.h"
#include "file.h"

#define NUM_LINE_MIXING_PARAMS 10

/* Workspace method: Doxygen documentation will be auto-generated */
void line_mixing_o2Init(// WS Output:
                        ArrayOfVector& line_mixing_o2,
                        ArrayOfArrayOfIndex& line_mixing_o2_lut,
                        // WS Input:
                        const ArrayOfArrayOfSpeciesTag& abs_species,
                        // Verbosity object:
                        const Verbosity&)
{
    line_mixing_o2.resize(0);
    line_mixing_o2_lut.resize(abs_species.nelem());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void line_mixing_o2Read(// WS Output:
                        ArrayOfVector& line_mixing_o2,
                        ArrayOfArrayOfIndex& line_mixing_o2_lut,
                        // WS Input:
                        const ArrayOfArrayOfLineRecord& abs_lines_per_species,
                        const ArrayOfArrayOfSpeciesTag& abs_species,
                        const String& filename,
                        const Verbosity& verbosity)
{
    CREATE_OUT2;
    CREATE_OUT3;

    if (abs_species.nelem() != line_mixing_o2_lut.nelem())
        throw runtime_error( "*line_mixing_o2_lut* doesn't match *abs_species*.\n"
                             "Make sure to call line_mixing_o2Init first." );
    ifstream ifs;
    String line;
    QuantumNumberRecord qnr;
    Array<QuantumNumberRecord> aqnr;
    Rational r;
    istringstream is;

    Index linenr = 0;

    Index species_o2_index = species_index_from_species_name("O2");
    // Find index of O2 in abs_lines_per_species
    Index abs_o2_index = -1;
    for (Index sp = 0; sp < abs_species.nelem(); sp++)
    {
        // FIXME: OLE: Temporarily we just look for an O2 tag. This will
        // change in the near future once we have a special linemixing tag
        if (abs_species[sp].nelem()
            && abs_species[sp][0].Species() == species_o2_index)
        {
            abs_o2_index = sp;
            break;
        }
    }

    if (abs_o2_index == -1)
        throw runtime_error("Can't find O2 in *abs_species*");

    open_input_file(ifs, filename);

    Vector params;
    params.resize(NUM_LINE_MIXING_PARAMS);

    line_mixing_o2.resize(0);

    // First we read the whole line mixing file and also
    // temporarily store the quantum numbers from the file
    while (!ifs.eof())
    {
        getline(ifs, line);
        linenr++;

        line.trim();
        // Skip empty lines and comments
        if (!line.nelem()
            || (line.nelem() && line[0] == '#'))
            continue;

        is.clear();
        is.str(line);
        try {
            is >> r; qnr.SetLower(QN_v1, r);
            is >> r; qnr.SetUpper(QN_N, r);
            is >> r; qnr.SetLower(QN_N, r);
            is >> r; qnr.SetUpper(QN_J, r);
            is >> r; qnr.SetLower(QN_J, r);
            aqnr.push_back(qnr);

            params = NAN;
            for (Index i = 0; i < NUM_LINE_MIXING_PARAMS; i++)
                is >> params[i];

            line_mixing_o2.push_back(params);
            
        } catch (runtime_error e) {
            ostringstream os;

            os << "Error parsing line mixing file in line " << linenr << endl;
            os << e.what();
            throw runtime_error(os.str());
        }
    }

    ArrayOfIndex matches;
    line_mixing_o2_lut[abs_o2_index].resize(abs_lines_per_species[abs_o2_index].nelem());
    line_mixing_o2_lut[abs_o2_index] = -1;

    // Now we use the quantum numbers to match the line mixing
    // data to lines in abs_lines_per_species
    Index nmatches = 0;
    for (Index i = 0; i < line_mixing_o2.nelem(); i++)
    {
        find_matching_lines(matches,
                            abs_lines_per_species[abs_o2_index],
                            species_o2_index,
                            -1,
                            aqnr[i]);

        if (!matches.nelem())
        {
            out3 << "  Found no matching lines for\n" << aqnr[i] << "\n";
        }
        else if (matches.nelem() == 1)
        {
            out3 << "  Found matching line for\n" << aqnr[i] << "\n";
            line_mixing_o2_lut[abs_o2_index][matches[0]] = i;
            nmatches++;
        }
        else
        {
            ostringstream os;
            os << "  Found multiple lines for\n" << aqnr[i] << endl
            << "  Matching lines are: " << endl;
            for (Index m = 0; m < matches.nelem(); m++)
                os << "  " << abs_lines_per_species[abs_o2_index][matches[m]] << endl
                << "  " << abs_lines_per_species[abs_o2_index][matches[m]].QuantumNumbers() << endl;
            throw runtime_error(os.str());
        }
    }

    out2 << "  Matched " << nmatches << " lines out of " << line_mixing_o2.nelem() << "\n";
    out2 << "  abs_lines_per_species contains " << abs_lines_per_species[abs_o2_index].nelem()
    << " lines.\n";
}


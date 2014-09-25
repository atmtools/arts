/* Copyright 2013, The ARTS Developers.
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "arts.h"
#include "absorption.h"
#include "file.h"
#include "linemixingrecord.h"


/* Workspace method: Doxygen documentation will be auto-generated */
void line_mixing_dataInit(// WS Output:
                          ArrayOfArrayOfLineMixingRecord& line_mixing_data,
                          ArrayOfArrayOfIndex& line_mixing_data_lut,
                          // WS Input:
                          const ArrayOfArrayOfSpeciesTag& abs_species,
                          // Verbosity object:
                          const Verbosity&)
{
    line_mixing_data.resize(abs_species.nelem());
    line_mixing_data_lut.resize(abs_species.nelem());
}


/* Workspace method: Doxygen documentation will be auto-generated */
void line_mixing_dataMatch(// WS Output:
                           //ArrayOfArrayOfLineMixingRecord& line_mixing_data,
                           //ArrayOfArrayOfIndex& line_mixing_data_lut,
                           // WS Output/Input:
                           ArrayOfArrayOfLineRecord& abs_lines_per_species,
                           // WS Input:
                           const ArrayOfArrayOfSpeciesTag& abs_species,
                           const String& species_tag,
                           const String& line_mixing_tag,
                           const ArrayOfLineMixingRecord& line_mixing_records,
                           const Verbosity& verbosity)
{
    CREATE_OUT2;
    CREATE_OUT3;
/*
    if (abs_species.nelem() != line_mixing_data.nelem())
        throw std::runtime_error( "*line_mixing_data* doesn't match *abs_species*.\n"
                             "Make sure to call line_mixing_dataInit first." );
    if (abs_species.nelem() != line_mixing_data_lut.nelem())
        throw std::runtime_error( "*line_mixing_data_lut* doesn't match *abs_species*.\n"
                            "Make sure to call line_mixing_dataInit first." );*/
    //std::cout<<"help\n";
    LineMixingData line_mixing_data_holder;
    line_mixing_data_holder.StorageTag2SetType(line_mixing_tag);

    // Find index of species_tag in abs_species
    SpeciesTag this_species( species_tag );
    Index species_index = -1;
    for (Index i = 0; species_index == -1 && i < abs_species.nelem(); i++)
        // We only look at the first SpeciesTag in each group because
        // line mixing tags can not be combined with other tags
        if (this_species == abs_species[i][0])
            species_index = i;

    if (species_index == -1)
    {
        ostringstream os;
        os << "Can't find tag \"" << species_tag << "\" in *abs_species*.";
        throw std::runtime_error(os.str());
    }
ArrayOfArrayOfLineMixingRecord line_mixing_data(abs_species.nelem());//FIXME: Added for convenience to match old use case
    line_mixing_data[species_index] = line_mixing_records;

    ArrayOfIndex matches;

    // Now we use the quantum numbers to match the line mixing
    // data to lines in abs_lines_per_species
    Index nmatches = 0;
    for (Index i = 0; i < line_mixing_data[species_index].nelem(); i++)
    {
        const LineMixingRecord& this_lmr = line_mixing_data[species_index][i];
        
        find_matching_lines(matches,
                            abs_lines_per_species[species_index],
                            this_lmr.Species(),
                            this_lmr.Isotopologue(),
                            this_lmr.Quantum());
        
        line_mixing_data_holder.SetDataFromVectorWithKnownType(line_mixing_data[species_index][i].Data());
        
        if (!matches.nelem())
        {
            out3 << "  Found no matching lines for\n" << this_lmr.Quantum() << "\n";
        }
        else if (matches.nelem() == 1)
        {
            out3 << "  Found matching line for\n" << this_lmr.Quantum() << "\n";
            abs_lines_per_species[species_index][matches[0]].SetLineMixingData(line_mixing_data_holder);
            nmatches++;
        }
        else
        {
            ostringstream os;
            os << "  Found multiple lines for\n" << this_lmr.Quantum() << std::endl
            << "  Matching lines are: " << std::endl;
            for (Index m = 0; m < matches.nelem(); m++)
                os << "  " << abs_lines_per_species[species_index][matches[m]] << std::endl
                << "  " << abs_lines_per_species[species_index][matches[m]].QuantumNumbers()
                << std::endl;
            throw std::runtime_error(os.str());
        }
    }

    out2 << "  Matched " << nmatches << " lines out of " << line_mixing_data[species_index].nelem()
         << "\n";
    out2 << "  abs_lines_per_species contains " << abs_lines_per_species[species_index].nelem()
         << " lines for " << species_tag << ".\n";
}


/* Workspace method: Doxygen documentation will be auto-generated */
void ArrayOfLineMixingRecordReadAscii(// Generic Output:
                                      ArrayOfLineMixingRecord& line_mixing_records,
                                      // Generic Input:
                                      const String& filename,
                                      const Verbosity& verbosity)
{
    ifstream ifs;
    open_input_file(ifs, filename);

    line_mixing_records.resize(0);

    // Read the line mixing file
    Index linenr = 0;
    String line;
    istringstream is;
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
            String species_string;
            SpeciesTag line_species;

            is >> species_string;
            line_species = SpeciesTag(species_string);

            LineMixingRecord lmr(line_species.Species(), line_species.Isotopologue());

            Rational r;
            is >> r; lmr.Quantum().SetLower(QN_v1, r);
            lmr.Quantum().SetUpper(QN_v1, r);
            is >> r; lmr.Quantum().SetUpper(QN_N,  r);
            is >> r; lmr.Quantum().SetLower(QN_N,  r);
            is >> r; lmr.Quantum().SetUpper(QN_J,  r);
            is >> r; lmr.Quantum().SetLower(QN_J,  r);

            vector<Numeric> temp_mixing_data;
            String s;
            char *c;
            while (is.good() && is)
            {
                is >> s;
                s.trim();
                if (s.nelem())
                {
                    temp_mixing_data.push_back(strtod(s.c_str(), &c));
                    if (c != s.c_str() + s.nelem())
                        throw std::runtime_error(line);
                }
            }

            lmr.Data() = temp_mixing_data;
            line_mixing_records.push_back(lmr);
        } catch (runtime_error e) {
            ostringstream os;

            os << "Error parsing line mixing file in line " << linenr << std::endl;
            os << e.what();
            throw std::runtime_error(os.str());
        }
    }

    CREATE_OUT2;
    out2 << "  Read " << line_mixing_records.nelem() << " lines from " << filename << ".\n";
}

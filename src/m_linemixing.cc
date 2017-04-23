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
#include "complex.h"
#include "lin_alg.h"
#include "global_data.h"
#include "Faddeeva.hh"
#include "linescaling.h"
#include "jacobian.h"
#include <Eigen/Eigenvalues>


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


void SetBandIdentifiersAuto(ArrayOfQuantumIdentifier& band_identifiers,
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const Verbosity& verbosity)
{
  CREATE_OUT3;
  
  out3 << "Will set band_identifiers to as many bands as is included by default in this setup.\n"
           "To extend these bands, please edit this function. Currently, supported species and \n"
           "bands are only for CO2 and O2-66, though not all bands are guaranteed to be available.\n";
  
  band_identifiers.resize(0);
  
  static const SpeciesTag CO2_626 = SpeciesTag("CO2-626");
  static const SpeciesTag CO2_627 = SpeciesTag("CO2-627");
  static const SpeciesTag CO2_628 = SpeciesTag("CO2-628");
  static const SpeciesTag CO2_636 = SpeciesTag("CO2-636");
  static const SpeciesTag CO2_638 = SpeciesTag("CO2-638");
  static const SpeciesTag CO2_828 = SpeciesTag("CO2-828");
  static const SpeciesTag O2_66 = SpeciesTag("O2-66");
  bool o2_66_done, co2_626_done, co2_638_done, co2_636_done, co2_627_done, co2_628_done, co2_828_done;
  o2_66_done = co2_626_done = co2_638_done = co2_636_done = co2_627_done = co2_628_done = co2_828_done = false;
  
  //for( const ArrayOfSpeciesTag& spec : abs_species)
  for( Index ispec = 0; ispec < abs_species.nelem(); ispec++)
  {
    const ArrayOfSpeciesTag& spec = abs_species[ispec];
    //for( const SpeciesTag& tag : spec )
    for( Index iabsorber = 0; iabsorber < spec.nelem(); iabsorber++ )
    {
      const SpeciesTag& tag = spec[iabsorber];
      if(tag.LineMixing() == SpeciesTag::LINE_MIXING_ON)
      {
        if(CO2_626.Species() == tag.Species())
        {
          if( CO2_626.Isotopologue() == tag.Isotopologue()  and not co2_626_done)
          {
            co2_626_done = true;
            QuantumIdentifier qi;
            Index nbands=0;
            
            //Frequency is around 471 cm-1
            qi.SetFromStringForCO2Band("20003", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 479 cm-1
            qi.SetFromStringForCO2Band("13302", "12201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 544 cm-1
            qi.SetFromStringForCO2Band("11102", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 568 cm-1
            qi.SetFromStringForCO2Band("13302", "04401", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 594 cm-1
            qi.SetFromStringForCO2Band("20002", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 597 cm-1
            qi.SetFromStringForCO2Band("11102", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 608 cm-1
            qi.SetFromStringForCO2Band("10012", "01111", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 615 cm-1
            qi.SetFromStringForCO2Band("20003", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 618 cm-1
            qi.SetFromStringForCO2Band("10002", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 647 cm-1
            qi.SetFromStringForCO2Band("11102", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 654 cm-1
            qi.SetFromStringForCO2Band("01111", "00011", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 655 cm-1
            qi.SetFromStringForCO2Band("02211", "01111", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 667 cm-1
            qi.SetFromStringForCO2Band("01101", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 667 cm-1
            qi.SetFromStringForCO2Band("02201", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 668 cm-1
            qi.SetFromStringForCO2Band("03301", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 668 cm-1
            qi.SetFromStringForCO2Band("04401", "03301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 681 cm-1
            qi.SetFromStringForCO2Band("13301", "12201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 683 cm-1
            qi.SetFromStringForCO2Band("12201", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 688 cm-1
            qi.SetFromStringForCO2Band("11101", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 703 cm-1
            qi.SetFromStringForCO2Band("21101", "20001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 710 cm-1
            qi.SetFromStringForCO2Band("10011", "01111", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 720 cm-1
            qi.SetFromStringForCO2Band("20001", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 720 cm-1
            qi.SetFromStringForCO2Band("10001", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 738 cm-1
            qi.SetFromStringForCO2Band("20002", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 739 cm-1
            qi.SetFromStringForCO2Band("21101", "12201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 741 cm-1
            qi.SetFromStringForCO2Band("11101", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 757 cm-1
            qi.SetFromStringForCO2Band("12201", "03301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 770 cm-1
            qi.SetFromStringForCO2Band("13301", "04401", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 791 cm-1
            qi.SetFromStringForCO2Band("11101", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 828 cm-1
            qi.SetFromStringForCO2Band("12201", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 829 cm-1
            qi.SetFromStringForCO2Band("21101", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 864 cm-1
            qi.SetFromStringForCO2Band("20001", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 898 cm-1
            qi.SetFromStringForCO2Band("02211", "12201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 917 cm-1
            qi.SetFromStringForCO2Band("10011", "20001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 927 cm-1
            qi.SetFromStringForCO2Band("01111", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 941 cm-1
            qi.SetFromStringForCO2Band("10012", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 952 cm-1
            qi.SetFromStringForCO2Band("21101", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 960 cm-1
            qi.SetFromStringForCO2Band("00011", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1043 cm-1
            qi.SetFromStringForCO2Band("10011", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1063 cm-1
            qi.SetFromStringForCO2Band("00011", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1064 cm-1
            qi.SetFromStringForCO2Band("10012", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1071 cm-1
            qi.SetFromStringForCO2Band("01111", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1880 cm-1
            qi.SetFromStringForCO2Band("20003", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1905 cm-1
            qi.SetFromStringForCO2Band("13302", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1932 cm-1
            qi.SetFromStringForCO2Band("11102", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2003 cm-1
            qi.SetFromStringForCO2Band("20002", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2076 cm-1
            qi.SetFromStringForCO2Band("11101", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2093 cm-1
            qi.SetFromStringForCO2Band("12201", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2107 cm-1
            qi.SetFromStringForCO2Band("13301", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2112 cm-1
            qi.SetFromStringForCO2Band("21101", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2129 cm-1
            qi.SetFromStringForCO2Band("20001", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2165 cm-1
            qi.SetFromStringForCO2Band("21101", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2170 cm-1
            qi.SetFromStringForCO2Band("11112", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2180 cm-1
            qi.SetFromStringForCO2Band("20012", "20001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2182 cm-1
            qi.SetFromStringForCO2Band("20013", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2224 cm-1
            qi.SetFromStringForCO2Band("10012", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2288 cm-1
            qi.SetFromStringForCO2Band("13311", "13301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2290 cm-1
            qi.SetFromStringForCO2Band("13312", "13302", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2299 cm-1
            qi.SetFromStringForCO2Band("04411", "04401", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2299 cm-1
            qi.SetFromStringForCO2Band("02221", "02211", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2301 cm-1
            qi.SetFromStringForCO2Band("12211", "12201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2301 cm-1
            qi.SetFromStringForCO2Band("10021", "10011", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2302 cm-1
            qi.SetFromStringForCO2Band("10022", "10012", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2302 cm-1
            qi.SetFromStringForCO2Band("20011", "20001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2305 cm-1
            qi.SetFromStringForCO2Band("20013", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2306 cm-1
            qi.SetFromStringForCO2Band("20012", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2311 cm-1
            qi.SetFromStringForCO2Band("03311", "03301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2311 cm-1
            qi.SetFromStringForCO2Band("01121", "01111", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2313 cm-1
            qi.SetFromStringForCO2Band("11111", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2315 cm-1
            qi.SetFromStringForCO2Band("11112", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2324 cm-1
            qi.SetFromStringForCO2Band("02211", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2324 cm-1
            qi.SetFromStringForCO2Band("00021", "00011", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2326 cm-1
            qi.SetFromStringForCO2Band("10011", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2327 cm-1
            qi.SetFromStringForCO2Band("10012", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2336 cm-1
            qi.SetFromStringForCO2Band("01111", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2349 cm-1
            qi.SetFromStringForCO2Band("00011", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2428 cm-1
            qi.SetFromStringForCO2Band("20011", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2429 cm-1
            qi.SetFromStringForCO2Band("10011", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2429 cm-1
            qi.SetFromStringForCO2Band("20012", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2458 cm-1
            qi.SetFromStringForCO2Band("11111", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3465 cm-1
            qi.SetFromStringForCO2Band("20013", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3500 cm-1
            qi.SetFromStringForCO2Band("21101", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3527 cm-1
            qi.SetFromStringForCO2Band("30014", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3527 cm-1
            qi.SetFromStringForCO2Band("13312", "03301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3550 cm-1
            qi.SetFromStringForCO2Band("30012", "20001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3556 cm-1
            qi.SetFromStringForCO2Band("30013", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3566 cm-1
            qi.SetFromStringForCO2Band("10022", "00011", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3568 cm-1
            qi.SetFromStringForCO2Band("20013", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3580 cm-1
            qi.SetFromStringForCO2Band("11112", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3589 cm-1
            qi.SetFromStringForCO2Band("20012", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3612 cm-1
            qi.SetFromStringForCO2Band("10012", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3667 cm-1
            qi.SetFromStringForCO2Band("10021", "00011", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3676 cm-1
            qi.SetFromStringForCO2Band("30012", "20002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3679 cm-1
            qi.SetFromStringForCO2Band("30013", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3692 cm-1
            qi.SetFromStringForCO2Band("20012", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3705 cm-1
            qi.SetFromStringForCO2Band("30011", "20001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3711 cm-1
            qi.SetFromStringForCO2Band("20011", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3714 cm-1
            qi.SetFromStringForCO2Band("10011", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3723 cm-1
            qi.SetFromStringForCO2Band("11111", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3726 cm-1
            qi.SetFromStringForCO2Band("12211", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3727 cm-1
            qi.SetFromStringForCO2Band("13311", "03301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3799 cm-1
            qi.SetFromStringForCO2Band("30012", "20003", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3814 cm-1
            qi.SetFromStringForCO2Band("20011", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4005 cm-1
            qi.SetFromStringForCO2Band("00021", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4687 cm-1
            qi.SetFromStringForCO2Band("30014", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4786 cm-1
            qi.SetFromStringForCO2Band("31113", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4790 cm-1
            qi.SetFromStringForCO2Band("30014", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4839 cm-1
            qi.SetFromStringForCO2Band("30013", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4853 cm-1
            qi.SetFromStringForCO2Band("20013", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4931 cm-1
            qi.SetFromStringForCO2Band("31113", "11102", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4942 cm-1
            qi.SetFromStringForCO2Band("30013", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4946 cm-1
            qi.SetFromStringForCO2Band("31112", "11101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4959 cm-1
            qi.SetFromStringForCO2Band("30012", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4977 cm-1
            qi.SetFromStringForCO2Band("20012", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5062 cm-1
            qi.SetFromStringForCO2Band("30012", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5099 cm-1
            qi.SetFromStringForCO2Band("20011", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5114 cm-1
            qi.SetFromStringForCO2Band("30011", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5217 cm-1
            qi.SetFromStringForCO2Band("30011", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5247 cm-1
            qi.SetFromStringForCO2Band("10022", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5291 cm-1
            qi.SetFromStringForCO2Band("02221", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5315 cm-1
            qi.SetFromStringForCO2Band("01121", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5349 cm-1
            qi.SetFromStringForCO2Band("10012", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5584 cm-1
            qi.SetFromStringForCO2Band("00031", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5687 cm-1
            qi.SetFromStringForCO2Band("00031", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6075 cm-1
            qi.SetFromStringForCO2Band("30014", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6196 cm-1
            qi.SetFromStringForCO2Band("31113", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6227 cm-1
            qi.SetFromStringForCO2Band("30013", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6347 cm-1
            qi.SetFromStringForCO2Band("30012", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6356 cm-1
            qi.SetFromStringForCO2Band("31112", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6503 cm-1
            qi.SetFromStringForCO2Band("30011", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6860 cm-1
            qi.SetFromStringForCO2Band("03331", "03301", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6897 cm-1
            qi.SetFromStringForCO2Band("02231", "02201", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6905 cm-1
            qi.SetFromStringForCO2Band("10031", "10001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6907 cm-1
            qi.SetFromStringForCO2Band("10032", "10002", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6935 cm-1
            qi.SetFromStringForCO2Band("01131", "01101", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6972 cm-1
            qi.SetFromStringForCO2Band("00031", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 8192 cm-1
            qi.SetFromStringForCO2Band("10032", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 8293 cm-1
            qi.SetFromStringForCO2Band("10031", "00001", "626"); band_identifiers.push_back(qi); nbands++;
            
            out3 << "Set " << nbands << " for CO2-626";
          }
          else if( CO2_636.Isotopologue() == tag.Isotopologue() and not co2_636_done)
          {
            co2_636_done = true;
            QuantumIdentifier qi;
            Index nbands=0;
            
            //Frequency is around 617 cm-1
            qi.SetFromStringForCO2Band("10002", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 636 cm-1
            qi.SetFromStringForCO2Band("01111", "00011", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 648 cm-1
            qi.SetFromStringForCO2Band("01101", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 648 cm-1
            qi.SetFromStringForCO2Band("02201", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 649 cm-1
            qi.SetFromStringForCO2Band("03301", "02201", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 667 cm-1
            qi.SetFromStringForCO2Band("11101", "10001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 721 cm-1
            qi.SetFromStringForCO2Band("10001", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 739 cm-1
            qi.SetFromStringForCO2Band("11101", "02201", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 771 cm-1
            qi.SetFromStringForCO2Band("11101", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 883 cm-1
            qi.SetFromStringForCO2Band("01111", "11101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 913 cm-1
            qi.SetFromStringForCO2Band("00011", "10001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1017 cm-1
            qi.SetFromStringForCO2Band("00011", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2037 cm-1
            qi.SetFromStringForCO2Band("11101", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2157 cm-1
            qi.SetFromStringForCO2Band("10012", "10001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2248 cm-1
            qi.SetFromStringForCO2Band("01121", "01111", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2250 cm-1
            qi.SetFromStringForCO2Band("11111", "11101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2260 cm-1
            qi.SetFromStringForCO2Band("02211", "02201", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2260 cm-1
            qi.SetFromStringForCO2Band("00021", "00011", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2261 cm-1
            qi.SetFromStringForCO2Band("10012", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2262 cm-1
            qi.SetFromStringForCO2Band("10011", "10001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2271 cm-1
            qi.SetFromStringForCO2Band("01111", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2283 cm-1
            qi.SetFromStringForCO2Band("00011", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2367 cm-1
            qi.SetFromStringForCO2Band("10011", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2498 cm-1
            qi.SetFromStringForCO2Band("11112", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3450 cm-1
            qi.SetFromStringForCO2Band("13312", "03301", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3460 cm-1
            qi.SetFromStringForCO2Band("21113", "11102", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3473 cm-1
            qi.SetFromStringForCO2Band("12212", "02201", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3482 cm-1
            qi.SetFromStringForCO2Band("20013", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3482 cm-1
            qi.SetFromStringForCO2Band("21112", "11101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3497 cm-1
            qi.SetFromStringForCO2Band("30001", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3517 cm-1
            qi.SetFromStringForCO2Band("20012", "10001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3527 cm-1
            qi.SetFromStringForCO2Band("10012", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3621 cm-1
            qi.SetFromStringForCO2Band("20011", "10001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3621 cm-1
            qi.SetFromStringForCO2Band("20012", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3623 cm-1
            qi.SetFromStringForCO2Band("21112", "11102", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3625 cm-1
            qi.SetFromStringForCO2Band("21111", "11101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3632 cm-1
            qi.SetFromStringForCO2Band("10011", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3639 cm-1
            qi.SetFromStringForCO2Band("11111", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3641 cm-1
            qi.SetFromStringForCO2Band("12211", "02201", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3641 cm-1
            qi.SetFromStringForCO2Band("13311", "03301", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3725 cm-1
            qi.SetFromStringForCO2Band("20011", "10002", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4748 cm-1
            qi.SetFromStringForCO2Band("20013", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4887 cm-1
            qi.SetFromStringForCO2Band("20012", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 4991 cm-1
            qi.SetFromStringForCO2Band("20011", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5013 cm-1
            qi.SetFromStringForCO2Band("21111", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 5168 cm-1
            qi.SetFromStringForCO2Band("01121", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6745 cm-1
            qi.SetFromStringForCO2Band("01131", "01101", "636"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 6780 cm-1
            qi.SetFromStringForCO2Band("00031", "00001", "636"); band_identifiers.push_back(qi); nbands++;
            
            out3 << "Set " << nbands << " for CO2-636";
          }
          else if( CO2_628.Isotopologue() == tag.Isotopologue() and not co2_628_done)
          {
            co2_628_done = true;
            QuantumIdentifier qi;
            Index nbands=0;
            
            //Frequency is around 597 cm-1
            qi.SetFromStringForCO2Band("10002", "01101", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 662 cm-1
            qi.SetFromStringForCO2Band("01101", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 662 cm-1
            qi.SetFromStringForCO2Band("02201", "01101", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 703 cm-1
            qi.SetFromStringForCO2Band("10001", "01101", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 966 cm-1
            qi.SetFromStringForCO2Band("00011", "10001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1072 cm-1
            qi.SetFromStringForCO2Band("00011", "10002", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1259 cm-1
            qi.SetFromStringForCO2Band("10002", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1365 cm-1
            qi.SetFromStringForCO2Band("10001", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2205 cm-1
            qi.SetFromStringForCO2Band("10012", "10001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2307 cm-1
            qi.SetFromStringForCO2Band("02211", "02201", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2309 cm-1
            qi.SetFromStringForCO2Band("10011", "10001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2311 cm-1
            qi.SetFromStringForCO2Band("10012", "10002", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2319 cm-1
            qi.SetFromStringForCO2Band("01111", "01101", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2332 cm-1
            qi.SetFromStringForCO2Band("00011", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2415 cm-1
            qi.SetFromStringForCO2Band("10011", "10002", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2614 cm-1
            qi.SetFromStringForCO2Band("20002", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2618 cm-1
            qi.SetFromStringForCO2Band("21102", "01101", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2757 cm-1
            qi.SetFromStringForCO2Band("20001", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3571 cm-1
            qi.SetFromStringForCO2Band("10012", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3675 cm-1
            qi.SetFromStringForCO2Band("10011", "00001", "628"); band_identifiers.push_back(qi); nbands++;
            
            out3 << "Set " << nbands << " for CO2-628";
          }
          else if( CO2_828.Isotopologue() == tag.Isotopologue() and not co2_828_done)
          {
            co2_828_done = true;
            QuantumIdentifier qi;
            Index nbands=0;
            
            //Frequency is around 2301 cm-1
            qi.SetFromStringForCO2Band("01111", "01101", "828"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2314 cm-1
            qi.SetFromStringForCO2Band("00011", "00001", "828"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3525 cm-1
            qi.SetFromStringForCO2Band("10012", "00001", "828"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3638 cm-1
            qi.SetFromStringForCO2Band("10011", "00001", "828"); band_identifiers.push_back(qi); nbands++;
            
            out3 << "Set " << nbands << " for CO2-828";
          }
          else if( CO2_627.Isotopologue() == tag.Isotopologue() and not co2_627_done)
          {
            co2_627_done = true;
            QuantumIdentifier qi;
            Index nbands=0;
            
            //Frequency is around 664 cm-1
            qi.SetFromStringForCO2Band("01101", "00001", "627"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 711 cm-1
            qi.SetFromStringForCO2Band("10001", "01101", "627"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 963 cm-1
            qi.SetFromStringForCO2Band("00011", "10001", "627"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 1376 cm-1
            qi.SetFromStringForCO2Band("10001", "00001", "627"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2340 cm-1
            qi.SetFromStringForCO2Band("00011", "00001", "627"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 2641 cm-1
            qi.SetFromStringForCO2Band("20002", "00001", "627"); band_identifiers.push_back(qi); nbands++;
            
            out3 << "Set " << nbands << " for CO2-627";
          }
          else if( CO2_638.Isotopologue() == tag.Isotopologue() and not co2_638_done)
          {
            co2_638_done = true;
            QuantumIdentifier qi;
            Index nbands=0;
            
            //Frequency is around 2265 cm-1
            qi.SetFromStringForCO2Band("00011", "00001", "638"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3490 cm-1
            qi.SetFromStringForCO2Band("10012", "00001", "638"); band_identifiers.push_back(qi); nbands++;
            
            //Frequency is around 3587 cm-1
            qi.SetFromStringForCO2Band("10011", "00001", "638"); band_identifiers.push_back(qi); nbands++;
            
            out3 << "Set " << nbands << " for CO2-638";
          }
        }
        else if(O2_66.Species() == tag.Species() and O2_66.Isotopologue() == tag.Isotopologue() and not o2_66_done)
        {
          o2_66_done = true;
          QuantumIdentifier qi;
          // The main band in the 60 GHz to 1.5 THz range.  LM mostly at 60 GHz, though the rest falls into this band.
          qi.SetFromString("O2-66 TR UP v1 0/1 LO v1 0/1"); band_identifiers.push_back(qi);
          // The secondary band in the 60 GHz to 1.5 THz range.  LM mostly at 60 GHz, though the rest falls into this band.
          qi.SetFromString("O2-66 TR UP v1 1/1 LO v1 1/1"); band_identifiers.push_back(qi);
        }
        else 
        {
          throw std::runtime_error("Unsupported species (edit m_linemixing.cc to fix or manually enter band_identifiers).\n");
        }
      }
    }
  }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_bandFromband_identifiers( ArrayOfArrayOfLineRecord&       abs_lines_per_band,
                                             ArrayOfArrayOfSpeciesTag&       abs_species_per_band,
                                             ArrayOfArrayOfLineRecord&       abs_lines_per_species,
                                             const ArrayOfArrayOfSpeciesTag& abs_species,
                                             const ArrayOfQuantumIdentifier& band_identifiers,
                                             const Verbosity&                verbosity)
{
  CREATE_OUT3;
  out3<<"Sets line mixing tag for provided bands. " <<
  "Requires \"*-LM-*\" tag in abs_species.\n";
  
  if(abs_lines_per_species.nelem() not_eq abs_species.nelem())
      throw std::runtime_error("Mismatching abs_species and abs_lines_per_species");
  
  abs_lines_per_band.resize(band_identifiers.nelem());
  abs_species_per_band.resize(band_identifiers.nelem());
  
  // This is a preallocated finding of a band
  LineMixingData lmd_byband;
  lmd_byband.SetByBand(); 
  
  #pragma omp parallel for        \
  if (!arts_omp_in_parallel())    
  for (Index qi = 0; qi < band_identifiers.nelem(); qi++)
  {
    const QuantumIdentifier& band_id = band_identifiers[qi];
    
    // Two variables that are used inside the loop
    ArrayOfIndex matches;
    ArrayOfQuantumMatchInfo match_info;
    
    for (Index s = 0; s < abs_lines_per_species.nelem(); s++)
    {
      // Skip this species if qi is not part of the species represented by this abs_lines
      if(abs_species[s][0].Species() not_eq band_id.Species() or 
         abs_species[s][0].LineMixing() == SpeciesTag::LINE_MIXING_OFF)
        continue;
      
      ArrayOfLineRecord& species_lines = abs_lines_per_species[s];
      
      // Copy this
      abs_species_per_band[qi] = abs_species[s];
      
      // Run internal mathcing routine
      match_lines_by_quantum_identifier(matches, match_info, band_id, species_lines);
      
      // Use info about mathced lines to tag the relevant parameter
      for (Index i = 0; i < matches.nelem(); i++)
      {
        QuantumMatchInfo& qm = match_info[i];
        
        LineRecord& lr = species_lines[matches[i]];
        
        // If any of the levels match partially or fully set the right quantum number
        if(qm.Upper()==QMI_NONE or qm.Lower()==QMI_NONE)
        {
          continue;
        }
        else // we will accept this match if both levels are at least partially matched
        {
          abs_lines_per_band[qi].push_back(lr);
          lr.SetLineMixingData(lmd_byband);
        }
      }
    }
  }
  
  for (Index qi = 0; qi < band_identifiers.nelem(); qi++)
  {
    Index nlines = abs_lines_per_band[qi].nelem();
    out3 << "Found " << nlines << " lines of the band: " << band_identifiers[qi] << "\n";
    if(nlines)
      out3 << "\tfrequency range: " << round(abs_lines_per_band[qi][0].F()/1e9) << " to " 
                                   << round(abs_lines_per_band[qi][nlines-1].F()/1e9) << " GHz\n";
    else 
      out3 << "\tfrequency range: NA\n";
  }
  
}

void calculate_xsec_from_relmat(ArrayOfMatrix& xsec,
                                ArrayOfArrayOfMatrix& dxsec_dx,
                                const PropmatPartialsData& ppd,
                                const ConstMatrixView Wmat,
                                const ConstMatrixView Wmat_dt,
                                const ConstVectorView f0,
                                const ConstVectorView f_grid,
                                const ConstVectorView d0,
                                const ConstVectorView d0_dt,
                                const ConstVectorView rhoT,
                                const ConstVectorView rhoT_dt,
                                const ConstVectorView psf,
                                const ConstVectorView psf_dt,
                                const Numeric&  T,
                                const Numeric&  isotopologue_mass,
                                const Numeric&  isotopologue_ratio,
                                const Index& this_species,
                                const Index& this_level,
                                const Index&    n )
{ 
    // Physical constants
    extern const Numeric BOLTZMAN_CONST;
    extern const Numeric AVOGADROS_NUMB;
    extern const Numeric SPEED_OF_LIGHT;
    extern const Numeric PI;
    extern const Numeric PLANCK_CONST;
    
    // internal constant
    const Index nf  = f_grid.nelem();
    const Numeric kT = BOLTZMAN_CONST * T;
    const Numeric doppler_const = sqrt( 2.0 * BOLTZMAN_CONST * AVOGADROS_NUMB / isotopologue_mass ) / SPEED_OF_LIGHT;
    static const Numeric invSqrtPI = 1.0/sqrt(PI);
    
    // Setup for the matrix to be diagonalized
    ComplexMatrix F0plusiPW(n, n), F0plusiPW_dt;
    {
      if(ppd.do_temperature())
      {
        F0plusiPW_dt.resize(n, n);
        
        for(Index i = 0; i < n; i++)
        {
          for(Index j = 0; j < n; j++)
          {
            if(j == i)
            {
              F0plusiPW(i, i) = Complex(f0[i], Wmat(i, i));
              F0plusiPW_dt(i, i) = Complex(f0[i], Wmat_dt(i, i));
            }
            else
            {
              F0plusiPW(i, j) = Complex(0.0, Wmat(i, j));
              F0plusiPW_dt(i, j) = Complex(0.0, Wmat_dt(i, j));
            }
          }
        }
      }
      else 
      {
        for(Index i = 0; i < n; i++)
        {
          for(Index j = 0; j < n; j++)
          {
            if(j == i)
            {
              F0plusiPW(i, i) = Complex(f0[i], Wmat(i, i));
            }
            else
            {
              F0plusiPW(i, j) = Complex(0.0, Wmat(i, j));
            }
          }
        }
      }
    }
    
    // z_eigs is set to contain the eigenvalues, 
    // E is the matrix of eigenvectors and invE is its inverse
    // the idea is that E*diag(z_eig)*invE is F0plusiPW
    ComplexVector z_eigs(n),z_eigs_dt;
    ComplexMatrix E(n, n), invE(n, n), E_dt, invE_dt;
    {
      diagonalize(E, z_eigs, F0plusiPW);
      inv(invE, E);
      if(ppd.do_temperature())
      {
        z_eigs_dt.resize(n);
        E_dt.resize(n, n);
        invE_dt.resize(n, n);
        diagonalize(E_dt, z_eigs_dt, F0plusiPW_dt);
        inv(invE_dt, E_dt);
      }
    }
    
    // Conjunction of eigenvalues and the Doppler broadening
    Vector sigma(n), sigma_dt;
    {
      if(ppd.do_temperature())
      {
        //dsigma_dt.resize(n);
        sigma_dt.resize(n);
        for(Index if1 = 0; if1 < n; if1++)
        {
          sigma[if1] = f0[if1] * doppler_const * sqrt(T);
          sigma_dt[if1] = f0[if1] * doppler_const * sqrt(T + ppd.Temperature_Perturbation()) ;
          
          z_eigs[if1] = conj(z_eigs[if1]); 
          z_eigs_dt[if1] = conj(z_eigs_dt[if1]);
        }
      }
      else
      {
        for(Index if1 = 0; if1 < n; if1++)
        {
          sigma[if1]= f0[if1] * doppler_const * sqrt(T);
          z_eigs[if1] = conj(z_eigs[if1]); 
        }
      }
    }
    
    // Equivalent line strength
    ComplexVector equivS0(n, 0.0), equivS0_dt;
    {
      if(ppd.do_temperature())
      {
        equivS0_dt.resize(n);
        for(Index k = 0; k < n; k++)
        {
          Complex z1 = 0.0;
          Complex z2 = 0.0;
          Complex z1_dt = 0.0;
          Complex z2_dt = 0.0;
          for(Index j = 0; j < n; j++)
          {
            z1 += d0[j] * E(j, k);
            z2 += rhoT[j] * d0[j] * invE(k, j);
            z1_dt += d0_dt[j] * E_dt(j, k);
            z2_dt += rhoT_dt[j] * d0_dt[j] * invE_dt(k, j);
            
          }
          equivS0[k] = conj(z1 * z2) * isotopologue_ratio;
          equivS0_dt[k] = conj(z1_dt * z2_dt) * isotopologue_ratio;
        }
      }
      else
      {
        for(Index k = 0; k < n; k++)
        {
          Complex z1 = 0.0;
          Complex z2 = 0.0;
          for(Index j = 0; j < n; j++)
          {
            z1 += d0[j] * E(j, k);
            z2 += rhoT[j] * d0[j] * invE(k, j);
            
          }
          equivS0[k] = conj(z1 * z2) * isotopologue_ratio;
        }
      }
    }
    
    // Perform the computations to get at the cross-sections
    for(Index iv = 0; iv < nf; iv++)
    {
      const Numeric s0_freqfac =  f_grid[iv] * (1 - exp(- PLANCK_CONST * f_grid[iv] / kT)); // Adapted from Niro's code --- why are they changing from f0 to f?
      Numeric s0_freqfac_dt = 0.0;
      if(ppd.do_temperature())
      {
        s0_freqfac_dt = f_grid[iv] * (1 - exp(- PLANCK_CONST * f_grid[iv] / (BOLTZMAN_CONST * (T+ppd.Temperature_Perturbation()))));
      }
      for(Index if0 = 0; if0 < n; if0++)
      {
        const Numeric ls_normfac = invSqrtPI / sigma[if0];
        const Complex z = (f_grid[iv] -(z_eigs[if0] + psf[if0])) / sigma[if0];
        const Complex w = Faddeeva::w(z) * ls_normfac;
        const Numeric xsec_this = (equivS0[if0] * w).real() * s0_freqfac;
        
        #pragma omp atomic
        xsec[this_species](iv, this_level) += xsec_this;
        
        Complex dw_df = Complex(0.0, 0.0);
        if(ppd.do_frequency())
          dw_df = 2 * (z * w - ls_normfac);
        
        for(Index iq = 0; iq < ppd.nelem(); iq++)
        {
          if(ppd(iq) == JQT_frequency  or ppd(iq) == JQT_wind_magnitude  or 
             ppd(iq) == JQT_wind_u  or ppd(iq) == JQT_wind_v  or ppd(iq) == JQT_wind_w)
          {
            #pragma omp atomic
            dxsec_dx[this_species][iq](iv, this_level) += real((dw_df + 1.0/f_grid[iv] * w) * s0_freqfac * equivS0[if0]);
          }
          else if(ppd(iq) == JQT_temperature)
          {
            const Numeric ls_normfac_dt = invSqrtPI / sigma_dt[if0];
            const Complex z_dt = (f_grid[iv] -(z_eigs_dt[if0] + psf_dt[if0])) / sigma_dt[if0];
            const Complex w_dt = Faddeeva::w(z_dt) * ls_normfac_dt;
            const Numeric xsec_this_dt = (equivS0_dt[if0] * w_dt).real() * s0_freqfac_dt;
            #pragma omp atomic
            dxsec_dx[this_species][iq](iv, this_level) += (xsec_this_dt - xsec_this) / ppd.Temperature_Perturbation();
          }
        }
      }
    }
}


static const Index hartman_tran_type = 0;
static const Index linear_type = 1;


void SetRelaxationMatrixCalcType( ArrayOfIndex& relmat_type_per_band,
                                  const ArrayOfArrayOfLineRecord&  abs_lines_per_band,
                                  const Index& type,
                                  const Verbosity& verbosity)
{
  CREATE_OUT2;
  if(type not_eq hartman_tran_type and type not_eq linear_type)
    throw std::runtime_error("Not a supported type.  Check documentation for supported types");
  
  const Index& n = abs_lines_per_band.nelem();
  Index i;
  relmat_type_per_band.resize(n);
  for(i=0; i<n; i++)
    relmat_type_per_band[i] = type;
  
  if(type == hartman_tran_type)
    out2 << "Hartmann-Tran type line mixing selected for all lines\n";
  else if(type == linear_type)
    out2 << "Mendaza linear type line mixing selected for all lines\n";
  else
    out2 << "Unknown type line mixing selected  for all lines (please fix by adding type)\n";
}


void SetRelaxationMatrixCalcType( ArrayOfIndex& relmat_type_per_band,
                                  const ArrayOfArrayOfLineRecord&  abs_lines_per_band,
                                  const ArrayOfIndex& type,
                                  const Verbosity& verbosity)
{
  CREATE_OUT2;
  const Index& n = type.nelem();
  
  if(n == 1)
    SetRelaxationMatrixCalcType(relmat_type_per_band, abs_lines_per_band, type[0], verbosity);
  else if(n not_eq type.nelem())
    throw std::runtime_error("Mismatching length of type and abs_lines_per_band");
  else
  {
    Index i;
    relmat_type_per_band.resize(n);
    for(i=0; i<n; i++)
    {
      if(type[i] not_eq hartman_tran_type or type[i] not_eq linear_type)
        throw std::runtime_error("Not a supported type.  Check documentation for supported types");
      relmat_type_per_band[i] = type[i];
    }
    
    out2 << "Mix of line mixing types set for each band\n";
    
  }
}

#ifdef ENABLE_RELMAT
extern "C"
{
    // This is the interfaces between the Fortran code that calculates W and ARTS
    extern void arts_relmat_interface__hartmann_and_niro_type(
        long   *nlines,
        double *fmin,
        double *fmax,
        long   *M,
        long   *I,
        double *v,
        double *S,
        double *gamma_air,
        double *E_double_prime,
        double *n_air,
        long   *upper,
        long   *lower,
        long   *g_prime,
        long   *g_double_prime,
        double *temperature,
        double *pressure,
        double *partition_function_t,
        double *partition_function_t0,
        double *isotopologue_mass,
        long   *number_of_perturbers,
        long   *molecule_code_perturber,
        long   *iso_code_perturber,
        double *perturber_mass,
        double *vmr,
        //output+input
        long   *debug_in__error_out,
        //outputs
        double *W,
        double *dipole,
        double *rhoT
    );
    
    extern void arts_relmat_interface__linear_type(
      long   *nlines,
      double *fmin,
      double *fmax,
      long   *M,
      long   *I,
      double *v,
      double *S,
      double *gamma_air,
      double *E_double_prime,
      double *n_air,
      long   *upper,
      long   *lower,
      long   *g_prime,
      long   *g_double_prime,
      double *temperature,
      double *pressure,
      double *partition_function_t,
      double *partition_function_t0,
      double *isotopologue_mass,
      long   *number_of_perturbers,
      long   *molecule_code_perturber,
      long   *iso_code_perturber,
      double *perturber_mass,
      double *vmr,
      //output+input
      long   *debug_in__error_out,
      //outputs
      double *W,
      double *dipole,
      double *rhoT
    );
}


void abs_xsec_per_speciesAddLineMixedBands( // WS Output:
                                            ArrayOfMatrix&                   abs_xsec_per_species,
                                            ArrayOfArrayOfMatrix&            dabs_xsec_per_species_dx,
                                            ArrayOfArrayOfMatrix&            relmat_per_band,
                                            // WS Input:                     
                                            const ArrayOfArrayOfLineRecord&  abs_lines_per_band,
                                            const ArrayOfArrayOfSpeciesTag&  abs_species_per_band,
                                            const ArrayOfArrayOfSpeciesTag&  abs_species,
                                            const SpeciesAuxData&            isotopologue_ratios,
                                            const SpeciesAuxData&            partition_functions,
                                            const ArrayOfRetrievalQuantity&  jacobian_quantities,
                                            const Vector&                    f_grid,
                                            const Vector&                    abs_p,
                                            const Vector&                    abs_t,
                                            const Numeric&                   lm_p_lim,
                                            const ArrayOfIndex&              relmat_type_per_band,
                                            const Index&                     write_relmat_per_band,
                                            const Index&                     debug,
                                            const Verbosity&)
{
  using global_data::species_data;
  using global_data::SpeciesMap;
  
  // Physical constants
  extern const Numeric SPEED_OF_LIGHT;
  extern const Numeric BOLTZMAN_CONST;
  extern const Numeric AVOGADROS_NUMB;
  extern const Numeric PI;
  extern const Numeric ATM2PA;
  
  // HITRAN to ARTS constants
  const Numeric doppler_const = sqrt( 2.0 * BOLTZMAN_CONST * AVOGADROS_NUMB ) / SPEED_OF_LIGHT;
  static const Numeric w2Hz               = SPEED_OF_LIGHT *1E2;
  static const Numeric lower_energy_const = wavenumber_to_joule(1.0);
  static const Numeric I0_hi2arts         = 1E-2 * SPEED_OF_LIGHT;
  static const Numeric gamma_hi2arts      = w2Hz / ATM2PA;
  static const Numeric invSqrtPI = 1.0/sqrt(PI);
  
  // size of atmosphere and input/output
  const Index nps         = abs_p.nelem();
  const Index nts         = abs_t.nelem();
  const Index nf          = f_grid.nelem();
  const Index nspecies    = abs_species.nelem();
  
  // Relmat constants
  const Numeric relmat_T0 = 296.0;
  
  // These should be identical
  if(nps not_eq nts)
    throw std::runtime_error("Different lengths of atmospheric input than expected.");
  
  if(nspecies == 0 or nspecies not_eq abs_xsec_per_species.nelem())
    throw std::runtime_error("Absorption species and abs_xsec_per_species are not from same source.");
  else
    for(Index i=0; i < nspecies; i++)
    {
      if(abs_xsec_per_species[i].ncols() not_eq nts)
      {
        throw std::runtime_error("Unexpected size of xsec matrix not matching abs_t and abs_p length.");
      }
      else if(abs_xsec_per_species[i].nrows() not_eq nf)
      {
        throw std::runtime_error("Unexpected size of xsec matrix not matching f_grid length.");
      }
      if(relmat_type_per_band.nelem() not_eq abs_lines_per_band.nelem())
      {
        throw std::runtime_error("Mismatching relmat_type_per_band and abs_lines_per_band.\n");
      }
    }
      
      PropmatPartialsData ppd(jacobian_quantities);
    ppd.supportsRelaxationMatrix();
  // It should support temperature and wind, and maybe even the line mixing parameters as an error vector
  
  const Index nbands = abs_lines_per_band.nelem();
  
  if(write_relmat_per_band)
  { 
    relmat_per_band.resize(nps);
    for(Index ip=0; ip<nps; ip++)
    {
      relmat_per_band[ip].resize(nbands);
    }
  }
  
  if(nbands not_eq abs_species_per_band.nelem())
    throw std::runtime_error("Error in definition of the bands.  Mismatching length of *_per_band arrays.");
  
  // Make constant input not constant and convert to wavenumber
  Vector v(nf);
  for(Index ifs=0; ifs<nf; ifs++)
    v[ifs] = f_grid[ifs] / w2Hz; 
  
  // Setting up thermal bath:  only in simplistic air for now
  // This means: 21% O2 and 79% N2
  //
  long    number_of_perturbers = 2;
  long   *molecule_code_perturber = new long[number_of_perturbers];
  long   *iso_code_perturber = new long[number_of_perturbers];
  double *perturber_mass = new double[number_of_perturbers];
  Vector  vmr(number_of_perturbers);
  
  // Setup air as background gas for now...
  { 
    vmr[0] = 0.2095;
    const SpeciesRecord& O2 = species_data[SpeciesMap.find("O2")->second];
    const IsotopologueRecord& O2_66 = O2.Isotopologue()[0];
    const Index O2_66_hitran_tag = O2_66.HitranTag();
    iso_code_perturber[0] = O2_66_hitran_tag%10;
    molecule_code_perturber[0] = (O2_66_hitran_tag)/10;
    perturber_mass[0] = O2_66.Mass();
    
    vmr[1] = 1.0-0.2095;
    const SpeciesRecord& N2 = species_data[SpeciesMap.find("N2")->second];
    const IsotopologueRecord& N2_44 = N2.Isotopologue()[0];
    const Index N2_44_hitran_tag = N2_44.HitranTag();
    iso_code_perturber[1] = N2_44_hitran_tag%10;
    molecule_code_perturber[1] = (N2_44_hitran_tag)/10;
    perturber_mass[1] = N2_44.Mass();
  }
  
  #pragma omp parallel for                    \
  if (!arts_omp_in_parallel() && nbands not_eq 1)
  for(Index iband=0;iband<nbands;iband++)
  {
    // band pointer
    const ArrayOfLineRecord& this_band = abs_lines_per_band[iband];
    
    long nlines = (long) this_band.nelem();
    
    // Worth doing anything?
    if(nlines==0) { continue; }
    
    // Send in frequency range
    Numeric fmin, fmax;
    
    // To store the xsec matrix we need to know where
    Index this_species=-1;
    
    // Allocation of band specific inputs (types: to be used as fortran input)
    long   M = (this_band[0].IsotopologueData().HitranTag())/10;
    long   I = (this_band[0].IsotopologueData().HitranTag())%10;
    long   *upper          = new long[4*nlines];
    long   *lower          = new long[4*nlines];
    long   *g_prime        = new long[nlines];
    long   *g_double_prime = new long[nlines];
    Vector v0(nlines);
    Vector S(nlines);
    Vector gamma_air(nlines);
    Vector delta_air(nlines);
    Vector E_double_prime(nlines);
    Vector n_air(nlines);
    Numeric mass, iso_ratio;
    
    for( long iline=0; iline<nlines; iline++ )
    {
      // Line data
      const LineRecord& this_line = this_band[iline];
      const IsotopologueRecord& this_iso = this_line.IsotopologueData();
      
      // For first line do something special with mass and name and such
      if( iline==0 )
      {
        // Isotopologue values
        mass  = this_iso.Mass();
        String iso_name = this_iso.Name();
        int isona;
        extract(isona,iso_name,iso_name.nelem());
        
        const ArrayOfSpeciesTag& band_tags=abs_species_per_band[iband];
        
        // Finds the first species in abs_species_per_band that matches abs_species
        bool this_one=false;
        for(Index ispecies=0; ispecies<nspecies; ispecies++)
        {
          const ArrayOfSpeciesTag& species_tags=abs_species[ispecies];
          const Index nbandtags = band_tags.nelem();
          const Index nspeciestags = species_tags.nelem();
          
          // Test if there is
          if(nbandtags not_eq nspeciestags) { break; }
          for(Index itags=0; itags<nspeciestags; itags++)
          {
            if(band_tags[itags]==species_tags[itags])
            {
              this_one = true;
              break;
            }
          }
          
          if(this_one)
          {
            this_species = ispecies;
            break;
          }
          
        }
        
        if(!this_one)
          throw std::runtime_error("abs_species and abs_species_per_band disagrees"
          " on absorption species");
        
        iso_ratio =  isotopologue_ratios.getParam(abs_lines_per_band[iband][iline].Species(),
                                                  abs_lines_per_band[iband][iline].Isotopologue())[0].data[0];
        
      }
      else 
      {
        if( mass not_eq this_iso.Mass() )
        {
          throw std::runtime_error("There are lines of different Isotopologues"
          " in abs_lines_per_band,");
        }
      }
      
      // Pressure broadening at relmat temperatures
      if(not this_line.PressureBroadening().isAirBroadening())
      {
        std::ostringstream os;
        os << "Line is not air broadening type but only air broadening types are suported.\n";
        os << "Its type is " << this_line.PressureBroadening().Type2StorageTag();
        throw std::runtime_error(os.str());
      }
      
      // Ensure that temperatures are sufficiently close
      if( 1e-4 < abs(this_line.Ti0()-relmat_T0))
      {
        std::ostringstream os;
        os << "Line is not of same standard temperature as relmat is expecting.\n";
        os << "Expecting: " << relmat_T0 <<" K.  Getting: " << this_line.Ti0() << " K.";
        throw std::runtime_error(os.str());
      }
      
      gamma_air[iline] = this_line.PressureBroadening().AirBroadeningAgam() / gamma_hi2arts;
      delta_air[iline] = this_line.PressureBroadening().AirBroadeningPsf();
      n_air[iline] = this_line.Nair();
      
      // Line information converted to relmat format --- i.e. to HITRAN format
      v0[iline]             = this_line.F()/w2Hz; 
      S[iline]              = this_line.I0()/I0_hi2arts;
      E_double_prime[iline] = this_line.Elow()/lower_energy_const;
      g_prime[iline]        = (long) this_line.G_upper(); // NB:  Numeric to long... why?
      g_double_prime[iline] = (long) this_line.G_lower(); // NB:  Numeric to long... why?
      
      // Quantum numbers converted to relmat format, again Numeric/Rational to long... why?
      Rational a;
      
      // l2 is for molecules like CO2
      {
        a = this_line.QuantumNumbers().Lower()[QN_l2];
        a.Simplify();
        if(a.isUndefined())
          lower[0+4*iline] = -1;
        else  if(a.Denom()==1)
          lower[0+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in l2.");
        a = this_line.QuantumNumbers().Upper()[QN_l2];
        a.Simplify();
        if(a.isUndefined())
          upper[0+4*iline] = -1;
        else  if(a.Denom()==1)
          upper[0+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in l2.");
      }
      
      // J is universally important for linear molecules
      {
        a = this_line.QuantumNumbers().Lower()[QN_J];
        a.Simplify();
        if(a.isUndefined())
          lower[1+4*iline] = -1;
        else  if(a.Denom()==1)
          lower[1+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in J.");
        a = this_line.QuantumNumbers().Upper()[QN_J];
        a.Simplify();
        if(a.isUndefined())
          upper[1+4*iline] = -1;
        else  if(a.Denom()==1)
          upper[1+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in J.");
      }
      
      // N is important for molecules with magnetic dipoles
      {
        a = this_line.QuantumNumbers().Lower()[QN_N];
        a.Simplify();
        if(a.isUndefined())
          lower[2+4*iline] = -1;
        else  if(a.Denom()==1)
          lower[2+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in N.");
        a = this_line.QuantumNumbers().Upper()[QN_N];
        a.Simplify();
        if(a.isUndefined())
          upper[2+4*iline] = -1;
        else  if(a.Denom()==1)
          upper[2+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in N.");
      }
      
      // S is important for molecules with magnetic dipoles
      {
        a = this_line.QuantumNumbers().Lower()[QN_S];
        a.Simplify();
        if(a.isUndefined())
          lower[3+4*iline] = -1;
        else  if(a.Denom()==1)
          lower[3+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in S.");
        a = this_line.QuantumNumbers().Upper()[QN_S];
        a.Simplify();
        if(a.isUndefined())
          upper[3+4*iline] = -1;
        else  if(a.Denom()==1)
          upper[3+4*iline] = (long) a.toIndex();
        else 
          throw std::runtime_error("Half quantum numbers not supported in S.");
      }
      
      // Set fmax and fmin.  Why do I need this again?
      if(iline==0)
      {
        fmin = v0[0]-1.0;
        fmax = v0[0]+1.0;
      }
      else 
      {
        if(fmin > v0[iline])
          fmin = v0[iline]-1.0;
        if(fmax < v0[iline])
          fmax = v0[iline]+1.0;
      }
      
    }
    
    Vector f0 = v0;
    f0 *= w2Hz;
    
    #pragma omp parallel for                    \
    if (!arts_omp_in_parallel() && nps not_eq 1)
      for(Index ip=0; ip<nps; ip++)
      {
        
      // Information on the lines will be here after relmat is done
      Matrix W(nlines,nlines);
      Vector dipole(nlines);
      Vector rhoT(nlines);
      
      Matrix W_dt;
      Vector dipole_dt;
      Vector rhoT_dt;
      if(ppd.do_temperature())
      {
        W_dt.resize(nlines,nlines);
        dipole_dt.resize(nlines);
        rhoT_dt.resize(nlines);
      }
      
      Vector psf(nlines), psf_dt;
      if(ppd.do_temperature())
        psf_dt.resize(nlines);
      
      if(abs_p[ip] > lm_p_lim)
      {
        
        // Find the partition function
        Numeric QT0;
        Numeric QT;
        
        // Get partition function information
        partition_function( QT0,
                            QT,
                            abs_lines_per_band[iband][0].Ti0(),
                            abs_t[ip],
                            partition_functions.getParamType(abs_lines_per_band[iband][0].Species(), 
                                                             abs_lines_per_band[iband][0].Isotopologue()),
                            partition_functions.getParam(abs_lines_per_band[iband][0].Species(), 
                                                         abs_lines_per_band[iband][0].Isotopologue()),
                            false);
        
        // Cannot be constants for Fortran's sake
        Numeric t;
        t = abs_t[ip];
        Numeric p;
        p = abs_p[ip]/ATM2PA; // HITRAN pressure unit is in atmospheres
        
        // Error handling variable is one or zero initially.  If 
        Index error_handling = debug?1:0;
        
        // Calling Teresa's code
        if(relmat_type_per_band[iband] == hartman_tran_type)
        {
          arts_relmat_interface__hartmann_and_niro_type(
            &nlines, &fmin, &fmax,
            &M, &I, v0.get_c_array(), S.get_c_array(),
            gamma_air.get_c_array(),E_double_prime.get_c_array(),n_air.get_c_array(),
            upper, lower,
            g_prime, g_double_prime,
            &t, &p, &QT, &QT0, &mass,
            &number_of_perturbers, molecule_code_perturber, 
            iso_code_perturber, perturber_mass, vmr.get_c_array(), &error_handling,
            W.get_c_array(), dipole.get_c_array(), rhoT.get_c_array() );
        }
        else if(relmat_type_per_band[iband] == linear_type)
        {
          arts_relmat_interface__linear_type(
            &nlines, &fmin, &fmax,
            &M, &I, v0.get_c_array(), S.get_c_array(),
            gamma_air.get_c_array(),E_double_prime.get_c_array(),n_air.get_c_array(),
            upper, lower,
            g_prime, g_double_prime,
            &t, &p, &QT, &QT0, &mass,
            &number_of_perturbers, molecule_code_perturber, 
            iso_code_perturber, perturber_mass, vmr.get_c_array(), &error_handling,
            W.get_c_array(), dipole.get_c_array(), rhoT.get_c_array() );
        }
        else 
        {
          throw std::runtime_error("Unsupported relaxation matrix type encountered.\n");
        }
        
        if(error_handling == 1) // If it is still one here, then the program has found some error
        {
           throw std:: runtime_error("Fatal error encountered in relmat calculations.  Check your input for sanity.\n"
           "To check what relmat is doing, activate the debug flag of this code.");
        }
        
        // Convert to SI-units
        W *= w2Hz / 2.0;
        dipole /= 100.0; // sqrt(I0_hi2arts / w2Hz) = 1/100;
        
        // The temperature derivatives are for now only possible to do with perturbations
        if(ppd.do_temperature())
        {
          // Do not write debug information in this section
          Index e_tmp = 0;
          
          // Perturbed temperature
          Numeric t_dt = t + ppd.Temperature_Perturbation();
          
          Numeric QT_dt;
          
          // Perturbed partition functions
          partition_function( QT0,
                              QT_dt,
                              abs_lines_per_band[iband][0].Ti0(),
                              t_dt,
                              partition_functions.getParamType(abs_lines_per_band[iband][0].Species(), 
                                                                abs_lines_per_band[iband][0].Isotopologue()),
                              partition_functions.getParam(abs_lines_per_band[iband][0].Species(), 
                                                            abs_lines_per_band[iband][0].Isotopologue()),
                              false);
          
          
          if(relmat_type_per_band[iband] == hartman_tran_type)
          {
            arts_relmat_interface__hartmann_and_niro_type(
              &nlines, &fmin, &fmax,
              &M, &I, v0.get_c_array(), S.get_c_array(),
              gamma_air.get_c_array(),E_double_prime.get_c_array(),n_air.get_c_array(),
              upper, lower,
              g_prime, g_double_prime,
              &t_dt, &p, &QT_dt, &QT0, &mass,
              &number_of_perturbers, molecule_code_perturber, 
              iso_code_perturber, perturber_mass, vmr.get_c_array(), &e_tmp,
              W_dt.get_c_array(), dipole_dt.get_c_array(), rhoT_dt.get_c_array() );
          }
          else if(relmat_type_per_band[iband] == linear_type)
          {
            arts_relmat_interface__linear_type(
              &nlines, &fmin, &fmax,
              &M, &I, v0.get_c_array(), S.get_c_array(),
              gamma_air.get_c_array(),E_double_prime.get_c_array(),n_air.get_c_array(),
              upper, lower,
              g_prime, g_double_prime,
              &t_dt, &p, &QT_dt, &QT0, &mass,
              &number_of_perturbers, molecule_code_perturber, 
              iso_code_perturber, perturber_mass, vmr.get_c_array(), &e_tmp,
              W_dt.get_c_array(), dipole_dt.get_c_array(), rhoT_dt.get_c_array() );
          }
          
          if(e_tmp == 1) // If it is still one here, then the program has found some error
            throw std:: runtime_error("Fatal error encountered in relmat calculations.  Check your input for sanity.\n"
            "To check what relmat is doing, activate the debug flag of this code.  You passed normal calculations but failed in\n"
            "the calculations of the partial derivatives...");
          
          // Convert to SI-units
          W_dt *=  w2Hz / 2.0;
          dipole_dt /= 100.0;
        }
        
        // Use the provided pressure shift  NOTE: this might be a bad idea
        if(ppd.do_temperature())
        {
          const Numeric t_dt = abs_t[ip] + ppd.Temperature_Perturbation();
          for(Index ii = 0; ii < nlines; ii++)
          {
            psf[ii] = delta_air[ii] * abs_p[ip] * pow ((relmat_T0/abs_t[ip]),(Numeric)0.25+(Numeric)1.5*n_air[ii]);
            psf_dt[ii] = delta_air[ii] * abs_p[ip] * pow ((relmat_T0/t_dt),(Numeric)0.25+(Numeric)1.5*n_air[ii]);
          }
          
        }
        else
        {
          for(Index ii = 0; ii < nlines; ii++)
          {
            psf[ii] = delta_air[ii] * abs_p[ip] * pow ((relmat_T0/abs_t[ip]),(Numeric)0.25+(Numeric)1.5*n_air[ii]);
          }
        }
        
        // NOTE: For now, _dt means that the values at a different temperature are sent into the function for numeric
        //       partial derivation at a later stage.
        
        // Using Rodrigues etal method
        if(not error_handling)
          calculate_xsec_from_relmat( abs_xsec_per_species,
                                      dabs_xsec_per_species_dx,
                                      ppd,
                                      W,
                                      W_dt,
                                      f0,
                                      f_grid,
                                      dipole,
                                      dipole_dt,
                                      rhoT,
                                      rhoT_dt,
                                      psf,
                                      psf_dt,
                                      abs_t[ip],
                                      mass,
                                      iso_ratio,
                                      this_species,
                                      ip,
                                      nlines);
        
        if(write_relmat_per_band not_eq 0)
        {
          relmat_per_band[ip][iband] = W;
        }
        
      }
      else 
      {
        for( long iline=0; iline<nlines; iline++ )
        {
          
          Numeric QT = -1, QT0 = -1, part_ratio, K1, K2, tmp1;
          static const Vector v_tmp(0);
          
          GetLineScalingData(QT, QT0, part_ratio, K1, K2, tmp1, tmp1, 
                             partition_functions.getParamType(abs_lines_per_band[iband][0].Species(), 
                                                              abs_lines_per_band[iband][0].Isotopologue()),
                             partition_functions.getParam(abs_lines_per_band[iband][0].Species(), 
                                                          abs_lines_per_band[iband][0].Isotopologue()), 
                             abs_t[ip],
                             abs_lines_per_band[iband][iline].Ti0(), 
                             abs_lines_per_band[iband][iline].F(),
                             abs_lines_per_band[iband][iline].Elow(),
                             false,
                             -1.0, -1.0, -1, -1, v_tmp);
          
          W = 0;
          //assert(abs_lines_per_band[iband][iline].PressureBroadening().isAirBroadening());
          abs_lines_per_band[iband][iline].PressureBroadening().GetAirBroadening(W(iline, iline),  psf[iline],
                                           abs_lines_per_band[iband][iline].Ti0()/abs_t[ip],
                                           abs_p[ip], 0.0);
            
          const Numeric S0 = abs_lines_per_band[iband][iline].I0() * part_ratio * K1 * K2 * iso_ratio;
          const Numeric F0_pseudo = abs_lines_per_band[iband][iline].F() + psf[iline];
          const Numeric sigma = abs_lines_per_band[iband][iline].F() * doppler_const * sqrt(abs_t[ip] / mass);
          const Numeric ls_normfac = invSqrtPI / sigma;
          
          const Complex z_tmp = Complex(-F0_pseudo,W(iline, iline))/sigma;
          
          Numeric ratio_dS0dT = 0.0, df0dT, dgdT;
          if(ppd.do_temperature())
          {
            Numeric dQT = -1, dpart_ratio, dK2;
            GetLineScalingData_dT(dQT, dK2, dpart_ratio, tmp1, tmp1, tmp1, QT, tmp1,
                                  partition_functions.getParamType(abs_lines_per_band[iband][0].Species(), 
                                                                   abs_lines_per_band[iband][0].Isotopologue()),
                                  partition_functions.getParam(abs_lines_per_band[iband][0].Species(), 
                                                               abs_lines_per_band[iband][0].Isotopologue()), 
                                  abs_t[ip],
                                  abs_lines_per_band[iband][iline].Ti0(), 
                                  ppd.Temperature_Perturbation(),
                                  abs_lines_per_band[iband][iline].F(),
                                  false,
                                  -1.0, -1.0, -1, -1, v_tmp );
            ratio_dS0dT = abs_lines_per_band[iband][iline].Elow() / BOLTZMAN_CONST / abs_t[ip] / abs_t[ip] + // K1-ratio 
                          dK2/K2 +  // K2-ratio
                          dQT; // Already as a ratio
                          
                          abs_lines_per_band[iband][iline].PressureBroadening().GetAirBroadening_dT(dgdT,  df0dT,
                                           abs_t[ip], abs_lines_per_band[iband][iline].Ti0(),
                                           abs_p[ip], 0.0);
          }
          
          for(Index iv = 0; iv < f_grid.nelem(); iv++)
          {
            const Complex z = z_tmp + f_grid[iv]/sigma;
            const Complex w = Faddeeva::w(z) * ls_normfac;
            
            const Numeric xsec_this = S0 * real(w);
            
            #pragma omp atomic
            abs_xsec_per_species[this_species](iv, ip) += xsec_this;
            
            Complex dwdz = 0.0;
            if(ppd.do_frequency() or ppd.do_temperature())
              dwdz = 2.0 * (z * w - ls_normfac);
            
            for(Index iq = 0; iq < ppd.nelem(); iq++)
            {
              if(ppd(iq) == JQT_frequency  or ppd(iq) == JQT_wind_magnitude  or 
                ppd(iq) == JQT_wind_u  or ppd(iq) == JQT_wind_v  or ppd(iq) == JQT_wind_w)
              {
                #pragma omp atomic
                dabs_xsec_per_species_dx[this_species][iq](iv, ip) += S0 * real(dwdz * z_tmp);
              }
              else if(ppd(iq) == JQT_temperature)
              {
                const Complex dzdT = - z / (2.0 * abs_t[ip]) + Complex(-df0dT, dgdT)/sigma;
                dabs_xsec_per_species_dx[this_species][iq](iv, ip) += 
                real( invSqrtPI * (S0 * dwdz * dzdT - 1/2/abs_t[ip] * xsec_this) ) + ratio_dS0dT * xsec_this;
              }
              
            }
          }
          
        }
        
        if(write_relmat_per_band not_eq 0)
        {
          relmat_per_band[ip][iband](joker, joker) = W;
        }
        
      }
    }
    
    delete[] g_prime;
    delete[] g_double_prime;
    delete[] upper;
    delete[] lower;
    
  }
  
  delete[] iso_code_perturber;
  delete[] molecule_code_perturber;
    
}

#else
void abs_xsec_per_speciesAddLineMixedBands(  // WS Output:
                                             ArrayOfMatrix&                   abs_xsec_per_species,
                                             ArrayOfArrayOfMatrix&            dabs_xsec_per_species_dx,
                                             ArrayOfArrayOfMatrix&            relmat_per_band,
                                             // WS Input:                     
                                             const ArrayOfArrayOfLineRecord&  abs_lines_per_band,
                                             const ArrayOfArrayOfSpeciesTag&  abs_species_per_band,
                                             const ArrayOfArrayOfSpeciesTag&  abs_species,
                                             const SpeciesAuxData&            isotopologue_ratios,
                                             const SpeciesAuxData&            partition_functions,
                                             const ArrayOfRetrievalQuantity&  jacobian_quantities,
                                             const Vector&                    f_grid,
                                             const Vector&                    abs_p,
                                             const Vector&                    abs_t,
                                             const Numeric&                   lm_p_lim,
                                             const ArrayOfIndex&              relmat_type_per_band,
                                             const Index&                     write_relmat_per_band,
                                             const Index&                     debug,
                                             const Verbosity&)
{
   throw std::runtime_error("This version of ARTS was compiled without external line mixing support.");
}

#endif //ENABLE_RELMAT

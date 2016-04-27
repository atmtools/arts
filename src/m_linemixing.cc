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

void abs_lines_per_bandInit(ArrayOfArrayOfLineRecord& abs_lines_per_band,
                            ArrayOfArrayOfMatrix&     relmat_per_band,
                            const Verbosity&)
{
    // These are now initialized and zero-sized.  All additional functions operating on these variables should work in unison
    abs_lines_per_band.resize(0);
    relmat_per_band.resize(0);
}

void abs_lines_per_bandLineMixingAppendCO2( ArrayOfArrayOfLineRecord& abs_lines_per_band,
                                            ArrayOfArrayOfMatrix&     relmat_per_band,
                                            const String&             bandinfo_file,
                                            const Numeric&            rel_str,
                                            const Numeric&            fmin,
                                            const Numeric&            fmax,
                                            const Verbosity&          verbosity)
{
    
    extern const Numeric HZ2CM;
    
    CREATE_OUT2;
    out2<< "Appending line mixing band following:\n\t"<<bandinfo_file<<"\n";
    
    // This will contain all the read file
    std::ifstream band_stream;
    
    // Store relative pathing for all the other files
    open_input_file(band_stream, bandinfo_file);
    ArrayOfString bandpath;
    bandinfo_file.split(bandpath, "/");
    String bandfolder="";
    for(Index i=0;i<bandpath.nelem()-1;i++)
    {
        bandfolder.append(bandpath[i]);
        bandfolder.append("/");
    }
    
    // This will contain the information of one line
    String line;
    
    // Since fortran uses "D" for long floats and C uses "E" for all floats...
    const String c ="E", f = "D";
    
    // Setup test of relative strength
    Numeric first_str=0;
    bool    first_loop=true;
    
    while(true)
    {
        // Get line and check if we are done with the file yet
        std::getline(band_stream, line);
        if(line.nelem()==0)
            break;
        else if(line.nelem()!=72)
            throw std::runtime_error("The band info is bad. Check the file.\n");
        
        // The filename should be connected to the first X characters
        String filename;
        filename = line.erase(0,13);
        filename = "S" + filename + ".dat";
        
        // The maximum cross-section is in the next X numbers [HITRAN unit line strength?]
        Numeric max_xsec;
        extract(max_xsec, line, 12);
        
        // Test relative xsec strength --- see documentation in the end
        if(first_loop)
        {
            first_loop=false;
            first_str=max_xsec;
        }
        else
        {
            if(max_xsec/first_str<rel_str)
                continue;
        }
        
        // The minimum frequency is in the next X numbers [Frequency]
        Numeric min_freq;
        extract(min_freq, line, 13);
        min_freq /= HZ2CM;
        
        // The maximum frequency is in the next X numbers [Frequency]
        Numeric max_freq;
        extract(max_freq, line, 13);
        max_freq /= HZ2CM;
        
        // Test the frequency range before continuing
        if(min_freq>fmax||max_freq<fmin)
            continue;
        
        // The wfit file
        String wfilename="WTfitXY.dat"; 
        const Index  ch1_int = (Index)filename[3]-48, ch2_int = (Index)filename[8]-48;
        
        // If XY is not 1 or 0 apart then skip
        if( abs( ch1_int - ch2_int ) > 1 )
        {
            continue;
        }
        else if(ch1_int<=ch2_int&&ch1_int<=5)
        {
            wfilename[5] = filename[3];
            wfilename[6] = filename[8];
        }
        else
        {
            continue;
        }
        
        // The number of lines per band (-1 is no lines of type)
        Index PJ_max, QJ_max, RJ_max;
        extract(PJ_max, line, 12);
        extract(QJ_max, line, 4);
        extract(RJ_max, line, 4);
        
        // The file paths for hitran-like and for relmat
        String path_hitran = bandfolder, path_relmat = bandfolder;
        path_hitran.append(filename);
        path_relmat.append(wfilename);
        
        // HITRAN reading for all lines in the file
        ArrayOfLineRecord this_band;
        std::ifstream hitran_stream;
        open_input_file(hitran_stream, path_hitran);
        while(! hitran_stream.eof())
        {
            LineRecord          one_line;
            if(one_line.ReadFromHitranModifiedStream(hitran_stream, verbosity)) // NOTE: still add extra line information?
                break;
            else
                one_line.ARTSCAT5FromARTSCAT3();
            this_band.push_back(one_line);
        }
        hitran_stream.close();
        
        // Create and set relmats.
        const Index nlines = this_band.nelem();
        ArrayOfMatrix relmats(4);
        for(Index i=0;i<4;i++)
        {
                relmats[i].resize(nlines,nlines);
                relmats[i] = 0.0;
        }
        
        // Find the interesting quantum numbers in simple manner
        ArrayOfIndex JUPPER_list(nlines),JLOWER_list(nlines);
        for(Index i=0;i<nlines;i++)
        {
            LineRecord& lr = this_band[i];
            JUPPER_list[i] = lr.Upper_J().toIndex();
            JLOWER_list[i] = lr.Lower_J().toIndex();
        }
        
        // Relmat routine --- ignore data that is none-existent (strange folder structures)
        std::ifstream relmat_stream;
        try
        {
            open_input_file(relmat_stream, path_relmat);
        }
        catch (const runtime_error& error)
        {
            continue;
        }
        
        while(true)
        {
            
            // old line is used so use line anew
            String orig_line, this_line="";
            std::getline(relmat_stream, orig_line);
            if(relmat_stream.eof())
                break;
            
            // c and fortran long float conversion
            ArrayOfString tmp;
            orig_line.split(tmp,f);
            for(Index i=0;i<tmp.nelem()-1;i++)
            {
                this_line.append(tmp[i]);
                this_line.append(c);
            }
            this_line.append(tmp[tmp.nelem()-1]);
            
            
            // Format of relmat data
            Numeric W0, W0_T, some_value, rel_err;
            Index   Jupper,Jlower,Jupper_p,Jlower_p;
            std::istringstream icecream(this_line);
            icecream >> W0 >> W0_T >> some_value >> rel_err >> Jupper >> Jlower >> Jupper_p >> Jlower_p;
            
            Numeric scale;
            if(Jupper&&Jupper_p)
                scale = ((Numeric)(Jupper*(Jupper+1))/((Numeric)(Jupper_p*(Jupper_p+1))));
            else if( Jupper )
                scale = (Numeric)(Jupper*(Jupper+1));
            else if (Jupper_p)
                scale = 1.0;
            else 
                scale = 0.0; // no scale?
            
            // Position of relmat data in remats
            for(Index i=0;i<nlines;i++)
            {
                // Add numerics for diagonal here
                if(Jupper==JUPPER_list[i])
                {
                    if(Jlower==JLOWER_list[i])
                    {
                        for(Index j=0;j<nlines;j++)
                        {
                            if(j==i)
                                continue; // This case is artificial and handled elsewhere
                            
                            if(Jupper_p==JUPPER_list[j])
                            {
                                if(Jlower_p==JLOWER_list[j])
                                {
                                    relmats[0](i,j)=W0;
                                    relmats[0](j,i)=W0*scale; //if 
                                    
                                    relmats[1](i,j)=W0_T;
                                    relmats[1](j,i)=W0_T; // scale? 1/W0_T?
                                    
                                    relmats[2](i,j)=some_value;
                                    relmats[2](j,i)=some_value; // scale?
                                    
                                    relmats[3](i,j)=rel_err;
                                    relmats[3](j,i)=rel_err; // scale?
                                }
                            }
                        }
                    }
                }
            }
            
        }
        relmat_stream.close();
        
        // End of one line
        abs_lines_per_band.push_back(this_band);
        relmat_per_band.push_back(relmats);
    }
    band_stream.close();
}

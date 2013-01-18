/* Copyright (C) 2012 Oliver Lemke <olemke@core-dump.info> and Stefan 
   Buehler <sbuehler@ltu.se>. 

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

/*!
  \file   m_cia.cc
  \author Oliver Lemke <olemke@core-dump.info>
  \author Stefan Buehler
  \date   2012-12-04

  \brief  Workspace methods for HITRAN CIA data.

*/

#include "arts.h"
#include "absorption.h"
#include "file.h"
#include "cia.h"
#include "messages.h"
#include "physics_funcs.h"


extern const Numeric SPEED_OF_LIGHT;

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddCIA(// WS Output:
                                ArrayOfMatrix& abs_xsec_per_species,
                                // WS Input:
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Vector& f_grid,
                                const Vector& abs_p,
                                const Vector& abs_t,
                                const Matrix& abs_vmrs,
                                const ArrayOfCIARecord& abs_cia_data,
                                // Verbosity object:
                                const Verbosity& verbosity)
{
    CREATE_OUTS;

  {
    // Check that all parameters that should have the number of tag
    // groups as a dimension are consistent:
    const Index n_tgs    = abs_species.nelem();
    const Index n_xsec   = abs_xsec_per_species.nelem();
    const Index nr_vmrs  = abs_vmrs.nrows();
//    const Index n_cia    = abs_cia_data.nelem();

    if (n_tgs != n_xsec  ||
        n_tgs != nr_vmrs) //  ||
//        n_tgs != n_cia)
      {
        ostringstream os;
        os << "The following variables must all have the same dimension:\n"
           << "abs_species:          " << n_tgs << "\n"
           << "abs_xsec_per_species: " << n_xsec << "\n"
           << "abs_vmrs.nrows:       " << nr_vmrs << "\n";
//           << "abs_cia_data:         " << n_cia;
        throw runtime_error(os.str());
      }
  }

  {
    // Check that all parameters that should have the the dimension of p_grid
    // are consistent:
    const Index n_p      = abs_p.nelem();
    const Index n_t      = abs_t.nelem();
    const Index nc_vmrs  = abs_vmrs.ncols();
    
    if (n_p != n_t  ||
        n_p != nc_vmrs)
      {
        ostringstream os;
        os << "The following variables must all have the same dimension:\n"
        << "abs_p:          " << n_p << "\n"
        << "abs_t:          " << n_t << "\n"
        << "abs_vmrs.ncols: " << nc_vmrs;
        throw runtime_error(os.str());
      }
  }
    
    // Allocate a vector with dimension frequencies for constructing our
    // cross-sections before adding them (more efficient to allocate this here
    // outside of the loops)
    Vector xsec_temp(f_grid.nelem());
    
    // Loop over CIA data sets.
    // Index i loops through the outer array (different tag groups),
    // index s through the inner array (different tags within each goup).
    for (Index i = 0; i < abs_species.nelem(); i++)
        for (Index s = 0; s < abs_species[i].nelem(); s++)
          {
            const SpeciesTag& this_species = abs_species[i][s];

            // Check if this is a CIA tag
            if (this_species.CIASecond() == -1)
                continue;
              
            // Get convenient references of this CIA data record and this
            // absorption cross-section record:
            Index this_cia_index = cia_get_index(abs_cia_data,
                                                 this_species.Species(),
                                                 this_species.CIASecond());

            assert(this_cia_index != -1);

            const CIARecord& this_cia = abs_cia_data[this_cia_index];
            Matrix&          this_xsec = abs_xsec_per_species[i];
            
            {
              // Some nice output to out2:
              ostringstream os;
              os << "  CIA Species found: "
                 << this_cia.MoleculeName(0) << "-CIA-"
                 << this_cia.MoleculeName(1) << "\n";
              out2 << os.str();
            }
            
            // Check that the dimension of this_xsec is
            // consistent with abs_p and f_grid.
            if (this_xsec.nrows()!=f_grid.nelem()) {
                ostringstream os;
                os << "Wrong dimension of abs_xsec_per_species.nrows for species "
                   << i << ":\n"
                   << "should match f_grid (" << f_grid.nelem() << ") but is "
                   << this_xsec.nrows() << ".";
                throw runtime_error(os.str());
            }
            if (this_xsec.ncols()!=abs_p.nelem()) {
                ostringstream os;
                os << "Wrong dimension of abs_xsec_per_species.ncols for species "
                << i << ":\n"
                << "should match abs_p (" << abs_p.nelem() << ") but is "
                << this_xsec.ncols() << ".";
                throw runtime_error(os.str());
            }
            
            // Find out index of VMR for the second CIA species.
            // (The index for the first species is simply i.)
            Index i_sec = find_first_species_tg(abs_species, this_cia.Species(1));

            // Loop over pressure:
            for (Index ip = 0; ip < abs_p.nelem(); ip++)
              {
                // Get the binary absorption cross sections from the CIA data:
                this_cia.Extract(xsec_temp, f_grid, abs_t[ip], this_species.CIADataset());
                
                // We have to multiply with the number density of the second CIA species.
                // We do not have to multiply with the first, since we still
                // want to return a (unary) absorption cross-section, not an
                // absorption coefficient.
                
                // Calculate number density from pressure and temperature.
                const Numeric n = abs_vmrs(i_sec,ip) * number_density(abs_p[ip],abs_t[ip]);
                xsec_temp *= n;
                
                // Add to result variable:
                this_xsec(joker,ip) += xsec_temp;
              }
            
          }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cia_dataReadFromCIA(// WS Output:
                            ArrayOfCIARecord& abs_cia_data,
                            // WS Input:
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const String& catalogpath,
                            const Verbosity& verbosity)
{
    ArrayOfString subfolders;
    subfolders.push_back("Main-Folder/");
    subfolders.push_back("Alternate-Folder/");

    abs_cia_data.resize(0);

    // Loop species tag groups to find CIA tags.
    // Index sp loops through the tag groups, index iso through the tags within
    // each group. Despite the name, iso does not denote the isotope!
    for (Index sp = 0; sp < abs_species.nelem(); sp++)
    {
        for (Index iso = 0; iso < abs_species[sp].nelem(); iso++)
        {
            if (abs_species[sp][iso].Type() != SpeciesTag::TYPE_CIA)
                continue;
            
            ostringstream cia_name;

            Index cia_index = cia_get_index(abs_cia_data,
                                            abs_species[sp][iso].Species(),
                                            abs_species[sp][iso].CIASecond());

            // If cia_index is not -1, we have already read this datafile earlier
            if (cia_index != -1)
                continue;

            cia_name
            << species_name_from_species_index(abs_species[sp][iso].Species()) << "-"
            << species_name_from_species_index(abs_species[sp][iso].CIASecond());

            if (cia_name)
            {
                ArrayOfString checked_dirs;

                bool found = false;
                for (Index dir = 0; !found && dir < subfolders.nelem(); dir++)
                {
                    ArrayOfString files;
                    checked_dirs.push_back(catalogpath + "/"
                                           + subfolders[dir]
                                           + cia_name.str() + "/");
                    try {
                        list_directory(files, *(checked_dirs.end()-1));
                    }
                    catch (runtime_error e) {
                        continue;
                    }

                    for (Index i = files.nelem()-1; i >= 0; i--)
                    {
                        if (files[i].find(cia_name.str()) != 0
                            || files[i].rfind(".cia") != files[i].length() - 4)
                        {
                            files.erase(files.begin() + i);
                        }
                    }
                    if (files.nelem())
                    {
                        CIARecord ciar;

                        found = true;
                        String catfile = *(checked_dirs.end()-1) + files[0];

                        ciar.SetSpecies(abs_species[sp][iso].Species(),
                                        abs_species[sp][iso].CIASecond());
                        ciar.ReadFromCIA(catfile, verbosity);

                        abs_cia_data.push_back(ciar);
                    }
                }

                if (!found)
                {
                    ostringstream os;
                    os << "Error: No data file found for CIA species "
                    << cia_name.str() << endl
                    << "Looked in directories: " << checked_dirs;

                    throw runtime_error(os.str());
                }
            }
        }
    }
}


/* Workspace method: Doxygen documentation will be auto-generated */
void CIAInfo(// Generic Input:
             const String& catalogpath,
             const ArrayOfString& cia_tags,
             const Verbosity& verbosity)
{
    CREATE_OUT1;

    ArrayOfArrayOfSpeciesTag species_tags;

    for (Index i = 0; i < cia_tags.nelem(); i++)
    {
        ArrayOfSpeciesTag this_species_tag;

        ArrayOfString species_names;

        cia_tags[i].split(species_names, "-");

        if (species_names.nelem() != 2)
        {
            ostringstream os;
            os << "ERROR: Cannot parse CIA tag: " << cia_tags[i];
            throw runtime_error(os.str());
        }

        this_species_tag.push_back(SpeciesTag(species_names[0]
                                              + "-CIA-"
                                              + species_names[1]
                                              + "-0"));

        species_tags.push_back(this_species_tag);
    }

    ArrayOfCIARecord cia_data;

    abs_cia_dataReadFromCIA(cia_data, species_tags, catalogpath, verbosity);

    out1 << "CIA tag; Spectral range [cm-1]; Temp range [K]; # of sets\n";
    for (Index i = 0; i < cia_data.nelem(); i++)
        for (Index j = 0; j < cia_data[i].DatasetCount(); j++)
        {
            Vector temp_grid = cia_data[i].TemperatureGrid(j);
            Vector freq_grid = cia_data[i].FrequencyGrid(j);

            ostringstream os;
            os << setprecision(2) << std::fixed << "  "
            << cia_data[i].MoleculeName(0) << "-CIA-" << cia_data[i].MoleculeName(1)
            << "-" << j
            << "; " << freq_grid[0] / 100. / SPEED_OF_LIGHT
            << " - " << freq_grid[freq_grid.nelem()-1] / 100. / SPEED_OF_LIGHT
            << std::fixed
            << "; " << temp_grid[0] << " - " << temp_grid[temp_grid.nelem()-1]
            << "; " << temp_grid.nelem()
            << "\n";

            out1 << os.str();
        }
}

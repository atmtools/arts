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

  \brief  Workspace methods for CIA catalog.

*/

#include "arts.h"
#include "absorption.h"
#include "file.h"
#include "cia.h"

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cia_dataInit(// WS Output:
                     ArrayOfArrayOfCIARecord& abs_cia_data,
                     const Verbosity&)
{
    abs_cia_data.resize(0);
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddCIA(// WS Output:
                                ArrayOfMatrix& abs_xsec_per_species,
                                // WS Input:
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Vector& f_grid,
                                const Vector& abs_p,
                                const Vector& abs_t,
                                const Matrix& abs_vmrs,
                                const ArrayOfArrayOfCIARecord& abs_cia_data,
                                // Verbosity object:
                                const Verbosity& verbosity)
{
  // FIXME
}


/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cia_dataReadFromCIA(// WS Output:
                            ArrayOfArrayOfCIARecord& abs_cia_data,
                            // WS Input:
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const String& catalogpath,
                            const Verbosity& verbosity)
{
    ArrayOfString subfolders;
    subfolders.push_back("Main-Folder/");
    subfolders.push_back("Alternate-Folder/");
    
    abs_cia_data.resize(abs_species.nelem());

    for (Index sp = 0; sp < abs_species.nelem(); sp++)
    {
        for (Index iso = 0; iso < abs_species[sp].nelem(); iso++)
        {
            if (abs_species[sp][iso].Cia() == -1)
                continue;
            
            ostringstream cia_name;

            cia_name
            << species_name_from_species_index(abs_species[sp][iso].Species()) << "-"
            << species_name_from_species_index(abs_species[sp][iso].Cia());

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
                                        abs_species[sp][iso].Cia());
                        ciar.ReadFromCIA(catfile, verbosity);

                        abs_cia_data[sp].push_back(ciar);
                    }
                }

                if (!found)
                {
                    ostringstream os;
                    os << "Error: No catalog file found for CIA species "
                    << cia_name.str() << endl
                    << "Looked in directories: " << checked_dirs;

                    throw runtime_error(os.str());
                }
            }
        }
    }
}

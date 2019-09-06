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

#include "absorption.h"
#include "arts.h"
#include "auto_md.h"
#include "cia.h"
#include "file.h"
#include "messages.h"
#include "physics_funcs.h"
#include "xml_io.h"

extern const Numeric SPEED_OF_LIGHT;

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddCIA(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const Matrix& abs_vmrs,
    const ArrayOfCIARecord& abs_cia_data,
    // WS Generic Input:
    const Numeric& T_extrapolfac,
    const Index& robust,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUTS;

  {
    // Check that all parameters that should have the number of tag
    // groups as a dimension are consistent:
    const Index n_tgs = abs_species.nelem();
    const Index n_xsec = abs_xsec_per_species.nelem();
    const Index nr_vmrs = abs_vmrs.nrows();
    //    const Index n_cia    = abs_cia_data.nelem();

    if (n_tgs != n_xsec || n_tgs != nr_vmrs)  //  ||
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

  // Jacobian overhead START
  /* NOTE:  The calculations below are inefficient and could 
            be made much better by using interp in Extract to
            return the derivatives as well. */
  const bool do_jac = supports_CIA(jacobian_quantities);
  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  Vector dfreq, dabs_t;
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);

  if (do_freq_jac) {
    dfreq = f_grid;
    dfreq += df;
  }
  if (do_temp_jac) {
    dabs_t.resize(abs_t.nelem());
    dabs_t = abs_t;
    dabs_t += dt;
  }
  // Jacobian overhead END

  // Useful if there is no Jacobian to calculate
  ArrayOfMatrix empty;

  {
    // Check that all parameters that should have the the dimension of p_grid
    // are consistent:
    const Index n_p = abs_p.nelem();
    const Index n_t = abs_t.nelem();
    const Index nc_vmrs = abs_vmrs.ncols();

    if (n_p != n_t || n_p != nc_vmrs) {
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

  // Jacobian vectors START
  Vector dxsec_temp_dT;
  Vector dxsec_temp_dF;
  if (do_freq_jac) dxsec_temp_dF.resize(f_grid.nelem());
  if (do_temp_jac) dxsec_temp_dT.resize(f_grid.nelem());
  // Jacobian vectors END

  // Loop over CIA data sets.
  // Index ii loops through the outer array (different tag groups),
  // index s through the inner array (different tags within each goup).
  for (Index ii = 0; ii < abs_species_active.nelem(); ii++) {
    const Index i = abs_species_active[ii];

    for (Index s = 0; s < abs_species[i].nelem(); s++) {
      const SpeciesTag& this_species = abs_species[i][s];

      // Check if this is a CIA tag
      if (this_species.Type() != SpeciesTag::TYPE_CIA) continue;

      // Get convenient references of this CIA data record and this
      // absorption cross-section record:
      Index this_cia_index = cia_get_index(
          abs_cia_data, this_species.Species(), this_species.CIASecond());

      assert(this_cia_index != -1);

      const CIARecord& this_cia = abs_cia_data[this_cia_index];
      Matrix& this_xsec = abs_xsec_per_species[i];

      if (out2.sufficient_priority()) {
        // Some nice output to out2:
        out2 << "  CIA Species found: " + this_species.Name() + "\n";
      }

      // Check that the dimension of this_xsec is
      // consistent with abs_p and f_grid.
      if (this_xsec.nrows() != f_grid.nelem()) {
        ostringstream os;
        os << "Wrong dimension of abs_xsec_per_species.nrows for species " << i
           << ":\n"
           << "should match f_grid (" << f_grid.nelem() << ") but is "
           << this_xsec.nrows() << ".";
        throw runtime_error(os.str());
      }
      if (this_xsec.ncols() != abs_p.nelem()) {
        ostringstream os;
        os << "Wrong dimension of abs_xsec_per_species.ncols for species " << i
           << ":\n"
           << "should match abs_p (" << abs_p.nelem() << ") but is "
           << this_xsec.ncols() << ".";
        throw runtime_error(os.str());
      }

      // Find out index of VMR for the second CIA species.
      // (The index for the first species is simply i.)
      Index i_sec = find_first_species_tg(abs_species, this_cia.Species(1));

      // Catch the case that the VMR for the second species does not exist:
      if (i_sec < 0) {
        ostringstream os;
        os << "VMR profile for second species in CIA species pair does not exist.\n"
           << "Tag " << this_species.Name() << " needs a VMR profile of "
           << this_cia.MoleculeName(1) << "!";
        throw runtime_error(os.str());
      }

      // Loop over pressure:
      for (Index ip = 0; ip < abs_p.nelem(); ip++) {
        // Get the binary absorption cross sections from the CIA data:
        try {
          this_cia.Extract(xsec_temp,
                           f_grid,
                           abs_t[ip],
                           this_species.CIADataset(),
                           T_extrapolfac,
                           robust,
                           verbosity);
          if (do_freq_jac)
            this_cia.Extract(dxsec_temp_dF,
                             dfreq,
                             abs_t[ip],
                             this_species.CIADataset(),
                             T_extrapolfac,
                             robust,
                             verbosity);
          if (do_temp_jac)
            this_cia.Extract(dxsec_temp_dT,
                             f_grid,
                             dabs_t[ip],
                             this_species.CIADataset(),
                             T_extrapolfac,
                             robust,
                             verbosity);
        } catch (const std::runtime_error& e) {
          ostringstream os;
          os << "Problem with CIA species " << this_species.Name() << ":\n"
             << e.what();
          throw runtime_error(os.str());
        }

        // We have to multiply with the number density of the second CIA species.
        // We do not have to multiply with the first, since we still
        // want to return a (unary) absorption cross-section, not an
        // absorption coefficient.

        // Calculate number density from pressure and temperature.
        const Numeric n =
            abs_vmrs(i_sec, ip) * number_density(abs_p[ip], abs_t[ip]);

        if (!do_jac) {
          xsec_temp *= n;
          // Add to result variable:
          this_xsec(joker, ip) += xsec_temp;
        } else {
          const Numeric dn_dT =
              abs_vmrs(i_sec, ip) * dnumber_density_dt(abs_p[ip], abs_t[ip]);

          for (Index iv = 0; iv < xsec_temp.nelem(); iv++) {
            this_xsec(iv, ip) += n * xsec_temp[iv];
            for (Index iq = 0; iq < jacobian_quantities_position.nelem();
                 iq++) {
              if (is_frequency_parameter(
                      jacobian_quantities[jacobian_quantities_position[iq]]))
                dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                    n * (dxsec_temp_dF[iv] - xsec_temp[iv]) / df;
              else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                       JacPropMatType::Temperature)
                dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                    n * (dxsec_temp_dT[iv] - xsec_temp[iv]) / dt +
                    xsec_temp[iv] * dn_dT;
              else if (species_match(jacobian_quantities
                                         [jacobian_quantities_position[iq]],
                                     this_species.BathSpecies()))
                dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                    number_density(abs_p[ip], abs_t[ip]) * xsec_temp[iv];
            }
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void CIARecordReadFromFile(  // WS GOutput:
    CIARecord& cia_record,
    // WS Generic Input:
    const String& species_tag,
    const String& filename,
    // Verbosity object:
    const Verbosity& verbosity) {
  SpeciesTag species(species_tag);

  if (species.Type() != SpeciesTag::TYPE_CIA) {
    ostringstream os;
    os << "Invalid species tag " << species_tag << ".\n"
       << "This is not recognized as a CIA type.\n";
    throw std::runtime_error(os.str());
  }

  cia_record.SetSpecies(species.Species(), species.CIASecond());
  cia_record.ReadFromCIA(filename, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cia_dataAddCIARecord(  // WS Output:
    ArrayOfCIARecord& abs_cia_data,
    // WS GInput:
    const CIARecord& cia_record,
    const Index& clobber,
    // WS Input:
    const Verbosity&) {
  Index cia_index =
      cia_get_index(abs_cia_data, cia_record.Species(0), cia_record.Species(1));
  if (cia_index == -1)
    abs_cia_data.push_back(cia_record);
  else if (clobber)
    abs_cia_data[cia_index] = cia_record;
  else
    abs_cia_data[cia_index].AppendDataset(cia_record);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cia_dataReadFromCIA(  // WS Output:
    ArrayOfCIARecord& abs_cia_data,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& catalogpath,
    const Verbosity& verbosity) {
  ArrayOfString subfolders;
  subfolders.push_back("Main-Folder/");
  subfolders.push_back("Alternate-Folder/");

  abs_cia_data.resize(0);

  // Loop species tag groups to find CIA tags.
  // Index sp loops through the tag groups, index iso through the tags within
  // each group. Despite the name, iso does not denote the isotope!
  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    for (Index iso = 0; iso < abs_species[sp].nelem(); iso++) {
      if (abs_species[sp][iso].Type() != SpeciesTag::TYPE_CIA) continue;

      ArrayOfString cia_names;

      Index cia_index = cia_get_index(abs_cia_data,
                                      abs_species[sp][iso].Species(),
                                      abs_species[sp][iso].CIASecond());

      // If cia_index is not -1, we have already read this datafile earlier
      if (cia_index != -1) continue;

      cia_names.push_back(
          species_name_from_species_index(abs_species[sp][iso].Species()) +
          "-" +
          species_name_from_species_index(abs_species[sp][iso].CIASecond()));

      cia_names.push_back(
          species_name_from_species_index(abs_species[sp][iso].CIASecond()) +
          "-" +
          species_name_from_species_index(abs_species[sp][iso].Species()));

      ArrayOfString checked_dirs;

      bool found = false;
      for (Index fname = 0; !found && fname < cia_names.nelem(); fname++) {
        String cia_name = cia_names[fname];

        for (Index dir = 0; !found && dir < subfolders.nelem(); dir++) {
          ArrayOfString files;
          checked_dirs.push_back(catalogpath + "/" + subfolders[dir] +
                                 cia_name + "/");
          try {
            list_directory(files, *(checked_dirs.end() - 1));
          } catch (const std::runtime_error& e) {
            continue;
          }

          for (Index i = files.nelem() - 1; i >= 0; i--) {
            if (files[i].find(cia_name) != 0 ||
                files[i].rfind(".cia") != files[i].length() - 4) {
              files.erase(files.begin() + i);
            }
          }
          if (files.nelem()) {
            CIARecord ciar;

            found = true;
            String catfile = *(checked_dirs.end() - 1) + files[0];

            ciar.SetSpecies(abs_species[sp][iso].Species(),
                            abs_species[sp][iso].CIASecond());
            ciar.ReadFromCIA(catfile, verbosity);

            abs_cia_data.push_back(ciar);
          }
        }
      }

      if (!found) {
        ostringstream os;
        os << "Error: No data file found for CIA species " << cia_names[0]
           << endl
           << "Looked in directories: " << checked_dirs;

        throw runtime_error(os.str());
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cia_dataReadFromXML(  // WS Output:
    ArrayOfCIARecord& abs_cia_data,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& filename,
    const Verbosity& verbosity) {
  xml_read_from_file(filename, abs_cia_data, verbosity);

  // Check that all CIA tags from abs_species are present in the
  // XML file

  vector<String> missing_tags;

  // Loop species tag groups to find CIA tags.
  // Index sp loops through the tag groups, index iso through the tags within
  // each group. Despite the name, iso does not denote the isotope!
  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    for (Index iso = 0; iso < abs_species[sp].nelem(); iso++) {
      if (abs_species[sp][iso].Type() != SpeciesTag::TYPE_CIA) continue;

      Index cia_index = cia_get_index(abs_cia_data,
                                      abs_species[sp][iso].Species(),
                                      abs_species[sp][iso].CIASecond());

      // If cia_index is -1, this CIA tag was not present in the input file
      if (cia_index == -1) {
        missing_tags.push_back(
            species_name_from_species_index(abs_species[sp][iso].Species()) +
            "-" +
            species_name_from_species_index(abs_species[sp][iso].CIASecond()));
      }
    }
  }

  if (missing_tags.size()) {
    ostringstream os;
    bool first = true;

    os << "Error: The following CIA tag(s) are missing in input file: ";
    for (size_t i = 0; i < missing_tags.size(); i++) {
      if (!first)
        os << ", ";
      else
        first = false;
      os << missing_tags[i];
    }
    throw runtime_error(os.str());
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void CIAInfo(  // Generic Input:
    const String& catalogpath,
    const ArrayOfString& cia_tags,
    const Verbosity& verbosity) {
  CREATE_OUT1;

  ArrayOfArrayOfSpeciesTag species_tags;

  for (Index i = 0; i < cia_tags.nelem(); i++) {
    ArrayOfSpeciesTag this_species_tag;

    ArrayOfString species_names;

    cia_tags[i].split(species_names, "-");

    if (species_names.nelem() != 2) {
      ostringstream os;
      os << "ERROR: Cannot parse CIA tag: " << cia_tags[i];
      throw runtime_error(os.str());
    }

    this_species_tag.push_back(
        SpeciesTag(species_names[0] + "-CIA-" + species_names[1] + "-0"));

    species_tags.push_back(this_species_tag);
  }

  ArrayOfCIARecord cia_data;

  abs_cia_dataReadFromCIA(cia_data, species_tags, catalogpath, verbosity);

  Print(cia_data, 1, verbosity);
}

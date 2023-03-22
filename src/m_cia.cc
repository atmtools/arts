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

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include "absorption.h"
#include "arts.h"
#include "arts_constants.h"
#include "auto_md.h"
#include "cia.h"
#include "debug.h"
#include "file.h"
#include "messages.h"
#include "physics_funcs.h"
#include "species.h"
#include "species_tags.h"
#include "xml_io.h"

inline constexpr Numeric SPEED_OF_LIGHT=Constant::speed_of_light;

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

    ARTS_USER_ERROR_IF (n_tgs != n_xsec || n_tgs != nr_vmrs,
      "The following variables must all have the same dimension:\n"
      "abs_species:          ", n_tgs, "\n"
      "abs_xsec_per_species: ", n_xsec, "\n"
      "abs_vmrs.nrows:       ", nr_vmrs, "\n")
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

    ARTS_USER_ERROR_IF (n_p != n_t || n_p != nc_vmrs,
        "The following variables must all have the same dimension:\n"
        "abs_p:          ", n_p, "\n"
        "abs_t:          ", n_t, "\n"
        "abs_vmrs.ncols: ", nc_vmrs)
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
      if (this_species.Type() != Species::TagType::Cia) continue;

      // Get convenient references of this CIA data record and this
      // absorption cross-section record:
      Index this_cia_index = cia_get_index(
          abs_cia_data, this_species.Spec(), this_species.cia_2nd_species);

      ARTS_ASSERT(this_cia_index != -1);

      const CIARecord& this_cia = abs_cia_data[this_cia_index];
      Matrix& this_xsec = abs_xsec_per_species[i];

      if (out2.sufficient_priority()) {
        // Some nice output to out2:
        out2 << "  CIA Species found: " + this_species.Name() + "\n";
      }

      // Check that the dimension of this_xsec is
      // consistent with abs_p and f_grid.
      ARTS_USER_ERROR_IF (this_xsec.nrows() != f_grid.nelem(),
          "Wrong dimension of abs_xsec_per_species.nrows for species ", i,
          ":\n"
          "should match f_grid (", f_grid.nelem(), ") but is ",
          this_xsec.nrows(), ".")
      ARTS_USER_ERROR_IF (this_xsec.ncols() != abs_p.nelem(),
          "Wrong dimension of abs_xsec_per_species.ncols for species ", i,
          ":\n"
          "should match abs_p (", abs_p.nelem(), ") but is ",
          this_xsec.ncols(), ".")

      // Find out index of VMR for the second CIA species.
      // (The index for the first species is simply i.)
      Index i_sec = find_first_species(abs_species, this_cia.Species(1));

      // Catch the case that the VMR for the second species does not exist:
      ARTS_USER_ERROR_IF (i_sec < 0,
          "VMR profile for second species in CIA species pair does not exist.\n"
          "Tag ", this_species.Name(), " needs a VMR profile of ",
          this_cia.MoleculeName(1), "!")

      // Loop over pressure:
      for (Index ip = 0; ip < abs_p.nelem(); ip++) {
        // Get the binary absorption cross sections from the CIA data:
        try {
          this_cia.Extract(xsec_temp,
                           f_grid,
                           abs_t[ip],
                           T_extrapolfac,
                           robust,
                           verbosity);
          if (do_freq_jac)
            this_cia.Extract(dxsec_temp_dF,
                             dfreq,
                             abs_t[ip],
                             T_extrapolfac,
                             robust,
                             verbosity);
          if (do_temp_jac)
            this_cia.Extract(dxsec_temp_dT,
                             f_grid,
                             dabs_t[ip],
                             T_extrapolfac,
                             robust,
                             verbosity);
        } catch (const std::runtime_error& e) {
          ARTS_USER_ERROR ("Problem with CIA species ",
                           this_species.Name(), ":\n", e.what())
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
            for (Index iq = 0; iq < jacobian_quantities.nelem();
                 iq++) {
              const auto& deriv = jacobian_quantities[iq];
            
            if (not deriv.propmattype()) continue;
              
            if (is_frequency_parameter(deriv))
                dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                    n * (dxsec_temp_dF[iv] - xsec_temp[iv]) / df;
            else if (deriv == Jacobian::Atm::Temperature)
              dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                  n * (dxsec_temp_dT[iv] - xsec_temp[iv]) / dt +
                  xsec_temp[iv] * dn_dT;
                  else if (species_match(deriv, this_species.cia_2nd_species))
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
void propmat_clearskyAddCIA(  // WS Output:
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfSpeciesTag& select_abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const AtmPoint& atm_point,
    const ArrayOfCIARecord& abs_cia_data,
    // WS Generic Input:
    const Numeric& T_extrapolfac,
    const Index& ignore_errors,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUTS;

  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();
  const Index ns = abs_species.nelem();

  // Possible things that can go wrong in this code (excluding line parameters)
  check_abs_species(abs_species);
  ARTS_USER_ERROR_IF(propmat_clearsky.NumberOfFrequencies() not_eq nf,
                     "*f_grid* must match *propmat_clearsky*")
  ARTS_USER_ERROR_IF(
      not nq and (nq not_eq dpropmat_clearsky_dx.nelem()),
      "*dpropmat_clearsky_dx* must match derived form of *jacobian_quantities*")
  ARTS_USER_ERROR_IF(
      not nq and bad_propmat(dpropmat_clearsky_dx, f_grid),
      "*dpropmat_clearsky_dx* must have frequency dim same as *f_grid*")
  ARTS_USER_ERROR_IF(any_negative(f_grid),
                     "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF(not is_increasing(f_grid),
                     "Must be sorted and increasing.")
  ARTS_USER_ERROR_IF(atm_point.temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(atm_point.pressure <= 0, "Non-positive pressure")

  // Jacobian overhead START
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);

  const bool do_jac = supports_CIA(
      jacobian_quantities);  // Throws runtime error if line parameters are wanted since we cannot know if the line is in the Continuum...
  const bool do_wind_jac = do_wind_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);

  Vector dfreq;
  Vector dabs_t{atm_point.temperature + dt};

  if (do_wind_jac) {
    dfreq.resize(f_grid.nelem());
    for (Index iv = 0; iv < f_grid.nelem(); iv++) dfreq[iv] = f_grid[iv] + df;
  }

  Vector dxsec_temp_dF, dxsec_temp_dT;

  if (do_jac) {
    if (do_wind_jac) {
      ARTS_USER_ERROR_IF(
          !std::isnormal(df), "df must be >0 and not NaN or Inf: ", df)
      dxsec_temp_dF.resize(f_grid.nelem());
    }
    if (do_temp_jac) {
      ARTS_USER_ERROR_IF(
          !std::isnormal(dt), "dt must be >0 and not NaN or Inf: ", dt)
      dxsec_temp_dT.resize(f_grid.nelem());
    }
  }
  // Jacobian overhead END

  // Useful if there is no Jacobian to calculate
  ArrayOfMatrix empty;

  // Allocate a vector with dimension frequencies for constructing our
  // cross-sections before adding them (more efficient to allocate this here
  // outside of the loops)
  Vector xsec_temp(f_grid.nelem());

  // Loop over CIA data sets.
  // Index ii loops through the outer array (different tag groups),
  // index s through the inner array (different tags within each goup).
  for (Index ispecies = 0; ispecies < ns; ispecies++) {
    if (select_abs_species.nelem() and
        select_abs_species not_eq abs_species[ispecies])
      continue;

    // Skip it if there are no species or there is Zeeman requested
    if (not abs_species[ispecies].nelem() or abs_species[ispecies].Zeeman())
      continue;

    // Go through the tags in the current tag group to see if they
    // are continuum tags:
    for (Index s = 0; s < abs_species[ispecies].nelem(); ++s) {
      const SpeciesTag& this_species = abs_species[ispecies][s];

      // Check if this is a CIA tag
      if (this_species.Type() != Species::TagType::Cia) continue;

      // Get convenient references of this CIA data record and this
      // absorption cross-section record:
      Index this_cia_index = cia_get_index(
          abs_cia_data, this_species.Spec(), this_species.cia_2nd_species);

      ARTS_ASSERT(this_cia_index != -1);

      const CIARecord& this_cia = abs_cia_data[this_cia_index];

      if (out2.sufficient_priority()) {
        // Some nice output to out2:
        out2 << "  CIA Species found: " + this_species.Name() + "\n";
      }

      // Find out index of VMR for the second CIA species.
      // (The index for the first species is simply i.)
      Index i_sec = find_first_species(abs_species, this_cia.Species(1));

      // Catch the case that the VMR for the second species does not exist:
      ARTS_USER_ERROR_IF(
          i_sec < 0,
          "VMR profile for second species in CIA species pair does not exist.\n"
          "Tag ",
          this_species.Name(),
          " needs a VMR profile of ",
          this_cia.MoleculeName(1),
          "!")

      // We have to multiply with the number density of the second CIA species.
      // We do not have to multiply with the first, since we still
      // want to return a (unary) absorption cross-section, not an
      // absorption coefficient.
      const Numeric nd_sec =
          number_density(atm_point.pressure, atm_point.temperature) * atm_point[this_cia.Species(1)];
      // Get the binary absorption cross sections from the CIA data:

      try {
        this_cia.Extract(xsec_temp,
                         f_grid,
                         atm_point.temperature,
                         T_extrapolfac,
                         ignore_errors,
                         verbosity);
        if (do_wind_jac) {
          this_cia.Extract(dxsec_temp_dF,
                           dfreq,
                           atm_point.temperature,
                           T_extrapolfac,
                           ignore_errors,
                           verbosity);
        }
        if (do_temp_jac) {
          this_cia.Extract(dxsec_temp_dT,
                           f_grid,
                           atm_point.temperature + dt,
                           T_extrapolfac,
                           ignore_errors,
                           verbosity);
        }
      } catch (const std::runtime_error& e) {
        ARTS_USER_ERROR(
            "Problem with CIA species ", this_species.Name(), ":\n", e.what())
      }

      if (!do_jac) {
        xsec_temp *= nd_sec * number_density(atm_point.pressure, atm_point.temperature) *
                     atm_point[this_cia.Species(0)];
        propmat_clearsky.Kjj() += xsec_temp;
      } else {  // The Jacobian block
        const Numeric nd = number_density(atm_point.pressure, atm_point.temperature);
        const Numeric dnd_dt =
            dnumber_density_dt(atm_point.pressure, atm_point.temperature);
        const Numeric dnd_dt_sec =
            dnumber_density_dt(atm_point.pressure, atm_point.temperature) * atm_point[this_cia.Species(1)];
        for (Index iv = 0; iv < f_grid.nelem(); iv++) {
          propmat_clearsky.Kjj()[iv] +=
              nd_sec * xsec_temp[iv] * nd * atm_point[this_cia.Species(0)];
          for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
            const auto& deriv = jacobian_quantities[iq];

            if (not deriv.propmattype()) continue;

            if (is_wind_parameter(deriv)) {
              dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                  nd_sec * (dxsec_temp_dF[iv] - xsec_temp[iv]) / df * nd *
                  atm_point[this_cia.Species(1)];
            } else if (deriv == Jacobian::Atm::Temperature) {
              dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                  ((nd_sec * (dxsec_temp_dT[iv] - xsec_temp[iv]) / dt +
                    xsec_temp[iv] * dnd_dt_sec) *
                       nd +
                   xsec_temp[iv] * nd_sec * dnd_dt) *
                  atm_point[this_cia.Species(0)];
            } else if (deriv == abs_species[ispecies]) {
              dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                  nd_sec * xsec_temp[iv] * nd;
            } else if (species_match(deriv, this_species.Spec())) {
              dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                  nd * nd_sec * xsec_temp[iv] *
                  (this_species.cia_2nd_species == this_species.Spec() ? 2.0
                                                                       : 1.0);
            } else if (species_match(deriv, this_species.cia_2nd_species)) {
              dpropmat_clearsky_dx[iq].Kjj()[iv] +=
                  nd * nd * xsec_temp[iv] * atm_point[this_cia.Species(0)] *
                  (this_species.cia_2nd_species == this_species.Spec() ? 2.0
                                                                       : 1.0);
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

  ARTS_USER_ERROR_IF (species.Type() != Species::TagType::Cia,
      "Invalid species tag ", species_tag, ".\n"
      "This is not recognized as a CIA type.\n")

  cia_record.SetSpecies(species.Spec(), species.cia_2nd_species);
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
      if (abs_species[sp][iso].Type() != Species::TagType::Cia) continue;

      ArrayOfString cia_names;

      Index cia_index = cia_get_index(abs_cia_data,
                                      abs_species[sp][iso].Spec(),
                                      abs_species[sp][iso].cia_2nd_species);

      // If cia_index is not -1, we have already read this datafile earlier
      if (cia_index != -1) continue;

      cia_names.push_back(String(Species::toShortName(abs_species[sp][iso].Spec())) +
          "-" +
          String(Species::toShortName(abs_species[sp][iso].cia_2nd_species)));

      cia_names.push_back(
          String(Species::toShortName(abs_species[sp][iso].cia_2nd_species)) +
          "-" +
          String(Species::toShortName(abs_species[sp][iso].Spec())));

      ArrayOfString checked_dirs;

      bool found = false;
      for (Index fname = 0; !found && fname < cia_names.nelem(); fname++) {
        String cia_name = cia_names[fname];

        for (Index dir = 0; !found && dir < subfolders.nelem(); dir++) {
          ArrayOfString files;
          checked_dirs.push_back(catalogpath + "/" + subfolders[dir] +
                                 cia_name + "/");
          try {
            files = list_directory(*(checked_dirs.end() - 1));
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

            ciar.SetSpecies(abs_species[sp][iso].Spec(),
                            abs_species[sp][iso].cia_2nd_species);
            ciar.ReadFromCIA(catfile, verbosity);

            abs_cia_data.push_back(ciar);
          }
        }
      }

      ARTS_USER_ERROR_IF (!found,
          "Error: No data file found for CIA species ", cia_names[0],
           "\n"
           "Looked in directories: ", checked_dirs)
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
      if (abs_species[sp][iso].Type() != Species::TagType::Cia) continue;

      Index cia_index = cia_get_index(abs_cia_data,
                                      abs_species[sp][iso].Spec(),
                                      abs_species[sp][iso].cia_2nd_species);

      // If cia_index is -1, this CIA tag was not present in the input file
      if (cia_index == -1) {
        missing_tags.push_back(
            String(Species::toShortName(abs_species[sp][iso].Spec())) +
            "-" +
            String(Species::toShortName(abs_species[sp][iso].cia_2nd_species)));
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
    ARTS_USER_ERROR (os.str());
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

    ARTS_USER_ERROR_IF (species_names.nelem() != 2,
      "ERROR: Cannot parse CIA tag: ", cia_tags[i])

    this_species_tag.push_back(
        SpeciesTag(species_names[0] + "-CIA-" + species_names[1] + "-0"));

    species_tags.push_back(this_species_tag);
  }

  ArrayOfCIARecord cia_data;

  abs_cia_dataReadFromCIA(cia_data, species_tags, catalogpath, verbosity);

  Print(cia_data, 1, verbosity);
}

void abs_cia_dataReadSpeciesSplitCatalog(
    ArrayOfCIARecord& abs_cia_data,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const String& basename,
    const Index& robust,
    const Verbosity& verbosity) {
  ArrayOfString names{};
  for (auto& spec : abs_species) {
    for (auto& tag : spec) {
      if (tag.type == Species::TagType::Cia) {
        names.emplace_back(
            var_string(Species::toShortName(tag.Spec()),
                       "-CIA-",
                       Species::toShortName(tag.cia_2nd_species)));
      }
    }
  }

  names.erase(std::unique(names.begin(), names.end()), names.end());

  const std::filesystem::path inpath{basename.c_str()};
  const bool is_dir = basename.back() == '/' or std::filesystem::is_directory(inpath);
  
  for (auto& name : names) {
    auto fil{inpath};
    if (is_dir) {
      fil /= name + ".xml";
    } else {
      fil += "." + name + ".xml";
    }

    xml_read_from_file(fil.c_str(), abs_cia_data.emplace_back(), verbosity);

    ARTS_USER_ERROR_IF(robust == 0 and abs_cia_data.back().DatasetCount() == 0,
                       "Cannot find any data for ", std::quoted(name),
                       " in file at ", fil)
  }
}

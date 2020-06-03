
/* Copyright (C) 2000-2012
   Stefan Buehler   <sbuehler@ltu.se>
   Patrick Eriksson <patrick.eriksson@chalmers.se>
   Axel von Engeln  <engeln@uni-bremen.de>
   Thomas Kuhn      <tkuhn@uni-bremen.de>

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

//

/**
   \file   m_abs.cc

   Stuff related to the calculation of absorption coefficients.

   \author Stefan Buehler
   \date   2001-03-12
*/
#include <algorithm>
#include <cmath>
#include "absorption.h"
#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "legacy_continua.h"
#include "file.h"
#include "global_data.h"
#include "jacobian.h"
#include "m_xml.h"
#include "math_funcs.h"
#include "matpackI.h"
#include "messages.h"
#include "montecarlo.h"
#include "optproperties.h"
#include "parameters.h"
#include "physics_funcs.h"
#include "rte.h"
#include "xml_io.h"

#ifdef ENABLE_NETCDF
#include <netcdf.h>
#include "nc_io.h"
#endif

extern const Numeric ELECTRON_CHARGE;
extern const Numeric ELECTRON_MASS;
extern const Numeric PI;
extern const Numeric SPEED_OF_LIGHT;
extern const Numeric VACUUM_PERMITTIVITY;

/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromRteScalars(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Matrix& abs_vmrs,
    // WS Input:
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const Vector& rtp_vmr,
    const Verbosity&) {
  // Prepare abs_p:
  abs_p.resize(1);
  abs_p = rtp_pressure;

  // Prepare abs_t:
  abs_t.resize(1);
  abs_t = rtp_temperature;

  // Prepare abs_vmrs:
  abs_vmrs.resize(rtp_vmr.nelem(), 1);
  abs_vmrs = rtp_vmr;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_lines_per_speciesCreateFromLines(  // WS Output:
    ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    // WS Input:
    const ArrayOfAbsorptionLines& abs_lines,
    const ArrayOfArrayOfSpeciesTag& tgs,
    const Verbosity&) {
  // Size is set but inner size will now change from the original definition of species tags...
  abs_lines_per_species.resize(tgs.nelem());

  // The inner arrays need to be emptied, because they may contain lines
  // from a previous calculation
  for (auto &tg : abs_lines_per_species)
    tg.resize(0);

  // Take copies because we have to support frequency ranges, so might have to delete
  for (AbsorptionLines lines: abs_lines) {
    
    // Skip empty lines
    if (lines.NumLines() == 0) continue;
    
    // Loop all the tags
    for (Index i=0; i<tgs.nelem() and lines.NumLines(); i++) {
      for (auto& this_tag: tgs[i]) {
        // Test species first, we might leave as we leave
        if (this_tag.Species() not_eq lines.Species()) break;
        
        // Test isotopologue, we have to hit the end of the list for no isotopologue or the exact value
        if (this_tag.Isotopologue() not_eq SpeciesDataOfBand(lines).Isotopologue().nelem() and
            this_tag.Isotopologue() not_eq lines.Isotopologue())
          continue;
        
        // If there is a frequency range, we have to check so that only selected lines are included
        if (this_tag.Lf() >= 0 or this_tag.Uf() >= 0) {
          const Numeric low = (this_tag.Lf() >= 0) ? this_tag.Lf() : std::numeric_limits<Numeric>::lowest();
          const Numeric upp = (this_tag.Uf() >= 0) ? this_tag.Uf() : std::numeric_limits<Numeric>::max();
          
          // Fill up a copy of the line record to match with the wished frequency criteria
          AbsorptionLines these_lines = createEmptyCopy(lines);
          for (Index k=lines.NumLines()-1; k>=0; k--)
            if (low <= lines.F0(k) and upp >= lines.F0(k))
              these_lines.AppendSingleLine(lines.PopLine(k));
          
          // Append these lines after sorting them if there are any of them
          if (these_lines.NumLines()) {
            these_lines.ReverseLines();
            abs_lines_per_species[i].push_back(these_lines);
          }
          
          // If this means we have deleted all lines, then we leave
          if (lines.NumLines() == 0)
            goto leave_inner_loop;
        }
        else {
          abs_lines_per_species[i].push_back(lines);
          goto leave_inner_loop;
        }
      }
    }
    leave_inner_loop: {}
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAllInScenario(  // WS Output:
    ArrayOfArrayOfSpeciesTag& tgs,
    Index& propmat_clearsky_agenda_checked,
    Index& abs_xsec_agenda_checked,
    // Control Parameters:
    const String& basename,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  // Invalidate agenda check flags
  propmat_clearsky_agenda_checked = false;
  abs_xsec_agenda_checked = false;

  // Species lookup data:
  using global_data::species_data;

  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for (Index i = 0; i < species_data.nelem(); ++i) {
    const String specname = species_data[i].Name();

    String filename = basename;
    if (basename.length() && basename[basename.length() - 1] != '/')
      filename += ".";
    filename += specname;

    try {
      find_xml_file(filename, verbosity);
      // Add to included list:
      included.push_back(specname);

      // Convert name of species to a SpeciesTag object:
      SpeciesTag this_tag(specname);

      // Create Array of SpeciesTags with length 1
      // (our tag group has only one tag):
      ArrayOfSpeciesTag this_group(1);
      this_group[0] = this_tag;

      // Add this tag group to tgs:
      tgs.push_back(this_group);
    } catch (const std::runtime_error& e) {
      // The file for the species could not be found.
      excluded.push_back(specname);
    }
  }

  // Some nice output:
  out2 << "  Included Species (" << included.nelem() << "):\n";
  for (Index i = 0; i < included.nelem(); ++i)
    out2 << "     " << included[i] << "\n";

  out2 << "  Excluded Species (" << excluded.nelem() << "):\n";
  for (Index i = 0; i < excluded.nelem(); ++i)
    out2 << "     " << excluded[i] << "\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_speciesDefineAll(  // WS Output:
    ArrayOfArrayOfSpeciesTag& abs_species,
    Index& propmat_clearsky_agenda_checked,
    Index& abs_xsec_agenda_checked,
    // Control Parameters:
    const Verbosity& verbosity) {
  // Species lookup data:
  using global_data::species_data;

  // We want to make lists of all species
  ArrayOfString specs(0);
  for (auto& spec: species_data) {
    specs.push_back(spec.Name());
  }

  // Set the values
  abs_speciesSet(abs_species, abs_xsec_agenda_checked, propmat_clearsky_agenda_checked, specs, verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void AbsInputFromAtmFields(  // WS Output:
    Vector& abs_p,
    Vector& abs_t,
    Matrix& abs_vmrs,
    // WS Input:
    const Index& atmosphere_dim,
    const Vector& p_grid,
    const Tensor3& t_field,
    const Tensor4& vmr_field,
    const Verbosity&) {
  // First, make sure that we really have a 1D atmosphere:
  if (1 != atmosphere_dim) {
    ostringstream os;
    os << "Atmospheric dimension must be 1D, but atmosphere_dim is "
       << atmosphere_dim << ".";
    throw runtime_error(os.str());
  }

  abs_p = p_grid;
  abs_t = t_field(joker, 0, 0);
  abs_vmrs = vmr_field(joker, joker, 0, 0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_coefCalcFromXsec(  // WS Output:
    Matrix& abs_coef,
    Matrix& src_coef,
    ArrayOfMatrix& dabs_coef_dx,
    ArrayOfMatrix& dsrc_coef_dx,
    ArrayOfMatrix& abs_coef_per_species,
    ArrayOfMatrix& src_coef_per_species,
    // WS Input:
    const ArrayOfMatrix& abs_xsec_per_species,
    const ArrayOfMatrix& src_xsec_per_species,
    const ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    const ArrayOfArrayOfMatrix& dsrc_xsec_per_species_dx,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Matrix& abs_vmrs,
    const Vector& abs_p,
    const Vector& abs_t,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Check that abs_vmrs and abs_xsec_per_species really have compatible
  // dimensions. In abs_vmrs there should be one row for each tg:
  if (abs_vmrs.nrows() != abs_xsec_per_species.nelem()) {
    ostringstream os;
    os << "Variable abs_vmrs must have compatible dimension to abs_xsec_per_species.\n"
       << "abs_vmrs.nrows() = " << abs_vmrs.nrows() << "\n"
       << "abs_xsec_per_species.nelem() = " << abs_xsec_per_species.nelem();
    throw runtime_error(os.str());
  }

  // Check that number of altitudes are compatible. We only check the
  // first element, this is possilble because within arts all elements
  // are on the same altitude grid.
  if (abs_vmrs.ncols() != abs_xsec_per_species[0].ncols()) {
    ostringstream os;
    os << "Variable abs_vmrs must have same numbers of altitudes as abs_xsec_per_species.\n"
       << "abs_vmrs.ncols() = " << abs_vmrs.ncols() << "\n"
       << "abs_xsec_per_species[0].ncols() = "
       << abs_xsec_per_species[0].ncols();
    throw runtime_error(os.str());
  }

  // Check dimensions of abs_p and abs_t:
  chk_size("abs_p", abs_p, abs_vmrs.ncols());
  chk_size("abs_t", abs_t, abs_vmrs.ncols());

  // Will these calculations deal with nlte?
  const bool do_src =
      (src_xsec_per_species.nelem() == abs_xsec_per_species.nelem())
          ? not src_xsec_per_species[0].empty() ? true : false
          : false;
  const Index src_rows = do_src ? src_xsec_per_species[0].nrows() : 0;

  // Initialize abs_coef and abs_coef_per_species. The array dimension of abs_coef_per_species
  // is the same as that of abs_xsec_per_species. The dimension of abs_coef should
  // be equal to one of the abs_xsec_per_species enries.

  abs_coef.resize(abs_xsec_per_species[0].nrows(),
                  abs_xsec_per_species[0].ncols());
  abs_coef = 0;
  src_coef.resize(src_rows, src_xsec_per_species[0].ncols());
  src_coef = 0;

  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);
  dabs_coef_dx.resize(jacobian_quantities_position.nelem());
  dsrc_coef_dx.resize(do_src ? jacobian_quantities_position.nelem() : 0);

  for (Index ii = 0; ii < jacobian_quantities_position.nelem(); ii++) {
    dabs_coef_dx[ii].resize(abs_xsec_per_species[0].nrows(),
                            abs_xsec_per_species[0].ncols());
    dabs_coef_dx[ii] = 0.0;
    if (do_src) {
      dsrc_coef_dx[ii].resize(src_rows, src_xsec_per_species[0].ncols());
      dsrc_coef_dx[ii] = 0.0;
    }
  }

  abs_coef_per_species.resize(abs_xsec_per_species.nelem());
  src_coef_per_species.resize(src_xsec_per_species.nelem());

  out3
      << "  Computing abs_coef and abs_coef_per_species from abs_xsec_per_species.\n";
  // Loop through all tag groups
  for (Index i = 0; i < abs_xsec_per_species.nelem(); ++i) {
    out3 << "  Tag group " << i << "\n";

    // Make this element of abs_xsec_per_species the right size:
    abs_coef_per_species[i].resize(abs_xsec_per_species[i].nrows(),
                                   abs_xsec_per_species[i].ncols());
    abs_coef_per_species[i] = 0;  // Initialize all elements to 0.

    if (do_src) {
      src_coef_per_species[i].resize(src_rows, src_xsec_per_species[i].ncols());
      src_coef_per_species[i] = 0;  // Initialize all elements to 0.
    }

    // Loop through all altitudes
    for (Index j = 0; j < abs_xsec_per_species[i].ncols(); j++) {
      // Calculate total number density from pressure and temperature.
      const Numeric n = number_density(abs_p[j], abs_t[j]);
      const Numeric dn_dT = dnumber_density_dt(abs_p[j], abs_t[j]);
      // Wasted calculations when Jacobians are not calculated...
      // Though this is called seldom enough that it this fine?  value is -1/t*n

      // Loop through all frequencies
      for (Index k = 0; k < abs_xsec_per_species[i].nrows(); k++) {
        abs_coef_per_species[i](k, j) =
            abs_xsec_per_species[i](k, j) * n * abs_vmrs(i, j);
        if (do_src)
          src_coef_per_species[i](k, j) =
              src_xsec_per_species[i](k, j) * n * abs_vmrs(i, j);

        for (Index iq = 0; iq < jacobian_quantities_position.nelem(); iq++) {
          if (jacobian_quantities[jacobian_quantities_position[iq]] ==
              Jacobian::Atm::Temperature) {
            dabs_coef_dx[iq](k, j) +=
                (dabs_xsec_per_species_dx[i][iq](k, j) * n +
                 abs_xsec_per_species[i](k, j) * dn_dT) *
                abs_vmrs(i, j);
            if (do_src)
              dsrc_coef_dx[iq](k, j) +=
                  (dsrc_xsec_per_species_dx[i][iq](k, j) * n +
                   src_xsec_per_species[i](k, j) * dn_dT) *
                  abs_vmrs(i, j);
          } else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                     Jacobian::Special::VMR) {
            bool seco = false, main = false;
            for (const auto& s : abs_species[i]) {
              if (species_match(
                      jacobian_quantities[jacobian_quantities_position[iq]],
                      s.BathSpecies()) or
                  s.Type() not_eq SpeciesTag::TYPE_CIA)
                seco = true;
              if (species_iso_match(
                      jacobian_quantities[jacobian_quantities_position[iq]],
                      s.Species(),
                      s.Isotopologue()))
                main = true;
            }
            if (main and seco) {
              dabs_coef_dx[iq](k, j) +=
                  (dabs_xsec_per_species_dx[i][iq](k, j) * abs_vmrs(i, j) +
                   abs_xsec_per_species[i](k, j)) *
                  n;
              if (do_src)
                dsrc_coef_dx[iq](k, j) +=
                    (dsrc_xsec_per_species_dx[i][iq](k, j) * abs_vmrs(i, j) +
                     src_xsec_per_species[i](k, j)) *
                    n;
            } else if (main) {
              dabs_coef_dx[iq](k, j) += abs_xsec_per_species[i](k, j) * n;
              if (do_src)
                dsrc_coef_dx[iq](k, j) += src_xsec_per_species[i](k, j) * n;
            } else if (seco) {
              dabs_coef_dx[iq](k, j) +=
                  dabs_xsec_per_species_dx[i][iq](k, j) * abs_vmrs(i, j) * n;
              if (do_src)
                dsrc_coef_dx[iq](k, j) +=
                    dsrc_xsec_per_species_dx[i][iq](k, j) * abs_vmrs(i, j) * n;
            }
          } else {
            dabs_coef_dx[iq](k, j) +=
                dabs_xsec_per_species_dx[i][iq](k, j) * n * abs_vmrs(i, j);
            if (do_src)
              dsrc_coef_dx[iq](k, j) +=
                  dsrc_xsec_per_species_dx[i][iq](k, j) * n * abs_vmrs(i, j);
          }
        }
      }
    }

    // Add up to the total absorption:
    abs_coef += abs_coef_per_species[i];  // In Matpack you can use the +=
        // operator to do elementwise addition.
    if (do_src) src_coef += src_coef_per_species[i];
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesInit(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfMatrix& src_xsec_per_species,
    ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    ArrayOfArrayOfMatrix& dsrc_xsec_per_species_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& tgs,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Index& abs_xsec_agenda_checked,
    const Index& nlte_do,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  if (!abs_xsec_agenda_checked)
    throw runtime_error(
        "You must call *abs_xsec_agenda_checkedCalc* before calling this method.");

  // We need to check that abs_species_active doesn't have more elements than
  // abs_species (abs_xsec_agenda_checkedCalc doesn't know abs_species_active.
  // Usually we come here through an agenda call, where abs_species_active has
  // been properly created somewhere internally. But we might get here by
  // direct call, and then need to be safe!).
  if (tgs.nelem() < abs_species_active.nelem()) {
    ostringstream os;
    os << "abs_species_active (n=" << abs_species_active.nelem()
       << ") not allowed to have more elements than abs_species (n="
       << tgs.nelem() << ")!\n";
    throw runtime_error(os.str());
  }

  // Initialize abs_xsec_per_species. The array dimension of abs_xsec_per_species
  // is the same as that of abs_lines_per_species.
  abs_xsec_per_species.resize(tgs.nelem());
  src_xsec_per_species.resize(tgs.nelem());
  
  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);

  dabs_xsec_per_species_dx.resize(tgs.nelem());
  dsrc_xsec_per_species_dx.resize(tgs.nelem());

  // Loop abs_xsec_per_species and make each matrix the right size,
  // initializing to zero.
  // But skip inactive species, loop only over the active ones.
  for (Index ii = 0; ii < abs_species_active.nelem(); ++ii) {
    const Index i = abs_species_active[ii];
    // Check that abs_species_active index is not higher than the number
    // of species
    if (i >= tgs.nelem()) {
      ostringstream os;
      os << "*abs_species_active* contains an invalid species index.\n"
         << "Species index must be between 0 and " << tgs.nelem() - 1;
      throw std::runtime_error(os.str());
    }
    // Make this element of abs_xsec_per_species the right size:
    abs_xsec_per_species[i].resize(f_grid.nelem(), abs_p.nelem());
    abs_xsec_per_species[i] = 0;  // Matpack can set all elements like this.
    if (nlte_do) {
      src_xsec_per_species[i].resize(f_grid.nelem(), abs_p.nelem());
      src_xsec_per_species[i] = 0;
    } else {
      src_xsec_per_species[i].resize(0, 0);
    }

    dabs_xsec_per_species_dx[ii] =
        ArrayOfMatrix(jacobian_quantities_position.nelem(),
                      Matrix(f_grid.nelem(), abs_p.nelem(), 0.0));
    if (nlte_do)
      dsrc_xsec_per_species_dx[ii] =
          ArrayOfMatrix(jacobian_quantities_position.nelem(),
                        Matrix(f_grid.nelem(), abs_p.nelem(), 0.0));
    else
      dsrc_xsec_per_species_dx[ii] =
          ArrayOfMatrix(jacobian_quantities_position.nelem(),
                        Matrix(0, 0));
  }

  ostringstream os;
  os << "  Initialized abs_xsec_per_species.\n"
     << "  Number of frequencies        : " << f_grid.nelem() << "\n"
     << "  Number of pressure levels    : " << abs_p.nelem() << "\n";
  out3 << os.str();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddConts(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& tgs,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const Matrix& abs_vmrs,
    const ArrayOfString& abs_cont_names,
    const ArrayOfVector& abs_cont_parameters,
    const ArrayOfString& abs_cont_models,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Needed for some continua, and set here from abs_vmrs:
  Vector abs_h2o, abs_n2, abs_o2;

  // Check that all paramters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs = tgs.nelem();
    const Index n_xsec = abs_xsec_per_species.nelem();
    const Index n_vmrs = abs_vmrs.nrows();

    if (n_tgs != n_xsec || n_tgs != n_vmrs) {
      ostringstream os;
      os << "The following variables must all have the same dimension:\n"
         << "tgs:          " << tgs.nelem() << "\n"
         << "abs_xsec_per_species:  " << abs_xsec_per_species.nelem() << "\n"
         << "abs_vmrs.nrows():      " << abs_vmrs.nrows();
      throw runtime_error(os.str());
    }
  }

  // Jacobian overhead START
  /* NOTE: if any of the functions inside continuum tags could 
           be made to give partial derivatives, then that would 
           speed things up.  Also be aware that line specific
           parameters cannot be retrieved while using these 
           models. */
  const bool do_jac = supports_continuum(
      jacobian_quantities);  // Throws runtime error if line parameters are wanted since we cannot know if the line is in the Continuum...
  const bool do_freq_jac = do_frequency_jacobian(jacobian_quantities);
  const bool do_temp_jac = do_temperature_jacobian(jacobian_quantities);
  Vector dfreq, dabs_t;
  const Numeric df = frequency_perturbation(jacobian_quantities);
  const Numeric dt = temperature_perturbation(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);

  if (do_freq_jac) {
    dfreq.resize(f_grid.nelem());
    for (Index iv = 0; iv < f_grid.nelem(); iv++) dfreq[iv] = f_grid[iv] + df;
  }
  if (do_temp_jac) {
    dabs_t.resize(abs_t.nelem());
    for (Index it = 0; it < abs_t.nelem(); it++) dabs_t[it] = abs_t[it] + dt;
  }

  Matrix jacs_df, jacs_dt, normal;
  if (do_jac) {
    if (do_freq_jac) jacs_df.resize(f_grid.nelem(), abs_p.nelem());
    if (do_temp_jac) jacs_dt.resize(f_grid.nelem(), abs_p.nelem());
    normal.resize(f_grid.nelem(), abs_p.nelem());
  }
  // Jacobian overhead END

  // Check, that dimensions of abs_cont_names and
  // abs_cont_parameters are consistent...
  if (abs_cont_names.nelem() != abs_cont_parameters.nelem()) {
    ostringstream os;

    for (Index i = 0; i < abs_cont_names.nelem(); ++i)
      os << "abs_xsec_per_speciesAddConts: " << i
         << " name : " << abs_cont_names[i] << "\n";

    for (Index i = 0; i < abs_cont_parameters.nelem(); ++i)
      os << "abs_xsec_per_speciesAddConts: " << i
         << " param: " << abs_cont_parameters[i] << "\n";

    for (Index i = 0; i < abs_cont_models.nelem(); ++i)
      os << "abs_xsec_per_speciesAddConts: " << i
         << " option: " << abs_cont_models[i] << "\n";

    os << "The following variables must have the same dimension:\n"
       << "abs_cont_names:      " << abs_cont_names.nelem() << "\n"
       << "abs_cont_parameters: " << abs_cont_parameters.nelem();

    throw runtime_error(os.str());
  }

  // Check that abs_p, abs_t, and abs_vmrs have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error.

  if (abs_t.nelem() != abs_p.nelem()) {
    ostringstream os;
    os << "Variable abs_t must have the same dimension as abs_p.\n"
       << "abs_t.nelem() = " << abs_t.nelem() << '\n'
       << "abs_p.nelem() = " << abs_p.nelem();
    throw runtime_error(os.str());
  }

  if (abs_vmrs.ncols() != abs_p.nelem()) {
    ostringstream os;
    os << "Variable dimension abs_vmrs.ncols() must\n"
       << "be the same as abs_p.nelem().\n"
       << "abs_vmrs.ncols() = " << abs_vmrs.ncols() << '\n'
       << "abs_p.nelem() = " << abs_p.nelem();
    throw runtime_error(os.str());
  }

  // We set abs_h2o, abs_n2, and abs_o2 later, because we only want to
  // do it if the parameters are really needed.

  out3 << "  Calculating continuum spectra.\n";

  // Loop tag groups:
  for (Index ii = 0; ii < abs_species_active.nelem(); ++ii) {
    const Index i = abs_species_active[ii];

    using global_data::species_data;

    // Go through the tags in the current tag group to see if they
    // are continuum tags:
    for (Index s = 0; s < tgs[i].nelem(); ++s) {
      // Continuum tags in the sense that we talk about here
      // (including complete absorption models) are marked by a special type.
      if (tgs[i][s].Type() == SpeciesTag::TYPE_PREDEF) {
        // We have identified a continuum tag!

        // Get only the continuum name. The full tag name is something like:
        // H2O-HITRAN96Self-*-*. We want only the `H2O-HITRAN96Self' part:
        const String name = species_data[tgs[i][s].Species()].Name() + "-" +
                            species_data[tgs[i][s].Species()]
                                .Isotopologue()[tgs[i][s].Isotopologue()]
                                .Name();

        if (name == "O2-MPM2020") continue;
                                
        // Check, if we have parameters for this model. For
        // this, the model name must be listed in
        // abs_cont_names.
        const Index n =
            find(abs_cont_names.begin(), abs_cont_names.end(), name) -
            abs_cont_names.begin();

        // n==abs_cont_names.nelem() indicates that
        // the name was not found.
        if (n == abs_cont_names.nelem()) {
          ostringstream os;
          os << "Cannot find model " << name << " in abs_cont_names.";
          throw runtime_error(os.str());
        }

        // Ok, the tag specifies a valid continuum model and
        // we have continuum parameters.

        if (out3.sufficient_priority()) {
          ostringstream os;
          os << "  Adding " << name << " to tag group " << i << ".\n";
          out3 << os.str();
        }

        // find the options for this continuum tag from the input array
        // of options. The actual field of the array is n:
        const String ContOption = abs_cont_models[n];

        // Set abs_h2o, abs_n2, and abs_o2 from the first matching species.
        set_vmr_from_first_species(abs_h2o, "H2O", tgs, abs_vmrs);
        set_vmr_from_first_species(abs_n2, "N2", tgs, abs_vmrs);
        set_vmr_from_first_species(abs_o2, "O2", tgs, abs_vmrs);

        // Add the continuum for this tag. The parameters in
        // this call should be clear. The vmr is in
        // abs_vmrs(i,Range(joker)). The other vmr variables,
        // abs_h2o, abs_n2, and abs_o2 contains the real vmr of H2O,
        // N2, nad O2, which are needed as additional information for
        // certain continua:
        // abs_h2o for
        //   O2-PWR88, O2-PWR93, O2-PWR98,
        //   O2-MPM85, O2-MPM87, O2-MPM89, O2-MPM92, O2-MPM93,
        //   O2-TRE05,
        //   O2-SelfContStandardType, O2-SelfContMPM93, O2-SelfContPWR93,
        //   N2-SelfContMPM93, N2-DryContATM01,
        //   N2-CIArotCKDMT252, N2-CIAfunCKDMT252
        // abs_n2 for
        //   H2O-SelfContCKD24, H2O-ForeignContCKD24,
        //   O2-v0v0CKDMT100,
        //   CO2-ForeignContPWR93, CO2-ForeignContHo66
        // abs_o2 for
        //   N2-CIArotCKDMT252, N2-CIAfunCKDMT252
        if (!do_jac)
          xsec_continuum_tag(abs_xsec_per_species[i],
                             name,
                             abs_cont_parameters[n],
                             abs_cont_models[n],
                             f_grid,
                             abs_p,
                             abs_t,
                             abs_n2,
                             abs_h2o,
                             abs_o2,
                             abs_vmrs(i, Range(joker)),
                             verbosity);
        else  // The Jacobian block
        {
          // Needs a reseted block here...
          for (Index iv = 0; iv < f_grid.nelem(); iv++) {
            for (Index ip = 0; ip < abs_p.nelem(); ip++) {
              if (do_freq_jac) jacs_df(iv, ip) = 0.0;
              if (do_temp_jac) jacs_dt(iv, ip) = 0.0;
              normal(iv, ip) = 0.0;
            }
          }

          // Normal calculations
          xsec_continuum_tag(normal,
                             name,
                             abs_cont_parameters[n],
                             abs_cont_models[n],
                             f_grid,
                             abs_p,
                             abs_t,
                             abs_n2,
                             abs_h2o,
                             abs_o2,
                             abs_vmrs(i, Range(joker)),
                             verbosity);

          // Frequency calculations
          if (do_freq_jac)
            xsec_continuum_tag(jacs_df,
                               name,
                               abs_cont_parameters[n],
                               abs_cont_models[n],
                               dfreq,
                               abs_p,
                               abs_t,
                               abs_n2,
                               abs_h2o,
                               abs_o2,
                               abs_vmrs(i, Range(joker)),
                               verbosity);

          //Temperature calculations
          if (do_temp_jac)
            xsec_continuum_tag(jacs_dt,
                               name,
                               abs_cont_parameters[n],
                               abs_cont_models[n],
                               f_grid,
                               abs_p,
                               dabs_t,
                               abs_n2,
                               abs_h2o,
                               abs_o2,
                               abs_vmrs(i, Range(joker)),
                               verbosity);
          for (Index iv = 0; iv < f_grid.nelem(); iv++) {
            for (Index ip = 0; ip < abs_p.nelem(); ip++) {
              abs_xsec_per_species[i](iv, ip) += normal(iv, ip);
              for (Index iq = 0; iq < jacobian_quantities_position.nelem();
                   iq++) {
                if (is_frequency_parameter(
                        jacobian_quantities[jacobian_quantities_position[iq]]))
                  dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                      (jacs_df(iv, ip) - normal(iv, ip)) * (1. / df);
                else if (jacobian_quantities
                             [jacobian_quantities_position[iq]] ==
                         Jacobian::Atm::Temperature)
                  dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                      (jacs_dt(iv, ip) - normal(iv, ip)) * (1. / dt);
              }
            }
          }
        }
        // Calling this function with a row of Matrix abs_vmrs
        // is possible because it uses Views.
      }
    }
  }
}

//======================================================================
//             Methods related to continua
//======================================================================

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cont_descriptionInit(  // WS Output:
    ArrayOfString& abs_cont_names,
    ArrayOfString& abs_cont_options,
    ArrayOfVector& abs_cont_parameters,
    const Verbosity& verbosity) {
  CREATE_OUT2;

  abs_cont_names.resize(0);
  abs_cont_options.resize(0);
  abs_cont_parameters.resize(0);
  out2 << "  Initialized abs_cont_names \n"
          "              abs_cont_models\n"
          "              abs_cont_parameters.\n";
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_cont_descriptionAppend(  // WS Output:
    ArrayOfString& abs_cont_names,
    ArrayOfString& abs_cont_models,
    ArrayOfVector& abs_cont_parameters,
    // Control Parameters:
    const String& tagname,
    const String& model,
    const Vector& userparameters,
    const Verbosity&) {
  // First we have to check that name matches a continuum species tag.
  check_continuum_model(tagname);

  //cout << "   + tagname:    " << tagname << "\n";
  //cout << "   + model:      " << model << "\n";
  //cout << "   + parameters: " << userparameters << "\n";

  // Add name and parameters to the apropriate variables:
  abs_cont_names.push_back(tagname);
  abs_cont_models.push_back(model);
  abs_cont_parameters.push_back(userparameters);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void nlte_sourceFromTemperatureAndSrcCoefPerSpecies(  // WS Output:
    ArrayOfStokesVector& nlte_source,
    ArrayOfStokesVector& dnlte_dx_source,
    ArrayOfStokesVector& nlte_dsource_dx,
    // WS Input:
    const ArrayOfMatrix& src_coef_per_species,
    const ArrayOfMatrix& dsrc_coef_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const Numeric& rtp_temperature,
    const Verbosity&) {
  // nlte_source has format
  // [ abs_species, f_grid, stokes_dim ].
  // src_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].

  Index n_species = src_coef_per_species.nelem();  // # species

  if (not n_species) {
    ostringstream os;
    os << "Must have at least one species.";
    throw runtime_error(os.str());
  }

  Index n_f = src_coef_per_species[0].nrows();  // # frequencies

  // # pressures must be 1:
  if (1 not_eq src_coef_per_species[0].ncols()) {
    ostringstream os;
    os << "Must have exactly one pressure.";
    throw runtime_error(os.str());
  }

  // Check species dimension of propmat_clearsky
  if (nlte_source.nelem() not_eq n_species) {
    ostringstream os;
    os << "Species dimension of propmat_clearsky does not\n"
       << "match src_coef_per_species.";
    throw std::runtime_error(os.str());
  }

  // Check frequency dimension of propmat_clearsky
  if (nlte_source[0].NumberOfFrequencies() not_eq n_f) {
    ostringstream os;
    os << "Frequency dimension of propmat_clearsky does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error(os.str());
  }

  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);
  Vector B(n_f);

  for (Index iv = 0; iv < n_f; iv++)
    B[iv] = planck(f_grid[iv], rtp_temperature);

  StokesVector sv(n_f, nlte_source[0].StokesDimensions());
  for (Index si = 0; si < n_species; ++si) {
    sv.Kjj() = src_coef_per_species[si](joker, 0);
    sv *= B;
    nlte_source[si].Kjj() += sv.Kjj();
  }

  // Jacobian
  for (Index ii = 0; ii < jacobian_quantities_position.nelem(); ii++) {
    if (jacobian_quantities[jacobian_quantities_position[ii]] ==
        Jacobian::Atm::Temperature) {
      Vector dB(n_f);
      for (Index iv = 0; iv < n_f; iv++)
        dB[iv] = dplanck_dt(f_grid[iv], rtp_temperature);

      for (Index si = 0; si < n_species; ++si) {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        nlte_dsource_dx[ii].Kjj() += sv.Kjj();
      }

      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_dx_source[ii].Kjj() += sv.Kjj();
    } else if (is_frequency_parameter(
                   jacobian_quantities[jacobian_quantities_position[ii]])) {
      Vector dB(n_f);
      for (Index iv = 0; iv < n_f; iv++)
        dB[iv] = dplanck_df(f_grid[iv], rtp_temperature);

      for (Index si = 0; si < n_species; ++si) {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        nlte_dsource_dx[ii].Kjj() += sv.Kjj();
      }

      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_dx_source[ii].Kjj() += sv.Kjj();
    } else {
      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_dx_source[ii].Kjj() += sv.Kjj();
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFromAbsCoefPerSpecies(  // WS Output:
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const ArrayOfMatrix& abs_coef_per_species,
    const ArrayOfMatrix& dabs_coef_dx,
    const Verbosity&) {
  // propmat_clearsky has format
  // [ abs_species, f_grid, stokes_dim, stokes_dim ].
  // abs_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].

  Index n_species = abs_coef_per_species.nelem();  // # species

  if (0 == n_species) {
    ostringstream os;
    os << "Must have at least one species.";
    throw runtime_error(os.str());
  }

  Index n_f = abs_coef_per_species[0].nrows();  // # frequencies

  // # pressures must be 1:
  if (1 not_eq abs_coef_per_species[0].ncols()) {
    ostringstream os;
    os << "Must have exactly one pressure.";
    throw runtime_error(os.str());
  }

  // Check species dimension of propmat_clearsky
  if (propmat_clearsky.nelem() not_eq n_species) {
    ostringstream os;
    os << "Species dimension of propmat_clearsky does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error(os.str());
  }

  // Check frequency dimension of propmat_clearsky
  if (propmat_clearsky[0].NumberOfFrequencies() not_eq n_f) {
    ostringstream os;
    os << "Frequency dimension of propmat_clearsky does not\n"
       << "match abs_coef_per_species.";
    throw runtime_error(os.str());
  }

  // Loop species and stokes dimensions, and add to propmat_clearsky:
  for (Index si = 0; si < n_species; ++si)
    propmat_clearsky[si].Kjj() += abs_coef_per_species[si](joker, 0);

  for (Index iqn = 0; iqn < dabs_coef_dx.nelem(); iqn++) {
    if (dabs_coef_dx[iqn].nrows() == n_f) {
      if (dabs_coef_dx[iqn].ncols() == 1) {
        dpropmat_clearsky_dx[iqn].Kjj() += dabs_coef_dx[iqn](joker, 0);
      } else
        throw std::runtime_error("Must have exactly one pressure.");
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyInit(  //WS Output
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfStokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_dx_source,
    ArrayOfStokesVector& nlte_dsource_dx,
    //WS Input
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const Index& stokes_dim,
    const Index& propmat_clearsky_agenda_checked,
    const Index& nlte_do,
    const Verbosity&) {
  if (!propmat_clearsky_agenda_checked)
    throw runtime_error(
        "You must call *propmat_clearsky_agenda_checkedCalc* before calling this method.");

  Index nf = f_grid.nelem();

  if (not abs_species.nelem())
    throw std::runtime_error("abs_species.nelem() = 0");

  if (not nf) throw runtime_error("nf = 0");

  if (not stokes_dim) throw runtime_error("stokes_dim = 0");

  const Index nq = equivalent_propmattype_indexes(jacobian_quantities).nelem();

  propmat_clearsky = ArrayOfPropagationMatrix(
      abs_species.nelem(), PropagationMatrix(nf, stokes_dim));
  dpropmat_clearsky_dx =
      ArrayOfPropagationMatrix(nq, PropagationMatrix(nf, stokes_dim));

  nlte_source = nlte_do ? ArrayOfStokesVector(abs_species.nelem(),
                                              StokesVector(nf, stokes_dim))
                        : ArrayOfStokesVector(0);
  dnlte_dx_source = nlte_do
                        ? ArrayOfStokesVector(nq, StokesVector(nf, stokes_dim))
                        : ArrayOfStokesVector(0);
  nlte_dsource_dx = nlte_do
                        ? ArrayOfStokesVector(nq, StokesVector(nf, stokes_dim))
                        : ArrayOfStokesVector(0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFaraday(
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& rtp_vmr,
    const Vector& rtp_los,
    const Vector& rtp_mag,
    const Verbosity&) {
  // All the physical constants joined into one static constant:
  // (abs as e defined as negative)
  static const Numeric FRconst =
      abs(ELECTRON_CHARGE * ELECTRON_CHARGE * ELECTRON_CHARGE /
          (8 * PI * PI * SPEED_OF_LIGHT * VACUUM_PERMITTIVITY * ELECTRON_MASS *
           ELECTRON_MASS));

  if (stokes_dim < 3)
    throw runtime_error(
        "To include Faraday rotation, stokes_dim >= 3 is required.");

  if (atmosphere_dim == 1 && rtp_los.nelem() < 1) {
    ostringstream os;
    os << "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
       << "(at least zenith angle component for atmosphere_dim==1),\n"
       << "but it is not.\n";
    throw runtime_error(os.str());
  } else if (atmosphere_dim > 1 && rtp_los.nelem() < 2) {
    ostringstream os;
    os << "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
       << "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
       << "but it is not.\n";
    throw runtime_error(os.str());
  }

  const bool do_jac = supports_faraday(jacobian_quantities);
  const bool do_magn_jac = do_magnetic_jacobian(jacobian_quantities);
  const Numeric dmag = magnetic_field_perturbation(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);

  Index ife = -1;
  for (Index sp = 0; sp < abs_species.nelem() && ife < 0; sp++) {
    if (abs_species[sp][0].Type() == SpeciesTag::TYPE_FREE_ELECTRONS) {
      ife = sp;
    }
  }

  if (ife < 0) {
    throw runtime_error(
        "Free electrons not found in *abs_species* and "
        "Faraday rotation can not be calculated.");
  } else {
    const Numeric ne = rtp_vmr[ife];

    if (ne != 0 && (rtp_mag[0] != 0 || rtp_mag[1] != 0 || rtp_mag[2] != 0)) {
      // Include remaining terms, beside /f^2
      const Numeric c1 =
          2 * FRconst * ne *
          dotprod_with_los(
              rtp_los, rtp_mag[0], rtp_mag[1], rtp_mag[2], atmosphere_dim);

      Numeric dc1_u = 0.0, dc1_v = 0.0, dc1_w = 0.0;
      if (do_magn_jac) {
        dc1_u = (2 * FRconst * ne *
                     dotprod_with_los(rtp_los,
                                      rtp_mag[0] + dmag,
                                      rtp_mag[1],
                                      rtp_mag[2],
                                      atmosphere_dim) -
                 c1) /
                dmag;
        dc1_v = (2 * FRconst * ne *
                     dotprod_with_los(rtp_los,
                                      rtp_mag[0],
                                      rtp_mag[1] + dmag,
                                      rtp_mag[2],
                                      atmosphere_dim) -
                 c1) /
                dmag;
        dc1_w = (2 * FRconst * ne *
                     dotprod_with_los(rtp_los,
                                      rtp_mag[0],
                                      rtp_mag[1],
                                      rtp_mag[2] + dmag,
                                      atmosphere_dim) -
                 c1) /
                dmag;
      }

      if (not do_jac) {
        for (Index iv = 0; iv < f_grid.nelem(); iv++) {
          const Numeric r = c1 / (f_grid[iv] * f_grid[iv]);
          propmat_clearsky[ife].SetFaraday(r, iv);
        }
      } else {
        for (Index iv = 0; iv < f_grid.nelem(); iv++) {
          const Numeric f2 = f_grid[iv] * f_grid[iv];
          const Numeric r = c1 / f2;
          propmat_clearsky[ife].SetFaraday(r, iv);

          // The Jacobian loop
          for (Index iq = 0; iq < jacobian_quantities_position.nelem(); iq++) {
            if (is_frequency_parameter(
                    jacobian_quantities[jacobian_quantities_position[iq]]))
              dpropmat_clearsky_dx[iq].AddFaraday(-2.0 * r / f_grid[iv], iv);
            else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                     Jacobian::Atm::MagneticU)
              dpropmat_clearsky_dx[iq].AddFaraday(dc1_u / f2, iv);
            else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                     Jacobian::Atm::MagneticV)
              dpropmat_clearsky_dx[iq].AddFaraday(dc1_v / f2, iv);
            else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                     Jacobian::Atm::MagneticW)
              dpropmat_clearsky_dx[iq].AddFaraday(dc1_w / f2, iv);
            else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                     Jacobian::Atm::Electrons)
              dpropmat_clearsky_dx[iq].AddFaraday(r / ne, iv);
          }
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddParticles(
    // WS Output:
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const Index& stokes_dim,
    const Index& atmosphere_dim,
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& rtp_vmr,
    const Vector& rtp_los,
    const Numeric& rtp_temperature,
    const ArrayOfArrayOfSingleScatteringData& scat_data,
    const Index& scat_data_checked,
    const Index& use_abs_as_ext,
    // Verbosity object:
    const Verbosity& verbosity) {
  CREATE_OUT1;

  // (i)yCalc only checks scat_data_checked if cloudbox is on. It is off here,
  // though, i.e. we need to check it here explicitly. (Also, cloudboxOff sets
  // scat_data_checked=0 as it does not check it and as we ususally don't need
  // scat_data for clearsky cases, hence don't want to check them by
  // scat_data_checkedCalc in that case. This approach seems to be the more
  // handy compared to cloudboxOff setting scat_data_checked=1 without checking
  // it assuming we won't use it anyways.)
  if (scat_data_checked != 1)
    throw runtime_error(
        "The scat_data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).");

  const Index ns = TotalNumberOfElements(scat_data);
  Index np = 0;
  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    if (abs_species[sp][0].Type() == SpeciesTag::TYPE_PARTICLES) {
      np++;
    }
  }

  if (np == 0) {
    ostringstream os;
    os << "For applying propmat_clearskyAddParticles, *abs_species* needs to"
       << "contain species 'particles', but it does not.\n";
    throw runtime_error(os.str());
  }

  if (ns != np) {
    ostringstream os;
    os << "Number of 'particles' entries in abs_species and of elements in\n"
       << "*scat_data* needs to be identical. But you have " << np
       << " 'particles' entries\n"
       << "and " << ns << " *scat_data* elements.\n";
    throw runtime_error(os.str());
  }

  if (atmosphere_dim == 1 && rtp_los.nelem() < 1) {
    ostringstream os;
    os << "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
       << "(at least zenith angle component for atmosphere_dim==1),\n"
       << "but it is not.\n";
    throw runtime_error(os.str());
  } else if (atmosphere_dim > 1 && rtp_los.nelem() < 2) {
    ostringstream os;
    os << "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
       << "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
       << "but it is not.\n";
    throw runtime_error(os.str());
  }

  // Use for rescaling vmr of particulates
  Numeric rtp_vmr_sum = 0.0;

  // Tests and setup partial derivatives
  const bool do_jac = supports_particles(jacobian_quantities);
  const bool do_jac_temperature = do_temperature_jacobian(jacobian_quantities);
  const bool do_jac_frequencies = do_frequency_jacobian(jacobian_quantities);
  const ArrayOfIndex jacobian_quantities_position =
      equivalent_propmattype_indexes(jacobian_quantities);
  const Numeric dT = temperature_perturbation(jacobian_quantities);

  const Index na = abs_species.nelem();
  Vector rtp_los_back;
  mirror_los(rtp_los_back, rtp_los, atmosphere_dim);

  // 170918 JM: along with transition to use of new-type (aka
  // pre-f_grid-interpolated) scat_data, freq perturbation switched off. Typical
  // clear-sky freq perturbations yield insignificant effects in particle
  // properties. Hence, this feature is neglected here.
  if (do_jac_frequencies) {
    out1 << "WARNING:\n"
         << "Frequency perturbation not available for absorbing particles.\n";
  }

  // creating temporary output containers
  ArrayOfArrayOfTensor5 ext_mat_Nse;
  ArrayOfArrayOfTensor4 abs_vec_Nse;
  ArrayOfArrayOfIndex ptypes_Nse;
  Matrix t_ok;

  // preparing input in format needed
  Vector T_array;
  if (do_jac_temperature) {
    T_array.resize(2);
    T_array = rtp_temperature;
    T_array[1] += dT;
  } else {
    T_array.resize(1);
    T_array = rtp_temperature;
  }
  Matrix dir_array(1, 2);
  dir_array(0, joker) = rtp_los_back;

  // ext/abs per scat element for all freqs at once
  opt_prop_NScatElems(ext_mat_Nse,
                      abs_vec_Nse,
                      ptypes_Nse,
                      t_ok,
                      scat_data,
                      stokes_dim,
                      T_array,
                      dir_array,
                      -1);

  const Index nf = abs_vec_Nse[0][0].nbooks();
  Tensor3 tmp(nf, stokes_dim, stokes_dim);

  // loop over the scat_data and link them with correct vmr_field entry according
  // to the position of the particle type entries in abs_species.
  Index sp = 0;
  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
      // forward to next particle entry in abs_species
      while (sp < na && abs_species[sp][0].Type() != SpeciesTag::TYPE_PARTICLES)
        sp++;

      // running beyond number of abs_species entries when looking for next
      // particle entry. shouldn't happen, though.
      assert(sp < na);
      if (rtp_vmr[sp] < 0.) {
        ostringstream os;
        os << "Negative absorbing particle 'vmr' (aka number density)"
           << " encountered:\n"
           << "scat species #" << i_ss << ", scat elem #" << i_se
           << " (vmr_field entry #" << sp << ")\n";
        throw runtime_error(os.str());
      } else if (rtp_vmr[sp] > 0.) {
        if (t_ok(i_se_flat, 0) < 0.) {
          ostringstream os;
          os << "Temperature interpolation error:\n"
             << "scat species #" << i_ss << ", scat elem #" << i_se << "\n";
          throw runtime_error(os.str());
        } else {
          if (use_abs_as_ext) {
            if (nf > 1)
              for (Index iv = 0; iv < f_grid.nelem(); iv++)
                propmat_clearsky[sp].AddAbsorptionVectorAtPosition(
                    abs_vec_Nse[i_ss][i_se](iv, 0, 0, joker), iv);
            else
              for (Index iv = 0; iv < f_grid.nelem(); iv++)
                propmat_clearsky[sp].AddAbsorptionVectorAtPosition(
                    abs_vec_Nse[i_ss][i_se](0, 0, 0, joker), iv);
          } else {
            if (nf > 1)
              for (Index iv = 0; iv < f_grid.nelem(); iv++)
                propmat_clearsky[sp].SetAtPosition(
                    ext_mat_Nse[i_ss][i_se](iv, 0, 0, joker, joker), iv);
            else
              for (Index iv = 0; iv < f_grid.nelem(); iv++)
                propmat_clearsky[sp].SetAtPosition(
                    ext_mat_Nse[i_ss][i_se](0, 0, 0, joker, joker), iv);
          }
          propmat_clearsky[sp] *= rtp_vmr[sp];
        }

        // For temperature derivatives (so we don't need to check it in jac loop)
        if (do_jac_temperature) {
          if (t_ok(i_se_flat, 1) < 0.) {
            ostringstream os;
            os << "Temperature interpolation error (in perturbation):\n"
               << "scat species #" << i_ss << ", scat elem #" << i_se << "\n";
            throw runtime_error(os.str());
          }
        }

        // For number density derivatives
        if (do_jac) rtp_vmr_sum += rtp_vmr[sp];

        for (Index iq = 0; iq < jacobian_quantities_position.nelem(); iq++) {
          
          if (jacobian_quantities[jacobian_quantities_position[iq]] ==
              Jacobian::Atm::Temperature) {
            if (use_abs_as_ext) {
              tmp(joker, joker, 0) =
                  abs_vec_Nse[i_ss][i_se](joker, 1, 0, joker);
              tmp(joker, joker, 0) -=
                  abs_vec_Nse[i_ss][i_se](joker, 0, 0, joker);
            } else {
              tmp = ext_mat_Nse[i_ss][i_se](joker, 1, 0, joker, joker);
              tmp -= ext_mat_Nse[i_ss][i_se](joker, 0, 0, joker, joker);
            }

            tmp *= rtp_vmr[sp];
            tmp /= dT;

            if (nf > 1)
              for (Index iv = 0; iv < f_grid.nelem(); iv++)
                if (use_abs_as_ext)
                  dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(
                      tmp(iv, joker, 0), iv);
                else
                  dpropmat_clearsky_dx[iq].AddAtPosition(tmp(iv, joker, joker),
                                                         iv);
            else
              for (Index iv = 0; iv < f_grid.nelem(); iv++)
                if (use_abs_as_ext)
                  dpropmat_clearsky_dx[iq].AddAbsorptionVectorAtPosition(
                      tmp(0, joker, 0), iv);
                else
                  dpropmat_clearsky_dx[iq].AddAtPosition(tmp(0, joker, joker),
                                                         iv);
          }

          else if (jacobian_quantities[jacobian_quantities_position[iq]] ==
                   Jacobian::Atm::Particulates) {
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              dpropmat_clearsky_dx[iq].AddAtPosition(propmat_clearsky[sp], iv);
          }
        }
      }
      sp++;
      i_se_flat++;
    }
  }

  //checking that no further 'particle' entry left after all scat_data entries
  //are processes. this is basically not necessary. but checking it anyway to
  //really be safe. remove later, when more extensively tested.
  while (sp < na) {
    assert(abs_species[sp][0].Type() != SpeciesTag::TYPE_PARTICLES);
    sp++;
  }

  if (rtp_vmr_sum != 0.0) {
    for (Index iq = 0; iq < jacobian_quantities_position.nelem(); iq++) {
      if (jacobian_quantities[jacobian_quantities_position[iq]] ==
          Jacobian::Atm::Particulates) {
        dpropmat_clearsky_dx[iq] /= rtp_vmr_sum;
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddOnTheFly(  // Workspace reference:
    Workspace& ws,
    // WS Output:
    ArrayOfPropagationMatrix& propmat_clearsky,
    ArrayOfStokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_dx_source,
    ArrayOfStokesVector& nlte_dsource_dx,
    // WS Input:
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Agenda& abs_xsec_agenda,
    // Verbosity object:
    const Verbosity& verbosity) {

  // Output of AbsInputFromRteScalars:
  Vector abs_p;
  Vector abs_t;
  Matrix abs_vmrs;
  // Output of abs_h2oSet:
  Vector abs_h2o;
  // Output of abs_coefCalc:
  Matrix abs_coef, src_coef;
  ArrayOfMatrix abs_coef_per_species, src_coef_per_species, dabs_coef_dx,
      dsrc_coef_dx;
      
  AbsInputFromRteScalars(abs_p,
                         abs_t,
                         abs_vmrs,
                         rtp_pressure,
                         rtp_temperature,
                         rtp_vmr,
                         verbosity);

  // Absorption cross sections per tag group.
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfMatrix src_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx, dsrc_xsec_per_species_dx;

  // Make all species active:
  ArrayOfIndex abs_species_active(abs_species.nelem());
  for (Index i = 0; i < abs_species.nelem(); ++i) abs_species_active[i] = i;

  // Call agenda to calculate absorption:
  abs_xsec_agendaExecute(ws,
                         abs_xsec_per_species,
                         src_xsec_per_species,
                         dabs_xsec_per_species_dx,
                         dsrc_xsec_per_species_dx,
                         abs_species,
                         jacobian_quantities,
                         abs_species_active,
                         f_grid,
                         abs_p,
                         abs_t,
                         rtp_nlte,
                         abs_vmrs,
                         abs_xsec_agenda);
  
  // Calculate absorption coefficients from cross sections:
  abs_coefCalcFromXsec(abs_coef,
                       src_coef,
                       dabs_coef_dx,
                       dsrc_coef_dx,
                       abs_coef_per_species,
                       src_coef_per_species,
                       abs_xsec_per_species,
                       src_xsec_per_species,
                       dabs_xsec_per_species_dx,
                       dsrc_xsec_per_species_dx,
                       abs_species,
                       jacobian_quantities,
                       abs_vmrs,
                       abs_p,
                       abs_t,
                       verbosity);

  // Now add abs_coef_per_species to propmat_clearsky:
  propmat_clearskyAddFromAbsCoefPerSpecies(propmat_clearsky,
                                           dpropmat_clearsky_dx,
                                           abs_coef_per_species,
                                           dabs_coef_dx,
                                           verbosity);

  // Now turn nlte_source from absorption into a proper source function
  if (not nlte_source.empty())
    nlte_sourceFromTemperatureAndSrcCoefPerSpecies(nlte_source,
                                                   dnlte_dx_source,
                                                   nlte_dsource_dx,
                                                   src_coef_per_species,
                                                   dsrc_coef_dx,
                                                   jacobian_quantities,
                                                   f_grid,
                                                   rtp_temperature,
                                                   verbosity);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyZero(ArrayOfPropagationMatrix& propmat_clearsky,
                          const Vector& f_grid,
                          const Index& stokes_dim,
                          const Verbosity&) {
  propmat_clearsky.resize(1);
  propmat_clearsky[0] = PropagationMatrix(f_grid.nelem(), stokes_dim);
  propmat_clearsky[0].SetZero();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyForceNegativeToZero(
    ArrayOfPropagationMatrix& propmat_clearsky, const Verbosity&) {
  for (auto& pm : propmat_clearsky)
    for (Index i = 0; i < pm.NumberOfFrequencies(); i++)
      if (pm.Kjj()[i] < 0.0) pm.SetAtPosition(0.0, i);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromBuiltin(SpeciesAuxData& isotopologue_ratios,
                                        const Verbosity&) {
  fillSpeciesAuxDataWithIsotopologueRatiosFromSpeciesData(isotopologue_ratios);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void partition_functionsInitFromBuiltin(SpeciesAuxData& partition_functions,
                                        const Verbosity&) {
  fillSpeciesAuxDataWithPartitionFunctionsFromSpeciesData(partition_functions);
}

#ifdef ENABLE_NETCDF
/* Workspace method: Doxygen documentation will be auto-generated */
/* Included by Claudia Emde, 20100707 */
void WriteMolTau(  //WS Input
    const Vector& f_grid,
    const Tensor3& z_field,
    const Tensor7& propmat_clearsky_field,
    const Index& atmosphere_dim,
    //Keyword
    const String& filename,
    const Verbosity&) {
  int retval, ncid;
  int nlev_dimid, nlyr_dimid, nwvl_dimid, stokes_dimid, none_dimid;
  int dimids[4];
  int wvlmin_varid, wvlmax_varid, z_varid, wvl_varid, tau_varid;

  if (atmosphere_dim != 1)
    throw runtime_error("WriteMolTau can only be used for atmosphere_dim=1");

#pragma omp critical(netcdf__critical_region)
  {
    // Open file
    if ((retval = nc_create(filename.c_str(), NC_CLOBBER, &ncid)))
      nca_error(retval, "nc_create");

    // Define dimensions
    if ((retval = nc_def_dim(ncid, "nlev", (int)z_field.npages(), &nlev_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval =
             nc_def_dim(ncid, "nlyr", (int)z_field.npages() - 1, &nlyr_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid, "nwvl", (int)f_grid.nelem(), &nwvl_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid, "none", 1, &none_dimid)))
      nca_error(retval, "nc_def_dim");

    if ((retval = nc_def_dim(ncid,
                             "nstk",
                             (int)propmat_clearsky_field.nbooks(),
                             &stokes_dimid)))
      nca_error(retval, "nc_def_dim");

    // Define variables
    if ((retval = nc_def_var(
             ncid, "wvlmin", NC_DOUBLE, 1, &none_dimid, &wvlmin_varid)))
      nca_error(retval, "nc_def_var wvlmin");

    if ((retval = nc_def_var(
             ncid, "wvlmax", NC_DOUBLE, 1, &none_dimid, &wvlmax_varid)))
      nca_error(retval, "nc_def_var wvlmax");

    if ((retval = nc_def_var(ncid, "z", NC_DOUBLE, 1, &nlev_dimid, &z_varid)))
      nca_error(retval, "nc_def_var z");

    if ((retval =
             nc_def_var(ncid, "wvl", NC_DOUBLE, 1, &nwvl_dimid, &wvl_varid)))
      nca_error(retval, "nc_def_var wvl");

    dimids[0] = nlyr_dimid;
    dimids[1] = nwvl_dimid;
    dimids[2] = stokes_dimid;
    dimids[3] = stokes_dimid;

    if ((retval =
             nc_def_var(ncid, "tau", NC_DOUBLE, 4, &dimids[0], &tau_varid)))
      nca_error(retval, "nc_def_var tau");

    // Units
    if ((retval = nc_put_att_text(ncid, wvlmin_varid, "units", 2, "nm")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, wvlmax_varid, "units", 2, "nm")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, z_varid, "units", 2, "km")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, wvl_varid, "units", 2, "nm")))
      nca_error(retval, "nc_put_att_text");

    if ((retval = nc_put_att_text(ncid, tau_varid, "units", 1, "-")))
      nca_error(retval, "nc_put_att_text");

    // End define mode. This tells netCDF we are done defining
    // metadata.
    if ((retval = nc_enddef(ncid))) nca_error(retval, "nc_enddef");

    // Assign data
    double wvlmin[1];
    wvlmin[0] = SPEED_OF_LIGHT / f_grid[f_grid.nelem() - 1] * 1e9;
    if ((retval = nc_put_var_double(ncid, wvlmin_varid, &wvlmin[0])))
      nca_error(retval, "nc_put_var");

    double wvlmax[1];
    wvlmax[0] = SPEED_OF_LIGHT / f_grid[0] * 1e9;
    if ((retval = nc_put_var_double(ncid, wvlmax_varid, &wvlmax[0])))
      nca_error(retval, "nc_put_var");

    double z[z_field.npages()];
    for (int iz = 0; iz < z_field.npages(); iz++)
      z[iz] = z_field(z_field.npages() - 1 - iz, 0, 0) * 1e-3;

    if ((retval = nc_put_var_double(ncid, z_varid, &z[0])))
      nca_error(retval, "nc_put_var");

    double wvl[f_grid.nelem()];
    for (int iv = 0; iv < f_grid.nelem(); iv++)
      wvl[iv] = SPEED_OF_LIGHT / f_grid[f_grid.nelem() - 1 - iv] * 1e9;

    if ((retval = nc_put_var_double(ncid, wvl_varid, &wvl[0])))
      nca_error(retval, "nc_put_var");

    const Index zfnp = z_field.npages() - 1;
    const Index fgne = f_grid.nelem();
    const Index amfnb = propmat_clearsky_field.nbooks();

    Tensor4 tau(zfnp, fgne, amfnb, amfnb, 0.);

    // Calculate average tau for layers
    for (int is = 0; is < propmat_clearsky_field.nlibraries(); is++)
      for (int iz = 0; iz < zfnp; iz++)
        for (int iv = 0; iv < fgne; iv++)
          for (int is1 = 0; is1 < amfnb; is1++)
            for (int is2 = 0; is2 < amfnb; is2++)
              // sum up all species
              tau(iz, iv, is1, is2) +=
                  0.5 *
                  (propmat_clearsky_field(is,
                                          f_grid.nelem() - 1 - iv,
                                          is1,
                                          is2,
                                          z_field.npages() - 1 - iz,
                                          0,
                                          0) +
                   propmat_clearsky_field(is,
                                          f_grid.nelem() - 1 - iv,
                                          is1,
                                          is2,
                                          z_field.npages() - 2 - iz,
                                          0,
                                          0)) *
                  (z_field(z_field.npages() - 1 - iz, 0, 0) -
                   z_field(z_field.npages() - 2 - iz, 0, 0));

    if ((retval = nc_put_var_double(ncid, tau_varid, tau.get_c_array())))
      nca_error(retval, "nc_put_var");

    // Close the file
    if ((retval = nc_close(ncid))) nca_error(retval, "nc_close");
  }
}

#else

void WriteMolTau(  //WS Input
    const Vector& f_grid _U_,
    const Tensor3& z_field _U_,
    const Tensor7& propmat_clearsky_field _U_,
    const Index& atmosphere_dim _U_,
    //Keyword
    const String& filename _U_,
    const Verbosity&) {
  throw runtime_error(
      "The workspace method WriteMolTau is not available"
      "because ARTS was compiled without NetCDF support.");
}

#endif /* ENABLE_NETCDF */

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLines(
    // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfMatrix& src_xsec_per_species,
    ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    ArrayOfArrayOfMatrix& dsrc_xsec_per_species_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const EnergyLevelMap& abs_nlte,
    const Matrix& abs_vmrs,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesAuxData& isotopologue_ratios,
    const SpeciesAuxData& partition_functions,
    const Index& lbl_checked,
    const Verbosity&) {
  if (not abs_lines_per_species.nelem()) return;
  
  if (not lbl_checked)
    throw std::runtime_error("Please set lbl_checked true to use this function");

  // Check that all temperatures are above 0 K
  if (min(abs_t) < 0) {
    std::ostringstream os;
    os << "Temperature must be at least 0 K. But you request an absorption\n"
       << "calculation at " << min(abs_t) << " K!";
    throw std::runtime_error(os.str());
  }

  // Check that all parameters that should have the number of tag
  // groups as a dimension are consistent.
  {
    const Index n_tgs = abs_species.nelem();
    const Index n_xsec = abs_xsec_per_species.nelem();
    const Index n_vmrs = abs_vmrs.nrows();
    const Index n_lines = abs_lines_per_species.nelem();

    if (n_tgs not_eq n_xsec or n_tgs not_eq n_vmrs or n_tgs not_eq n_lines) {
      std::ostringstream os;
      os << "The following variables must all have the same dimension:\n"
         << "abs_species:           " << abs_species.nelem() << '\n'
         << "abs_xsec_per_species:  " << abs_xsec_per_species.nelem() << '\n'
         << "abs_vmrs:              " << abs_vmrs.nrows() << '\n'
         << "abs_lines_per_species: " << abs_lines_per_species.nelem() << '\n';
      throw std::runtime_error(os.str());
    }
  }

  // Meta variables that explain the calculations required
  const ArrayOfIndex jac_pos = equivalent_propmattype_indexes(jacobian_quantities);

  // Skipping uninteresting data
  static Matrix dummy1(0, 0);
  static ArrayOfMatrix dummy2(0);

  // Call xsec_species for each tag group.
  for (Index ii = 0; ii < abs_species_active.nelem(); ++ii) {
    const Index i = abs_species_active[ii];
    
    if (not abs_species[i].nelem() or is_zeeman(abs_species[i]))
      continue;
    
    for (auto& lines: abs_lines_per_species[i]) {
      xsec_species(
          abs_xsec_per_species[i],
          src_xsec_per_species[i],
          dummy1,
          dabs_xsec_per_species_dx[i],
          dsrc_xsec_per_species_dx[i],
          dummy2,
          jacobian_quantities,
          jac_pos,
          f_grid,
          abs_p,
          abs_t,
          abs_nlte,
          abs_vmrs,
          abs_species,
          lines,
          isotopologue_ratios.getIsotopologueRatio(lines.QuantumIdentity()),
          partition_functions.getParamType(lines.QuantumIdentity()),
          partition_functions.getParam(lines.QuantumIdentity()));
    }
  }  // End of species for loop.
}

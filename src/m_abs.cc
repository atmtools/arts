
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
#include "absorption.h"
#include "array.h"
#include "arts.h"
#include "auto_md.h"
#include "check_input.h"
#include "depr.h"
#include "file.h"
#include "global_data.h"
#include "hitran_species.h"
#include "jacobian.h"
#include "legacy_continua.h"
#include "lineshape.h"
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
#include <algorithm>
#include <cmath>

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
        // Test isotopologue, we have to hit the end of the list for no isotopologue or the exact value
        if (not same_or_joker(this_tag.Isotopologue(), lines.Isotopologue()))
          continue;
        
        // If there is a frequency range, we have to check so that only selected lines are included
        if (this_tag.lower_freq >= 0 or this_tag.upper_freq >= 0) {
          const Numeric low = (this_tag.lower_freq >= 0) ? this_tag.lower_freq : std::numeric_limits<Numeric>::lowest();
          const Numeric upp = (this_tag.upper_freq >= 0) ? this_tag.upper_freq : std::numeric_limits<Numeric>::max();
          
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

  // We want to make lists of included and excluded species:
  ArrayOfString included(0), excluded(0);

  tgs.resize(0);

  for (Index i = 0; i < Index(Species::Species::FINAL); ++i) {
    const String specname = Species::toShortName(Species::Species(i));

    String filename = basename;
    if (basename.length() && basename[basename.length() - 1] != '/')
      filename += ".";
    filename += specname;

    try {
      find_xml_file(filename, verbosity);
      // Add to included list:
      included.push_back(specname);
      
      // Add this tag group to tgs:
      tgs.emplace_back(ArrayOfSpeciesTag(specname));
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

  // We want to make lists of all species
  ArrayOfString specs(0);
  for (Index i = 0; i < Index(Species::Species::FINAL); ++i) {
    if (Species::Species(i) not_eq Species::Species::Bath) {
      specs.emplace_back(Species::toShortName(Species::Species(i)));
    }
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
  ARTS_USER_ERROR_IF (1 != atmosphere_dim,
      "Atmospheric dimension must be 1D, but atmosphere_dim is ",
      atmosphere_dim, ".")

  abs_p = p_grid;
  abs_t = t_field(joker, 0, 0);
  abs_vmrs = vmr_field(joker, joker, 0, 0);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_coefCalcFromXsec(  // WS Output:
    Matrix& abs_coef,
    ArrayOfMatrix& dabs_coef_dx,
    ArrayOfMatrix& abs_coef_per_species,
    // WS Input:
    const ArrayOfMatrix& abs_xsec_per_species,
    const ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Matrix& abs_vmrs,
    const Vector& abs_p,
    const Vector& abs_t,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  // Check that abs_vmrs and abs_xsec_per_species really have compatible
  // dimensions. In abs_vmrs there should be one row for each tg:
  ARTS_USER_ERROR_IF(abs_vmrs.nrows() != abs_xsec_per_species.nelem(),
    "Variable abs_vmrs must have compatible dimension to abs_xsec_per_species.\n"
    "abs_vmrs.nrows() = ", abs_vmrs.nrows(), "\n"
    "abs_xsec_per_species.nelem() = ", abs_xsec_per_species.nelem())

  // Check that number of altitudes are compatible. We only check the
  // first element, this is possilble because within arts all elements
  // are on the same altitude grid.
  ARTS_USER_ERROR_IF(abs_vmrs.ncols() != abs_xsec_per_species[0].ncols(),
    "Variable abs_vmrs must have same numbers of altitudes as abs_xsec_per_species.\n"
    "abs_vmrs.ncols() = ", abs_vmrs.ncols(), "\n"
    "abs_xsec_per_species[0].ncols() = ", abs_xsec_per_species[0].ncols());

  // Check dimensions of abs_p and abs_t:
  chk_size("abs_p", abs_p, abs_vmrs.ncols());
  chk_size("abs_t", abs_t, abs_vmrs.ncols());

  // Initialize abs_coef and abs_coef_per_species. The array dimension of abs_coef_per_species
  // is the same as that of abs_xsec_per_species. The dimension of abs_coef should
  // be equal to one of the abs_xsec_per_species enries.

  abs_coef.resize(abs_xsec_per_species[0].nrows(),
                  abs_xsec_per_species[0].ncols());
  abs_coef = 0;
  
  dabs_coef_dx.resize(jacobian_quantities.nelem());

  for (Index ii = 0; ii < jacobian_quantities.nelem(); ii++) {
    const auto& deriv = jacobian_quantities[ii];
    
    if (not deriv.propmattype()) continue;
    
    dabs_coef_dx[ii].resize(abs_xsec_per_species[0].nrows(),
                            abs_xsec_per_species[0].ncols());
    dabs_coef_dx[ii] = 0.0;
  }

  abs_coef_per_species.resize(abs_xsec_per_species.nelem());

  out3
      << "  Computing abs_coef and abs_coef_per_species from abs_xsec_per_species.\n";
  // Loop through all tag groups
  for (Index i = 0; i < abs_xsec_per_species.nelem(); ++i) {
    out3 << "  Tag group " << i << "\n";

    // Make this element of abs_xsec_per_species the right size:
    abs_coef_per_species[i].resize(abs_xsec_per_species[i].nrows(),
                                   abs_xsec_per_species[i].ncols());
    abs_coef_per_species[i] = 0;  // Initialize all elements to 0.

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

        for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
          const auto& deriv = jacobian_quantities[iq];
          
          if (not deriv.propmattype()) continue;
          
          if (deriv == Jacobian::Atm::Temperature) {
            dabs_coef_dx[iq](k, j) +=
                (dabs_xsec_per_species_dx[i][iq](k, j) * n +
                 abs_xsec_per_species[i](k, j) * dn_dT) *
                abs_vmrs(i, j);
          } else if (deriv == Jacobian::Line::VMR) {
            bool seco = false, main = false;
            for (const auto& s : abs_species[i]) {
              if (species_match(
                      deriv, s.cia_2nd_species) or
                  s.type not_eq Species::TagType::Cia)
                seco = true;
              if (species_iso_match(
                      deriv,
                      s.Isotopologue()))
                main = true;
            }
            if (main and seco) {
              dabs_coef_dx[iq](k, j) +=
                  (dabs_xsec_per_species_dx[i][iq](k, j) * abs_vmrs(i, j) +
                   abs_xsec_per_species[i](k, j)) *
                  n;
            } else if (main) {
              dabs_coef_dx[iq](k, j) += abs_xsec_per_species[i](k, j) * n;
            } else if (seco) {
              dabs_coef_dx[iq](k, j) +=
                  dabs_xsec_per_species_dx[i][iq](k, j) * abs_vmrs(i, j) * n;
            }
          } else if (deriv == Jacobian::Special::ArrayOfSpeciesTagVMR) {
            dabs_coef_dx[iq](k, j) += abs_xsec_per_species[i](k, j) * n;
          } else {
            dabs_coef_dx[iq](k, j) +=
                dabs_xsec_per_species_dx[i][iq](k, j) * n * abs_vmrs(i, j);
          }
        }
      }
    }

    // Add up to the total absorption:
    abs_coef += abs_coef_per_species[i];  // In Matpack you can use the +=
        // operator to do elementwise addition.
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesInit(  // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfArrayOfMatrix& dabs_xsec_per_species_dx,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& tgs,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Index& abs_xsec_agenda_checked,
    const Verbosity& verbosity) {
  CREATE_OUT3;

  ARTS_USER_ERROR_IF (!abs_xsec_agenda_checked,
        "You must call *abs_xsec_agenda_checkedCalc* before calling this method.");
  
  // Sizes
  const Index nq = jacobian_quantities.nelem();
  const Index nf = f_grid.nelem();
  const Index np = abs_p.nelem();
  const Index ns = tgs.nelem();

  // We need to check that abs_species_active doesn't have more elements than
  // abs_species (abs_xsec_agenda_checkedCalc doesn't know abs_species_active.
  // Usually we come here through an agenda call, where abs_species_active has
  // been properly created somewhere internally. But we might get here by
  // direct call, and then need to be safe!).
  ARTS_USER_ERROR_IF (ns < abs_species_active.nelem(),
    "abs_species_active (n=", abs_species_active.nelem(),
    ") not allowed to have more elements than abs_species (n=",
    ns, ")!\n")
  
  // Make elements the right size if they are not already the right size
  if (abs_xsec_per_species.nelem() not_eq ns) abs_xsec_per_species.resize(ns);
  if (dabs_xsec_per_species_dx.nelem() not_eq ns) dabs_xsec_per_species_dx.resize(ns);
  
  // Loop abs_xsec_per_species and make each matrix the right size,
  // initializing to zero.
  // But skip inactive species, loop only over the active ones.
  for (auto& i: abs_species_active) {
    ARTS_USER_ERROR_IF (i >= ns,
      "*abs_species_active* contains an invalid species index.\n"
      "Species index must be between 0 and ", ns - 1)
    
    // Make elements the right size if they are not already the right size, then reset them
    if (abs_xsec_per_species[i].nrows() == nf and abs_xsec_per_species[i].ncols() == np) {
      abs_xsec_per_species[i] = 0.0;
    } else {
      abs_xsec_per_species[i] = Matrix(nf, np, 0.0);
    }
    
    // Make elements the right size if they are not already the right size, then reset them
    if (dabs_xsec_per_species_dx[i].nelem() not_eq nq) {
      dabs_xsec_per_species_dx[i] = ArrayOfMatrix(nq, Matrix(nf, np, 0.0));
    } else {
      for (Index j=0; j<nq; j++) {
        if (dabs_xsec_per_species_dx[i][j].nrows() == nf and dabs_xsec_per_species_dx[i][j].ncols() == np) {
          dabs_xsec_per_species_dx[i][j] = 0.0;
        } else {
          dabs_xsec_per_species_dx[i][j] = Matrix(nf, np, 0.0);
        }
      }
    }
  }

  ostringstream os;
  os << "  Initialized abs_xsec_per_species.\n"
     << "  Number of frequencies        : " << nf << "\n"
     << "  Number of pressure levels    : " << np << "\n";
  out3 << os.str();
}


String continua_model_error_message(const ArrayOfString& abs_cont_names,
                                    const ArrayOfVector& abs_cont_parameters,
                                    const ArrayOfString& abs_cont_models) {
  std::ostringstream os;
  
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
  
  return os.str();
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
  ARTS_USER_ERROR_IF (tgs.nelem() != abs_xsec_per_species.nelem() || tgs.nelem() != abs_vmrs.nrows(),
    "The following variables must all have the same dimension:\n"
    "tgs:          ", tgs.nelem(), "\n"
    "abs_xsec_per_species:  ", abs_xsec_per_species.nelem(), "\n"
    "abs_vmrs.nrows():      ", abs_vmrs.nrows())

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
  ARTS_USER_ERROR_IF (abs_cont_names.nelem() != abs_cont_parameters.nelem(),
    continua_model_error_message(abs_cont_names, abs_cont_parameters, abs_cont_models))

  // Check that abs_p, abs_t, and abs_vmrs have the same
  // dimension. This could be a user error, so we throw a
  // runtime_error.

  ARTS_USER_ERROR_IF (abs_t.nelem() != abs_p.nelem(),
    "Variable abs_t must have the same dimension as abs_p.\n"
    "abs_t.nelem() = ", abs_t.nelem(), '\n',
    "abs_p.nelem() = ", abs_p.nelem())

  ARTS_USER_ERROR_IF (abs_vmrs.ncols() != abs_p.nelem(),
    "Variable dimension abs_vmrs.ncols() must\n"
    "be the same as abs_p.nelem().\n"
    "abs_vmrs.ncols() = ", abs_vmrs.ncols(), '\n',
    "abs_p.nelem() = ", abs_p.nelem())

  // We set abs_h2o, abs_n2, and abs_o2 later, because we only want to
  // do it if the parameters are really needed.

  out3 << "  Calculating continuum spectra.\n";

  // Loop tag groups:
  for (Index ii = 0; ii < abs_species_active.nelem(); ++ii) {
    const Index i = abs_species_active[ii];

    // Go through the tags in the current tag group to see if they
    // are continuum tags:
    for (Index s = 0; s < tgs[i].nelem(); ++s) {
      // Continuum tags in the sense that we talk about here
      // (including complete absorption models) are marked by a special type.
      if (tgs[i][s].Type() == Species::TagType::Predefined) {
        // We have identified a continuum tag!

        // Get only the continuum name. The full tag name is something like:
        // H2O-HITRAN96Self-*-*. We want only the `H2O-HITRAN96Self' part:
        const String name = tgs[i][s].Isotopologue().FullName();

        if (name == "O2-MPM2020") continue;
                                
        // Check, if we have parameters for this model. For
        // this, the model name must be listed in
        // abs_cont_names.
        const Index n =
            find(abs_cont_names.begin(), abs_cont_names.end(), name) -
            abs_cont_names.begin();

        // n==abs_cont_names.nelem() indicates that
        // the name was not found.
        ARTS_USER_ERROR_IF (n == abs_cont_names.nelem(),
          "Cannot find model ", name, " in abs_cont_names.")

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
              for (Index iq = 0; iq < jacobian_quantities.nelem();
                   iq++) {
                const auto& deriv = jacobian_quantities[iq];
                
                if (not deriv.propmattype()) continue;
                
                if (is_frequency_parameter(deriv))
                  dabs_xsec_per_species_dx[i][iq](iv, ip) +=
                      (jacs_df(iv, ip) - normal(iv, ip)) * (1. / df);
                else if (deriv == Jacobian::Atm::Temperature)
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
    StokesVector& nlte_source,
    ArrayOfStokesVector& dnlte_source_dx,
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

  ARTS_USER_ERROR_IF (not n_species,
    "Must have at least one species.")

  Index n_f = src_coef_per_species[0].nrows();  // # frequencies

  // # pressures must be 1:
  ARTS_USER_ERROR_IF (1 not_eq src_coef_per_species[0].ncols(),
    "Must have exactly one pressure.")

  // Check frequency dimension of propmat_clearsky
  ARTS_USER_ERROR_IF (nlte_source.NumberOfFrequencies() not_eq n_f,
                      "Frequency dimension of nlte_source does not\n"
    "match abs_coef_per_species.")
  
  const Vector B = planck(f_grid, rtp_temperature);

  StokesVector sv(n_f, nlte_source.StokesDimensions());
  for (Index si = 0; si < n_species; ++si) {
    sv.Kjj() = src_coef_per_species[si](joker, 0);
    sv *= B;
    nlte_source.Kjj() += sv.Kjj();
  }

  // Jacobian
  for (Index ii = 0; ii < jacobian_quantities.nelem(); ii++) {
    const auto& deriv = jacobian_quantities[ii];
    
    if (not deriv.propmattype()) continue;
    
    if (deriv == Jacobian::Atm::Temperature) {
      const Vector dB = dplanck_dt(f_grid, rtp_temperature);

      for (Index si = 0; si < n_species; ++si) {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        dnlte_source_dx[ii].Kjj() += sv.Kjj();
      }

      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_source_dx[ii].Kjj() += sv.Kjj();
    } else if (is_frequency_parameter(deriv)) {
      const Vector dB = dplanck_df(f_grid, rtp_temperature);

      for (Index si = 0; si < n_species; ++si) {
        sv.Kjj() = src_coef_per_species[si](joker, 0);
        sv *= dB;
        dnlte_source_dx[ii].Kjj() += sv.Kjj();
      }

      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_source_dx[ii].Kjj() += sv.Kjj();
    } else {
      sv.Kjj() = dsrc_coef_dx[ii](joker, 0);
      sv *= B;
      dnlte_source_dx[ii].Kjj() += sv.Kjj();
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFromAbsCoefPerSpecies(  // WS Output:
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const ArrayOfMatrix& abs_coef_per_species,
    const ArrayOfMatrix& dabs_coef_dx,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfSpeciesTag& abs_species) {
  // propmat_clearsky has format
  // [ abs_species, f_grid, stokes_dim, stokes_dim ].
  // abs_coef_per_species has format ArrayOfMatrix (over species),
  // where for each species the matrix has format [f_grid, abs_p].

  Index n_species = abs_coef_per_species.nelem();  // # species

  ARTS_USER_ERROR_IF (0 == n_species,
    "Must have at least one species.")

  Index n_f = abs_coef_per_species[0].nrows();  // # frequencies

  // # pressures must be 1:
  ARTS_USER_ERROR_IF (1 not_eq abs_coef_per_species[0].ncols(),
    "Must have exactly one pressure.")

  // Check frequency dimension of propmat_clearsky
  ARTS_USER_ERROR_IF (propmat_clearsky.NumberOfFrequencies() not_eq n_f,
    "Frequency dimension of propmat_clearsky does not\n"
    "match abs_coef_per_species.")

  // Loop species and stokes dimensions, and add to propmat_clearsky:
  for (Index si = 0; si < n_species; ++si)
    propmat_clearsky.Kjj() += abs_coef_per_species[si](joker, 0);

  for (Index iqn = 0; iqn < dabs_coef_dx.nelem(); iqn++) {
    bool found_special = false;
    for (Index isp=0; isp<abs_species.nelem(); isp++) {
      if (jacobian_quantities[iqn] == abs_species[isp]) {
        dpropmat_clearsky_dx[iqn].Kjj() += abs_coef_per_species[isp](joker, 0);
        found_special = true;
        break;
      }
    }
    
    if (not found_special) {
      if (dabs_coef_dx[iqn].nrows() == n_f) {
        if (dabs_coef_dx[iqn].ncols() == 1) {
          dpropmat_clearsky_dx[iqn].Kjj() += dabs_coef_dx[iqn](joker, 0);
        } else {
          ARTS_USER_ERROR("Must have exactly one pressure.");
        }
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyInit(  //WS Output
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    //WS Input
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Vector& f_grid,
    const Index& stokes_dim,
    const Index& propmat_clearsky_agenda_checked,
    const Verbosity&) {
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();
  
  ARTS_USER_ERROR_IF (!propmat_clearsky_agenda_checked,
        "You must call *propmat_clearsky_agenda_checkedCalc* before calling this method.")

  ARTS_USER_ERROR_IF (not nf, "No frequencies");

  ARTS_USER_ERROR_IF (stokes_dim < 1 or stokes_dim > 4, "stokes_dim not in [1, 2, 3, 4]");
  
  // Set size of propmat_clearsky or reset it's values
  if (propmat_clearsky.StokesDimensions() == stokes_dim and propmat_clearsky.NumberOfFrequencies() == nf) {
    propmat_clearsky.SetZero();
  } else {
    propmat_clearsky = PropagationMatrix(nf, stokes_dim);
  }
  
  // Set size of dpropmat_clearsky_dx or reset it's values
  if (dpropmat_clearsky_dx.nelem() not_eq nq) {
    dpropmat_clearsky_dx = ArrayOfPropagationMatrix(nq, PropagationMatrix(nf, stokes_dim));
  } else {
    for (auto& pm: dpropmat_clearsky_dx) {
      if (pm.StokesDimensions() == stokes_dim and pm.NumberOfFrequencies() == nf) {
        pm.SetZero();
      } else {
        pm = PropagationMatrix(nf, stokes_dim);
      }
    }
  }
  
  // Set size of nlte_source or reset it's values
  if (nlte_source.StokesDimensions() == stokes_dim and nlte_source.NumberOfFrequencies() == nf) {
    nlte_source.SetZero();
  } else {
    nlte_source = StokesVector(nf, stokes_dim);
  }
  
  // Set size of dnlte_source_dx or reset it's values
  if (dnlte_source_dx.nelem() not_eq nq) {
    dnlte_source_dx = ArrayOfStokesVector(nq, StokesVector(nf, stokes_dim));
  } else {
    for (auto& pm: dnlte_source_dx) {
      if (pm.StokesDimensions() == stokes_dim and pm.NumberOfFrequencies() == nf) {
        pm.SetZero();
      } else {
        pm = StokesVector(nf, stokes_dim);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddFaraday(
    PropagationMatrix& propmat_clearsky,
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

  ARTS_USER_ERROR_IF (stokes_dim < 3,
        "To include Faraday rotation, stokes_dim >= 3 is required.")
  ARTS_USER_ERROR_IF (atmosphere_dim == 1 && rtp_los.nelem() < 1,
    "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
    "(at least zenith angle component for atmosphere_dim==1),\n"
    "but it is not.\n")
  ARTS_USER_ERROR_IF (atmosphere_dim > 1 && rtp_los.nelem() < 2,
    "For applying propmat_clearskyAddFaraday, los needs to be specified\n"
    "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
    "but it is not.\n")

  const bool do_magn_jac = do_magnetic_jacobian(jacobian_quantities);
  const Numeric dmag = magnetic_field_perturbation(jacobian_quantities);

  Index ife = -1;
  for (Index sp = 0; sp < abs_species.nelem() && ife < 0; sp++) {
    if (abs_species[sp].FreeElectrons()) {
      ife = sp;
    }
  }

  ARTS_USER_ERROR_IF (ife < 0,
    "Free electrons not found in *abs_species* and "
    "Faraday rotation can not be calculated.");

  const Numeric ne = rtp_vmr[ife];

  if (ne != 0 && (rtp_mag[0] != 0 || rtp_mag[1] != 0 || rtp_mag[2] != 0)) {
    // Include remaining terms, beside /f^2
    const Numeric c1 =
        2 * FRconst *
        dotprod_with_los(
            rtp_los, rtp_mag[0], rtp_mag[1], rtp_mag[2], atmosphere_dim);

    Numeric dc1_u = 0.0, dc1_v = 0.0, dc1_w = 0.0;
    if (do_magn_jac) {
      dc1_u = (2 * FRconst *
                    dotprod_with_los(rtp_los,
                                    rtp_mag[0] + dmag,
                                    rtp_mag[1],
                                    rtp_mag[2],
                                    atmosphere_dim) -
                c1) /
              dmag;
      dc1_v = (2 * FRconst *
                    dotprod_with_los(rtp_los,
                                    rtp_mag[0],
                                    rtp_mag[1] + dmag,
                                    rtp_mag[2],
                                    atmosphere_dim) -
                c1) /
              dmag;
      dc1_w = (2 * FRconst *
                    dotprod_with_los(rtp_los,
                                    rtp_mag[0],
                                    rtp_mag[1],
                                    rtp_mag[2] + dmag,
                                    atmosphere_dim) -
                c1) /
              dmag;
    }

    for (Index iv = 0; iv < f_grid.nelem(); iv++) {
      const Numeric f2 = f_grid[iv] * f_grid[iv];
      const Numeric r = ne * c1 / f2;
      propmat_clearsky.AddFaraday(r, iv);

      // The Jacobian loop
      for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
        if (is_frequency_parameter(jacobian_quantities[iq]))
          dpropmat_clearsky_dx[iq].AddFaraday(-2.0 * ne * r / f_grid[iv], iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::MagneticU)
          dpropmat_clearsky_dx[iq].AddFaraday(ne * dc1_u / f2, iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::MagneticV)
          dpropmat_clearsky_dx[iq].AddFaraday(ne * dc1_v / f2, iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::MagneticW)
          dpropmat_clearsky_dx[iq].AddFaraday(ne * dc1_w / f2, iv);
        else if (jacobian_quantities[iq] == Jacobian::Atm::Electrons)
          dpropmat_clearsky_dx[iq].AddFaraday(r, iv);
        else if (jacobian_quantities[iq] == abs_species[ife])
          dpropmat_clearsky_dx[iq].AddFaraday(r, iv);
      }
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddParticles(
    // WS Output:
    PropagationMatrix& propmat_clearsky,
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
  ARTS_USER_ERROR_IF (scat_data_checked != 1,
        "The scat_data must be flagged to have "
        "passed a consistency check (scat_data_checked=1).")

  const Index ns = TotalNumberOfElements(scat_data);
  Index np = 0;
  for (Index sp = 0; sp < abs_species.nelem(); sp++) {
    if (abs_species[sp].Particles()) {
      np++;
    }
  }

  ARTS_USER_ERROR_IF (np == 0,
    "For applying propmat_clearskyAddParticles, *abs_species* needs to"
    "contain species 'particles', but it does not.\n")

  ARTS_USER_ERROR_IF (ns != np,
    "Number of 'particles' entries in abs_species and of elements in\n"
    "*scat_data* needs to be identical. But you have " , np,
    " 'particles' entries\n"
    "and ", ns, " *scat_data* elements.\n")

  ARTS_USER_ERROR_IF (atmosphere_dim == 1 && rtp_los.nelem() < 1,
    "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
    "(at least zenith angle component for atmosphere_dim==1),\n"
    "but it is not.\n")
  ARTS_USER_ERROR_IF (atmosphere_dim > 1 && rtp_los.nelem() < 2,
    "For applying *propmat_clearskyAddParticles*, *rtp_los* needs to be specified\n"
    "(both zenith and azimuth angle components for atmosphere_dim>1),\n"
    "but it is not.\n")

  // Use for rescaling vmr of particulates
  Numeric rtp_vmr_sum = 0.0;

  // Tests and setup partial derivatives
  const bool do_jac_temperature = do_temperature_jacobian(jacobian_quantities);
  const bool do_jac_frequencies = do_frequency_jacobian(jacobian_quantities);
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

  // Internal computations necessary since it relies on zero start
  PropagationMatrix internal_propmat(propmat_clearsky.NumberOfFrequencies(), propmat_clearsky.StokesDimensions());
  
  // loop over the scat_data and link them with correct vmr_field entry according
  // to the position of the particle type entries in abs_species.
  Index sp = 0;
  Index i_se_flat = 0;
  for (Index i_ss = 0; i_ss < scat_data.nelem(); i_ss++) {
    for (Index i_se = 0; i_se < scat_data[i_ss].nelem(); i_se++) {
      // forward to next particle entry in abs_species
      while (sp < na && not abs_species[sp].Particles())
        sp++;
      internal_propmat.SetZero();

      // running beyond number of abs_species entries when looking for next
      // particle entry. shouldn't happen, though.
      ARTS_ASSERT(sp < na);
      ARTS_USER_ERROR_IF (rtp_vmr[sp] < 0.,
        "Negative absorbing particle 'vmr' (aka number density)"
        " encountered:\n"
        "scat species #", i_ss, ", scat elem #", i_se,
        " (vmr_field entry #", sp, ")\n")
      
      if (rtp_vmr[sp] > 0.) {
        ARTS_USER_ERROR_IF (t_ok(i_se_flat, 0) < 0.,
          "Temperature interpolation error:\n"
          "scat species #", i_ss, ", scat elem #", i_se, "\n")
        if (use_abs_as_ext) {
          if (nf > 1)
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.AddAbsorptionVectorAtPosition(
                  abs_vec_Nse[i_ss][i_se](iv, 0, 0, joker), iv);
          else
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.AddAbsorptionVectorAtPosition(
                  abs_vec_Nse[i_ss][i_se](0, 0, 0, joker), iv);
        } else {
          if (nf > 1)
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.SetAtPosition(
                  ext_mat_Nse[i_ss][i_se](iv, 0, 0, joker, joker), iv);
          else
            for (Index iv = 0; iv < f_grid.nelem(); iv++)
              internal_propmat.SetAtPosition(
                  ext_mat_Nse[i_ss][i_se](0, 0, 0, joker, joker), iv);
        }
        propmat_clearsky += rtp_vmr[sp] * internal_propmat;
      }

      // For temperature derivatives (so we don't need to check it in jac loop)
      if (do_jac_temperature) {
        ARTS_USER_ERROR_IF (t_ok(i_se_flat, 1) < 0.,
            "Temperature interpolation error (in perturbation):\n"
            "scat species #", i_ss, ", scat elem #", i_se, "\n")
      }

      // For number density derivatives
      if (jacobian_quantities.nelem()) rtp_vmr_sum += rtp_vmr[sp];

      for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
        const auto& deriv = jacobian_quantities[iq];
        
        if (not deriv.propmattype()) continue;
        
        if (deriv == Jacobian::Atm::Temperature) {
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

        else if (deriv == Jacobian::Atm::Particulates) {
          for (Index iv = 0; iv < f_grid.nelem(); iv++)
            dpropmat_clearsky_dx[iq].AddAtPosition(internal_propmat, iv);
        }
        
        else if (deriv == abs_species[sp]) {
          dpropmat_clearsky_dx[iq] += internal_propmat;
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
    ARTS_ASSERT(abs_species[sp][0].Type() != Species::TagType::Particles);
    sp++;
  }

  if (rtp_vmr_sum != 0.0) {
    for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {
      const auto& deriv = jacobian_quantities[iq];
      
      if (not deriv.propmattype()) continue;
      
      if (deriv == Jacobian::Atm::Particulates) {
        dpropmat_clearsky_dx[iq] /= rtp_vmr_sum;
      }
    }
  }
}


void sparse_f_gridFromFrequencyGrid(Vector& sparse_f_grid,
                                    const Vector& f_grid,
                                    const Numeric& sparse_df,
                                    const String& speedup_option,
                                    // Verbosity object:
                                    const Verbosity&)
{
  // Return empty for nothing
  if (not f_grid.nelem()) {
    sparse_f_grid.resize(0);
    return;
  };
  
  switch (Options::toLblSpeedupOrThrow(speedup_option)) {
    case Options::LblSpeedup::LinearIndependent:
      sparse_f_grid = LineShape::linear_sparse_f_grid(f_grid, sparse_df);
      ARTS_ASSERT(LineShape::good_linear_sparse_f_grid(f_grid, sparse_f_grid))
      break;
    case Options::LblSpeedup::QuadraticIndependent:
      sparse_f_grid = LineShape::triple_sparse_f_grid(f_grid, sparse_df);
      break;
    case Options::LblSpeedup::None:
      sparse_f_grid.resize(0);
      break;
    case Options::LblSpeedup::FINAL: { /* Leave last */ }
  }
}

Vector create_sparse_f_grid_internal(const Vector& f_grid,
                                     const Numeric& sparse_df,
                                     const String& speedup_option,
                                     // Verbosity object:
                                     const Verbosity& verbosity)
{
  Vector sparse_f_grid;
  sparse_f_gridFromFrequencyGrid(sparse_f_grid, f_grid, sparse_df, speedup_option, verbosity);
  return sparse_f_grid;
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddLines(  // Workspace reference:
    // WS Output:
    PropagationMatrix& propmat_clearsky,
    StokesVector& nlte_source,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    ArrayOfStokesVector& dnlte_source_dx,
    // WS Input:
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
    const EnergyLevelMap& rtp_nlte,
    const Vector& rtp_vmr,
    const Index& nlte_do,
    const Index& lbl_checked,
    // WS User Generic inputs
    const Numeric& sparse_df,
    const Numeric& sparse_lim,
    const String& speedup_option,
    const ArrayOfSpeciesTag& select_speciestags,
    // Verbosity object:
    const Verbosity& verbosity) {
  
  // Size of problem
  const Index nf = f_grid.nelem();
  const Index nq = jacobian_quantities.nelem();
  const Index ns = abs_species.nelem();
  
  // Possible things that can go wrong in this code (excluding line parameters)
  ARTS_USER_ERROR_IF(not lbl_checked, "Must check LBL calculations")
  check_abs_species(abs_species);
  ARTS_USER_ERROR_IF(rtp_vmr.nelem() not_eq abs_species.nelem(),
                     "*rtp_vmr* must match *abs_species*")
  ARTS_USER_ERROR_IF(propmat_clearsky.NumberOfFrequencies() not_eq nf,
                     "*f_grid* must match *propmat_clearsky*")
  ARTS_USER_ERROR_IF(nlte_source.NumberOfFrequencies() not_eq nf,
                     "*f_grid* must match *nlte_source*")
  ARTS_USER_ERROR_IF(not nq and (nq not_eq dpropmat_clearsky_dx.nelem()),
                     "*dpropmat_clearsky_dx* must match derived form of *jacobian_quantities*")
  ARTS_USER_ERROR_IF(not nq and bad_propmat(dpropmat_clearsky_dx, f_grid),
                     "*dpropmat_clearsky_dx* must have frequency dim same as *f_grid*")
  ARTS_USER_ERROR_IF(nlte_do and (nq not_eq dnlte_source_dx.nelem()),
                     "*dnlte_source_dx* must match derived form of *jacobian_quantities* when non-LTE is on")
  ARTS_USER_ERROR_IF(nlte_do and bad_propmat(dnlte_source_dx, f_grid),
                     "*dnlte_source_dx* must have frequency dim same as *f_grid* when non-LTE is on")
  ARTS_USER_ERROR_IF(any_negative(f_grid), "Negative frequency (at least one value).")
  ARTS_USER_ERROR_IF(not is_increasing(f_grid), "Must be sorted and increasing.")
  ARTS_USER_ERROR_IF(any_negative(rtp_vmr), "Negative VMR (at least one value).")
  ARTS_USER_ERROR_IF(any_negative(rtp_nlte.Data()), "Negative NLTE (at least one value).")
  ARTS_USER_ERROR_IF(rtp_temperature <= 0, "Non-positive temperature")
  ARTS_USER_ERROR_IF(rtp_pressure <= 0, "Non-positive pressure")
  ARTS_USER_ERROR_IF(sparse_lim > 0 and sparse_df > sparse_lim, 
                    "If sparse grids are to be used, the limit must be larger than the grid-spacing.\n"
                    "The limit is ", sparse_lim, " Hz and the grid_spacing is ", sparse_df, " Hz")
  
  if (not nf) return;
  
  // Deal with sparse computational grid
  const Vector f_grid_sparse = create_sparse_f_grid_internal(f_grid, sparse_df, speedup_option, verbosity);
  const Options::LblSpeedup speedup_type = f_grid_sparse.nelem() ? Options::toLblSpeedupOrThrow(speedup_option) : Options::LblSpeedup::None;
  ARTS_USER_ERROR_IF(sparse_lim <= 0 and speedup_type not_eq Options::LblSpeedup::None,
                     "Must have a sparse limit if you set speedup_option")
  
  // Calculations data
  LineShape::ComputeData com(f_grid, jacobian_quantities, nlte_do);
  LineShape::ComputeData sparse_com(f_grid_sparse, jacobian_quantities, nlte_do);
  
  for (Index ispecies = 0; ispecies < ns; ispecies++) {
    if (select_speciestags.nelem() and select_speciestags not_eq abs_species[ispecies]) continue;
    
    // Skip it if there are no species or there is Zeeman requested
    if (not abs_species[ispecies].nelem() or abs_species[ispecies].Zeeman() or not abs_lines_per_species[ispecies].nelem())
      continue;
    
    for (auto& band : abs_lines_per_species[ispecies]) {
      LineShape::compute(com, sparse_com, band, jacobian_quantities, rtp_nlte, band.BroadeningSpeciesVMR(rtp_vmr, abs_species), abs_species[ispecies], rtp_vmr[ispecies],
                          isotopologue_ratios[band.Isotopologue()], rtp_pressure, rtp_temperature, 0, sparse_lim,
                          false, Zeeman::Polarization::Pi, speedup_type);
      
    }
  }
  
  switch (speedup_type) {
    case Options::LblSpeedup::LinearIndependent: com.interp_add_even(sparse_com); break;
    case Options::LblSpeedup::QuadraticIndependent: com.interp_add_triplequad(sparse_com); break;
    case Options::LblSpeedup::None: /* Do nothing */ break;
    case Options::LblSpeedup::FINAL: { /* Leave last */ }
  }
    
  // Sum up the propagation matrix
  propmat_clearsky.Kjj() += com.F.real();
  
  // Sum up the Jacobian
  for (Index j=0; j<nq; j++) {
    if (not jacobian_quantities[j].propmattype()) continue;
    dpropmat_clearsky_dx[j].Kjj() += com.dF.real()(joker, j);
  }
  
  if (nlte_do) {
    // Sum up the source vector
    nlte_source.Kjj() += com.N.real();
    
    // Sum up the Jacobian
    for (Index j=0; j<nq; j++) {
      if (not jacobian_quantities[j].propmattype()) continue;
      dnlte_source_dx[j].Kjj() += com.dN.real()(joker, j);
    }
  }
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyAddXsecAgenda(  // Workspace reference:
    Workspace& ws,
    // WS Output:
    PropagationMatrix& propmat_clearsky,
    ArrayOfPropagationMatrix& dpropmat_clearsky_dx,
    // WS Input:
    const Vector& f_grid,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const Numeric& rtp_pressure,
    const Numeric& rtp_temperature,
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
  Matrix abs_coef;
  ArrayOfMatrix abs_coef_per_species, dabs_coef_dx;
      
  AbsInputFromRteScalars(abs_p,
                         abs_t,
                         abs_vmrs,
                         rtp_pressure,
                         rtp_temperature,
                         rtp_vmr,
                         verbosity);

  // Absorption cross sections per tag group.
  ArrayOfMatrix abs_xsec_per_species;
  ArrayOfArrayOfMatrix dabs_xsec_per_species_dx;

  // Make all species active:
  ArrayOfIndex abs_species_active(abs_species.nelem());
  for (Index i = 0; i < abs_species.nelem(); ++i) abs_species_active[i] = i;

  // Call agenda to calculate absorption:
  abs_xsec_agendaExecute(ws,
                         abs_xsec_per_species,
                         dabs_xsec_per_species_dx,
                         abs_species,
                         jacobian_quantities,
                         abs_species_active,
                         f_grid,
                         abs_p,
                         abs_t,
                         abs_vmrs,
                         abs_xsec_agenda);
  
  // Calculate absorption coefficients from cross sections:
  abs_coefCalcFromXsec(abs_coef,
                       dabs_coef_dx,
                       abs_coef_per_species,
                       abs_xsec_per_species,
                       dabs_xsec_per_species_dx,
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
                                           jacobian_quantities, abs_species);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyZero(PropagationMatrix& propmat_clearsky,
                          const Vector& f_grid,
                          const Index& stokes_dim,
                          const Verbosity&) {
  propmat_clearsky = PropagationMatrix(f_grid.nelem(), stokes_dim);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void propmat_clearskyForceNegativeToZero(
    PropagationMatrix& propmat_clearsky, const Verbosity&) {
  for (Index i = 0; i < propmat_clearsky.NumberOfFrequencies(); i++)
    if (propmat_clearsky.Kjj()[i] < 0.0) propmat_clearsky.SetAtPosition(0.0, i);
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromBuiltin(SpeciesIsotopologueRatios& isotopologue_ratios,
                                       const Verbosity&) {
  isotopologue_ratios = Species::isotopologue_ratiosInitFromBuiltin();
}

/* Workspace method: Doxygen documentation will be auto-generated */
void isotopologue_ratiosInitFromHitran(SpeciesIsotopologueRatios& isotopologue_ratios,
                                       const String& option,
                                       const Verbosity&) {
  isotopologue_ratios = Hitran::isotopologue_ratios(Hitran::toTypeOrThrow(option));
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

  ARTS_USER_ERROR_IF (atmosphere_dim != 1,
    "WriteMolTau can only be used for atmosphere_dim=1")

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
  ARTS_USER_ERROR_IF(true,
      "The workspace method WriteMolTau is not available"
      "because ARTS was compiled without NetCDF support.");
}

#endif /* ENABLE_NETCDF */

/* Workspace method: Doxygen documentation will be auto-generated */
void abs_xsec_per_speciesAddLines(
    // WS Output:
    ArrayOfMatrix& abs_xsec_per_species,
    ArrayOfArrayOfMatrix&,
    // WS Input:
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfRetrievalQuantity& jacobian_quantities,
    const ArrayOfIndex& abs_species_active,
    const Vector& f_grid,
    const Vector& abs_p,
    const Vector& abs_t,
    const Matrix& abs_vmrs,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
    const SpeciesIsotopologueRatios& isotopologue_ratios,
    const Index& lbl_checked,
    const Verbosity&) {
  DEPRECATED_FUNCTION("abs_xsec_per_speciesAddLines", "2021-07-13",
             "This function is no longer up to date.  It only exists to satisfy "
             "lookup table calculations before these are updated.\n"
             "Once the lookup table calculations are up-to-date, this function "
             "is fully replaced with propmat_clearskyAddLines, with better functionality\n" )
  
  ARTS_USER_ERROR_IF(not lbl_checked, "Must check LBL calculations")
  ARTS_USER_ERROR_IF(jacobian_quantities.nelem(), "There's a hard deprecation of derivatives using old style lbl-calculations with derivatives, switch to propmat_clearskyAddLines")
  ARTS_USER_ERROR_IF (abs_species.nelem() not_eq abs_xsec_per_species.nelem() or
                      abs_species.nelem() not_eq abs_vmrs.nrows() or
                      abs_species.nelem() not_eq abs_lines_per_species.nelem(),
    "The following variables must all have the same dimension:\n"
    "abs_species:           ", abs_species.nelem(), '\n',
    "abs_xsec_per_species:  ", abs_xsec_per_species.nelem(), '\n',
    "abs_vmrs:              ", abs_vmrs.nrows(), '\n',
    "abs_lines_per_species: ", abs_lines_per_species.nelem(), '\n')
  ARTS_USER_ERROR_IF (min(abs_t) < 0,
    "Temperature must be at least 0 K. But you request an absorption\n"
    "calculation at ", min(abs_t), " K!")
  
  // Size of problem
  const Index np = abs_p.nelem();
  
  // Calculations data
  LineShape::ComputeData com(f_grid, jacobian_quantities, false);
  LineShape::ComputeData sparse_com(Vector(0), jacobian_quantities, false);
  constexpr Options::LblSpeedup speedup_type = Options::LblSpeedup::None;
  const EnergyLevelMap rtp_nlte;
  
  for (Index ip=0; ip<np; ip++) {
    for (Index ispecies: abs_species_active) {
      // Skip it if there are no species or there is Zeeman requested
      if (not abs_species[ispecies].nelem() or abs_species[ispecies].Zeeman() or not abs_lines_per_species[ispecies].nelem())
        continue;
      
      // Reset for legacy VMR jacobian
      com.reset();
      sparse_com.reset();
      
      for (auto& band : abs_lines_per_species[ispecies]) {
        LineShape::compute(com, sparse_com, band, jacobian_quantities, rtp_nlte, band.BroadeningSpeciesVMR(abs_vmrs(joker, ip), abs_species), abs_species[ispecies], abs_vmrs(ispecies, ip),
                           isotopologue_ratios[band.Isotopologue()], abs_p[ip], abs_t[ip], 0, 0, false, Zeeman::Polarization::Pi, speedup_type);
        
      }
      
      // Sum up the propagation matrix
      com.F /= number_density(abs_p[ip], abs_t[ip]) * abs_vmrs(ispecies, ip);
      abs_xsec_per_species[ispecies](joker, ip) += com.F.real();
    }
  }
}

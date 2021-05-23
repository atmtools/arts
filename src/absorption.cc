/* Copyright (C) 2000-2012
   Stefan Buehler  <sbuehler@ltu.se>
   Axel von Engeln <engeln@uni-bremen.de>

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

/**
  \file   absorption.cc

  Physical absorption routines. 

  The absorption workspace methods are
  in file m_abs.cc

  This is the file from arts-1-0, back-ported to arts-1-1.

  \author Stefan Buehler and Axel von Engeln
*/

#include "absorption.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <map>
#include "arts.h"
#include "auto_md.h"
#include "file.h"
#include "linescaling.h"
#include "lineshape.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"

#include "global_data.h"
#include "linefunctions.h"


void checkPartitionFunctions(const ArrayOfArrayOfSpeciesTag& abs_species) {
  for (auto& specs: abs_species) {
    for (auto& spec: specs) {
      if (spec.type == Species::TagType::Plain or spec.type == Species::TagType::Zeeman) {
        if (spec.Isotopologue().isotname not_eq Species::Joker) {
          ARTS_USER_ERROR_IF (not PartitionFunctions::has_partfun(spec.Isotopologue()),
            "Species: ", spec.Isotopologue().FullName(), " has no partition function\n",
            "You must recompile ARTS partition functions with data for this species to continue your calculations,\n"
            "or exclude the species from your computation setup")
        } else {
          auto all_isots = Species::isotopologues(spec.Isotopologue().spec);
          for (auto& isot: all_isots) {
            if (not Species::is_predefined_model(isot)) {
              ARTS_USER_ERROR_IF (not PartitionFunctions::has_partfun(isot),
                                  "Species: ", isot.FullName(), " has no partition function\n",
                                  "You must recompiler ARTS partition functions with data for this species to continue your calculations,\n"
                                  "or exclude the species from your computation setup.  Note that it is part a joker-species, so the\n"
                                  "troubling isotopologue is not directly defined in its name")
            }
          }
        }
      }
    }
  }
}

void checkIsotopologueRatios(const ArrayOfArrayOfSpeciesTag& abs_species,
                             const Species::IsotopologueRatios& isoratios) {
  // For the selected species, we check all isotopes by looping over the
  // species data. (Trying to check only the isotopes actually used gets
  // quite complicated, actually, so we do the simple thing here.)

  // Loop over the absorption species:
  for (auto& specs: abs_species) {
    for (auto& spec: specs) {
      const Index i = Species::find_species_index(spec.Isotopologue());
      ARTS_USER_ERROR_IF(i < 0 or i >= Species::IsotopologueRatios::maxsize, "Cannot find species: ", spec)
      ARTS_USER_ERROR_IF(std::isnan(isoratios[i]), "Chosen species: ", spec, " has no isotopologue ratio")
    }
  }
}

/** A little helper function to convert energy from units of
    wavenumber (cm^-1) to Joule (J). 

    This is used when reading HITRAN or JPL catalogue files, which
    have the lower state energy in cm^-1.

    \return Energy in J.
    \param[in]  e Energy in cm^-1.

    \author Stefan Buehler
    \date   2001-06-26 */
Numeric wavenumber_to_joule(Numeric e) {
  // Planck constant [Js]
  extern const Numeric PLANCK_CONST;

  // Speed of light [m/s]
  extern const Numeric SPEED_OF_LIGHT;

  // Constant to convert lower state energy from cm^-1 to J
  const Numeric lower_energy_const = PLANCK_CONST * SPEED_OF_LIGHT * 1E2;

  return e * lower_energy_const;
}

//!  set_abs_from_first_species.
/*!
 Returns vmr for the profile of the first tag group containing
 the given species.

 \author Oliver Lemke

 \param[out] vmr          Volume mixing ratio
 \param[in]  species_name Species Name
 \param[in]  abs_species  WS Input
 \param[in]  abs_vmrs     WS Input
 */
void set_vmr_from_first_species(Vector& vmr,
                                const String& species_name,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Matrix& abs_vmrs) {
  const Index index = find_first_species(abs_species, Species::fromShortName(species_name));

  vmr.resize(abs_vmrs.ncols());
  if (index < 0)
    vmr = -99;
  else
    vmr = abs_vmrs(index, Range(joker));
}

void xsec_species(Matrix& xsec,
                  Matrix& source,
                  Matrix& phase,
                  ArrayOfMatrix& dxsec_dx,
                  ArrayOfMatrix& dsource_dx,
                  ArrayOfMatrix& dphase_dx,
                  const ArrayOfRetrievalQuantity& jacobian_quantities,
                  const Vector& f_grid,
                  const Vector& abs_p,
                  const Vector& abs_t,
                  const EnergyLevelMap& abs_nlte,
                  const Matrix& abs_vmrs,
                  const ArrayOfArrayOfSpeciesTag& abs_species,
                  const AbsorptionLines& band,
                  const Numeric& isot_ratio) {
  // Size of problem
  const Index np = abs_p.nelem();      // number of pressure levels
  const Index nf = f_grid.nelem();     // number of Dirac frequencies
  const Index nl = band.NumLines();  // number of lines in the catalog
  const Index nj = jacobian_quantities.nelem();  // number of partial derivatives
  const Index nt = source.nrows();         // number of energy levels in NLTE

  // Type of problem
  const bool do_nonlte = nt;

  Linefunctions::InternalData scratch(nf, nj);
  Linefunctions::InternalData sum(nf, nj);
  
  // Test if the size of the problem is 0
  if (not np or not nf or not nl) return;
  
  // Constant for all lines
  const Numeric QT0 = single_partition_function(band.T0(), band.Isotopologue());

  ArrayOfString fail_msg;
  bool do_abort = false;

#pragma omp parallel for if (!arts_omp_in_parallel() && np > 1) \
    firstprivate(scratch, sum)
  for (Index ip = 0; ip < np; ip++) {
    if (do_abort) continue;
    try {
      // Constants for this level
      const Numeric& temperature = abs_t[ip];
      const Numeric& pressure = abs_p[ip];

      // Constants for this level
      const Numeric QT =
          single_partition_function(temperature, band.Isotopologue());
      const Numeric dQTdT = dsingle_partition_function_dT(
          temperature, band.Isotopologue());
      const Numeric DC =
          Linefunctions::DopplerConstant(temperature, band.SpeciesMass());
      const Numeric dDCdT = Linefunctions::dDopplerConstant_dT(temperature, DC);
      const Vector line_shape_vmr =
          band.BroadeningSpeciesVMR(abs_vmrs(joker, ip), abs_species);

      Linefunctions::set_cross_section_of_band(scratch,
                                               sum,
                                               f_grid,
                                               band,
                                               jacobian_quantities,
                                               line_shape_vmr,
                                               abs_nlte[ip],
                                               pressure,
                                               temperature,
                                               isot_ratio,
                                               0,
                                               DC,
                                               dDCdT,
                                               QT,
                                               dQTdT,
                                               QT0,
                                               false);

      // absorption cross-section
      MapToEigen(xsec).col(ip).noalias() += sum.F.real();
      for (Index j = 0; j < nj; j++) {
        if (not jacobian_quantities[j].propmattype()) continue;
        MapToEigen(dxsec_dx[j]).col(ip).noalias() += sum.dF.col(j).real();
      }

      // phase cross-section
      if (not phase.empty()) {
        MapToEigen(phase).col(ip).noalias() += sum.F.imag();
        for (Index j = 0; j < nj; j++) {
          if (not jacobian_quantities[j].propmattype()) continue;
          MapToEigen(dphase_dx[j]).col(ip).noalias() += sum.dF.col(j).imag();
        }
      }

      // source ratio cross-section
      if (do_nonlte) {
        MapToEigen(source).col(ip).noalias() += sum.N.real();
        for (Index j = 0; j < nj; j++) {
          if (not jacobian_quantities[j].propmattype()) continue;
          MapToEigen(dsource_dx[j]).col(ip).noalias() += sum.dN.col(j).real();
        }
      }
    } catch (const std::runtime_error& e) {
      ostringstream os;
      os << "Runtime-error in cross-section calculation at p_abs index " << ip
         << ": \n";
      os << e.what();
#pragma omp critical(xsec_species_cross_sections)
      {
        do_abort = true;
        fail_msg.push_back(os.str());
      }
    }
  }

  ARTS_USER_ERROR_IF (do_abort,
    "Error messages from failed cases:\n", fail_msg)
}

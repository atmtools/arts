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
#include "arts.h"
#include "arts_conversions.h"
#include "auto_md.h"
#include "file.h"
#include "linescaling.h"
#include "lineshape.h"
#include "logic.h"
#include "math_funcs.h"
#include "messages.h"
#include "partfun.h"
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <map>

#include "global_data.h"

void checkPartitionFunctions(const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species) {
  for (auto& abs_lines: abs_lines_per_species) {
    for (auto& band: abs_lines) {
      ARTS_USER_ERROR_IF (not PartitionFunctions::has_partfun(band.Isotopologue()),
                          "Species: ", band.Isotopologue().FullName(), " has no partition function\n",
                          "You must recompile ARTS partition functions with data for this species to continue your calculations,\n"
                          "or exclude the species from your computation setup")
    }
  }
}

void checkIsotopologueRatios(const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                             const Species::IsotopologueRatios& isoratios) {
  
  for (auto& abs_lines: abs_lines_per_species) {
    for (auto& band: abs_lines) {
      ARTS_USER_ERROR_IF (std::isnan(isoratios[band.Isotopologue()]),
                          "Species: ", band.Isotopologue().FullName(), " has no isotopologue ratios\n",
                          "You must add its isotopologue ratios to your included data or\n"
                          "exclude the species from your computation setup")
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
  return Conversion::kaycm2joule(e);
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

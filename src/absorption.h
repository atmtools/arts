/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@ltu.se>
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

/** \file
    Declarations required for the calculation of absorption coefficients.

    This is the file from arts-1-0, back-ported to arts-1-1.

    \author Stefan Buehler, Axel von Engeln
*/

#ifndef absorption_h
#define absorption_h

#include "species_tags.h"
#include "absorptionlines.h"
#include "array.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "matpack_data.h"
#include "mystring.h"
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>

/** Check that ARTS was compiled for all requested species tags */
void checkPartitionFunctions(const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species);

/** Check that isotopologue ratios for the given species are correctly defined. */
void checkIsotopologueRatios(const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species,
                             const Species::IsotopologueRatios& isoratios);

// A helper function for energy conversion:
Numeric wavenumber_to_joule(Numeric e);

//======================================================================
//             Functions related to species
//======================================================================

void set_vmr_from_first_species(Vector& vmr,
                                const String& species_name,
                                const ArrayOfArrayOfSpeciesTag& abs_species,
                                const Matrix& abs_vmrs);

#endif  // absorption_h

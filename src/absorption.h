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
#include "energylevelmap.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "matpack_data.h"
#include "messages.h"
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

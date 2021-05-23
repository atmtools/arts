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

#include <array>
#include <map>
#include <cmath>
#include <stdexcept>
#include <limits>
#include "species_tags.h"
#include "array.h"
#include "energylevelmap.h"
#include "gridded_fields.h"
#include "jacobian.h"
#include "matpackI.h"
#include "messages.h"
#include "mystring.h"
#include "absorptionlines.h"

/** Check that ARTS was compiled for all requested species tags */
void checkPartitionFunctions(const ArrayOfArrayOfSpeciesTag& abs_species);

namespace Species {
struct IsotopologueRatios {
  static constexpr Index maxsize = Index(Isotopologues.size());
  std::array<Numeric, maxsize> data;
  
  constexpr IsotopologueRatios() noexcept : data() {
    for (auto& x: data) x = std::numeric_limits<Numeric>::quiet_NaN();
  }
  
  constexpr Numeric operator[](const Index spec_ind) const ARTS_NOEXCEPT {
    ARTS_ASSERT(spec_ind < maxsize and spec_ind >= 0)
    return data[spec_ind];
  }
  
  constexpr Numeric operator[](const IsotopeRecord& ir) const {
    const Index spec_ind = find_species_index(ir);
    ARTS_USER_ERROR_IF(spec_ind < maxsize and spec_ind >= 0,
      "Cannot understand: ", ir.FullName(), " as a valid species")
    return data[spec_ind];
  }
  
  friend std::ostream& operator<<(std::ostream& os, const IsotopologueRatios& iso_rat) {
    for (size_t i=0; i<iso_rat.maxsize; i++) {
      if (i not_eq 0)
        os << '\n';
      os << Isotopologues[i].FullName() << ' ' << iso_rat.data[i];
    }
    return os;
  }
};
}

using SpeciesIsotopologueRatios = Species::IsotopologueRatios;

/** Check that isotopologue ratios for the given species are correctly defined. */
void checkIsotopologueRatios(const ArrayOfArrayOfSpeciesTag& abs_species,
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

/** Cross-section algorithm
 * 
 *  @param[in,out] xsec Cross section of one tag group. This is now the true attenuation cross section in units of m^2.
 *  @param[in,out] sourceCross section of one tag group. This is now the true source cross section in units of m^2.
 *  @param[in,out] phase Cross section of one tag group. This is now the true phase cross section in units of m^2.
 *  @param[in,out] dxsec Partial derivatives of xsec.
 *  @param[in,out] dsource Partial derivatives of source.
 *  @param[in,out] dphase Partial derivatives of phase.
 *  @param[in] jacobian_quantities As WSV
 *  @param[in] f_grid As WSV
 *  @param[in] abs_p As WSV
 *  @param[in] abs_t As WSV
 *  @param[in] abs_nlte As WSV
 *  @param[in] abs_vmrs As WSV
 *  @param[in] abs_species As WSV
 *  @param[in] band A single absorption band
 *  @param[in] isot_ratio Isotopologue ratio of this species
 *  @param[in] partfun_type Partition function type for this species
 *  @param[in] partfun_data Partition function model data for this species
 * 
 *  @author Richard Larsson
 *  @date   2019-10-10
 */
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
                  const Numeric& isot_ratio);

#endif  // absorption_h

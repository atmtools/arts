/* Copyright (C) 2018 Oliver Lemke <oliver.lemke@uni-hamburg.de>

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
  \file   hitran_xsec.h
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/

#include "hitran_xsec.h"

#include <algorithm>
#include <numeric>

#include "absorption.h"
#include "arts.h"
#include "check_input.h"
#include "debug.h"
#include "interpolation.h"
#include "physics_funcs.h"

String XsecRecord::SpeciesName() const {
  // The function species_name_from_species_index internally does an assertion
  // that the species with this index really exists.
  return Species::toShortName(mspecies);
}

void XsecRecord::SetVersion(const Index version) {
  if (version != 2) {
    ARTS_USER_ERROR("Invalid version, only 2 supported")
  }

  mversion = version;
}

void XsecRecord::Extract(VectorView result,
                         const Vector& f_grid,
                         const Numeric pressure,
                         const Numeric temperature,
                         const Verbosity& verbosity) const {
  CREATE_OUTS;

  const Index nf = f_grid.nelem();

  ARTS_ASSERT(result.nelem() == nf)

  result = 0.;

  const Index ndatasets = mfitcoeffs.nelem();
  for (Index this_dataset_i = 0; this_dataset_i < ndatasets; this_dataset_i++) {
    const Vector& data_f_grid = mfitcoeffs[this_dataset_i].get_numeric_grid(0);
    const Numeric data_fmin = data_f_grid[0];
    const Numeric data_fmax = data_f_grid[data_f_grid.nelem() - 1];

    if (out3.sufficient_priority()) {
      ostringstream os;
      os << "    f_grid:      " << f_grid[0] << " - " << f_grid[nf - 1]
         << " Hz\n"
         << "    data_f_grid: " << data_fmin << " - " << data_fmax << " Hz\n"
         << "    pressure: " << pressure << " K\n";
      out3 << os.str();
    }

    // We want to return result zero for all f_grid points that are outside the
    // data_f_grid, because xsec datasets are defined only where the absorption
    // was measured. So, we have to find out which part of f_grid is inside
    // data_f_grid.
    const Numeric* f_grid_begin = f_grid.get_c_array();
    const Numeric* f_grid_end = f_grid_begin + f_grid.nelem();
    const Index i_fstart = std::distance(
        f_grid_begin, std::lower_bound(f_grid_begin, f_grid_end, data_fmin));
    const Index i_fstop =
        std::distance(
            f_grid_begin,
            std::upper_bound(f_grid_begin + i_fstart, f_grid_end, data_fmax)) -
        1;

    // Ignore band if all frequencies are below or above data_f_grid:
    if (i_fstart == nf || i_fstop == -1) continue;

    const Index f_extent = i_fstop - i_fstart + 1;

    if (out3.sufficient_priority()) {
      ostringstream os;
      os << "    " << f_extent << " frequency extraction points starting at "
         << "frequency index " << i_fstart << ".\n";
      out3 << os.str();
    }

    // If f_extent is less than one, then the entire data_f_grid is between two
    // grid points of f_grid. (So that we do not have any f_grid points inside
    // data_f_grid.) Return also in this case.
    if (f_extent < 1) continue;

    // This is the part of f_grid that lies inside the xsec band
    const ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

    // Determine the index of the first grid points in the xsec band that lie
    // outside or on the boundary of the part of f_grid that lies inside the
    // xsec band. We can ignore the remaining data.
    const Numeric f_grid_fmin = f_grid[i_fstart];
    const Numeric f_grid_fmax = f_grid[i_fstop];

    const Numeric* data_f_grid_begin = data_f_grid.get_c_array();
    const Numeric* data_f_grid_end = data_f_grid_begin + data_f_grid.size() - 1;
    const Index i_data_fstart =
        std::distance(
            data_f_grid_begin,
            std::upper_bound(data_f_grid_begin, data_f_grid_end, f_grid_fmin)) -
        1;
    const Index i_data_fstop = std::distance(
        data_f_grid_begin,
        std::upper_bound(
            data_f_grid_begin + i_data_fstart, data_f_grid_end, f_grid_fmax));

    ARTS_ASSERT(i_data_fstart >= 0)
    ARTS_ASSERT(f_grid[i_fstart] > data_f_grid[i_data_fstart])
    ARTS_ASSERT(f_grid[i_fstart] < data_f_grid[i_data_fstart + 1])
    ARTS_ASSERT(f_grid[i_fstop] < data_f_grid[i_data_fstop])
    ARTS_ASSERT(f_grid[i_fstop] > data_f_grid[i_data_fstop - 1])

    // Extent for active data frequency vector:
    const Index data_f_extent = i_data_fstop - i_data_fstart + 1;

    // This is the part of the xsec dataset for which we have to do the
    // interpolation.
    const Range active_range(i_data_fstart, data_f_extent);
    const ConstVectorView data_f_grid_active = data_f_grid[active_range];

    Vector fit_result(data_f_grid.nelem());
    VectorView fit_result_active = fit_result[active_range];

    // We have to create a matching view on the result vector:
    VectorView result_active = result[Range(i_fstart, f_extent)];
    Vector xsec_interp(f_extent);

    CalcXsec(fit_result, this_dataset_i, pressure, temperature);

    RemoveNegativeXsec(fit_result);

    // Check if frequency is inside the range covered by the data:
    chk_interpolation_grids("Frequency interpolation for cross sections",
                            data_f_grid,
                            f_grid_active);

    {
      // Find frequency grid positions:
      ArrayOfGridPos f_gp(f_grid_active.nelem());
      gridpos(f_gp, data_f_grid_active, f_grid_active);

      Matrix itw(f_gp.nelem(), 2);
      interpweights(itw, f_gp);
      interp(xsec_interp, itw, fit_result_active, f_gp);
    }

    result_active += xsec_interp;
  }
}

void XsecRecord::CalcXsec(VectorView& xsec,
                          const Index dataset,
                          const Numeric pressure,
                          const Numeric temperature) const {
  for (Index i = 0; i < xsec.nelem(); i++) {
    const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
    xsec[i] = coeffs[P00] + coeffs[P10] * temperature + coeffs[P01] * pressure +
              coeffs[P20] * temperature * temperature;
  }
}

// void XsecRecord::CalcDT(VectorView& xsec_dt,
//                         const Index dataset,
//                         const Numeric temperature) const {
//   for (Index i = 0; i < xsec_dt.nelem(); i++) {
//     const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
//     xsec_dt[i] = coeffs[P10] + 2. * coeffs[P20] * temperature;
//   }
// }

// void XsecRecord::CalcDP(VectorView& xsec_dp,
//                         const Index dataset,
//                         const Numeric pressure) const {
//   for (Index i = 0; i < xsec_dp.nelem(); i++) {
//     const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
//     xsec_dp[i] = coeffs[P01] + 2. * coeffs[P02] * pressure;
//   }
// }

void XsecRecord::RemoveNegativeXsec(VectorView& xsec) const {
  Numeric sum_xsec{};
  Numeric sum_xsec_non_negative{};
  for (auto& x : xsec) {
    sum_xsec += x;
    if (x < 0.) x = 0.;
    sum_xsec_non_negative += x;
  }

  if (sum_xsec > 0. && sum_xsec != sum_xsec_non_negative) {
    xsec *= sum_xsec / sum_xsec_non_negative;
  }
}

/** Get the index in hitran_xsec_data for the given species.

 \param[in] hitran_xsec_data Hitran Xsec data array
 \param[in] species Species name

 \returns Index of this species in hitran_xsec_data. -1 if not found.
 */
Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            const Species::Species species) {
  for (Index i = 0; i < xsec_data.nelem(); i++)
    if (xsec_data[i].Species() == species) return i;

  return -1;
}

std::ostream& operator<<(std::ostream& os, const XsecRecord& xd) {
  os << "Species: " << xd.Species() << std::endl;
  return os;
}

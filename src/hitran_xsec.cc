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

#include "absorption.h"
#include "arts.h"
#include "check_input.h"
#include "interpolation.h"

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
                         const ConstVectorView& f_grid,
                         const Numeric pressure,
                         const Numeric temperature,
                         const Index extrapolate_p,
                         const Index extrapolate_t,
                         const Verbosity& verbosity) const {
  CREATE_OUTS;

  const Index nf = f_grid.nelem();

  // Assert that result vector has right size:
  ARTS_ASSERT(result.nelem() == nf)

  // Initialize result to zero (important for those frequencies outside the data grid).
  result = 0.;

  const Index ndatasets = mfitcoeffs.nelem();
  for (Index this_dataset_i = 0; this_dataset_i < ndatasets; this_dataset_i++) {
    const Vector& data_f_grid = mfitcoeffs[this_dataset_i].get_numeric_grid(0);
    const Numeric fmin = data_f_grid[0];
    const Numeric fmax = data_f_grid[data_f_grid.nelem() - 1];
    const Index data_nf = data_f_grid.nelem();

    if (out3.sufficient_priority()) {
      // Some detailed information to the most verbose output stream:
      ostringstream os;
      os << "    f_grid:      " << f_grid[0] << " - " << f_grid[nf - 1]
         << " Hz\n"
         << "    data_f_grid: " << fmin << " - " << fmax << " Hz\n"
         << "    pressure: " << pressure << " K\n";
      out3 << os.str();
    }

    // We want to return result zero for all f_grid points that are outside the
    // data_f_grid, because xsec datasets are defined only where the absorption
    // was measured. So, we have to find out which part of f_grid is inside
    // data_f_grid.
    Index i_fstart, i_fstop;

    for (i_fstart = 0; i_fstart < nf; ++i_fstart)
      if (f_grid[i_fstart] >= fmin) break;

    // Return immediately if all frequencies are below data_f_grid:
    if (i_fstart == nf) continue;

    for (i_fstop = nf - 1; i_fstop >= 0; --i_fstop)
      if (f_grid[i_fstop] <= fmax) break;

    // Return immediately if all frequencies are above data_f_grid:
    if (i_fstop == -1) continue;

    // Extent for active frequency vector:
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

    // This is the part of f_grid for which we have to do the interpolation.
    ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

    // We also need to determine the range in the xsec dataset that's inside
    // f_grid. We can ignore the remaining data.
    Index i_data_fstart, i_data_fstop;

    for (i_data_fstart = 0; i_data_fstart < data_nf; ++i_data_fstart)
      if (data_f_grid[i_data_fstart] >= fmin) break;

    for (i_data_fstop = data_nf - 1; i_data_fstop >= 0; --i_data_fstop)
      if (data_f_grid[i_data_fstop] <= fmax) break;

    // Extent for active data frequency vector:
    const Index data_f_extent = i_data_fstop - i_data_fstart + 1;

    // This is the part of f_grid for which we have to do the interpolation.
    ConstVectorView data_f_grid_active =
        data_f_grid[Range(i_data_fstart, data_f_extent)];

    // This is the part of the xsec dataset for which we have to do the
    // interpolation.
    Range active_range(i_data_fstart, data_f_extent);
    const ConstMatrixView coeffs_active =
        mfitcoeffs[this_dataset_i].data(active_range, joker);

    Vector fit_result(data_f_extent);
    Vector derivative;
    // Only allocate derivative if extrapolation is enabled
    if (extrapolate_p || extrapolate_t) {
      derivative.resize(data_f_extent);
    }

    // We have to create a matching view on the result vector:
    VectorView result_active = result[Range(i_fstart, f_extent)];
    Vector xsec_interp(f_extent);

    const Numeric min_p = mfitminpressures[this_dataset_i];
    const Numeric max_p = mfitmaxpressures[this_dataset_i];
    const Numeric min_t = mfitmintemperatures[this_dataset_i];
    const Numeric max_t = mfitmaxtemperatures[this_dataset_i];

    const Numeric active_temperature =
        temperature < min_t ? min_t
                            : (temperature > max_t ? max_t : temperature);

    const Numeric active_pressure =
        pressure < min_p ? min_p : (pressure > max_p ? max_p : pressure);

    CalcXsec(fit_result,
             this_dataset_i,
             active_range,
             active_pressure,
             active_temperature);

    if (extrapolate_p && (pressure < min_p || pressure > max_p)) {
      CalcDP(derivative,
             this_dataset_i,
             active_range,
             active_pressure,
             active_temperature);
      derivative *= pressure - active_pressure;
      fit_result += derivative;
    }

    if (extrapolate_t && (temperature < min_t || temperature > max_t)) {
      CalcDT(derivative,
             this_dataset_i,
             active_range,
             active_pressure,
             active_temperature);
      derivative *= temperature - active_temperature;
      fit_result += derivative;
    }

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
      interp(xsec_interp, itw, fit_result, f_gp);
    }

    result_active += xsec_interp;
  }
}

void XsecRecord::CalcXsec(VectorView& xsec,
                          const Index dataset,
                          const Range range,
                          const Numeric pressure,
                          const Numeric temperature) const {
  const Numeric logp = log10(pressure);
  for (Index i = range.get_start();
       i <= range.get_start() + range.get_extent() - 1;
       i++) {
    const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
    xsec[i] = coeffs[P00] + coeffs[P10] * temperature + coeffs[P01] * logp +
              coeffs[P20] * temperature * temperature +
              coeffs[P11] * temperature * logp + coeffs[P02] * logp * logp;
    xsec[i] *= xsec[i];
  }
}

void XsecRecord::CalcDT(VectorView& xsec_dt,
                        const Index dataset,
                        const Range range,
                        const Numeric pressure,
                        const Numeric temperature) const {
  const Numeric logp = log10(pressure);
  for (Index i = range.get_start();
       i <= range.get_start() + range.get_extent() - 1;
       i++) {
    const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
    xsec_dt[i] =
        2. *
        (coeffs[P10] + 2. * coeffs[P20] * temperature + coeffs[P11] * logp) *
        (coeffs[P00] + coeffs[P10] * temperature +
         coeffs[P20] * temperature * temperature + coeffs[P01] * logp +
         coeffs[P11] * temperature * logp + coeffs[P02] * logp * logp);
  }
}

void XsecRecord::CalcDP(VectorView& xsec_dp,
                        const Index dataset,
                        const Range range,
                        const Numeric pressure,
                        const Numeric temperature) const {
  const Numeric logp = log10(pressure);
  const Numeric plog = pressure * log(10);

  for (Index i = range.get_start();
       i <= range.get_start() + range.get_extent() - 1;
       i++) {
    const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
    xsec_dp[i] = 2. *
                 (coeffs[P01] / plog + coeffs[P11] * temperature / plog +
                  2. * coeffs[P02] * logp / plog) *
                 (coeffs[P00] + coeffs[P10] * temperature +
                  coeffs[P20] * temperature * temperature + coeffs[P01] * logp +
                  coeffs[P11] * temperature * logp + coeffs[P02] * logp * logp);
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

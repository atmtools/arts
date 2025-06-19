/*!
  \file   xsec_fit.h
  \author Oliver Lemke <oliver.lemke@uni-hamburg.de>
  \date   2018-01-08

  \brief  Methods and classes for HITRAN absorption cross section data.
*/

#include "xsec_fit.h"

#include <algorithm>
#include <memory>

#include "debug.h"
#include "interp.h"

const SpeciesEnum& XsecRecord::Species() const { return mspecies; };

/** Set species name */
void XsecRecord::SetSpecies(const SpeciesEnum species) { mspecies = species; };

const Vector& XsecRecord::FitMinPressures() const { return mfitminpressures; };

/** Get maximum pressures from fit */
const Vector& XsecRecord::FitMaxPressures() const { return mfitmaxpressures; };

/** Get mininum temperatures from fit */
const Vector& XsecRecord::FitMinTemperatures() const {
  return mfitmintemperatures;
};

/** Get maximum temperatures */
const Vector& XsecRecord::FitMaxTemperatures() const {
  return mfitmaxtemperatures;
};

/** Get coefficients */
const ArrayOfGriddedField1Named& XsecRecord::FitCoeffs() const {
  return mfitcoeffs;
};

Vector& XsecRecord::FitMinPressures() { return mfitminpressures; };

/** Get maximum pressures from fit */
Vector& XsecRecord::FitMaxPressures() { return mfitmaxpressures; };

/** Get mininum temperatures from fit */
Vector& XsecRecord::FitMinTemperatures() { return mfitmintemperatures; };

/** Get maximum temperatures */
Vector& XsecRecord::FitMaxTemperatures() { return mfitmaxtemperatures; };

/** Get coefficients */
ArrayOfGriddedField1Named& XsecRecord::FitCoeffs() { return mfitcoeffs; };

namespace {
void RemoveNegativeXsec(Vector& xsec) {
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
}  // namespace

String XsecRecord::SpeciesName() const {
  // The function species_name_from_species_index internally does an assertion
  // that the species with this index really exists.
  return String{toString<1>(mspecies)};
}

void XsecRecord::SetVersion(const Index version) {
  if (version != mversion) {
    ARTS_USER_ERROR(
        "Invalid version {}, only version {} supported", version, mversion)
  }
}

void XsecRecord::Extract(VectorView result,
                         const Vector& f_grid,
                         const Numeric pressure,
                         const Numeric temperature) const {
  const Size nf = f_grid.size();

  assert(result.size() == nf);

  result = 0.;

  const Size ndatasets = mfitcoeffs.size();
  for (Size this_dataset_i = 0; this_dataset_i < ndatasets; this_dataset_i++) {
    const Vector& data_f_grid = mfitcoeffs[this_dataset_i].grid<0>();
    const Numeric data_fmin   = data_f_grid[0];
    const Numeric data_fmax   = data_f_grid[data_f_grid.size() - 1];

    // We want to return result zero for all f_grid points that are outside the
    // data_f_grid, because xsec datasets are defined only where the absorption
    // was measured. So, we have to find out which part of f_grid is inside
    // data_f_grid.
    const Numeric* f_grid_begin = f_grid.data_handle();
    const Numeric* f_grid_end   = f_grid_begin + f_grid.size();
    const Index i_fstart        = std::distance(
        f_grid_begin, std::lower_bound(f_grid_begin, f_grid_end, data_fmin));
    const Index i_fstop =
        std::distance(
            f_grid_begin,
            std::upper_bound(f_grid_begin + i_fstart, f_grid_end, data_fmax)) -
        1;

    // Ignore band if all frequencies are below or above data_f_grid:
    if (static_cast<Size>(i_fstart) == nf || i_fstop == -1) continue;

    const Index f_extent = i_fstop - i_fstart + 1;

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

    const Numeric* data_f_grid_begin = data_f_grid.data_handle();
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

    assert(i_data_fstart >= 0);
    assert(f_grid[i_fstart] >= data_f_grid[i_data_fstart]);
    assert(f_grid[i_fstart] < data_f_grid[i_data_fstart + 1]);
    assert(f_grid[i_fstop] <= data_f_grid[i_data_fstop]);
    assert(f_grid[i_fstop] > data_f_grid[i_data_fstop - 1]);

    // Extent for active data frequency vector:
    const Index data_f_extent = i_data_fstop - i_data_fstart + 1;

    // This is the part of the xsec dataset for which we have to do the
    // interpolation.
    const Range active_range(i_data_fstart, data_f_extent);
    const ConstVectorView data_f_grid_active = data_f_grid[active_range];

    Vector fit_result(data_f_grid.size());
    VectorView fit_result_active = fit_result[active_range];

    // We have to create a matching view on the result vector:
    VectorView result_active = result[Range(i_fstart, f_extent)];
    Vector xsec_interp(f_extent);

    CalcXsec(fit_result, this_dataset_i, pressure, temperature);

    RemoveNegativeXsec(fit_result);

    {
      const auto f_gp =
          my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<1>>(
              data_f_grid_active, f_grid_active, 0.5, "Frequency");
      const auto f_itw = interpweights(f_gp);
      // Find frequency grid positions:
      my_interp::reinterp(xsec_interp, fit_result_active, f_itw, f_gp);
    }

    result_active += xsec_interp;
  }
}

void XsecRecord::CalcXsec(VectorView xsec,
                          const Index dataset,
                          const Numeric pressure,
                          const Numeric temperature) const {
  for (Size i = 0; i < xsec.size(); i++) {
    const ConstVectorView coeffs = mfitcoeffs[dataset].data[i, joker];
    xsec[i] = coeffs[P00] + coeffs[P10] * temperature + coeffs[P01] * pressure +
              coeffs[P20] * temperature * temperature;
  }
}

// void XsecRecord::CalcDT(VectorView xsec_dt,
//                         const Index dataset,
//                         const Numeric temperature) const {
//   for (Index i = 0; i < xsec_dt.size(); i++) {
//     const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
//     xsec_dt[i] = coeffs[P10] + 2. * coeffs[P20] * temperature;
//   }
// }

// void XsecRecord::CalcDP(VectorView xsec_dp,
//                         const Index dataset,
//                         const Numeric pressure) const {
//   for (Index i = 0; i < xsec_dp.size(); i++) {
//     const ConstVectorView coeffs = mfitcoeffs[dataset].data(i, joker);
//     xsec_dp[i] = coeffs[P01] + 2. * coeffs[P02] * pressure;
//   }
// }

/** Get the index in xsec_fit_data for the given species.

 \param[in] xsec_fit_data Hitran Xsec data array
 \param[in] species Species name

 \returns Index of this species in xsec_fit_data. -1 if not found.
 */
Index hitran_xsec_get_index(const ArrayOfXsecRecord& xsec_data,
                            const SpeciesEnum species) {
  for (Size i = 0; i < xsec_data.size(); i++)
    if (xsec_data[i].Species() == species) return i;

  return -1;
}

/** Get the index in xsec_fit_data for the given species.

 \param[in] xsec_fit_data Hitran Xsec data array
 \param[in] species Species name

 \returns Correct CIA record or nullptr if not found.
 */
XsecRecord* hitran_xsec_get_data(
    const std::shared_ptr<std::vector<XsecRecord>>& xsec_data,
    const SpeciesEnum species) {
  for (auto& xsec : *xsec_data) {
    if (xsec.Species() == species) return &xsec;
  }
  return nullptr;
}

std::ostream& operator<<(std::ostream& os, const XsecRecord& xd) {
  os << "Species: " << xd.Species() << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const ArrayOfXsecRecord& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}

/* Copyright (C) 2012 Stefan Buehler <sbuehler@ltu.se>
 
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
 Implementation file for work with HITRAN collision induced absorption (CIA).
 
 The CIA data are part of the HITRAN distribution. They are described in
 Richard, C., I. E. Gordon, L. S. Rothman, M. Abel, L. Frommhold, M. Gustafsson,
 J.-M. Hartmann, C. Hermans, W. J. Lafferty, G. S. Orton, K.M. Smith, and H. Tran (2012),
 New section of the HITRAN database: Collision-induced absorption (CIA),
 J. Quant. Spectrosc. Radiat. Transfer, 113, 1276-1285, doi:10.1016/j.jqsrt.2011.11.004.
 
 \author Stefan Buehler
 \date   2012-11-30
 */

#include "cia.h"
#include <cmath>
#include "abs_species_tags.h"
#include "absorption.h"
#include "file.h"
#include "interpolation_poly.h"

extern const Numeric SPEED_OF_LIGHT;

/** Interpolate CIA data.
 
 Interpolate CIA data to given frequency vector and given scalar temperature.
 Uses third order interpolation in both coordinates, if grid length allows,
 otherwise lower order or no interpolation.
 
 \param[out] result     CIA value for given frequency grid and temperature.
 \param[in] f_grid      Frequency grid.
 \param[in] temperature Scalar temperature.
 \param[in] cia_data    The CIA dataset to interpolate.
 \param[in] robust      Set to 1 to suppress runtime errors (and return NAN values instead).
 \param[in] verbosity   Standard verbosity object.
 */
void cia_interpolation(VectorView result,
                       ConstVectorView f_grid,
                       const Numeric& temperature,
                       const GriddedField2& cia_data,
                       const Numeric& T_extrapolfac,
                       const Index& robust,
                       const Verbosity& verbosity) {
  CREATE_OUTS;

  const Index nf = f_grid.nelem();

  // Assert that result vector has right size:
  assert(result.nelem() == nf);

  // Get data grids:
  ConstVectorView data_f_grid = cia_data.get_numeric_grid(0);
  ConstVectorView data_T_grid = cia_data.get_numeric_grid(1);

  if (out3.sufficient_priority()) {
    // Some detailed information to the most verbose output stream:
    ostringstream os;
    os << "    f_grid:      " << f_grid[0] << " - " << f_grid[nf - 1] << " Hz\n"
       << "    data_f_grid: " << data_f_grid[0] << " - "
       << data_f_grid[data_f_grid.nelem() - 1] << " Hz\n"
       << "    temperature: " << temperature << " K\n"
       << "    data_T_grid: " << data_T_grid[0] << " - "
       << data_T_grid[data_T_grid.nelem() - 1] << " K\n";
    out3 << os.str();
  }

  // Initialize result to zero (important for those frequencies outside the data grid).
  result = 0;

  // We want to return result zero for all f_grid points that are outside the
  // data_f_grid, because some CIA datasets are defined only where the absorption
  // is not zero. So, we have to find out which part of f_grid is inside
  // data_f_grid.
  Index i_fstart, i_fstop;

  for (i_fstart = 0; i_fstart < nf; ++i_fstart)
    if (f_grid[i_fstart] >= data_f_grid[0]) break;

  // Return directly if all frequencies are below data_f_grid:
  if (i_fstart == nf) return;

  for (i_fstop = nf - 1; i_fstop >= 0; --i_fstop)
    if (f_grid[i_fstop] <= data_f_grid[data_f_grid.nelem() - 1]) break;

  // Return directly if all frequencies are above data_f_grid:
  if (i_fstop == -1) return;

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
  if (f_extent < 1) return;

  // This is the part of f_grid for which we have to do the interpolation.
  ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

  // We have to create a matching view on the result vector:
  VectorView result_active = result[Range(i_fstart, f_extent)];

  // Decide on interpolation orders:
  const Index f_order = 3;

  // The frequency grid has to have enough points for this interpolation
  // order, otherwise throw a runtime error.
  if (data_f_grid.nelem() < f_order + 1) {
    ostringstream os;
    os << "Not enough frequency grid points in CIA data.\n"
       << "You have only " << data_f_grid.nelem() << " grid points.\n"
       << "But need at least " << f_order + 1 << ".";
    throw runtime_error(os.str());
  }

  // For T we have to be adaptive, since sometimes there is only one T in
  // the data
  Index T_order;
  switch (data_T_grid.nelem()) {
    case 1:
      T_order = 0;
      break;
    case 2:
      T_order = 1;
      break;
    case 3:
      T_order = 2;
      break;
    default:
      T_order = 3;
      break;
  }

  // Check if frequency is inside the range covered by the data:
  chk_interpolation_grids("Frequency interpolation for CIA continuum",
                          data_f_grid,
                          f_grid_active,
                          f_order);

  // Check if temperature is inside the range covered by the data:
  if (T_order > 0) {
    try {
      chk_interpolation_grids("Temperature interpolation for CIA continuum",
                              data_T_grid,
                              temperature,
                              T_order,
                              T_extrapolfac);
    } catch (const std::runtime_error& e) {
      //            cout << "Gotcha!\n";
      if (robust) {
        // Just return NANs, but continue.
        result_active = NAN;
        return;
      } else {
        // Re-throw the exception.
        throw runtime_error(e.what());
      }
    }
  }

  // Find frequency grid positions:
  ArrayOfGridPosPoly f_gp(f_grid_active.nelem()), T_gp(1);
  gridpos_poly(f_gp, data_f_grid, f_grid_active, f_order);

  // Do the rest of the interpolation.
  if (T_order == 0) {
    // No temperature interpolation in this case, just a frequency interpolation.

    Matrix itw(f_gp.nelem(), f_order + 1);
    interpweights(itw, f_gp);
    interp(result_active, itw, cia_data.data(joker, 0), f_gp);
  } else {
    // Temperature and frequency interpolation.

    // Find temperature grid position:
    gridpos_poly(T_gp, data_T_grid, temperature, T_order, T_extrapolfac);

    // Calculate combined interpolation weights:
    Tensor3 itw(f_gp.nelem(), T_gp.nelem(), (f_order + 1) * (T_order + 1));
    interpweights(itw, f_gp, T_gp);

    // Make a matrix view of the results vector:
    MatrixView result_matrix(result_active);

    //        cout << "result_matrix r/c: " << result_matrix.nrows() << " / "
    //             << result_matrix.ncols() << endl;
    //        cout << "result_matrix before: " << result_matrix << endl;

    //        cout << "cia_data: " << cia_data.data << endl;

    // Actual interpolation:
    interp(result_matrix, itw, cia_data.data, f_gp, T_gp);

    //        cout << "result_matrix after: " << result_matrix << endl;
  }

  //    cout << "result_active before: " << result_active << endl;

  // Set negative values to zero. (These could happen due to overshooting
  // of the higher order interpolation.)
  for (Index i = 0; i < result_active.nelem(); ++i)
    if (result_active[i] < 0) result_active[i] = 0;

  //    cout << "result_active after: " << result_active << endl;
}

/** Get the index in cia_data for the two given species.
 
 \param[in] cia_data CIA data array
 \param[in] sp1 First species index
 \param[in] sp2 Second species index

 \returns Index of this species combination in cia_data. -1 if not found.
 */
Index cia_get_index(const ArrayOfCIARecord& cia_data,
                    const Index sp1,
                    const Index sp2) {
  for (Index i = 0; i < cia_data.nelem(); i++)
    if ((cia_data[i].Species(0) == sp1 && cia_data[i].Species(1) == sp2) ||
        (cia_data[i].Species(0) == sp2 && cia_data[i].Species(1) == sp1))
      return i;

  return -1;
}

// Documentation in header file.
void CIARecord::Extract(VectorView result,
                        ConstVectorView f_grid,
                        const Numeric& temperature,
                        const Index& dataset,
                        const Numeric& T_extrapolfac,
                        const Index& robust,
                        const Verbosity& verbosity) const {
  // If there is more than one dataset available for this species pair,
  // we have to decide on which one to use. The rest is done by helper function
  // cia_interpolate.

  // Make sure dataset index exists
  if (dataset >= mdata.nelem()) {
    ostringstream os;
    os << "There are only " << mdata.nelem() << " datasets in this CIA file.\n"
       << "But you are trying to use dataset " << dataset
       << ". (Zero-based indexing.)";
    throw runtime_error(os.str());
  }

  // Get a handle on this dataset:
  const GriddedField2& this_cia = mdata[dataset];

  cia_interpolation(
      result, f_grid, temperature, this_cia, T_extrapolfac, robust, verbosity);
}

// Documentation in header file.
String CIARecord::MoleculeName(const Index i) const {
  // Assert that i is 0 or 1:
  assert(i >= 0);
  assert(i <= 1);

  // The function species_name_from_species_index internally does an assertion
  // that the species with this index really exists.
  return species_name_from_species_index(mspecies[i]);
}

// Documentation in header file.
void CIARecord::SetMoleculeName(const Index i, const String& name) {
  // Assert that i is 0 or 1:
  assert(i >= 0);
  assert(i <= 1);

  // Find out the species index for name:
  Index spec_ind = species_index_from_species_name(name);

  // Function species_index_from_species_name returns -1 if the species does
  // not exist. Check this:
  if (spec_ind < 0) {
    ostringstream os;
    os << "Species does not exist in ARTS: " << name;
    throw runtime_error(os.str());
  }

  // Assign species:
  mspecies[i] = spec_ind;
}

//! Read CIA catalog file.
/*!
 Reads the given CIA catalog file into this CIARecord.

 \param[in]  filename  Path of catalog file to read.
 \param[in]  verbosity.
 \return os
 */
void CIARecord::ReadFromCIA(const String& filename,
                            const Verbosity& verbosity) {
  CREATE_OUT2;

  ifstream is;

  out2 << "  Reading file: " << filename << "\n";
  open_input_file(is, filename);

  // Number of points for spectral range in current dataset
  Index npoints = -1;

  // Min/max wave numbers for this dataset
  Numeric wave_min = -1.;
  Numeric wave_max = -1.;

  // Current dataset index
  Index ndataset = -1;

  // Current set in dataset
  Index nset = -1;

  // Frequency, temp and cia values for current dataset
  Vector freq;
  ArrayOfNumeric temp;
  ArrayOfVector cia;

  // Keep track of current line in file
  Index nline = 0;

  mdata.resize(0);
  istringstream istr;

  while (is) {
    String line;

    // Extract needed information from dataset header line
    //////////////////////////////////////////////////////
    nline++;
    getline(is, line);
    if (is.eof()) continue;

    if (line.nelem() < 100) {
      ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << endl
         << "Header line unexpectedly short: " << endl
         << line;

      throw runtime_error(os.str());
    }

    if (is.bad()) {
      ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << endl;

      throw runtime_error(os.str());
    }

    line.erase(0, 20);

    // Data for current set
    Index set_npoints;
    Numeric set_temp = NAN;
    Numeric set_wave_min = NAN;
    Numeric set_wave_max = NAN;

    istr.str(line);
    istr.clear();
    istr >> set_wave_min >> set_wave_max >> set_npoints >> set_temp;

    if (!istr || std::isnan(set_temp) || std::isnan(set_wave_min) ||
        std::isnan(set_wave_max)) {
      ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << endl;

      throw runtime_error(os.str());
    }

    // If the min/max wave numbers of this set are different from the
    // previous one, a new dataset starts
    if (npoints == -1 || wave_min != set_wave_min || wave_max != set_wave_max) {
      if (ndataset != -1) AppendDataset(freq, temp, cia);

      npoints = set_npoints;
      ndataset++;
      wave_min = set_wave_min;
      wave_max = set_wave_max;
      nset = 0;
      freq.resize(set_npoints);
      temp.resize(0);
      cia.resize(0);
    }
    if (npoints != set_npoints) {
      ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << endl
         << "Inconsistent number of data points. Expected " << npoints
         << ", got " << set_npoints;

      throw runtime_error(os.str());
    }

    temp.push_back(set_temp);
    cia.push_back(Vector(npoints));

    // Read dataset
    ////////////////////
    for (Index i = 0; i < npoints; i++) {
      Numeric w, c;

      nline++;
      getline(is, line);

      w = c = NAN;
      istr.str(line);
      istr.clear();
      istr >> w >> c;

      if (std::isnan(w) || std::isnan(c) || is.bad() || istr.bad()) {
        ostringstream os;
        os << "Error in line " << nline << " reading CIA catalog file "
           << filename << ":" << endl
           << line;

        throw runtime_error(os.str());
      }

      // Convert wavenumbers to Herz:
      freq[i] = 100. * w * SPEED_OF_LIGHT;

      // Convert binary absorption cross-sections from
      // cm^5 molec^(-2) to m^5 molec^(-2):
      cia[nset][i] = c / 1e10;
    }

    nset++;
  }

  if (is.bad()) {
    ostringstream os;
    os << "Error in line " << nline << " reading CIA catalog file " << filename
       << endl;

    throw runtime_error(os.str());
  }

  AppendDataset(freq, temp, cia);

  //    // For debugging
  //    for(Index i = 0; i < mdata.nelem(); i++)
  //    {
  //        cout << i << " ";
  //        cout << mdata[i].get_numeric_grid(0).nelem() << " ";
  //        cout << mdata[i].get_numeric_grid(1).nelem() << endl;
  //    }
}

/** Append data dataset to mdata. */
void CIARecord::AppendDataset(const Vector& freq,
                              const ArrayOfNumeric& temp,
                              const ArrayOfVector& cia) {
  GriddedField2 dataset;
  dataset.resize(freq.nelem(), temp.nelem());
  dataset.set_grid(0, freq);
  dataset.set_grid_name(0, "Frequency");

  Vector temp_t;
  temp_t = temp;
  dataset.set_grid(1, temp_t);
  dataset.set_grid_name(1, "Temperature");

  for (Index t = 0; t < temp.nelem(); t++) dataset.data(joker, t) = cia[t];
  mdata.push_back(dataset);
}

/** Append other CIARecord to this. */
void CIARecord::AppendDataset(const CIARecord& c2) {
  for (Index ii = 0; ii < c2.DatasetCount(); ii++) {
    mdata.push_back(c2.Dataset(ii));
  }
}

//! Output operator for CIARecord
/*!
 Outputs the grids for the given CIARecord.

 \param[in,out]  os  Output stream.
 \param[in]      cr  CIARecord.
 \return os
 */
ostream& operator<<(ostream& os, const CIARecord& /* cr */) {
  os << "CIARecord output operator not yet implemented." << endl;
  return os;
}

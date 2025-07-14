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
#include <memory>

#include "arts_constants.h"
#include "double_imanip.h"
#include "file.h"
#include "interp.h"

inline constexpr Numeric SPEED_OF_LIGHT = Constant::speed_of_light;

CIARecord::CIARecord(ArrayOfGriddedField2 data,
                     SpeciesEnum spec1,
                     SpeciesEnum spec2)
    : mdata(std::move(data)), mspecies({spec1, spec2}) {}

std::array<SpeciesEnum, 2> CIARecord::TwoSpecies() const { return mspecies; }

std::array<SpeciesEnum, 2>& CIARecord::TwoSpecies() { return mspecies; }

Numeric CIARecord::Extract(const Numeric& frequency,
                           const Numeric& temperature,
                           const Numeric& T_extrapolfac,
                           const Index& robust) const {
  Vector result(1);
  const Vector freqvec(1, frequency);

  Extract(result, freqvec, temperature, T_extrapolfac, robust);

  return result[0];
}

Index CIARecord::DatasetCount() const { return mdata.size(); }

const ArrayOfGriddedField2& CIARecord::Data() const { return mdata; }

/** Return CIA data.
   */
ArrayOfGriddedField2& CIARecord::Data() { return mdata; }

/** Set CIA species.
     \param[in] first CIA Species.
     \param[in] second CIA Species.
     */
void CIARecord::SetSpecies(const SpeciesEnum first, const SpeciesEnum second) {
  mspecies[0] = first;
  mspecies[1] = second;
}

SpeciesEnum CIARecord::Species(const Index i) const {
  // Assert that i is 0 or 1:
  assert(i >= 0);
  assert(i <= 1);

  return mspecies[i];
}

ConstVectorView CIARecord::FrequencyGrid(Size dataset) const {
  assert(dataset < mdata.size());

  return mdata[dataset].grid<0>();
}

ConstVectorView CIARecord::TemperatureGrid(Size dataset) const {
  assert(dataset < mdata.size());

  return mdata[dataset].grid<1>();
}

const GriddedField2& CIARecord::Dataset(Size dataset) const {
  assert(dataset < mdata.size());

  return mdata[dataset];
}

std::ostream& operator<<(std::ostream& os, const ArrayOfCIARecord& x) {
  for (auto& a : x) os << a << '\n';
  return os;
}

/** Interpolate CIA data.
 
 Interpolate CIA data to given frequency vector and given scalar temperature.
 Uses third order interpolation in both coordinates, if grid length allows,
 otherwise lower order or no interpolation.
 
 \param[out] result     CIA value for given frequency grid and temperature.
 \param[in] f_grid      Frequency grid.
 \param[in] temperature Scalar temperature.
 \param[in] cia_data    The CIA dataset to interpolate.
 \param[in] robust      Set to 1 to suppress runtime errors (and return NAN values instead).
 */
void cia_interpolation(VectorView result,
                       const ConstVectorView& f_grid,
                       const Numeric& temperature,
                       const GriddedField2& cia_data,
                       const Numeric& T_extrapolfac,
                       const Index& robust) try {
  const Index nf = f_grid.size();

  // Assert that result vector has right size:
  assert(result.size() == static_cast<Size>(nf));

  // Get data grids:
  ConstVectorView data_f_grid = cia_data.grid<0>();
  ConstVectorView data_T_grid = cia_data.grid<1>();

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
    if (f_grid[i_fstop] <= data_f_grid[data_f_grid.size() - 1]) break;

  // Return directly if all frequencies are above data_f_grid:
  if (i_fstop == -1) return;

  // Extent for active frequency vector:
  const Index f_extent = i_fstop - i_fstart + 1;

  // If f_extent is less than one, then the entire data_f_grid is between two
  // grid points of f_grid. (So that we do not have any f_grid points inside
  // data_f_grid.) Return also in this case.
  if (f_extent < 1) return;

  // This is the part of f_grid for which we have to do the interpolation.
  ConstVectorView f_grid_active = f_grid[Range(i_fstart, f_extent)];

  // We have to create a matching view on the result vector:
  VectorView result_active = result[Range(i_fstart, f_extent)];

  // Decide on interpolation orders:
  constexpr Index f_order = 3;

  // The frequency grid has to have enough points for this interpolation
  // order, otherwise throw a runtime error.
  if (data_f_grid.size() < f_order + 1) {
    std::ostringstream os;
    os << "Not enough frequency grid points in CIA data.\n"
       << "You have only " << data_f_grid.size() << " grid points.\n"
       << "But need at least " << f_order + 1 << ".";
    throw std::runtime_error(os.str());
  }

  // For T we have to be adaptive, since sometimes there is only one T in
  // the data
  Index T_order;
  switch (data_T_grid.size()) {
    case 1:  T_order = 0; break;
    case 2:  T_order = 1; break;
    case 3:  T_order = 2; break;
    default: T_order = 3; break;
  }

  // Find frequency grid positions:
  const auto f_lag = my_interp::lagrange_interpolation_list<
      FixedLagrangeInterpolation<f_order>>(
      f_grid_active, data_f_grid, 0.5, "Frequency");

  // Do the rest of the interpolation.
  if (T_order == 0) {
    // No temperature interpolation in this case, just a frequency interpolation.
    result_active =
        reinterp(cia_data.data[joker, 0], interpweights(f_lag), f_lag);
  } else {
    // Temperature and frequency interpolation.
    const auto Tnew = matpack::cdata_t<Numeric, 1>{temperature};
    if (T_order == 1) {
      const auto T_lag =
          my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<1>>(
              Tnew, data_T_grid, T_extrapolfac, "Temperature");
      result_active =
          reinterp(cia_data.data, interpweights(f_lag, T_lag), f_lag, T_lag)
              .reshape(f_lag.size());
    } else if (T_order == 2) {
      const auto T_lag =
          my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<2>>(
              Tnew, data_T_grid, T_extrapolfac, "Temperature");
      result_active =
          reinterp(cia_data.data, interpweights(f_lag, T_lag), f_lag, T_lag)
              .reshape(f_lag.size());
    } else if (T_order == 3) {
      const auto T_lag =
          my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<3>>(
              Tnew, data_T_grid, T_extrapolfac, "Temperature");
      result_active =
          reinterp(cia_data.data, interpweights(f_lag, T_lag), f_lag, T_lag)
              .reshape(f_lag.size());
    } else {
      throw std::runtime_error(
          "Cannot have this T_order, you must update the code...");
    }
  }

  //    cout << "result_active before: " << result_active << '\n';

  // Set negative values to zero. (These could happen due to overshooting
  // of the higher order interpolation.)
  for (Size i = 0; i < result_active.size(); ++i)
    if (result_active[i] < 0) result_active[i] = 0;

  //    cout << "result_active after: " << result_active << '\n';
} catch (const std::exception&) {
  if (robust) {
    result = NAN;
  } else {
    throw;
  }
}

/** Get the index in cia_data for the two given species.
 
 \param[in] cia_data CIA data array
 \param[in] sp1 First species index
 \param[in] sp2 Second species index

 \returns Index of this species combination in cia_data. -1 if not found.
 */
Index cia_get_index(const ArrayOfCIARecord& cia_data,
                    const SpeciesEnum sp1,
                    const SpeciesEnum sp2) {
  for (Size i = 0; i < cia_data.size(); i++)
    if ((cia_data[i].Species(0) == sp1 && cia_data[i].Species(1) == sp2) ||
        (cia_data[i].Species(0) == sp2 && cia_data[i].Species(1) == sp1))
      return i;

  return -1;
}

/** Get the correct cia recrod
 
 \param[in] cia_data CIA data array
 \param[in] sp1 First species index
 \param[in] sp2 Second species index

 \returns Correct CIA record or nullptr if not found.
 */
CIARecord* cia_get_data(const std::shared_ptr<std::vector<CIARecord>>& cia_data,
                        const SpeciesEnum sp1,
                        const SpeciesEnum sp2) {
  for (auto& data : *cia_data)
    if ((data.Species(0) == sp1 && data.Species(1) == sp2) ||
        (data.Species(0) == sp2 && data.Species(1) == sp1))
      return &data;

  return nullptr;
}

// Documentation in header file.
void CIARecord::Extract(VectorView res,
                        const ConstVectorView& f_grid,
                        const Numeric& temperature,
                        const Numeric& T_extrapolfac,
                        const Index& robust) const {
  res = 0;

  Vector result(res.size());
  for (auto& this_cia : mdata) {
    cia_interpolation(
        result, f_grid, temperature, this_cia, T_extrapolfac, robust);
    res += result;
  }
}

// Documentation in header file.
String CIARecord::MoleculeName(const Index i) const {
  // Assert that i is 0 or 1:
  assert(i == 0 or i == 1);

  // The function species_name_from_species_index internally does an assertion
  // that the species with this index really exists.
  return String{toString<1>(mspecies[i])};
}

// Documentation in header file.
void CIARecord::SetMoleculeName(const Index i, const String& name) {
  // Assert that i is 0 or 1:
  assert(i >= 0);
  assert(i <= 1);

  // Find out the species index for name:
  SpeciesEnum spec_ind = to<SpeciesEnum>(name);

  // Function species_index_from_species_name returns -1 if the species does
  // not exist. Check this:
  ARTS_USER_ERROR_IF(
      not good_enum(spec_ind), "Species does not exist in ARTS: {}", name)

  // Assign species:
  mspecies[i] = spec_ind;
}

//! Read CIA catalog file.
/*!
 Reads the given CIA catalog file into this CIARecord.

 \param[in]  filename  Path of catalog file to read.
 \return os
 */
void CIARecord::ReadFromCIA(const String& filename) {
  std::ifstream is;

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
  std::vector<Numeric> temp;
  ArrayOfVector cia;

  // Keep track of current line in file
  Index nline = 0;

  mdata.resize(0);
  std::istringstream istr;

  while (is) {
    String line;

    // Extract needed information from dataset header line
    //////////////////////////////////////////////////////
    nline++;
    getline(is, line);
    if (is.eof()) continue;

    if (line.size() < 100) {
      std::ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << '\n'
         << "Header line unexpectedly short: " << '\n'
         << line;

      throw std::runtime_error(os.str());
    }

    if (is.bad()) {
      std::ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << '\n';

      throw std::runtime_error(os.str());
    }

    line.erase(0, 20);

    // Data for current set
    Index set_npoints;
    Numeric set_temp     = NAN;
    Numeric set_wave_min = NAN;
    Numeric set_wave_max = NAN;

    istr.str(line);
    istr.clear();
    istr >> double_imanip() >> set_wave_min >> set_wave_max;
    istr >> set_npoints;
    istr >> double_imanip() >> set_temp;

    if (!istr || std::isnan(set_temp) || std::isnan(set_wave_min) ||
        std::isnan(set_wave_max)) {
      std::ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << '\n';

      throw std::runtime_error(os.str());
    }

    // If the min/max wave numbers of this set are different from the
    // previous one, a new dataset starts
    if (npoints == -1 || wave_min != set_wave_min || wave_max != set_wave_max) {
      if (ndataset != -1) AppendDataset(freq, Vector{temp}, cia);

      npoints = set_npoints;
      ndataset++;
      wave_min = set_wave_min;
      wave_max = set_wave_max;
      nset     = 0;
      freq.resize(set_npoints);
      temp.resize(0);
      cia.resize(0);
    }
    if (npoints != set_npoints) {
      std::ostringstream os;
      os << "Error in line " << nline << " reading CIA catalog file "
         << filename << '\n'
         << "Inconsistent number of data points. Expected " << npoints
         << ", got " << set_npoints;

      throw std::runtime_error(os.str());
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
        std::ostringstream os;
        os << "Error in line " << nline << " reading CIA catalog file "
           << filename << ":" << '\n'
           << line;

        throw std::runtime_error(os.str());
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
    std::ostringstream os;
    os << "Error in line " << nline << " reading CIA catalog file " << filename
       << '\n';

    throw std::runtime_error(os.str());
  }

  AppendDataset(freq, Vector{temp}, cia);

  //    // For debugging
  //    for(Index i = 0; i < mdata.size(); i++)
  //    {
  //        cout << i << " ";
  //        cout << mdata[i].get_numeric_grid(0).size() << " ";
  //        cout << mdata[i].get_numeric_grid(1).size() << '\n';
  //    }
}

/** Append data dataset to mdata. */
void CIARecord::AppendDataset(const Vector& freq,
                              const Vector& temp,
                              const ArrayOfVector& cia) {
  GriddedField2 dataset;
  dataset.resize(freq.size(), temp.size());
  dataset.grid<0>()     = freq;
  dataset.gridname<0>() = "Frequency";

  Vector temp_t;
  temp_t                = temp;
  dataset.grid<1>()     = temp_t;
  dataset.gridname<1>() = "Temperature";

  for (Size t = 0; t < temp.size(); t++) dataset.data[joker, t] = cia[t];
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
std::ostream& operator<<(std::ostream& os, const CIARecord& /* cr */) {
  os << "CIARecord output operator not yet implemented." << '\n';
  return os;
}

void xml_io_stream<CIARecord>::write(std::ostream& os,
                                     const CIARecord& x,
                                     bofstream* pbofs,
                                     std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.TwoSpecies(), pbofs);
  xml_write_to_stream(os, x.Data(), pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<CIARecord>::read(std::istream& is,
                                    CIARecord& x,
                                    bifstream* pbifs) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.TwoSpecies(), pbifs);
  xml_read_from_stream(is, x.Data(), pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  throw std::runtime_error(
      std::format("Error reading {}:\n{}", type_name.data(), e.what()));
}

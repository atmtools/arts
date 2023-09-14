/**
 * @file   m_raw.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to generating y-data from raw data
 */


#include <arts_omp.h>
#include <artstime.h>
#include "raw.h"


void yColdAtmHot(Vector& y,
                 const Vector& cold,
                 const Vector& atm,
                 const Vector& hot,
                 const Numeric& cold_temp,
                 const Numeric& hot_temp,
                 const Index& calib)
{
  ARTS_USER_ERROR_IF (cold.size() not_eq atm.size() or atm.size() not_eq hot.size(),
                      "Length of vectors must be correct");
   
  y.resize(atm.size());
  if (calib) {
    y = Raw::Calibration::calibration(cold, atm, hot, cold_temp, hot_temp);
  } else {
    y = Raw::Calibration::systemtemp(cold, hot, cold_temp, hot_temp);
  }
}


void ybatchColdAtmHotAtmCycle(ArrayOfVector& ybatch,
                              ArrayOfTime& sensor_time,
                              const ArrayOfVector& level0_data,
                              const ArrayOfTime& level0_time,
                              const Vector& cold_temp,
                              const Vector& hot_temp,
                              const Index& first_c_index)
{
  ARTS_USER_ERROR_IF(level0_data.size() not_eq static_cast<Size>(cold_temp.size()) or
                     level0_data.size() not_eq static_cast<Size>(hot_temp.size()),
                     "Length of vectors must be correct");
  ARTS_USER_ERROR_IF (level0_time.size() not_eq level0_data.size() and
                      level0_time.size() not_eq 0,
                      "Bad level0_time length, must be empty of same as level0_data");
  
  ybatch = Raw::Calibration::caha(level0_data, cold_temp, hot_temp, first_c_index);
  
  // Fix time using the same method as CAHA
  if (level0_time.size()) {
    sensor_time.resize(ybatch.size());
    // First position as described by method is at H if this index is too large
    const Index pos = first_c_index - ((first_c_index > 1) ? 2 : 0);
    for (Size i=0; i<sensor_time.size(); i++) {
      sensor_time[i] = level0_time[pos + 2*i];
    }
  }
}


void ybatchTimeAveraging(ArrayOfVector& ybatch,
                         ArrayOfTime& sensor_time,
                         const String& time_step,
                         const Index& disregard_first,
                         const Index& disregard_last)
{
  // Size of problem
  const Size n=sensor_time.size();
  ARTS_USER_ERROR_IF (sensor_time.size() not_eq n,
                      "Time vector length must match input data length");
  
  // Time is not decreasing
  ARTS_USER_ERROR_IF (not std::is_sorted(sensor_time.cbegin(), sensor_time.cend()),
                      "Time vector cannot decrease");
  
  // Find the limits of the range
  auto lims = time_steps(sensor_time, time_stepper_selection(time_step));
  
  // Output variables
  ArrayOfVector ybatch_out;
  ArrayOfTime sensor_time_out;
  
  if (static_cast<Size>(lims.front()) == n) {
    ybatch_out.resize(0);
    sensor_time_out.resize(0);
  } else {
    
    // Frequency grids
    const Index k = ybatch[0].size();
    ARTS_USER_ERROR_IF (not std::all_of(ybatch.cbegin(), ybatch.cend(),
                                        [k](auto& x){return x.size() == k;}),
                        "Bad frequency grid size in input data; expects all equal");
    
    // Allocate output
    const Index m = lims.size() - 1;
    ARTS_USER_ERROR_IF (m < 0,
                        "Must include last if time step covers all of the range");
    ybatch_out = ArrayOfVector(m, Vector(k));
    sensor_time_out = ArrayOfTime(m);
    
    // Fill output
    #pragma omp parallel for if (not arts_omp_in_parallel()) schedule(guided)
    for (Index i=0; i<m; i++) {
      sensor_time_out[i] = mean_time(sensor_time, lims[i], lims[i+1]);
      Raw::Average::nanavg(ybatch_out[i], ybatch, lims[i], lims[i+1]);
    }
  }
  
  if (disregard_first) {
    ybatch_out.erase(ybatch_out.begin());
    sensor_time_out.erase(sensor_time_out.begin());
  }
  
  if (disregard_last) {
    ybatch_out.pop_back();
    sensor_time_out.pop_back();
  }
  
  ybatch = ybatch_out;
  sensor_time = sensor_time_out;
}


void ybatchTroposphericCorrectionNaiveMedianForward(ArrayOfVector& ybatch_corr,
                                                    ArrayOfVector& ybatch,
                                                    const ArrayOfIndex& range,
                                                    const Vector& trop_temp,
                                                    const Numeric& targ_temp)
{
  // Size of problem
  const Index n=ybatch.size();
  const Index m=n ? ybatch[0].size() : 0;
  
  ARTS_USER_ERROR_IF (m == 0, "A frequency range is required")
  
  ARTS_USER_ERROR_IF (std::any_of(ybatch.begin(), ybatch.end(),
                                  [m](auto& y){return y.size() not_eq m;}),
                      "Bad input size, all of ybatch must match itself");
  
  ARTS_USER_ERROR_IF (trop_temp.size() not_eq n and trop_temp.size() not_eq 1,
                      "Bad input size, trop_temp must match ybatch or be 1-long,\n"
                      "trop_temp: [", trop_temp, "]\ntrop_temp.size(): ",
                      trop_temp.size(), "\nybatch.size(): ", n);
  
  // This algorithm stores partial transmission and median and tropospheric temperature in the correction terms
  ybatch_corr = ArrayOfVector(n, Vector(3));
  
  // Compute tropospheric correction
  for (Index i=0; i<n; i++) {
    ybatch_corr[i][2] = trop_temp.size() > 1 ? trop_temp[i] : trop_temp[0];
    ybatch_corr[i][0] = Raw::Average::nanmedian(ybatch[i], range);
    ybatch_corr[i][1] = std::exp(- std::log((ybatch_corr[i][2] - ybatch_corr[i][0])  / (ybatch_corr[i][2] - targ_temp)));
  }
  
  // Apply correction
  for (Index i=0; i<n; i++) {
    ybatch[i] *= ybatch_corr[i][1];
    ybatch[i] += ybatch_corr[i][2] * (1 - ybatch_corr[i][1]);
  }
}


void ybatchTroposphericCorrectionNaiveMedianInverse(ArrayOfVector& ybatch,
                                                    const ArrayOfVector& ybatch_corr)
{
  // Size of problem
  const Size n=ybatch.size();
  ARTS_USER_ERROR_IF ((std::any_of(ybatch_corr.begin(), ybatch_corr.end(),
                                   [](auto& corr){return corr.size() not_eq 3;})) or ybatch_corr.size() not_eq n,
                      "Bad input size, all of ybatch_corr must match ybatch and have three elements each");
  
  // Apply inverse of correction
  for (Size i=0; i<n; i++) {
    ybatch[i] -= ybatch_corr[i][2] * (1 - ybatch_corr[i][1]);
    ybatch[i] /= ybatch_corr[i][1];
  }
}

void yMaskOutsideMedianRange(Vector& y,
                             const Numeric& dx) {
  const Numeric median = Raw::Average::nanmedian(y);
  auto mask = Raw::Mask::out_of_bounds(y, median-dx, median+dx);
  
  for (Index i=0; i<y.size(); i++) {
    y[i] = mask[i] ? std::numeric_limits<Numeric>::quiet_NaN() : y[i];
  }
}

void ybatchMaskOutsideMedianRange(ArrayOfVector& ybatch,
                                  const Numeric& dx) {
  for (Vector& y: ybatch) yMaskOutsideMedianRange(y, dx);
}

void yDoublingMeanFocus(Vector& f_grid,
                        Vector& y,
                        const Numeric& f0,
                        const Numeric& df) {
  if (f_grid.size() not_eq y.size()) {
    throw std::runtime_error("f_grid and y must have the same size");
  }
  
  if (not is_increasing(f_grid)) {
    throw std::runtime_error("f_grid must be sorted and ever increasing");
  }
  
  if (f_grid.size() < 2) {
    throw std::runtime_error("Must have at least 2 frequency grid points");
  }
  
  // Sets defaults by method description
  const Numeric F0 = f0 > 0 ? f0 : mean(f_grid);
  const Numeric DF = df > 0 ? df : 10 * (f_grid[1] - f_grid[0]);
  
  if (f_grid[0] <= F0 or F0 >= last(f_grid)) {
    throw std::runtime_error("Needs F0 in the range of f_grid");
  }
  
  auto red = Raw::Reduce::focus_doublescaling(f_grid, F0, DF);
  y = Raw::Reduce::nanfocus(y, red);
  f_grid = Raw::Reduce::nanfocus(f_grid, red);
}

void ybatchDoublingMeanFocus(Vector& f_grid,
                             ArrayOfVector& ybatch,
                             const Numeric& f0,
                             const Numeric& df) {
  for (Vector& y: ybatch) {
    if (f_grid.size() not_eq y.size()) {
      throw std::runtime_error("f_grid and all of ybatch must have the same size");
    }
  }
  
  if (not is_increasing(f_grid)) {
    throw std::runtime_error("f_grid must be sorted and ever increasing");
  }
  
  if (f_grid.size() < 2) {
    throw std::runtime_error("Must have at least 2 frequency grid points");
  }
  
  // Sets defaults by method description
  const Numeric F0 = f0 > 0 ? f0 : mean(f_grid);
  const Numeric DF = df > 0 ? df : 10 * (f_grid[1] - f_grid[0]);
  
  if (last(f_grid) <= F0 or f_grid[0] >= F0) {
    throw std::runtime_error("Needs F0 in the range of f_grid");
  }
  
  auto red = Raw::Reduce::focus_doublescaling(f_grid, F0, DF);
  for (Vector& y: ybatch) y = Raw::Reduce::nanfocus(y, red);
  f_grid = Raw::Reduce::nanfocus(f_grid, red);
}

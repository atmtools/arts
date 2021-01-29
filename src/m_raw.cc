/* Copyright (C) 2020
 * Richard Larsson <larsson@mps.mpg.de>
 * 
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 2, or (at your option) any
 * later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307,
 * USA. */

/**
 * @file   m_raw.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to generating y-data from raw data
 */


#include "artstime.h"
#include "raw.h"


void yColdAtmHot(Vector& y,
                 const Vector& cold,
                 const Vector& atm,
                 const Vector& hot,
                 const Numeric& cold_temp,
                 const Numeric& hot_temp,
                 const Index& calib,
                 const Verbosity&)
{
  ARTS_USER_ERROR_IF (cold.nelem() not_eq atm.nelem() or atm.nelem() not_eq hot.nelem(),
                      "Length of vectors must be correct");
   
  y.resize(atm.nelem());
  if (calib) {
    y = Raw::Calibration::calibration(cold, atm, hot, cold_temp, hot_temp);
  } else {
    y = Raw::Calibration::systemtemp(cold, hot, cold_temp, hot_temp);
  }
}


void ybatchCAHA(ArrayOfVector& ybatch,
                ArrayOfTime& time_grid,
                const ArrayOfVector& rawdata,
                const ArrayOfTime& timedata,
                const Vector& cold_temp,
                const Vector& hot_temp,
                const Index& c_offset,
                const Verbosity&)
{
  ARTS_USER_ERROR_IF(rawdata.nelem() not_eq cold_temp.nelem() or
                     rawdata.nelem() not_eq hot_temp.nelem(),
                     "Length of vectors must be correct");
  ARTS_USER_ERROR_IF (timedata.nelem() not_eq rawdata.nelem() and timedata.nelem() not_eq 0,
                      "Bad timedata length, must be empty of same as data");
  
  ybatch = Raw::Calibration::caha(rawdata, cold_temp, hot_temp, c_offset);
  
  // Fix time using the same method as CAHA
  if (timedata.nelem()) {
    time_grid.resize(ybatch.nelem());
    const Index pos = c_offset - ((c_offset > 1) ? 2 : 0);
    for (Index i=0; i<time_grid.nelem(); i++) {
      time_grid[i] = timedata[pos + 2*i];
    }
  }
}


void ybatchTimeAveraging(ArrayOfVector& ybatch,
                         ArrayOfTime& time_grid,
                         const String& time_step,
                         const Index& disregard_first,
                         const Index& disregard_last,
                         const Verbosity&)
{
  // Size of problem
  const Index n=time_grid.nelem();
  ARTS_USER_ERROR_IF (time_grid.nelem() not_eq n,
                      "Time vector length must match input data length");
  
  // Time is not decreasing
  ARTS_USER_ERROR_IF (not std::is_sorted(time_grid.cbegin(), time_grid.cend()),
                      "Time vector cannot decrease");
  
  // Find the limits of the range
  auto lims = time_steps(time_grid, time_stepper_selection(time_step));
  
  // Output variables
  ArrayOfVector ybatch_out;
  ArrayOfTime time_grid_out;
  
  if (lims.front() == n) {
    ybatch_out.resize(0);
    time_grid_out.resize(0);
  } else {
    
    // Frequency grids
    const Index k = ybatch[0].nelem();
    ARTS_USER_ERROR_IF (not std::all_of(ybatch.cbegin(), ybatch.cend(),
                                        [k](auto& x){return x.nelem() == k;}),
                        "Bad frequency grid size in input data; expects all equal");
    
    // Allocate output
    const Index m = lims.nelem() - 1;
    ARTS_USER_ERROR_IF (m < 0,
                        "Must include last if time step covers all of the range");
    ybatch_out = ArrayOfVector(m, Vector(k));
    time_grid_out = ArrayOfTime(m);
    
    // Fill output
    #pragma omp parallel for if (not arts_omp_in_parallel()) schedule(guided)
    for (Index i=0; i<m; i++) {
      time_grid_out[i] = mean_time(time_grid, lims[i], lims[i+1]);
      Raw::Average::nanavg(ybatch_out[i], ybatch, lims[i], lims[i+1]);
    }
  }
  
  if (disregard_first) {
    ybatch_out.erase(ybatch_out.begin());
    time_grid_out.erase(time_grid_out.begin());
  }
  
  if (disregard_last) {
    ybatch_out.pop_back();
    time_grid_out.pop_back();
  }
  
  ybatch = ybatch_out;
  time_grid = time_grid_out;
}


void ybatchTroposphericCorrectionNaiveMedianForward(ArrayOfVector& ybatch_corr,
                                                    ArrayOfVector& ybatch,
                                                    const ArrayOfIndex& range,
                                                    const Vector& trop_temp,
                                                    const Numeric& targ_temp,
                                                    const Verbosity&)
{
  // Size of problem
  const Index n=ybatch.nelem();
  
  const Index m=n?ybatch[0].nelem():0;
  
  ARTS_USER_ERROR_IF (std::any_of(ybatch.begin(), ybatch.end(),
                                  [m](auto& y){return y.nelem() not_eq m;}),
                      "Bad input size, all of ybatch must match itself");
  ARTS_USER_ERROR_IF (trop_temp.nelem() not_eq n,
                      "Bad input size, trop_temp must match ybatch");
  
  // This algorithm stores partial transmission and median and tropospheric temperature in the correction terms
  ybatch_corr = ArrayOfVector(n, Vector(3));
  
  // Compute tropospheric correction
  for (Index i=0; i<n; i++) {
    ybatch_corr[i][2] = trop_temp.nelem() > 1 ? trop_temp[i] : trop_temp[0];
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
                                                    const ArrayOfVector& ybatch_corr,
                                                    const Verbosity&)
{
  // Size of problem
  const Index n=ybatch.nelem();
  ARTS_USER_ERROR_IF ((std::any_of(ybatch_corr.begin(), ybatch_corr.end(),
                                   [](auto& corr){return corr.nelem() not_eq 3;})) or ybatch_corr.nelem() not_eq n,
                      "Bad input size, all of ybatch_corr must match ybatch and have three elements each");
  
  // Apply inverse of correction
  for (Index i=0; i<n; i++) {
    ybatch[i] -= ybatch_corr[i][2] * (1 - ybatch_corr[i][1]);
    ybatch[i] /= ybatch_corr[i][1];
  }
}

void yMaskOutsideMedianRange(Vector& y,
                             const Numeric& dx,
                             const Verbosity&) {
  const Numeric median = Raw::Average::nanmedian(y);
  auto mask = Raw::Mask::out_of_bounds(y, median-dx, median+dx);
  
  for (Index i=0; i<y.nelem(); i++) {
    y[i] = mask[i] ? std::numeric_limits<Numeric>::quiet_NaN() : y[i];
  }
}

void ybatchMaskOutsideMedianRange(ArrayOfVector& ybatch,
                                  const Numeric& dx,
                                  const Verbosity& verbosity) {
  for (Vector& y: ybatch) yMaskOutsideMedianRange(y, dx, verbosity);
}

void yDoublingMeanFocus(Vector& f_grid,
                        Vector& y,
                        const Numeric& f0,
                        const Numeric& df,
                        const Verbosity&) {
  if (f_grid.nelem() not_eq y.nelem()) {
    throw std::runtime_error("f_grid and y must have the same size");
  }
  
  if (not is_increasing(f_grid)) {
    throw std::runtime_error("f_grid must be sorted and ever increasing");
  }
  
  if (f_grid.nelem() < 2) {
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
                             const Numeric& df,
                             const Verbosity&) {
  for (Vector& y: ybatch) {
    if (f_grid.nelem() not_eq y.nelem()) {
      throw std::runtime_error("f_grid and all of ybatch must have the same size");
    }
  }
  
  if (not is_increasing(f_grid)) {
    throw std::runtime_error("f_grid must be sorted and ever increasing");
  }
  
  if (f_grid.nelem() < 2) {
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

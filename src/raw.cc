/* Copyright (C) 2020
 * Richard Larsson <ric.larsson@gmail.com>
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
 * @file   raw.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to generating y-data from raw data
*/

#include <algorithm>

#include "constants.h"
#include "raw.h"

template <typename T>
constexpr bool isnormal_or_zero(T x) noexcept {return std::isnormal(x) or x == 0;}

namespace Raw {
namespace Calibration {
Vector calibration(const Vector& pc, const Vector& pa, const Vector& ph, Numeric tc, Numeric th) noexcept
{
  const Index N = Index(pc.size());
  Vector calib(N);
  for (Index i=0; i<N; i++) calib[i] = calibration(pc[i], pa[i], ph[i], tc, th);
  return calib;
}

Vector systemtemp(const Vector& pc, const Vector& ph, Numeric tc, Numeric th) noexcept
{
  const Index N = Index(pc.size());
  Vector systemp(N);
  for (Index i=0; i<N; i++) systemp[i] = systemtemp(Numeric(pc[i]), Numeric(ph[i]), tc, th);
  return systemp;
}

ArrayOfVector caha(const ArrayOfVector& data, const Vector& tcvec, const Vector& thvec, const Index first_c_index)
{
  assert(first_c_index >= 0);
  
  /* Must compute the positions due to c-Index offset
   * 
   * The start position should be on C or H
   * if first_c_index > 1 then start at H two steps early
   * 
   * The cold position is always a multiple of 4 from
   * first_c_index.  In case first_c_index is > 3, we
   * must take the modulus of the position
   * 
   * The hot position is either two indices before or
   * two indices behind the cold position, whichever
   * is lowest and positive and below 4 */
  const Index start = first_c_index - ((first_c_index > 1) ? 2 : 0);
  const Index cpos = first_c_index % 4;
  const Index hpos = cpos + ((cpos < 2) ? 2 : -2);
  
  // Initialize pointers and values of interest
  const Vector * pc = nullptr;
  const Vector * pa = nullptr;
  const Vector * ph = nullptr;
  Numeric tc=0, th=0;
  
  // Output vector
  ArrayOfVector out;
  
  // For all data entries
  for (Index i=start; i<data.nelem(); i++) {
    // Cycle is CAHA → 4 steps
    const Index pos = i % 4;
    
    // Set the values when we are at the correct position
    if (pos == cpos) {
      pc = &data[i];
      tc = tcvec[i];
    } else if (pos == hpos) {
      ph = &data[i];
      th = thvec[i];
    } else {
      pa = &data[i];
    }
    
    /* Add new measurement if we are at C or H since we have
     *  then completed a new calibration cycle 
     * 
     * Note that nullptr is false */
    if ((pos == cpos or pos == hpos) and pc and pa and ph) {
      out.emplace_back(calibration(*pc, *pa, *ph, tc, th));
    }
  }
  
  return out;
}
}  // Calibration

namespace Average {
VectorView avg(VectorView y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{ 
  // True end
  const Index end = end_tmp >= 0 ? end_tmp : 1 + ys.nelem() + end_tmp;
  
  // Scale
  const Numeric scale = 1 / Numeric(end - start);
  
  // Reset
  y = 0;
  
  // Compute
  for (Index k=start; k<end; k++)
    for (Index i=0; i<y.nelem(); i++)
      y[i] += ys[k][i] * scale;
  return y;
}

VectorView nanavg(VectorView y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{ 
  // True end
  const Index end = end_tmp >= 0 ? end_tmp : 1 + ys.nelem() + end_tmp;
  
  // Reset y
  y = 0;
  
  // Compute the averages ignoring all NaN
  for (Index i=0; i<y.nelem(); i++) {
    Index numnormal=0;
    for (Index k=start; k<end; k++) {
      if (isnormal_or_zero(ys[k][i])) {
        y[i] += ys[k][i];
        numnormal++;
      }
    }
    
    // Compute the average or set to NaN
    if (numnormal) {
      y[i] /= Numeric(numnormal);
    } else {
      y[i] = std::numeric_limits<Numeric>::quiet_NaN();
    }
  }
  
  return y;
}

VectorView var(VectorView var, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{
  // True end
  const Index end = end_tmp >= 0 ? end_tmp : 1 + ys.nelem() + end_tmp;
  
  // Scale
  const Numeric scale = 1 / Numeric(end - start);
  
  // Reset
  var = 0;
  
  // Compute
  for (Index k=start; k<end; k++)
    for (Index i=0; i<y.nelem(); i++)
      var[i] += Constant::pow2(ys[k][i] - y[i]) * scale;
  return var;
}

VectorView nanvar(VectorView var, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{
  // True end
  const Index end = end_tmp >= 0 ? end_tmp : 1 + ys.nelem() + end_tmp;
  
  // Reset
  var = 0;
  
  // Compute
  for (Index i=0; i<y.nelem(); i++) {
    Index numnormal=0;
    for (Index k=start; k<end; k++) {
      if (isnormal_or_zero(ys[k][i])) {
        var[i] += Constant::pow2(ys[k][i] - y[i]);
        numnormal++;
      }
    }
    
    // Compute the average or set to NaN
    if (numnormal > 0) {
      var[i] /= Numeric(numnormal);
    } else {
      var[i] = std::numeric_limits<Numeric>::quiet_NaN();
    }
  }
  return var;
}

VectorView std(VectorView std, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{
  var(std, y, ys, start, end_tmp);
  std::transform(std.begin(), std.end(), std.begin(), [](const auto& x){return std::sqrt(x);});
  return std;
}

VectorView nanstd(VectorView std, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{
  nanvar(std, y, ys, start, end_tmp);
  std::transform(std.begin(), std.end(), std.begin(), [](const auto& x){return std::sqrt(x);});
  return std;
}

MatrixView cov(MatrixView cov, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{
  // True end
  const Index end = end_tmp >= 0 ? end_tmp : 1 + ys.nelem() + end_tmp;
  
  // Scale
  const Numeric scale = 1 / Numeric(end - start - 1);
  
  // Reset
  cov = 0;
  
  // Compute
  for (Index k=start; k<end; k++)
    for (Index i=0; i<y.nelem(); i++)
      for (Index j=0; j<y.nelem(); j++)
        cov(i, j) += (ys[k][i] - y[i]) * (ys[k][j] - y[j]) * scale;
  return cov;
}

Numeric median(const ConstVectorView v, const ArrayOfIndex& pos)
{
  // Size of problem
  const Index n = pos.nelem() ? pos.nelem() : v.nelem();
  
  // Create a vector to sort
  ArrayOfNumeric calc(n);
  for (Index i=0; i<n; i++)
    if (pos.nelem())
      calc[i] = v[pos[i]];
    else
      calc[i] = v[i];
  
  // Sort the vector
  std::sort(calc.begin(), calc.end());
  
  // Return the median
  if (n % 2)
    return calc[n/2];
  else
    return (calc[(n-1)/2] + calc[n/2]) / 2;
}

Numeric nanmedian(const ConstVectorView v, const ArrayOfIndex& pos)
{
  // Size of problem
  const Index n = pos.nelem() ? pos.nelem() : v.nelem();
  
  // Create a vector to sort
  ArrayOfNumeric calc(n);
  for (Index i=0; i<n; i++) {
    if (pos.nelem()) {
      calc[i] = v[pos[i]];
    } else {
      calc[i] = v[i];
    }
  }
    
  // Sort the vector
  std::sort(calc.begin(), calc.end());
  
  // Remove non-normal
  auto newend = std::remove_if(calc.begin(), calc.end(), [](auto x){return not isnormal_or_zero(x);});
  const Index numnormal = n - Index(calc.end() - newend);
  
  // Return the median
  if (numnormal) {
    if (numnormal % 2) {
      return calc[numnormal/2];
    } else {
      return (calc[(numnormal-1)/2] + calc[numnormal/2]) / 2;
    }
  } else {
    return std::numeric_limits<Numeric>::quiet_NaN();
  }
}
}  // Average

namespace Reduce {
ArrayOfIndex focus_doublescaling(const Vector& x, const Numeric x0, const Numeric dx)
{
  ArrayOfIndex scale(x.size());
  
  // Scaling function, rounds down an Index
  auto scaling = [](Numeric X, Numeric X0, Numeric DX) {
    const Index p = 1 + Index(std::abs(X - X0) / DX);  // floor + 1 is close enough
    Index this_scale = 0;
    while ((p >> this_scale) > 1) this_scale++;
    return 1 << this_scale;
  };
  
  for (Index i=0; i<x.size(); i++) scale[i] = scaling(x[i], x0, dx);
  return scale;
}

std::pair<Index, Index> find_first_and_last_1(const ArrayOfIndex& x)
{
  Index first=x.nelem() - 1;
  Index last=0;
  for (Index i=0; i<x.nelem(); i++) {
    if (x[i] == 1) {
      last = std::max(i, last);
      first = std::min(i, first);
    }
  }
  
  assert(first <= last);
  return {first, last};
}

Vector focus(const Vector& x, const ArrayOfIndex& scaling)
{
  const Index N = x.size();
  Index first_i, last_i;
  
  // Find the first and last one in scaling
  const auto [firstone, lastone] = find_first_and_last_1(scaling);
  
  // All the central numbers
  std::vector<Numeric> out;
  for (Index i=firstone; i<=lastone; i++) out.emplace_back(x[i]);
  
  // All the last numbers appended
  first_i = lastone + 1;
  while (first_i not_eq N) {
    last_i = first_i + scaling[first_i];
    if (last_i > N) last_i = N;
    out.emplace_back(mean(x[Range(first_i, last_i-first_i)]));
    first_i = last_i;
  }
  
  // All the first numbers prepended
  std::vector<Numeric> fx(0);
  last_i = firstone;
  while (last_i > 0) {
    first_i = last_i - scaling[last_i];
    if (first_i < 0) first_i = 0;
    out.insert(out.begin(), mean(x[Range(first_i, last_i-first_i)]));
    last_i = first_i;
  }
  
  return out;
}

Vector nanfocus(const Vector& x, const ArrayOfIndex& scaling)
{
  const Index N = x.size();
  Index first_i, last_i;
  
  // Find the first and last one in scaling
  const auto [firstone, lastone] = find_first_and_last_1(scaling);
  
  // All the central numbers
  std::vector<Numeric> out;
  for (Index i=firstone; i<=lastone; i++) {
    if (isnormal_or_zero(x[i])) {
      out.emplace_back(x[i]);
    } else {
      out.emplace_back(std::numeric_limits<Numeric>::quiet_NaN());
    }
  }
  
  // All the last numbers appended
  first_i = lastone + 1;
  while (first_i not_eq N) {
    last_i = first_i + scaling[first_i];
    if (last_i > N) last_i = N;
    out.emplace_back(nanmean(x[Range(first_i, last_i-first_i)]));
    first_i = last_i;
  }
  
  // All the first numbers prepended
  std::vector<Numeric> fx(0);
  last_i = firstone;
  while (last_i > 0) {
    first_i = last_i - scaling[last_i];
    if (first_i < 0) first_i = 0;
    out.insert(out.begin(), nanmean(x[Range(first_i, last_i-first_i)]));
    last_i = first_i;
  }
  
  return out;
}
}  // Reduce

namespace Mask {
std::vector<bool> out_of_bounds(const Vector& x, const Numeric xmin, const Numeric xmax)
{
  std::vector<bool> mask(x.size(), false);
  for (Index i=0; i<x.size(); i++) mask[i] = (xmin > x[i]) or (xmax < x[i]) or (not isnormal_or_zero(x[i]));
  return mask;
}

VectorView mask(VectorView x, const std::vector<bool>& masking)
{
  const Index N = x.nelem();
  for (Index i=0; i<N; i++)
    if (masking[i])
      x[i] = std::numeric_limits<Numeric>::quiet_NaN();
  return x;
}
}  // Mask

namespace Correction {
  Numeric naive_tropospheric_singletau_median(const ConstVectorView bt, const Numeric trop_bt, const Numeric target_bt) {
    return - std::log((trop_bt - Average::nanmedian(bt)) / (trop_bt - target_bt));
}

VectorView naive_tropospheric(VectorView bt, const Numeric tau, const Numeric trop_bt) {
  const Numeric t = std::exp(-tau);
  if (std::isnormal(t)) {
    bt /= t;
    bt += trop_bt * std::expm1(-tau) / t;
  }
  return bt;
}
}  // Correction
}  // Raw

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
 * @file   raw.cc
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to generating y-data from raw data
*/

#include <algorithm>

#include "constants.h"
#include "raw.h"

void linalg::avg(VectorView y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
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
}

void linalg::var(VectorView var, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
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
}

void linalg::std(VectorView std, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
{
  var(std, y, ys, start, end_tmp);
  std::transform(std.begin(), std.end(), std.begin(), [](auto& x){return std::sqrt(x);});
}

void linalg::cov(MatrixView cov, const Vector& y, const ArrayOfVector& ys, const Index start, const Index end_tmp)
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
}

Numeric linalg::median(const ConstVectorView v, const ArrayOfIndex& pos)
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

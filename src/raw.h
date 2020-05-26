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
 * @file   raw.h
 * @author Richard Larsson
 * @date   2020-04-13
 * 
 * @brief  Stuff related to generating y-data from raw data
*/

#ifndef RAW_H
#define RAW_H

#include "covariance_matrix.h"
#include "matpackI.h"

/** Computes the linear calibration formula
 * 
 * \f[ T_a = T_c + \frac{(T_h - T_c)(P_a - P_c)}{P_h - P_c} \f]
 * 
 * @param[in] pc \f$ P_c \f$, linear power of cold load; cold count
 * @param[in] pa \f$ P_a \f$, linear power of atmospheric measurement; atmospheric count
 * @param[in] ph \f$ P_h \f$, linear power of got load; hot count
 * @param[in] tc \f$ T_c \f$, cold load temperature
 * @param[in] th \f$ T_h \f$, hot load temperature
 * @return \f$ T_a \f$, atmospheric temperature
 */
constexpr Numeric calibration(Numeric pc, Numeric pa, Numeric ph, Numeric tc, Numeric th) noexcept
{
  return tc + (th - tc) * (pa - pc) / (ph - pc);
}

/** Computes the linear receiver temperature formula
 * 
 * \f[ T_r = \frac{T_hP_c - T_cP_h}{P_h - P_c} \f]
 * 
 * @param[in] pc \f$ P_c \f$, linear power of cold load; cold count
 * @param[in] ph \f$ P_h \f$, linear power of got load; hot count
 * @param[in] tc \f$ T_c \f$, cold load temperature
 * @param[in] th \f$ T_h \f$, hot load temperature
 * @return \f$ T_r \f$, receiver temperature
 */
constexpr Numeric systemtemp(Numeric pc, Numeric ph, Numeric tc, Numeric th) noexcept
{
  return (th*pc - tc*ph) / (ph - pc);
}

namespace linalg {
/** Compute the average of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] y N-dimensional Vector; Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 */
void avg(VectorView y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the standard deviation of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] std N-dimensional Vector; Standard deviations
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 */
void std(VectorView std, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the variance of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] var N-dimensional Vector; Variance
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 */
void var(VectorView var, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the covariance matrix of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] cov N-by-N-dimensional Matrix; Covariance matrix
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 */
void cov(MatrixView cov, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Get the median of the vector in the range
 * 
 * If pos is empty, use all numbers
 * 
 * @param[in] v A Vector of any length
 * @param[in] pos Positions to use in the vector
 */
Numeric median(const ConstVectorView v, const ArrayOfIndex& pos=ArrayOfIndex{});
};  // linalg

#endif  // RAW_H

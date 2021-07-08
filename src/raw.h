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

namespace Raw {
namespace Calibration {
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
constexpr Numeric calibration(Numeric pc, Numeric pa, Numeric ph, Numeric tc, Numeric th) noexcept {
  return tc + (th - tc) * (pa - pc) / (ph - pc);
}

/** Computes the linear calibration formula elementwise
 * 
 * \f[ \vec{T}_a = T_c + \frac{(T_h - T_c)(\vec{P}_a - \vec{P}_c)}{\vec{P}_h - \vec{P}_c} \f]
 * 
 * @param[in] pc \f$ \vec{P}_c \f$, linear power of cold load; cold count
 * @param[in] pa \f$ \vec{P}_a \f$, linear power of atmospheric measurement; atmospheric count
 * @param[in] ph \f$ \vec{P}_h \f$, linear power of got load; hot count
 * @param[in] tc \f$ T_c \f$, cold load temperature
 * @param[in] th \f$ T_h \f$, hot load temperature
 * @return \f$ \vec{T}_a \f$, atmospheric temperature
 */
Vector calibration(const Vector& pc, const Vector& pa, const Vector& ph, Numeric tc, Numeric th) noexcept;

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
constexpr Numeric systemtemp(Numeric pc, Numeric ph, Numeric tc, Numeric th) noexcept {
  return (th*pc - tc*ph) / (ph - pc);
}

/** Computes the linear receiver temperature formula elementwise
 * 
 * \f[ \vec{T}_r = \frac{T_h \vec{P}_c - T_c \vec{P}_h}{\vec{P}_h - \vec{P}_c} \f]
 * 
 * @param[in] pc \f$ \vec{P}_c \f$, linear power of cold load; cold count
 * @param[in] ph \f$ \vec{P}_h \f$, linear power of got load; hot count
 * @param[in] tc \f$ T_c \f$, cold load temperature
 * @param[in] th \f$ T_h \f$, hot load temperature
 * @return \f$ \vec{T}_r \f$, receiver temperature
 */
Vector systemtemp(const Vector& pc, const Vector& ph, Numeric tc, Numeric th) noexcept;

/** Calibrate the data by CAHA
 * 
 * CAHA stands for cold-atm-hot-atm, and is the presumed
 * observation scheme for the calibration.  The calibration
 * will happen on the atmospheric measurement between any two
 * pairs of cold and hot (cold then hot; hot then cold)
 * 
 * @param[in] data A list or array or vector of measurements
 * @param[in] tcvec A Vector of cold load temperatures at corresponding data entry
 * @param[in] thvec A Vector of hot load temperatures at corresponding data entry
 * @param[in] first_c_index The position of the first C in data [must be positive]
 * @return Calibrated values [use calibration() if calib, else use systemtemp()]
 */
ArrayOfVector caha(const ArrayOfVector& data, const Vector& tcvec, const Vector& thvec, const Index first_c_index);
}  // Calibration

namespace Average {
/** Compute the average of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] y N-dimensional Vector; Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return y
 */
VectorView avg(VectorView y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the average of the ranged ys ignoring all non-normal values
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * If this range contains only NaN, the value is set to NaN
 * 
 * @param[in,out] y N-dimensional Vector; Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return y
 */
VectorView nanavg(VectorView y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the standard deviation of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] std N-dimensional Vector; Standard deviations
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return std
 */
VectorView std(VectorView std, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the standard deviation of the ranged ys ignoring all non-normal values
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * If this range contains only NaN, the value is set to NaN. Same if y contains NaN
 * 
 * @param[in,out] std N-dimensional Vector; Standard deviations
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return std
 */
VectorView nanstd(VectorView std, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the variance of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] var N-dimensional Vector; Variance
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return var
 */
VectorView var(VectorView var, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the variance of the ranged ys ignoring all non-normal values
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * If this range contains only NaN, the value is set to NaN. Same if y contains NaN
 * 
 * @param[in,out] var N-dimensional Vector; Variance
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return var
 */
VectorView nanvar(VectorView var, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Compute the covariance matrix of the ranged ys
 * 
 * The range of ys is [start, end) if end is positive or [start, ys.nelem()+end] otherwise
 * 
 * @param[in,out] cov N-by-N-dimensional Matrix; Covariance matrix
 * @param[in] y Averages of y in given range
 * @param[in] ys list of N-dimensional Vector(s) to average
 * @param[in] start First index in ys
 * @param[in] end Past last index representation
 * @return cov
 */
MatrixView cov(MatrixView cov, const Vector& y, const ArrayOfVector& ys, const Index start=0, const Index end=-1);

/** Get the median of the vector in the range
 * 
 * If pos is empty, use all numbers
 * 
 * @param[in] v A Vector of any length
 * @param[in] pos Positions to use in the vector
 * @return The median
 */
Numeric median(const ConstVectorView v, const ArrayOfIndex& pos=ArrayOfIndex{});

/** Get the median of the vector in the range ignoring all non-normal values
 * 
 * If pos is empty, use all numbers
 * 
 * @param[in] v A Vector of any length
 * @param[in] pos Positions to use in the vector
 * @return The median ignoring all non-normal values (or NaN)
 */
Numeric nanmedian(const ConstVectorView v, const ArrayOfIndex& pos=ArrayOfIndex{});
};  // Average

namespace Reduce {
/** Returns the relative position scale for each value in x
 * 
 * This is meant to be used with the focus routine to retain mostly values
 * close to x0, as a reduction mechanism
 * 
 * This particular routine will compute the focus scaling by doubling the points
 * so that for (positive values only, mirroring applies at absolute values):
 * 
 * x0 + 0 * dx to x0 + 1 * dx will have 1-long steps
 * x0 + 1 * dx to x0 + 2 * dx will have 2-long steps
 * x0 + 2 * dx to x0 + 4 * dx will have 4-long steps
 * x0 + 4 * dx to x0 + 8 * dx will have 8-long steps
 * 
 * This doubling continues to the end (on both side around x0).  The focus routine
 * might still keep more noisy data at the end of the output focus vector when there
 * are too few data points to have exactly a doubling of the numbers of steps
 * 
 * Note: the rounding is simplistic so that x0 + N * dx belongs to the 'next' step
 * 
 * @param[in] x The grid of x-values [strictly increasing by the same amount]
 * @param[in] x0 A central position so that x[0] < x0 < x[x.size() - 1]
 * @param[in] dx A stepsize so that can be used to compute a relative step
 * @return Scaling of some positions (same size at x)
 */
ArrayOfIndex focus_doublescaling(const Vector& x, const Numeric x0, const Numeric dx);

/** Returns the re-averaged values of x according to scaling
 * 
 * This refocus the input Vector x to reduce its size
 * 
 * The scaling should be such that it contains how many points you want the
 * rescaled position to contain at most.  Note that if there is not enough
 * points left, all the remaining points are returned averaged to what they
 * can be averaged to.
 * 
 * Note that this cannot defocus so there must be atleast 1 1 in scaling
 * representing retained focus.  Also note that all 1s must be continuous
 * in scaling (if they are all non-1s inbetween are ignored).  Lastly note
 * that 0s are not allowed in scaling (if there are, the function will loop
 * forever).
 * 
 * @param[in] x A vector to refocus
 * @param[in] scaling The scaling
 * @return A rescaled and focused x
 */
Vector focus(const Vector& x, const ArrayOfIndex& scaling);

/** Returns the re-averaged values of x according to scaling ignoring all non-normal values
 * 
 * This refocus the input Vector x to reduce its size
 * 
 * The scaling should be such that it contains how many points you want the
 * rescaled position to contain at most.  Note that if there is not enough
 * points left, all the remaining points are returned averaged to what they
 * can be averaged to.
 * 
 * Note that this cannot defocus so there must be atleast 1 1 in scaling
 * representing retained focus.  Also note that all 1s must be continuous
 * in scaling (if they are all non-1s inbetween are ignored).  Lastly note
 * that 0s are not allowed in scaling (if there are, the function will loop
 * forever).
 * 
 * If there is only NaN in a range then that value is NaN in the output
 * (retaining original sizes)
 * 
 * @param[in,out] x A vector to refocus
 * @param[in] scaling The scaling
 * @return A rescaled and focused x
 */
Vector nanfocus(const Vector& x, const ArrayOfIndex& scaling);
} // Reduce

namespace Mask {
/** Masks values that are out of bounds in x
 * 
 * Acceptable values are xmin <= x[i] <= xmax
 * 
 * All other i are marked as masked (true)
 * 
 * @param[in] x A vector to mask
 * @param[in] xmin The minimum acceptable value of x
 * @param[in] xmax The maximum acceptable value of x
 * @return A mask (true is masked, false is unmasked)
 */
std::vector<bool> out_of_bounds(const Vector& x, const Numeric xmin, const Numeric xmax);

/** Masks all true entries of masking in x by turning them into NaN
 * 
 * @param[in,out] x A Vector to mask
 * @param[in] masking The mask
 * @return x masked (note, x is masked in-place anyways)
 */
VectorView mask(VectorView x, const std::vector<bool>& masking);
}  // Mask

namespace Correction {
/** Naive tropospheric correction parameterization
 * 
 * Tropospheric correction comes from the idea seen in Raw::Correction::naive_tropospheric
 * 
 * The "naive" part of this implementation comes from assuming
 * 
 * \f[ \overline{\tau(f)} := - \ln\left(\frac{\overline{T_{b, trop}(f)} - \overline{T_{b, 0}(f)}}{\overline{T_{b, trop}(f)} - \overline{T_{b, target}(f)}} \right), \f]
 * 
 * where \f$ \overline{T_{b, 0}(f)} \f$ is the nan-median of input brightness temperature at the surface,
 * \f$ \overline{T_{b, trop}(f)} \f$ is infinite brightness temperature of the troposphere, and
 * \f$ \overline{T_{b, target}(f)} \f$ is brightness temperature of the tropopause.
 * Here \f$ f \f$ is the frequency dimension of bt, from which only the nan-median is extracted
 *
 * @param[in] bt The measurement vector \f$ T_{b, 0}(f) \f$
 * @param[in] trop_bt \f$ \overline{T_{b, trop}(f)} \f$
 * @param[in] target_bt \f$ \overline{T_{b, target}(f)} \f$
 * @return \f$ \overline{\tau(f)} \f$ as equations above
 */
Numeric naive_tropospheric_singletau_median(const ConstVectorView bt, const Numeric trop_bt, const Numeric target_bt);

/** Apply tropospheric correction
 * 
 * Tropospheric correction comes here from
 * \f[ T_{b, 1}(f) = \frac{T_{b, 0}(f) - \overline{T_{b, trop}(f)}(1 - e^{-\overline{\tau(f)}})}{e^{-\overline{\tau(f)}}}, \f]
 * 
 * where
 * \f$ T_{b, 1}(f) \f$ is the corrected radiation to some higher altitude,
 * \f$ T_{b, 0}(f) \f$ is the measured brightness temperature at surface altitude,
 * \f$ \overline{T_{b, trop}(f)} \f$ is the averaged brightness temperature emission of the troposphere and 
 * \f$ \overline{\tau(f)} \f$ is the averaged finite opacity of the troposphere.
 * Here \f$ f \f$ is the frequency dimension of bt
 * 
 * This method assumes \f$ \overline{T_{b, trop}(f)} \f$ and  \f$ \overline{\tau(f)} \f$ are 
 * constant for all frequencies.  They have a frequency dependency, but we naively
 * ignore it, thus giving the method its name
 * 
 * Note that \f$ \overline{\tau(f)} \f$ so large that \f$ e^{-\overline{\tau(f)}} := 0 \f$ 
 * are ignored and the original input is retained for such an input
 * 
 * @param[in,out] bt \f$ T_{b, 0}(f) \f$
 * @param[in] tau \f$ \overline{\tau(f)} \f$
 * @param[in] trop_bt \f$ \overline{T_{b, trop}(f)} \f$
 * @return \f$ T_{b, 1}(f) \f$
 */
VectorView naive_tropospheric(VectorView bt, const Numeric tau, const Numeric trop_bt);
}  // Correction
} // Raw
#endif  // RAW_H

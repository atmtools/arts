/* Copyright (C) 2019
   Richard Larsson
                            
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

/**
 * @file radiation_field.h
 * @author Richard Larsson
 * @date 2019-09-04
 * 
 * @brief Radiation field calculations
 */

#ifndef radiation_field_h
#define radiation_field_h

#include "matpackI.h"
#include "mystring.h"
#include "rte.h"

/** Throws an error if integration values are bad
 * 
 * @param[in] error_msg Error message to print
 * @param[in] value_that_should_be_unity Check value
 */
void error_in_integrate(const String& error_msg,
                        const Numeric& value_that_should_be_unity);

/** Convolved integration over the line
 * 
 * f must be sorted and in same order as F
 * 
 * @param[in] F Line shape
 * @param[in] f Frequency grod
 * @return Numeric Integrated line shape
 */
Numeric test_integrate_convolved(const Eigen::Ref<Eigen::VectorXcd> F,
                                 const Vector& f);

/** Integrate cos(za)
 * 
 * @param[in] cosza cos of zenith angle
 * @param[in] sorted_index Order of zenith angles
 * @return Numeric Integration
 */
Numeric test_integrate_zenith(const Vector& cosza,
                              const Array<Index>& sorted_index);

/** Integrate the convolved strength of the absorption of the line
 * 
 * f must be sorted and in same order as F and I
 * 
 * F must be normalized
 * 
 * @param[in] I Intensity vector
 * @param[in] F Line shape normalized
 * @param[in] f Frequency grid
 * @return Numeric Integrated absorption
 */
Numeric integrate_convolved(const RadiationVector& I,
                            const Eigen::VectorXcd& F,
                            const Vector& f);

/** Integrate the convolved strength ofthe transmission of the line
 * 
 * f must be sorted and in same order as F and T
 * 
 * F must be normalized
 * 
 * Only consider [0, 0] position of T
 * 
 * @param[in] T Transmission matrix
 * @param[in] F Line shape normalized
 * @param[in] f Frequency grid
 * @return Numeric Integrated absorption
 */
Numeric integrate_convolved(const TransmissionMatrix& T,
                            const Eigen::VectorXcd& F,
                            const Vector& f);

/** Integrate the convolved source function over zenith angles
 * 
 * @param[in] j Source vector
 * @param[in] cosza cos of zenith angle
 * @param[in] sorted_index Order of zenith angles
 * @return Numeric Integrade source function
 */
Numeric integrate_zenith(const VectorView j,
                         const Vector& cosza,
                         const Array<Index>& sorted_index);

/** Get a position from grid pos
 * 
 * Assumes the grid pos has been set to 
 * extend over the path
 * 
 * @param[in] gp Grid position
 * @return Index Position
 */
Index grid_index_from_gp(const GridPos& gp);

/** Get sorting of zenith angles in field of ppath
 * 
 * @param[in] sorted_index Order of zenith angles
 * @param[in] cosza cos of zenith angle
 * @param[in] ppath_field As WSV
 */
void sorted_index_of_ppath_field(ArrayOfArrayOfIndex& sorted_index,
                                 ArrayOfVector& cosza,
                                 const ArrayOfPpath& ppath_field);

#endif  // radiation_field_h

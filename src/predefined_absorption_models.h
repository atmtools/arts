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

/*!
 * @file   fullmodel.h
 * @author Richard Larsson
 * @date   2020-01-29
 * 
 * @brief  Full absorption models of various kinds
 */


#ifndef fullmodel_h
#define fullmodel_h

#include <Faddeeva/Faddeeva.hh>
#include "linefunctions.h"

namespace Absorption {
namespace PredefinedModel {

/** Adds Makarov MPM2020 O2 absorption lines to the absorption matrix
 * 
 * Adds negative values outside of bounds.  Works for Earth only
 * 
 * Water is 10% more effective at broadening than dry air so must
 * be included.  Does not deal with Zeeman effect and ignores negative
 * absorption far from the line center (this can be fixed but requires
 * clear-cut motivation)
 * 
 * @param[in,out] xsec Cross-section of oxygen (size: [f x p])
 * @param[in,out] dxsec Cross-section derivatives of oxygen (size: [jacs_pos][f x p])
 * @param[in]     f Frequency grid of computations (size: [f])
 * @param[in]     p Pressure grid of computations (size: [p])
 * @param[in]     t Temperature grid of computations (size: [p])
 * @param[in]     water_vmr Water volume mixing ratio (size: [p])
 * @param[in]     jacs The Jacobian descriptions (size: [greater than max(jacs_pos)])
 * @param[in]     jacs_pos The Jacobian matrix positions in dxsec (size: [jacs_pos])
 */
void makarov2020_o2_lines_mpm(Matrix& xsec,
                              ArrayOfMatrix& dxsec,
                              const Vector& f,
                              const Vector& p,
                              const Vector& t,
                              const Vector& water_vmr,
                              const ArrayOfRetrievalQuantity& jacs,
                              const ArrayOfIndex& jacs_pos);

};  //PredefinedModel 
};  //Absorption

#endif  // fullmodel_h

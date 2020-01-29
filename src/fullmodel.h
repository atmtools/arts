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

struct makarov2020_o2_lines_control {
  Index pos_in_xsec;
  Index frequency;
  Index pressure;
  Index temperature;
  Index f0;
  Index intens;
  Index a2;
  Index gamma;
  Index y0;
  Index y1;
  Index g0;
  Index g1;
  Index dv0;
  Index dv1;
  Index x;
};

void makarov2020_o2_lines(ArrayOfMatrix& xsec,
                          ArrayOfArrayOfMatrix& dxsec,
                          const Vector& f,
                          const Vector& p,
                          const Vector& t,
                          const makarov2020_o2_lines_control& ctrl);

#endif  // fullmodel_h

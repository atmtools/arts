/* Copyright (C) 2003-2012 Nikolay Koulev <nkoulev@uni-bremen.de>
  
   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.
  
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
  
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

/*!
  \file   geomag_calc.h
  \author Nikolay Koulev <nkoulev@uni-bremen.de>
  \date   Mon Jul 28 11:35:24 2003
  
 \brief The header file for the functions in geomag_calc.cc


*/

#ifndef geomag_calc_h
#define geomag_calc_h

#include "matpack.h"

void magfield_nk(   // Output
    Numeric& B_r,   // radial component of the geomagnetic field
    Numeric& B_th,  // colatitudinal component of the geomagnetic field
    Numeric& B_ph,  // longitudinal component of the geomagnetic field

    // Input
    const Numeric r,      // radial distance to the point
    const Numeric theta,  // geocentric colatitude of the point
    const Numeric phi,    // longitude of the point
    // All coordinates - geocentric!

    const Index Ny  // number of elapsed years after an epoch year, J - [0,4]
);
#endif

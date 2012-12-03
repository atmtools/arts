/* Copyright (C) 2012 Stefan Buehler <sbuehler@ltu.se>
 
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

/** \file
 Implementation file for work with HITRAN collision induced absorption (CIA).
 
 The CIA data are part of the HITRAN distribution. They are described in
 Richard, C., I. E. Gordon, L. S. Rothman, M. Abel, L. Frommhold, M. Gustafsson,
 J.-M. Hartmann, C. Hermans, W. J. Lafferty, G. S. Orton, K.M. Smith, and H. Tran (2012),
 New section of the HITRAN database: Collision-induced absorption (CIA),
 J. Quant. Spectrosc. Radiat. Transfer, 113, 1276-1285, doi:10.1016/j.jqsrt.2011.11.004.
 
 \author Stefan Buehler
 \date   2012-11-30
 */

#include "cia.h"

/** Interpolate CIA data.
 
 Interpolate CIA data to given frequency vector and given scalar temperature.
 Uses second order interpolation in both coordinates, if grid length allows,
 otherwise first order or no interpolation.
 
 /retval result CIA value for given frequency grid and temperature.
 /param frequency Frequency grid
 /param temperature Scalar temparature
 /param cia_data The CIA dataset to interpolate */
void cia_interpolation(VectorView result,
                       ConstVectorView frequency,
                       const Numeric& temperature,
                       const GriddedField2& cia_data)
{
    // Get data grids:
    ConstVectorView f_grid = cia_data.get_numeric_grid(1);
    ConstVectorView T_grid = cia_data.get_numeric_grid(2);

    // Decide on interpolation orders:
    const Index f_order = 2;

    // For T we have to be adaptive, since sometimes there is only one T in
    // the data
    Index T_order;
    switch (T_grid.nelem()) {
        case 1:
            T_order = 0;
            break;
        case 2:
            T_order = 1;
            break;
        default:
            T_order = 2;
            break;
    }
    
    // Check if frequency is inside the range covered by the data:
    chk_interpolation_grids("Frequency interpolation for CIA continuum",
                                f_grid,
                                frequency,
                                f_order);
    

    // Check if temperature is inside the range covered by the data:
    if (T_order > 0) {
        chk_interpolation_grids("Temperature interpolation for CIA continuum",
                                T_grid,
                                temperature,
                                T_order);
    }

}

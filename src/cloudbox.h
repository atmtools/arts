/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
   USA. 
*/

 /*!
     \file   cloudbox.h
     \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
     \date   Thu May  23 14:34:05 2002
     
     \brief  Internal cloudbox functions.
     
   */

#ifndef cloudbox_h
#define cloudbox_h

#include "matpackVII.h"
#include "interpolation.h"

void cloudbox_getOutgoing1D(// Output:
                            MatrixView i_out,
                            //WS Input:
                            ConstTensor7View scat_i_p,
                            const GridPos& rte_gp_p,
                            ConstVectorView rte_los,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& stokes_dim,
                            ConstVectorView scat_za_grid,
                            ConstVectorView f_grid
                            );

void cloudbox_getOutgoingCubic1D(// Output:
                            MatrixView i_out,
                            //WS Input:
                            ConstTensor7View scat_i_p,
                            const GridPos& rte_gp_p,
                            ConstVectorView rte_los,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& stokes_dim,
                            ConstVectorView scat_za_grid,
                            ConstVectorView f_grid
                            );

void cloudbox_getOutgoing3D(// Output:
                            MatrixView   i_out,
                            // Input:
                            ConstTensor7View scat_i_p,
                            ConstTensor7View scat_i_lat,
                            ConstTensor7View scat_i_lon,
                            const GridPos& rte_gp_p,
                            const GridPos& rte_gp_lat,
                            const GridPos& rte_gp_lon,
                            ConstVectorView rte_los,
                            const ArrayOfIndex& cloudbox_limits,
                            const Index& stokes_dim,
                            ConstVectorView scat_za_grid,
                            ConstVectorView scat_aa_grid,
                            ConstVectorView f_grid);

void cloudbox_boundary_check(// Output:
                             bool& on_p_bd,
                             bool& on_lat_bd,
                             bool& on_lon_bd,
                             GridPos& cloud_gp_p,
                             GridPos& cloud_gp_lat,
                             GridPos& cloud_gp_lon,
                             // Input:
                             const ArrayOfIndex& cloudbox_limits
                             );

#endif //cloudbox_h


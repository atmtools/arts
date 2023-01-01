/* Copyright (C) 2021 Patrick Eriksson <patrick.eriksson@chalmers.se>

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
 * @file   ppath.cc
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2023-01-01
 *
 * @brief  Functions releated to calculation of propagation paths.
 *
 * The term propagation path is here shortened to ppath.
 */

/*===========================================================================
  === External declarations
  ===========================================================================*/


#include <algorithm>

#include "geodeticZZZ.h"
#include "ppathZZZ.h"
#include "variousZZZ.h"

bool ppath_l2toa_from_above(Numeric& l2toa,
                            ConstVectorView rte_pos,
                            ConstVectorView rte_los,
                            ConstVectorView ecef,
                            ConstVectorView decef,
                            const Vector& refellipsoid,
                            const Numeric& z_toa)
{
  // Cases that are inside atmosphere
  if (rte_pos[0] < z_toa || (rte_pos[0] == z_toa && rte_los[0] > 90)) {
    l2toa = 0;
    return false;

  // Outside of the atmosphere 
  } else {
    if (rte_los[0] > 90) {
      // Returns negative if no intersection
      l2toa = intersection_altitude(ecef, decef, refellipsoid, z_toa);
    }
    return true;
  }
}

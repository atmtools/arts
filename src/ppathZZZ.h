/* Copyright (C) 2021 Patrick Eriksson <Patrick.Eriksson@chalmers.se>

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
 * @file   ppath.h
 * @author Patrick Eriksson <patrick.eriksson@chalmers.se>
 * @date   2023-01-01
 * 
 * @brief  Functions releated to calculation of propagation paths.
 */

#ifndef ppathZZZ_h
#define ppathZZZ_h

#include "agenda_class.h"
#include "gridded_fields.h"
#include "interpolation.h"


/*===========================================================================
  === Functions
  ===========================================================================*/


/** Locates rte_pos with respect to the top of the atmosphere

   @param[out]  l2toa         Length to TOA. Set to zero if sensor is
                              inside  of atmosphere, and to -1 if TOA 
                              not reached from above.
   @param[in]   rte_pos       As the WSV with the same name.
   @param[in]   rte_los       As the WSV with the same name.
   @param[in]   ecef          rte_pos in ECEF.
   @param[in]   decef         rte_los in ECEF.
   @param[in]   refellipsoid  As the WSV with same name.
   @param[in]   z_toa         Top-of-the-atmosphere altitude

   @return  True if rte_pos is above z_toa, false otherwise

   @author Patrick Eriksson
   @date   2023-01-01
 */
bool ppath_l2toa_from_above(Numeric& l2toa,
                            ConstVectorView rte_pos,
                            ConstVectorView rte_los,
                            ConstVectorView ecef,
                            ConstVectorView decef,
                            const Vector& refellipsoid,
                            const Numeric& z_toa);

#endif  // ppathZZZ_h

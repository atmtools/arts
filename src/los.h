/* Copyright (C) 2000 Patrick Eriksson <patrick@rss.chalmers.se>

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

/*-----------------------------------------------------------------------
FILE:      los.h

INCLUDES:  The declaration of the Los data type, moved here from workspace.h

FUNCTIONS: None

HISTORY:   10.06.2000 Created by Stefan Buehler
-----------------------------------------------------------------------*/

#ifndef los_h
#define los_h

#include "vecmat.h"

/** The line of sight (LOS). The LOS structure has the fields:
    \begin{verbatim}
       ARRAYofVECTOR  p;
       VECTOR         l_step;
       ARRAY<int>     ground;
       ARRAY<int>     start;
       ARRAY<int>     stop;
    where 
       p        The pressures along LOS
       l_step   The geometrical length along LOS between the points.
       start    start index for the iteration
       stop     stop index for the iteration
       ground   O if no intersection with the ground. Else, GROUND
                gives the index for the ground.  
    \end{verbatim}

    The LOS is defined in equal long geometrical steps along the path.
    This step length (L_STEP) can vary between the viewing angles.

    Spectra are calculated in the following way (by RTE_ITERATE in m_los):
    \begin{enumerate}
    \item Iteration from START down to 1 or GROUND
    \item If GROUND, including the effect of the ground reflection.
    \item Iteration from 1 or GROUND-1 to STOP
    \end{enumerate}

    The START and STOP variables make it possible to use a possible symmetry
    for 1D calculations. For example, for limb sounding from space, START
    and STOP are both set to the length of P. The GROUND variable is for
    1D calculations either 0 or 1.

    For cases without symmetry (upward looking and 2D), STOP is always 1
    and corresponds to the point closest to the sensor. Accordingly, START
    corresponds to the point of LOS furthest away from the sensor.

    The GROUND variable is used both as a flag to indicate ground 
    intersections of the LOS, and a variable to give the position of the
    ground. As mentioned, for 1D cases, the ground is always placed at 
    index 1. For 2D cases, GROUND gives the index for the ground point, 
    that is, the point of LOS with index GROUND corresponds to the ground 
    level.

    @author Patrick Eriksson 07.06.00 */
struct Los {
  ARRAYofVECTOR  p;
  VECTOR         l_step;
  ARRAY<int>     ground;
  ARRAY<int>     start;
  ARRAY<int>     stop;
};

#endif  // los_h

/* Copyright (C) 2000 Stefan Buehler <sbuehler@uni-bremen.de>
                      Patrick Eriksson <patrick@rss.chalmers.se>

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
FILE:      workspace.h

INCLUDES:  The declaration of the workspace data type. You have to
           edit this file if you want to add a new workspace variable.

FUNCTIONS: None

CLASSES:   Workspace

HISTORY:   ??.??.???? Created by Stefan Buehler
           10.06.2000 Cleaned up by Stefan Buehler
-----------------------------------------------------------------------*/

#ifndef workspace_h
#define workspace_h

#include "arts.h"
#include "absorption.h"
#include "los.h"
#include "wsv_group.h"
#include "wsv_aux.h"


/** The declaration of the (great) workspace.
    @author Stefan Buehler */
class WorkSpace {
public:
  //                       < Spectroscopy Stuff >
  ARRAYofLineRecord        lines;
  ARRAYofARRAYofLineRecord lines_per_tg;
  TagGroups           	   tag_groups;
  //                  	   < 1D Input Atmosphere Stuff >
  MATRIX              	   raw_ptz_1d;
  ARRAYofMATRIX       	   raw_vmrs_1d;
  //                  	   < General Absorption Stuff >
  VECTOR              	   p_abs;
  VECTOR              	   f_abs;
  //                  	   < 2D Absorption Stuff >
  ARRAYofVECTOR       	   t_abs_2d;
  ARRAYofVECTOR       	   z_abs_2d;
  ARRAYofMATRIX       	   vmrs_2d;
  ARRAYofMATRIX       	   abs_2d;
  //                  	   < 1D Absorption Stuff >
  VECTOR              	   t_abs;
  VECTOR              	   z_abs;
  ARRAYofVECTOR       	   vmrs;
  MATRIX              	   abs;
  ARRAYofMATRIX       	   abs_per_tg;
  //                  	   < RT Stuff >
  VECTOR              	   view1;
  Numeric             	   z_plat;
  Numeric             	   l_step;
  int                 	   refr;
  Numeric             	   l_step_refr;
  VECTOR              	   refr_index;
  Numeric             	   z_ground;
  Numeric             	   t_ground;
  VECTOR              	   e_ground;
  Los                 	   los; 
  ARRAYofMATRIX       	   source;
  ARRAYofMATRIX       	   trans;
  VECTOR              	   y_space;
  VECTOR              	   y;
  //                  	   < WF stuff >
  ARRAYofMATRIX       	   klos;
  VECTOR              	   k_grid;
  MATRIX                   k;
};



#endif  // workshpace_h

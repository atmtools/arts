/* Copyright (C) 2000, 2001 Patrick Eriksson <patrick@rss.chalmers.se>

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



////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   ppath.h
  \author Patrick Eriksson <Patrick.Eriksson@rss.chalmers.se>
  \date   Sun May  5  2002
  
  \brief  Propagation path structure and functions.
  
   This file contains the definition of the Ppath structure and the
   functions in ppath.c that are of interest elsewhere.
*/


#ifndef ppath_h
#define ppath_h


#include <math.h>
#include "arts.h"
#include "check_input.h"
#include "math_funcs.h"
#include "matpackI.h"
#include "mystring.h"



////////////////////////////////////////////////////////////////////////////
//   The Ppath structure
////////////////////////////////////////////////////////////////////////////

struct Ppath {
  Index     dim;
  Index     np;
  Matrix    pos;
  Vector    z;
  Vector    l_step;
  Matrix    gridindex;
  Matrix    los;
  String    background;
  Index     ground;
  Index     i_ground;
  Vector    tan_pos;
  Index     symmetry;
  Index     i_symmetry;
};



////////////////////////////////////////////////////////////////////////////
//   Functions from ppath.cc
////////////////////////////////////////////////////////////////////////////



#endif  // ppath_h

/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

/*!
  \file   workspace_aux.cc
  \brief  Auxiliary material for the workspace

  This file contains auxiliary material for the workspace, which used
  to be in workspace.cc. The reason for the separation is that the
  stuff here hardly ever should be changed, whereas workspace.cc has
  to be edited each time a new variable is added.

  \author Stefan Buehler
  \date 2000-06-10 */

#include "arts.h"
#include <map>
#include "matpackI.h"
#include "array.h"
#include "auto_wsv_groups.h"
#include "wsv_aux.h"
#include "agenda_record.h"

/*! The map assiciated with wsv_data. */
map<String, Index> WsvMap;

void define_wsv_map()
{
  extern const Array<WsvRecord> wsv_data;
  extern map<String, Index> WsvMap;

  for ( Index i=0 ; i<wsv_data.nelem() ; ++i )
    {
      WsvMap[wsv_data[i].Name()] = i;
    }
}



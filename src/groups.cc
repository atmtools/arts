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

/*!
  \file   groups.cc
  \brief  Defines workspace variable groups.

  If you want to add new workspace variable groups you have to do it
  in this file. This is used by the program make_wsv_group_h to
  generate the header file wsv_group.h

  \author Stefan Buehler
  \date   2000-08-04 */

#include "arts.h"
#include "vecmat.h"

/*! The names associated with Wsv groups as strings.
  See function define_wsv_group_names for more information. */
ARRAY<string> wsv_group_names;


//! Define the array of workspace variable group names.
/*!
  This defines the global variable wsv_group_names. It is used in two
  different programs:

  1. In arts.

  2. In make_wsv_group_h.

  \author Stefan Buehler
  \date   2000-08-04
*/
void define_wsv_group_names()
{

  //--------------------< Build the group names array >--------------------
  // Initialize to empty, just in case.
  wsv_group_names.clear();

  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     Must be consistent with the enum in wsv_group.h! Later on,
     these enums could also be generated automatically, but that would
     have to take place before the wsv_data is defined, since that
     needs these enums.
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
  wsv_group_names.push_back("string");
  wsv_group_names.push_back("int");
  wsv_group_names.push_back("Numeric");
  wsv_group_names.push_back("VECTOR");
  wsv_group_names.push_back("MATRIX");
  wsv_group_names.push_back("ARRAYofMATRIX");
  wsv_group_names.push_back("ARRAYofVECTOR");
  wsv_group_names.push_back("Los");
  wsv_group_names.push_back("ARRAYofLineRecord");
  wsv_group_names.push_back("ARRAYofARRAYofLineRecord");
  wsv_group_names.push_back("TagGroups");
  wsv_group_names.push_back("SPARSEMATRIX");
}

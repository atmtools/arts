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

/*! The names associated with Wsv groups as Strings.
  See function define_wsv_group_names for more information. */
Array<String> wsv_group_names;


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
  resize(wsv_group_names,0);

  wsv_group_names.push_back("String");
  wsv_group_names.push_back("size_t");
  wsv_group_names.push_back("int");
  wsv_group_names.push_back("Numeric");
  wsv_group_names.push_back("Vector");
  wsv_group_names.push_back("Matrix");
  wsv_group_names.push_back("SYMMETRIC");
  wsv_group_names.push_back("ArrayofString");
  wsv_group_names.push_back("Arrayofsizet");
  wsv_group_names.push_back("ArrayofVector");
  wsv_group_names.push_back("ArrayofMatrix");
  wsv_group_names.push_back("LOS");
  wsv_group_names.push_back("ArrayofLineRecord");
  wsv_group_names.push_back("ArrayofArrayofLineRecord");
  wsv_group_names.push_back("TagGroups");
  wsv_group_names.push_back("Hmatrix");
  wsv_group_names.push_back("ArrayofLineshapeSpec");
}


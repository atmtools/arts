/* Copyright (C) 2000, 2001 Stefan Buehler <sbuehler@uni-bremen.de>
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
#include "array.h"
#include "mystring.h"
// #include "supergeneric.h"
// #include "ppath.h"
// #include "gas_abs_lookup.h"

/*! The names associated with Wsv groups as Strings.
  See function define_wsv_group_names for more information. */
ArrayOfString wsv_group_names;


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
  wsv_group_names.resize(0);

  wsv_group_names.push_back("Any");
  wsv_group_names.push_back("Index");
  wsv_group_names.push_back("Numeric");
  wsv_group_names.push_back("String");
  wsv_group_names.push_back("Vector");
  wsv_group_names.push_back("Matrix");
  wsv_group_names.push_back("Sparse");
  wsv_group_names.push_back("Tensor3");
  wsv_group_names.push_back("Tensor4");
  wsv_group_names.push_back("Tensor5");
  wsv_group_names.push_back("Tensor6");
  wsv_group_names.push_back("Tensor7");
  wsv_group_names.push_back("Timer");
  wsv_group_names.push_back("ArrayOfIndex");
  wsv_group_names.push_back("ArrayOfString");
  wsv_group_names.push_back("ArrayOfVector");
  wsv_group_names.push_back("ArrayOfMatrix");
  wsv_group_names.push_back("ArrayOfArrayOfMatrix");
  wsv_group_names.push_back("ArrayOfTensor3");
  wsv_group_names.push_back("ArrayOfArrayOfTensor3");
  wsv_group_names.push_back("ArrayOfTensor6");
  wsv_group_names.push_back("ArrayOfTensor7");
  wsv_group_names.push_back("ArrayOfArrayOfTensor6");
//   wsv_group_names.push_back("ArrayOfLineRecord");
//   wsv_group_names.push_back("ArrayOfArrayOfLineRecord");
  wsv_group_names.push_back("ArrayOfArrayOfSpeciesTag");
  wsv_group_names.push_back("Ppath");
  wsv_group_names.push_back("Agenda");
  wsv_group_names.push_back("GridPos");
  wsv_group_names.push_back("ArrayOfArrayOfArrayOfArrayOfGridPos");
  wsv_group_names.push_back("GasAbsLookup");
  wsv_group_names.push_back("SingleScatteringData");
  wsv_group_names.push_back("ArrayOfSingleScatteringData");
  wsv_group_names.push_back("GriddedField3");
  wsv_group_names.push_back("ArrayOfGriddedField3");
  wsv_group_names.push_back("ArrayOfRetrievalQuantity");
}


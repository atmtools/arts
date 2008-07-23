/* Copyright (C) 2000-2008
   Stefan Buehler <sbuehler@ltu.se>
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

#include <map>
#include "arts.h"
#include "array.h"
#include "mystring.h"


/*! The names associated with Wsv groups as Strings.
  See function define_wsv_group_names for more information. */
ArrayOfString wsv_group_names;
map<String, Index> WsvGroupMap;

/*! Groups that can be used as keywords */
ArrayOfIndex valid_keyword_groups;


void define_valid_keyword_groups()
{
  valid_keyword_groups.resize(0);
  valid_keyword_groups.push_back(get_wsv_group_id("String"));
  valid_keyword_groups.push_back(get_wsv_group_id("Index"));
  valid_keyword_groups.push_back(get_wsv_group_id("Numeric"));
  valid_keyword_groups.push_back(get_wsv_group_id("ArrayOfString"));
  valid_keyword_groups.push_back(get_wsv_group_id("ArrayOfIndex"));
  valid_keyword_groups.push_back(get_wsv_group_id("Vector"));
}


void define_wsv_group_map()
{
  for ( Index i=0 ; i<wsv_group_names.nelem() ; ++i )
    {
      WsvGroupMap[wsv_group_names[i]] = i;
    }
}


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
  wsv_group_names.push_back("ArrayOfArrayOfIndex");
  wsv_group_names.push_back("ArrayOfString");
  wsv_group_names.push_back("ArrayOfVector");
  wsv_group_names.push_back("ArrayOfMatrix");
  wsv_group_names.push_back("ArrayOfArrayOfMatrix");
  wsv_group_names.push_back("ArrayOfTensor3");
  wsv_group_names.push_back("ArrayOfArrayOfTensor3");
  wsv_group_names.push_back("ArrayOfTensor4");
  wsv_group_names.push_back("ArrayOfTensor6");
  wsv_group_names.push_back("ArrayOfTensor7");
  wsv_group_names.push_back("ArrayOfArrayOfTensor6");
  wsv_group_names.push_back("ArrayOfLineRecord");
  wsv_group_names.push_back("ArrayOfArrayOfLineRecord");
  wsv_group_names.push_back("ArrayOfLineshapeSpec");
  wsv_group_names.push_back("ArrayOfArrayOfSpeciesTag");
  wsv_group_names.push_back("Ppath");
  wsv_group_names.push_back("ArrayOfPpath");
  wsv_group_names.push_back("Agenda");
  wsv_group_names.push_back("GridPos");
  wsv_group_names.push_back("ArrayOfArrayOfArrayOfArrayOfGridPos");
  wsv_group_names.push_back("GasAbsLookup");
  wsv_group_names.push_back("SingleScatteringData");
  wsv_group_names.push_back("ArrayOfSingleScatteringData");
  wsv_group_names.push_back("GField1");
  wsv_group_names.push_back("GField2");
  wsv_group_names.push_back("GField3");
  wsv_group_names.push_back("GField4");
  wsv_group_names.push_back("ArrayOfGField1");
  wsv_group_names.push_back("ArrayOfGField2");
  wsv_group_names.push_back("ArrayOfGField3");
  wsv_group_names.push_back("ArrayOfGField4");
  wsv_group_names.push_back("ArrayOfArrayOfGField1");
  wsv_group_names.push_back("ArrayOfArrayOfGField3");
  wsv_group_names.push_back("ArrayOfRetrievalQuantity");
  wsv_group_names.push_back("MCAntenna");
  wsv_group_names.push_back("SLIData2");

  define_wsv_group_map();
  define_valid_keyword_groups();
}


bool is_valid_keyword_group(const Index group)
{
  for (Index i = 0; i < valid_keyword_groups.nelem(); i++)
    {
      if (valid_keyword_groups[i] == group)
        return true;
    }

  return false;
}

Index get_wsv_group_id(const String& name)
{
  map<String, Index>::const_iterator it = WsvGroupMap.find (name);
  if (it == WsvGroupMap.end())
    return -1;
  else
    return it->second;
}


/* Copyright (C) 2000-2012
   Stefan Buehler <sbuehler@ltu.se>
   Patrick Eriksson <patrick.eriksson@chalmers.se>

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
#include "array.h"
#include "arts.h"
#include "mystring.h"
#include "wsv_aux.h"

/*! The names associated with Wsv groups as Strings.
  See function define_wsv_group_names for more information. */
namespace global_data {
ArrayOfString wsv_group_names;
map<String, Index> WsvGroupMap;
}  // namespace global_data

/*! Groups that can be used as keywords */
ArrayOfIndex valid_keyword_groups;

void define_valid_keyword_groups() {
  valid_keyword_groups.resize(0);
  valid_keyword_groups.push_back(get_wsv_group_id("String"));
  valid_keyword_groups.push_back(get_wsv_group_id("Index"));
  valid_keyword_groups.push_back(get_wsv_group_id("Numeric"));
  valid_keyword_groups.push_back(get_wsv_group_id("ArrayOfString"));
  valid_keyword_groups.push_back(get_wsv_group_id("ArrayOfIndex"));
  valid_keyword_groups.push_back(get_wsv_group_id("Vector"));
}

void define_wsv_group_map() {
  using global_data::wsv_group_names;
  using global_data::WsvGroupMap;
  for (Index i = 0; i < wsv_group_names.nelem(); ++i) {
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
void define_wsv_group_names() {
  using global_data::wsv_group_names;

  //--------------------< Build the group names array >--------------------
  // Initialize to empty, just in case.
  wsv_group_names.resize(0);
  
  wsv_group_names.push_back("AbsorptionLines");
  wsv_group_names.push_back("Agenda");
  wsv_group_names.push_back("Any");
  wsv_group_names.push_back("ArrayOfAbsorptionLines");
  wsv_group_names.push_back("ArrayOfArrayOfAbsorptionLines");
  wsv_group_names.push_back("ArrayOfAgenda");
  wsv_group_names.push_back("ArrayOfArrayOfGriddedField1");
  wsv_group_names.push_back("ArrayOfArrayOfGriddedField2");
  wsv_group_names.push_back("ArrayOfArrayOfGriddedField3");
  wsv_group_names.push_back("ArrayOfArrayOfIndex");
  wsv_group_names.push_back("ArrayOfArrayOfMatrix");
  wsv_group_names.push_back("ArrayOfPpath");
  wsv_group_names.push_back("ArrayOfArrayOfPropagationMatrix");
  wsv_group_names.push_back("ArrayOfArrayOfRadiationVector");
  wsv_group_names.push_back("ArrayOfArrayOfScatteringMetaData");
  wsv_group_names.push_back("ArrayOfArrayOfSingleScatteringData");
  wsv_group_names.push_back("ArrayOfArrayOfSpeciesTag");
  wsv_group_names.push_back("ArrayOfArrayOfStokesVector");
  wsv_group_names.push_back("ArrayOfArrayOfString");
  wsv_group_names.push_back("ArrayOfArrayOfTensor3");
  wsv_group_names.push_back("ArrayOfArrayOfTensor6");
  wsv_group_names.push_back("ArrayOfArrayOfTransmissionMatrix");
  wsv_group_names.push_back("ArrayOfArrayOfVector");
  wsv_group_names.push_back("ArrayOfCIARecord");
  wsv_group_names.push_back("ArrayOfGriddedField1");
  wsv_group_names.push_back("ArrayOfGriddedField2");
  wsv_group_names.push_back("ArrayOfGriddedField3");
  wsv_group_names.push_back("ArrayOfGriddedField4");
  wsv_group_names.push_back("ArrayOfIndex");
  wsv_group_names.push_back("ArrayOfMatrix");
  wsv_group_names.push_back("ArrayOfPropagationMatrix");
  wsv_group_names.push_back("ArrayOfQuantumIdentifier");
  wsv_group_names.push_back("ArrayOfRadiationVector");
  wsv_group_names.push_back("ArrayOfRetrievalQuantity");
  wsv_group_names.push_back("ArrayOfScatteringMetaData");
  wsv_group_names.push_back("ArrayOfSingleScatteringData");
  wsv_group_names.push_back("ArrayOfSparse");
  wsv_group_names.push_back("ArrayOfStokesVector");
  wsv_group_names.push_back("ArrayOfString");
  wsv_group_names.push_back("ArrayOfTelsemAtlas");
  wsv_group_names.push_back("ArrayOfTensor3");
  wsv_group_names.push_back("ArrayOfTensor4");
  wsv_group_names.push_back("ArrayOfTensor5");
  wsv_group_names.push_back("ArrayOfTensor6");
  wsv_group_names.push_back("ArrayOfTensor7");
  wsv_group_names.push_back("ArrayOfTransmissionMatrix");
  wsv_group_names.push_back("ArrayOfVector");
  wsv_group_names.push_back("ArrayOfXsecRecord");
  wsv_group_names.push_back("CIARecord");
  wsv_group_names.push_back("CovarianceMatrix");
  wsv_group_names.push_back("EnergyLevelMap");
  wsv_group_names.push_back("GasAbsLookup");
  wsv_group_names.push_back("GridPos");
  wsv_group_names.push_back("GriddedField1");
  wsv_group_names.push_back("GriddedField2");
  wsv_group_names.push_back("GriddedField3");
  wsv_group_names.push_back("GriddedField4");
  wsv_group_names.push_back("GriddedField5");
  wsv_group_names.push_back("GriddedField6");
  wsv_group_names.push_back("Index");
  wsv_group_names.push_back("MCAntenna");
  wsv_group_names.push_back("Matrix");
  wsv_group_names.push_back("Numeric");
  wsv_group_names.push_back("Ppath");
  wsv_group_names.push_back("PropagationMatrix");
  wsv_group_names.push_back("QuantumIdentifier");
  wsv_group_names.push_back("RadiationVector");
  wsv_group_names.push_back("Rational");
  wsv_group_names.push_back("ScatteringMetaData");
  wsv_group_names.push_back("SingleScatteringData");
  wsv_group_names.push_back("Sparse");
  wsv_group_names.push_back("SpeciesAuxData");
  wsv_group_names.push_back("StokesVector");
  wsv_group_names.push_back("String");
  wsv_group_names.push_back("TelsemAtlas");
  wsv_group_names.push_back("Tensor3");
  wsv_group_names.push_back("Tensor4");
  wsv_group_names.push_back("Tensor5");
  wsv_group_names.push_back("Tensor6");
  wsv_group_names.push_back("Tensor7");
  wsv_group_names.push_back("Timer");
  wsv_group_names.push_back("TessemNN");
  wsv_group_names.push_back("TransmissionMatrix");
  wsv_group_names.push_back("Vector");
  wsv_group_names.push_back("Verbosity");

  define_wsv_group_map();
  define_valid_keyword_groups();
}

bool is_valid_keyword_group(const Index group) {
  for (Index i = 0; i < valid_keyword_groups.nelem(); i++) {
    if (valid_keyword_groups[i] == group) return true;
  }

  return false;
}

void get_wsv_group_ids(ArrayOfIndex& ids, String name) {
  ids.resize(0);

  Index pos = 0;
  while (pos < name.nelem()) {
    switch (name[pos]) {
      case ' ':
      case '\r':
      case '\t':
      case '#':
        name.erase(pos, 1);
        break;
      default:
        pos++;
    }
  }

  pos = 0;
  Index prev = 0;
  while (pos < name.nelem()) {
    while (pos < name.nelem() && name[pos] != ',') pos++;
    Index id = get_wsv_group_id(name.substr(prev, pos - prev));
    if (id == -1) {
      ids.resize(0);
      return;
    }
    ids.push_back(id);
    pos++;
    prev = pos;
  }
}

bool is_agenda_group_id(const Index group) {
  return (group == get_wsv_group_id("Agenda") ||
          group == get_wsv_group_id("ArrayOfAgenda"));
}

Index get_wsv_group_id(const String& name) {
  using global_data::WsvGroupMap;
  map<String, Index>::const_iterator it = WsvGroupMap.find(name);
  if (it == WsvGroupMap.end())
    return -1;
  else
    return it->second;
}

String get_array_groups_as_string(bool basetype_is_group,
                                  bool return_basetype_only) {
  using global_data::wsv_group_names;
  String arraygroups;

  bool first = true;
  for (Index i = 0; i < wsv_group_names.nelem(); i++) {
    if (wsv_group_names[i].substr(0, String("ArrayOf").length()) == "ArrayOf") {
      const String basetype = wsv_group_names[i].substr(
          String("ArrayOf").length(), wsv_group_names[i].length());
      bool basetype_exists = (get_wsv_group_id(basetype) != -1);

      if (return_basetype_only) {
        // Return only the basetype of the array,
        // skip arrays whose basetype is not a WSV group
        if (basetype_exists) {
          if (!first)
            arraygroups += ", ";
          else
            first = false;
          arraygroups += basetype;
        }
      } else {
        if (!basetype_is_group || (basetype_is_group && basetype_exists)) {
          if (!first)
            arraygroups += ", ";
          else
            first = false;
          arraygroups += wsv_group_names[i];
        }
      }
    }
  }
  return arraygroups;
}

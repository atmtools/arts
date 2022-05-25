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

#include "array.h"
#include "arts.h"
#include "groups.h"
#include "mystring.h"
#include "wsv_aux.h"
#include <map>

/*! The names associated with Wsv groups as Strings.
  See function define_wsv_groups for more information. */
namespace global_data {
ArrayOfGroupRecord wsv_groups;
map<String, Index> WsvGroupMap;
}  // namespace global_data

void define_wsv_group_map() {
  using global_data::wsv_groups;
  using global_data::WsvGroupMap;
  for (Index i = 0; i < wsv_groups.nelem(); ++i) {
    WsvGroupMap[wsv_groups[i]] = i;
  }
}

//! Define the array of workspace variable group names.
/*!
  This defines the global variable wsv_groups. It is used in two
  different programs:

  1. In arts.

  2. In make_wsv_group_h.

  \author Stefan Buehler
  \date   2000-08-04
*/
void define_wsv_groups() {
  using global_data::wsv_groups;

  //--------------------< Build the group names array >--------------------
  // Initialize to empty, just in case.
  wsv_groups.resize(0);
  
  wsv_groups.emplace_back("AbsorptionLines");
  wsv_groups.emplace_back("Agenda");
  wsv_groups.emplace_back("Any");
  wsv_groups.emplace_back("ArrayOfAbsorptionLines");
  wsv_groups.emplace_back("ArrayOfArrayOfAbsorptionLines");
  wsv_groups.emplace_back("ArrayOfAgenda");
  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField1");
  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField2");
  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField3");
  wsv_groups.emplace_back("ArrayOfArrayOfIndex");
  wsv_groups.emplace_back("ArrayOfArrayOfMatrix");
  wsv_groups.emplace_back("ArrayOfPpath");
  wsv_groups.emplace_back("ArrayOfArrayOfPropagationMatrix");
  wsv_groups.emplace_back("ArrayOfArrayOfRadiationVector");
  wsv_groups.emplace_back("ArrayOfArrayOfScatteringMetaData");
  wsv_groups.emplace_back("ArrayOfArrayOfSingleScatteringData");
  wsv_groups.emplace_back("ArrayOfArrayOfSpeciesTag");
  wsv_groups.emplace_back("ArrayOfArrayOfStokesVector");
  wsv_groups.emplace_back("ArrayOfArrayOfString");
  wsv_groups.emplace_back("ArrayOfArrayOfTensor3");
  wsv_groups.emplace_back("ArrayOfArrayOfTensor6");
  wsv_groups.emplace_back("ArrayOfArrayOfTime");
  wsv_groups.emplace_back("ArrayOfArrayOfTransmissionMatrix");
  wsv_groups.emplace_back("ArrayOfArrayOfVector");
  wsv_groups.emplace_back("ArrayOfCIARecord");
  wsv_groups.emplace_back("ArrayOfGriddedField1");
  wsv_groups.emplace_back("ArrayOfGriddedField2");
  wsv_groups.emplace_back("ArrayOfGriddedField3");
  wsv_groups.emplace_back("ArrayOfGriddedField4");
  wsv_groups.emplace_back("ArrayOfIndex");
  wsv_groups.emplace_back("ArrayOfJacobianTarget");
  wsv_groups.emplace_back("ArrayOfMatrix");
  wsv_groups.emplace_back("ArrayOfPropagationMatrix");
  wsv_groups.emplace_back("ArrayOfQuantumIdentifier");
  wsv_groups.emplace_back("ArrayOfRadiationVector");
  wsv_groups.emplace_back("ArrayOfRetrievalQuantity");
  wsv_groups.emplace_back("ArrayOfScatteringMetaData");
  wsv_groups.emplace_back("ArrayOfSingleScatteringData");
  wsv_groups.emplace_back("ArrayOfSpeciesTag");
  wsv_groups.emplace_back("ArrayOfSparse");
  wsv_groups.emplace_back("ArrayOfStokesVector");
  wsv_groups.emplace_back("ArrayOfString");
  wsv_groups.emplace_back("ArrayOfTelsemAtlas");
  wsv_groups.emplace_back("ArrayOfTensor3");
  wsv_groups.emplace_back("ArrayOfTensor4");
  wsv_groups.emplace_back("ArrayOfTensor5");
  wsv_groups.emplace_back("ArrayOfTensor6");
  wsv_groups.emplace_back("ArrayOfTensor7");
  wsv_groups.emplace_back("ArrayOfTime");
  wsv_groups.emplace_back("ArrayOfTransmissionMatrix");
  wsv_groups.emplace_back("ArrayOfVector");
  wsv_groups.emplace_back("ArrayOfXsecRecord");
  wsv_groups.emplace_back("CIARecord");
  wsv_groups.emplace_back("CallbackFunction");
  wsv_groups.emplace_back("CovarianceMatrix");
  wsv_groups.emplace_back("EnergyLevelMap");
  wsv_groups.emplace_back("GasAbsLookup");
  wsv_groups.emplace_back("GridPos");
  wsv_groups.emplace_back("GriddedField1");
  wsv_groups.emplace_back("GriddedField2");
  wsv_groups.emplace_back("GriddedField3");
  wsv_groups.emplace_back("GriddedField4");
  wsv_groups.emplace_back("GriddedField5");
  wsv_groups.emplace_back("GriddedField6");
  wsv_groups.emplace_back("HitranRelaxationMatrixData");
  wsv_groups.emplace_back("Index");
  wsv_groups.emplace_back("JacobianTarget");
  wsv_groups.emplace_back("MapOfErrorCorrectedSuddenData");
  wsv_groups.emplace_back("MCAntenna");
  wsv_groups.emplace_back("Matrix");
  wsv_groups.emplace_back("Numeric");
  wsv_groups.emplace_back("Ppath");
  wsv_groups.emplace_back("PropagationMatrix");
  wsv_groups.emplace_back("QuantumIdentifier");
  wsv_groups.emplace_back("RadiationVector");
  wsv_groups.emplace_back("Rational");
  wsv_groups.emplace_back("ScatteringMetaData");
  wsv_groups.emplace_back("SingleScatteringData");
  wsv_groups.emplace_back("Sparse");
  wsv_groups.emplace_back("SpeciesIsotopologueRatios");
  wsv_groups.emplace_back("StokesVector");
  wsv_groups.emplace_back("String");
  wsv_groups.emplace_back("TelsemAtlas");
  wsv_groups.emplace_back("Tensor3");
  wsv_groups.emplace_back("Tensor4");
  wsv_groups.emplace_back("Tensor5");
  wsv_groups.emplace_back("Tensor6");
  wsv_groups.emplace_back("Tensor7");
  wsv_groups.emplace_back("Timer");
  wsv_groups.emplace_back("Time");
  wsv_groups.emplace_back("TessemNN");
  wsv_groups.emplace_back("TransmissionMatrix");
  wsv_groups.emplace_back("Vector", R"--(A 1 dimensional array of numbers

Vector is used in a lot of math requiring vector semantics,
or just as a way to transport numbers between functions
)--");
  wsv_groups.emplace_back("Verbosity");

  define_wsv_group_map();
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
  auto it = WsvGroupMap.find(name);
  if (it == WsvGroupMap.end())
    return -1;
  return it->second;
}

String get_array_groups_as_string(bool basetype_is_group,
                                  bool return_basetype_only) {
  using global_data::wsv_groups;
  String arraygroups;

  bool first = true;
  for (Index i = 0; i < wsv_groups.nelem(); i++) {
    if (wsv_groups[i].name.substr(0, String("ArrayOf").length()) == "ArrayOf") {
      const String basetype = wsv_groups[i].name.substr(
          String("ArrayOf").length(), wsv_groups[i].name.length());
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
          arraygroups += wsv_groups[i];
        }
      }
    }
  }
  return arraygroups;
}

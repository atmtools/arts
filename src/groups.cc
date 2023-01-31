/*!
  \file   groups.cc
  \brief  Defines workspace variable groups.

  If you want to add new workspace variable groups you have to do it
  in this file. This is used by the program make_wsv_group_h to
  generate the header file wsv_group.h

  \author Stefan Buehler
  \date   2000-08-04 */

#include "groups.h"

#include <map>

#include "array.h"
#include "arts.h"
#include "mystring.h"

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

  wsv_groups.emplace_back(
      "AbsorptionLines",
      "Contains line-by-line absorption information for a number of related absorption lines");

  wsv_groups.emplace_back(
      "Agenda", "Describes a set of function calls and variable definitions");

  wsv_groups.emplace_back(
      "Any",
      "Meta type for when methods can take any argument (avoid manual use)");

  wsv_groups.emplace_back("ArrayOfAbsorptionLines",
                          "A list of *AbsorptionLines*");

  wsv_groups.emplace_back("ArrayOfArrayOfAbsorptionLines",
                          "A list of *ArrayOfAbsorptionLines*");

  wsv_groups.emplace_back("ArrayOfAgenda", "A list of *Agenda*");

  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField1",
                          "A list of *ArrayOfGriddedField1*");

  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField2",
                          "A list of *ArrayOfGriddedField2*");

  wsv_groups.emplace_back("ArrayOfArrayOfGriddedField3",
                          "A list of *ArrayOfGriddedField3*");

  wsv_groups.emplace_back("ArrayOfArrayOfIndex", "A list of *ArrayOfIndex*");

  wsv_groups.emplace_back("ArrayOfArrayOfMatrix", "A list of *ArrayOfMatrix*");

  wsv_groups.emplace_back("ArrayOfPpath", "A list of *Ppath*");

  wsv_groups.emplace_back("ArrayOfArrayOfPropagationMatrix",
                          "A list of *ArrayOfPropagationMatrix*");

  wsv_groups.emplace_back("ArrayOfArrayOfRadiationVector",
                          "A list of *ArrayOfRadiationVector*");

  wsv_groups.emplace_back("ArrayOfArrayOfScatteringMetaData",
                          "A list of *ArrayOfScatteringMetaData*");

  wsv_groups.emplace_back("ArrayOfArrayOfSingleScatteringData",
                          "A list of *ArrayOfSingleScatteringData*");

  wsv_groups.emplace_back("ArrayOfArrayOfSpeciesTag",
                          "A list of *ArrayOfSpeciesTag*");

  wsv_groups.emplace_back("ArrayOfArrayOfStokesVector",
                          "A list of *ArrayOfStokesVector*");

  wsv_groups.emplace_back("ArrayOfArrayOfString", "A list of *ArrayOfString*");

  wsv_groups.emplace_back("ArrayOfArrayOfTensor3",
                          "A list of *ArrayOfTensor3*");

  wsv_groups.emplace_back("ArrayOfArrayOfTensor6",
                          "A list of *ArrayOfTensor6*");

  wsv_groups.emplace_back("ArrayOfArrayOfTime", "A list of *ArrayOfTime*");

  wsv_groups.emplace_back("ArrayOfArrayOfTransmissionMatrix",
                          "A list of *ArrayOfTransmissionMatrix*");

  wsv_groups.emplace_back("ArrayOfArrayOfVector", "A list of *ArrayOfVector*");

  wsv_groups.emplace_back("ArrayOfCIARecord", "A list of *CIARecord*");

  wsv_groups.emplace_back("ArrayOfGriddedField1", "A list of *GriddedField1*");

  wsv_groups.emplace_back("ArrayOfGriddedField2", "A list of *GriddedField2*");

  wsv_groups.emplace_back("ArrayOfGriddedField3", "A list of *GriddedField3*");

  wsv_groups.emplace_back("ArrayOfGriddedField4", "A list of *GriddedField4*");

  wsv_groups.emplace_back("ArrayOfIndex", "A list of *Index*");

  wsv_groups.emplace_back("ArrayOfJacobianTarget",
                          "A list of *JacobianTarget*");

  wsv_groups.emplace_back("ArrayOfMatrix", "A list of *Matrix*");

  wsv_groups.emplace_back("ArrayOfPropagationMatrix",
                          "A list of *PropagationMatrix*");

  wsv_groups.emplace_back("ArrayOfQuantumIdentifier",
                          "A list of *QuantumIdentifier*");

  wsv_groups.emplace_back("ArrayOfRadiationVector",
                          "A list of *RadiationVector*");

  wsv_groups.emplace_back("ArrayOfRetrievalQuantity",
                          "A list of retrieval quantitities");

  wsv_groups.emplace_back("ArrayOfScatteringMetaData",
                          "A list of *ScatteringMetaData*");

  wsv_groups.emplace_back("ArrayOfSingleScatteringData",
                          "A list of *SingleScatteringData*");

  wsv_groups.emplace_back("ArrayOfSpeciesTag", R"--(A list of species tags

These tags include the species and a lot of optional information
about the isotopologue, the absorption scheme, and the frequency limits)--");

  wsv_groups.emplace_back("ArrayOfSparse", "A list of *Sparse*");

  wsv_groups.emplace_back("ArrayOfSun", "A list of sun");

  wsv_groups.emplace_back("ArrayOfStokesVector", "A list of *StokesVector*");

  wsv_groups.emplace_back("ArrayOfString", "A list of *String*");

  wsv_groups.emplace_back("ArrayOfTelsemAtlas", "A list of *TelsemAtlas*");

  wsv_groups.emplace_back("ArrayOfTensor3", "A list of *Tensor3*");

  wsv_groups.emplace_back("ArrayOfTensor4", "A list of *Tensor4*");

  wsv_groups.emplace_back("ArrayOfTensor5", "A list of *Tensor5*");

  wsv_groups.emplace_back("ArrayOfTensor6", "A list of *Tensor6*");

  wsv_groups.emplace_back("ArrayOfTensor7", "A list of *Tensor7*");

  wsv_groups.emplace_back("ArrayOfTime", "A list of *Time*");

  wsv_groups.emplace_back("ArrayOfTransmissionMatrix",
                          "A list of *TransmissionMatrix*");

  wsv_groups.emplace_back("ArrayOfVector", "A list of *Vector*");

  wsv_groups.emplace_back("ArrayOfXsecRecord",
                          R"--(A list of cross-section records

These cross-section records contains information about the valid temperature and
pressure ranges as well as well as the fitting coefficients used to compute
and interpolate the cross-section to other temperatures and pressures)--");

  wsv_groups.emplace_back("AtmField", R"--(An atmospheric field)--");

  wsv_groups.emplace_back("AtmPoint", R"--(An atmospheric point)--");

  wsv_groups.emplace_back(
      "CIARecord",
      R"--(Contains information to compute collision induced absorption for a pair of species

Holds an the record data in a gridded field with grids of temperature and frequency in
units of m^5 molec^(-2) )--");

  wsv_groups.emplace_back("CallbackFunction",
                          "Used to inject custom code into *Agenda*");

  wsv_groups.emplace_back("CovarianceMatrix", "Contains the covariance matrix");

  wsv_groups.emplace_back("EnergyLevelMap",
                          R"--(Maps data based on energy levels

Used for keeping track of non-local thermodynamic equilibrium data)--");

  wsv_groups.emplace_back("GasAbsLookup", R"--(An absorption lookup table

This class holds an absorption lookup table, as well as all
information that is necessary to use the table to extract
absorption)--");

  wsv_groups.emplace_back("GridPos", "A position in a grid");

  wsv_groups.emplace_back("GriddedField1",
                          R"--(A 1 dimensional gridded set of *Numeric* data

The grid is 1 *Vector* or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField2",
                          R"--(A 2 dimensional gridded set *Numeric* data

The grid is a combination of 2 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField3",
                          R"--(A 3 dimensional gridded set of *Numeric* data

The grid is a combination of 3 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField4",
                          R"--(A 4 dimensional gridded set of *Numeric* data

The grid is a combination of 4 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField5",
                          R"--(A 5 dimensional gridded set  of *Numeric* data

The grid is a combination of 5 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("GriddedField6",
                          R"--(A 6 dimensional gridded set of *Numeric* data

The grid is a combination of 6 *Vector* and/or *ArrayOfString*

Both the data and the grid may be named)--");

  wsv_groups.emplace_back("HitranRelaxationMatrixData",
                          "Wraps data required to use Hitran line mixing");

  wsv_groups.emplace_back("Index", "A 64 bit signed integer type");

  wsv_groups.emplace_back(
      "JacobianTarget", "A single target if a partial derivative computation");

  wsv_groups.emplace_back(
      "MapOfErrorCorrectedSuddenData",
      R"--(A map of data required for computing the error-corrected-sudden relaxation matrix

This map contains a list of an underlying data type.  This underlying data type contains a
*QuantumIdentifier* and a list of species dependent computational data for various components
required to compute the relaxation matrix

If there is no identifier or species avaialable, default values that approximates a diagonal
relaxation matrix are set)--");

  wsv_groups.emplace_back("MCAntenna", "An Antenna object used by *MCGeneral*");

  wsv_groups.emplace_back("Matrix", "A 2 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Numeric", "IEEE 754 binary64 floating point number");

  wsv_groups.emplace_back("Ppath", "Describes a propagation path");

  wsv_groups.emplace_back("PredefinedModelData",
                          R"--(Contains any data required for a predefined model)--");

  wsv_groups.emplace_back("PropagationMatrix",
                          R"--(The propagation matrix data is help by this type

This type is related to *StokesVector*

The data type is *Tensor4* in units of [1/m]

The dimensionality is kept as:

Number of frequencies as *Index* (usually from *f_grid*)
Number of zenith angles as *Index* 
Number of azimuth angles as *Index* 
The Stokes dimension as *Index* (usually from *stokes_dim*)

An individual propagation matrix (i.e., for a given frequency, zenith,
and azimuth angle) follows certain symmetries depending on the Stokes
dimension

For Stokes dimension 4:

K11  K12  K13  K14
K12  K11  K23  K24
K13 -K23  K11  K34
K14 -K24 -K34  K11

For Stokes dimension 3:

K11  K12  K13
K12  K11  K23
K13 -K23  K11

For Stokes dimension 2:

K11  K12
K12  K11

For Stokes dimension 1:

K11

The propagation matrix make use of these symmetries to computate the matrix inverses and exponents
required to turn the data into a *TransmissionMatrix* (with information about the distance))--");

  wsv_groups.emplace_back("QuantumIdentifier",
                          R"--(An ID for an absorption species state

It contains information about the species and a set of quantum numbers
and can thus be used to identify one of the following:
1) a species
2) an isotopologue of a species
3) an absorption band of an isotopologue
4) an absorption line of an isotopologue
5) the energy level of absorption band(s) of an isotopologue
6) the energy level of absorption line(s) of an isotopologue)--");

  wsv_groups.emplace_back(
      "RadiationVector",
      R"--(Contains the radiation vector as a function of frequency

This type is related to *TransmissionMatrix*

The stokes dimensionality translates directly to the size of the vector

Internally, this holds an efficiently packed list of these vectors

This is often used in combination with *TransmissionMatrix* to compute the radiative
transfer through the atmosphere

It holds information about the radiance, unlike its cousin *StokesVector*, which holds information
about the vector absorption/emission)--");

  wsv_groups.emplace_back("Rational",
                          "Holds a rational number as two *Index* n / d");

  wsv_groups.emplace_back("ScatteringMetaData",
                          "Holds meta data about the scattering");

  wsv_groups.emplace_back("SingleScatteringData",
                          "Holds single scattering data");

  wsv_groups.emplace_back("Sparse", "A sparse version of *Matrix*");

  wsv_groups.emplace_back(
      "SpeciesIsotopologueRatios",
      "Contains a list of isotopologue ratios for all defined species");

  wsv_groups.emplace_back("StokesVector", R"--(A stokes vector

This type is related to *PropagationMatrix*

The data type is *Tensor4* in units of [1/m]

The dimensionality is kept as:

Number of frequencies as *Index* (usually from *f_grid*)
Number of zenith angles as *Index* 
Number of azimuth angles as *Index* 
The Stokes dimension as *Index* (usually from *stokes_dim*)

This is often used to compute the source emission with the help of a *PropagationMatrix*)--");

  wsv_groups.emplace_back("String", "Basic string type");

  wsv_groups.emplace_back("TelsemAtlas", R"--(A telsem atlas

Represents a Telsem2 atlas containing land surface microwave emissivities.
Since the Atlas contains emissivities only for land surfaces, the data is
stored in a sparse format.
 
The emissivities are represented on an equal area grid and numbered
sequentially starting with the first latitude band at -90 degrees and
moving up to 90 degrees.

The correspondance array contains the data indices for each cellnumber
if it is contained in the Atlas and NAN otherwise.)--");

  wsv_groups.emplace_back("Tensor3", "A 3 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor4", "A 4 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor5", "A 5 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor6", "A 6 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Tensor7", "A 7 dimensional array of *Numeric*");

  wsv_groups.emplace_back("Timer", "Represents a clock");

  wsv_groups.emplace_back("Time", R"(Represents a time stamp in the format:
"YEAR-MONTH-DAY HOUR:MINUTE:SECOND", e.g., "2023-03-06 14:32:35.35"

Note that most direct user input of a Time accepts a string as above to
represent the time stamp.
)");

  wsv_groups.emplace_back(
      "TessemNN", "Data required by TESSEM to calculate surface emissivity");

  wsv_groups.emplace_back(
      "TransmissionMatrix",
      R"--(Contains the transmission matrix as a function of frequency

This type is related to *RadiationVector*

The stokes dimensionality squared translates directly to the size of the matrix

Internally, this holds an efficiently packed list of these matrices

This is often used in combination with *RadiationVector* to compute the radiative
transfer through the atmosphere

The transmission matrix is often computed from the combination of two *PropagationMatrix*
at different atmospheric path points (using the distance between these points)

It holds information about the polarized transmission, unlike its cousin *PropagationMatrix*,
which holds information about the polarized absorption)--");

  wsv_groups.emplace_back("Vector", "A 1 dimensional array of *Numeric*");

  wsv_groups.emplace_back(
      "Verbosity",
      "Controls the screen, agenda, and file verbosity level (i.e. the level of information printed)");

  wsv_groups.emplace_back("VibrationalEnergyLevels", "A map of vibrational energy levels for NLTE calculations");

  std::sort(wsv_groups.begin(), wsv_groups.end(), [](auto& a, auto& b) {
    return a.name < b.name;
  });

  define_wsv_group_map();
}

Index get_wsv_group_id(const String& name) {
  using global_data::WsvGroupMap;
  auto it = WsvGroupMap.find(name);
  if (it == WsvGroupMap.end()) return -1;
  return it->second;
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

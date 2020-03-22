/* Copyright (C) 2020 Richard Larsson <larsson@mps.mpg.de>

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License as
   published by the Free Software Foundation; either version 2 of the
   License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA. */

////////////////////////////////////////////////////////////////////////////
//   File description
////////////////////////////////////////////////////////////////////////////
/*!
  \file   arts_api_classes.cc
  \author Richard Larsson <larsson@mps.mpg.de>
  \date   2020-03-12

  \brief This file contains all declarations of the ARTS C API class interface.
*/

#include "arts_api_classes.h"

#include "absorption.h"
#include "absorptionlines.h"
#include "energylevelmap.h"
#include "global_data.h"
#include "lineshapemodel.h"
#include "quantum.h"
#include "xml_io.h"
#include "xml_io_types.h"
#include "zeemandata.h"


#define BasicInterfaceCAPI(TYPE)                            \
void * create##TYPE()                                       \
{                                                           \
    return new TYPE;                                        \
}                                                           \
                                                            \
void delete##TYPE(void * data)                              \
{                                                           \
    delete static_cast<TYPE *>(data);                       \
}                                                           \
                                                            \
void print##TYPE(void * data)                               \
{                                                           \
  std::cout << (*static_cast<TYPE *>(data)) << std::endl;   \
}


#define GetterSetterCAPI(TYPE, VALUE, BASETYPE)     \
BASETYPE get##VALUE##TYPE(void * data)              \
{                                                   \
    return static_cast<TYPE *>(data) -> VALUE ();   \
}                                                   \
void set##VALUE##TYPE(void * data, BASETYPE newval) \
{                                                   \
    static_cast<TYPE *>(data) -> VALUE (newval);    \
}


#define EnumGetterSetterCAPI(TYPE, VALUE, ENUM)                       \
Index get##VALUE##TYPE(void * data)                                   \
{                                                                     \
  return Index(static_cast<TYPE *>(data) -> VALUE ());                \
}                                                                     \
Index set##VALUE##TYPE(void * data, Index newval)                     \
{                                                                     \
  if (static_cast<TYPE *>(data) -> validIndexFor##VALUE (newval)) {   \
    static_cast<TYPE *>(data) -> VALUE (ENUM(newval));                \
    return EXIT_SUCCESS;                                              \
  } else {                                                            \
    return EXIT_FAILURE;                                              \
  }                                                                   \
}                                                                     \
Index string2index##VALUE##TYPE(void * data, char * newval)           \
{                                                                     \
  return Index(static_cast<TYPE *>(data) -> string2##VALUE(newval));  \
}

#define VoidGetterCAPI(TYPE, VALUE)                 \
void * get##VALUE##TYPE(void * data)                \
{                                                   \
    return &static_cast<TYPE *>(data) -> VALUE();   \
}


#define BasicInputOutputCAPI(TYPE)                                                                                      \
Index xmlread##TYPE(void * data, char * filepath)                                                                       \
{                                                                                                                       \
    try {                                                                                                               \
        xml_read_from_file(filepath, *static_cast<TYPE *>(data), Verbosity());                                          \
        return EXIT_SUCCESS;                                                                                            \
    } catch (std::runtime_error& e) {                                                                                   \
        return EXIT_FAILURE;                                                                                            \
    }                                                                                                                   \
}                                                                                                                       \
                                                                                                                        \
Index xmlsave##TYPE(void * data, char * filepath, Index filetype, Index clobber)                                        \
{                                                                                                                       \
    try {                                                                                                               \
        xml_write_to_file(filepath, *static_cast<const TYPE *>(data), FileType(filetype), not clobber, Verbosity());    \
        return EXIT_SUCCESS;                                                                                            \
    } catch (std::runtime_error& e) {                                                                                   \
        return EXIT_FAILURE;                                                                                            \
    }                                                                                                                   \
}


#define VoidArrayCAPI(TYPE)                           \
Index size##TYPE(void * data)                         \
{                                                     \
  return static_cast<TYPE *>(data) -> size();         \
}                                                     \
void resize##TYPE(Index n, void * data)               \
{                                                     \
  static_cast<TYPE *>(data) -> resize(n);             \
}                                                     \
void * getelem##TYPE(Index i, void * data)            \
{                                                     \
  return &static_cast<TYPE *>(data) -> operator[](i); \
}


#define VoidArrayElemCAPI(TYPE, ELEM)                   \
Index size##ELEM##TYPE(void * data)                     \
{                                                       \
  return static_cast<TYPE *>(data) -> ELEM().size();    \
}                                                       \
void resize##ELEM##TYPE(Index n, void * data)           \
{                                                       \
  return static_cast<TYPE *>(data) -> ELEM().resize(n); \
}                                                       \
void * getelem##ELEM##TYPE(Index i, void * data)        \
{                                                       \
  return &static_cast<TYPE *>(data) -> ELEM()[i];       \
}


// Index
BasicInterfaceCAPI(Index)
BasicInputOutputCAPI(Index)
Index getIndex(void * data) { return *static_cast<Index *>(data); }
void setIndex(void * data, Index newval) { *static_cast<Index *>(data) = newval; }
VoidArrayCAPI(ArrayOfIndex)
BasicInterfaceCAPI(ArrayOfIndex)
BasicInputOutputCAPI(ArrayOfIndex)
VoidArrayCAPI(ArrayOfArrayOfIndex)
BasicInterfaceCAPI(ArrayOfArrayOfIndex)
BasicInputOutputCAPI(ArrayOfArrayOfIndex)


// Numeric
BasicInterfaceCAPI(Numeric)
BasicInputOutputCAPI(Numeric)
Numeric getNumeric(void * data) { return *static_cast<Numeric *>(data); }
void setNumeric(void * data, Numeric newval) { *static_cast<Numeric *>(data) = newval; }


// ZeemanModel
BasicInterfaceCAPI(ZeemanModel)
GetterSetterCAPI(ZeemanModel, gu, Numeric)
GetterSetterCAPI(ZeemanModel, gl, Numeric)


// Rational
BasicInterfaceCAPI(Rational)
BasicInputOutputCAPI(Rational)
GetterSetterCAPI(Rational, Nom, Index)
GetterSetterCAPI(Rational, Denom, Index)
void simplifyRational(void * data) { static_cast<Rational *>(data)->simplify_in_place(); }


// LineShape::ModelParameters
void printLineShapeModelParameters(void * data) { std::cout << (*static_cast<LineShape::ModelParameters *>(data)) << std::endl; }
Index getLineShapeModelParametersType(char * data) { try { return Index(LineShape::string2temperaturemodel(data)); } catch (std::runtime_error& e) { return -1; } }


// LineShape::SingleSpeciesModel
BasicInterfaceCAPI(LineShapeSingleSpeciesModel)
VoidGetterCAPI(LineShapeSingleSpeciesModel, G0)
VoidGetterCAPI(LineShapeSingleSpeciesModel, D0)
VoidGetterCAPI(LineShapeSingleSpeciesModel, G2)
VoidGetterCAPI(LineShapeSingleSpeciesModel, D2)
VoidGetterCAPI(LineShapeSingleSpeciesModel, FVC)
VoidGetterCAPI(LineShapeSingleSpeciesModel, ETA)
VoidGetterCAPI(LineShapeSingleSpeciesModel, Y)
VoidGetterCAPI(LineShapeSingleSpeciesModel, G)
VoidGetterCAPI(LineShapeSingleSpeciesModel, DV)


// LineShape::Model
BasicInterfaceCAPI(LineShapeModel)
VoidArrayCAPI(LineShapeModel)


// Absorption::SingleLine
BasicInterfaceCAPI(AbsorptionSingleLine)
GetterSetterCAPI(AbsorptionSingleLine, F0, Numeric)
GetterSetterCAPI(AbsorptionSingleLine, I0, Numeric)
GetterSetterCAPI(AbsorptionSingleLine, E0, Numeric)
GetterSetterCAPI(AbsorptionSingleLine, g_low, Numeric)
GetterSetterCAPI(AbsorptionSingleLine, g_upp, Numeric)
GetterSetterCAPI(AbsorptionSingleLine, A, Numeric)
VoidGetterCAPI(AbsorptionSingleLine, Zeeman)
VoidGetterCAPI(AbsorptionSingleLine, LineShape)
VoidArrayElemCAPI(AbsorptionSingleLine, LowerQuantumNumbers)
VoidArrayElemCAPI(AbsorptionSingleLine, UpperQuantumNumbers)


// QuantumNumbers
BasicInterfaceCAPI(QuantumNumbers)
void * getelemQuantumNumbers(Index i, void * data) { return &static_cast<QuantumNumbers *>(data)->operator[](i); }
Index sizeQuantumNumbers() { return Index(QuantumNumberType::FINAL_ENTRY); }
Index string2quantumnumbersindex(char * str) { return Index(string2quantumnumbertype(str)); }


// QuantumIdentifier
BasicInterfaceCAPI(QuantumIdentifier)
BasicInputOutputCAPI(QuantumIdentifier)
EnumGetterSetterCAPI(QuantumIdentifier, Type, QuantumIdentifier::QType)
GetterSetterCAPI(QuantumIdentifier, Species, Index)
GetterSetterCAPI(QuantumIdentifier, Isotopologue, Index)
VoidGetterCAPI(QuantumIdentifier, EnergyLevelQuantumNumbers)
VoidGetterCAPI(QuantumIdentifier, LowerQuantumNumbers)
VoidGetterCAPI(QuantumIdentifier, UpperQuantumNumbers)

// ArrayOfQuantumIdentifier
BasicInterfaceCAPI(ArrayOfQuantumIdentifier)
BasicInputOutputCAPI(ArrayOfQuantumIdentifier)
VoidArrayCAPI(ArrayOfQuantumIdentifier)

// SpeciesTag
BasicInterfaceCAPI(SpeciesTag)
BasicInputOutputCAPI(SpeciesTag)
GetterSetterCAPI(SpeciesTag, Species, Index)
GetterSetterCAPI(SpeciesTag, Isotopologue, Index)
GetterSetterCAPI(SpeciesTag, Uf, Numeric)
GetterSetterCAPI(SpeciesTag, Lf, Numeric)
GetterSetterCAPI(SpeciesTag, CIASecond, Index)
GetterSetterCAPI(SpeciesTag, CIADataset, Index)
EnumGetterSetterCAPI(SpeciesTag, Type, Index)
VoidArrayCAPI(ArrayOfSpeciesTag)
BasicInterfaceCAPI(ArrayOfSpeciesTag)
BasicInputOutputCAPI(ArrayOfSpeciesTag)
VoidArrayCAPI(ArrayOfArrayOfSpeciesTag)
BasicInterfaceCAPI(ArrayOfArrayOfSpeciesTag)
BasicInputOutputCAPI(ArrayOfArrayOfSpeciesTag)

Index setSpeciesTag(void * data, char * newdata)
{
  try {
    *static_cast<SpeciesTag *>(data) = SpeciesTag(newdata);
    return EXIT_SUCCESS;
  } catch(std::exception& e) {
    return EXIT_FAILURE;
  }
}

Index validSpecies(Index spec)
{
  if (spec >= 0 and spec < global_data::species_data.nelem())
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

Index validAllIsotopologues(Index spec, Index isot)
{
  auto& species = global_data::species_data[spec];
  if (isot == species.Isotopologue().nelem())
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

Index validIsotopologue(Index spec, Index isot)
{
  auto& species = global_data::species_data[spec];
  if (isot >= 0 and isot < species.Isotopologue().nelem() and not species.Isotopologue()[isot].isContinuum())
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}

Index validContinuum(Index spec, Index isot)
{
  auto& species = global_data::species_data[spec];
  if (isot >= 0 and isot < species.Isotopologue().nelem() and species.Isotopologue()[isot].isContinuum())
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}


// AbsorptionLines
BasicInterfaceCAPI(AbsorptionLines)
BasicInputOutputCAPI(AbsorptionLines)
GetterSetterCAPI(AbsorptionLines, Self, bool)
GetterSetterCAPI(AbsorptionLines, Bath, bool)
EnumGetterSetterCAPI(AbsorptionLines, Cutoff, Absorption::CutoffType)
EnumGetterSetterCAPI(AbsorptionLines, LineShapeType, LineShape::Type)
EnumGetterSetterCAPI(AbsorptionLines, Mirroring, Absorption::MirroringType)
EnumGetterSetterCAPI(AbsorptionLines, Population, Absorption::PopulationType)
EnumGetterSetterCAPI(AbsorptionLines, Normalization, Absorption::NormalizationType)
GetterSetterCAPI(AbsorptionLines, T0, Numeric)
GetterSetterCAPI(AbsorptionLines, CutoffFreqValue, Numeric)
GetterSetterCAPI(AbsorptionLines, LinemixingLimit, Numeric)
VoidGetterCAPI(AbsorptionLines, QuantumIdentity)
VoidArrayElemCAPI(AbsorptionLines, LocalQuanta)
VoidGetterCAPI(AbsorptionLines, BroadeningSpecies)
VoidArrayElemCAPI(AbsorptionLines, AllLines)
VoidArrayCAPI(ArrayOfAbsorptionLines)
BasicInterfaceCAPI(ArrayOfAbsorptionLines)
BasicInputOutputCAPI(ArrayOfAbsorptionLines)
VoidArrayCAPI(ArrayOfArrayOfAbsorptionLines)
BasicInterfaceCAPI(ArrayOfArrayOfAbsorptionLines)
BasicInputOutputCAPI(ArrayOfArrayOfAbsorptionLines)
void printmetaAbsorptionLines(void * data) { std::cout << static_cast<AbsorptionLines *>(data) -> MetaData() << std::endl; }
Index isAbsorptionLinesOK(void * data) { return Index(static_cast<AbsorptionLines *>(data) -> OK()); }


// EnergyLevelMap
BasicInterfaceCAPI(EnergyLevelMap)
BasicInputOutputCAPI(EnergyLevelMap)
EnumGetterSetterCAPI(EnergyLevelMap, Type, EnergyLevelMapType)
VoidGetterCAPI(EnergyLevelMap, Levels)
VoidGetterCAPI(EnergyLevelMap, Energies)
VoidGetterCAPI(EnergyLevelMap, Data)


// generic
Index string2filetypeindex(char * data) { try { return Index(string2filetype(data)); } catch (std::runtime_error& e) { return -1; } }


#undef BasicInterfaceCAPI
#undef GetterSetterCAPI
#undef EnumGetterSetterCAPI
#undef VoidGetterCAPI
#undef BasicInputOutputCAPI
#undef VoidArrayCAPI
#undef VoidArrayElemCAPI

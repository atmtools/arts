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
#include "covariance_matrix.h"
#include "energylevelmap.h"
#include "global_data.h"
#include "lineshapemodel.h"
#include "quantum.h"
#include "supergeneric.h"
#include "xml_io.h"
#include "xml_io_types.h"
#include "zeemandata.h"

#ifdef TIME_SUPPORT
#include <unistd.h>
#endif

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

#define VoidStructGetterCAPI(TYPE, VALUE)     \
void * get##VALUE##TYPE(void * data)          \
{                                             \
  return &static_cast<TYPE *>(data) -> VALUE; \
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
  static_cast<TYPE *>(data) -> ELEM().resize(n);        \
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
bool getOKEnergyLevelMap(void * data) {return static_cast<EnergyLevelMap *>(data) -> OK();}


// Vector
BasicInterfaceCAPI(Vector)
BasicInputOutputCAPI(Vector)
VoidArrayCAPI(ArrayOfVector)
BasicInterfaceCAPI(ArrayOfVector)
BasicInputOutputCAPI(ArrayOfVector)
VoidArrayCAPI(ArrayOfArrayOfVector)
BasicInterfaceCAPI(ArrayOfArrayOfVector)
BasicInputOutputCAPI(ArrayOfArrayOfVector)
void resizeVector(Index n, void * data) {static_cast<Vector *>(data) -> resize(n);}
Index nelemVector(void * data) {return static_cast<Vector *>(data) -> nelem();}
Numeric * getDataVector(void * data) {return static_cast<Vector *>(data) -> get_c_array();}


// Matrix
BasicInterfaceCAPI(Matrix)
BasicInputOutputCAPI(Matrix)
VoidArrayCAPI(ArrayOfMatrix)
BasicInterfaceCAPI(ArrayOfMatrix)
BasicInputOutputCAPI(ArrayOfMatrix)
VoidArrayCAPI(ArrayOfArrayOfMatrix)
BasicInterfaceCAPI(ArrayOfArrayOfMatrix)
BasicInputOutputCAPI(ArrayOfArrayOfMatrix)
void resizeMatrix(Index nrows, Index ncols, void * data) {static_cast<Matrix *>(data) -> resize(nrows, ncols);}
Index rowsMatrix(void * data) {return static_cast<Matrix *>(data) -> nrows();}
Index colsMatrix(void * data) {return static_cast<Matrix *>(data) -> ncols();}
Numeric * getDataMatrix(void * data) {return static_cast<Matrix *>(data) -> get_c_array();}


// Tensor3
BasicInterfaceCAPI(Tensor3)
BasicInputOutputCAPI(Tensor3)
VoidArrayCAPI(ArrayOfTensor3)
BasicInterfaceCAPI(ArrayOfTensor3)
BasicInputOutputCAPI(ArrayOfTensor3)
VoidArrayCAPI(ArrayOfArrayOfTensor3)
BasicInterfaceCAPI(ArrayOfArrayOfTensor3)
BasicInputOutputCAPI(ArrayOfArrayOfTensor3)
void resizeTensor3(Index npages, Index nrows, Index ncols, void * data) {static_cast<Tensor3 *>(data) -> resize(npages, nrows, ncols);}
Index pagesTensor3(void * data) {return static_cast<Tensor3 *>(data) -> npages();}
Index rowsTensor3(void * data) {return static_cast<Tensor3 *>(data) -> nrows();}
Index colsTensor3(void * data) {return static_cast<Tensor3 *>(data) -> ncols();}
Numeric * getDataTensor3(void * data) {return static_cast<Tensor3 *>(data) -> get_c_array();}


// Tensor4
BasicInterfaceCAPI(Tensor4)
BasicInputOutputCAPI(Tensor4)
VoidArrayCAPI(ArrayOfTensor4)
BasicInterfaceCAPI(ArrayOfTensor4)
BasicInputOutputCAPI(ArrayOfTensor4)
// VoidArrayCAPI(ArrayOfArrayOfTensor4)
// BasicInterfaceCAPI(ArrayOfArrayOfTensor4)
// BasicInputOutputCAPI(ArrayOfArrayOfTensor4)
void resizeTensor4(Index nbooks, Index npages, Index nrows, Index ncols, void * data) {static_cast<Tensor4 *>(data) -> resize(nbooks, npages, nrows, ncols);}
Index booksTensor4(void * data) {return static_cast<Tensor4 *>(data) -> nbooks();}
Index pagesTensor4(void * data) {return static_cast<Tensor4 *>(data) -> npages();}
Index rowsTensor4(void * data) {return static_cast<Tensor4 *>(data) -> nrows();}
Index colsTensor4(void * data) {return static_cast<Tensor4 *>(data) -> ncols();}
Numeric * getDataTensor4(void * data) {return static_cast<Tensor4 *>(data) -> get_c_array();}


// Tensor5
BasicInterfaceCAPI(Tensor5)
BasicInputOutputCAPI(Tensor5)
VoidArrayCAPI(ArrayOfTensor5)
BasicInterfaceCAPI(ArrayOfTensor5)
BasicInputOutputCAPI(ArrayOfTensor5)
// VoidArrayCAPI(ArrayOfArrayOfTensor5)
// BasicInterfaceCAPI(ArrayOfArrayOfTensor5)
// BasicInputOutputCAPI(ArrayOfArrayOfTensor5)
void resizeTensor5(Index nshelves, Index nbooks, Index npages, Index nrows, Index ncols, void * data) {static_cast<Tensor5 *>(data) -> resize(nshelves, nbooks, npages, nrows, ncols);}
Index shelvesTensor5(void * data) {return static_cast<Tensor5 *>(data) -> nshelves();}
Index booksTensor5(void * data) {return static_cast<Tensor5 *>(data) -> nbooks();}
Index pagesTensor5(void * data) {return static_cast<Tensor5 *>(data) -> npages();}
Index rowsTensor5(void * data) {return static_cast<Tensor5 *>(data) -> nrows();}
Index colsTensor5(void * data) {return static_cast<Tensor5 *>(data) -> ncols();}
Numeric * getDataTensor5(void * data) {return static_cast<Tensor5 *>(data) -> get_c_array();}


// Tensor6
BasicInterfaceCAPI(Tensor6)
BasicInputOutputCAPI(Tensor6)
VoidArrayCAPI(ArrayOfTensor6)
BasicInterfaceCAPI(ArrayOfTensor6)
BasicInputOutputCAPI(ArrayOfTensor6)
VoidArrayCAPI(ArrayOfArrayOfTensor6)
BasicInterfaceCAPI(ArrayOfArrayOfTensor6)
BasicInputOutputCAPI(ArrayOfArrayOfTensor6)
void resizeTensor6(Index nvitrines, Index nshelves, Index nbooks, Index npages, Index nrows, Index ncols, void * data) {static_cast<Tensor6 *>(data) -> resize(nvitrines, nshelves, nbooks, npages, nrows, ncols);}
Index vitrinesTensor6(void * data) {return static_cast<Tensor6 *>(data) -> nvitrines();}
Index shelvesTensor6(void * data) {return static_cast<Tensor6 *>(data) -> nshelves();}
Index booksTensor6(void * data) {return static_cast<Tensor6 *>(data) -> nbooks();}
Index pagesTensor6(void * data) {return static_cast<Tensor6 *>(data) -> npages();}
Index rowsTensor6(void * data) {return static_cast<Tensor6 *>(data) -> nrows();}
Index colsTensor6(void * data) {return static_cast<Tensor6 *>(data) -> ncols();}
Numeric * getDataTensor6(void * data) {return static_cast<Tensor6 *>(data) -> get_c_array();}


// Tensor7
BasicInterfaceCAPI(Tensor7)
BasicInputOutputCAPI(Tensor7)
VoidArrayCAPI(ArrayOfTensor7)
BasicInterfaceCAPI(ArrayOfTensor7)
BasicInputOutputCAPI(ArrayOfTensor7)
// VoidArrayCAPI(ArrayOfArrayOfTensor7)
// BasicInterfaceCAPI(ArrayOfArrayOfTensor7)
// BasicInputOutputCAPI(ArrayOfArrayOfTensor7)
void resizeTensor7(Index nlibraries, Index nvitrines, Index nshelves, Index nbooks, Index npages, Index nrows, Index ncols, void * data) {static_cast<Tensor7 *>(data) -> resize(nlibraries, nvitrines, nshelves, nbooks, npages, nrows, ncols);}
Index librariesTensor7(void * data) {return static_cast<Tensor7 *>(data) -> nlibraries();}
Index vitrinesTensor7(void * data) {return static_cast<Tensor7 *>(data) -> nvitrines();}
Index shelvesTensor7(void * data) {return static_cast<Tensor7 *>(data) -> nshelves();}
Index booksTensor7(void * data) {return static_cast<Tensor7 *>(data) -> nbooks();}
Index pagesTensor7(void * data) {return static_cast<Tensor7 *>(data) -> npages();}
Index rowsTensor7(void * data) {return static_cast<Tensor7 *>(data) -> nrows();}
Index colsTensor7(void * data) {return static_cast<Tensor7 *>(data) -> ncols();}
Numeric * getDataTensor7(void * data) {return static_cast<Tensor7 *>(data) -> get_c_array();}


// PropagationMatrix
BasicInterfaceCAPI(PropagationMatrix)
BasicInputOutputCAPI(PropagationMatrix)
VoidGetterCAPI(PropagationMatrix, Data)
VoidArrayCAPI(ArrayOfPropagationMatrix)
BasicInterfaceCAPI(ArrayOfPropagationMatrix)
BasicInputOutputCAPI(ArrayOfPropagationMatrix)
VoidArrayCAPI(ArrayOfArrayOfPropagationMatrix)
BasicInterfaceCAPI(ArrayOfArrayOfPropagationMatrix)
BasicInputOutputCAPI(ArrayOfArrayOfPropagationMatrix)
Index stokesPropagationMatrix(void * data) {return static_cast<PropagationMatrix *>(data) -> StokesDimensions();}
Index frequenciesPropagationMatrix(void * data) {return static_cast<PropagationMatrix *>(data) -> NumberOfFrequencies();}
Index zenithsPropagationMatrix(void * data) {return static_cast<PropagationMatrix *>(data) -> NumberOfZenithAngles();}
Index azimuthsPropagationMatrix(void * data) {return static_cast<PropagationMatrix *>(data) -> NumberOfAzimuthAngles();}
Index setPropagationMatrix(void * data, Index f, Index s, Index z, Index a, Numeric v)
{
  if (s >= 0 and s < 5 and f >= 0 and z >= 0 and a >= 0) {
    static_cast<PropagationMatrix *>(data) -> operator=(PropagationMatrix(f, s, z, a, v));
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
bool getOKPropagationMatrix(void * data) {return static_cast<PropagationMatrix *>(data) -> OK();}


// StokesVector
BasicInterfaceCAPI(StokesVector)
BasicInputOutputCAPI(StokesVector)
VoidGetterCAPI(StokesVector, Data)
VoidArrayCAPI(ArrayOfStokesVector)
BasicInterfaceCAPI(ArrayOfStokesVector)
BasicInputOutputCAPI(ArrayOfStokesVector)
VoidArrayCAPI(ArrayOfArrayOfStokesVector)
BasicInterfaceCAPI(ArrayOfArrayOfStokesVector)
BasicInputOutputCAPI(ArrayOfArrayOfStokesVector)
Index stokesStokesVector(void * data) {return static_cast<StokesVector *>(data) -> StokesDimensions();}
Index frequenciesStokesVector(void * data) {return static_cast<StokesVector *>(data) -> NumberOfFrequencies();}
Index zenithsStokesVector(void * data) {return static_cast<StokesVector *>(data) -> NumberOfZenithAngles();}
Index azimuthsStokesVector(void * data) {return static_cast<StokesVector *>(data) -> NumberOfAzimuthAngles();}
Index setStokesVector(void * data, Index f, Index s, Index z, Index a, Numeric v)
{
  if (s >= 0 and s < 5 and f >= 0 and z >= 0 and a >= 0) {
    static_cast<StokesVector *>(data) -> operator=(StokesVector(f, s, z, a, v));
    return EXIT_SUCCESS;
  } else {
    return EXIT_FAILURE;
  }
}
bool getOKStokesVector(void * data) {return static_cast<StokesVector *>(data) -> OK();}


// String
BasicInterfaceCAPI(String)
BasicInputOutputCAPI(String)
VoidArrayCAPI(ArrayOfString)
BasicInterfaceCAPI(ArrayOfString)
BasicInputOutputCAPI(ArrayOfString)
VoidArrayCAPI(ArrayOfArrayOfString)
BasicInterfaceCAPI(ArrayOfArrayOfString)
BasicInputOutputCAPI(ArrayOfArrayOfString)
void setString(void * data, char * newdata) {static_cast<String *>(data) -> operator=(String(newdata));}
char * getString(void * data) {return const_cast<char *>(static_cast<String *>(data) -> data());}


// GridPos
BasicInterfaceCAPI(GridPos)
BasicInputOutputCAPI(GridPos)
VoidStructGetterCAPI(GridPos, idx)
VoidStructGetterCAPI(GridPos, fd)
VoidArrayCAPI(ArrayOfGridPos)
BasicInterfaceCAPI(ArrayOfGridPos)
BasicInputOutputCAPI(ArrayOfGridPos)


// GridPosPoly
BasicInterfaceCAPI(GridPosPoly)
VoidStructGetterCAPI(GridPosPoly, idx)
VoidStructGetterCAPI(GridPosPoly, w)
VoidArrayCAPI(ArrayOfGridPosPoly)
BasicInterfaceCAPI(ArrayOfGridPosPoly)
Index xmlreadArrayOfGridPosPoly(void *, char *) {return 1;}
Index xmlsaveArrayOfGridPosPoly(void *, char *, Index, Index) {return 1;}


// Ppath
// BasicInterfaceCAPI(Ppath)
void * createPpath() {return new Ppath;}
void deletePpath(void * data) {delete static_cast<Ppath *>(data);}
void printPpath(void *) {std::cout << std::endl;}
BasicInputOutputCAPI(Ppath)
VoidStructGetterCAPI(Ppath, dim)
VoidStructGetterCAPI(Ppath, np)
VoidStructGetterCAPI(Ppath, constant)
VoidStructGetterCAPI(Ppath, background)
VoidStructGetterCAPI(Ppath, start_pos)
VoidStructGetterCAPI(Ppath, start_los)
VoidStructGetterCAPI(Ppath, start_lstep)
VoidStructGetterCAPI(Ppath, pos)
VoidStructGetterCAPI(Ppath, los)
VoidStructGetterCAPI(Ppath, r)
VoidStructGetterCAPI(Ppath, lstep)
VoidStructGetterCAPI(Ppath, end_pos)
VoidStructGetterCAPI(Ppath, end_los)
VoidStructGetterCAPI(Ppath, end_lstep)
VoidStructGetterCAPI(Ppath, nreal)
VoidStructGetterCAPI(Ppath, ngroup)
VoidStructGetterCAPI(Ppath, gp_p)
VoidStructGetterCAPI(Ppath, gp_lat)
VoidStructGetterCAPI(Ppath, gp_lon)
VoidArrayCAPI(ArrayOfPpath)
void * createArrayOfPpath() {return new ArrayOfPpath;}
void deleteArrayOfPpath(void * data) {delete static_cast<ArrayOfPpath *>(data);}
void printArrayOfPpath(void *) {std::cout << std::endl;}
BasicInputOutputCAPI(ArrayOfPpath)


// TransmissionMatrix
BasicInterfaceCAPI(TransmissionMatrix)
BasicInputOutputCAPI(TransmissionMatrix)
VoidArrayCAPI(ArrayOfTransmissionMatrix)
BasicInterfaceCAPI(ArrayOfTransmissionMatrix)
BasicInputOutputCAPI(ArrayOfTransmissionMatrix)
VoidArrayCAPI(ArrayOfArrayOfTransmissionMatrix)
BasicInterfaceCAPI(ArrayOfArrayOfTransmissionMatrix)
BasicInputOutputCAPI(ArrayOfArrayOfTransmissionMatrix)
Numeric * getMat1TransmissionMatrix(Index i, void * data) {return static_cast<TransmissionMatrix *>(data) -> Mat1(i).data();}
Numeric * getMat2TransmissionMatrix(Index i, void * data) {return static_cast<TransmissionMatrix *>(data) -> Mat2(i).data();}
Numeric * getMat3TransmissionMatrix(Index i, void * data) {return static_cast<TransmissionMatrix *>(data) -> Mat3(i).data();}
Numeric * getMat4TransmissionMatrix(Index i, void * data) {return static_cast<TransmissionMatrix *>(data) -> Mat4(i).data();}
void setTransmissionMatrix(void * data, Index stokes, Index freqs) {static_cast<TransmissionMatrix *>(data) -> operator=(TransmissionMatrix(freqs, stokes));}
Index getStokesDimTransmissionMatrix(void * data) {return static_cast<TransmissionMatrix *>(data) -> StokesDim();}
Index getFrequenciesTransmissionMatrix(void * data) {return static_cast<TransmissionMatrix *>(data) -> Frequencies();}


// RadiationVector
BasicInterfaceCAPI(RadiationVector)
BasicInputOutputCAPI(RadiationVector)
VoidArrayCAPI(ArrayOfRadiationVector)
BasicInterfaceCAPI(ArrayOfRadiationVector)
BasicInputOutputCAPI(ArrayOfRadiationVector)
VoidArrayCAPI(ArrayOfArrayOfRadiationVector)
BasicInterfaceCAPI(ArrayOfArrayOfRadiationVector)
BasicInputOutputCAPI(ArrayOfArrayOfRadiationVector)
Numeric * getVec1RadiationVector(Index i, void * data) {return static_cast<RadiationVector *>(data) -> Vec1(i).data();}
Numeric * getVec2RadiationVector(Index i, void * data) {return static_cast<RadiationVector *>(data) -> Vec2(i).data();}
Numeric * getVec3RadiationVector(Index i, void * data) {return static_cast<RadiationVector *>(data) -> Vec3(i).data();}
Numeric * getVec4RadiationVector(Index i, void * data) {return static_cast<RadiationVector *>(data) -> Vec4(i).data();}
void setRadiationVector(void * data, Index stokes, Index freqs) {static_cast<RadiationVector *>(data) -> operator=(RadiationVector(freqs, stokes));}
Index getStokesDimRadiationVector(void * data) {return static_cast<RadiationVector *>(data) -> StokesDim();}
Index getFrequenciesRadiationVector(void * data) {return static_cast<RadiationVector *>(data) -> Frequencies();}


// GriddedField1
BasicInterfaceCAPI(GriddedField1)
BasicInputOutputCAPI(GriddedField1)
VoidArrayCAPI(ArrayOfGriddedField1)
BasicInterfaceCAPI(ArrayOfGriddedField1)
BasicInputOutputCAPI(ArrayOfGriddedField1)
VoidArrayCAPI(ArrayOfArrayOfGriddedField1)
BasicInterfaceCAPI(ArrayOfArrayOfGriddedField1)
BasicInputOutputCAPI(ArrayOfArrayOfGriddedField1)
Index get_dimGriddedField1(void * data) {return static_cast<GriddedField1 *>(data) -> get_dim();}
Index get_grid_typeIndexGriddedField1(Index i, void * data) {return Index(static_cast<GriddedField1 *>(data) -> get_grid_type(i));}
Index get_grid_sizeGriddedField1(Index i, void * data) {return static_cast<GriddedField1 *>(data) -> get_grid_size(i);}
char * get_nameGriddedField1(void * data) {return const_cast<char *>(static_cast<GriddedField1 *>(data) -> get_name().data());}
void set_nameGriddedField1(void * data, char * newdata) {static_cast<GriddedField1 *>(data) -> set_name(newdata);}
char * get_grid_nameGriddedField1(Index i, void * data) {return const_cast<char *>(static_cast<GriddedField1 *>(data) -> get_grid_name(i).data());}
void set_grid_nameGriddedField1(Index i, void * data, char * newdata) {static_cast<GriddedField1 *>(data) -> set_grid_name(i, newdata);}
void * get_numeric_gridGriddedField1(Index i, void * data) {return &static_cast<GriddedField1 *>(data) -> get_numeric_grid(i);}
void * get_string_gridGriddedField1(Index i, void * data) {return &static_cast<GriddedField1 *>(data) -> get_string_grid(i);}
void set_gridGriddedField1(Index i, void * data, void * newdata, bool NumericType) 
{
  if (NumericType)
    static_cast<GriddedField1 *>(data) -> set_grid(i, *static_cast<Vector *>(newdata));
  else
    static_cast<GriddedField1 *>(data) -> set_grid(i, *static_cast<ArrayOfString *>(newdata));
}
void * dataGriddedField1(void * data) {return &static_cast<GriddedField1 *>(data) -> data;}
bool checksizeGriddedField1(void * data) {return static_cast<GriddedField1 *>(data) -> checksize();}


// GriddedField2
BasicInterfaceCAPI(GriddedField2)
BasicInputOutputCAPI(GriddedField2)
VoidArrayCAPI(ArrayOfGriddedField2)
BasicInterfaceCAPI(ArrayOfGriddedField2)
BasicInputOutputCAPI(ArrayOfGriddedField2)
VoidArrayCAPI(ArrayOfArrayOfGriddedField2)
BasicInterfaceCAPI(ArrayOfArrayOfGriddedField2)
BasicInputOutputCAPI(ArrayOfArrayOfGriddedField2)
Index get_dimGriddedField2(void * data) {return static_cast<GriddedField2 *>(data) -> get_dim();}
Index get_grid_typeIndexGriddedField2(Index i, void * data) {return Index(static_cast<GriddedField2 *>(data) -> get_grid_type(i));}
Index get_grid_sizeGriddedField2(Index i, void * data) {return static_cast<GriddedField2 *>(data) -> get_grid_size(i);}
char * get_nameGriddedField2(void * data) {return const_cast<char *>(static_cast<GriddedField2 *>(data) -> get_name().data());}
void set_nameGriddedField2(void * data, char * newdata) {static_cast<GriddedField2 *>(data) -> set_name(newdata);}
char * get_grid_nameGriddedField2(Index i, void * data) {return const_cast<char *>(static_cast<GriddedField2 *>(data) -> get_grid_name(i).data());}
void set_grid_nameGriddedField2(Index i, void * data, char * newdata) {static_cast<GriddedField2 *>(data) -> set_grid_name(i, newdata);}
void * get_numeric_gridGriddedField2(Index i, void * data) {return &static_cast<GriddedField2 *>(data) -> get_numeric_grid(i);}
void * get_string_gridGriddedField2(Index i, void * data) {return &static_cast<GriddedField2 *>(data) -> get_string_grid(i);}
void set_gridGriddedField2(Index i, void * data, void * newdata, bool NumericType) 
{
  if (NumericType)
    static_cast<GriddedField2 *>(data) -> set_grid(i, *static_cast<Vector *>(newdata));
  else
    static_cast<GriddedField2 *>(data) -> set_grid(i, *static_cast<ArrayOfString *>(newdata));
}
void * dataGriddedField2(void * data) {return &static_cast<GriddedField2 *>(data) -> data;}
bool checksizeGriddedField2(void * data) {return static_cast<GriddedField2 *>(data) -> checksize();}


// GriddedField3
BasicInterfaceCAPI(GriddedField3)
BasicInputOutputCAPI(GriddedField3)
VoidArrayCAPI(ArrayOfGriddedField3)
BasicInterfaceCAPI(ArrayOfGriddedField3)
BasicInputOutputCAPI(ArrayOfGriddedField3)
VoidArrayCAPI(ArrayOfArrayOfGriddedField3)
BasicInterfaceCAPI(ArrayOfArrayOfGriddedField3)
BasicInputOutputCAPI(ArrayOfArrayOfGriddedField3)
Index get_dimGriddedField3(void * data) {return static_cast<GriddedField3 *>(data) -> get_dim();}
Index get_grid_typeIndexGriddedField3(Index i, void * data) {return Index(static_cast<GriddedField3 *>(data) -> get_grid_type(i));}
Index get_grid_sizeGriddedField3(Index i, void * data) {return static_cast<GriddedField3 *>(data) -> get_grid_size(i);}
char * get_nameGriddedField3(void * data) {return const_cast<char *>(static_cast<GriddedField3 *>(data) -> get_name().data());}
void set_nameGriddedField3(void * data, char * newdata) {static_cast<GriddedField3 *>(data) -> set_name(newdata);}
char * get_grid_nameGriddedField3(Index i, void * data) {return const_cast<char *>(static_cast<GriddedField3 *>(data) -> get_grid_name(i).data());}
void set_grid_nameGriddedField3(Index i, void * data, char * newdata) {static_cast<GriddedField3 *>(data) -> set_grid_name(i, newdata);}
void * get_numeric_gridGriddedField3(Index i, void * data) {return &static_cast<GriddedField3 *>(data) -> get_numeric_grid(i);}
void * get_string_gridGriddedField3(Index i, void * data) {return &static_cast<GriddedField3 *>(data) -> get_string_grid(i);}
void set_gridGriddedField3(Index i, void * data, void * newdata, bool NumericType) 
{
  if (NumericType)
    static_cast<GriddedField3 *>(data) -> set_grid(i, *static_cast<Vector *>(newdata));
  else
    static_cast<GriddedField3 *>(data) -> set_grid(i, *static_cast<ArrayOfString *>(newdata));
}
void * dataGriddedField3(void * data) {return &static_cast<GriddedField3 *>(data) -> data;}
bool checksizeGriddedField3(void * data) {return static_cast<GriddedField3 *>(data) -> checksize();}


// GriddedField4
BasicInterfaceCAPI(GriddedField4)
BasicInputOutputCAPI(GriddedField4)
VoidArrayCAPI(ArrayOfGriddedField4)
BasicInterfaceCAPI(ArrayOfGriddedField4)
BasicInputOutputCAPI(ArrayOfGriddedField4)
Index get_dimGriddedField4(void * data) {return static_cast<GriddedField4 *>(data) -> get_dim();}
Index get_grid_typeIndexGriddedField4(Index i, void * data) {return Index(static_cast<GriddedField4 *>(data) -> get_grid_type(i));}
Index get_grid_sizeGriddedField4(Index i, void * data) {return static_cast<GriddedField4 *>(data) -> get_grid_size(i);}
char * get_nameGriddedField4(void * data) {return const_cast<char *>(static_cast<GriddedField4 *>(data) -> get_name().data());}
void set_nameGriddedField4(void * data, char * newdata) {static_cast<GriddedField4 *>(data) -> set_name(newdata);}
char * get_grid_nameGriddedField4(Index i, void * data) {return const_cast<char *>(static_cast<GriddedField4 *>(data) -> get_grid_name(i).data());}
void set_grid_nameGriddedField4(Index i, void * data, char * newdata) {static_cast<GriddedField4 *>(data) -> set_grid_name(i, newdata);}
void * get_numeric_gridGriddedField4(Index i, void * data) {return &static_cast<GriddedField4 *>(data) -> get_numeric_grid(i);}
void * get_string_gridGriddedField4(Index i, void * data) {return &static_cast<GriddedField4 *>(data) -> get_string_grid(i);}
void set_gridGriddedField4(Index i, void * data, void * newdata, bool NumericType) 
{
  if (NumericType)
    static_cast<GriddedField4 *>(data) -> set_grid(i, *static_cast<Vector *>(newdata));
  else
    static_cast<GriddedField4 *>(data) -> set_grid(i, *static_cast<ArrayOfString *>(newdata));
}
void * dataGriddedField4(void * data) {return &static_cast<GriddedField4 *>(data) -> data;}
bool checksizeGriddedField4(void * data) {return static_cast<GriddedField4 *>(data) -> checksize();}


// GriddedField5
BasicInterfaceCAPI(GriddedField5)
BasicInputOutputCAPI(GriddedField5)
Index get_dimGriddedField5(void * data) {return static_cast<GriddedField5 *>(data) -> get_dim();}
Index get_grid_typeIndexGriddedField5(Index i, void * data) {return Index(static_cast<GriddedField5 *>(data) -> get_grid_type(i));}
Index get_grid_sizeGriddedField5(Index i, void * data) {return static_cast<GriddedField5 *>(data) -> get_grid_size(i);}
char * get_nameGriddedField5(void * data) {return const_cast<char *>(static_cast<GriddedField5 *>(data) -> get_name().data());}
void set_nameGriddedField5(void * data, char * newdata) {static_cast<GriddedField5 *>(data) -> set_name(newdata);}
char * get_grid_nameGriddedField5(Index i, void * data) {return const_cast<char *>(static_cast<GriddedField5 *>(data) -> get_grid_name(i).data());}
void set_grid_nameGriddedField5(Index i, void * data, char * newdata) {static_cast<GriddedField5 *>(data) -> set_grid_name(i, newdata);}
void * get_numeric_gridGriddedField5(Index i, void * data) {return &static_cast<GriddedField5 *>(data) -> get_numeric_grid(i);}
void * get_string_gridGriddedField5(Index i, void * data) {return &static_cast<GriddedField5 *>(data) -> get_string_grid(i);}
void set_gridGriddedField5(Index i, void * data, void * newdata, bool NumericType) 
{
  if (NumericType)
    static_cast<GriddedField5 *>(data) -> set_grid(i, *static_cast<Vector *>(newdata));
  else
    static_cast<GriddedField5 *>(data) -> set_grid(i, *static_cast<ArrayOfString *>(newdata));
}
void * dataGriddedField5(void * data) {return &static_cast<GriddedField5 *>(data) -> data;}
bool checksizeGriddedField5(void * data) {return static_cast<GriddedField5 *>(data) -> checksize();}


// GriddedField6
BasicInterfaceCAPI(GriddedField6)
BasicInputOutputCAPI(GriddedField6)
Index get_dimGriddedField6(void * data) {return static_cast<GriddedField6 *>(data) -> get_dim();}
Index get_grid_typeIndexGriddedField6(Index i, void * data) {return Index(static_cast<GriddedField6 *>(data) -> get_grid_type(i));}
Index get_grid_sizeGriddedField6(Index i, void * data) {return static_cast<GriddedField6 *>(data) -> get_grid_size(i);}
char * get_nameGriddedField6(void * data) {return const_cast<char *>(static_cast<GriddedField6 *>(data) -> get_name().data());}
void set_nameGriddedField6(void * data, char * newdata) {static_cast<GriddedField6 *>(data) -> set_name(newdata);}
char * get_grid_nameGriddedField6(Index i, void * data) {return const_cast<char *>(static_cast<GriddedField6 *>(data) -> get_grid_name(i).data());}
void set_grid_nameGriddedField6(Index i, void * data, char * newdata) {static_cast<GriddedField6 *>(data) -> set_grid_name(i, newdata);}
void * get_numeric_gridGriddedField6(Index i, void * data) {return &static_cast<GriddedField6 *>(data) -> get_numeric_grid(i);}
void * get_string_gridGriddedField6(Index i, void * data) {return &static_cast<GriddedField6 *>(data) -> get_string_grid(i);}
void set_gridGriddedField6(Index i, void * data, void * newdata, bool NumericType) 
{
  if (NumericType)
    static_cast<GriddedField6 *>(data) -> set_grid(i, *static_cast<Vector *>(newdata));
  else
    static_cast<GriddedField6 *>(data) -> set_grid(i, *static_cast<ArrayOfString *>(newdata));
}
void * dataGriddedField6(void * data) {return &static_cast<GriddedField6 *>(data) -> data;}
bool checksizeGriddedField6(void * data) {return static_cast<GriddedField6 *>(data) -> checksize();}


// SpeciesAuxData
BasicInterfaceCAPI(SpeciesAuxData)
BasicInputOutputCAPI(SpeciesAuxData)
void initSpeciesAuxData(void * data) {static_cast<SpeciesAuxData *>(data) -> InitFromSpeciesData();}
bool validindexSpeciesAuxData(void * data, Index s, Index i) {return static_cast<SpeciesAuxData *>(data) -> validIndex(s, i);}
void * getDataSpeciesAuxData(void * data, Index s, Index i) {return &static_cast<SpeciesAuxData *>(data) -> Data(s, i);}
Index setTypeFromIndexSpeciesAuxData(void * data, Index s, Index i, Index t) {return static_cast<SpeciesAuxData *>(data) -> setParamType(s, i, t);}
Index getTypeSpeciesAuxData(void * data, Index s, Index i) {return Index(static_cast<SpeciesAuxData *>(data) -> getParamType(s, i));}


// CIARecord
BasicInterfaceCAPI(CIARecord)
BasicInputOutputCAPI(CIARecord)
VoidGetterCAPI(CIARecord, Data)
VoidArrayCAPI(ArrayOfCIARecord)
BasicInterfaceCAPI(ArrayOfCIARecord)
BasicInputOutputCAPI(ArrayOfCIARecord)
Index getSpecies1CIARecord(void * data) {return static_cast<CIARecord *>(data) -> Species(0);}
Index getSpecies2CIARecord(void * data) {return static_cast<CIARecord *>(data) -> Species(1);}
void setSpeciesCIARecord(void * data, Index newval1, Index newval2) {return static_cast<CIARecord *>(data) -> SetSpecies(newval1, newval2);}


// Verbosity
BasicInterfaceCAPI(Verbosity)
BasicInputOutputCAPI(Verbosity)
Index getAgendaVerbosity(void * data) {return static_cast<Verbosity *>(data) -> get_agenda_verbosity();}
Index getScreenVerbosity(void * data) {return static_cast<Verbosity *>(data) -> get_screen_verbosity();}
Index getFileVerbosity(void * data) {return static_cast<Verbosity *>(data) -> get_file_verbosity();}
bool getMainVerbosity(void * data) {return static_cast<Verbosity *>(data) -> is_main_agenda();}
void setVerbosity(void * data, Index a, Index s, Index f, bool m) {
  auto x = static_cast<Verbosity *>(data);
  x -> set_agenda_verbosity(a);
  x -> set_screen_verbosity(s);
  x -> set_file_verbosity(f);
  x -> set_main_agenda(m);
}


// TessemNN
void * createTessemNN() {return new TessemNN;}
void deleteTessemNN(void * data) {delete static_cast<TessemNN *>(data);}
void printTessemNN(void *) {std::cout << std::endl;}
BasicInputOutputCAPI(TessemNN)
VoidStructGetterCAPI(TessemNN, nb_inputs)
VoidStructGetterCAPI(TessemNN, nb_outputs)
VoidStructGetterCAPI(TessemNN, nb_cache)
VoidStructGetterCAPI(TessemNN, b1)
VoidStructGetterCAPI(TessemNN, b2)
VoidStructGetterCAPI(TessemNN, w1)
VoidStructGetterCAPI(TessemNN, w2)
VoidStructGetterCAPI(TessemNN, x_min)
VoidStructGetterCAPI(TessemNN, x_max)
VoidStructGetterCAPI(TessemNN, y_min)
VoidStructGetterCAPI(TessemNN, y_max)


// SingleScatteringData
BasicInterfaceCAPI(SingleScatteringData)
BasicInputOutputCAPI(SingleScatteringData)
VoidStructGetterCAPI(SingleScatteringData, ptype)
VoidStructGetterCAPI(SingleScatteringData, description)
VoidStructGetterCAPI(SingleScatteringData, f_grid)
VoidStructGetterCAPI(SingleScatteringData, T_grid)
VoidStructGetterCAPI(SingleScatteringData, za_grid)
VoidStructGetterCAPI(SingleScatteringData, aa_grid)
VoidStructGetterCAPI(SingleScatteringData, pha_mat_data)
VoidStructGetterCAPI(SingleScatteringData, ext_mat_data)
VoidStructGetterCAPI(SingleScatteringData, abs_vec_data)
VoidArrayCAPI(ArrayOfSingleScatteringData)
BasicInterfaceCAPI(ArrayOfSingleScatteringData)
BasicInputOutputCAPI(ArrayOfSingleScatteringData)
VoidArrayCAPI(ArrayOfArrayOfSingleScatteringData)
BasicInterfaceCAPI(ArrayOfArrayOfSingleScatteringData)
BasicInputOutputCAPI(ArrayOfArrayOfSingleScatteringData)


// ScatteringMetaData
BasicInterfaceCAPI(ScatteringMetaData)
BasicInputOutputCAPI(ScatteringMetaData)
VoidStructGetterCAPI(ScatteringMetaData, description)
VoidStructGetterCAPI(ScatteringMetaData, source)
VoidStructGetterCAPI(ScatteringMetaData, refr_index)
VoidStructGetterCAPI(ScatteringMetaData, mass)
VoidStructGetterCAPI(ScatteringMetaData, diameter_max)
VoidStructGetterCAPI(ScatteringMetaData, diameter_volume_equ)
VoidStructGetterCAPI(ScatteringMetaData, diameter_area_equ_aerodynamical)
VoidArrayCAPI(ArrayOfScatteringMetaData)
BasicInterfaceCAPI(ArrayOfScatteringMetaData)
BasicInputOutputCAPI(ArrayOfScatteringMetaData)
VoidArrayCAPI(ArrayOfArrayOfScatteringMetaData)
BasicInterfaceCAPI(ArrayOfArrayOfScatteringMetaData)
BasicInputOutputCAPI(ArrayOfArrayOfScatteringMetaData)


// Timer
void * createTimer() {return new Timer;}
void deleteTimer(void * data) {delete static_cast<Timer *>(data);}
void printTimer(void *) {std::cout << std::endl;}
BasicInputOutputCAPI(Timer)
bool getrunningTimer(void * data) {return static_cast<Timer *>(data) -> running;}
bool getfinishedTimer(void * data) {return static_cast<Timer *>(data) -> finished;}
Index getcputime_start_utimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_start.tms_utime;
#else
  return 0;
#endif
}
Index getcputime_start_stimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_start.tms_stime;
#else
  return 0;
#endif
}
Index getcputime_start_cutimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_start.tms_cutime;
#else
  return 0;
#endif
}
Index getcputime_start_cstimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_start.tms_cstime;
#else
  return 0;
#endif
}
Index getrealtime_startTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> realtime_start;
#else
  return 0;
#endif
}
Index getcputime_end_utimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_end.tms_utime;
#else
  return 0;
#endif
}
Index getcputime_end_stimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_end.tms_stime;
#else
  return 0;
#endif
}
Index getcputime_end_cutimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_end.tms_cutime;
#else
  return 0;
#endif
}
Index getcputime_end_cstimeTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> cputime_end.tms_cstime;
#else
  return 0;
#endif
}
Index getrealtime_endTimer(void * data)
{
#ifdef TIME_SUPPORT
  return static_cast<Timer *>(data) -> realtime_end;
#else
  return 0;
#endif
}
void setrunningTimer(void * data, bool newdata) {static_cast<Timer *>(data) -> running = newdata;}
void setfinishedTimer(void * data, bool newdata) {static_cast<Timer *>(data) -> finished = newdata;}
void setcputime_start_utimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_start.tms_utime = clock_t(newdata);
#endif
}
void setcputime_start_stimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_start.tms_stime = clock_t(newdata);
#endif
}
void setcputime_start_cutimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_start.tms_cutime = clock_t(newdata);
#endif
}
void setcputime_start_cstimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_start.tms_cstime = clock_t(newdata);
#endif
}
void setrealtime_startTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> realtime_start = clock_t(newdata);
#endif
}
void setcputime_end_utimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_end.tms_utime = clock_t(newdata);
#endif
}
void setcputime_end_stimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_end.tms_stime = clock_t(newdata);
#endif
}
void setcputime_end_cutimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_end.tms_cutime = clock_t(newdata);
#endif
}
void setcputime_end_cstimeTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> cputime_end.tms_cstime = clock_t(newdata);
#endif
}
void setrealtime_endTimer(void * data, Index newdata)
{
#ifdef TIME_SUPPORT
  static_cast<Timer *>(data) -> realtime_end = clock_t(newdata);
#endif
}
bool supportTimer()
{
#ifdef TIME_SUPPORT
  return true;
#else
  return false;
#endif
}
Index tickTimer()
{
#ifdef TIME_SUPPORT
  return sysconf(_SC_CLK_TCK);
#else
  return 0;
#endif
}


// TelsemAtlas
BasicInterfaceCAPI(TelsemAtlas)
BasicInputOutputCAPI(TelsemAtlas)
VoidGetterCAPI(TelsemAtlas, DataCount)
VoidGetterCAPI(TelsemAtlas, ChannelCount)
VoidGetterCAPI(TelsemAtlas, Name)
VoidGetterCAPI(TelsemAtlas, Month)
VoidGetterCAPI(TelsemAtlas, Lat)
VoidGetterCAPI(TelsemAtlas, Cells)
VoidGetterCAPI(TelsemAtlas, FirstCells)
VoidGetterCAPI(TelsemAtlas, Emis)
VoidGetterCAPI(TelsemAtlas, Emis_err)
VoidGetterCAPI(TelsemAtlas, Correlations)
VoidGetterCAPI(TelsemAtlas, Classes1)
VoidGetterCAPI(TelsemAtlas, Classes2)
VoidGetterCAPI(TelsemAtlas, Cellnumber)
VoidGetterCAPI(TelsemAtlas, Correspondance)
VoidArrayCAPI(ArrayOfTelsemAtlas)
BasicInterfaceCAPI(ArrayOfTelsemAtlas)
BasicInputOutputCAPI(ArrayOfTelsemAtlas)
Numeric getA0_K0TelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A0_K0(i);}
Numeric getA0_K1TelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A0_K1(i);}
Numeric getA0_K2TelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A0_K2(i);}
Numeric getA0_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A0_EVEH(i);}
Numeric getA1_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A1_EVEH(i);}
Numeric getA2_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A2_EVEH(i);}
Numeric getA3_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> A3_EVEH(i);}
Numeric getB0_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> B0_EVEH(i);}
Numeric getB1_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> B1_EVEH(i);}
Numeric getB2_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> B2_EVEH(i);}
Numeric getB3_EVEHTelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> B3_EVEH(i);}
Numeric getRAPPORT43_32TelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> RAPPORT43_32(i);}
Numeric getRAPPORT54_43TelsemAtlas(Index i, void* data) {return static_cast<TelsemAtlas *>(data) -> RAPPORT54_43(i);}


// MCAntenna
BasicInterfaceCAPI(MCAntenna)
BasicInputOutputCAPI(MCAntenna)
VoidGetterCAPI(MCAntenna, saa)
VoidGetterCAPI(MCAntenna, sza)
VoidGetterCAPI(MCAntenna, aag)
VoidGetterCAPI(MCAntenna, zag)
VoidGetterCAPI(MCAntenna, G)
Index getTypeMCAntenna(void * data) {return Index(static_cast<MCAntenna *>(data) -> Type());}
Index setTypeMCAntenna(void * data, Index newval) {return Index(static_cast<MCAntenna *>(data) -> Type(AntennaType(newval)));}


// GasAbsLookup
BasicInterfaceCAPI(GasAbsLookup)
BasicInputOutputCAPI(GasAbsLookup)
VoidGetterCAPI(GasAbsLookup, Species)
VoidGetterCAPI(GasAbsLookup, NonLinearSpecies)
VoidGetterCAPI(GasAbsLookup, Fgrid)
VoidGetterCAPI(GasAbsLookup, FGPDefault)
VoidGetterCAPI(GasAbsLookup, Pgrid)
VoidGetterCAPI(GasAbsLookup, LogPgrid)
VoidGetterCAPI(GasAbsLookup, VMRs)
VoidGetterCAPI(GasAbsLookup, Tref)
VoidGetterCAPI(GasAbsLookup, Tpert)
VoidGetterCAPI(GasAbsLookup, NLSPert)
VoidGetterCAPI(GasAbsLookup, Xsec)


// XsecRecord
BasicInterfaceCAPI(XsecRecord)
BasicInputOutputCAPI(XsecRecord)
VoidGetterCAPI(XsecRecord, Coeffs)
VoidGetterCAPI(XsecRecord, RefPressure)
VoidGetterCAPI(XsecRecord, RefTemperature)
VoidGetterCAPI(XsecRecord, Fgrids)
VoidGetterCAPI(XsecRecord, Xsecs)
VoidGetterCAPI(XsecRecord, TemperatureSlope)
VoidGetterCAPI(XsecRecord, TemperatureIntersect)
VoidArrayCAPI(ArrayOfXsecRecord)
BasicInterfaceCAPI(ArrayOfXsecRecord)
BasicInputOutputCAPI(ArrayOfXsecRecord)
Index getSpeciesXsecRecord(void * data) {return static_cast<XsecRecord *>(data) -> Species();}
void setSpeciesXsecRecord(void * data, Index newdata) {static_cast<XsecRecord *>(data) -> SetSpecies(newdata);}


// Sparse
BasicInterfaceCAPI(Sparse)
BasicInputOutputCAPI(Sparse)
VoidArrayCAPI(ArrayOfSparse)
BasicInterfaceCAPI(ArrayOfSparse)
BasicInputOutputCAPI(ArrayOfSparse)
void resizeSparse(Index nrows, Index ncols, void * data) {static_cast<Sparse *>(data) -> resize(nrows, ncols);}
Index rowsSparse(void * data) {return static_cast<Sparse *>(data) -> nrows();}
Index colsSparse(void * data) {return static_cast<Sparse *>(data) -> ncols();}
Index sizeSparse(void * data) {return static_cast<Sparse *>(data) -> nnz();}
int * rowsptrSparse(void * data) {return static_cast<Sparse *>(data) -> get_row_start_pointer();}
int * colsptrSparse(void * data) {return static_cast<Sparse *>(data) -> get_column_index_pointer();}
Numeric * getDataSparse(void * data) {return static_cast<Sparse *>(data) -> get_element_pointer();}
void setDataSparse(void * data, Index r, Index c, Numeric v) {(static_cast<Sparse *>(data) -> rw(r, c)) = v;}


// CovarianceMatrix
BasicInterfaceCAPI(CovarianceMatrix)
BasicInputOutputCAPI(CovarianceMatrix)
Index sizeget_blocksCovarianceMatrix(void * data) { return static_cast<CovarianceMatrix *>(data) -> get_blocks().size(); }
void resizeget_blocksCovarianceMatrix(Index n, void * data) { static_cast<CovarianceMatrix *>(data) -> get_blocks() = std::vector<Block>(n, Block(Range(joker), Range(joker), {0, 0}, std::make_shared<Matrix>(Matrix()))); }
void * getelemget_blocksCovarianceMatrix(Index i, void * data) { return &static_cast<CovarianceMatrix *>(data) -> get_blocks()[i]; }
Index sizeget_inverse_blocksCovarianceMatrix(void * data) { return static_cast<CovarianceMatrix *>(data) -> get_inverse_blocks().size(); }
void resizeget_inverse_blocksCovarianceMatrix(Index n, void * data) { static_cast<CovarianceMatrix *>(data) -> get_inverse_blocks() = std::vector<Block>(n, Block(Range(joker), Range(joker), {0, 0}, std::make_shared<Matrix>(Matrix()))); }
void * getelemget_inverse_blocksCovarianceMatrix(Index i, void * data) { return &static_cast<CovarianceMatrix *>(data) -> get_inverse_blocks()[i]; }


// Any
void * createAny() {return new Any;}
void deleteAny(void * data) {delete static_cast<Any *>(data);}
void printAny(void *) {std::cout << std::endl;}
Index xmlreadAny(void *, char *) {return 1;}
Index xmlsaveAny(void *, char *, Index, Index) {return 1;}


// Agenda
BasicInterfaceCAPI(Agenda)
BasicInputOutputCAPI(Agenda)
VoidArrayCAPI(ArrayOfAgenda)
BasicInterfaceCAPI(ArrayOfAgenda)
BasicInputOutputCAPI(ArrayOfAgenda)


// RetrievalQuantity
BasicInterfaceCAPI(RetrievalQuantity)
BasicInputOutputCAPI(RetrievalQuantity)
VoidGetterCAPI(RetrievalQuantity, MainTag)
VoidGetterCAPI(RetrievalQuantity, SubTag)
VoidGetterCAPI(RetrievalQuantity, SubSubTag)
VoidGetterCAPI(RetrievalQuantity, Mode)
VoidGetterCAPI(RetrievalQuantity, Analytical)
VoidGetterCAPI(RetrievalQuantity, Perturbation)
VoidGetterCAPI(RetrievalQuantity, Grids)
VoidGetterCAPI(RetrievalQuantity, QuantumIdentity)
VoidGetterCAPI(RetrievalQuantity, TransformationFunc)
VoidGetterCAPI(RetrievalQuantity, TFuncParameters)
VoidGetterCAPI(RetrievalQuantity, Transformation)
VoidGetterCAPI(RetrievalQuantity, Offset)
GetterSetterCAPI(RetrievalQuantity, Integration, bool)
VoidArrayCAPI(ArrayOfRetrievalQuantity)
BasicInterfaceCAPI(ArrayOfRetrievalQuantity)
BasicInputOutputCAPI(ArrayOfRetrievalQuantity)
Index getTypeRetrievalQuantity(void * data) {return Index(static_cast<RetrievalQuantity *>(data) -> Proptype());}
Index setTypeRetrievalQuantity(void * data, Index newval) {return Index(static_cast<RetrievalQuantity *>(data) -> Proptype(JacPropMatType(newval)));}


// Range
void * createRange() {return new Range(joker);}
void deleteRange(void * data) {delete static_cast<Range *>(data);}
void printRange(void * data) {std::cout << (*static_cast<Range *>(data)) << std::endl;}
Index get_startRange(void * data) {return static_cast<Range *>(data) -> get_start();}
Index get_strideRange(void * data) {return static_cast<Range *>(data) -> get_stride();}
Index get_extentRange(void * data) {return static_cast<Range *>(data) -> get_extent();}
void setRange(void * data, Index start, Index extent, Index stride)
{
  if (extent >= 0)
    static_cast<Range *>(data) -> operator=(Range(start, extent, stride));
  else
    static_cast<Range *>(data) -> operator=(Range(start, joker, stride));
}


// Block
void * createBlock() {return new Block(Range(joker), Range(joker), {0, 0}, std::make_shared<Matrix>(Matrix()));}
void deleteBlock(void * data) {delete static_cast<Block *>(data);}
void printBlock(void *) {std::cout << std::endl;}
VoidGetterCAPI(Block, get_row_range)
VoidGetterCAPI(Block, get_column_range)
VoidGetterCAPI(Block, get_dense)
VoidGetterCAPI(Block, get_sparse)
Index get_matrix_typeBlock(void * data) {return Index(static_cast<Block *>(data) -> get_matrix_type());}
Index get_index1Block(void * data) {return static_cast<Block *>(data) -> get_indices().first;}
Index get_index2Block(void * data) {return static_cast<Block *>(data) -> get_indices().second;}
void set_indicesBlock(void * data, Index i1, Index i2) {static_cast<Block *>(data) -> set_indices(i1, i2);}
void set_matrixBlock(void * data, void * newdata, bool dense)
{
  auto x = static_cast<Block *>(data);
  if (dense) {
    x -> operator=(Block(x -> get_row_range(), x -> get_column_range(), x -> get_indices(), std::make_shared<Matrix>(*static_cast<Matrix *>(newdata))));
  } else {
    x -> operator=(Block(x -> get_row_range(), x -> get_column_range(), x -> get_indices(), std::make_shared<Sparse>(*static_cast<Sparse *>(newdata))));
  }
}


// Time
BasicInterfaceCAPI(Time)
BasicInputOutputCAPI(Time)
GetterSetterCAPI(Time, Seconds, Numeric)
VoidArrayCAPI(ArrayOfTime)
BasicInterfaceCAPI(ArrayOfTime)
BasicInputOutputCAPI(ArrayOfTime)


// generic
Index string2filetypeindex(char * data) { try { return Index(string2filetype(data)); } catch (std::runtime_error& e) { return -1; } }


#undef BasicInterfaceCAPI
#undef GetterSetterCAPI
#undef EnumGetterSetterCAPI
#undef VoidGetterCAPI
#undef VoidStructGetterCAPI
#undef BasicInputOutputCAPI
#undef VoidArrayCAPI
#undef VoidArrayElemCAPI

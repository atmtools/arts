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
BasicInputOutputCAPI(GridPos)
VoidArrayCAPI(ArrayOfGridPos)
BasicInterfaceCAPI(ArrayOfGridPos)
BasicInputOutputCAPI(ArrayOfGridPos)
void printGridPos(void * data) {std::cout << (*static_cast<GridPos *>(data)) << std::endl;}


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

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
  \file   arts_api_classes.h
  \author Richard Larsson <larsson@mps.mpg.de>
  \date   2020-03-12

  \brief This file contains all declarations of the ARTS C API class interface.
*/

#ifndef _ARTS_ARTS_API_CLASS_H_
#define _ARTS_ARTS_API_CLASS_H_

#include "matpack.h"

#ifndef DLL_PUBLIC
#define DLL_PUBLIC __attribute__((visibility("default")))
#define REMOVE_DLL_PUBLIC 1
#else
#define REMOVE_DLL_PUBLIC 0
#endif

#define BasicInterfaceCAPI(TYPE)        \
__attribute__((visibility("default")))  \
void * create##TYPE();                  \
__attribute__((visibility("default")))  \
void delete##TYPE(void * data);         \
__attribute__((visibility("default")))  \
void print##TYPE(void * data);


#define GetterSetterCAPI(TYPE, VALUE, BASETYPE) \
__attribute__((visibility("default")))          \
BASETYPE get##VALUE##TYPE(void * data);         \
__attribute__((visibility("default")))          \
void set##VALUE##TYPE(void * data, BASETYPE newval);


#define EnumGetterSetterCAPI(TYPE, VALUE, ENUM)     \
__attribute__((visibility("default")))              \
Index get##VALUE##TYPE(void * data);                \
__attribute__((visibility("default")))              \
Index set##VALUE##TYPE(void * data, Index newval);  \
__attribute__((visibility("default")))              \
Index string2index##VALUE##TYPE(void * data, char * newval);


#define VoidGetterCAPI(TYPE, VALUE)     \
__attribute__((visibility("default")))  \
void * get##VALUE##TYPE(void * data);


#define VoidStructGetterCAPI(TYPE, VALUE) \
__attribute__((visibility("default")))    \
void * get##VALUE##TYPE(void * data);


#define BasicInputOutputCAPI(TYPE)                  \
__attribute__((visibility("default")))              \
Index xmlread##TYPE(void * data, char * filepath);  \
__attribute__((visibility("default")))              \
Index xmlsave##TYPE(void * data, char * filepath, Index filetype, Index clobber);


#define VoidArrayCAPI(TYPE)                 \
__attribute__((visibility("default")))      \
Index size##TYPE(void * data);              \
__attribute__((visibility("default")))      \
void resize##TYPE(Index n, void * data);    \
__attribute__((visibility("default")))      \
void * getelem##TYPE(Index i, void * data);


#define VoidArrayElemCAPI(TYPE, ELEM)           \
__attribute__((visibility("default")))          \
Index size##ELEM##TYPE(void * data);            \
__attribute__((visibility("default")))          \
void resize##ELEM##TYPE(Index n, void * data);  \
__attribute__((visibility("default")))          \
void * getelem##ELEM##TYPE(Index i, void * data);


extern "C" {
    // Index
    BasicInterfaceCAPI(Index)
    BasicInputOutputCAPI(Index)
    VoidArrayCAPI(ArrayOfIndex)
    BasicInterfaceCAPI(ArrayOfIndex)
    BasicInputOutputCAPI(ArrayOfIndex)
    VoidArrayCAPI(ArrayOfArrayOfIndex)
    BasicInterfaceCAPI(ArrayOfArrayOfIndex)
    BasicInputOutputCAPI(ArrayOfArrayOfIndex)
    DLL_PUBLIC Index getIndex(void *);
    DLL_PUBLIC void setIndex(void *, Index);
    
    // Numeric
    BasicInterfaceCAPI(Numeric)
    BasicInputOutputCAPI(Numeric)
    DLL_PUBLIC Numeric getNumeric(void * data);
    DLL_PUBLIC void setNumeric(void * data, Numeric newval);
  
    // ZeemanModel
    BasicInterfaceCAPI(ZeemanModel)
    GetterSetterCAPI(ZeemanModel, gu, Numeric)
    GetterSetterCAPI(ZeemanModel, gl, Numeric)

    // Rational
    BasicInterfaceCAPI(Rational)
    BasicInputOutputCAPI(Rational)
    GetterSetterCAPI(Rational, Nom, Index)
    GetterSetterCAPI(Rational, Denom, Index)
    DLL_PUBLIC void simplifyRational(void *);
  
    // LineShape::ModelParameters
    DLL_PUBLIC void printLineShapeModelParameters(void *);
    DLL_PUBLIC Index getLineShapeModelParametersType(char *);
  
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
    DLL_PUBLIC void * getelemQuantumNumbers(Index, void *);
    DLL_PUBLIC Index sizeQuantumNumbers();
    DLL_PUBLIC Index string2quantumnumbersindex(char *);
  
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
    DLL_PUBLIC Index setSpeciesTag(void *, char *);
    DLL_PUBLIC Index validSpecies(Index);
    DLL_PUBLIC Index validAllIsotopologues(Index, Index);
    DLL_PUBLIC Index validIsotopologue(Index, Index);
    DLL_PUBLIC Index validContinuum(Index, Index);
  
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
    DLL_PUBLIC void printmetaAbsorptionLines(void *);
    DLL_PUBLIC Index isAbsorptionLinesOK(void *);
    
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
    DLL_PUBLIC void resizeVector(Index, void *);
    DLL_PUBLIC Index nelemVector(void *);
    DLL_PUBLIC Numeric * getDataVector(void *);
    
    // Matrix
    BasicInterfaceCAPI(Matrix)
    BasicInputOutputCAPI(Matrix)
    VoidArrayCAPI(ArrayOfMatrix)
    BasicInterfaceCAPI(ArrayOfMatrix)
    BasicInputOutputCAPI(ArrayOfMatrix)
    VoidArrayCAPI(ArrayOfArrayOfMatrix)
    BasicInterfaceCAPI(ArrayOfArrayOfMatrix)
    BasicInputOutputCAPI(ArrayOfArrayOfMatrix)
    DLL_PUBLIC void resizeMatrix(Index, Index, void *);
    DLL_PUBLIC Index rowsMatrix(void *);
    DLL_PUBLIC Index colsMatrix(void *);
    DLL_PUBLIC Numeric * getDataMatrix(void *);
    
    // Tensor3
    BasicInterfaceCAPI(Tensor3)
    BasicInputOutputCAPI(Tensor3)
    VoidArrayCAPI(ArrayOfTensor3)
    BasicInterfaceCAPI(ArrayOfTensor3)
    BasicInputOutputCAPI(ArrayOfTensor3)
    VoidArrayCAPI(ArrayOfArrayOfTensor3)
    BasicInterfaceCAPI(ArrayOfArrayOfTensor3)
    BasicInputOutputCAPI(ArrayOfArrayOfTensor3)
    DLL_PUBLIC void resizeTensor3(Index, Index, Index, void *);
    DLL_PUBLIC Index pagesTensor3(void *);
    DLL_PUBLIC Index rowsTensor3(void *);
    DLL_PUBLIC Index colsTensor3(void *);
    DLL_PUBLIC Numeric * getDataTensor3(void *);
    
    // Tensor4
    BasicInterfaceCAPI(Tensor4)
    BasicInputOutputCAPI(Tensor4)
    VoidArrayCAPI(ArrayOfTensor4)
    BasicInterfaceCAPI(ArrayOfTensor4)
    BasicInputOutputCAPI(ArrayOfTensor4)
//     VoidArrayCAPI(ArrayOfArrayOfTensor4)
//     BasicInterfaceCAPI(ArrayOfArrayOfTensor4)
//     BasicInputOutputCAPI(ArrayOfArrayOfTensor4)
    DLL_PUBLIC void resizeTensor4(Index, Index, Index, Index, void *);
    DLL_PUBLIC Index booksTensor4(void *);
    DLL_PUBLIC Index pagesTensor4(void *);
    DLL_PUBLIC Index rowsTensor4(void *);
    DLL_PUBLIC Index colsTensor4(void *);
    DLL_PUBLIC Numeric * getDataTensor4(void *);
    
    // Tensor5
    BasicInterfaceCAPI(Tensor5)
    BasicInputOutputCAPI(Tensor5)
    VoidArrayCAPI(ArrayOfTensor5)
    BasicInterfaceCAPI(ArrayOfTensor5)
    BasicInputOutputCAPI(ArrayOfTensor5)
//     VoidArrayCAPI(ArrayOfArrayOfTensor5)
//     BasicInterfaceCAPI(ArrayOfArrayOfTensor5)
//     BasicInputOutputCAPI(ArrayOfArrayOfTensor5)
    DLL_PUBLIC void resizeTensor5(Index, Index, Index, Index, Index, void *);
    DLL_PUBLIC Index shelvesTensor5(void *);
    DLL_PUBLIC Index booksTensor5(void *);
    DLL_PUBLIC Index pagesTensor5(void *);
    DLL_PUBLIC Index rowsTensor5(void *);
    DLL_PUBLIC Index colsTensor5(void *);
    DLL_PUBLIC Numeric * getDataTensor5(void *);
    
    // Tensor6
    BasicInterfaceCAPI(Tensor6)
    BasicInputOutputCAPI(Tensor6)
    VoidArrayCAPI(ArrayOfTensor6)
    BasicInterfaceCAPI(ArrayOfTensor6)
    BasicInputOutputCAPI(ArrayOfTensor6)
    VoidArrayCAPI(ArrayOfArrayOfTensor6)
    BasicInterfaceCAPI(ArrayOfArrayOfTensor6)
    BasicInputOutputCAPI(ArrayOfArrayOfTensor6)
    DLL_PUBLIC void resizeTensor6(Index, Index, Index, Index, Index, Index, void *);
    DLL_PUBLIC Index vitrinesTensor6(void *);
    DLL_PUBLIC Index shelvesTensor6(void *);
    DLL_PUBLIC Index booksTensor6(void *);
    DLL_PUBLIC Index pagesTensor6(void *);
    DLL_PUBLIC Index rowsTensor6(void *);
    DLL_PUBLIC Index colsTensor6(void *);
    DLL_PUBLIC Numeric * getDataTensor6(void *);
    
    // Tensor7
    BasicInterfaceCAPI(Tensor7)
    BasicInputOutputCAPI(Tensor7)
    VoidArrayCAPI(ArrayOfTensor7)
    BasicInterfaceCAPI(ArrayOfTensor7)
    BasicInputOutputCAPI(ArrayOfTensor7)
//     VoidArrayCAPI(ArrayOfArrayOfTensor7)
//     BasicInterfaceCAPI(ArrayOfArrayOfTensor7)
//     BasicInputOutputCAPI(ArrayOfArrayOfTensor7)
    DLL_PUBLIC void resizeTensor7(Index, Index, Index, Index, Index, Index, Index, void *);
    DLL_PUBLIC Index librariesTensor7(void *);
    DLL_PUBLIC Index vitrinesTensor7(void *);
    DLL_PUBLIC Index shelvesTensor7(void *);
    DLL_PUBLIC Index booksTensor7(void *);
    DLL_PUBLIC Index pagesTensor7(void *);
    DLL_PUBLIC Index rowsTensor7(void *);
    DLL_PUBLIC Index colsTensor7(void *);
    DLL_PUBLIC Numeric * getDataTensor7(void *);
    
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
    DLL_PUBLIC Index stokesPropagationMatrix(void *);
    DLL_PUBLIC Index frequenciesPropagationMatrix(void *);
    DLL_PUBLIC Index zenithsPropagationMatrix(void *);
    DLL_PUBLIC Index azimuthsPropagationMatrix(void *);
    DLL_PUBLIC Index setPropagationMatrix(void *, Index, Index, Index, Index, Numeric);
    
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
    DLL_PUBLIC Index stokesStokesVector(void *);
    DLL_PUBLIC Index frequenciesStokesVector(void *);
    DLL_PUBLIC Index zenithsStokesVector(void *);
    DLL_PUBLIC Index azimuthsStokesVector(void *);
    DLL_PUBLIC Index setStokesVector(void *, Index, Index, Index, Index, Numeric);
    
    // String
    BasicInterfaceCAPI(String)
    BasicInputOutputCAPI(String)
    VoidArrayCAPI(ArrayOfString)
    BasicInterfaceCAPI(ArrayOfString)
    BasicInputOutputCAPI(ArrayOfString)
    VoidArrayCAPI(ArrayOfArrayOfString)
    BasicInterfaceCAPI(ArrayOfArrayOfString)
    BasicInputOutputCAPI(ArrayOfArrayOfString)
    DLL_PUBLIC void setString(void * data, char * newdata);
    DLL_PUBLIC char * getString(void * data);
    
    // GridPos
    BasicInputOutputCAPI(GridPos)
    VoidArrayCAPI(ArrayOfGridPos)
    BasicInterfaceCAPI(ArrayOfGridPos)
    BasicInputOutputCAPI(ArrayOfGridPos)
    DLL_PUBLIC void printGridPos(void *);
    
    // Ppath
    BasicInterfaceCAPI(Ppath)
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
    BasicInterfaceCAPI(ArrayOfPpath)
    BasicInputOutputCAPI(ArrayOfPpath)
  
    // generic
    DLL_PUBLIC Index string2filetypeindex(char *);
}


#undef BasicInterfaceCAPI
#undef GetterSetterCAPI
#undef EnumGetterSetterCAPI
#undef VoidGetterCAPI
#undef VoidStructGetterCAPI
#undef BasicInputOutputCAPI
#undef VoidArrayCAPI
#undef VoidArrayElemCAPI


#if REMOVE_DLL_PUBLIC
#undef DLL_PUBLIC
#endif
#undef REMOVE_DLL_PUBLIC

#endif // _ARTS_ARTS_API_CLASS_H_


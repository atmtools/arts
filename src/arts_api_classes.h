/* Copyright (C) 2020 Richard Larsson <ric.larsson@gmail.com>

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
  \author Richard Larsson <ric.larsson@gmail.com>
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


#define StringEnumGetterSetterCAPI(TYPE, ELEM)          \
__attribute__((visibility("default")))                  \
bool enumset ## ELEM ## TYPE (void * data, char * in);  \
__attribute__((visibility("default")))                  \
void * enumget ## ELEM ## TYPE (void * data);


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

#define StringEnumPointersCAPI(ENUM)            \
__attribute__((visibility("default")))          \
void * get##ENUM##String(void *);               \
__attribute__((visibility("default")))          \
int set##ENUM##String(void *, char *);

#define EnumMacroInterfaceCAPI(TYPE)                                      \
  __attribute__((visibility("default"))) void *get##TYPE(void *data);    \
  __attribute__((visibility("default"))) int set##TYPE(void *data, char *val);

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
    VoidArrayCAPI(ArrayOfNumeric)
    BasicInterfaceCAPI(ArrayOfNumeric)
    BasicInputOutputCAPI(ArrayOfNumeric)
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
    
    // LineShapeTemperatureModel
    BasicInterfaceCAPI(LineShapeTemperatureModel)
    StringEnumPointersCAPI(LineShapeTemperatureModel)
  
    // LineShape::ModelParameters
    BasicInterfaceCAPI(LineShapeModelParameters)
    VoidStructGetterCAPI(LineShapeModelParameters, type)
    VoidStructGetterCAPI(LineShapeModelParameters, X0)
    VoidStructGetterCAPI(LineShapeModelParameters, X1)
    VoidStructGetterCAPI(LineShapeModelParameters, X2)
    VoidStructGetterCAPI(LineShapeModelParameters, X3)
  
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
    
    // LineShapeType
    BasicInterfaceCAPI(LineShapeType)
    StringEnumPointersCAPI(LineShapeType)
  
    // LineShape::Model
    BasicInterfaceCAPI(LineShapeModel)
    VoidArrayCAPI(LineShapeModel)
  
    // Absorption::SingleLine
    BasicInterfaceCAPI(AbsorptionSingleLine)
    VoidStructGetterCAPI(AbsorptionSingleLine, F0)
    VoidStructGetterCAPI(AbsorptionSingleLine, I0)
    VoidStructGetterCAPI(AbsorptionSingleLine, E0)
    VoidStructGetterCAPI(AbsorptionSingleLine, glow)
    VoidStructGetterCAPI(AbsorptionSingleLine, gupp)
    VoidStructGetterCAPI(AbsorptionSingleLine, A)
    VoidStructGetterCAPI(AbsorptionSingleLine, zeeman)
    VoidStructGetterCAPI(AbsorptionSingleLine, lineshape)
    VoidStructGetterCAPI(AbsorptionSingleLine, localquanta)
    VoidArrayCAPI(ArrayOfAbsorptionSingleLine)
    BasicInputOutputCAPI(ArrayOfAbsorptionSingleLine)
    BasicInterfaceCAPI(ArrayOfAbsorptionSingleLine)
    
    // QuantumNumbers
    BasicInterfaceCAPI(QuantumNumberType)
    EnumMacroInterfaceCAPI(QuantumNumberType)

    // QuantumNumberValue
    BasicInterfaceCAPI(QuantumNumberValue)
    VoidStructGetterCAPI(QuantumNumberValue, type)
    DLL_PUBLIC int setqnQuantumNumberValue(void * data, char * val);
    DLL_PUBLIC void * getuppQuantumNumberValue(void * data, bool str);
    DLL_PUBLIC void * getlowQuantumNumberValue(void * data, bool str);

    // QuantumNumberValueList
    BasicInterfaceCAPI(QuantumNumberValueList)
    DLL_PUBLIC void * getQuantumNumberValueList(void *, char *);
    DLL_PUBLIC int setQuantumNumberValueList(void *, char *, char *);
    DLL_PUBLIC void * getQuantumNumberValueListString(void * data);
    DLL_PUBLIC int setQuantumNumberValueListString(void * data, char * value);

    // QuantumNumberLocalState
    BasicInterfaceCAPI(QuantumNumberLocalState)
    VoidStructGetterCAPI(QuantumNumberLocalState, val)
    
    // Species::Species
    DLL_PUBLIC void * createSpecies();
    DLL_PUBLIC void deleteSpecies(void *);
    DLL_PUBLIC void printSpecies(void *);
    DLL_PUBLIC int setSpeciesLongName(void *, char *);
    DLL_PUBLIC int setSpeciesShortName(void *, char *);
    DLL_PUBLIC void * getSpeciesLongName(void *);
    DLL_PUBLIC void * getSpeciesShortName(void *);
    VoidArrayCAPI(ArrayOfSpecies)
    BasicInputOutputCAPI(ArrayOfSpecies)
    BasicInterfaceCAPI(ArrayOfSpecies)
    
    // SpeciesIsotopeRecord
    BasicInterfaceCAPI(SpeciesIsotopeRecord)
    DLL_PUBLIC Index getIndexSpeciesIsotopeRecordFromNames(char * spec, char * isot);
    DLL_PUBLIC Index getIndexSpeciesIsotopeRecordFromData(void *);
    DLL_PUBLIC int setSpeciesIsotopeRecordToIndex(void *, Index);
    DLL_PUBLIC void * getSpeciesSpeciesIsotopeRecord(void *);
    DLL_PUBLIC void * getIsotnameSpeciesIsotopeRecord(void *);
    DLL_PUBLIC Numeric getMassSpeciesIsotopeRecord(void *);
    DLL_PUBLIC Index getGSpeciesIsotopeRecord(void *);
    DLL_PUBLIC Index nelemSpeciesIsotopeRecordDefined();
    
    // SpeciesIsotopologueRatios
    BasicInterfaceCAPI(SpeciesIsotopologueRatios)
    BasicInputOutputCAPI(SpeciesIsotopologueRatios)
    DLL_PUBLIC Numeric * getdataSpeciesIsotopologueRatios(void *);
  
    // QuantumIdentifier
    BasicInterfaceCAPI(QuantumIdentifier)
    BasicInputOutputCAPI(QuantumIdentifier)
    VoidStructGetterCAPI(QuantumIdentifier, isotopologue_index)
    VoidStructGetterCAPI(QuantumIdentifier, val)
    DLL_PUBLIC Index fromstringQuantumIdentifier(void *, char *);
    
    // ArrayOfQuantumIdentifier
    BasicInterfaceCAPI(ArrayOfQuantumIdentifier)
    BasicInputOutputCAPI(ArrayOfQuantumIdentifier)
    VoidArrayCAPI(ArrayOfQuantumIdentifier)
  
    // SpeciesTagType
    BasicInterfaceCAPI(SpeciesTagType)
    StringEnumPointersCAPI(SpeciesTagType)
    
    // SpeciesTag
    BasicInterfaceCAPI(SpeciesTag)
    BasicInputOutputCAPI(SpeciesTag)
    VoidStructGetterCAPI(SpeciesTag, spec_ind)
    VoidStructGetterCAPI(SpeciesTag, lower_freq)
    VoidStructGetterCAPI(SpeciesTag, upper_freq)
    VoidStructGetterCAPI(SpeciesTag, type)
    VoidStructGetterCAPI(SpeciesTag, cia_2nd_species)
    VoidStructGetterCAPI(SpeciesTag, cia_dataset_index)
    VoidArrayCAPI(ArrayOfSpeciesTag)
    BasicInterfaceCAPI(ArrayOfSpeciesTag)
    BasicInputOutputCAPI(ArrayOfSpeciesTag)
    VoidArrayCAPI(ArrayOfArrayOfSpeciesTag)
    BasicInterfaceCAPI(ArrayOfArrayOfSpeciesTag)
    BasicInputOutputCAPI(ArrayOfArrayOfSpeciesTag)
    DLL_PUBLIC void * getNameSpeciesTag(void *);
    DLL_PUBLIC Index setSpeciesTag(void *, char *);
    
    // AbsorptionNormalizationType
    BasicInterfaceCAPI(AbsorptionNormalizationType)
    StringEnumPointersCAPI(AbsorptionNormalizationType)
    
    // AbsorptionPopulationType
    BasicInterfaceCAPI(AbsorptionPopulationType)
    StringEnumPointersCAPI(AbsorptionPopulationType)
    
    // AbsorptionMirroringType
    BasicInterfaceCAPI(AbsorptionMirroringType)
    StringEnumPointersCAPI(AbsorptionMirroringType)
    
    // AbsorptionCutoffType
    BasicInterfaceCAPI(AbsorptionCutoffType)
    StringEnumPointersCAPI(AbsorptionCutoffType)
    
    // AbsorptionLines
    BasicInterfaceCAPI(AbsorptionLines)
    BasicInputOutputCAPI(AbsorptionLines)
    VoidStructGetterCAPI(AbsorptionLines, selfbroadening)
    VoidStructGetterCAPI(AbsorptionLines, bathbroadening)
    VoidStructGetterCAPI(AbsorptionLines, cutoff)
    VoidStructGetterCAPI(AbsorptionLines, lineshapetype)
    VoidStructGetterCAPI(AbsorptionLines, mirroring)
    VoidStructGetterCAPI(AbsorptionLines, population)
    VoidStructGetterCAPI(AbsorptionLines, normalization)
    VoidStructGetterCAPI(AbsorptionLines, T0)
    VoidStructGetterCAPI(AbsorptionLines, cutofffreq)
    VoidStructGetterCAPI(AbsorptionLines, linemixinglimit)
    VoidStructGetterCAPI(AbsorptionLines, quantumidentity)
    VoidStructGetterCAPI(AbsorptionLines, broadeningspecies)
    VoidStructGetterCAPI(AbsorptionLines, lines)
    VoidArrayCAPI(ArrayOfAbsorptionLines)
    BasicInterfaceCAPI(ArrayOfAbsorptionLines)
    BasicInputOutputCAPI(ArrayOfAbsorptionLines)
    VoidArrayCAPI(ArrayOfArrayOfAbsorptionLines)
    BasicInterfaceCAPI(ArrayOfArrayOfAbsorptionLines)
    BasicInputOutputCAPI(ArrayOfArrayOfAbsorptionLines)
    DLL_PUBLIC void printmetaAbsorptionLines(void *);
    DLL_PUBLIC Index isAbsorptionLinesOK(void *);
    DLL_PUBLIC void * getSpeciesNameAbsorptionLines(void *);

    // EnergyLevelMap
    BasicInterfaceCAPI(EnergyLevelMap)
    BasicInputOutputCAPI(EnergyLevelMap)
    EnumGetterSetterCAPI(EnergyLevelMap, Type, EnergyLevelMapType)
    VoidGetterCAPI(EnergyLevelMap, Levels)
    VoidGetterCAPI(EnergyLevelMap, Energies)
    VoidGetterCAPI(EnergyLevelMap, Data)
    DLL_PUBLIC bool getOKEnergyLevelMap(void *);
    
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
    DLL_PUBLIC bool getOKPropagationMatrix(void *);
    
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
    DLL_PUBLIC bool getOKStokesVector(void *);
    
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
    BasicInterfaceCAPI(GridPos)
    BasicInputOutputCAPI(GridPos)
    VoidStructGetterCAPI(GridPos, idx)
    VoidStructGetterCAPI(GridPos, fd)
    VoidArrayCAPI(ArrayOfGridPos)
    BasicInterfaceCAPI(ArrayOfGridPos)
    BasicInputOutputCAPI(ArrayOfGridPos)
    
    // LagrangeInterpolation
    BasicInterfaceCAPI(LagrangeInterpolation)
    VoidStructGetterCAPI(LagrangeInterpolation, pos)
    VoidStructGetterCAPI(LagrangeInterpolation, lx)
    VoidStructGetterCAPI(LagrangeInterpolation, dlx)
    VoidArrayCAPI(ArrayOfLagrangeInterpolation)
    BasicInterfaceCAPI(ArrayOfLagrangeInterpolation)
    BasicInputOutputCAPI(ArrayOfLagrangeInterpolation)
    
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
    
    // TransmissionMatrix
    BasicInterfaceCAPI(TransmissionMatrix)
    BasicInputOutputCAPI(TransmissionMatrix)
    VoidArrayCAPI(ArrayOfTransmissionMatrix)
    BasicInterfaceCAPI(ArrayOfTransmissionMatrix)
    BasicInputOutputCAPI(ArrayOfTransmissionMatrix)
    VoidArrayCAPI(ArrayOfArrayOfTransmissionMatrix)
    BasicInterfaceCAPI(ArrayOfArrayOfTransmissionMatrix)
    BasicInputOutputCAPI(ArrayOfArrayOfTransmissionMatrix)
    DLL_PUBLIC Numeric * getMat1TransmissionMatrix(Index, void *);
    DLL_PUBLIC Numeric * getMat2TransmissionMatrix(Index, void *);
    DLL_PUBLIC Numeric * getMat3TransmissionMatrix(Index, void *);
    DLL_PUBLIC Numeric * getMat4TransmissionMatrix(Index, void *);
    DLL_PUBLIC void setTransmissionMatrix(void *, Index, Index);
    DLL_PUBLIC Index getStokesDimTransmissionMatrix(void *);
    DLL_PUBLIC Index getFrequenciesTransmissionMatrix(void *);
    
    // RadiationVector
    BasicInterfaceCAPI(RadiationVector)
    BasicInputOutputCAPI(RadiationVector)
    VoidArrayCAPI(ArrayOfRadiationVector)
    BasicInterfaceCAPI(ArrayOfRadiationVector)
    BasicInputOutputCAPI(ArrayOfRadiationVector)
    VoidArrayCAPI(ArrayOfArrayOfRadiationVector)
    BasicInterfaceCAPI(ArrayOfArrayOfRadiationVector)
    BasicInputOutputCAPI(ArrayOfArrayOfRadiationVector)
    DLL_PUBLIC Numeric * getVec1RadiationVector(Index, void *);
    DLL_PUBLIC Numeric * getVec2RadiationVector(Index, void *);
    DLL_PUBLIC Numeric * getVec3RadiationVector(Index, void *);
    DLL_PUBLIC Numeric * getVec4RadiationVector(Index, void *);
    DLL_PUBLIC void setRadiationVector(void *, Index, Index);
    DLL_PUBLIC Index getStokesDimRadiationVector(void *);
    DLL_PUBLIC Index getFrequenciesRadiationVector(void *);
    
    // GriddedField1
    BasicInterfaceCAPI(GriddedField1)
    BasicInputOutputCAPI(GriddedField1)
    VoidArrayCAPI(ArrayOfGriddedField1)
    BasicInterfaceCAPI(ArrayOfGriddedField1)
    BasicInputOutputCAPI(ArrayOfGriddedField1)
    VoidArrayCAPI(ArrayOfArrayOfGriddedField1)
    BasicInterfaceCAPI(ArrayOfArrayOfGriddedField1)
    BasicInputOutputCAPI(ArrayOfArrayOfGriddedField1)
    DLL_PUBLIC Index get_dimGriddedField1(void * data);
    DLL_PUBLIC Index get_grid_typeIndexGriddedField1(Index i, void * data);
    DLL_PUBLIC Index get_grid_sizeGriddedField1(Index i, void * data);
    DLL_PUBLIC char * get_nameGriddedField1(void * data);
    DLL_PUBLIC void set_nameGriddedField1(void * data, char * newdata);
    DLL_PUBLIC char * get_grid_nameGriddedField1(Index i, void * data);
    DLL_PUBLIC void set_grid_nameGriddedField1(Index i, void * data, char * newdata);
    DLL_PUBLIC void * get_numeric_gridGriddedField1(Index i, void * data);
    DLL_PUBLIC void * get_string_gridGriddedField1(Index i, void * data);
    DLL_PUBLIC void set_gridGriddedField1(Index i, void * data, void * newdata, bool NumericType);
    DLL_PUBLIC void * dataGriddedField1(void * data);
    DLL_PUBLIC bool checksizeGriddedField1(void * data);
    
    // GriddedField2
    BasicInterfaceCAPI(GriddedField2)
    BasicInputOutputCAPI(GriddedField2)
    VoidArrayCAPI(ArrayOfGriddedField2)
    BasicInterfaceCAPI(ArrayOfGriddedField2)
    BasicInputOutputCAPI(ArrayOfGriddedField2)
    VoidArrayCAPI(ArrayOfArrayOfGriddedField2)
    BasicInterfaceCAPI(ArrayOfArrayOfGriddedField2)
    BasicInputOutputCAPI(ArrayOfArrayOfGriddedField2)
    DLL_PUBLIC Index get_dimGriddedField2(void * data);
    DLL_PUBLIC Index get_grid_typeIndexGriddedField2(Index i, void * data);
    DLL_PUBLIC Index get_grid_sizeGriddedField2(Index i, void * data);
    DLL_PUBLIC char * get_nameGriddedField2(void * data);
    DLL_PUBLIC void set_nameGriddedField2(void * data, char * newdata);
    DLL_PUBLIC char * get_grid_nameGriddedField2(Index i, void * data);
    DLL_PUBLIC void set_grid_nameGriddedField2(Index i, void * data, char * newdata);
    DLL_PUBLIC void * get_numeric_gridGriddedField2(Index i, void * data);
    DLL_PUBLIC void * get_string_gridGriddedField2(Index i, void * data);
    DLL_PUBLIC void set_gridGriddedField2(Index i, void * data, void * newdata, bool NumericType);
    DLL_PUBLIC void * dataGriddedField2(void * data);
    DLL_PUBLIC bool checksizeGriddedField2(void * data);
    
    // GriddedField3
    BasicInterfaceCAPI(GriddedField3)
    BasicInputOutputCAPI(GriddedField3)
    VoidArrayCAPI(ArrayOfGriddedField3)
    BasicInterfaceCAPI(ArrayOfGriddedField3)
    BasicInputOutputCAPI(ArrayOfGriddedField3)
    VoidArrayCAPI(ArrayOfArrayOfGriddedField3)
    BasicInterfaceCAPI(ArrayOfArrayOfGriddedField3)
    BasicInputOutputCAPI(ArrayOfArrayOfGriddedField3)
    DLL_PUBLIC Index get_dimGriddedField3(void * data);
    DLL_PUBLIC Index get_grid_typeIndexGriddedField3(Index i, void * data);
    DLL_PUBLIC Index get_grid_sizeGriddedField3(Index i, void * data);
    DLL_PUBLIC char * get_nameGriddedField3(void * data);
    DLL_PUBLIC void set_nameGriddedField3(void * data, char * newdata);
    DLL_PUBLIC char * get_grid_nameGriddedField3(Index i, void * data);
    DLL_PUBLIC void set_grid_nameGriddedField3(Index i, void * data, char * newdata);
    DLL_PUBLIC void * get_numeric_gridGriddedField3(Index i, void * data);
    DLL_PUBLIC void * get_string_gridGriddedField3(Index i, void * data);
    DLL_PUBLIC void set_gridGriddedField3(Index i, void * data, void * newdata, bool NumericType);
    DLL_PUBLIC void * dataGriddedField3(void * data);
    DLL_PUBLIC bool checksizeGriddedField3(void * data);
    
    // GriddedField4
    BasicInterfaceCAPI(GriddedField4)
    BasicInputOutputCAPI(GriddedField4)
    VoidArrayCAPI(ArrayOfGriddedField4)
    BasicInterfaceCAPI(ArrayOfGriddedField4)
    BasicInputOutputCAPI(ArrayOfGriddedField4)
    DLL_PUBLIC Index get_dimGriddedField4(void * data);
    DLL_PUBLIC Index get_grid_typeIndexGriddedField4(Index i, void * data);
    DLL_PUBLIC Index get_grid_sizeGriddedField4(Index i, void * data);
    DLL_PUBLIC char * get_nameGriddedField4(void * data);
    DLL_PUBLIC void set_nameGriddedField4(void * data, char * newdata);
    DLL_PUBLIC char * get_grid_nameGriddedField4(Index i, void * data);
    DLL_PUBLIC void set_grid_nameGriddedField4(Index i, void * data, char * newdata);
    DLL_PUBLIC void * get_numeric_gridGriddedField4(Index i, void * data);
    DLL_PUBLIC void * get_string_gridGriddedField4(Index i, void * data);
    DLL_PUBLIC void set_gridGriddedField4(Index i, void * data, void * newdata, bool NumericType);
    DLL_PUBLIC void * dataGriddedField4(void * data);
    DLL_PUBLIC bool checksizeGriddedField4(void * data);
    
    // GriddedField5
    BasicInterfaceCAPI(GriddedField5)
    BasicInputOutputCAPI(GriddedField5)
    DLL_PUBLIC Index get_dimGriddedField5(void * data);
    DLL_PUBLIC Index get_grid_typeIndexGriddedField5(Index i, void * data);
    DLL_PUBLIC Index get_grid_sizeGriddedField5(Index i, void * data);
    DLL_PUBLIC char * get_nameGriddedField5(void * data);
    DLL_PUBLIC void set_nameGriddedField5(void * data, char * newdata);
    DLL_PUBLIC char * get_grid_nameGriddedField5(Index i, void * data);
    DLL_PUBLIC void set_grid_nameGriddedField5(Index i, void * data, char * newdata);
    DLL_PUBLIC void * get_numeric_gridGriddedField5(Index i, void * data);
    DLL_PUBLIC void * get_string_gridGriddedField5(Index i, void * data);
    DLL_PUBLIC void set_gridGriddedField5(Index i, void * data, void * newdata, bool NumericType);
    DLL_PUBLIC void * dataGriddedField5(void * data);
    DLL_PUBLIC bool checksizeGriddedField5(void * data);
    
    // GriddedField6
    BasicInterfaceCAPI(GriddedField6)
    BasicInputOutputCAPI(GriddedField6)
    DLL_PUBLIC Index get_dimGriddedField6(void * data);
    DLL_PUBLIC Index get_grid_typeIndexGriddedField6(Index i, void * data);
    DLL_PUBLIC Index get_grid_sizeGriddedField6(Index i, void * data);
    DLL_PUBLIC char * get_nameGriddedField6(void * data);
    DLL_PUBLIC void set_nameGriddedField6(void * data, char * newdata);
    DLL_PUBLIC char * get_grid_nameGriddedField6(Index i, void * data);
    DLL_PUBLIC void set_grid_nameGriddedField6(Index i, void * data, char * newdata);
    DLL_PUBLIC void * get_numeric_gridGriddedField6(Index i, void * data);
    DLL_PUBLIC void * get_string_gridGriddedField6(Index i, void * data);
    DLL_PUBLIC void set_gridGriddedField6(Index i, void * data, void * newdata, bool NumericType);
    DLL_PUBLIC void * dataGriddedField6(void * data);
    DLL_PUBLIC bool checksizeGriddedField6(void * data);
    
    // CIARecord
    BasicInterfaceCAPI(CIARecord)
    BasicInputOutputCAPI(CIARecord)
    VoidGetterCAPI(CIARecord, Data)
    VoidArrayCAPI(ArrayOfCIARecord)
    BasicInterfaceCAPI(ArrayOfCIARecord)
    BasicInputOutputCAPI(ArrayOfCIARecord)
    DLL_PUBLIC void * getSpecies1CIARecord(void *);
    DLL_PUBLIC void * getSpecies2CIARecord(void *);
    DLL_PUBLIC void setSpeciesCIARecord(void *, void *, void *);
    
    // Verbosity
    BasicInterfaceCAPI(Verbosity)
    BasicInputOutputCAPI(Verbosity)
    DLL_PUBLIC Index getAgendaVerbosity(void *);
    DLL_PUBLIC Index getScreenVerbosity(void *);
    DLL_PUBLIC Index getFileVerbosity(void *);
    DLL_PUBLIC bool getMainVerbosity(void *);
    DLL_PUBLIC void setVerbosity(void *, Index, Index, Index, bool);
    
    // TessemNN
    BasicInterfaceCAPI(TessemNN)
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
    BasicInterfaceCAPI(Timer)
    BasicInputOutputCAPI(Timer)
    DLL_PUBLIC bool getrunningTimer(void *);
    DLL_PUBLIC bool getfinishedTimer(void *);
    DLL_PUBLIC Index getcputime_start_utimeTimer(void *);
    DLL_PUBLIC Index getcputime_start_stimeTimer(void *);
    DLL_PUBLIC Index getcputime_start_cutimeTimer(void *);
    DLL_PUBLIC Index getcputime_start_cstimeTimer(void *);
    DLL_PUBLIC Index getrealtime_startTimer(void *);
    DLL_PUBLIC Index getcputime_end_utimeTimer(void *);
    DLL_PUBLIC Index getcputime_end_stimeTimer(void *);
    DLL_PUBLIC Index getcputime_end_cutimeTimer(void *);
    DLL_PUBLIC Index getcputime_end_cstimeTimer(void *);
    DLL_PUBLIC Index getrealtime_endTimer(void *);
    DLL_PUBLIC void setrunningTimer(void *, bool);
    DLL_PUBLIC void setfinishedTimer(void *, bool);
    DLL_PUBLIC void setcputime_start_utimeTimer(void *, Index);
    DLL_PUBLIC void setcputime_start_stimeTimer(void *, Index);
    DLL_PUBLIC void setcputime_start_cutimeTimer(void *, Index);
    DLL_PUBLIC void setcputime_start_cstimeTimer(void *, Index);
    DLL_PUBLIC void setrealtime_startTimer(void *, Index);
    DLL_PUBLIC void setcputime_end_utimeTimer(void *, Index);
    DLL_PUBLIC void setcputime_end_stimeTimer(void *, Index);
    DLL_PUBLIC void setcputime_end_cutimeTimer(void *, Index);
    DLL_PUBLIC void setcputime_end_cstimeTimer(void *, Index);
    DLL_PUBLIC void setrealtime_endTimer(void *, Index);
    DLL_PUBLIC bool supportTimer();
    DLL_PUBLIC Index tickTimer();
    
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
    DLL_PUBLIC Numeric getA0_K0TelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getA0_K1TelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getA0_K2TelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getA0_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getA1_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getA2_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getA3_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getB0_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getB1_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getB2_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getB3_EVEHTelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getRAPPORT43_32TelsemAtlas(Index, void*);
    DLL_PUBLIC Numeric getRAPPORT54_43TelsemAtlas(Index, void*);
    
    // MCAntenna
    BasicInterfaceCAPI(MCAntenna)
    BasicInputOutputCAPI(MCAntenna)
    VoidGetterCAPI(MCAntenna, saa)
    VoidGetterCAPI(MCAntenna, sza)
    VoidGetterCAPI(MCAntenna, aag)
    VoidGetterCAPI(MCAntenna, zag)
    VoidGetterCAPI(MCAntenna, G)
    DLL_PUBLIC Index getTypeMCAntenna(void *);
    DLL_PUBLIC Index setTypeMCAntenna(void *, Index);
    
    // GasAbsLookup
    BasicInterfaceCAPI(GasAbsLookup)
    BasicInputOutputCAPI(GasAbsLookup)
    VoidGetterCAPI(GasAbsLookup, Species)
    VoidGetterCAPI(GasAbsLookup, NonLinearSpecies)
    VoidGetterCAPI(GasAbsLookup, Fgrid)
    VoidGetterCAPI(GasAbsLookup, FLAGDefault)
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
    DLL_PUBLIC Index getVersionXsecRecord(void *);
    DLL_PUBLIC void setVersionXsecRecord(void *, Index);
    DLL_PUBLIC void * getSpeciesXsecRecord(void *);
    DLL_PUBLIC void setSpeciesXsecRecord(void *, void *);
    VoidGetterCAPI(XsecRecord, FitMinPressures)
    VoidGetterCAPI(XsecRecord, FitMaxPressures)
    VoidGetterCAPI(XsecRecord, FitMinTemperatures)
    VoidGetterCAPI(XsecRecord, FitMaxTemperatures)
    VoidGetterCAPI(XsecRecord, FitCoeffs)
    VoidArrayCAPI(ArrayOfXsecRecord)
    BasicInterfaceCAPI(ArrayOfXsecRecord)
    BasicInputOutputCAPI(ArrayOfXsecRecord)
    
    // Sparse
    BasicInterfaceCAPI(Sparse)
    BasicInputOutputCAPI(Sparse)
    VoidArrayCAPI(ArrayOfSparse)
    BasicInterfaceCAPI(ArrayOfSparse)
    BasicInputOutputCAPI(ArrayOfSparse)
    DLL_PUBLIC void resizeSparse(Index, Index, void *);
    DLL_PUBLIC Index rowsSparse(void *);
    DLL_PUBLIC Index colsSparse(void *);
    DLL_PUBLIC Index sizeSparse(void *);
    DLL_PUBLIC int * rowsptrSparse(void *);
    DLL_PUBLIC int * colsptrSparse(void *);
    DLL_PUBLIC Numeric * getDataSparse(void *);
    DLL_PUBLIC void setDataSparse(void *, Index, Index, Numeric);
    
    // CovarianceMatrix
    BasicInterfaceCAPI(CovarianceMatrix)
    BasicInputOutputCAPI(CovarianceMatrix)
    VoidArrayElemCAPI(CovarianceMatrix, get_blocks)
    VoidArrayElemCAPI(CovarianceMatrix, get_inverse_blocks)
    
    // Any
    BasicInterfaceCAPI(Any)
    BasicInputOutputCAPI(Any)
    
    // Agenda
    BasicInterfaceCAPI(Agenda)
    BasicInputOutputCAPI(Agenda)
    VoidArrayCAPI(ArrayOfAgenda)
    BasicInterfaceCAPI(ArrayOfAgenda)
    BasicInputOutputCAPI(ArrayOfAgenda)
    
    // Jacobian::Target
    BasicInterfaceCAPI(JacobianTarget)
    BasicInputOutputCAPI(JacobianTarget)
    VoidGetterCAPI(JacobianTarget, Perturbation)
    VoidGetterCAPI(JacobianTarget, QuantumIdentity)
    VoidGetterCAPI(JacobianTarget, StringKey)
    VoidGetterCAPI(JacobianTarget, SpeciesList)
    StringEnumGetterSetterCAPI(JacobianTarget, TargetType)
    StringEnumGetterSetterCAPI(JacobianTarget, TargetSubType)
    VoidArrayCAPI(ArrayOfJacobianTarget)
    BasicInterfaceCAPI(ArrayOfJacobianTarget)
    BasicInputOutputCAPI(ArrayOfJacobianTarget)
    
    // RetrievalQuantity
    BasicInterfaceCAPI(RetrievalQuantity)
    BasicInputOutputCAPI(RetrievalQuantity)
    VoidGetterCAPI(RetrievalQuantity, SubTag)
    VoidGetterCAPI(RetrievalQuantity, SubSubTag)
    VoidGetterCAPI(RetrievalQuantity, Mode)
    VoidGetterCAPI(RetrievalQuantity, Grids)
    VoidGetterCAPI(RetrievalQuantity, Target)
    VoidGetterCAPI(RetrievalQuantity, TransformationFunc)
    VoidGetterCAPI(RetrievalQuantity, TFuncParameters)
    VoidGetterCAPI(RetrievalQuantity, Transformation)
    VoidGetterCAPI(RetrievalQuantity, Offset)
    GetterSetterCAPI(RetrievalQuantity, Integration, bool)
    VoidArrayCAPI(ArrayOfRetrievalQuantity)
    BasicInterfaceCAPI(ArrayOfRetrievalQuantity)
    BasicInputOutputCAPI(ArrayOfRetrievalQuantity)
    
    // Range
    BasicInterfaceCAPI(Range)
    DLL_PUBLIC Index get_startRange(void *);
    DLL_PUBLIC Index get_strideRange(void *);
    DLL_PUBLIC Index get_extentRange(void *);
    DLL_PUBLIC void setRange(void *, Index, Index, Index);
    
    // Block
    BasicInterfaceCAPI(Block)
    VoidGetterCAPI(Block, get_row_range)
    VoidGetterCAPI(Block, get_column_range)
    VoidGetterCAPI(Block, get_dense)
    VoidGetterCAPI(Block, get_sparse)
    DLL_PUBLIC Index get_matrix_typeBlock(void *);
    DLL_PUBLIC Index get_index1Block(void *);
    DLL_PUBLIC Index get_index2Block(void *);
    DLL_PUBLIC void set_indicesBlock(void *, Index, Index);
    DLL_PUBLIC void set_matrixBlock(void *, void *, bool);
    
    // Time
    BasicInterfaceCAPI(Time)
    BasicInputOutputCAPI(Time)
    GetterSetterCAPI(Time, Seconds, Numeric)
    VoidArrayCAPI(ArrayOfTime)
    BasicInterfaceCAPI(ArrayOfTime)
    BasicInputOutputCAPI(ArrayOfTime)
    VoidArrayCAPI(ArrayOfArrayOfTime)
    BasicInterfaceCAPI(ArrayOfArrayOfTime)
    BasicInputOutputCAPI(ArrayOfArrayOfTime)
    DLL_PUBLIC Index setTimeFromString(void *, char *);
    DLL_PUBLIC void setTime(void *, void *);
    DLL_PUBLIC bool equalTime(void *, void *);
    DLL_PUBLIC bool lessTime(void *, void *);
    
    // HitranRelaxationMatrixData
    BasicInterfaceCAPI(HitranRelaxationMatrixData)
    BasicInputOutputCAPI(HitranRelaxationMatrixData)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0rr)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0rr)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0rq)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0rq)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0rp)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0rp)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0qr)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0qr)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0qq)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0qq)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0qp)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0qp)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0pr)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0pr)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0pq)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0pq)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, W0pp)
    VoidStructGetterCAPI(HitranRelaxationMatrixData, B0pp)
    
    // PartitionFunctionsType
    BasicInterfaceCAPI(PartitionFunctionsType)
    StringEnumPointersCAPI(PartitionFunctionsType)
    
    // PartitionFunctionsData
    BasicInterfaceCAPI(PartitionFunctionsData)
    BasicInputOutputCAPI(PartitionFunctionsData)
    VoidStructGetterCAPI(PartitionFunctionsData, type)
    VoidStructGetterCAPI(PartitionFunctionsData, data)
    
    // SpeciesErrorCorrectedSuddenData
    BasicInterfaceCAPI(SpeciesErrorCorrectedSuddenData)
    VoidStructGetterCAPI(SpeciesErrorCorrectedSuddenData, spec)
    VoidStructGetterCAPI(SpeciesErrorCorrectedSuddenData, scaling)
    VoidStructGetterCAPI(SpeciesErrorCorrectedSuddenData, beta)
    VoidStructGetterCAPI(SpeciesErrorCorrectedSuddenData, lambda)
    VoidStructGetterCAPI(SpeciesErrorCorrectedSuddenData, collisional_distance)
    VoidStructGetterCAPI(SpeciesErrorCorrectedSuddenData, mass)
    
    // ErrorCorrectedSuddenData
    BasicInterfaceCAPI(ErrorCorrectedSuddenData)
    VoidStructGetterCAPI(ErrorCorrectedSuddenData, id)
    DLL_PUBLIC void * getErrorCorrectedSuddenData(void *, void *);
    DLL_PUBLIC Index getnelemErrorCorrectedSuddenData(void *);
    DLL_PUBLIC void * getSpeciesErrorCorrectedSuddenDataAtErrorCorrectedSuddenData(void *, Index);
    
    // MapOfErrorCorrectedSuddenData
    BasicInterfaceCAPI(MapOfErrorCorrectedSuddenData)
    BasicInputOutputCAPI(MapOfErrorCorrectedSuddenData)
    DLL_PUBLIC void * getMapOfErrorCorrectedSuddenData(void *, void *);
    DLL_PUBLIC Index getnelemMapOfErrorCorrectedSuddenData(void *);
    DLL_PUBLIC void * getErrorCorrectedSuddenDataAtMapOfErrorCorrectedSuddenData(void *, Index);
    
    // generic
    DLL_PUBLIC Index string2filetypeindex(char *);
    DLL_PUBLIC void * get_list_of_all_workspace_classes();
}


#undef BasicInterfaceCAPI
#undef GetterSetterCAPI
#undef EnumGetterSetterCAPI
#undef StringEnumGetterSetterCAPI
#undef VoidGetterCAPI
#undef VoidStructGetterCAPI
#undef BasicInputOutputCAPI
#undef VoidArrayCAPI
#undef VoidArrayElemCAPI
#undef StringEnumPointersCAPI
#undef EnumMacroInterfaceCAPI

#if REMOVE_DLL_PUBLIC
#undef DLL_PUBLIC
#endif
#undef REMOVE_DLL_PUBLIC

#endif // _ARTS_ARTS_API_CLASS_H_


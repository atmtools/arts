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

extern "C" {
  // Zeeman::Model
  DLL_PUBLIC void * createZeemanModel();
  DLL_PUBLIC void deleteZeemanModel(void *);
  DLL_PUBLIC Numeric getZeemanModelGU(void *);
  DLL_PUBLIC Numeric getZeemanModelGL(void *);
  DLL_PUBLIC void setZeemanModelGU(void *, Numeric);
  DLL_PUBLIC void setZeemanModelGL(void *, Numeric);
  DLL_PUBLIC void printZeemanModel(void *);
  
  // Rational
  DLL_PUBLIC void * createRational();
  DLL_PUBLIC void deleteRational(void *);
  DLL_PUBLIC Index getRationalNom(void *);
  DLL_PUBLIC Index getRationalDenom(void *);
  DLL_PUBLIC void setRationalNom(void *, Index);
  DLL_PUBLIC void setRationalDenom(void *, Index);
  DLL_PUBLIC void printRational(void *);
  DLL_PUBLIC void simplifyRational(void *);
  
  // LineShape::ModelParameters
  DLL_PUBLIC void printLineShapeModelParameters(void *);
  DLL_PUBLIC Index getLineShapeModelParametersType(char *);
  
  // LineShape::SingleSpeciesModel
  DLL_PUBLIC void * createLineShapeSingleSpeciesModel();
  DLL_PUBLIC void deleteLineShapeSingleSpeciesModel(void *);
  DLL_PUBLIC void printLineShapeSingleSpeciesModel(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelG0(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelD0(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelG2(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelD2(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelFVC(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelETA(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelY(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelG(void *);
  DLL_PUBLIC void * getLineShapeSingleSpeciesModelDV(void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelG0(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelD0(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelG2(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelD2(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelFVC(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelETA(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelY(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelG(void *, void *);
  DLL_PUBLIC void setLineShapeSingleSpeciesModelDV(void *, void *);
  
  // LineShape::Model
  DLL_PUBLIC void * createLineShapeModel();
  DLL_PUBLIC void deleteLineShapeModel(void *);
  DLL_PUBLIC void printLineShapeModel(void *);
  DLL_PUBLIC void resizeLineShapeModel(Index, void *);
  DLL_PUBLIC Index sizeLineShapeModel(void *);
  DLL_PUBLIC void * getLineShapeModelSingleSpeciesModel(Index, void *);
  
  // Absorption::SingleLine
  DLL_PUBLIC void * createAbsorptionSingleLine();
  DLL_PUBLIC void deleteAbsorptionSingleLine(void *);
  DLL_PUBLIC void printAbsorptionSingleLine(void *);
  DLL_PUBLIC Numeric getAbsorptionSingleLineF0(void *);
  DLL_PUBLIC Numeric getAbsorptionSingleLineI0(void *);
  DLL_PUBLIC Numeric getAbsorptionSingleLineE0(void *);
  DLL_PUBLIC Numeric getAbsorptionSingleLineGL(void *);
  DLL_PUBLIC Numeric getAbsorptionSingleLineGU(void *);
  DLL_PUBLIC Numeric getAbsorptionSingleLineA(void *);
  DLL_PUBLIC void setAbsorptionSingleLineF0(void *, Numeric);
  DLL_PUBLIC void setAbsorptionSingleLineI0(void *, Numeric);
  DLL_PUBLIC void setAbsorptionSingleLineE0(void *, Numeric);
  DLL_PUBLIC void setAbsorptionSingleLineGL(void *, Numeric);
  DLL_PUBLIC void setAbsorptionSingleLineGU(void *, Numeric);
  DLL_PUBLIC void setAbsorptionSingleLineA(void *, Numeric);
  DLL_PUBLIC void * getAbsorptionSingleLineZeemanModel(void *);
  DLL_PUBLIC void * getAbsorptionSingleLineLineShapeModel(void *);
  DLL_PUBLIC Index sizeAbsorptionSingleLineLowerQuantas(void *);
  DLL_PUBLIC Index sizeAbsorptionSingleLineUpperQuantas(void *);
  DLL_PUBLIC void resizeAbsorptionSingleLineLowerQuantas(Index, void *);
  DLL_PUBLIC void resizeAbsorptionSingleLineUpperQuantas(Index, void *);
  DLL_PUBLIC void * getAbsorptionSingleLineLowerQuanta(Index, void *);
  DLL_PUBLIC void * getAbsorptionSingleLineUpperQuanta(Index, void *);
  
  // QuantumNumbers
  DLL_PUBLIC void * createQuantumNumbers();
  DLL_PUBLIC void deleteQuantumNumbers(void *);
  DLL_PUBLIC void printQuantumNumbers(void *);
  DLL_PUBLIC void * getQuantumNumbersNumber(Index, void *);
  
  // QuantumIdentifier
  DLL_PUBLIC void * createQuantumIdentifier();
  DLL_PUBLIC void deleteQuantumIdentifier(void *);
  DLL_PUBLIC void printQuantumIdentifier(void *);
  DLL_PUBLIC Index getQuantumIdentifierType(void *);
  DLL_PUBLIC Index setQuantumIdentifierTypeFromString(void *, char *);
  DLL_PUBLIC Index setQuantumIdentifierTypeFromIndex(void *, Index);
  DLL_PUBLIC Index getQuantumIdentifierSpecies(void *);
  DLL_PUBLIC void setQuantumIdentifierSpecies(void *, Index);
  DLL_PUBLIC Index getQuantumIdentifierIsotopologue(void *);
  DLL_PUBLIC void setQuantumIdentifierIsotopologue(void *, Index);
  DLL_PUBLIC void * getQuantumIdentifierLevelQuantumNumbers(void *);
  DLL_PUBLIC void * getQuantumIdentifierLowerQuantumNumbers(void *);
  DLL_PUBLIC void * getQuantumIdentifierUpperQuantumNumbers(void *);
  
  // SpeciesTag
  DLL_PUBLIC void * createSpeciesTag();
  DLL_PUBLIC void deleteSpeciesTag(void *);
  DLL_PUBLIC void printSpeciesTag(void *);
  DLL_PUBLIC Index setSpeciesTag(void *, char *);
  DLL_PUBLIC Index getSpeciesTagSpecies(void *);
  DLL_PUBLIC Index getSpeciesTagIsotopologue(void *);
  DLL_PUBLIC Numeric getSpeciesTagLowerFrequency(void *);
  DLL_PUBLIC Numeric getSpeciesTagUpperFrequency(void *);
  DLL_PUBLIC Index getSpeciesTagType(void *);
  DLL_PUBLIC Index getSpeciesTagCIASecond(void *);
  DLL_PUBLIC Index getSpeciesTagCIADataset(void *);
  DLL_PUBLIC void setSpeciesTagSpecies(void *, Index);
  DLL_PUBLIC void setSpeciesTagIsotopologue(void *, Index);
  DLL_PUBLIC void setSpeciesTagLowerFrequency(void *, Numeric);
  DLL_PUBLIC void setSpeciesTagUpperFrequency(void *, Numeric);
  DLL_PUBLIC void setSpeciesTagType(void *, Index);
  DLL_PUBLIC void setSpeciesTagCIASecond(void *, Index);
  DLL_PUBLIC void setSpeciesTagCIADataset(void *, Index);
  
  // Absorption::Lines
  DLL_PUBLIC void * createAbsorptionLines();
  DLL_PUBLIC void deleteAbsorptionLines(void *);
  DLL_PUBLIC void printAbsorptionLines(void *);
  DLL_PUBLIC Index xmlreadAbsorptionLines(void *, char *);
  DLL_PUBLIC Index xmlsaveAbsorptionLines(void *, char *, Index, Index);
  DLL_PUBLIC void printmetaAbsorptionLines(void *);
  DLL_PUBLIC bool getAbsorptionLinesSelfBroadening(void *);
  DLL_PUBLIC void setAbsorptionLinesSelfBroadening(void *, bool);
  DLL_PUBLIC bool getAbsorptionLinesBathBroadening(void *);
  DLL_PUBLIC void setAbsorptionLinesBathBroadening(void *, bool);
  DLL_PUBLIC Index getAbsorptionLinesCutoffType(void *);
  DLL_PUBLIC Index setAbsorptionLinesCutoffType(void *, char *);
  DLL_PUBLIC void setAbsorptionLinesCutoffTypeByIndex(void *, Index);
  DLL_PUBLIC Index getAbsorptionLinesMirroringType(void *);
  DLL_PUBLIC Index setAbsorptionLinesMirroringType(void *, char *);
  DLL_PUBLIC void setAbsorptionLinesMirroringTypeByIndex(void *, Index);
  DLL_PUBLIC Index getAbsorptionLinesPopulationType(void *);
  DLL_PUBLIC Index setAbsorptionLinesPopulationType(void *, char *);
  DLL_PUBLIC void setAbsorptionLinesPopulationTypeByIndex(void *, Index);
  DLL_PUBLIC Index getAbsorptionLinesNormalizationType(void *);
  DLL_PUBLIC Index setAbsorptionLinesNormalizationType(void *, char *);
  DLL_PUBLIC void setAbsorptionLinesNormalizationTypeByIndex(void *, Index);
  DLL_PUBLIC Index getAbsorptionLinesLineShapeType(void *);
  DLL_PUBLIC Index setAbsorptionLinesLineShapeType(void *, char *);
  DLL_PUBLIC void setAbsorptionLinesLineShapeTypeByIndex(void *, Index);
  DLL_PUBLIC Numeric getAbsorptionLinesT0(void *);
  DLL_PUBLIC void setAbsorptionLinesT0(void *, Numeric);
  DLL_PUBLIC Numeric getAbsorptionLinesCutoffFrequency(void *);
  DLL_PUBLIC void setAbsorptionLinesCutoffFrequency(void *, Numeric);
  DLL_PUBLIC Numeric getAbsorptionLinesLinemixingLimit(void *);
  DLL_PUBLIC void setAbsorptionLinesLinemixingLimit(void *, Numeric);
  DLL_PUBLIC void * getAbsorptionLinesQuantumIdentifier(void *);
  DLL_PUBLIC void resizeAbsorptionLinesLocalQuantumNumber(Index, void *);
  DLL_PUBLIC Index getAbsorptionLinesLocalQuantumNumber(Index, void *);
  DLL_PUBLIC void setAbsorptionLinesLocalQuantumNumber(Index, void *, Index);
  DLL_PUBLIC Index getAbsorptionLinesLocalQuantumNumberCount(void *);
  DLL_PUBLIC void resizeAbsorptionLinesSpeciesTag(Index, void *);
  DLL_PUBLIC void * getAbsorptionLinesSpeciesTag(Index, void *);
  DLL_PUBLIC Index getAbsorptionLinesSpeciesTagCount(void *);
  DLL_PUBLIC Index isAbsorptionLinesOK(void *);
  DLL_PUBLIC void resizeAbsorptionLinesSingleLine(Index, void *);
  DLL_PUBLIC void * getAbsorptionLinesSingleLine(Index, void *);
  DLL_PUBLIC Index getAbsorptionLinesSingleLineCount(void *);
  
  // generic
  DLL_PUBLIC Index validSpecies(Index);
  DLL_PUBLIC Index validIsotopologue(Index, Index);
  DLL_PUBLIC Index getQuantumNumbersMaxNumber();
  DLL_PUBLIC Index string2quantumnumbersindex(char *);
  DLL_PUBLIC Index string2filetypeindex(char *);
}

#if REMOVE_DLL_PUBLIC
#undef DLL_PUBLIC
#endif
#undef REMOVE_DLL_PUBLIC

#endif // _ARTS_ARTS_API_CLASS_H_


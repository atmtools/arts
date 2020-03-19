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
#include "global_data.h"
#include "lineshapemodel.h"
#include "quantum.h"
#include "xml_io.h"
#include "xml_io_types.h"
#include "zeemandata.h"
  
void * createZeemanModel()
{
  return new Zeeman::Model;
}

void deleteZeemanModel(void * data)
{
  delete static_cast<Zeeman::Model *>(data);
}

void printZeemanModel(void * data)
{
  std::cout << (*static_cast<Zeeman::Model *>(data)) << std::endl;
}

Numeric getZeemanModelGU(void * data)
{
  return static_cast<Zeeman::Model *>(data)->gu();
}

Numeric getZeemanModelGL(void * data)
{
  return static_cast<Zeeman::Model *>(data)->gl();
}

void setZeemanModelGU(void * data, Numeric newdata)
{
  static_cast<Zeeman::Model *>(data)->gu() = newdata;
}

void setZeemanModelGL(void * data, Numeric newdata)
{
  static_cast<Zeeman::Model *>(data)->gl() = newdata;
}


void * createRational()
{
  return new Rational;
}

void deleteRational(void * data)
{
  delete static_cast<Rational *>(data);
}

void printRational(void * data)
{
  std::cout << (*static_cast<Rational *>(data)) << std::endl;
}

Index getRationalNom(void * data)
{
  return static_cast<Rational *>(data)->Nom();
}

Index getRationalDenom(void * data)
{
  return static_cast<Rational *>(data)->Denom();
}

void setRationalNom(void * data, Index newdata)
{
  static_cast<Rational *>(data)->Nom() = newdata;
}

void setRationalDenom(void * data, Index newdata)
{
  static_cast<Rational *>(data)->Denom() = newdata;
}

void simplifyRational(void * data)
{
  static_cast<Rational *>(data)->simplify_in_place();
}


void printLineShapeModelParameters(void * data)
{
  std::cout << (*static_cast<LineShape::ModelParameters *>(data)) << std::endl;
}

Index getLineShapeModelParametersType(char * data)
{
  try {
    return Index(LineShape::string2temperaturemodel(data));
  } catch (std::runtime_error& e) {
    return -1;
  }
}


void * createLineShapeSingleSpeciesModel()
{
  return new LineShape::SingleSpeciesModel;
}

void deleteLineShapeSingleSpeciesModel(void * data)
{
  delete static_cast<LineShape::SingleSpeciesModel *>(data);
}

void printLineShapeSingleSpeciesModel(void * data)
{
  std::cout << (*static_cast<LineShape::SingleSpeciesModel *>(data)) << std::endl;
}

void * getLineShapeSingleSpeciesModelG0(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->G0();
}

void setLineShapeSingleSpeciesModelG0(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->G0() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelD0(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->D0();
}

void setLineShapeSingleSpeciesModelD0(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->D0() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelG2(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->G2();
}

void setLineShapeSingleSpeciesModelG2(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->G2() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelD2(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->D2();
}

void setLineShapeSingleSpeciesModelD2(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->D2() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelFVC(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->FVC();
}

void setLineShapeSingleSpeciesModelFVC(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->FVC() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelETA(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->ETA();
}

void setLineShapeSingleSpeciesModelETA(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->ETA() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelY(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->Y();
}

void setLineShapeSingleSpeciesModelY(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->Y() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelG(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->G();
}

void setLineShapeSingleSpeciesModelG(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->G() = *static_cast<LineShape::ModelParameters *>(newdata);
}

void * getLineShapeSingleSpeciesModelDV(void * data)
{
  return &static_cast<LineShape::SingleSpeciesModel *>(data)->DV();
}

void setLineShapeSingleSpeciesModelDV(void * data, void * newdata)
{
  static_cast<LineShape::SingleSpeciesModel *>(data)->DV() = *static_cast<LineShape::ModelParameters *>(newdata);
}


void * createLineShapeModel()
{
  return new LineShape::Model;
}

void deleteLineShapeModel(void * data)
{
  delete static_cast<LineShape::Model *>(data);
}

void printLineShapeModel(void * data)
{
  std::cout << (*static_cast<LineShape::Model *>(data)) << std::endl;
}

Index sizeLineShapeModel(void * data)
{
  return static_cast<LineShape::Model *>(data)->nelem();
}

void resizeLineShapeModel(Index n, void * data)
{
  static_cast<LineShape::Model *>(data)->resize(n);
}

void * getLineShapeModelSingleSpeciesModel(Index i, void * data)
{
  return &static_cast<LineShape::Model *>(data)->Data()[i];
}


void * createAbsorptionSingleLine()
{
  return new Absorption::SingleLine;
}

void deleteAbsorptionSingleLine(void * data)
{
  delete static_cast<Absorption::SingleLine *>(data);
}

void printAbsorptionSingleLine(void * data)
{
  std::cout << (*static_cast<Absorption::SingleLine *>(data)) << std::endl;
}

Numeric getAbsorptionSingleLineF0(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->F0();
}

Numeric getAbsorptionSingleLineI0(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->I0();
}

Numeric getAbsorptionSingleLineE0(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->E0();
}

Numeric getAbsorptionSingleLineGL(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->g_low();
}

Numeric getAbsorptionSingleLineGU(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->g_upp();
}

Numeric getAbsorptionSingleLineA(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->A();
}

void setAbsorptionSingleLineF0(void * data, Numeric newdata)
{
  static_cast<Absorption::SingleLine *>(data)->F0() = newdata;
}

void setAbsorptionSingleLineI0(void * data, Numeric newdata)
{
  static_cast<Absorption::SingleLine *>(data)->I0() = newdata;
}

void setAbsorptionSingleLineE0(void * data, Numeric newdata)
{
  static_cast<Absorption::SingleLine *>(data)->E0() = newdata;
}

void setAbsorptionSingleLineGL(void * data, Numeric newdata)
{
  static_cast<Absorption::SingleLine *>(data)->g_low() = newdata;
}

void setAbsorptionSingleLineGU(void * data, Numeric newdata)
{
  static_cast<Absorption::SingleLine *>(data)->g_upp() = newdata;
}

void setAbsorptionSingleLineA(void * data, Numeric newdata)
{
  static_cast<Absorption::SingleLine *>(data)->A() = newdata;
}

void * getAbsorptionSingleLineZeemanModel(void * data)
{
  return &static_cast<Absorption::SingleLine *>(data)->Zeeman();
}

void * getAbsorptionSingleLineLineShapeModel(void * data)
{
  return &static_cast<Absorption::SingleLine *>(data)->LineShape();
}

Index sizeAbsorptionSingleLineLowerQuantas(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->LowerQuantumNumbers().size();
}

Index sizeAbsorptionSingleLineUpperQuantas(void * data)
{
  return static_cast<Absorption::SingleLine *>(data)->UpperQuantumNumbers().size();
}

void resizeAbsorptionSingleLineLowerQuantas(Index n, void * data)
{
  static_cast<Absorption::SingleLine *>(data)->LowerQuantumNumbers().resize(n);
}

void resizeAbsorptionSingleLineUpperQuantas(Index n, void * data)
{
  static_cast<Absorption::SingleLine *>(data)->UpperQuantumNumbers().resize(n);
}

void * getAbsorptionSingleLineLowerQuanta(Index i, void * data)
{
  return &static_cast<Absorption::SingleLine *>(data)->LowerQuantumNumbers()[i];
}

void * getAbsorptionSingleLineUpperQuanta(Index i, void * data)
{
  return &static_cast<Absorption::SingleLine *>(data)->UpperQuantumNumbers()[i];
}


void * createQuantumNumbers()
{
  return new QuantumNumbers;
}

void deleteQuantumNumbers(void * data)
{
  delete static_cast<QuantumNumbers *>(data);
}

void printQuantumNumbers(void * data)
{
  std::cout << (*static_cast<QuantumNumbers *>(data)) << std::endl;
}

Index getQuantumNumbersMaxNumber()
{
    return Index(QuantumNumberType::FINAL_ENTRY);
}

void * getQuantumNumbersNumber(Index i, void * data)
{
    return &static_cast<QuantumNumbers *>(data)->operator[](i);
}

Index string2quantumnumbersindex(char * str)
{
    return Index(string2quantumnumbertype(str));
}


void * createQuantumIdentifier()
{
  return new QuantumIdentifier;
}

void deleteQuantumIdentifier(void * data)
{
  delete static_cast<QuantumIdentifier *>(data);
}

void printQuantumIdentifier(void * data)
{
  std::cout << (*static_cast<QuantumIdentifier *>(data)) << std::endl;
}

Index getQuantumIdentifierType(void * data)
{
    return Index(static_cast<QuantumIdentifier *>(data)->Type());
}

Index setQuantumIdentifierTypeFromIndex(void * data, Index ind)
{
    if (QuantumIdentifier::ENERGY_LEVEL == ind) {
        static_cast<QuantumIdentifier *>(data)->SetEnergyLevel(static_cast<QuantumIdentifier *>(data)->EnergyLevelQuantumNumbers());
        return EXIT_SUCCESS;
    } else if (QuantumIdentifier::TRANSITION == ind) {
        static_cast<QuantumIdentifier *>(data)->SetTransition();
        return EXIT_SUCCESS;
    } 
    else if (QuantumIdentifier::ALL == ind) {
        static_cast<QuantumIdentifier *>(data)->SetAll();
        return EXIT_SUCCESS;
    } 
    else if (QuantumIdentifier::NONE == ind) {
        static_cast<QuantumIdentifier *>(data)->SetNone();
        return EXIT_SUCCESS;
    } else
        return EXIT_FAILURE;
}

Index setQuantumIdentifierTypeFromString(void * data, char * str)
{
    if (std::string("ENERGY_LEVEL") == str) {
        static_cast<QuantumIdentifier *>(data)->SetEnergyLevel(static_cast<QuantumIdentifier *>(data)->EnergyLevelQuantumNumbers());
        return EXIT_SUCCESS;
    } else if (std::string("TRANSITION") == str) {
        static_cast<QuantumIdentifier *>(data)->SetTransition();
        return EXIT_SUCCESS;
    } 
    else if (std::string("ALL") == str) {
        static_cast<QuantumIdentifier *>(data)->SetAll();
        return EXIT_SUCCESS;
    } 
    else if (std::string("NONE") == str) {
        static_cast<QuantumIdentifier *>(data)->SetNone();
        return EXIT_SUCCESS;
    } else
        return EXIT_FAILURE;
}

Index getQuantumIdentifierSpecies(void * data)
{
    return static_cast<QuantumIdentifier *>(data)->Species();
}

Index validSpecies(Index spec)
{
    if (spec >= 0 and spec < global_data::species_data.nelem())
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}

void setQuantumIdentifierSpecies(void * data, Index spec)
{
    static_cast<QuantumIdentifier *>(data)->Species() = spec;
}

Index getQuantumIdentifierIsotopologue(void * data)
{
    return static_cast<QuantumIdentifier *>(data)->Isotopologue();
}

Index validIsotopologue(Index spec, Index isot)
{
    auto& species = global_data::species_data[spec];
    if (isot >= 0 and isot < species.Isotopologue().nelem()) {
        if (not species.Isotopologue()[isot].isContinuum()) {
            return EXIT_SUCCESS;
        }
    }
    
    return EXIT_FAILURE;
}

void setQuantumIdentifierIsotopologue(void * data, Index isot)
{
    static_cast<QuantumIdentifier *>(data)->Isotopologue() = isot;
}

void * getQuantumIdentifierLevelQuantumNumbers(void * data)
{
    return &static_cast<QuantumIdentifier *>(data) ->EnergyLevelQuantumNumbers();
}

void * getQuantumIdentifierLowerQuantumNumbers(void * data)
{
    return &static_cast<QuantumIdentifier *>(data) ->LowerQuantumNumbers();
}

void * getQuantumIdentifierUpperQuantumNumbers(void * data)
{
    return &static_cast<QuantumIdentifier *>(data) ->UpperQuantumNumbers();
}


void * createSpeciesTag()
{
  return new SpeciesTag;
}

void deleteSpeciesTag(void * data)
{
    delete static_cast<SpeciesTag *>(data);
}

void printSpeciesTag(void * data)
{
    std::cout << (*static_cast<SpeciesTag *>(data)) << std::endl;
}

Index setSpeciesTag(void * data, char * newdata)
{
    try {
      *static_cast<SpeciesTag *>(data) = SpeciesTag(newdata);
      return EXIT_SUCCESS;
    } catch(std::exception& e) {
      return EXIT_FAILURE;
    }
}

Index getSpeciesTagSpecies(void * data)
{
    return static_cast<SpeciesTag *>(data)->Species();
}

Index getSpeciesTagIsotopologue(void * data)
{
    return static_cast<SpeciesTag *>(data)->Isotopologue();
}

Numeric getSpeciesTagLowerFrequency(void * data)
{
    return static_cast<SpeciesTag *>(data)->Lf();
}

Numeric getSpeciesTagUpperFrequency(void * data)
{
    return static_cast<SpeciesTag *>(data)->Uf();
}

Index getSpeciesTagType(void * data)
{
    return static_cast<SpeciesTag *>(data)->Type();
}

Index getSpeciesTagCIASecond(void * data)
{
    return static_cast<SpeciesTag *>(data)->CIASecond();
}

Index getSpeciesTagCIADataset(void * data)
{
    return static_cast<SpeciesTag *>(data)->CIADataset();
}

void setSpeciesTagSpecies(void * data, Index newdata)
{
  static_cast<SpeciesTag *>(data)->Species(newdata);
}

void setSpeciesTagIsotopologue(void * data, Index newdata)
{
  static_cast<SpeciesTag *>(data)->Isotopologue(newdata);
}

void setSpeciesTagLowerFrequency(void * data, Numeric newdata)
{
  static_cast<SpeciesTag *>(data)->Lf(newdata);
}

void setSpeciesTagUpperFrequency(void * data, Numeric newdata)
{
  static_cast<SpeciesTag *>(data)->Uf(newdata);
}

void setSpeciesTagType(void * data, Index newdata)
{
  static_cast<SpeciesTag *>(data)->Type(newdata);
}

void setSpeciesTagCIASecond(void * data, Index newdata)
{
  static_cast<SpeciesTag *>(data)->CIASecond(newdata);
}

void setSpeciesTagCIADataset(void * data, Index newdata)
{
  static_cast<SpeciesTag *>(data)->CIADataset(newdata);
}


void * createAbsorptionLines()
{
  return new AbsorptionLines;
}

void deleteAbsorptionLines(void * data)
{
    delete static_cast<AbsorptionLines *>(data);
}

void printAbsorptionLines(void * data)
{
    std::cout << (*static_cast<AbsorptionLines *>(data)) << std::endl;
}

Index xmlreadAbsorptionLines(void * data, char * filepath)
{
    try {
        xml_read_from_file(filepath, *static_cast<AbsorptionLines *>(data), Verbosity());
        return EXIT_SUCCESS;
    } catch (std::runtime_error& e) {
        return EXIT_FAILURE;
    }
}

Index xmlsaveAbsorptionLines(void * data, char * filepath, Index filetype, Index clobber)
{
    try {
        xml_write_to_file(filepath, *static_cast<const AbsorptionLines *>(data), FileType(filetype), not clobber, Verbosity());
        return EXIT_SUCCESS;
    } catch (std::runtime_error& e) {
        return EXIT_FAILURE;
    }
}

void printmetaAbsorptionLines(void * data)
{
    std::cout << static_cast<AbsorptionLines *>(data)->MetaData() << std::endl;
}

bool getAbsorptionLinesSelfBroadening(void * data)
{
    return static_cast<AbsorptionLines *>(data)->Self();
}

void setAbsorptionLinesSelfBroadening(void * data, bool newdata)
{
    static_cast<AbsorptionLines *>(data)->Self() = newdata;
}

bool getAbsorptionLinesBathBroadening(void * data)
{
    return static_cast<AbsorptionLines *>(data)->Bath();
}

void setAbsorptionLinesBathBroadening(void * data, bool newdata)
{
    static_cast<AbsorptionLines *>(data)->Bath() = newdata;
}

Index getAbsorptionLinesCutoffType(void * data)
{
    return Index(static_cast<AbsorptionLines *>(data)->Cutoff());
}

Index setAbsorptionLinesCutoffType(void * data, char * newdata)
{
    try {
      static_cast<AbsorptionLines *>(data)->Cutoff(Absorption::string2cutofftype(newdata));
      return EXIT_SUCCESS;
    } catch(const std::exception& e) {
      return EXIT_FAILURE;
    }
}

void setAbsorptionLinesCutoffTypeByIndex(void * data, Index newdata)
{
  static_cast<AbsorptionLines *>(data)->Cutoff(Absorption::CutoffType(newdata));
}

Index getAbsorptionLinesMirroringType(void * data)
{
    return Index(static_cast<AbsorptionLines *>(data)->Mirroring());
}

Index setAbsorptionLinesMirroringType(void * data, char * newdata)
{
    try {
      static_cast<AbsorptionLines *>(data)->Mirroring(Absorption::string2mirroringtype(newdata));
      return EXIT_SUCCESS;
    } catch(const std::exception& e) {
      return EXIT_FAILURE;
    }
}

void setAbsorptionLinesMirroringTypeByIndex(void * data, Index newdata)
{
  static_cast<AbsorptionLines *>(data)->Mirroring(Absorption::MirroringType(newdata));
}

Index getAbsorptionLinesPopulationType(void * data)
{
    return Index(static_cast<AbsorptionLines *>(data)->Population());
}

Index setAbsorptionLinesPopulationType(void * data, char * newdata)
{
    try {
      static_cast<AbsorptionLines *>(data)->Population(Absorption::string2populationtype(newdata));
      return EXIT_SUCCESS;
    } catch(const std::exception& e) {
      return EXIT_FAILURE;
    }
}

void setAbsorptionLinesPopulationTypeByIndex(void * data, Index newdata)
{
  static_cast<AbsorptionLines *>(data)->Population(Absorption::PopulationType(newdata));
}

Index getAbsorptionLinesNormalizationType(void * data)
{
    return Index(static_cast<AbsorptionLines *>(data)->Normalization());
}

Index setAbsorptionLinesNormalizationType(void * data, char * newdata)
{
    try {
      static_cast<AbsorptionLines *>(data)->Normalization(Absorption::string2normalizationtype(newdata));
      return EXIT_SUCCESS;
    } catch(const std::exception& e) {
      return EXIT_FAILURE;
    }
}

void setAbsorptionLinesNormalizationTypeByIndex(void * data, Index newdata)
{
  static_cast<AbsorptionLines *>(data)->Normalization(Absorption::NormalizationType(newdata));
}

Index getAbsorptionLinesLineShapeType(void * data)
{
    return Index(static_cast<AbsorptionLines *>(data)->LineShapeType());
}

Index setAbsorptionLinesLineShapeType(void * data, char * newdata)
{
    try {
      static_cast<AbsorptionLines *>(data)->LineShapeType(LineShape::string2shapetype(newdata));
      return EXIT_SUCCESS;
    } catch(const std::exception& e) {
      return EXIT_FAILURE;
    }
}

void setAbsorptionLinesLineShapeTypeByIndex(void * data, Index newdata)
{
  static_cast<AbsorptionLines *>(data)->LineShapeType(LineShape::Type(newdata));
}

Numeric getAbsorptionLinesT0(void * data)
{
    return static_cast<AbsorptionLines *>(data)->T0();
}

void setAbsorptionLinesT0(void * data, Numeric newdata)
{
    static_cast<AbsorptionLines *>(data)->T0(newdata);
}

Numeric getAbsorptionLinesCutoffFrequency(void * data)
{
    return static_cast<AbsorptionLines *>(data)->CutoffFreqValue();
}

void setAbsorptionLinesCutoffFrequency(void * data, Numeric newdata)
{
  static_cast<AbsorptionLines *>(data)->CutoffFreqValue(newdata);
}

Numeric getAbsorptionLinesLinemixingLimit(void * data)
{
    return static_cast<AbsorptionLines *>(data)->LinemixingLimit();
}

void setAbsorptionLinesLinemixingLimit(void * data, Numeric newdata)
{
    static_cast<AbsorptionLines *>(data)->LinemixingLimit(newdata);
}

void * getAbsorptionLinesQuantumIdentifier(void * data)
{
    return &static_cast<AbsorptionLines *>(data)->QuantumIdentity();
}

void resizeAbsorptionLinesLocalQuantumNumber(Index n, void * data)
{
  static_cast<AbsorptionLines *>(data)->LocalQuanta().resize(n);
}

Index getAbsorptionLinesLocalQuantumNumber(Index i, void * data)
{
    return Index(static_cast<AbsorptionLines *>(data)->LocalQuanta()[i]);
}

void setAbsorptionLinesLocalQuantumNumber(Index i, void * data, Index newdata)
{
    static_cast<AbsorptionLines *>(data)->LocalQuanta()[i] = QuantumNumberType(newdata);
}

Index getAbsorptionLinesLocalQuantumNumberCount(void * data)
{
    return static_cast<AbsorptionLines *>(data)->NumLocalQuanta();
}

void resizeAbsorptionLinesSpeciesTag(Index n, void * data)
{
    static_cast<AbsorptionLines *>(data)->BroadeningSpecies().resize(n);
}

void * getAbsorptionLinesSpeciesTag(Index i, void * data)
{
    return &static_cast<AbsorptionLines *>(data)->BroadeningSpecies()[i];
}

Index getAbsorptionLinesSpeciesTagCount(void * data)
{
    return static_cast<AbsorptionLines *>(data)->NumBroadeners();
}

void resizeAbsorptionLinesSingleLine(Index n, void * data)
{
  static_cast<AbsorptionLines *>(data)->AllLines().resize(n);
}

void * getAbsorptionLinesSingleLine(Index i, void * data)
{
    return &static_cast<AbsorptionLines *>(data)->Line(i);
}

Index getAbsorptionLinesSingleLineCount(void * data)
{
    return static_cast<AbsorptionLines *>(data)->NumLines();
}

Index isAbsorptionLinesOK(void * data)
{
  if (static_cast<AbsorptionLines *>(data) -> OK())
    return 1;
  else
    return 0;
}


Index string2filetypeindex(char * data)
{
    try {
        return Index(string2filetype(data));
    } catch (std::runtime_error& e) {
        return -1;
    }
}


#define WorkspaceArrayInterfaceCAPI_CCFILE(BASETYPE)                                                                                               \
void * createArrayOf##BASETYPE()                                                                                                    \
{                                                                                                                                   \
    return new ArrayOf##BASETYPE;                                                                                                   \
}                                                                                                                                   \
                                                                                                                                    \
void deleteArrayOf##BASETYPE(void * data)                                                                                           \
{                                                                                                                                   \
    delete static_cast<ArrayOf##BASETYPE *>(data);                                                                                  \
}                                                                                                                                   \
                                                                                                                                    \
void printArrayOf##BASETYPE(void * data)                                                                                            \
{                                                                                                                                   \
  std::cout << (*static_cast<ArrayOf##BASETYPE *>(data)) << std::endl;                                                              \
}                                                                                                                                   \
                                                                                                                                    \
Index sizeArrayOf##BASETYPE(void * data)                                                                                            \
{                                                                                                                                   \
    return static_cast<ArrayOf##BASETYPE *>(data) -> nelem();                                                                       \
}                                                                                                                                   \
                                                                                                                                    \
void resizeArrayOf##BASETYPE(Index n, void * data)                                                                                  \
{                                                                                                                                   \
    static_cast<ArrayOf##BASETYPE *>(data) -> resize(n);                                                                            \
}                                                                                                                                   \
                                                                                                                                    \
void * getelemArrayOf##BASETYPE(Index i, void * data)                                                                               \
{                                                                                                                                   \
    if (i >= 0 and i < static_cast<ArrayOf##BASETYPE *>(data) -> nelem())                                                           \
        return &static_cast<ArrayOf##BASETYPE *>(data) -> operator[](i);                                                            \
    else                                                                                                                            \
        return nullptr;                                                                                                             \
}                                                                                                                                   \
                                                                                                                                    \
Index xmlreadArrayOf##BASETYPE(void * data, char * filepath)                                                                        \
{                                                                                                                                   \
    try {                                                                                                                           \
        xml_read_from_file(filepath, *static_cast<ArrayOf##BASETYPE *>(data), Verbosity());                                         \
        return EXIT_SUCCESS;                                                                                                        \
    } catch (std::runtime_error& e) {                                                                                               \
        return EXIT_FAILURE;                                                                                                        \
    }                                                                                                                               \
}                                                                                                                                   \
                                                                                                                                    \
Index xmlsaveArrayOf##BASETYPE(void * data, char * filepath, Index filetype, Index clobber)                                         \
{                                                                                                                                   \
    try {                                                                                                                           \
        xml_write_to_file(filepath, *static_cast<const ArrayOf##BASETYPE *>(data), FileType(filetype), not clobber, Verbosity());   \
        return EXIT_SUCCESS;                                                                                                        \
    } catch (std::runtime_error& e) {                                                                                               \
        return EXIT_FAILURE;                                                                                                        \
    }                                                                                                                               \
}

WorkspaceArrayInterfaceCAPI_CCFILE(AbsorptionLines)
WorkspaceArrayInterfaceCAPI_CCFILE(ArrayOfAbsorptionLines)

#undef WorkspaceArrayInterfaceCAPI_CCFILE

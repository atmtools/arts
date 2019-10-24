/* Copyright (C) 2019
 Richard Larsson <larsson@mps.mpg.de>

 
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

/** Contains the absorption namespace
 * @file   absorptionlines.h
 * @author Richard Larsson
 * @date   2019-09-07
 * 
 * @brief  Contains the absorption lines implementation
 * 
 * This namespace contains classes to deal with absorption lines
 **/

#ifndef absorptionlines_h
#define absorptionlines_h

#include <vector>
#include "lineshapemodel.h"
#include "matpack.h"
#include "quantum.h"
#include "zeemandata.h"

/** Namespace to contain things required for absorption calculations */
namespace Absorption {
/** Describes the type of mirroring line effects
 * 
 * Each type but None has to have an implemented effect
 */
enum class MirroringType {
  None,             // No mirroring
  Lorentz,          // Mirror, but use Lorentz line shape
  SameAsLineShape,  // Mirror using the same line shape
  Manual,           // Mirror by having a line in the array of line record with negative F0
};  // MirroringType

inline MirroringType string2mirroringtype(const String& in) {
  if (in == "None")
    return MirroringType::None;
  else if (in == "Lorentz")
    return MirroringType::Lorentz;
  else if (in == "Same")
    return MirroringType::SameAsLineShape;
  else if (in == "Manual")
    return MirroringType::Manual;
  else
    throw std::runtime_error("Cannot recognize the mirroring type");
}

inline String mirroringtype2string(MirroringType in) {
  if (in == MirroringType::None)
    return "None";
  else if (in == MirroringType::Lorentz)
    return "Lorentz";
  else if (in == MirroringType::SameAsLineShape)
    return "Same";
  else if (in == MirroringType::Manual)
    return "Manual";
  std::terminate();
}

inline String mirroringtype2metadatastring(MirroringType in) {
  if (in == MirroringType::None)
    return "These lines are not mirrored at 0 Hz.\n";
  else if (in == MirroringType::Lorentz)
    return "These lines are mirrored around 0 Hz using the Lorentz line shape for the f0<0 mirrors.\n";
  else if (in == MirroringType::SameAsLineShape)
    return "These line are mirrored around 0 Hz using the original line shape for the f0<0 mirrors.\n";
  else if (in == MirroringType::Manual)
    return "There are manual line entires in the catalog at -f0 to mirror this line.\n";
  std::terminate();
}

/** Describes the type of normalization line effects
 *
 * Each type but None has to have an implemented effect
 */
enum class NormalizationType {
  None,  // Do not renormalize the line shape
  VVH,   // Renormalize with Van Vleck and Huber specifications
  VVW,   // Renormalize with Van Vleck and Weiskopf specifications
  RosenkranzQuadratic,  // Renormalize using Rosenkranz's quadratic specifications
};  // LineNormalizationType

inline NormalizationType string2normalizationtype(const String& in) {
  if (in == "None")
    return NormalizationType::None;
  else if (in == "VVH")
    return NormalizationType::VVH;
  else if (in == "VVW")
    return NormalizationType::VVW;
  else if (in == "RQ")
    return NormalizationType::RosenkranzQuadratic;
  else
    throw std::runtime_error("Cannot recognize the normalization type");
}

inline String normalizationtype2string(NormalizationType in) {
  if (in == NormalizationType::None)
    return "None";
  else if (in == NormalizationType::VVH)
    return "VVH";
  else if (in == NormalizationType::VVW)
    return "VVW";
  else if (in == NormalizationType::RosenkranzQuadratic)
    return "RQ";
  std::terminate();
}

inline String normalizationtype2metadatastring(NormalizationType in) {
  if (in == NormalizationType::None)
    return "No re-normalization in the far wing will be applied.\n";
  else if (in == NormalizationType::VVH)
    return "van Vleck and Huber far-wing renormalization will be applied, "
      "i.e. F ~ (f tanh(hf/2kT))/(f0 tanh(hf0/2kT))\n";
  else if (in == NormalizationType::VVW)
    return "van Vleck and Weisskopf far-wing renormalization will be applied, "
      "i.e. F ~ (f/f0)^2\n";
  else if (in == NormalizationType::RosenkranzQuadratic)
    return "Rosenkranz quadratic far-wing renormalization will be applied, "
      "i.e. F ~ hf0/2kT sinh(hf0/2kT) (f/f0)^2\n";
  std::terminate();
}

/** Describes the type of population level counter
 *
 * The types here might require that different data is available at runtime absorption calculations
 */
enum class PopulationType {
  ByLTE,                          // Assume line is in LTE
  ByRelmatMendazaLTE,             // Assume line is in LTE but requires Relaxation matrix calculations - follows Mendaza method
  ByRelmatHartmannLTE,            // Assume line is in LTE but requires Relaxation matrix calculations - follows Hartmann method
  ByNLTEVibrationalTemperatures,  // Assume line is in NLTE described by vibrational temperatures
  ByNLTEPopulationDistribution,   // Assume line is in NLTE and the upper-to-lower ratio is known
};  // PopulationType

inline PopulationType string2populationtype(const String& in) {
  if (in == "LTE")
    return PopulationType::ByLTE;
  if (in == "MendazaRelmat")
    return PopulationType::ByRelmatMendazaLTE;
  if (in == "HartmannRelmat")
    return PopulationType::ByRelmatHartmannLTE;
  else if (in == "NLTE-VibrationalTemperatures")
    return PopulationType::ByNLTEVibrationalTemperatures;
  else if (in == "NLTE")
    return PopulationType::ByNLTEPopulationDistribution;
  else
    throw std::runtime_error("Cannot recognize the population type");
}

inline String populationtype2string(PopulationType in) {
  switch (in) {
    case PopulationType::ByLTE:
      return "LTE";
    case PopulationType::ByRelmatMendazaLTE:
      return "MendazaRelmat";
    case PopulationType::ByRelmatHartmannLTE:
      return "HartmannRelmat";
    case PopulationType::ByNLTEVibrationalTemperatures:
      return "NLTE-VibrationalTemperatures";
    case PopulationType::ByNLTEPopulationDistribution:
      return "NLTE";
  } std::terminate();
}

inline String populationtype2metadatastring(PopulationType in) {
  switch (in) {
    case PopulationType::ByLTE:
      return "The lines are considered as in pure LTE.\n";
    case PopulationType::ByRelmatMendazaLTE:
      return "The lines requires Relaxation matrix calculations in LTE - Mendaza method.\n";
    case PopulationType::ByRelmatHartmannLTE:
      return "The lines requires Relaxation matrix calculations in LTE - Hartmann method.\n";
    case PopulationType::ByNLTEVibrationalTemperatures:
      return "The lines are considered as in NLTE by vibrational temperatures.\n";
    case PopulationType::ByNLTEPopulationDistribution:
      return "The lines are considered as in pure NLTE.\n";
  } std::terminate();
}

inline bool relaxationtype_relmat(PopulationType in) {
  return in == PopulationType::ByRelmatMendazaLTE or
         in == PopulationType::ByRelmatHartmannLTE;
}

/** Describes the type of cutoff calculations */
enum class CutoffType {
  None,                // No cutoff frequency at all
  LineByLineOffset,    // The cutoff frequency is at SingleLine::F0 plus the cutoff frequency
  BandFixedFrequency,  // The curoff frequency is the cutoff frequency for all SingleLine(s)
};  // LineCutoffType

inline CutoffType string2cutofftype(const String& in) {
  if (in == "None")
    return CutoffType::None;
  else if (in == "ByLine")
    return CutoffType::LineByLineOffset;
  else if (in == "ByBand")
    return CutoffType::BandFixedFrequency;
  else
    throw std::runtime_error("Cannot recognize the cutoff type");
}

inline String cutofftype2string(CutoffType in) {
  if (in == CutoffType::None)
    return "None";
  else if (in == CutoffType::LineByLineOffset)
    return "ByLine";
  else if (in == CutoffType::BandFixedFrequency)
    return "ByBand";
  std::terminate();
}

inline String cutofftype2metadatastring(CutoffType in, Numeric cutoff) {
  std::ostringstream os;
  if (in == CutoffType::None)
    os << "No cut-off will be applied.\n";
  else if (in == CutoffType::LineByLineOffset)
    os << "The lines will be cut-off " << cutoff << " Hz from the line center.\n";
  else if (in == CutoffType::BandFixedFrequency)
    os << "All lines are cut-off at " << cutoff << " Hz.\n";
  return os.str();
}

/** Computations and data for a single absorption line */
class SingleLine {
private:
  /** Central frequency */
  Numeric mF0;
  
  /** Reference intensity */
  Numeric mI0;
  
  /** Lower state energy level */
  Numeric mE0;
  
  /** Lower level statistical weight */
  Numeric mglow;
  
  /** Upper level statistical weight */
  Numeric mgupp;
  
  /** Einstein spontaneous emission coefficient */
  Numeric mA;
  
  /** Zeeman model */
  Zeeman::Model mzeeman;
  
  /** Line shape model */
  LineShape::Model2 mlineshape;
  
  /** Lower level quantum numbers */
  std::vector<Rational> mlowerquanta;
  
  /** Upper level quantum numbers */
  std::vector<Rational> mupperquanta;

public:
  /** Default initialization 
   * 
   * @param[in] F0 Central frequency
   * @param[in] I0 Reference line strength at external T0
   * @param[in] E0 Lower energy level
   * @param[in] glow Lower level statistical weight
   * @param[in] gupp Upper level statistical weight
   * @param[in] A Einstein spontaneous emission coefficient
   * @param[in] zeeman Zeeman model
   * @param[in] lineshape Line shape model
   * @param[in] lowerquanta Lower quantum numbers
   * @param[in] upperquanta Upper quantum numbers
   */
  SingleLine(Numeric F0=0,
             Numeric I0=0,
             Numeric E0=0,
             Numeric glow=0,
             Numeric gupp=0,
             Numeric A=0,
             Zeeman::Model zeeman=Zeeman::Model(),
             const LineShape::Model2& lineshape=LineShape::Model2(),
             const std::vector<Rational>& lowerquanta={},
             const std::vector<Rational>& upperquanta={}) :
             mF0(F0),
             mI0(I0),
             mE0(E0),
             mglow(glow),
             mgupp(gupp),
             mA(A),
             mzeeman(zeeman),
             mlineshape(lineshape),
             mlowerquanta(lowerquanta),
             mupperquanta(upperquanta) {}
  
  /** Initialization for constant sizes
   * 
   * @param[in] nbroadeners Number of broadening species
   * @param[in] nquanta Number of local quantum numbers
   */
  SingleLine(size_t nbroadeners, size_t nquanta, const LineShape::Model2& metamodel) :
  mlineshape(metamodel), mlowerquanta(nquanta), mupperquanta(nquanta) {
    if(Index(nbroadeners) not_eq mlineshape.nelem())
      throw std::runtime_error("Mismatch between broadeners and model");
  }
  
  //////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Counts
  //////////////////////////////////////////////////////////////////
  
  /** Number of lineshape elements */
  Index LineShapeElems() const noexcept {return mlineshape.nelem();}
  
  /** Number of lower quantum numbers */
  Index LowerQuantumElems() const noexcept {return mlowerquanta.size();}
  
  /** Number of upper quantum numbers */
  Index UpperQuantumElems() const noexcept {return mupperquanta.size();}
  
  //////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////// Constant access
  //////////////////////////////////////////////////////////////////
  
  /** Central frequency */
  Numeric F0() const noexcept {return mF0;}
  
  /** Lower level energy */
  Numeric E0() const noexcept {return mE0;}
  
  /** Reference line strength */
  Numeric I0() const noexcept {return mI0;}
  
  /** Einstein spontaneous emission */
  Numeric A() const noexcept {return mA;}
  
  /** Lower level statistical weight */
  Numeric g_low() const noexcept {return mglow;}
  
  /** Upper level statistical weight */
  Numeric g_upp() const noexcept {return mgupp;}
  
  /** Zeeman model */
  Zeeman::Model Zeeman() const noexcept {return mzeeman;}
  
  /** Line shape model */
  const LineShape::Model2& LineShape() const noexcept {return mlineshape;}
  
  /** Lower level quantum numbers */
  const std::vector<Rational>& LowerQuantumNumbers() const noexcept {return mlowerquanta;}
  
  /** Upper level quantum numbers */
  const std::vector<Rational>& UpperQuantumNumbers() const noexcept {return mupperquanta;}
  
  //////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////// Reference access
  //////////////////////////////////////////////////////////////////
  
  /** Central frequency */
  Numeric& F0() noexcept {return mF0;}
  
  /** Lower level energy */
  Numeric& E0() noexcept {return mE0;}
  
  /** Reference line strength */
  Numeric& I0() noexcept {return mI0;}
  
  /** Einstein spontaneous emission */
  Numeric& A() noexcept {return mA;}
  
  /** Lower level statistical weight */
  Numeric& g_low() noexcept {return mglow;}
  
  /** Upper level statistical weight */
  Numeric& g_upp() noexcept {return mgupp;}
  
  /** Zeeman model */
  Zeeman::Model& Zeeman() noexcept {return mzeeman;}
  
  /** Line shape model */
  LineShape::Model2& LineShape() noexcept {return mlineshape;}
  
  /** Lower level quantum numbers */
  std::vector<Rational>& LowerQuantumNumbers() noexcept {return mlowerquanta;}
  
  /** Upper level quantum numbers */
  std::vector<Rational>& UpperQuantumNumbers() noexcept {return mupperquanta;}
  
  //////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////// Special access
  //////////////////////////////////////////////////////////////////
  
  /** Lower quantum number */
  Rational LowerQuantumNumber(size_t i) const noexcept {return mlowerquanta[i];}
  
  /** Upper quantum number */
  Rational UpperQuantumNumber(size_t i) const noexcept {return mupperquanta[i];}
  
  /** Lower quantum number */
  Rational& LowerQuantumNumber(size_t i) noexcept {return mlowerquanta[i];}
  
  /** Upper quantum number */
  Rational& UpperQuantumNumber(size_t i) noexcept {return mupperquanta[i];}
  
  /** Checks if the quantum numbers are the same of the two lines */
  bool SameQuantumNumbers(const SingleLine& sl) const noexcept;
  
  //////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////// Special settings
  //////////////////////////////////////////////////////////////////
  
  /** Set Zeeman effect by automatic detection
   * 
   * Will fail if the available and provided quantum numbers are bad
   * 
   * @param[in] qid Copy of the global identifier to fill by local numbers
   * @param[in] keys List of quantum number keys in this line's local quantum number lists
   */
  void SetAutomaticZeeman(QuantumIdentifier qid, const std::vector<QuantumNumberType>& keys) {
    for(size_t i=0; i<keys.size(); i++) {
      qid.LowerQuantumNumber(keys[i]) = mlowerquanta[i];
      qid.UpperQuantumNumber(keys[i]) = mupperquanta[i];
    }
    
    mzeeman = Zeeman::Model(qid);
  }
  
  /** Set the line mixing model to 2nd order
   * 
   * @param[in] d Data in 2nd order format
   */
  void SetLineMixing2SecondOrderData(const Vector& d) {
    mlineshape.SetLineMixingModel(
      LineShape::LegacyLineMixingData::vector2modellm(
        d, LineShape::LegacyLineMixingData::TypeLM::LM_2NDORDER)
      .Data()[0]);
  }
  
  /** Set the line mixing model to AER kind
   * 
   * @param[in] d Data in AER format
   */
  void SetLineMixing2AER(const Vector& d) {
    const LineShape::ModelParameters Y = {LineShape::TemperatureModel::LM_AER, d[4], d[5], d[6], d[7]};
    const LineShape::ModelParameters G = {LineShape::TemperatureModel::LM_AER, d[8], d[9], d[10], d[11]};
    for (auto& sm : mlineshape.Data()) {
      sm.Y() = Y;
      sm.G() = G;
    }
  }
};  // SingleLine

std::ostream& operator<<(std::ostream&, const SingleLine&);

std::istream& operator>>(std::istream&, SingleLine&);

/** Single line reading output */
struct SingleLineExternal {
  bool bad=true;
  bool selfbroadening=false;
  bool bathbroadening=false;
  CutoffType cutoff=CutoffType::None;
  MirroringType mirroring=MirroringType::None;
  PopulationType population=PopulationType::ByLTE;
  NormalizationType normalization=NormalizationType::None;
  LineShape::Type lineshapetype=LineShape::Type::DP;
  Numeric T0=0;
  Numeric cutofffreq=0;
  Numeric linemixinglimit=-1;
  QuantumIdentifier quantumidentity=QuantumIdentifier(QuantumIdentifier::TRANSITION, -1, -1);
  ArrayOfSpeciesTag species;
  SingleLine line;
};

class Lines {
private:
  /** Does the line broadening have self broadening */
  bool mselfbroadening;
  
  /** Does the line broadening have bath broadening */
  bool mbathbroadening;
  
  /** cutoff type, by band or by line */
  CutoffType mcutoff;
  
  /** Mirroring type */
  MirroringType mmirroring;
  
  /** Line population distribution */
  PopulationType mpopulation;
  
  /** Line normalization type */
  NormalizationType mnormalization;

  /** Type of line shape */
  LineShape::Type mlineshapetype;
  
  /** Reference temperature for all parameters of the lines */
  Numeric mT0;
  
  /** cutoff frequency */
  Numeric mcutofffreq;
  
  /** linemixing limit */
  Numeric mlinemixinglimit;
  
  /** Catalog ID */
  QuantumIdentifier mquantumidentity;
  
  /** List of local quantum numbers, these must be defined */
  std::vector<QuantumNumberType> mlocalquanta;
  
  /** A list of broadening species */
  ArrayOfSpeciesTag mbroadeningspecies;
  
  /** A list of individual lines */
  std::vector<SingleLine> mlines;
  
public:
  /** Default initialization
   * 
   * @param[in] selfbroadening Do self broadening
   * @param[in] bathbroadening Do bath broadening
   * @param[in] cutoff Type of cutoff frequency
   * @param[in] mirroring Type of mirroring
   * @param[in] population Type of line strengths distributions
   * @param[in] normalization Type of normalization
   * @param[in] lineshapetype Type of line shape
   * @param[in] T0 Reference temperature
   * @param[in] cutofffreq Cutoff frequency
   * @param[in] linemixinglimit Line mixing limit
   * @param[in] quantumidentity Identity of global lines
   * @param[in] localquanta List of local quantum number(s)
   * @param[in] broadeningspecies List of broadening species
   * @param[in] lines List of SingleLine(s)
   */
  Lines(bool selfbroadening=false,
        bool bathbroadening=false,
        CutoffType cutoff=CutoffType::None,
        MirroringType mirroring=MirroringType::None,
        PopulationType population=PopulationType::ByLTE,
        NormalizationType normalization=NormalizationType::None,
        LineShape::Type lineshapetype=LineShape::Type::DP,
        Numeric T0=296,
        Numeric cutofffreq=-1,
        Numeric linemixinglimit=-1,
        const QuantumIdentifier& quantumidentity=QuantumIdentifier(),
        const std::vector<QuantumNumberType>& localquanta={},
        const ArrayOfSpeciesTag& broadeningspecies={},
        const std::vector<SingleLine>& lines={}) :
        mselfbroadening(selfbroadening),
        mbathbroadening(bathbroadening),
        mcutoff(cutoff),
        mmirroring(mirroring),
        mpopulation(population),
        mnormalization(normalization),
        mlineshapetype(lineshapetype),
        mT0(T0),
        mcutofffreq(cutofffreq),
        mlinemixinglimit(linemixinglimit),
        mquantumidentity(quantumidentity),
        mlocalquanta(localquanta),
        mbroadeningspecies(broadeningspecies),
        mlines(lines) {};
  
  /** XML-tag initialization
   * 
   * @param[in] selfbroadening Do self broadening
   * @param[in] bathbroadening Do bath broadening
   * @param[in] nlines Number of SingleLine(s) to initiate as empty
   * @param[in] cutoff Type of cutoff frequency
   * @param[in] mirroring Type of mirroring
   * @param[in] population Type of line strengths distributions
   * @param[in] normalization Type of normalization
   * @param[in] lineshapetype Type of line shape
   * @param[in] T0 Reference temperature
   * @param[in] cutofffreq Cutoff frequency
   * @param[in] linemixinglimit Line mixing limit
   * @param[in] quantumidentity Identity of global lines
   * @param[in] localquanta List of local quantum number(s)
   * @param[in] broadeningspecies List of broadening species
   * @param[in] metamodel A line shape model with defined shapes
   */
  Lines(bool selfbroadening,
        bool bathbroadening,
        size_t nlines,
        CutoffType cutoff,
        MirroringType mirroring,
        PopulationType population,
        NormalizationType normalization,
        LineShape::Type lineshapetype,
        Numeric T0,
        Numeric cutofffreq,
        Numeric linemixinglimit,
        const QuantumIdentifier& quantumidentity,
        const std::vector<QuantumNumberType>& localquanta,
        const ArrayOfSpeciesTag& broadeningspecies,
        const LineShape::Model2& metamodel) :
        mselfbroadening(selfbroadening),
        mbathbroadening(bathbroadening),
        mcutoff(cutoff),
        mmirroring(mirroring),
        mpopulation(population),
        mnormalization(normalization),
        mlineshapetype(lineshapetype),
        mT0(T0),
        mcutofffreq(cutofffreq),
        mlinemixinglimit(linemixinglimit),
        mquantumidentity(quantumidentity),
        mlocalquanta(localquanta),
        mbroadeningspecies(broadeningspecies),
        mlines(nlines,
               SingleLine(broadeningspecies.size(),
               localquanta.size(), metamodel)) {};
  
  /** Appends a single line to the absorption lines
   * 
   * Useful for reading undefined number of lines and setting
   * their structures
   * 
   * Warning: caller must guarantee that the broadening species
   * and the quantum numbers of both levels have the correct
   * order and the correct size.  Only the sizes can be and are
   * tested.
   * 
   * @param[in] sl A single line
   */
  void AppendSingleLine(SingleLine&& sl) {
    if(NumBroadeners() not_eq sl.LowerQuantumElems() or
       NumBroadeners() not_eq sl.UpperQuantumElems())
      throw std::runtime_error("Error calling appending function, bad size of quantum numbers");
    
    if(NumLines() not_eq 0 and 
       sl.LineShapeElems() not_eq mlines[0].LineShapeElems())
      throw std::runtime_error("Error calling appending function, bad size of broadening species");
    
    mlines.push_back(std::move(sl));
  }
  
  /** Appends a single line to the absorption lines
   * 
   * Useful for reading undefined number of lines and setting
   * their structures
   * 
   * Warning: caller must guarantee that the broadening species
   * and the quantum numbers of both levels have the correct
   * order and the correct size.  Only the sizes can be and are
   * tested.
   * 
   * @param[in] sl A single line
   */
  void AppendSingleLine(const SingleLine& sl) {
    if(NumBroadeners() not_eq sl.LowerQuantumElems() or
       NumBroadeners() not_eq sl.UpperQuantumElems())
      throw std::runtime_error("Error calling appending function, bad size of quantum numbers");
    
    if(NumLines() not_eq 0 and 
       sl.LineShapeElems() not_eq mlines[0].LineShapeElems())
      throw std::runtime_error("Error calling appending function, bad size of broadening species");
    
    mlines.push_back(sl);
  }
  
  /** Checks if an external line matches this structure
   * 
   * @param[in] sle Full external lines
   * @param[in] quantumidentity Expected global quantum id of the line
   */
  bool MatchWithExternal(const SingleLineExternal& sle, const QuantumIdentifier& quantumidentity) const noexcept {
    if(sle.bad)
      return false;
    else if(sle.selfbroadening not_eq mselfbroadening)
      return false;
    else if(sle.bathbroadening not_eq mbathbroadening)
      return false;
    else if(sle.cutoff not_eq mcutoff)
      return false;
    else if(sle.mirroring not_eq mmirroring)
      return false;
    else if(sle.population not_eq mpopulation)
      return false;
    else if(sle.normalization not_eq mnormalization)
      return false;
    else if(sle.lineshapetype not_eq mlineshapetype)
      return false;
    else if(sle.T0 not_eq mT0)
      return false;
    else if(sle.cutofffreq not_eq mcutofffreq)
      return false;
    else if(sle.linemixinglimit not_eq mlinemixinglimit)
      return false;
    else if(quantumidentity not_eq mquantumidentity)
      return false;
    else if(not std::equal(sle.species.cbegin(), sle.species.cend(), mbroadeningspecies.cbegin(), mbroadeningspecies.cend()))
      return false;
    else if(NumLines() not_eq 0 and not sle.line.LineShape().Match(mlines[0].LineShape()))
      return false;
    else
      return true;
  }
  
  /** Checks if another line list matches this structure
   * 
   * @param[in] sle Full external lines
   * @param[in] quantumidentity Expected global quantum id of the line
   */
  bool Match(const Lines& l) const noexcept {
    if(l.mselfbroadening not_eq mselfbroadening)
      return false;
    else if(l.mbathbroadening not_eq mbathbroadening)
      return false;
    else if(l.mcutoff not_eq mcutoff)
      return false;
    else if(l.mmirroring not_eq mmirroring)
      return false;
    else if(l.mpopulation not_eq mpopulation)
      return false;
    else if(l.mnormalization not_eq mnormalization)
      return false;
    else if(l.mlineshapetype not_eq mlineshapetype)
      return false;
    else if(l.mT0 not_eq mT0)
      return false;
    else if(l.mcutofffreq not_eq mcutofffreq)
      return false;
    else if(l.mlinemixinglimit not_eq mlinemixinglimit)
      return false;
    else if(l.mquantumidentity not_eq mquantumidentity)
      return false;
    else if(not std::equal(l.mbroadeningspecies.cbegin(), l.mbroadeningspecies.cend(), mbroadeningspecies.cbegin(), mbroadeningspecies.cend()))
      return false;
    else if(not std::equal(l.mlocalquanta.cbegin(), l.mlocalquanta.cend(), mlocalquanta.cbegin(), mlocalquanta.cend()))
      return false;
    else if(NumLines() not_eq 0 and l.NumLines() not_eq 0 and not l.mlines[0].LineShape().Match(mlines[0].LineShape()))
      return false;
    else
      return true;
  }
  
  /** Sort inner line list by frequency */
  void sort_by_frequency() {
    std::sort(mlines.begin(), mlines.end(),
              [](const SingleLine& a, const SingleLine& b){return a.F0() < b.F0();});
  }
  
  /** Sort inner line list by Einstein coefficient */
  void sort_by_einstein() {
    std::sort(mlines.begin(), mlines.end(),
              [](const SingleLine& a, const SingleLine& b){return a.A() < b.A();});
  }
  
  /** Removes all global quantum numbers */
  void truncate_global_quantum_numbers() {
    mquantumidentity.SetTransition(QuantumNumbers(), QuantumNumbers());
  }
  
  /** Species Name */
  String SpeciesName() const noexcept;
  
  /** Upper quantum numbers string */
  String UpperQuantumNumbers() const noexcept;
  
  /** Lower quantum numbers string */
  String LowerQuantumNumbers() const noexcept;
  
  /** Meta data for the line shape if it exists */
  String LineShapeMetaData() const noexcept {
    return NumLines() ?
      LineShape::ModelShape2MetaData(mlines[0].LineShape()) :
      "";
  }
  
  /** Species Index */
  Index Species() const noexcept {return mquantumidentity.Species();}
  
  /** Isotopologue Index */
  Index Isotopologue() const noexcept {return mquantumidentity.Isotopologue();}
  
  /** Number of lines */
  Index NumLines() const noexcept {return Index(mlines.size());}
  
  /** Lines */
  const std::vector<SingleLine>& AllLines() const noexcept {return mlines;}
  
  /** Lines */
  std::vector<SingleLine>& AllLines() noexcept {return mlines;}
  
  /** Number of broadening species */
  Index NumBroadeners() const noexcept {return Index(mlocalquanta.size());}
  
  /** Remove quantum numbers that are not used by even a single line
   */
  void RemoveUnusedLocalQuantums();
  
  /** Remove quantum numbers at the given position from all lines 
   */
  void RemoveLocalQuantum(size_t);
  
  /** Quantum number lower level
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] qnt Quantum number type
   * @return Quantum number
   */
  Rational LowerQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept;
  
  /** Quantum number upper level
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] qnt Quantum number type
   * @return Quantum number
   */
  Rational UpperQuantumNumber(size_t k, QuantumNumberType qnt) const noexcept;
  
  /** Checks if all defined quantum numbers in qid are defined and equal to 
   * the k:th line's local quantum numbers
   * 
   * qid must be a transition or allow for all matches
   * 
   * @param[in] qid Energy level quantum identifier
   * @param[in] k Line number (less than NumLines())
   * @return true if they are
   */
  bool LocalIn(const QuantumIdentifier& qid, size_t k) const noexcept;
  
  /** Checks if all defined quantum numbers in qid are defined and equal to 
   * the k:th line's local quantum numbers
   * 
   * qid must be a transition or all-allowing type to return true
   * 
   * @param[in] qid Energy level quantum identifier
   * @param[in] k Line number (less than NumLines())
   * @return true if they are
   */
  bool InLocal(const QuantumIdentifier& qid, size_t k) const noexcept;
  
  /** Checks if all defined local quantum numbers of the k:th line are equal to 
   * qid's quantum numbers
   * 
   * qid must be a transition or all-allowing type to return true
   * 
   * @param[in] qid Energy level quantum identifier
   * @param[in] k Line number (less than NumLines())
   * @return true if they are
   */
  bool InLocalUpper(const QuantumIdentifier& qid, size_t k) const noexcept;
  
  /** Checks if all defined upper local quantum numbers of the k:th line are equal to 
   * qid's quantum numbers
   * 
   * qid is either checked for its upper level or as an energy level
   * 
   * @param[in] qid Energy level quantum identifier
   * @param[in] k Line number (less than NumLines())
   * @return true if they are
   */
  bool LocalUpperIn(const QuantumIdentifier& qid, size_t k) const noexcept;
  
  /** Checks if all defined lower local quantum numbers of the k:th line are equal to 
   * qid's quantum numbers
   * 
   * qid is either checked for its lower level or as an energy level
   * 
   * @param[in] qid Energy level quantum identifier
   * @param[in] k Line number (less than NumLines())
   * @return true if they are
   */
  bool LocalLowerIn(const QuantumIdentifier& qid, size_t k) const noexcept;
  
  /** Checks if all defined quantum numbers in qid are defined and equal to 
   * the k:th line's local lower quantum numbers
   * 
   * qid is either checked for its lower level or as an energy level
   * 
   * @param[in] qid Energy level quantum identifier
   * @param[in] k Line number (less than NumLines())
   * @return true if they are
   */
  bool InLocalLower(const QuantumIdentifier& qid, size_t k) const noexcept;
  
  /** Returns the number of Zeeman split lines
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] type Type of Zeeman polarization
   */
  Index ZeemanCount(size_t k, Zeeman::Polarization type) const noexcept {
    return Zeeman::nelem(UpperQuantumNumber(k, QuantumNumberType::J),
                         LowerQuantumNumber(k, QuantumNumberType::J),
                         type);
  }
  
  /** Returns the strength of a Zeeman split line
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] type Type of Zeeman polarization
   * @param[in] i Zeeman line count
   */
  Numeric ZeemanStrength(size_t k, Zeeman::Polarization type, Index i) const noexcept {
    return mlines[k].Zeeman().Strength(UpperQuantumNumber(k, QuantumNumberType::J),
                                       LowerQuantumNumber(k, QuantumNumberType::J),
                                       type, i);
  }
  
  /** Returns the splitting of a Zeeman split line
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] type Type of Zeeman polarization
   * @param[in] i Zeeman line count
   */
  Numeric ZeemanSplitting(size_t k, Zeeman::Polarization type, Index i) const noexcept {
    return mlines[k].Zeeman().Splitting(UpperQuantumNumber(k, QuantumNumberType::J),
                                        LowerQuantumNumber(k, QuantumNumberType::J),
                                        type, i);
  }
  
  /** Set Zeeman effect for all lines that have the correct quantum numbers */
  void SetAutomaticZeeman() noexcept {
    for(auto& line: mlines)
      line.SetAutomaticZeeman(mquantumidentity, mlocalquanta);
  }
  
  /** Central frequency
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Central frequency
   */
  Numeric F0(size_t k) const noexcept {return mlines[k].F0();}
  
  /** Central frequency
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Central frequency
   */
  Numeric& F0(size_t k) noexcept {return mlines[k].F0();}
  
  /** Mean frequency
   * 
   * @return Mean frequency
   */
  Numeric F_mean() const noexcept {
    Numeric x = 0, div = 1 / Numeric(NumLines());
    for (auto& line: mlines)
      x = std::fma(div, line.F0(), x);
    return x;
  }
  
  /** Lower level energy
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Lower level energy
   */
  Numeric E0(size_t k) const noexcept {return mlines[k].E0();}
  
  /** Lower level energy
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Lower level energy
   */
  Numeric& E0(size_t k) noexcept {return mlines[k].E0();}
  
  /** Reference line strength
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Reference line strength
   */
  Numeric I0(size_t k) const noexcept {return mlines[k].I0();}
  
  /** Reference line strength
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Reference line strength
   */
  Numeric& I0(size_t k) noexcept {return mlines[k].I0();}
  
  /** Einstein spontaneous emission
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Einstein spontaneous emission
   */
  Numeric A(size_t k) const noexcept {return mlines[k].A();}
  
  /** Einstein spontaneous emission
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Einstein spontaneous emission
   */
  Numeric& A(size_t k) noexcept {return mlines[k].A();}
  
  /** Lower level statistical weight
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Lower level statistical weight
   */
  Numeric g_low(size_t k) const noexcept {return mlines[k].g_low();}
  
  /** Lower level statistical weight
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Lower level statistical weight
   */
  Numeric& g_low(size_t k) noexcept {return mlines[k].g_low();}
  
  /** Upper level statistical weight
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Upper level statistical weight
   */
  Numeric g_upp(size_t k) const noexcept {return mlines[k].g_upp();}
  
  /** Upper level statistical weight
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Upper level statistical weight
   */
  Numeric& g_upp(size_t k) noexcept {return mlines[k].g_upp();}
  
  /** Returns mirroring style */
  MirroringType Mirroring() const noexcept {return mmirroring;}
  
  /** Returns normalization style */
  NormalizationType Normalization() const noexcept {return mnormalization;}
  
  /** Returns cutoff style */
  CutoffType Cutoff() const noexcept {return mcutoff;}
  
  /** Returns population style */
  PopulationType Population() const noexcept {return mpopulation;}
  
  /** Returns population style */
  LineShape::Type LineShapeType() const noexcept {return mlineshapetype;}
  
  /** Returns mirroring style */
  void Mirroring(MirroringType x) noexcept {mmirroring = x;}
  
  /** Returns normalization style */
  void Normalization(NormalizationType x) noexcept {mnormalization = x;}
  
  /** Returns cutoff style */
  void Cutoff(CutoffType x) noexcept {mcutoff = x;}
  
  /** Returns population style */
  void Population(PopulationType x) noexcept {mpopulation = x;}
  
  /** Returns population style */
  void LineShapeType(LineShape::Type x) noexcept {mlineshapetype = x;}
  
  /** Returns if the pressure should do line mixing
   * 
   * @param[in] P Atmospheric pressure
   * @return true if no limit or P less than limit
   */
  bool DoLineMixing(Numeric P) const noexcept {
    return mlinemixinglimit < 0 ? true : mlinemixinglimit > P;
  }
  
  /** Line shape parameters
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] vmrs Line broadener species's volume mixing ratio
   * @return Line shape parameters
   */
  LineShape::Output ShapeParameters(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept;
  
  /** Line shape parameters temperature derivatives
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] vmrs Line broadener's volume mixing ratio
   * @return Line shape parameters temperature derivatives
   */
  LineShape::Output ShapeParameters_dT(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept;
  
  /** Position among broadening species or -1
   * 
   * @param[in] A species index that might be among the broadener species
   * @return Position among broadening species or -1
   */
  Index LineShapePos(const Index& spec) const noexcept;
  
  /** Position among broadening species or -1
   * 
   * @param[in] An identity that might be among the broadener species
   * @return Position among broadening species or -1
   */
  Index LineShapePos(const QuantumIdentifier& qid) const noexcept {
    return LineShapePos(qid.Species());
  }
  
  /** Line shape parameters vmr derivative
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] vmr_qid Identity of species whose VMR derivative is requested
   * @return Line shape parameters vmr derivative
   */
  LineShape::Output ShapeParameters_dVMR(size_t k, Numeric T, Numeric P,
                                         const QuantumIdentifier& vmr_qid) const noexcept;
  
  /** Line shape parameter internal derivative
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] vmrs Line broadener's volume mixing ratio
   * @param[in] derivative Type of line shape derivative
   * @return Line shape parameter internal derivative
   */
  Numeric ShapeParameter_dInternal(size_t k, Numeric T, Numeric P,
                                   const Vector& vmrs,
                                   const RetrievalQuantity& derivative) const noexcept;
  
  /** Returns cutoff frequency or maximum value
   * 
   * @param[in] k Line number (less than NumLines())
   * @returns Cutoff frequency or 0
   */
  Numeric CutoffFreq(size_t k) const noexcept {
    switch(mcutoff) {
      case CutoffType::LineByLineOffset:
        return F0(k) + mcutofffreq;
      case CutoffType::BandFixedFrequency:
        return mcutofffreq;
      case CutoffType::None:
        return std::numeric_limits<Numeric>::max();
    }
    std::terminate();
  }
  
  /** Returns negative cutoff frequency or lowest value
   * 
   * @param[in] k Line number (less than NumLines())
   * @returns Negative cutoff frequency or the lowest value
   */
  Numeric CutoffFreqMinus(size_t k, Numeric fmean) const noexcept {
    switch(mcutoff) {
      case CutoffType::LineByLineOffset:
        return F0(k) - mcutofffreq;
      case CutoffType::BandFixedFrequency:
        return mcutofffreq - 2*fmean;
      case CutoffType::None:
        return std::numeric_limits<Numeric>::lowest();
    }
    std::terminate();
  }
  
  /** Returns reference temperature */
  Numeric T0() const noexcept {
    return mT0;
  }
  
  /** Sets reference temperature */
  void T0(Numeric x) noexcept {
    mT0 = x;
  }
  
  /** Returns internal cutoff frequency value */
  Numeric CutoffFreqValue() const noexcept {
    return mcutofffreq;
  }
  
  /** Sets internal cutoff frequency value */
  void CutoffFreqValue(Numeric x) noexcept {
    mcutofffreq = x;
  }
  
  /** Returns line mixing limit */
  Numeric LinemixingLimit() const noexcept {
    return mlinemixinglimit;
  }
  
  /** Sets line mixing limit */
  void LinemixingLimit(Numeric x) noexcept {
    mlinemixinglimit = x;
  }
  
  /** Returns local quantum numbers */
  const std::vector<QuantumNumberType>& LocalQuanta() const noexcept {
    return mlocalquanta;
  }
  
  /** Returns the broadening species */
  const ArrayOfSpeciesTag& BroadeningSpecies() const noexcept {
    return mbroadeningspecies;
  }
  
  /** Returns self broadening status */
  bool Self() const noexcept {
    return mselfbroadening;
  }
  
  /** Returns bath broadening status */
  bool Bath() const noexcept {
    return mbathbroadening;
  }
  
  /** Returns identity status */
  const QuantumIdentifier& QuantumIdentity() const noexcept {
    return mquantumidentity;
  }
  
  /** Returns identity status */
  QuantumIdentifier QuantumIdentityOfLine(Index k) const noexcept {
    QuantumIdentifier copy(mquantumidentity);
    for (size_t i=0; i<mlocalquanta.size(); i++) {
      copy.UpperQuantumNumber(mlocalquanta[i]) = mlines[k].UpperQuantumNumber(i);
      copy.LowerQuantumNumber(mlocalquanta[i]) = mlines[k].LowerQuantumNumber(i);
    }
    return copy;
  }
  
  /** Returns a printable statement about the lines */
  String MetaData() const noexcept;
  
  /** Removes a single line */
  void RemoveLine(Index) noexcept;
  
  /** Pops a single line */
  SingleLine PopLine(Index) noexcept;
  
  /** Returns a single line */
  SingleLine& Line(Index) noexcept;
  
  /** Returns a single line */
  const SingleLine& Line(Index) const noexcept;
  
  /** Reverses the order of the internal lines */
  void ReverseLines() noexcept;
  
  /** Mass of the molecule */
  Numeric SpeciesMass() const noexcept;
  
  /** Returns the VMRs of the broadening species
   * 
   * @param[in] atm_vmrs Atmospheric VMRs
   * @param[in] atm_spec Atmospheric Species
   * @return VMR list of the species
   */
  Vector BroadeningSpeciesVMR(const ConstVectorView, const ArrayOfArrayOfSpeciesTag&) const;
  
  /** Returns the VMR of the species
   * 
   * @param[in] atm_vmrs Atmospheric VMRs
   * @param[in] atm_spec Atmospheric Species
   * @return VMR of the species
   */
  Numeric SelfVMR(const ConstVectorView, const ArrayOfArrayOfSpeciesTag&) const;
};  // Lines

std::ostream& operator<<(std::ostream&, const Lines&);
std::istream& operator>>(std::istream&, Lines&);

/** Read from ARTSCAT-3
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromArtscat3Stream(istream& is);

/** Read from ARTSCAT-4
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromArtscat4Stream(istream& is);

/** Read from ARTSCAT-5
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromArtscat5Stream(istream& is);

/** Read from LBLRTM
 * 
 * LBLRTM follows the old HITRAN format from before 2004.  This
 * HITRAN format is as follows (directly from the HITRAN documentation):
 *
 * @verbatim
  Each line consists of 100
  bytes of ASCII text data, followed by a line feed (ASCII 10) and
  carriage return (ASCII 13) character, for a total of 102 bytes per line.
  Each line can be read using the following READ and FORMAT statement pair
  (for a FORTRAN sequential access read):

        READ(3,800) MO,ISO,V,S,R,AGAM,SGAM,E,N,d,V1,V2,Q1,Q2,IERF,IERS,
       *  IERH,IREFF,IREFS,IREFH
  800   FORMAT(I2,I1,F12.6,1P2E10.3,0P2F5.4,F10.4,F4.2,F8.6,2I3,2A9,3I1,3I2)

  Each item is defined below, with its format shown in parenthesis.

    MO  (I2)  = molecule number
    ISO (I1)  = isotopologue number (1 = most abundant, 2 = second, etc)
    V (F12.6) = frequency of transition in wavenumbers (cm-1)
    S (E10.3) = intensity in cm-1/(molec * cm-2) at 296 Kelvin
    R (E10.3) = transition probability squared in Debyes**2
    AGAM (F5.4) = air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
    SGAM (F5.4) = self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
    E (F10.4) = lower state energy in wavenumbers (cm-1)
    N (F4.2) = coefficient of temperature dependence of air-broadened halfwidth
    d (F8.6) = shift of transition due to pressure (cm-1)
    V1 (I3) = upper state global quanta index
    V2 (I3) = lower state global quanta index
    Q1 (A9) = upper state local quanta
    Q2 (A9) = lower state local quanta
    IERF (I1) = accuracy index for frequency reference
    IERS (I1) = accuracy index for intensity reference
    IERH (I1) = accuracy index for halfwidth reference
    IREFF (I2) = lookup index for frequency
    IREFS (I2) = lookup index for intensity
    IREFH (I2) = lookup index for halfwidth

  The molecule numbers are encoded as shown in the table below:

    0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
    6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
    12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
    18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
    24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
    30=  SF6   31=  H2S   32=HCOOH
 * @endverbatim
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromLBLRTMStream(istream& is);

/** Read from newer HITRAN
 *
 * The HITRAN format is as follows:
 *
 * @verbatim
  Each line consists of 160 ASCII characters, followed by a line feed (ASCII 10)
  and carriage return (ASCII 13) character, for a total of 162 bytes per line.

  Each item is defined below, with its Fortran format shown in parenthesis.

  (I2)     molecule number
  (I1)     isotopologue number (1 = most abundant, 2 = second, etc)
  (F12.6)  vacuum wavenumbers (cm-1)
  (E10.3)  intensity in cm-1/(molec * cm-2) at 296 Kelvin
  (E10.3)  Einstein-A coefficient (s-1)
  (F5.4)   air-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
  (F5.4)   self-broadened halfwidth (HWHM) in cm-1/atm at 296 Kelvin
  (F10.4)  lower state energy (cm-1)
  (F4.2)   coefficient of temperature dependence of air-broadened halfwidth
  (F8.6)   air-broadened pressure shift of line transition at 296 K (cm-1)
  (A15)    upper state global quanta
  (A15)    lower state global quanta
  (A15)    upper state local quanta
  (A15)    lower state local quanta
  (I1)     uncertainty index for wavenumber
  (I1)     uncertainty index for intensity
  (I1)     uncertainty index for air-broadened half-width
  (I1)     uncertainty index for self-broadened half-width
  (I1)     uncertainty index for temperature dependence
  (I1)     uncertainty index for pressure shift
  (I2)     index for table of references correspond. to wavenumber
  (I2)     index for table of references correspond. to intensity
  (I2)     index for table of references correspond. to air-broadened half-width
  (I2)     index for table of references correspond. to self-broadened half-width
  (I2)     index for table of references correspond. to temperature dependence
  (I2)     index for table of references correspond. to pressure shift
  (A1)     flag (*) for lines supplied with line-coupling algorithm
  (F7.1)   upper state statistical weight
  (F7.1)   lower state statistical weight

  The molecule numbers are encoded as shown in the table below:

    0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=    CO
    6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=   NH3
    12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=    HI
    18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=   HCN
    24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29=  COF2
    30=  SF6   31=  H2S   32=HCOOH   33=  HO2   34=    O   35=ClONO2
    36=  NO+   37= HOBr   38= C2H4
 * @endverbatim
 * 
 * @param[in] is Input stream
 * @param[in] fmin Lowest frequency to continue reading
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromHitran2004Stream(istream& is);

/** Read from Mytran2
 * The MYTRAN2
 * format is as follows (directly taken from the abs_my.c documentation):
 *
 * @verbatim
  The MYTRAN format is as follows (FORTRAN notation):
  FORMAT(I2,I1,F13.4,1PE10.3,0P2F5.2,F10.4,2F4.2,F8.6,F6.4,2I3,2A9,4I1,3I2)
  
  Each item is defined below, with its FORMAT String shown in
  parenthesis.
  
      MO  (I2)      = molecule number
      ISO (I1)      = isotopologue number (1 = most abundant, 2 = second, etc)
  *  F (F13.4)     = frequency of transition in MHz
  *  errf (F8.4)   = error in f in MHz
      S (E10.3)     = intensity in cm-1/(molec * cm-2) at 296 K
  *  AGAM (F5.4)   = air-broadened halfwidth (HWHM) in MHz/Torr at Tref
  *  SGAM (F5.4)   = self-broadened halfwidth (HWHM) in MHz/Torr at Tref
      E (F10.4)     = lower state energy in wavenumbers (cm-1)
      N (F4.2)      = coefficient of temperature dependence of 
                      air-broadened halfwidth
  *  N_self (F4.2) = coefficient of temperature dependence of 
                      self-broadened halfwidth
  *  Tref (F7.2)   = reference temperature for AGAM and SGAM 
  *  d (F9.7)      = shift of transition due to pressure (MHz/Torr)
      V1 (I3)       = upper state global quanta index
      V2 (I3)       = lower state global quanta index
      Q1 (A9)       = upper state local quanta
      Q2 (A9)       = lower state local quanta
      IERS (I1)     = accuracy index for S
      IERH (I1)     = accuracy index for AGAM
  *  IERN (I1)     = accuracy index for N

  
  The asterisks mark entries that are different from HITRAN.

  Note that AGAM and SGAM are for the temperature Tref, while S is
  still for 296 K!
  
  The molecule numbers are encoded as shown in the table below:
  
     0= Null    1=  H2O    2=  CO2    3=   O3    4=  N2O    5=   CO
     6=  CH4    7=   O2    8=   NO    9=  SO2   10=  NO2   11=  NH3
    12= HNO3   13=   OH   14=   HF   15=  HCl   16=  HBr   17=   HI
    18=  ClO   19=  OCS   20= H2CO   21= HOCl   22=   N2   23=  HCN
    24=CH3Cl   25= H2O2   26= C2H2   27= C2H6   28=  PH3   29= COF2
    30=  SF6   31=  H2S   32=HCOOH   33= HO2    34=    O   35= CLONO2
    36=  NO+   37= Null   38= Null   39= Null   40=H2O_L   41= Null
    42= Null   43= OCLO   44= Null   45= Null   46=BRO     47= Null
    48= H2SO4  49=CL2O2

  All molecule numbers are from HITRAN, except for species with id's
  greater or equals 40, which are not included in HITRAN.
  (E.g.: For BrO, iso=1 is Br-79-O,iso=2 is  Br-81-O.)
 * @endverbatim
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromMytran2Stream(istream& is);

/** Read from JPL
 * 
 *  The JPL format is as follows (directly taken from the JPL documentation):
 * 
 * @verbatim 
    The catalog line files are composed of 80-character lines, with one
    line entry per spectral line.  The format of each line is:

    \label{lfmt}
    \begin{tabular}{@{}lccccccccr@{}}
    FREQ, & ERR, & LGINT, & DR, & ELO, & GUP, & TAG, & QNFMT, & QN${'}$, & QN${''}$\\ 
    (F13.4, & F8.4, & F8.4, & I2, & F10.4, & I3, & I7, & I4, & 6I2, & 6I2)\\
    \end{tabular}

    \begin{tabular}{lp{4.5in}} 
    FREQ: & Frequency of the line in MHz.\\ 
    ERR: & Estimated or experimental error of FREQ in MHz.\\ 
    LGINT: &Base 10 logarithm of the integrated intensity 
    in units of \linebreak nm$^2$$\cdot$MHz at 300 K. (See Section 3 for 
    conversions to other units.)\\ 
    DR: & Degrees of freedom in the rotational partition 
    function (0 for atoms, 2 for linear molecules, and 3 for nonlinear 
    molecules).\\ 
    ELO: &Lower state energy in cm$^{-1}$ relative to the lowest energy 
    spin--rotation level in ground vibronic state.\\ 
    GUP: & Upper state degeneracy.\\ 
    TAG: & Species tag or molecular identifier. 
    A negative value flags that the line frequency has 
    been measured in the laboratory.  The absolute value of TAG is then the 
    species tag and ERR is the reported experimental error.  The three most 
    significant digits of the species tag are coded as the mass number of the 
    species, as explained above.\\ 
    QNFMT: &Identifies the format of the quantum numbers 
    given in the field QN. These quantum number formats are given in Section 5 
    and are different from those in the first two editions of the catalog.\\ 
    QN${'}$: & Quantum numbers for the upper state coded 
    according to QNFMT.\\ 
    QN${''}$: & Quantum numbers for the lower state.\\
    \end{tabular} 
 * @endverbatim
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromJplStream(istream& is);

/** Splits a list of lines into proper Lines
 * 
 * Ensures that all but SingleLine list in Lines is the same in a full
 * Lines
 * 
 * @param[in] lines A list of lines
 * @param[in] localquantas List of quantum numbers to be presumed local
 * @param[in] globalquantas List of quantum numbers to be presumed global
 * @return A list of properly ordered Lines
 */
std::vector<Lines> split_list_of_external_lines(const std::vector<SingleLineExternal>& external_lines,
                                                const std::vector<QuantumNumberType>& localquantas={},
                                                const std::vector<QuantumNumberType>& globalquantas={});

/** Creates a copy of the input lines structure
 * 
 * The output will have zero lines but be otherwise a copy of the input
 * 
 * @param[in] al Lines which structure is copied
 */
Lines createEmptyCopy(const Lines& al) noexcept;

/** Checks if the external quantum identifier match a line's ID
 * 
 * The check demands that all defined quantum numbers for the line
 * are the same as for the id
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool line_in_id(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Checks if the external quantum identifier is equal to a line's identifier
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool line_is_id(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Checks if the external quantum identifier match a line's ID
 * 
 * The check demands that all defined quantum numbers for the id
 * are the same as for the line
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool id_in_line(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Checks if the external quantum identifier match a line's ID
 * 
 * The check demands that all defined upper quantum numbers for the line
 * are the same as for the id
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool line_upper_in_id(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Checks if the external quantum identifier match a line's ID
 * 
 * The check demands that all defined lower quantum numbers for the line
 * are the same as for the id
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool line_lower_in_id(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Checks if the external quantum identifier match a line's ID
 * 
 * The check demands that all defined quantum numbers for the id
 * are the same as for the line's upper numbers
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool id_in_line_upper(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Checks if the external quantum identifier match a line's ID
 * 
 * The check demands that all defined quantum numbers for the id
 * are the same as for the line's lower numbers
 * 
 * @param[in] band The band of lines
 * @param[in] id An identifier
 * @param[in] line_index The local line
 */
bool id_in_line_lower(const Lines& band, const QuantumIdentifier& id, size_t line_index);

/** Number of lines */
inline Index nelem(const Lines& l) {return l.NumLines();}

/** Number of lines in list */
inline Index nelem(const Array<Lines>& l) {Index n=0; for (auto& x:l) n+=nelem(x); return n;}

/** Number of lines in lists */
inline Index nelem(const Array<Array<Lines>>& l) {Index n=0; for (auto& x:l) n+=nelem(x); return n;}

/** Compute the reduced rovibrational dipole moment
 * 
 * @param[in] Jf Final J
 * @param[in] Ji Initial J
 * @param[in] lf Final l2
 * @param[in] li Initial l2
 * @param[in] k Type of transition
 * @return As titled
 */
Numeric reduced_rovibrational_dipole(Rational Jf, Rational Ji, Rational lf, Rational li, Rational k = 1);

/** Compute the reduced magnetic quadrapole moment
 * 
 * @param[in] Jf Final J
 * @param[in] Ji Initial J
 * @param[in] N The quantum number (upper should be equal to lower)
 * @return As titled
 */
Numeric reduced_magnetic_quadrapole(Rational Jf, Rational Ji, Rational N);
};  // Absorption

typedef Absorption::Lines AbsorptionLines;
typedef Array<AbsorptionLines> ArrayOfAbsorptionLines;
typedef Array<ArrayOfAbsorptionLines> ArrayOfArrayOfAbsorptionLines;

std::ostream& operator<<(std::ostream&, const ArrayOfAbsorptionLines&);

std::ostream& operator<<(std::ostream&, const ArrayOfArrayOfAbsorptionLines&);

#endif  // absorptionlines_h

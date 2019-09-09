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
 * @file   lineshapemodel.h
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

/** Describes the type of population level counter
 *
 * The types here might require that different data is available at runtime absorption calculations
 */
enum class PopulationType {
  ByLTE,                      // Assume line is in LTE
  ByVibrationalTemperatures,  // Assume line is in NLTE described by vibrational temperatures
  ByPopulationDistribution,   // Assume line is in NLTE and the upper-to-lower ratio is known
};  // LinePopulationType

/** Describes the type of cutoff calculations */
enum class CutoffType {
  None,                // No cutoff frequency at all
  LineByLineOffset,    // The cutoff frequency is at SingleLine::F0 plus the cutoff frequency
  BandFixedFrequency,  // The curoff frequency is the cutoff frequency for all SingleLine(s)
};  // LineCutoffType

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
  SingleLine(size_t nbroadeners, size_t nquanta) noexcept :
  mlineshape(nbroadeners), mlowerquanta(nquanta), mupperquanta(nquanta) {}
  
  /** Number of lineshape elements */
  Index LineShapeElems() const noexcept {return mlineshape.nelem();}
  
  /** Number of lower quantum numbers */
  Index LowerQuantumCount() const noexcept {return mlowerquanta.size();}
  
  /** Lower quantum number */
  Rational LowerQuantumNumber(size_t i) const noexcept {return mlowerquanta[i];}
  
  /** Number of upper quantum numbers */
  Index UpperQuantumCount() const noexcept {return mupperquanta.size();}
  
  /** Upper quantum number */
  Rational UpperQuantumNumber(size_t i) const noexcept {return mupperquanta[i];}
  
  /** Central frequency */
  Numeric F0() const noexcept {return mF0;}
  
  /** Lower level energy */
  Numeric E0() const noexcept {return mE0;}
  
  /** Reference line strength */
  LineShape::Type mlineshapetype;
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
};  // SingleLine

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
  std::vector<SpeciesTag> mbroadeningspecies;
  
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
        const std::vector<SpeciesTag>& broadeningspecies={},
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
        const std::vector<SpeciesTag>& broadeningspecies) :
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
               localquanta.size())) {};
  
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
    if(NumBroadeners() not_eq sl.LowerQuantumCount() or
       NumBroadeners() not_eq sl.UpperQuantumCount())
      throw std::runtime_error("Error calling appending function, bad size of quantum numbers");
    
    if(NumLines() not_eq 0 and 
       sl.LineShapeElems() not_eq mlines[0].LineShapeElems())
      throw std::runtime_error("Error calling appending function, bad size of broadening species");
    
    mlines.push_back(std::move(sl));
  }
  
  /** Species Index */
  Index Species() const noexcept {return mquantumidentity.Species();}
  
  /** Isotopologue Index */
  Index Isotopologue() const noexcept {return mquantumidentity.Isotopologue();}
  
  /** Number of lines */
  Index NumLines() const noexcept {return Index(mlines.size());}
  
  /** Number of broadening species */
  Index NumBroadeners() const noexcept {return Index(mlocalquanta.size());}
  
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
  
  /** Checks if all defined quantum numbers in qid are equal to the lower levels
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] qid Energy level quantum identifier
   * @return If qid is in lower level
   */
  bool InLowerLevel(size_t k, const QuantumIdentifier& qid) const noexcept;
  
  /** Checks if all defined quantum numbers in qid are equal to the upper levels 
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] qid Energy level quantum identifier
   * @return If qid is in upper level
   */
  bool InUpperLevel(size_t k, const QuantumIdentifier& qid) const noexcept;
  
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
  
  /** Central frequency
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Central frequency
   */
  Numeric F0(size_t k) const noexcept {return mlines[k].F0();}
  
  /** Lower level energy
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Lower level energy
   */
  Numeric E0(size_t k) const noexcept {return mlines[k].E0();}
  
  /** Reference line strength
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Reference line strength
   */
  Numeric I0(size_t k) const noexcept {return mlines[k].I0();}
  
  /** Einstein spontaneous emission
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Einstein spontaneous emission
   */
  Numeric A(size_t k) const noexcept {return mlines[k].A();}
  
  /** Lower level statistical weight
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Lower level statistical weight
   */
  Numeric g_low(size_t k) const noexcept {return mlines[k].g_low();}
  
  /** Upper level statistical weight
   * 
   * @param[in] k Line number (less than NumLines())
   * @return Upper level statistical weight
   */
  Numeric g_upp(size_t k) const noexcept {return mlines[k].g_upp();}
  
  /** Returns mirroring style */
  MirroringType Mirroring() const noexcept {return mmirroring;}
  
  /** Returns normalization style */
  NormalizationType Normalization() const noexcept {return mnormalization;}
  
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
  
  /** Returns cutoff frequency or 0
   * 
   * @param[in] k Line number (less than NumLines())
   * @returns Cutoff frequency or 0
   */
  Numeric Cutoff(size_t k) const noexcept {
    switch(mcutoff) {
      case CutoffType::LineByLineOffset:
        return F0(k) + mcutofffreq;
      case CutoffType::BandFixedFrequency:
        return mcutofffreq;
      case CutoffType::None:
        return 0;
    }
    std::terminate();
  }
};  // Lines

struct SingleLineInternalRead {SingleLine line; bool OK;};

SingleLineInternalRead ReadFromArtscat3Stream(istream& is);

SingleLineInternalRead ReadFromArtscat4Stream(istream& is, const std::vector<QuantumNumberType>& localquantums={});

SingleLineInternalRead ReadFromArtscat5Stream(istream& is, const std::vector<QuantumNumberType>& localquantums={});

SingleLineInternalRead ReadFromLBLRTMStream(istream& is, const std::vector<QuantumNumberType>& localquantums={});

SingleLineInternalRead ReadFromHitran2004Stream(istream& is, const std::vector<QuantumNumberType>& localquantums={}, Numeric fmin=-1);

};  // Absorption

typedef Absorption::Lines AbsorptionLines;
typedef Array<AbsorptionLines> ArrayOfAbsorptionLines;
typedef Array<ArrayOfAbsorptionLines> ArrayOfArrayOfAbsorptionLines;

#endif  // absorptionlines_h

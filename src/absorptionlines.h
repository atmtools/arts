/* Copyright (C) 2019
 Richard Larsson <ric.larsson@gmail.com>

 
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

#include "bifstream.h"
#include "bofstream.h"
#include "enums.h"
#include "lineshapemodel.h"
#include "matpack.h"
#include "quantum_numbers.h"
#include "species_tags.h"
#include "zeemandata.h"
#include "jacobian.h"

#include <utility>
#include <vector>


/** Namespace to contain things required for absorption calculations */
namespace Absorption {
/** Describes the type of mirroring line effects
 * 
 * Each type but None has to have an implemented effect
 */
ENUMCLASS(MirroringType, char,
  None,             // No mirroring
  Lorentz,          // Mirror, but use Lorentz line shape
  SameAsLineShape,  // Mirror using the same line shape
  Manual            // Mirror by having a line in the array of line record with negative F0
)  // MirroringType

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr std::string_view mirroringtype2metadatastring(MirroringType in) noexcept {
  switch (in) {
    case MirroringType::None:
      return "These lines are not mirrored at 0 Hz.\n";
    case MirroringType::Lorentz:
      return "These lines are mirrored around 0 Hz using the Lorentz line shape.\n";
    case MirroringType::SameAsLineShape:
      return "These line are mirrored around 0 Hz using the original line shape.\n";
    case MirroringType::Manual:
      return "There are manual line entries in the catalog to mirror this line.\n";
    case MirroringType::FINAL: break;
  }
}
#pragma GCC diagnostic pop

/** Describes the type of normalization line effects
 *
 * Each type but None has to have an implemented effect
 */
ENUMCLASS(NormalizationType, char,
  None,                // Do not renormalize the line shape
  VVH,                 // Renormalize with Van Vleck and Huber specifications
  VVW,                 // Renormalize with Van Vleck and Weiskopf specifications
  RQ,                  // Renormalize using Rosenkranz's quadratic specifications
  SFS                  // Renormalize using simple frequency scaling of the line strength
)  // NormalizationType

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr std::string_view normalizationtype2metadatastring(NormalizationType in) {
  switch (in) {
    case NormalizationType::None:
      return "No re-normalization in the far wing will be applied.\n";
    case NormalizationType::VVH:
      return "van Vleck and Huber far-wing renormalization will be applied, "
        "i.e. F ~ (f tanh(hf/2kT))/(f0 tanh(hf0/2kT))\n";
    case NormalizationType::VVW:
      return "van Vleck and Weisskopf far-wing renormalization will be applied, "
        "i.e. F ~ (f/f0)^2\n";
    case NormalizationType::RQ:
      return "Rosenkranz quadratic far-wing renormalization will be applied, "
        "i.e. F ~ hf0/2kT sinh(hf0/2kT) (f/f0)^2\n";
    case NormalizationType::SFS:
      return "Simple frequency scaling of the far-wings will be applied, "
        "i.e. F ~ (f / f0) * ((1 - exp(- hf / kT)) / (1 - exp(- hf0 / kT)))\n";
    case NormalizationType::FINAL: break;
  }
}
#pragma GCC diagnostic pop

/** Describes the type of population level counter
 *
 * The types here might require that different data is available at runtime absorption calculations
 */
ENUMCLASS(PopulationType, char,
  LTE,                            // Assume band is in LTE
  NLTE,                           // Assume band is in NLTE and the upper-to-lower ratio is known
  VibTemps,                       // Assume band is in NLTE described by vibrational temperatures and LTE at other levels
  ByHITRANRosenkranzRelmat,       // Assume band needs to compute relaxation matrix to derive HITRAN Y-coefficients
  ByHITRANFullRelmat,             // Assume band needs to compute and directly use the relaxation matrix according to HITRAN
  ByMakarovFullRelmat,            // Assume band needs to compute and directly use the relaxation matrix according to Makarov et al 2020
  ByRovibLinearDipoleLineMixing   // Assume band needs to compute and directly use the relaxation matrix according to Hartmann, Boulet, Robert, 2008, 1st edition
)  // PopulationType

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr std::string_view populationtype2metadatastring(PopulationType in) {
  switch (in) {
    case PopulationType::LTE:
      return "The lines are considered as in pure LTE.\n";
    case PopulationType::ByMakarovFullRelmat:
      return "The lines requires relaxation matrix calculations in LTE - Makarov et al 2020 full method.\n";
    case PopulationType::ByRovibLinearDipoleLineMixing:
      return "The lines requires relaxation matrix calculations in LTE - Hartmann, Boulet, Robert, 2008, 1st edition method.\n";
    case PopulationType::ByHITRANFullRelmat:
      return "The lines requires relaxation matrix calculations in LTE - HITRAN full method.\n";
    case PopulationType::ByHITRANRosenkranzRelmat:
      return "The lines requires Relaxation matrix calculations in LTE - HITRAN Rosenkranz method.\n";
    case PopulationType::VibTemps:
      return "The lines are considered as in NLTE by vibrational temperatures.\n";
    case PopulationType::NLTE:
      return "The lines are considered as in pure NLTE.\n";
    case PopulationType::FINAL: return "There's an error";
  }
}
#pragma GCC diagnostic pop

constexpr bool relaxationtype_relmat(PopulationType in) noexcept {
  return in == PopulationType::ByHITRANFullRelmat or
         in == PopulationType::ByMakarovFullRelmat or
         in == PopulationType::ByHITRANRosenkranzRelmat or
         in == PopulationType::ByRovibLinearDipoleLineMixing;
}

/** Describes the type of cutoff calculations */
ENUMCLASS(CutoffType, char,
  None,                             // No cutoff frequency at all
  ByLine                            // The cutoff frequency is at SingleLine::F0 plus the cutoff frequency plus the speed independent pressure shift
)  // CutoffType

String cutofftype2metadatastring(CutoffType in, Numeric cutoff);

/** Computations and data for a single absorption line */
struct SingleLine {
  /** Central frequency */
  Numeric F0{};
  
  /** Reference intensity */
  Numeric I0{};
  
  /** Lower state energy level */
  Numeric E0{};
  
  /** Lower level statistical weight */
  Numeric glow{};
  
  /** Upper level statistical weight */
  Numeric gupp{};
  
  /** Einstein spontaneous emission coefficient */
  Numeric A{};
  
  /** Zeeman model */
  Zeeman::Model zeeman{};
  
  /** Line shape model */
  LineShape::Model lineshape{};
  
  /** Local quantum numbers */
  Quantum::Number::LocalState localquanta{};

  /** Default initialization 
   * 
   * @param[in] F0_ Central frequency
   * @param[in] I0_ Reference line strength at external T0
   * @param[in] E0_ Lower energy level
   * @param[in] glow_ Lower level statistical weight
   * @param[in] gupp_ Upper level statistical weight
   * @param[in] A_ Einstein spontaneous emission coefficient
   * @param[in] zeeman_ Zeeman model
   * @param[in] lineshape_ Line shape model
   * @param[in] localquanta_ Local quantum numbers
   */
  SingleLine(Numeric F0_=0,
             Numeric I0_=0,
             Numeric E0_=0,
             Numeric glow_=0,
             Numeric gupp_=0,
             Numeric A_=0,
             Zeeman::Model zeeman_=Zeeman::Model(),
             LineShape::Model lineshape_=LineShape::Model(),
             Quantum::Number::LocalState localquanta_={}) :
             F0(F0_),
             I0(I0_),
             E0(E0_),
             glow(glow_),
             gupp(gupp_),
             A(A_),
             zeeman(zeeman_),
             lineshape(std::move(lineshape_)),
             localquanta(std::move(localquanta_)) {}
  
  /** Initialization for constant sizes
   * 
   * @param metaquanta A quantum number state with the right sizes and access points
   * @param metamodel A line shape model with the right sizes and access points
   */
  SingleLine(Quantum::Number::LocalState metaquanta, LineShape::Model metamodel) :
  lineshape(std::move(metamodel)), localquanta(std::move(metaquanta)) {}
  
  //////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////// Counts
  //////////////////////////////////////////////////////////////////
  
  /** Number of lineshape elements */
  [[nodiscard]] Index LineShapeElems() const noexcept {return lineshape.nelem();}
  
  /** Number of lower quantum numbers */
  [[nodiscard]] Index LocalQuantumElems() const noexcept {return localquanta.val.nelem();}
  
  //////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////// Special settings
  //////////////////////////////////////////////////////////////////
  
  /** Set Zeeman effect by automatic detection
   * 
   * Will fail if the available and provided quantum numbers are bad
   * 
   * @param[in] qid Copy of the global identifier to fill by local numbers
   */
  void SetAutomaticZeeman(QuantumIdentifier qid);
  
  /** Set the line mixing model to 2nd order
   * 
   * @param[in] d Data in 2nd order format
   */
  void SetLineMixing2SecondOrderData(const Vector& d);
  
  /** Set the line mixing model to AER kind
   * 
   * @param[in] d Data in AER format
   */
  void SetLineMixing2AER(const Vector& d);
  
  /** Binary read for AbsorptionLines */
  bifstream& read(bifstream& bif);
  
  /** Binary write for AbsorptionLines */
  bofstream& write(bofstream& bof) const;
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
  PopulationType population=PopulationType::LTE;
  NormalizationType normalization=NormalizationType::None;
  LineShape::Type lineshapetype=LineShape::Type::DP;
  Numeric T0=0;
  Numeric cutofffreq=0;
  Numeric linemixinglimit=-1;
  QuantumIdentifier quantumidentity;
  ArrayOfSpecies species;
  SingleLine line;
};

struct Lines {
  static constexpr Index version = 2;

  /** Does the line broadening have self broadening */
  bool selfbroadening;
  
  /** Does the line broadening have bath broadening */
  bool bathbroadening;
  
  /** cutoff type, by band or by line */
  CutoffType cutoff;
  
  /** Mirroring type */
  MirroringType mirroring;
  
  /** Line population distribution */
  PopulationType population;
  
  /** Line normalization type */
  NormalizationType normalization;

  /** Type of line shape */
  LineShape::Type lineshapetype;
  
  /** Reference temperature for all parameters of the lines */
  Numeric T0;
  
  /** cutoff frequency */
  Numeric cutofffreq;
  
  /** linemixing limit */
  Numeric linemixinglimit;
  
  /** Catalog ID */
  QuantumIdentifier quantumidentity;
  
  /** A list of broadening species */
  ArrayOfSpecies broadeningspecies;
  
  /** A list of individual lines */
  Array<SingleLine> lines;
  
  /** Default initialization
   * 
   * @param[in] selfbroadening_ Do self broadening
   * @param[in] bathbroadening_ Do bath broadening
   * @param[in] cutoff_ Type of cutoff frequency
   * @param[in] mirroring_ Type of mirroring
   * @param[in] population_ Type of line strengths distributions
   * @param[in] normalization_ Type of normalization
   * @param[in] lineshapetype_ Type of line shape
   * @param[in] T0_ Reference temperature
   * @param[in] cutofffreq_ Cutoff frequency
   * @param[in] linemixinglimit_ Line mixing limit
   * @param[in] quantumidentity_ Identity of global lines
   * @param[in] broadeningspecies_ List of broadening species
   * @param[in] lines_ List of SingleLine(s)
   */
  Lines(bool selfbroadening_=false,
        bool bathbroadening_=false,
        CutoffType cutoff_=CutoffType::None,
        MirroringType mirroring_=MirroringType::None,
        PopulationType population_=PopulationType::LTE,
        NormalizationType normalization_=NormalizationType::None,
        LineShape::Type lineshapetype_=LineShape::Type::DP,
        Numeric T0_=296,
        Numeric cutofffreq_=-1,
        Numeric linemixinglimit_=-1,
        QuantumIdentifier quantumidentity_=QuantumIdentifier(),
        ArrayOfSpecies broadeningspecies_={},
        Array<SingleLine> lines_={}) :
        selfbroadening(selfbroadening_),
        bathbroadening(bathbroadening_),
        cutoff(cutoff_),
        mirroring(mirroring_),
        population(population_),
        normalization(normalization_),
        lineshapetype(lineshapetype_),
        T0(T0_),
        cutofffreq(cutofffreq_),
        linemixinglimit(linemixinglimit_),
        quantumidentity(std::move(quantumidentity_)),
        broadeningspecies(std::move(broadeningspecies_)),
        lines(std::move(lines_)) {
    if (selfbroadening) broadeningspecies.front() = quantumidentity.Species();
    if (bathbroadening) broadeningspecies.back() = Species::Species::Bath;
  }
  
  /** XML-tag initialization
   * 
   * @param[in] selfbroadening_ Do self broadening
   * @param[in] bathbroadening_ Do bath broadening
   * @param[in] nlines Number of SingleLine(s) to initiate as empty
   * @param[in] cutoff_ Type of cutoff frequency
   * @param[in] mirroring_ Type of mirroring
   * @param[in] population_ Type of line strengths distributions
   * @param[in] normalization_ Type of normalization
   * @param[in] lineshapetype_ Type of line shape
   * @param[in] T0_ Reference temperature
   * @param[in] cutofffreq_ Cutoff frequency
   * @param[in] linemixinglimit_ Line mixing limit
   * @param[in] quantumidentity_ Identity of global lines
   * @param[in] broadeningspecies_ List of broadening species
   * @param[in] metalocalquanta A local state with defined quantum numbers
   * @param[in] metamodel A line shape model with defined shapes
   */
  Lines(bool selfbroadening_,
        bool bathbroadening_,
        size_t nlines,
        CutoffType cutoff_,
        MirroringType mirroring_,
        PopulationType population_,
        NormalizationType normalization_,
        LineShape::Type lineshapetype_,
        Numeric T0_,
        Numeric cutofffreq_,
        Numeric linemixinglimit_,
        QuantumIdentifier  quantumidentity_,
        ArrayOfSpecies  broadeningspecies_,
        const Quantum::Number::LocalState& metalocalquanta,
        const LineShape::Model& metamodel) :
        selfbroadening(selfbroadening_),
        bathbroadening(bathbroadening_),
        cutoff(cutoff_),
        mirroring(mirroring_),
        population(population_),
        normalization(normalization_),
        lineshapetype(lineshapetype_),
        T0(T0_),
        cutofffreq(cutofffreq_),
        linemixinglimit(linemixinglimit_),
        quantumidentity(std::move(quantumidentity_)),
        broadeningspecies(std::move(broadeningspecies_)),
        lines(nlines, SingleLine(metalocalquanta, metamodel)) {
    if (selfbroadening) broadeningspecies.front() = quantumidentity.Species();
    if (bathbroadening) broadeningspecies.back() = Species::Species::Bath;
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
  void AppendSingleLine(SingleLine&& sl);
  
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
  void AppendSingleLine(const SingleLine& sl);
  
  /** Checks if an external line matches this structure
   * 
   * @param[in] sle Full external lines
   * @param[in] quantumidentity Expected global quantum id of the line
   */
  [[nodiscard]] bool MatchWithExternal(const SingleLineExternal& sle, const QuantumIdentifier& quantumidentity) const ARTS_NOEXCEPT;
  
  /** Checks if another line list matches this structure
   * 
   * @param[in] sle Full external lines
   * @param[in] quantumidentity Expected global quantum id of the line
   * @return first: match; second: nullable line shape
   */
  [[nodiscard]] std::pair<bool, bool> Match(const Lines& l) const noexcept;
  
  /** Sort inner line list by frequency */
  void sort_by_frequency();
  
  /** Sort inner line list by Einstein coefficient */
  void sort_by_einstein();
  
  /** Species Name */
  [[nodiscard]] String SpeciesName() const noexcept;
  
  /** Meta data for the line shape if it exists */
  [[nodiscard]] String LineShapeMetaData() const noexcept;
  
  /** Species Enum */
  [[nodiscard]] Species::Species Species() const noexcept {return quantumidentity.Species();}
  
  /** Isotopologue Index */
  [[nodiscard]] Species::IsotopeRecord Isotopologue() const noexcept {return quantumidentity.Isotopologue();}
  
  /** Number of lines */
  [[nodiscard]] Index NumLines() const noexcept {return Index(lines.size());}
  
  /** Make a common line shape if possible */
  void MakeLineShapeModelCommon();
  
  /** Number of broadening species */
  [[nodiscard]] Index NumBroadeners() const noexcept {return Index(broadeningspecies.nelem());}
  
  /** Number of broadening species */
  [[nodiscard]] Index NumLocalQuanta() const noexcept {
    return lines.size() ? lines.front().localquanta.val.nelem() : 0;}
  
  /** Returns the number of Zeeman split lines
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] type Type of Zeeman polarization
   */
  [[nodiscard]] Index ZeemanCount(size_t k, Zeeman::Polarization type) const noexcept;
  
  /** Returns the strength of a Zeeman split line
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] type Type of Zeeman polarization
   * @param[in] i Zeeman line count
   */
  [[nodiscard]] Numeric ZeemanStrength(size_t k, Zeeman::Polarization type, Index i) const noexcept;
  
  /** Returns the splitting of a Zeeman split line
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] type Type of Zeeman polarization
   * @param[in] i Zeeman line count
   */
  [[nodiscard]] Numeric ZeemanSplitting(size_t k, Zeeman::Polarization type, Index i) const noexcept;
  
  /** Set Zeeman effect for all lines that have the correct quantum numbers */
  void SetAutomaticZeeman() noexcept;
  
  /** Mean frequency by weight of line strength
   * 
   * @param[in] T Temperature at which to compute the line strength (T <= 0 means at T0 is used)
   * @return Mean frequency
   */
  [[nodiscard]] Numeric F_mean(Numeric T=0) const noexcept;
  
  /** Mean frequency by weight of line strengt
   * 
   * @param[in] wgts Weight of averaging
   * @return Mean frequency
   */
  [[nodiscard]] Numeric F_mean(const ConstVectorView& wgts) const noexcept;
  
  /** On-the-fly line mixing */
  [[nodiscard]] bool OnTheFlyLineMixing() const noexcept {
    return population == PopulationType::ByMakarovFullRelmat or
           population == PopulationType::ByRovibLinearDipoleLineMixing;
  }
  
  /** Returns if the pressure should do line mixing
   * 
   * @param[in] P Atmospheric pressure
   * @return true if no limit or P less than limit
   */
  [[nodiscard]] bool DoLineMixing(Numeric P) const noexcept {
    return linemixinglimit < 0 ? true : linemixinglimit > P;
  }

  [[nodiscard]] bool DoVmrDerivative(const QuantumIdentifier& qid) const noexcept {
    return qid.Isotopologue() == quantumidentity.Isotopologue() or 
    (qid.Isotopologue().joker() and qid.Species() == quantumidentity.Species()) or
    std::any_of(broadeningspecies.begin(), broadeningspecies.end(), [s=qid.Species()](auto& a){return a == s;});
  }

  /** @return Whether the band may require linemixing */
  [[nodiscard]] bool AnyLinemixing() const noexcept {
    for (auto& line : lines) {
      for (auto& shape : line.lineshape.Data()) {
        if (shape.Y().type not_eq LineShape::TemperatureModel::None or
            shape.G().type not_eq LineShape::TemperatureModel::None or
            shape.DV().type not_eq LineShape::TemperatureModel::None) {
          return true;
        }
      }
    }
    return false;
  }

  /** Line shape parameters
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] vmrs Line broadener species's volume mixing ratio
   * @return Line shape parameters
   */
  [[nodiscard]] LineShape::Output ShapeParameters(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept;
  
  /** Line shape parameters
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] m Line broadening species position
   * @return Line shape parameters
   */
  [[nodiscard]] LineShape::Output ShapeParameters(size_t k, Numeric T, Numeric P, size_t m) const noexcept;
  
  /** Line shape parameters temperature derivatives
   * 
   * @param[in] k Line number (less than NumLines())
   * @param[in] T Atmospheric temperature
   * @param[in] P Atmospheric pressure
   * @param[in] vmrs Line broadener's volume mixing ratio
   * @return Line shape parameters temperature derivatives
   */
  [[nodiscard]] LineShape::Output ShapeParameters_dT(size_t k, Numeric T, Numeric P, const Vector& vmrs) const noexcept;
  
  /** Position among broadening species or -1
   * 
   * @param[in] A species index that might be among the broadener species
   * @return Position among broadening species or -1
   */
  [[nodiscard]] Index LineShapePos(const Species::Species spec) const noexcept;
  
  /** Position among broadening species or -1
   * 
   * @param[in] An identity that might be among the broadener species
   * @return Position among broadening species or -1
   */
  [[nodiscard]] Index LineShapePos(const QuantumIdentifier& qid) const noexcept {
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
  [[nodiscard]] LineShape::Output ShapeParameters_dVMR(size_t k, Numeric T, Numeric P,
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
  [[nodiscard]] Numeric ShapeParameter_dInternal(size_t k, Numeric T, Numeric P,
                                   const Vector& vmrs,
                                   const RetrievalQuantity& derivative) const noexcept;
  
  /** Returns cutoff frequency or maximum value
   * 
   * @param[in] k Line number (less than NumLines())
   * @returns Cutoff frequency or 0
   */
  [[nodiscard]] Numeric CutoffFreq(size_t k, Numeric shift=0) const noexcept;
  
  /** Returns negative cutoff frequency or lowest value
   * 
   * @param[in] k Line number (less than NumLines())
   * @returns Negative cutoff frequency or the lowest value
   */
  [[nodiscard]] Numeric CutoffFreqMinus(size_t k, Numeric shift=0) const noexcept;
  
  /** Position of species if available or -1 else */
  [[nodiscard]] Index BroadeningSpeciesPosition(Species::Species spec) const noexcept {
    if (auto ptr = std::find(broadeningspecies.cbegin(),
      broadeningspecies.cend(), spec); ptr not_eq broadeningspecies.cend())
      return std::distance(broadeningspecies.cbegin(), ptr);
    return -1;
  }
  
  /** Returns a printable statement about the lines */
  [[nodiscard]] String MetaData() const;
  
  /** Removes a single line */
  void RemoveLine(Index) noexcept;
  
  /** Pops a single line */
  SingleLine PopLine(Index) noexcept;
  
  /** Reverses the order of the internal lines */
  void ReverseLines() noexcept;
  
  /** Mass of the molecule */
  [[nodiscard]] Numeric SpeciesMass() const noexcept;
  
  /** Returns the VMRs of the broadening species
   * 
   * @param[in] atm_vmrs Atmospheric VMRs
   * @param[in] atm_spec Atmospheric Species
   * @return VMR list of the species
   */
  [[nodiscard]] Vector BroadeningSpeciesVMR(const ConstVectorView&, const ArrayOfArrayOfSpeciesTag&) const;
  
  /** Returns the mass of the broadening species
   * 
   * @param[in] atm_vmrs Atmospheric VMRs
   * @param[in] atm_spec Atmospheric Species
   * @param[in] bath_mass Mass of Bath/Air (optional, will compute it if <=0)
   * @return Mass list of the species
   */
  [[nodiscard]] Vector BroadeningSpeciesMass(const ConstVectorView&, const ArrayOfArrayOfSpeciesTag&, const SpeciesIsotopologueRatios&, const Numeric& bath_mass=0) const;
  
  /** Returns the VMR of the species
   * 
   * @param[in] atm_vmrs Atmospheric VMRs
   * @param[in] atm_spec Atmospheric Species
   * @return VMR of the species
   */
  [[nodiscard]] Numeric SelfVMR(const ConstVectorView&, const ArrayOfArrayOfSpeciesTag&) const;
  
  /** Binary read for Lines */
  bifstream& read(bifstream& is);
  
  /** Binary write for Lines */
  bofstream& write(bofstream& os) const;
  
  [[nodiscard]] bool OK() const noexcept;
  
  [[nodiscard]] Numeric DopplerConstant(Numeric T) const noexcept;

  [[nodiscard]] QuantumIdentifier QuantumIdentityOfLine(Index k) const noexcept;
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
 * Beyond the HITRAN pre-2004 format, there is one more tag for line mixing
 * available in LBLRTM.  This is a sign at the end of the line to indicate that
 * the very next line gives line mixing information.
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
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromHitran2004Stream(istream& is);

/** Read from HITRAN online
 * 
 * The data format from online should be a .par line
 * followed by upper state quantum numbers and then
 * lower state quantum numbers.  See ReadFromHitran2004Stream
 * for the format of the .par-bit.  The quantum numbers are
 * parsed by name and should look as:
 * 
 * J=5.5;N1=2.5;parity=-;kronigParity=f [[tab]] J=6.5;N1=2.5;parity=-;kronigParity=f
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
*/ 
SingleLineExternal ReadFromHitranOnlineStream(istream& is);

/** Read from HITRAN before 2004
 * 
 * See ReadFromLBLRTMStream for details on format
 * 
 * @param[in] is Input stream
 * @return SingleLineExternal 
 */
SingleLineExternal ReadFromHitran2001Stream(istream& is);

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
std::vector<Lines> split_list_of_external_lines(std::vector<SingleLineExternal>& external_lines,
                                                const std::vector<QuantumNumberType>& localquantas={},
                                                const std::vector<QuantumNumberType>& globalquantas={});

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
Numeric reduced_rovibrational_dipole(Rational Jf, Rational Ji, Rational lf, Rational li, Rational k = Rational(1));

/** Compute the reduced magnetic quadrapole moment
 * 
 * @param[in] Jf Final J
 * @param[in] Ji Initial J
 * @param[in] N The quantum number (upper should be equal to lower)
 * @return As titled
 */
Numeric reduced_magnetic_quadrapole(Rational Jf, Rational Ji, Rational N);
} // namespace Absorption

using AbsorptionSingleLine = Absorption::SingleLine;
using ArrayOfAbsorptionSingleLine = Array<AbsorptionSingleLine>;
using AbsorptionLines = Absorption::Lines;
using ArrayOfAbsorptionLines = Array<AbsorptionLines>;
using ArrayOfArrayOfAbsorptionLines = Array<ArrayOfAbsorptionLines>;

std::ostream& operator<<(std::ostream&, const ArrayOfAbsorptionLines&);

std::ostream& operator<<(std::ostream&, const ArrayOfArrayOfAbsorptionLines&);

using AbsorptionNormalizationType = Absorption::NormalizationType;
using AbsorptionPopulationType = Absorption::PopulationType;
using AbsorptionMirroringType = Absorption::MirroringType;
using AbsorptionCutoffType = Absorption::CutoffType;

struct AbsorptionMirroringTagTypeStatus {
  bool None{false}, Lorentz{false}, SameAsLineShape{false}, Manual{false};
  AbsorptionMirroringTagTypeStatus(const ArrayOfArrayOfAbsorptionLines&);
  friend std::ostream& operator<<(std::ostream&, AbsorptionMirroringTagTypeStatus);
};

struct AbsorptionNormalizationTagTypeStatus {
  bool None{false}, VVH{false}, VVW{false}, RQ{false}, SFS{false};
  AbsorptionNormalizationTagTypeStatus(const ArrayOfArrayOfAbsorptionLines&);
  friend std::ostream& operator<<(std::ostream&, AbsorptionNormalizationTagTypeStatus);
};

struct AbsorptionPopulationTagTypeStatus {
  bool LTE{false}, NLTE{false}, VibTemps{false},
      ByHITRANRosenkranzRelmat{false}, ByHITRANFullRelmat{false},
      ByMakarovFullRelmat{false}, ByRovibLinearDipoleLineMixing{false};
  AbsorptionPopulationTagTypeStatus(const ArrayOfArrayOfAbsorptionLines&);
  friend std::ostream& operator<<(std::ostream&, AbsorptionPopulationTagTypeStatus);
};

struct AbsorptionCutoffTagTypeStatus {
  bool None{false}, ByLine{false};
  AbsorptionCutoffTagTypeStatus(const ArrayOfArrayOfAbsorptionLines&);
  friend std::ostream& operator<<(std::ostream&, AbsorptionCutoffTagTypeStatus);
};

struct AbsorptionLineShapeTagTypeStatus {
  bool DP{false}, LP{false}, VP{false}, SDVP{false}, HTP{false};
  AbsorptionLineShapeTagTypeStatus(const ArrayOfArrayOfAbsorptionLines&);
  friend std::ostream& operator<<(std::ostream&, AbsorptionLineShapeTagTypeStatus);
};

struct AbsorptionTagTypesStatus {
  AbsorptionMirroringTagTypeStatus mirroring;
  AbsorptionNormalizationTagTypeStatus normalization;
  AbsorptionPopulationTagTypeStatus population;
  AbsorptionCutoffTagTypeStatus cutoff;
  AbsorptionLineShapeTagTypeStatus lineshapetype;

  AbsorptionTagTypesStatus(const ArrayOfArrayOfAbsorptionLines& lines)
      : mirroring(lines),
        normalization(lines),
        population(lines),
        cutoff(lines),
        lineshapetype(lines) {}
  friend std::ostream& operator<<(std::ostream&, AbsorptionTagTypesStatus);
};

//! Helper struct for flat_index
struct AbsorptionSpeciesBandIndex {
  //! The species index in abs_species/abs_lines_per_species
  Index ispecies;

  //! The band index in abs_lines_per_species[ispecies]
  Index iband;
};

/** Get a flat index pair for species and band

  @param[in] i: Index smaller than the total number of bands but at least 0
  @param[in] abs_species: As WSV
  @param[in] abs_lines_per_species: As WSV
  @return A valid AbsorptionSpeciesBandIndex
*/
AbsorptionSpeciesBandIndex flat_index(
    Index i,
    const ArrayOfArrayOfSpeciesTag& abs_species,
    const ArrayOfArrayOfAbsorptionLines& abs_lines_per_species);

#endif  // absorptionlines_h

/* Copyright (C) 2018
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

/** Contains the line function data class
 * \file   linefunctiondata.h
 * 
 * \author Richard Larsson
 * \date   2018-09-19
 **/

#ifndef linefunctiondata_h
#define linefunctiondata_h

// For LEGACY only
#include "pressurebroadeningdata.h"
#include "linemixingdata.h"

// Actually needed
#include "abs_species_tags.h"
#include "jacobian.h"


// List of all variables that can be returned
typedef struct{
  /** Speed independent broadening */
  Numeric G0;
  
  /** Speed independent shifting */
  Numeric D0;
  
  /** Speed dependent broadening */
  Numeric G2;
  
  /** Speed dependent shifting */
  Numeric D2;
  
  /** Frequency of velocity changing collisions */
  Numeric FVC;
  
  /** The correlation parameter */
  Numeric ETA;
  
  /** First order line mixing parameter */
  Numeric Y;
  
  /** Second order line mixing parameter --- changes the strength of the line */
  Numeric G;
  
  /** Second order line mixing parameter --- changes the frequency of the line */
  Numeric DV;
} LineFunctionDataOutput;
std::ostream& operator<<(std::ostream& os, const LineFunctionDataOutput& v);

class LineFunctionData {
public:

  #define LineFunctionData_SelfBroadening String("SELF")
  #define LineFunctionData_BathBroadening String("AIR")
  #define LineFunctionData_NoLineMixing String("#")
  
  /** Line shape models */
  enum class LineShapeType {DP, LP, VP, SDVP, HTP, };
  
  /** Line mixing models */
  enum class LineMixingOrderType {None, LM1, LM2, ConstG, Interp, };
  
  /** Temperature models; T# for analytical, rest for interpolated */
  enum class TemperatureType {None, T0, T1, T2, T3, T4, T5, LM_AER };
  
  /** Doppler line shape model parameters */
  enum class DopplerParam    : Index {                          Size};
  
  /** Lorentz line shape model parameters */
  enum class LorentzParam    : Index {G0, D0,                   Size};
  
  /** Voigt line shape model parameters */
  enum class VoigtParam      : Index {G0, D0,                   Size};
  
  /** Speed-dependent Voigt line shape model parameters */
  enum class SpeedVoigtParam : Index {G0, D0, G2, D2,           Size};
  
  /** Hartmann-Tran line shape model parameters */
  enum class HTPParam        : Index {G0, D0, G2, D2, FVC, ETA, Size};
  
  /** No line mixing model */
  enum class NoLineMixingParam : Index {                        Size};
  
  /** First order line mixing model */
  enum class FirstOrderParam   : Index {Y,                      Size};
  
  /** Second-order line mixing model */
  enum class SecondOrderParam  : Index {Y, G, DV,               Size};
  
  /** Only compute G line mixing model */
  enum class ConstGParam       : Index {   G,                   Size};
  
  /** Interpolated line mixing model */
  enum class InterpParam       : Index {INTERPOLATED_VARIABLES, Size};
  
  /** Default constructor to be used with iostream */
  LineFunctionData() = default;
  
  /** Air-broadening constructor to keep up with old definitions of pure air-broadening
   * 
   * Will set up so the calculations are:
   * G0 = sgam * svmr * theta^nself + agam * (1 - avmr) * theta&nair
   * D0 = psf * theta^(1.5*nair + 0.25)
   * 
   * The rest are zero
   * 
   * \param sgam Self speed independent pressure broadening coefficient
   * \param nself Self speed independent pressure broadening exponent coefficient
   * \param agam Air speed independent pressure broadening coefficient
   * \param nair Air speed independent pressure broadening exponent exponent coefficient
   * \param psf Air speed independent pressure frequency shift coefficient
   */
  LineFunctionData(const Numeric& sgam, const Numeric& nself, const Numeric& agam, const Numeric& nair, const Numeric& psf) :
  do_line_in_standard_calculations(true), mself(true), mbath(true), mspecies(2), mtypes(2, {TemperatureType::T1, TemperatureType::T5}),
  mdata(2), merrors(2, Vector(4, 0)), mp(LineShapeType::VP), mlm(LineMixingOrderType::None) {
    mdata[0] = {sgam, nself, psf, nair};
    mdata[1] = {agam, nair, psf, nair};
  };
  
  /** General-broadening constructor to keep up with old definitions of all types of broadening.
   
   This constructor should be considered deprecated but still exists for historical reasons. 
   
   \param pb Pressure broadening data in old format
   \param lm Line mixing data in old format
   \param species Self species
   \param T0 Reference temperature of parameters in lm and pb
   */
  LineFunctionData(const PressureBroadeningData& pb, const LineMixingData& lm, const String& species, const Numeric& T0);
  
  /** Returns the number of line mixing elements requested */
  Index LineMixingTypeNelem() const noexcept {
    switch(mlm) {
      case LineMixingOrderType::None:   return Index(NoLineMixingParam::Size);
      case LineMixingOrderType::LM1:    return Index(FirstOrderParam::Size);
      case LineMixingOrderType::LM2:    return Index(SecondOrderParam::Size);
      case LineMixingOrderType::Interp: return Index(InterpParam::Size);
      case LineMixingOrderType::ConstG: return Index(ConstGParam::Size);;
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  /** Returns the storage string for this line mixing */
  String LineMixingType2String() const noexcept {
    switch(mlm) {
      case LineMixingOrderType::None:   return LineFunctionData_NoLineMixing;
      case LineMixingOrderType::LM1:    return "LM1";
      case LineMixingOrderType::LM2:    return "LM2";
      case LineMixingOrderType::Interp: return "INT";
      case LineMixingOrderType::ConstG: return "ConstG";
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  /** Sets line mixing type according to string describing type */
  void StringSetLineMixingType(const String& type) {
    if(type == LineFunctionData_NoLineMixing) mlm = LineMixingOrderType::None; 
    else if(type == String("LM1"))            mlm = LineMixingOrderType::LM1;
    else if(type == String("LM2"))            mlm = LineMixingOrderType::LM2;
    else if(type == String("INT"))            mlm = LineMixingOrderType::Interp;
    else if(type == String("ConstG"))         mlm = LineMixingOrderType::ConstG;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }
  
  /** Returns the number of line shape elements requested */
  Index LineShapeTypeNelem() const noexcept {
    switch(mp) {
      case LineShapeType::DP:   return Index(DopplerParam::Size);
      case LineShapeType::LP:   return Index(LorentzParam::Size);
      case LineShapeType::VP:   return Index(VoigtParam::Size);
      case LineShapeType::SDVP: return Index(SpeedVoigtParam::Size);
      case LineShapeType::HTP:  return Index(HTPParam::Size);
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  /** Returns the storage string for this line shape */
  String LineShapeType2String() const noexcept {
    switch(mp) {
      case LineShapeType::DP:   return "DP";
      case LineShapeType::LP:   return "LP";
      case LineShapeType::VP:   return "VP";
      case LineShapeType::SDVP: return "SDVP";
      case LineShapeType::HTP:  return "HTP";
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  /** Sets line shape type according to string describing type */
  void StringSetLineShapeType(const String& type) {
    if(type == String("DP"))        mp = LineShapeType::DP;
    else if(type == String("LP"))   mp = LineShapeType::LP;
    else if(type == String("VP"))   mp = LineShapeType::VP;
    else if(type == String("SDVP")) mp = LineShapeType::SDVP;
    else if(type == String("HTP"))  mp = LineShapeType::HTP;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }
  
  /** Returns the number of parameters for the select temperature model */
  Index TemperatureTypeNelem(const TemperatureType type) const noexcept {
    switch(type) {
      case TemperatureType::LM_AER: return 12;
      case TemperatureType::None:   return 0;
      case TemperatureType::T0:     return 1;
      case TemperatureType::T1:     return 2;
      case TemperatureType::T2:     return 3;
      case TemperatureType::T3:     return 2;
      case TemperatureType::T4:     return 3;
      case TemperatureType::T5:     return 2;
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  /** Returns the number of temperature parameters requested for the line shape model */
  Index LineShapeDataNelemForSpecies(Index ispec) const noexcept {
    Index count=0;
    for(auto i=0; i<LineShapeTypeNelem(); i++)
      count += TemperatureTypeNelem(mtypes[ispec][i]);
    return count;
  }
  
  /** Returns the storage string for the temperature model type */
  String TemperatureType2String(const TemperatureType type) const noexcept {
    switch(type) {
      case TemperatureType::LM_AER: return "LM_AER";
      case TemperatureType::None:   return LineFunctionData_NoLineMixing;
      case TemperatureType::T0:     return "T0";
      case TemperatureType::T1:     return "T1";
      case TemperatureType::T2:     return "T2";
      case TemperatureType::T3:     return "T3";
      case TemperatureType::T4:     return "T4";
      case TemperatureType::T5:     return "T5";
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  /** Sets the temperature model type from the storage string */
  void StringSetTemperatureType(const Index ispecies, const Index iparam, const String& type) {
    if(type == String("LM_AER"))                   mtypes[ispecies][iparam] = TemperatureType::LM_AER;
    else if(type == LineFunctionData_NoLineMixing) mtypes[ispecies][iparam] = TemperatureType::None;
    else if(type == String("T0"))                  mtypes[ispecies][iparam] = TemperatureType::T0;
    else if(type == String("T1"))                  mtypes[ispecies][iparam] = TemperatureType::T1;
    else if(type == String("T2"))                  mtypes[ispecies][iparam] = TemperatureType::T2;
    else if(type == String("T3"))                  mtypes[ispecies][iparam] = TemperatureType::T3;
    else if(type == String("T4"))                  mtypes[ispecies][iparam] = TemperatureType::T4;
    else if(type == String("T5"))                  mtypes[ispecies][iparam] = TemperatureType::T5;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }

  LineFunctionDataOutput GetParams(const Numeric& T0,
                                   const Numeric& T,
                                   const Numeric& P,
                                   const Numeric& self_vmr,
                                   const ConstVectorView& rtp_vmr,
                                   const ArrayOfArrayOfSpeciesTag& abs_species,
                                   const bool do_linemixing=true,
                                   const bool normalization=true) const noexcept;

  LineFunctionDataOutput GetTemperatureDerivs(const Numeric& T0,
                                              const Numeric& T,
                                              const Numeric& dT,
                                              const Numeric& P,
                                              const Numeric& self_vmr,
                                              const ConstVectorView& rtp_vmr,
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const bool do_linemixing=true,
                                              const bool normalization=true) const noexcept;

  LineFunctionDataOutput GetReferenceT0Derivs(const Numeric& T0,
                                              const Numeric& T,
                                              const Numeric& P,
                                              const Numeric& self_vmr,
                                              const ConstVectorView& rtp_vmr, 
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const RetrievalQuantity& rt, 
                                              const QuantumIdentifier& line_qi,
                                              const bool do_linemixing=true,
                                              const bool normalization=true) const noexcept;

  LineFunctionDataOutput GetVMRDerivs(const Numeric& T0,
                                      const Numeric& T,
                                      const Numeric& P,
                                      const Numeric& self_vmr,
                                      const ConstVectorView& rtp_vmr,
                                      const ArrayOfArrayOfSpeciesTag& abs_species,
                                      const QuantumIdentifier& vmr_qi, 
                                      const QuantumIdentifier& line_qi,
                                      const bool do_linemixing=true,
                                      const bool normalization=true) const noexcept;
              
  Numeric GetLineParamDeriv(const Numeric& T0,
                            const Numeric& T,
                            const Numeric& P,
                            const Numeric& self_vmr,
                            const ConstVectorView& rtp_vmr, 
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const RetrievalQuantity& rt, 
                            const QuantumIdentifier& line_qi,
                            const bool do_linemixing=true,
                            const bool normalization=true) const noexcept;
                            
  //! Access to mspecies
  const ArrayOfSpeciesTag& Species() const noexcept { return mspecies; }
  
  //! Access to mself
  bool SelfBroadening() const noexcept { return mself; }
  
  //! Access to mbath
  bool AirBroadening() const noexcept { return mbath; }
  
  //! Access to do_line_in_standard_calculations
  bool DoStandardCalcs() const noexcept { return do_line_in_standard_calculations; }
  
  //! Access to mp
  LineShapeType LineShape() const noexcept { return mp; }
  
  //! Access to mlm
  LineMixingOrderType LineMixing() const noexcept { return mlm; }
  
  friend std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd);
  friend std::istream& operator>>(std::istream& is, LineFunctionData& lfd);
  
  // Special functions translating old behavior, all DEPRECATED
  //! Access to self broadening if AirBroadening was used as constructor else UB
  Numeric SelfG0() const;
  
  //! Access to self broadening exponent if AirBroadening was used as constructor else UB
  Numeric SelfN() const;
  
  //! Access to air broadening if AirBroadening was used as constructor else UB
  Numeric AirG0() const;
  
  //! Access to air freq shift if AirBroadening was used as constructor else UB
  Numeric AirD0() const;
  
  //! Access to air broadening exponent if AirBroadening was used as constructor else UB
  Numeric AirN() const;
  
  //! Access to self broadening error if AirBroadening was used as constructor else UB
  Numeric dSelfG0() const;
  
  //! Access to self broadening exponent error if AirBroadening was used as constructor else UB
  Numeric dSelfN() const;
  
  //! Access to air broadening error if AirBroadening was used as constructor else UB
  Numeric dAirG0() const;
  
  //! Access to air freq shifting error if AirBroadening was used as constructor else UB
  Numeric dAirD0() const;
  
  //! Access to air broadening exponent error if AirBroadening was used as constructor else UB
  Numeric dAirN() const;
  
  //! Access to Planetary broadening data if following old method else UB
  Vector PlanetaryForeignG0() const;
  
  //! Access to Planetary freq shift data if following old method else UB
  Vector PlanetaryForeignD0() const;
  
  //! Access to Planetary broadening exponent data if following old method else UB
  Vector PlanetaryForeignN() const;
  
  LineFunctionDataOutput AirBroadening(const Numeric& theta, const Numeric& P, const Numeric& self_vmr) const;
  
  //! Change line mixing model by replacing the current one with the new data and new temperature dependencies  
  void ChangeLineMixing(const LineMixingOrderType lm, const Array<TemperatureType>& ts, const Vector& data);
  
  //! Change line mixing model to LM2 style from Tretyakov et al.
  void ChangeLineMixingfromSimpleLM2(const Vector& data) {ChangeLineMixing(LineMixingOrderType::LM2, {TemperatureType::T4, TemperatureType::T4, TemperatureType::T4}, data);}
  
  //! Change line mixing model to AER style
  void ChangeLineMixing2AER(const Vector& data) {ChangeLineMixing(LineMixingOrderType::Interp, {TemperatureType::LM_AER}, data);}
  
  //! Set do_line_in_standard_calculations
  void SetStandard(const bool standard) {do_line_in_standard_calculations=standard;}
  
  //! Set variable and coefficient of species to X.  Slowly.
  void Set(const Numeric& X, const String& species, const String& coefficient, const String& variable);
  
  //! Get variable and coefficient of species.  Slowly.
  Numeric Get(const String& species, const String& coefficient, const String& variable) const;
  
  //! Remove a species
  void Remove(Index);
  
  //! Remove self if it exists
  void RemoveSelf() {if(mself) Remove(0);}
  
  //! Keeps only air broadening but throws if mbath is false
  void KeepOnlyBath() {if(not mbath) throw std::runtime_error("Cannot comply because no air broadening exist!"); else while(mdata.nelem() > 1) Remove(0);}

private:
  /** Should ANY normal calculations be performed? */
  bool do_line_in_standard_calculations;
  
  /** Do we have self-broadening? */
  bool mself;
  
  /** Do we have air-broadening? */
  bool mbath;
  
  /** List of species.  Note that first is ill-formed if mself, and last is ill-formed if mbath */
  ArrayOfSpeciesTag mspecies;
  
  /** List of temperature types per species, outer size is species, inner size is number of line shape and line mixing parameters in that order */
  Array<Array<TemperatureType>> mtypes;
  
  /** List of data per species, the internal vector is the sum of all temperature parameters required for line shape and line mixing computations */
  ArrayOfVector mdata;
  
  /** List of errors per species, same size as mdata or empty */
  ArrayOfVector merrors;
  
  LineShapeType mp;
  
  LineMixingOrderType mlm;
  
  // Model has data
  bool ComputesParam(const String& type) const noexcept;
  Index IndexOfParam(const String& type) const noexcept;
};

std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd);
std::istream& operator>>(std::istream& is, LineFunctionData& lfd);

//! Select the derivative that will be used in Jacobian calculations --- also checks validity of var and coeff
JacPropMatType select_derivativeLineFunctionData(const String& var, const String& coeff);

//! {"X0", "X1", "X2"}
ArrayOfString all_coefficientsLineFunctionData();

//! {"G0", "D0", "G2", "D2", "ETA", "FVC", "Y", "G", "DV"}
ArrayOfString all_variablesLineFunctionData();

constexpr LineFunctionDataOutput NoLineFunctionDataOutput(){return {0, 0, 0, 0, 0, 0, 0, 0, 0};}

LineFunctionDataOutput mirroredOutput(LineFunctionDataOutput v) noexcept;

LineFunctionDataOutput si2cgs(LineFunctionDataOutput v);

LineFunctionDataOutput cgs2si(LineFunctionDataOutput v);

#endif // linefunctiondata_h

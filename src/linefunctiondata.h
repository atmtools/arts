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
typedef struct{Numeric G0, D0, G2, D2, FVC, ETA, Y, G, DV;} LineFunctionDataOutput;

class LineFunctionData {
public:

  
  #define LineFunctionData_SelfBroadening String("SELF")
  #define LineFunctionData_BathBroadening String("AIR")
  #define LineFunctionData_NoLineMixing String("#")
  
  // Line shape models
  enum class LineShapeType {DP, LP, VP, SDVP, HTP, };
  
  // Line mixing models
  enum class LineMixingOrderType {None, LM1, LM2, ConstG, Interp, };
  
  // Temperature dependencies
  enum class TemperatureType {None, T0, T1, T2, T3, T4, T5, LM_AER };
  
  // Local variables allowed in the temperature functions
  enum class VariableType {X0, X1, X2, };
  
  // Line shape models are defined by these variables (do leave Size at end)
  enum class DopplerParam    : Index {                          Size};
  enum class LorentzParam    : Index {G0, D0,                   Size};
  enum class VoigtParam      : Index {G0, D0,                   Size};
  enum class SpeedVoigtParam : Index {G0, D0, G2, D2,           Size};
  enum class HTPParam        : Index {G0, D0, G2, D2, FVC, ETA, Size};
  
  // Line mixing models are defined by these variables (do leave Size at end)
  enum class NoLineMixingParam : Index {                        Size};
  enum class FirstOrderParam   : Index {Y,                      Size};
  enum class SecondOrderParam  : Index {Y, G, DV,               Size};
  enum class ConstGParam       : Index {   G,                   Size};
  enum class InterpParam       : Index {INTERPOLATED_VARIABLES, Size};
  
  LineFunctionData() = default;
  
  LineFunctionData(const Numeric& sgam, const Numeric& nself, const Numeric& agam, const Numeric& nair, const Numeric& psf) :
  do_line_in_standard_calculations(true), mself(true), mbath(true), mspecies(2), mtypes(2, {TemperatureType::T1, TemperatureType::T5}),
  mdata(2), merrors(2, Vector(4, 0)), mp(LineShapeType::VP), mlm(LineMixingOrderType::None)
  {
    mdata[0] = {sgam, nself, psf, nair};
    mdata[1] = {agam, nair, psf, nair};
  };
  
  //! TRANSLATION MECHANISM TO BE DEPRECATED ASAP
  LineFunctionData(const PressureBroadeningData& pb, const LineMixingData& lm, const String& species, const Numeric& T0);
  
  // Identification calls
  Index LineMixingTypeNelem() const {
    switch(mlm) {
      case LineMixingOrderType::None:   return Index(NoLineMixingParam::Size);
      case LineMixingOrderType::LM1:    return Index(FirstOrderParam::Size);
      case LineMixingOrderType::LM2:    return Index(SecondOrderParam::Size);
      case LineMixingOrderType::Interp: return Index(InterpParam::Size);
      case LineMixingOrderType::ConstG: return Index(ConstGParam::Size);;
    }

    return -1;  // Should not reach
  }
  
  String LineMixingType2String() const {
    switch(mlm) {
      case LineMixingOrderType::None:   return LineFunctionData_NoLineMixing;
      case LineMixingOrderType::LM1:    return "LM1";
      case LineMixingOrderType::LM2:    return "LM2";
      case LineMixingOrderType::Interp: return "INT";
      case LineMixingOrderType::ConstG: return "ConstG";
    }

    return "-1";
  }
  
  void StringSetLineMixingType(const String& type) {
    if(type == LineFunctionData_NoLineMixing) mlm = LineMixingOrderType::None; 
    else if(type == String("LM1"))            mlm = LineMixingOrderType::LM1;
    else if(type == String("LM2"))            mlm = LineMixingOrderType::LM2;
    else if(type == String("INT"))            mlm = LineMixingOrderType::Interp;
    else if(type == String("ConstG"))         mlm = LineMixingOrderType::ConstG;
    else {
      ostringstream os;
      os << "Cannot recognize type " << type << " as a line mixing functionality\n" ;
      throw std::runtime_error(os.str());
    }

  }
  
  Index LineShapeTypeNelem() const {
    switch(mp) {
      case LineShapeType::DP:   return Index(DopplerParam::Size);
      case LineShapeType::LP:   return Index(LorentzParam::Size);
      case LineShapeType::VP:   return Index(VoigtParam::Size);
      case LineShapeType::SDVP: return Index(SpeedVoigtParam::Size);
      case LineShapeType::HTP:  return Index(HTPParam::Size);
    }

    return -1;  // Should not reach
  }
  
  String LineShapeType2String() const {
    switch(mp) {
      case LineShapeType::DP:   return "DP";
      case LineShapeType::LP:   return "LP";
      case LineShapeType::VP:   return "VP";
      case LineShapeType::SDVP: return "SDVP";
      case LineShapeType::HTP:  return "HTP";
    }

    return "-1";
  }
  
  void StringSetLineShapeType(const String& type) {
    if(type == String("DP"))        mp = LineShapeType::DP;
    else if(type == String("LP"))   mp = LineShapeType::LP;
    else if(type == String("VP"))   mp = LineShapeType::VP;
    else if(type == String("SDVP")) mp = LineShapeType::SDVP;
    else if(type == String("HTP"))  mp = LineShapeType::HTP;
    else {
      ostringstream os;
      os << "Cannot recognize type " << type << " as a line shape functionality\n" ;
      throw std::runtime_error(os.str());
    }

  }
  
  Index TemperatureTypeNelem(const TemperatureType type) const {
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

    return -1;
  }
  
  Index LineShapeDataNelemForSpecies(Index ispec) const {
    Index count=0;
    for(auto i=0; i<LineShapeTypeNelem(); i++)
      count += TemperatureTypeNelem(mtypes[ispec][i]);
    return count;
  }
  
  String TemperatureType2String(const TemperatureType type) const {
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

    return "-1";
  }
  
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
      ostringstream os;
      os << "Cannot recognize type " << type << " as a type of temperature type\n" ;
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
                                   const bool normalization=true) const;

  LineFunctionDataOutput GetTemperatureDerivs(const Numeric& T0,
                                              const Numeric& T,
                                              const Numeric& dT,
                                              const Numeric& P,
                                              const Numeric& self_vmr,
                                              const ConstVectorView& rtp_vmr,
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const bool do_linemixing=true,
                                              const bool normalization=true) const;

  LineFunctionDataOutput GetReferenceT0Derivs(const Numeric& T0,
                                              const Numeric& T,
                                              const Numeric& P,
                                              const Numeric& self_vmr,
                                              const ConstVectorView& rtp_vmr, 
                                              const ArrayOfArrayOfSpeciesTag& abs_species,
                                              const RetrievalQuantity& rt, 
                                              const QuantumIdentifier& line_qi,
                                              const bool do_linemixing=true,
                                              const bool normalization=true) const;

  LineFunctionDataOutput GetVMRDerivs(const Numeric& T0,
                                      const Numeric& T,
                                      const Numeric& P,
                                      const Numeric& self_vmr,
                                      const ConstVectorView& rtp_vmr,
                                      const ArrayOfArrayOfSpeciesTag& abs_species,
                                      const QuantumIdentifier& vmr_qi, 
                                      const QuantumIdentifier& line_qi,
                                      const bool do_linemixing=true,
                                      const bool normalization=true) const;
              
  Numeric GetLineParamDeriv(const Numeric& T0,
                            const Numeric& T,
                            const Numeric& P,
                            const Numeric& self_vmr,
                            const ConstVectorView& rtp_vmr, 
                            const ArrayOfArrayOfSpeciesTag& abs_species,
                            const RetrievalQuantity& rt, 
                            const QuantumIdentifier& line_qi,
                            const bool do_linemixing=true,
                            const bool normalization=true) const;
                            
  // Read data
  const ArrayOfSpeciesTag& Species() const { return mspecies; }
  bool SelfBroadening() const { return mself; }
  bool AirBroadening() const { return mbath; }
  bool DoStandardCalcs() const { return do_line_in_standard_calculations; }
  LineShapeType LineShape() const { return mp; }
  LineMixingOrderType LineMixing() const { return mlm; }
  
  friend std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd);
  friend std::istream& operator>>(std::istream& is, LineFunctionData& lfd);
  
  // Special functions translating old behavior, all DEPRECATED
  Numeric SelfG0() const;
  Numeric SelfN() const;
  Numeric AirG0() const;
  Numeric AirD0() const;
  Numeric AirN() const;
  Numeric dSelfG0() const;
  Numeric dSelfN() const;
  Numeric dAirG0() const;
  Numeric dAirD0() const;
  Numeric dAirN() const;
  Vector PlanetaryForeignG0() const;
  Vector PlanetaryForeignD0() const;
  Vector PlanetaryForeignN() const;
  LineFunctionDataOutput AirBroadening(const Numeric& theta, const Numeric& P, const Numeric& self_vmr) const;
  
  // Changes and alterations of internal data
  void ChangeLineMixing(const LineMixingOrderType lm, const Array<TemperatureType>& ts, const Vector& data);
  void ChangeLineMixingfromSimpleLM2(const Vector& data) {ChangeLineMixing(LineMixingOrderType::LM2, {TemperatureType::T4, TemperatureType::T4, TemperatureType::T4}, data);}
  void ChangeLineMixing2AER(const Vector& data) {ChangeLineMixing(LineMixingOrderType::Interp, {TemperatureType::LM_AER}, data);}
  void SetStandard(const bool standard) {do_line_in_standard_calculations=standard;}
  void Set(const Numeric& X, const String& species, const String& coefficient, const String& variable);
  Numeric Get(const String& species, const String& coefficient, const String& variable) const;
    void Remove(Index);
    void RemoveSelf() {if(mself) Remove(0);}
    void KeepOnlyBath() {if(not mbath) throw std::runtime_error("Cannot comply because no air broadening exist!"); else while(mdata.nelem() > 1) Remove(0);}

private:
  bool do_line_in_standard_calculations;
  bool mself;
  bool mbath;
  ArrayOfSpeciesTag mspecies;
  Array<Array<TemperatureType>> mtypes;
  ArrayOfVector mdata;
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

LineFunctionDataOutput mirroredOutput(const LineFunctionDataOutput& v) noexcept;

#endif // linefunctiondata_h

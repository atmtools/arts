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

/** Contains the pressure broadening data class
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

class LineFunctionData {
public:
  
  #define SELF "SELF"
  #define BATH "AIR"
  #define NO_LM "#"
  
  // Line shape models
  enum class LineShapeType {DP, LP, VP, SDVP, HTP, };
  
  // Line mixing models
  enum class LineMixingOrderType {None, LM1, LM2, Interp, };
  
  // Temperature dependencies
  enum class TemperatureType {None, T0, T1, T2, T3, T4, T5, LM_AER, };
  
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
  enum class InterpParam       : Index {INTERPOLATED_VARIABLES, Size};
  
  LineFunctionData() = default;
  
   //! TRANSLATION MECHANISM TO BE DEPRECATED ASAP
   LineFunctionData(const PressureBroadeningData& pb, const LineMixingData& lm, const String& species, const Numeric& T0);
  
  Index LineMixingTypeNelem() const {
    switch(mlm) {
      case LineMixingOrderType::None:   return Index(NoLineMixingParam::Size);
      case LineMixingOrderType::LM1:    return Index(FirstOrderParam::Size);
      case LineMixingOrderType::LM2:    return Index(SecondOrderParam::Size);
      case LineMixingOrderType::Interp: return Index(InterpParam::Size);
    }
    return -1;  // Should not reach
  }
  
  String LineMixingType2String() const {
    switch(mlm) {
      case LineMixingOrderType::None:   return NO_LM;
      case LineMixingOrderType::LM1:    return "LM1";
      case LineMixingOrderType::LM2:    return "LM2";
      case LineMixingOrderType::Interp: return "INT";
    }
    return "-1";
  }
  
  void StringSetLineMixingType(const String& type) {
    if(type == NO_LM)      mlm = LineMixingOrderType::None; 
    else if(type == "LM1") mlm = LineMixingOrderType::LM1;
    else if(type == "LM2") mlm = LineMixingOrderType::LM2;
    else if(type == "INT") mlm = LineMixingOrderType::Interp;
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
    if(type == "DP")        mp = LineShapeType::DP;
    else if(type == "LP")   mp = LineShapeType::LP;
    else if(type == "VP")   mp = LineShapeType::VP;
    else if(type == "SDVP") mp = LineShapeType::SDVP;
    else if(type == "HTP")  mp = LineShapeType::HTP;
    else {
      ostringstream os;
      os << "Cannot recognize type " << type << " as a line shape functionality\n" ;
      throw std::runtime_error(os.str());
    }
  }
  
  Index TemperatureTypeNelem(const TemperatureType type) const {
    switch(type) {
      case TemperatureType::LM_AER: return -1;
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
  
  String TemperatureType2String(const TemperatureType type) const {
    switch(type) {
      case TemperatureType::LM_AER: return "LM_AER";
      case TemperatureType::None:   return NO_LM;
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
    if(type == "LM_AER")   mtypes[ispecies][iparam] = TemperatureType::LM_AER;
    else if(type == NO_LM) mtypes[ispecies][iparam] = TemperatureType::None;
    else if(type == "T0")  mtypes[ispecies][iparam] = TemperatureType::T0;
    else if(type == "T1")  mtypes[ispecies][iparam] = TemperatureType::T1;
    else if(type == "T2")  mtypes[ispecies][iparam] = TemperatureType::T2;
    else if(type == "T3")  mtypes[ispecies][iparam] = TemperatureType::T3;
    else if(type == "T4")  mtypes[ispecies][iparam] = TemperatureType::T4;
    else if(type == "T5")  mtypes[ispecies][iparam] = TemperatureType::T5;
    else {
      ostringstream os;
      os << "Cannot recognize type " << type << " as a type of temperature type\n" ;
      throw std::runtime_error(os.str());
    }
  }
  
  void GetParams(Numeric& G0, Numeric& D0, 
                 Numeric& G2, Numeric& D2, 
                 Numeric& FVC, Numeric& ETA,
                 Numeric& Y, Numeric& G, Numeric& DV,
                 const Numeric& T0, const Numeric& T,
                 const Numeric& P, const Numeric& self_vmr,
                 const ConstVectorView rtp_vmr, const ArrayOfSpeciesTag& abs_species,
                 const bool normalization=true) const;
  
  void GetTemperatureDerivs(Numeric& dG0, Numeric& dD0, 
                            Numeric& dG2, Numeric& dD2, 
                            Numeric& dFVC, Numeric& dETA,
                            Numeric& dY, Numeric& dG, Numeric& dDV,
                            const Numeric& T0, const Numeric& T, const Numeric& dT,
                            const Numeric& P, const Numeric& self_vmr,
                            const ConstVectorView rtp_vmr, 
                            const ArrayOfSpeciesTag& abs_species,
                            const bool normalization=true) const;
  
  void GetReferenceT0Derivs(Numeric& dG0, Numeric& dD0, 
                            Numeric& dG2, Numeric& dD2, 
                            Numeric& dFVC, Numeric& dETA,
                            Numeric& dY, Numeric& dG, Numeric& dDV,
                            const Numeric& T0, const Numeric& T,
                            const Numeric& P, const Numeric& self_vmr,
                            const ConstVectorView rtp_vmr, 
                            const ArrayOfSpeciesTag& abs_species,
                            const RetrievalQuantity& rt, 
                            const QuantumIdentifier& line_qi,
                            const bool normalization=true) const;
  
  Numeric GetLineParamDeriv(const Numeric& T0, const Numeric& T,
                            const Numeric& P, const Numeric& self_vmr,
                            const ConstVectorView rtp_vmr, 
                            const ArrayOfSpeciesTag& abs_species,
                            const RetrievalQuantity& rt, 
                            const QuantumIdentifier& line_qi,
                            const bool normalization) const;

  // Model has data
  bool ComputesParam(const String& type) const;
                            
  // Read data
  const ArrayOfVector& Data() const { return mdata; }
  const ArrayOfSpeciesTag& Species() const { return mspecies; }
  bool SelfBroadening() const { return mself; }
  bool AirBroadening() const { return mbath; }
  bool DoLine() const { return do_line_in_standard_calculations; }
  LineShapeType LineShape() const { return mp; }
  LineMixingOrderType LineMixing() const { return mlm; }
  
  friend std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd);
  friend std::istream& operator>>(std::istream& is, LineFunctionData& lfd);
  
private:
  bool do_line_in_standard_calculations;
  bool mself;
  bool mbath;
  ArrayOfSpeciesTag mspecies;
  Array<Array<TemperatureType>> mtypes;
  ArrayOfVector mdata;
  LineShapeType mp;
  LineMixingOrderType mlm;
};

std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd);
std::istream& operator>>(std::istream& is, LineFunctionData& lfd);

#endif // linefunctiondata_h

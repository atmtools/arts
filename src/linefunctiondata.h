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

// Using some algorithms
#include <numeric>
#include "constants.h"

#ifndef NEWARTSCAT

// List of all variables that can be returned
struct LineFunctionDataOutput{
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
};
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
  mdata(2), mp(LineShapeType::VP), mlm(LineMixingOrderType::None) {
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
  
  LineShapeType mp;
  
  LineMixingOrderType mlm;
  
  // Model has data
  bool ComputesParam(const String& type) const noexcept;
  Index IndexOfParam(const String& type) const noexcept;
};

std::ostream& operator<<(std::ostream& os, const LineFunctionData& lfd);
std::istream& operator>>(std::istream& is, LineFunctionData& lfd);

constexpr LineFunctionDataOutput NoLineFunctionDataOutput(){return {0, 0, 0, 0, 0, 0, 0, 0, 0};}

LineFunctionDataOutput mirroredOutput(LineFunctionDataOutput v) noexcept;

LineFunctionDataOutput si2cgs(LineFunctionDataOutput v);

LineFunctionDataOutput cgs2si(LineFunctionDataOutput v);
#endif


//! Select the derivative that will be used in Jacobian calculations --- also checks validity of var and coeff
JacPropMatType select_derivativeLineShape(const String& var, const String& coeff);


//! {"X0", "X1", "X2"}
ArrayOfString all_coefficientsLineFunctionData();


//! {"G0", "D0", "G2", "D2", "ETA", "FVC", "Y", "G", "DV"}
ArrayOfString all_variablesLineFunctionData();


namespace LineShape {
  
  /** Temperature models */
  enum class TemperatureModel : Index {
    None,   // 0
    T0,     // Constant, X0
    T1,     // Standard, X0 * (T0/T) ^ X1
    T2,     // X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0));
    T3,     // X0 + X1 * (T - T0)
    T4,     // (X0 + X1 * (T0/T - 1)) * (T0/T)^X2;
    T5,     // X0 * (T0/T)^(0.25 + 1.5*X1)
    LM_AER  // Interpolation AER-style
    // ALWAYS ADD NEW AT THE END
  };
  
  inline String temperaturemodel2string(TemperatureModel type) noexcept {
    switch(type) {
      case TemperatureModel::None:   return "#";
      case TemperatureModel::T0:     return "T0";
      case TemperatureModel::T1:     return "T1";
      case TemperatureModel::T2:     return "T2";
      case TemperatureModel::T3:     return "T3";
      case TemperatureModel::T4:     return "T4";
      case TemperatureModel::T5:     return "T5";
      case TemperatureModel::LM_AER: return "LM_AER";
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  inline TemperatureModel string2temperaturemodel(const String& type) {
    if(type == "#")                   return TemperatureModel::None;
    else if(type == String("T0"))     return TemperatureModel::T0;
    else if(type == String("T1"))     return TemperatureModel::T1;
    else if(type == String("T2"))     return TemperatureModel::T2;
    else if(type == String("T3"))     return TemperatureModel::T3;
    else if(type == String("T4"))     return TemperatureModel::T4;
    else if(type == String("T5"))     return TemperatureModel::T5;
    else if(type == String("LM_AER")) return TemperatureModel::LM_AER;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  "
         << "See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }
  
  /** List of possible shape variables */
  enum class Variable {
    G0=0,   // Pressure broadening speed-independent
    D0=1,   // Pressure f-shifting speed-dependent
    G2=2,   // Pressure broadening speed-dependent
    D2=3,   // Pressure f-shifting speed-independent
    FVC=4,  // Frequency of velocity-changing collisions
    ETA=5,  // Correlation
    Y=6,    // First order line mixing coefficient
    G=7,    // Second order line mixing coefficient
    DV=8    // Second order line mixing f-shifting
    // ALWAYS ADD NEW AT THE END
  };
  
  inline std::ostream& operator<<(std::ostream& os, Variable v) {
    switch(v) {
      case Variable::G0:  os << "G0";  break;
      case Variable::D0:  os << "D0";  break;
      case Variable::G2:  os << "G2";  break;
      case Variable::D2:  os << "D2";  break;
      case Variable::FVC: os << "FVC"; break;
      case Variable::ETA: os << "ETA"; break;
      case Variable::Y:   os << "Y";   break;
      case Variable::G:   os << "G";   break;
      case Variable::DV:  os << "DV";  break;
    }
    return os;
  }
  
  inline String variable2string(Variable type) noexcept {
    #define VARIABLE2STRINGDEF(X) case Variable::X: return # X
    switch(type) {
      VARIABLE2STRINGDEF(G0);
      VARIABLE2STRINGDEF(D0);
      VARIABLE2STRINGDEF(G2);
      VARIABLE2STRINGDEF(D2);
      VARIABLE2STRINGDEF(FVC);
      VARIABLE2STRINGDEF(ETA);
      VARIABLE2STRINGDEF(Y);
      VARIABLE2STRINGDEF(G);
      VARIABLE2STRINGDEF(DV);
    }
    #undef VARIABLE2STRINGDEF
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  inline Variable string2variable(const String& type) {
    #define STRING2VARIABLEDEF(X) if(type == # X) return Variable::X
    STRING2VARIABLEDEF(G0);
    else STRING2VARIABLEDEF(D0);
    else STRING2VARIABLEDEF(G2);
    else STRING2VARIABLEDEF(D2);
    else STRING2VARIABLEDEF(FVC);
    else STRING2VARIABLEDEF(ETA);
    else STRING2VARIABLEDEF(Y);
    else STRING2VARIABLEDEF(G);
    else STRING2VARIABLEDEF(DV);
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  "
      << "See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }
  
  struct ModelParameters {
    TemperatureModel type;
    Numeric X0;
    Numeric X1;
    Numeric X2;
    // ALWAYS ADD NEW AT THE END
  };
  
  inline Numeric& SingleModelParameter(ModelParameters& mp, const String& type) {
    if(type =="X0")
      return mp.X0;
    else if(type =="X1")
      return mp.X1;
    else if(type =="X2")
      return mp.X2;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  "
      << "See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
    std::terminate();  
  }
  
  inline std::ostream& operator<<(std::ostream& os, const ModelParameters& mp) {
    os << temperaturemodel2string(mp.type) << ' ' << mp.X0 << ' ' << mp.X1 << ' ' << mp.X2;
    return os;
  }
  
  inline std::istream& operator>>(std::istream& is, ModelParameters& mp) {
    String tmp;
    is >> tmp >> mp.X0 >> mp.X1 >> mp.X2;
    mp.type = string2temperaturemodel(tmp);
    return is;
  }

  constexpr Index nmaxTempModelParams=3;
  constexpr Index nVars=9;
  constexpr Index nmaxInterpModels=12;

  class SingleSpeciesModel {
  private:
    std::array<ModelParameters, nVars> X;
    std::array<Numeric, nmaxInterpModels> V;

    Numeric special_linemixing_aer(Numeric T,
                                   Variable var) const noexcept
    {
      // Data starts at 4 for Y and 8 for G, and must be either
      const Index i = (var == Variable::Y) ? 4 : 8;
      
      if(T < V[1])
        return V[i+0] + (T-V[0]) * (V[i+1]-V[i+0]) / (V[1]-V[0]);
      else if(T > V[2])
        return V[i+2] + (T-V[2]) * (V[i+3]-V[i+2]) / (V[3]-V[2]);
      else
        return V[i+1] + (T-V[1]) * (V[i+2]-V[i+1]) / (V[2]-V[1]);
    }
    
    Numeric special_linemixing_aer_dT(Numeric T,
                                      Variable var) const noexcept
    {
      // Data starts at 4 for Y and 8 for G, and must be either
      const Index i = (var == Variable::Y) ? 4 : 8;
      
      if(T < V[1])
        return (V[i+1]-V[i+0]) / (V[1]-V[0]);
      else if(T > V[2])
        return (V[i+3]-V[i+2]) / (V[3]-V[2]);
      else
        return (V[i+2]-V[i+1]) / (V[2]-V[1]);
    }
    
  public:
    constexpr
    SingleSpeciesModel(ModelParameters G0  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters D0  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters G2  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters D2  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters FVC = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters ETA = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters Y   = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters G   = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters DV  = {TemperatureModel::None, 0, 0, 0},
                       std::array<Numeric, nmaxInterpModels> interp = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) :
    X({G0, D0, G2, D2, FVC, ETA, Y, G, DV}), V(interp)
    {}
    
    #define x0 X[Index(var)].X0
    #define x1 X[Index(var)].X1
    #define x2 X[Index(var)].X2
    
    // Standard compute feature
    Numeric compute(Numeric T, Numeric T0, Variable var) const noexcept {
      using std::pow;
      using std::log;
      
      switch(X[Index(var)].type) {
        case TemperatureModel::None:
          return 0;
        case TemperatureModel::T0:
          return x0;
        case TemperatureModel::T1:
          return x0 * pow(T0/T, x1);
        case TemperatureModel::T2:
          return x0 * pow(T0/T, x1) * (1 + x2 * log(T/T0));
        case TemperatureModel::T3:
          return x0 + x1 * (T - T0);
        case TemperatureModel::T4:
          return (x0 + x1 * (T0/T - 1.)) * pow(T0/T, x2);
        case TemperatureModel::T5:
          return x0 * pow(T0/T, 0.25 + 1.5*x1);
        case TemperatureModel::LM_AER:
          return special_linemixing_aer(T, var);
      }
      std::terminate();
    }
    
    // Standard compute feature for derivative wrt x0
    Numeric compute_dX0(Numeric T, Numeric T0, Variable var) const noexcept {
      using std::pow;
      using std::log;
      
      switch(X[Index(var)].type) {
        case TemperatureModel::None:
          return 0;
        case TemperatureModel::T0:
          return 1;
        case TemperatureModel::T1:
          return pow(T0/T, x1);
        case TemperatureModel::T2:
          return pow(T0/T, x1)*(x2*log(T/T0) + 1);
        case TemperatureModel::T3:
          return 1;
        case TemperatureModel::T4:
          return pow(T0/T, x2);
        case TemperatureModel::T5:
          return pow(T0/T, 1.5*x1 + 0.25);
        case TemperatureModel::LM_AER:
          return 0;
      }
      std::terminate();
    }
    
    // Standard compute feature for derivative wrt x1
    Numeric compute_dX1(Numeric T, Numeric T0, Variable var) const noexcept {
      using std::pow;
      using std::log;
      
      switch(X[Index(var)].type) {
        case TemperatureModel::None:
          return 0;
        case TemperatureModel::T0:
          return 0;
        case TemperatureModel::T1:
          return x0*pow(T0/T, x1)*log(T0/T);
        case TemperatureModel::T2:
          return x0*pow(T0/T,x1)*(x2*log(T/T0)+1.)*log(T0/T);
        case TemperatureModel::T3:
          return (T - T0); 
        case TemperatureModel::T4:
          return pow(T0/T, x2)*(T0/T - 1.);
        case TemperatureModel::T5:
          return 1.5*x0*pow(T0/T, 1.5*x1 + 0.25)*log(T0/T);
        case TemperatureModel::LM_AER:
          return 0;
      }
      std::terminate();
    }
    
    // Standard compute feature for derivative wrt x2
    Numeric compute_dX2(Numeric T, Numeric T0, Variable var) const noexcept {
      using std::pow;
      using std::log;
      
      switch(X[Index(var)].type) {
        case TemperatureModel::None:
          return 0;
        case TemperatureModel::T0:
          return 0;
        case TemperatureModel::T1:
          return 0;
        case TemperatureModel::T2:
          return x0*pow(T0/T, x1)*log(T/T0);
        case TemperatureModel::T3:
          return 0; 
        case TemperatureModel::T4:
          return pow(T0/T,x2)*(x0+x1*(T0/T-1))*log(T0/T);
        case TemperatureModel::T5:
          return 0;
        case TemperatureModel::LM_AER:
          return 0;
      }
      std::terminate();
    }
    
    // Standard compute feature for derivative wrt T
    Numeric compute_dT(Numeric T, Numeric T0, Variable var) const noexcept {
      using std::pow;
      using std::log;
      
      switch(X[Index(var)].type) {
        case TemperatureModel::None:
          return 0;
        case TemperatureModel::T0:
          return 0;
        case TemperatureModel::T1:
          return -x0*x1*pow(T0/T, x1)/T;
        case TemperatureModel::T2:
          return -x0*x1*pow(T0/T,x1)*(x2*log(T/T0)+1.)/T+x0*x2*pow(T0/T,x1)/T;
        case TemperatureModel::T3:
          return x1;
        case TemperatureModel::T4:
          return -x2*pow(T0/T,x2)*(x0+x1*(T0/T-1.))/T-T0*x1*pow(T0/T,x2)/pow(T,2);
        case TemperatureModel::T5:
          return -x0*pow(T0/T,1.5*x1 + 0.25)*(1.5*x1 + 0.25)/T;
        case TemperatureModel::LM_AER:
          return special_linemixing_aer_dT(T, var);
      }
      std::terminate();
    }
    
    // Standard compute feature for derivative wrt T0
    Numeric compute_dT0(Numeric T, Numeric T0, Variable var) const noexcept {
      using std::pow;
      using std::log;
      
      switch(X[Index(var)].type) {
        case TemperatureModel::None:
          return 0;
        case TemperatureModel::T0:
          return 0;
        case TemperatureModel::T1:
          return x0*x1*pow(T0/T, x1)/T0;
        case TemperatureModel::T2:
          return x0*x1*pow(T0/T,x1)*(x2*log(T/T0)+1.)/T0-x0*x2*pow(T0/T,x1)/T0;
        case TemperatureModel::T3:
          return -x1;
        case TemperatureModel::T4:
          return x2*pow(T0/T,x2)*(x0+x1*(T0/T-1.))/T0+x1*pow(T0/T,x2)/T;
        case TemperatureModel::T5:
          return x0*pow(T0/T,1.5*x1 + 0.25)*(1.5*x1 + 0.25)/T0;
        case TemperatureModel::LM_AER:
          return 0;
      }
      std::terminate();
    }

    #undef x0
    #undef x1
    #undef x2
    
    // Access normal data
    #define ACCESS_INTERNAL(VARPOS)                                               \
    ModelParameters& VARPOS()       noexcept {return X[Index(Variable::VARPOS)];} \
    ModelParameters  VARPOS() const noexcept {return X[Index(Variable::VARPOS)];}
    ACCESS_INTERNAL(G0); ACCESS_INTERNAL(G2);
    ACCESS_INTERNAL(D0); ACCESS_INTERNAL(D2);
    ACCESS_INTERNAL(FVC); ACCESS_INTERNAL(ETA);
    ACCESS_INTERNAL(Y); ACCESS_INTERNAL(G); ACCESS_INTERNAL(DV);
    #undef ACCESS_INTERNAL
    
    // Access to all data
    std::array<ModelParameters, nVars>&       Data()       noexcept {return X;}
    const std::array<ModelParameters, nVars>& Data() const noexcept {return X;}
    std::array<Numeric, nmaxInterpModels>&       Interp()       noexcept {return V;}
    const std::array<Numeric, nmaxInterpModels>& Interp() const noexcept {return V;}
    
    void Set(Variable var, const ModelParameters& x) noexcept {
      #define MODELPARAMCASESETTER(X) case Variable::X: X () = x; break
      switch(var) {
        MODELPARAMCASESETTER(G0);
        MODELPARAMCASESETTER(D0);
        MODELPARAMCASESETTER(G2);
        MODELPARAMCASESETTER(D2);
        MODELPARAMCASESETTER(FVC);
        MODELPARAMCASESETTER(ETA);
        MODELPARAMCASESETTER(Y);
        MODELPARAMCASESETTER(G);
        MODELPARAMCASESETTER(DV);
      }
      #undef MODELPARAMCASESETTER
    }
    
    ModelParameters Get(Variable var) const noexcept {
      #define MODELPARAMCASEGETTER(X) case Variable::X: return X ();
      switch(var) {
        MODELPARAMCASEGETTER(G0);
        MODELPARAMCASEGETTER(D0);
        MODELPARAMCASEGETTER(G2);
        MODELPARAMCASEGETTER(D2);
        MODELPARAMCASEGETTER(FVC);
        MODELPARAMCASEGETTER(ETA);
        MODELPARAMCASEGETTER(Y);
        MODELPARAMCASEGETTER(G);
        MODELPARAMCASEGETTER(DV);
      }
      #undef MODELPARAMCASEGETTER
      std::terminate();
    }
    
    // Read complete binary data... FIXME? adopt for versioning?  constexpr numbers above are assumed
    std::istream& read(std::istream& is) {
      is.read(reinterpret_cast<char*>(this), sizeof(SingleSpeciesModel));
      return is;
    }
    
    // Write complete binary data
    std::ostream& write(std::ostream& os) const {
      os.write(reinterpret_cast<const char*>(this), sizeof(SingleSpeciesModel));
      return os;
    }
  };
  
  inline std::ostream& operator<<(std::ostream& os, const SingleSpeciesModel& ssm) {
    for(const auto& mp: ssm.Data()) os << mp << ' ';
    for(const auto& num: ssm.Interp()) os << num << ' ';
    return os;
  }
  
  inline std::istream& operator>>(std::istream& is, SingleSpeciesModel& ssm) {
    for(auto& mp: ssm.Data()) is >> mp;
    for(auto& num: ssm.Interp()) is >> num;
    return is;
  }
  
  enum class Type {
    DP,   // Doppler
    LP,   // Lorentz
    VP,   // Voigt
    SDVP, // Speed-dependent Voigt
    HTP,  // Hartmann-Tran
  };
  
  inline String shapetype2string(Type type) noexcept {
    switch(type) {
      case Type::DP:   return "DP";
      case Type::LP:   return "LP";
      case Type::VP:   return "VP";
      case Type::SDVP: return "SDVP";
      case Type::HTP:  return "HTP";
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  inline Type string2shapetype(const String& type) {
    if(type == "DP")                return Type::DP;
    else if(type == String("LP"))   return Type::LP;
    else if(type == String("VP"))   return Type::VP;
    else if(type == String("SDVP")) return Type::SDVP;
    else if(type == String("HTP"))  return Type::HTP;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  "
      << "See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }
  
  struct Output {
    Numeric
    G0,   // Pressure broadening speed-independent
    D0,   // Pressure f-shifting speed-independent
    G2,   // Pressure broadening speed-dependent
    D2,   // Pressure f-shifting speed-dependent
    FVC,  // Frequency of velocity-changing collisions
    ETA,  // Correlation
    Y,    // First order line mixing coefficient
    G,    // Second order line mixing coefficient
    DV;   // Second order line mixing f-shifting
  };
  
  inline std::ostream& operator<<(std::ostream& os, Output x) {
    return os <<  "G0: " << x.G0
              << " D0: " << x.D0
              << " G2: " << x.G2
              << " D2: " << x.D2
              << " FVC: " << x.FVC
              << " ETA: " << x.ETA
              << " Y: " << x.Y
              << " G: " << x.G
              << " DV: " << x.DV;
  }
  
  constexpr Output mirroredOutput(Output x) noexcept {
    return {x.G0, -x.D0, x.G2, -x.D2, x.FVC, x.ETA, x.Y, x.G, -x.DV};
  }
  
  
  constexpr Output negativeOutput(Output x) noexcept {
    return {-x.G0, -x.D0, -x.G2, -x.D2, -x.FVC, -x.ETA, -x.Y, -x.G, -x.DV};
  }
  
  
  constexpr Output si2cgs(Output x) noexcept {
    using Conversion::freq2kaycm;
    return {freq2kaycm(x.G0), freq2kaycm(x.D0),
            freq2kaycm(x.D2), freq2kaycm(x.G2),
            freq2kaycm(x.FVC), x.ETA,
            x.Y, x.G, freq2kaycm(x.DV)};
  }

  static constexpr const char*const bath_broadening = "AIR";
  static constexpr const char*const self_broadening = "SELF";
  
  class Model {
  private:
    Type mtype;
    bool mself;
    bool mbath;
    ArrayOfSpeciesTag mspecies;
    std::vector<SingleSpeciesModel> mdata;
    
    bool OK() const noexcept {
      Index n = mdata.size();
      Index k = mspecies.nelem();
      Index m = Index(mself) + Index(mbath);
      bool needs_any = mtype not_eq Type::DP;
      if(n not_eq k or m > n or (needs_any and not n))
        return false;
      else
        return true;
    }
  public:
    // No line shape parameters means this is a Doppler line
    Model() noexcept : mtype(Type::DP), mself(false), mbath(false),
                       mspecies(0), mdata(0) {}
    
    // Standard HITRAN means this is a Voigt line (with interp it is lblrtm line mixing)
    Model(Numeric sgam, Numeric nself, Numeric agam, Numeric nair, Numeric psf,
          std::array<Numeric, nmaxInterpModels> interp = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) noexcept :
    mtype(Type::VP), mself(true), mbath(true), mspecies(2), mdata(2) {
      mdata.front().G0() = {TemperatureModel::T1, sgam, nself, 0};
      mdata.front().D0() = {TemperatureModel::T5, psf,  nair,  0};
      mdata.front().Interp() = interp;
      
      mdata.back(). G0() = {TemperatureModel::T1, agam, nair,  0};
      mdata.back(). D0() = {TemperatureModel::T5, psf,  nair,  0};
      mdata.back ().Interp() = interp;
    }
    
    Model(const SingleSpeciesModel& bath, Type type) noexcept :
    mtype(type), mself(false), mbath(true), mspecies(1), mdata(1, bath) {}
    
    Model(Type type, bool self, bool bath, const ArrayOfSpeciesTag& species,
          const std::vector<SingleSpeciesModel>& ssms) :
    mtype(type), mself(self), mbath(bath), mspecies(species), mdata(ssms) {
      if(not OK())
        throw std::runtime_error("Bad initialization with different sizes or bad types, see documentation for valid initializations of LineShape::Model");
    }
    
    // Check if this model has the same species as another model
    bool same_broadening_species(const Model& other) const noexcept {
      if(mself not_eq other.mself)
        return false;
      else if(mbath not_eq other.mbath)
        return false;
      else if(mspecies.nelem() not_eq other.mspecies.nelem())
        return false;
      else {
        for(Index i=Index(mself); i<mspecies.nelem()-Index(mbath); i++)
          if(mspecies[i].Species() not_eq other.mspecies[i].Species())
            return false;
        return true;
      }
    }
    
    // The VMR vector required based on atmospheric species and self-existence
    Vector vmrs(const ConstVectorView& atmospheric_vmrs,
                const ArrayOfArrayOfSpeciesTag& atmospheric_species,
                const QuantumIdentifier& self) const {
      assert(atmospheric_species.nelem() == atmospheric_vmrs.nelem());
      
      // Initialize list of VMRS to 0
      Vector line_vmrs(mspecies.nelem(), 0);
      const Index back = mspecies.nelem()-1;  // Last index
      
      if(mtype == Type::DP)
        return line_vmrs;
      
      // Loop species ignoring self and bath
      for(Index i=0; i<mspecies.nelem(); i++) {
        if(mbath and i == back) {}
        else {
          // Select target in-case this is self-broadening
          const auto target = (mself and i == 0) ? 
                                  self.Species() :
                           mspecies[i].Species() ;
          
          // Find species in list or do nothing at all
          Index this_species_index = -1;
          for(Index j=0; j<atmospheric_species.nelem(); j++)
            if(atmospheric_species[j][0].Species() == target)
              this_species_index = j;
          
          // Set to non-zero in-case species exists
          if(this_species_index not_eq -1)
            line_vmrs[i] = atmospheric_vmrs[this_species_index];
        }
      }
      
      // Renormalize, if bath-species exist this is automatic.
      if(mbath)
        line_vmrs[back] = 1.0 - line_vmrs.sum();
      else
        line_vmrs /= line_vmrs.sum();
      
      // The result must be non-zero, a real number, and finite
      if(not std::isnormal(line_vmrs.sum()))
        throw std::runtime_error("Bad VMRs, your atmosphere does not support the line of interest");
      
      return line_vmrs;
    }
    
    Index this_species(const Index& self_species) const noexcept {
      if(mself)
        return 0;
      for(Index i=mself; i<nelem()-mbath; i++)
        if(self_species == mspecies[i].Species())
          return i;
      return -1;
    }
    
    Index this_species(const QuantumIdentifier& self) const noexcept {
      return this_species(self.Species());
    }
    
    // All main calculations
    #define LSPC(XVAR, PVAR)                                                                              \
    Numeric XVAR (Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Vector& vmrs) const noexcept { \
      return PVAR * std::inner_product(mdata.cbegin(), mdata.cend(), vmrs.begin(), 0.0,                   \
      std::plus<Numeric>(), [=](const SingleSpeciesModel& x, Numeric vmr) -> Numeric                      \
      {return vmr * x.compute(T, T0, Variable::XVAR);});                                                  \
    }
    LSPC(G0, P) LSPC(G2, P) LSPC(D0, P) LSPC(D2, P) LSPC(ETA, 1) LSPC(FVC, P) LSPC(Y, P) LSPC(G, P*P) LSPC(DV, P*P)
    #undef LSPC
    
    // All VMR derivatives
    #define LSPC(XVAR, PVAR)                                                                                               \
    Numeric d ## XVAR ## _dVMR (Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Index deriv_pos) const noexcept { \
      if(deriv_pos not_eq -1) return PVAR * mdata[deriv_pos].compute(T, T0, Variable::XVAR);                               \
      else return 0;                                                                                                       \
    }
    LSPC(G0, P) LSPC(G2, P) LSPC(D0, P) LSPC(D2, P) LSPC(ETA, 1) LSPC(FVC, P) LSPC(Y, P) LSPC(G, P*P) LSPC(DV, P*P)
    #undef LSPC
    
    // All shape model derivatives
    #define LSPC(XVAR, PVAR)                                                                                          \
    Numeric d ## XVAR ## _dT (Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Vector& vmrs) const noexcept { \
      return PVAR * std::inner_product(mdata.cbegin(), mdata.cend(), vmrs.begin(), 0.0,                               \
      std::plus<Numeric>(), [=](const SingleSpeciesModel& x, Numeric vmr) -> Numeric                                  \
      {return vmr * x.compute_dT (T, T0, Variable::XVAR);});                                                          \
    }
    LSPC(G0, P) LSPC(G2, P) LSPC(D0, P) LSPC(D2, P) LSPC(ETA, 1) LSPC(FVC, P) LSPC(Y, P) LSPC(G, P*P) LSPC(DV, P*P)
    #undef LSPC
    
    // All shape model derivatives
    #define LSPDC(XVAR, DERIV, PVAR)                                                               \
    Numeric d ## XVAR ## DERIV (Numeric T, Numeric T0, Numeric P [[maybe_unused]],                 \
                                Index deriv_pos, const Vector& vmrs) const noexcept {              \
      if(deriv_pos not_eq -1)                                                                      \
        return vmrs[deriv_pos] * PVAR * mdata[deriv_pos].compute ## DERIV (T, T0, Variable::XVAR); \
      else return 0;                                                                     \
    }
    LSPDC(G0,  _dT0,  P  ) LSPDC(G0,  _dX0, P  ) LSPDC(G0,  _dX1, P  ) LSPDC(G0,  _dX2, P  )
    LSPDC(G2,  _dT0,  P  ) LSPDC(G2,  _dX0, P  ) LSPDC(G2,  _dX1, P  ) LSPDC(G2,  _dX2, P  )
    LSPDC(D0,  _dT0,  P  ) LSPDC(D0,  _dX0, P  ) LSPDC(D0,  _dX1, P  ) LSPDC(D0,  _dX2, P  )
    LSPDC(D2,  _dT0,  P  ) LSPDC(D2,  _dX0, P  ) LSPDC(D2,  _dX1, P  ) LSPDC(D2,  _dX2, P  )
    LSPDC(ETA, _dT0,  1  ) LSPDC(ETA, _dX0, 1  ) LSPDC(ETA, _dX1, 1  ) LSPDC(ETA, _dX2, 1  )
    LSPDC(FVC, _dT0,  P  ) LSPDC(FVC, _dX0, P  ) LSPDC(FVC, _dX1, P  ) LSPDC(FVC, _dX2, P  )
    LSPDC(Y,   _dT0,  P  ) LSPDC(Y,   _dX0, P  ) LSPDC(Y,   _dX1, P  ) LSPDC(Y,   _dX2, P  )
    LSPDC(G,   _dT0,  P*P) LSPDC(G,   _dX0, P*P) LSPDC(G,   _dX1, P*P) LSPDC(G,   _dX2, P*P)
    LSPDC(DV,  _dT0,  P*P) LSPDC(DV,  _dX0, P*P) LSPDC(DV,  _dX1, P*P) LSPDC(DV,  _dX2, P*P)
    #undef LSPDC
    
    Output GetParams(Numeric T, Numeric T0, Numeric P, const Vector& vmrs) const {
      return {G0(T, T0, P, vmrs), D0(T, T0, P, vmrs),
              G2(T, T0, P, vmrs), D2(T, T0, P, vmrs),
              FVC(T, T0, P, vmrs), ETA(T, T0, P, vmrs),
              Y(T, T0, P, vmrs), G(T, T0, P, vmrs), DV(T, T0, P, vmrs)};
    }
    
    Output GetTemperatureDerivs(Numeric T, Numeric T0, Numeric P, const Vector& vmrs) const {
      return {dG0_dT(T, T0, P, vmrs), dD0_dT(T, T0, P, vmrs),
              dG2_dT(T, T0, P, vmrs), dD2_dT(T, T0, P, vmrs),
              dFVC_dT(T, T0, P, vmrs), dETA_dT(T, T0, P, vmrs),
              dY_dT(T, T0, P, vmrs), dG_dT(T, T0, P, vmrs), dDV_dT(T, T0, P, vmrs)};
    }
    
    Output GetVMRDerivs(Numeric T, Numeric T0, Numeric P, const Index pos) const noexcept {
      return {dG0_dVMR(T, T0, P, pos), dD0_dVMR(T, T0, P, pos),
              dG2_dVMR(T, T0, P, pos), dD2_dVMR(T, T0, P, pos),
              dFVC_dVMR(T, T0, P, pos), dETA_dVMR(T, T0, P, pos),
              dY_dVMR(T, T0, P, pos), dG_dVMR(T, T0, P, pos), dDV_dVMR(T, T0, P, pos)};
    }
    
    Numeric GetInternalDeriv(Numeric T, Numeric T0, Numeric P, Index pos, const Vector& vmrs, JacPropMatType deriv) const {
      if(pos < 0)
        return 0;
      
      #define RETURNINTERNALDERIVATIVE(TYPE)                                                         \
        case JacPropMatType::LineShape ## TYPE ## X0: return d ## TYPE ## _dX0(T, T0, P, pos, vmrs); \
        case JacPropMatType::LineShape ## TYPE ## X1: return d ## TYPE ## _dX1(T, T0, P, pos, vmrs); \
        case JacPropMatType::LineShape ## TYPE ## X2: return d ## TYPE ## _dX2(T, T0, P, pos, vmrs)
      switch(deriv) {
        RETURNINTERNALDERIVATIVE(G0);
        RETURNINTERNALDERIVATIVE(D0);
        RETURNINTERNALDERIVATIVE(G2);
        RETURNINTERNALDERIVATIVE(D2);
        RETURNINTERNALDERIVATIVE(FVC);
        RETURNINTERNALDERIVATIVE(ETA);
        RETURNINTERNALDERIVATIVE(Y);
        RETURNINTERNALDERIVATIVE(G);
        RETURNINTERNALDERIVATIVE(DV);
        default: return 0;
      }
      #undef RETURNINTERNALDERIVATIVE
    }
    
    // Size and size manipulation
    Index nelem() const {return mspecies.nelem();}
    void resize(Index n) {mspecies.resize(n); mdata.resize(n);}
    void reserve(Index n) {mspecies.reserve(n); mdata.reserve(n);}
    
    friend inline std::istream& operator>>(std::istream& is, Model& m);
    friend inline std::ostream& operator<<(std::ostream& os, const Model& m);
    friend std::istream& from_artscat4(std::istream& is, Model& lsc, const QuantumIdentifier& qid);
    friend std::istream& from_linefunctiondata(std::istream& data, Model& lsc);;
    
    Type ModelType() const noexcept {return mtype;}
    bool Self() const noexcept {return mself;}
    bool Bath() const noexcept {return mbath;}
    const std::vector<SingleSpeciesModel>& Data() const noexcept {return mdata;}
    const ArrayOfSpeciesTag& Species() const noexcept {return mspecies;}
    
    std::vector<SingleSpeciesModel>& Data() noexcept {return mdata;}
    
    void Set(const ModelParameters& param, const String& spec, const Variable var);
    
    ModelParameters Get(const String& spec, const Variable var) const {
      bool self = spec == self_broadening;
      bool bath = spec == bath_broadening;
      if(mself and self)
        return mdata.front().Get(var);
      else if(self)
        throw std::runtime_error("No self species but trying to get self in line shape model");
      else if(mbath and bath)
        return mdata.back().Get(var);
      else if(bath)
        throw std::runtime_error("No bath species but trying to get bath in line shape model");
      else {
        const SpeciesTag sp(spec);
        for(Index i=Index(mself); i<nelem()-Index(mbath); i++)
          if(sp.Species() == mspecies[i].Species())
            return mdata[i].Get(var);
        std::ostringstream os;
        os << "No species of type " << spec << " found in line shape model\n";
        os << "Available species are: " << mspecies << "\n";
        throw std::runtime_error(os.str());
      }
    }
    
    void Remove(Index i) {
      mspecies.erase(mspecies.begin()+i);
      mdata.erase(mdata.begin()+i);
      if(i == 0 and mself) mself=false;
      else if(i == nelem() and mbath) mbath=false;
    }
    
    void SetLineMixingModel(SingleSpeciesModel x) {
      for(auto& ssm: mdata) {
        ssm.Y() = x.Y();
        ssm.G() = x.G();
        ssm.DV() = x.DV();
        if(x.Y().type == TemperatureModel::LM_AER or
           x.G().type == TemperatureModel::LM_AER)
          ssm.Interp() = x.Interp();
      }
    }
  };

  std::istream& from_artscat4(std::istream& is, Model& lsc, const QuantumIdentifier& qid);
  std::istream& from_linefunctiondata(std::istream& data, Model& lsc);
  std::istream& from_linemixingdata(std::istream& data, Model& lsc);
  std::istream& from_pressurebroadeningdata(std::istream& data, Model& lsc, const QuantumIdentifier& qid);
  
  inline std::ostream& operator<<(std::ostream& os, const Model& m) {
    os << shapetype2string(m.mtype) << ' ' << m.nelem() 
                                    << ' ' << m.mself
                                    << ' ' << m.mbath;
    for(Index i=0; i<m.nelem(); i++) {
      if(m.mself and i==0) os << ' ' << self_broadening;
      else if(m.mbath and i==m.nelem()-1) os << ' ' << bath_broadening;
      else os << ' ' << m.mspecies[i].SpeciesNameMain();
    }
    for(auto& d: m.mdata) os << ' ' << d;
    return os;
  }
  
  inline std::istream& operator>>(std::istream& is, Model& m) {
    String tmp;
    Index nelem;
    is >> tmp >> nelem >> m.mself >> m.mbath;
    m.mtype = string2shapetype(tmp);
    m.mspecies.resize(nelem);
    m.mdata.resize(nelem);
    for(Index i=0; i<m.nelem(); i++) {
      is >> tmp;
      if(m.mself and i==0) {}
      else if(m.mbath and i==m.nelem()-1) {}
      else m.mspecies[i] = SpeciesTag(tmp);
    }
    for(auto& d: m.mdata) is >> d;
    return is;
  }
  
  // Legacy dealing with reading old LineFunctionData
  namespace LegacyLineFunctionData {
    /** Length per variable */
    inline Index temperaturemodel2legacynelem(TemperatureModel type) noexcept {
      switch(type) {
        case TemperatureModel::None:   return  0;
        case TemperatureModel::T0:     return  1;
        case TemperatureModel::T1:     return  2;
        case TemperatureModel::T2:     return  3;
        case TemperatureModel::T3:     return  2;
        case TemperatureModel::T4:     return  3;
        case TemperatureModel::T5:     return  2;
        case TemperatureModel::LM_AER: return 12;
      }
      std::terminate();  // Not allowed to reach, fix higher level code
    }
    
    /** Line shape models */
    inline std::vector<Variable> lineshapetag2variablesvector(String type) {
      if(type == String("DP"))        return {};
      else if(type == String("LP"))   return {Variable::G0, Variable::D0};
      else if(type == String("VP"))   return {Variable::G0, Variable::D0};
      else if(type == String("SDVP")) return {Variable::G0, Variable::D0, Variable::G2, Variable::D2};
      else if(type == String("HTP"))  return {Variable::G0, Variable::D0, Variable::G2, Variable::D2, Variable::FVC, Variable::ETA};
      else {
        std::ostringstream os;
        os << "Type: " << type << ", is not accepted.  "
        << "See documentation for accepted types\n";
        throw std::runtime_error(os.str());
      }
    }
    
    /** Line mixing models */
    inline std::vector<Variable> linemixingtag2variablesvector(String type) {
      if(type == "#")           return {};
      else if(type == "LM1")    return {Variable::Y};
      else if(type == "LM2")    return {Variable::Y, Variable::G, Variable::DV};
      else if(type == "INT")    return {};
      else if(type == "ConstG") return {Variable::G};
      else {
        std::ostringstream os;
        os << "Type: " << type << ", is not accepted.  "
        << "See documentation for accepted types\n";
        throw std::runtime_error(os.str());
      }
    }
  };
  
  namespace LegacyLineMixingData {
    enum class TypeLM {
      LM_NONE,                          // Reserved for no line mixing
      LM_LBLRTM,                        // Reserved for LBLRTM line mixing
      LM_LBLRTM_O2NonResonant,          // Reserved for the non-resonant O2 line in LBLRTM
      LM_1STORDER,                      // Reserved for Tretyakov et al. 2005 1st order of line mixing
      LM_2NDORDER,                      // Reserved for Makarov et al. 2011 second order of line mixing
      LM_BYBAND                         // Reserved for Paris data of relaxation matrix line mixing for band
    };
    
    inline LegacyLineMixingData::TypeLM string2typelm(String type) {
      if(type == "NA") // The standard case
        return TypeLM::LM_NONE;
      else if(type == "LL") // The LBLRTM case
        return TypeLM::LM_LBLRTM;
      else if(type == "NR") // The LBLRTM O2 non-resonant case
        return TypeLM::LM_LBLRTM_O2NonResonant;
      else if(type == "L2") // The 2nd order case
        return TypeLM::LM_2NDORDER;
      else if(type == "L1") // The 2nd order case
        return TypeLM::LM_1STORDER;
      else if(type == "BB") // The band class
        return TypeLM::LM_BYBAND;
      else {
        std::ostringstream os;
        os << "Type: " << type << ", is not accepted.  "
        << "See documentation for accepted types\n";
        throw std::runtime_error(os.str());
      }
    }
    
    inline Index typelm2nelem(LegacyLineMixingData::TypeLM type) {
      switch(type) {
        case TypeLM::LM_NONE: // The standard case
          return 0;
        case TypeLM::LM_LBLRTM: // The LBLRTM case
          return 12;
        case TypeLM::LM_LBLRTM_O2NonResonant: // Nonresonant is just a tag
          return 1;
        case TypeLM::LM_2NDORDER: // The 2nd order case
          return 10;
        case TypeLM::LM_1STORDER: // The 2nd order case
          return 3;
        case TypeLM::LM_BYBAND: // The band class
          return 1;
      }
      std::terminate();
    }
    
    Model vector2modellm(Vector x, LegacyLineMixingData::TypeLM type);
  };
  
  namespace LegacyPressureBroadeningData {
    enum class TypePB {
      PB_NONE,                          // No pressure broadening
      PB_AIR_BROADENING,                // Air broadening and self broadening only
      PB_AIR_AND_WATER_BROADENING,      // Air, water, and self broadening
      PB_PLANETARY_BROADENING,          // Gas broadening as done for solar system planets
      // PB_SD_AIR_VOLUME,                 // HTP in air for SD limit NOT SUPPORTED
      // PB_HTP_AIR_VOLUME,                // HTP in air NOT SUPPORTED
      // PB_VOIGT_TEST_WATER,              // Voigt parameters for testing NOT SUPPORTED
      // PB_SD_TEST_WATER,                 // SD parameters for testing NOT SUPPORTED
      // PB_PURELY_FOR_TESTING             // Testing tag for new input structures --- can be changed by anyone... NOT SUPPORTED
    };
    
    inline LegacyPressureBroadeningData::TypePB string2typepb(String type) {
      if(type == "NA") // The none case
        return TypePB::PB_NONE;
      else if(type == "N2") // Air Broadening is N2 broadening mostly...
        return TypePB::PB_AIR_BROADENING;
      else if(type == "WA") // Water and Air Broadening
        return TypePB::PB_AIR_AND_WATER_BROADENING;
      else if(type == "AP") // Planetary broadening
        return TypePB::PB_PLANETARY_BROADENING;
      else {
        std::ostringstream os;
        os << "Type: " << type << ", is not accepted.  "
        << "See documentation for accepted types\n";
        throw std::runtime_error(os.str());
      }
    }
    
    inline bool self_listed(const QuantumIdentifier& qid, LegacyPressureBroadeningData::TypePB t) {
      if(t == TypePB::PB_PLANETARY_BROADENING and
        (qid.Species() == SpeciesTag(String("N2")).Species() or
         qid.Species() == SpeciesTag(String("O2")).Species() or
         qid.Species() == SpeciesTag(String("H2O")).Species() or
         qid.Species() == SpeciesTag(String("CO2")).Species() or
         qid.Species() == SpeciesTag(String("H2")).Species() or
         qid.Species() == SpeciesTag(String("He")).Species()))
        return true;
      else if(t == TypePB::PB_AIR_AND_WATER_BROADENING and
              qid.Species() == SpeciesTag(String("H2O")).Species())
        return true;
      else
        return false;
    }
    
    inline Index typepb2nelem(LegacyPressureBroadeningData::TypePB type) {
      switch(type) {
        case TypePB::PB_NONE:
          return 0;
        case TypePB::PB_AIR_BROADENING:
          return 10;
        case TypePB::PB_AIR_AND_WATER_BROADENING:
          return 9;
        case TypePB::PB_PLANETARY_BROADENING:
          return 20;
      }
      std::terminate();
    }
    
    Model vector2modelpb(Vector x, LegacyPressureBroadeningData::TypePB type, bool self_in_list);
  };
};

#endif // linefunctiondata_h

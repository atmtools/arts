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
    G2=1,   // Pressure broadening speed-dependent
    D0=2,   // Pressure f-shifting speed-independent
    D2=3,   // Pressure f-shifting speed-dependent
    ETA=4,  // Correlation
    FVC=5,  // Frequency of velocity-changing collisions
    Y=6,    // First order line mixing coefficient
    G=7,    // Second order line mixing coefficient
    DV=8    // Second order line mixing f-shifting
    // ALWAYS ADD NEW AT THE END
  };
  
  struct ModelParameters {
    TemperatureModel type;
    Numeric X0;
    Numeric X1;
    Numeric X2;
    // ALWAYS ADD NEW AT THE END
  };
  
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
                       ModelParameters G2  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters D0  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters D2  = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters ETA = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters FVC = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters Y   = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters G   = {TemperatureModel::None, 0, 0, 0},
                       ModelParameters DV  = {TemperatureModel::None, 0, 0, 0}) :
    X({G0, G2, D0, D2, ETA, FVC, Y, G, DV}), V({0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0})
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
    }

    #undef x0
    #undef x1
    #undef x2
    
    // Access normal data
    #define ACCESS_INTERNAL(VARPOS) \
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
  
  enum class ShapeType {
    DP,   // Doppler
    LP,   // Lorentz
    VP,   // Voigt
    SDVP, // Speed-dependent Voigt
    HTP,  // Hartmann-Tran
  };
  
  inline String shapetype2string(ShapeType type) noexcept {
    switch(type) {
      case ShapeType::DP:   return "DP";
      case ShapeType::LP:   return "LP";
      case ShapeType::VP:   return "VP";
      case ShapeType::SDVP: return "SDVP";
      case ShapeType::HTP:  return "HTP";
    }
    std::terminate();  // Not allowed to reach, fix higher level code
  }
  
  inline ShapeType string2shapetype(const String& type) {
    if(type == "DP")                return ShapeType::DP;
    else if(type == String("LP"))   return ShapeType::LP;
    else if(type == String("VP"))   return ShapeType::VP;
    else if(type == String("SDVP")) return ShapeType::SDVP;
    else if(type == String("HTP"))  return ShapeType::HTP;
    else {
      std::ostringstream os;
      os << "Type: " << type << ", is not accepted.  "
      << "See documentation for accepted types\n";
      throw std::runtime_error(os.str());
    }
  }

  class Model {
  private:
    ShapeType mshapetype;
    bool mself;
    bool mbath;
    ArrayOfSpeciesTag mspecies;
    std::vector<SingleSpeciesModel> mdata;
    
  public:
    // No line shape parameters means this is a Doppler line
    Model() noexcept : mshapetype(ShapeType::DP), mself(false), mbath(false),
                       mspecies(0), mdata(0) {}
    
    // Standard HITRAN means this is a Voigt line
    Model(Numeric sgam, Numeric nself, Numeric agam, Numeric nair, Numeric psf) noexcept :
    mshapetype(ShapeType::VP), mself(true), mbath(true), mspecies(2), mdata(2) {
      mdata.front().G0() = {TemperatureModel::T1, sgam, nself, 0};
      mdata.front().D0() = {TemperatureModel::T5, psf,  nair,  0};
      mdata.back(). G0() = {TemperatureModel::T1, agam, nair,  0};
      mdata.back(). D0() = {TemperatureModel::T5, psf,  nair,  0};
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
        for(auto& vmr: line_vmrs) vmr /= line_vmrs.sum();
      
      // The result must be non-zero, a real number, and finite
      if(not std::isnormal(line_vmrs.sum()))
        throw std::runtime_error("Bad VMRS, do the species exists or have any VMR?");
      
      return line_vmrs;
    }
    
    // All main calculations
    #define LSPC(XVAR, PVAR) \
    Numeric XVAR (Numeric T, Numeric T0, Numeric P  [[maybe_unused]], const Vector& vmrs) const noexcept { \
      return PVAR * std::inner_product(mdata.cbegin(), mdata.cend(), vmrs.begin(), 0.0, \
                                       std::plus<Numeric>(), [=](const SingleSpeciesModel& x, Numeric vmr) -> Numeric \
                                       {return vmr * x.compute(T, T0, Variable::XVAR);}); \
    }
    LSPC(G0, P) LSPC(G2, P) LSPC(D0, P) LSPC(D2, P) LSPC(ETA, 1) LSPC(FVC, P) LSPC(Y, P) LSPC(G, P*P) LSPC(DV, P*P)
    #undef LSPC
    
    // All VMR derivatives
    #define LSPC(XVAR, PVAR) \
    Numeric XVAR (Numeric T, Numeric T0, Numeric P  [[maybe_unused]], const Index deriv_pos) const noexcept { \
      if(deriv_pos not_eq -1) return PVAR * mdata[deriv_pos].compute(T, T0, Variable::XVAR); \
      else return 0; \
    }
    LSPC(G0, P) LSPC(G2, P) LSPC(D0, P) LSPC(D2, P) LSPC(ETA, 1) LSPC(FVC, P) LSPC(Y, P) LSPC(G, P*P) LSPC(DV, P*P)
    #undef LSPC
    
    // All shape model derivatives
    #define LSPDC(XVAR, DERIV, PVAR) \
    Numeric d ## XVAR ## DERIV (Numeric T, Numeric T0, Numeric P  [[maybe_unused]], const Vector& vmrs) const noexcept { \
      return PVAR * std::inner_product(mdata.cbegin(), mdata.cend(), vmrs.begin(), 0.0, \
      std::plus<Numeric>(), [=](const SingleSpeciesModel& x, Numeric vmr) -> Numeric \
      {return vmr * x.compute ## DERIV (T, T0, Variable::XVAR);}); \
    }
    LSPDC(G0,  _dT,  P  ) LSPDC(G0,  _dT0,  P  ) LSPDC(G0,  _dX0, P  ) LSPDC(G0,  _dX1, P  ) LSPDC(G0,  _dX2, P  )
    LSPDC(G2,  _dT,  P  ) LSPDC(G2,  _dT0,  P  ) LSPDC(G2,  _dX0, P  ) LSPDC(G2,  _dX1, P  ) LSPDC(G2,  _dX2, P  )
    LSPDC(D0,  _dT,  P  ) LSPDC(D0,  _dT0,  P  ) LSPDC(D0,  _dX0, P  ) LSPDC(D0,  _dX1, P  ) LSPDC(D0,  _dX2, P  )
    LSPDC(D2,  _dT,  P  ) LSPDC(D2,  _dT0,  P  ) LSPDC(D2,  _dX0, P  ) LSPDC(D2,  _dX1, P  ) LSPDC(D2,  _dX2, P  )
    LSPDC(ETA, _dT,  1  ) LSPDC(ETA, _dT0,  1  ) LSPDC(ETA, _dX0, 1  ) LSPDC(ETA, _dX1, 1  ) LSPDC(ETA, _dX2, 1  )
    LSPDC(FVC, _dT,  P  ) LSPDC(FVC, _dT0,  P  ) LSPDC(FVC, _dX0, P  ) LSPDC(FVC, _dX1, P  ) LSPDC(FVC, _dX2, P  )
    LSPDC(Y,   _dT,  P  ) LSPDC(Y,   _dT0,  P  ) LSPDC(Y,   _dX0, P  ) LSPDC(Y,   _dX1, P  ) LSPDC(Y,   _dX2, P  )
    LSPDC(G,   _dT,  P*P) LSPDC(G,   _dT0,  P*P) LSPDC(G,   _dX0, P*P) LSPDC(G,   _dX1, P*P) LSPDC(G,   _dX2, P*P)
    LSPDC(DV,  _dT,  P*P) LSPDC(DV,  _dT0,  P*P) LSPDC(DV,  _dX0, P*P) LSPDC(DV,  _dX1, P*P) LSPDC(DV,  _dX2, P*P)
    #undef LSPDC
    
    // Size and size manipulation
    Index nelem() const {return mspecies.nelem();}
    void resize(Index n) {mspecies.resize(n); mdata.resize(n);}
    void reserve(Index n) {mspecies.reserve(n); mdata.reserve(n);}
    
    friend inline std::istream& operator>>(std::istream& is, Model& m);
    friend inline std::ostream& operator<<(std::ostream& os, const Model& m);
    friend std::istream& from_artscat4(std::istream& is, Model& lsc, bool self_in_list);
    friend std::istream& from_linefunctiondata(std::istream& data, Model& lsc);;
  };

  std::istream& from_artscat4(std::istream& is, Model& lsc, bool self_in_list);
  std::istream& from_linefunctiondata(std::istream& data, Model& lsc);
  
  inline std::ostream& operator<<(std::ostream& os, const Model& m) {
    os << shapetype2string(m.mshapetype) << ' ' << m.nelem() << ' ' << m.mself << ' ' << m.mbath;
    for(auto& s: m.mspecies) os << ' ' << s.SpeciesNameMain();
    for(auto& d: m.mdata) os << ' ' << d;
    return os;
  }
  
  inline std::istream& operator>>(std::istream& is, Model& m) {
    String tmp;
    Index nelem;
    is >> tmp >> nelem >> m.mself >> m.mbath;
    m.mshapetype = string2shapetype(tmp);
    m.mspecies.resize(nelem);
    m.mdata.resize(nelem);
    for(auto& s: m.mspecies) {is >> tmp; s = SpeciesTag(tmp);}
    for(auto& d: m.mdata) is >> d;
    return is;
  }
  
  // Legacy dealing with reading old LineFunctionData
  namespace LegacyLineFunctionData {
    // Length per variable
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
};

#endif // linefunctiondata_h

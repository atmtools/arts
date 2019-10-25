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

/** Contains the line shape namespace
 * @file   lineshapemodel.h
 * @author Richard Larsson
 * @date   2018-09-19
 * 
 * @brief  Contains the line shape namespace
 * 
 * This namespace computes all line shape parameters
 * for any set of line shape we can use in ARTS.  Should
 * be extended for more use as seen fit.
 **/

#ifndef lineshapemodel_h
#define lineshapemodel_h

#include <numeric>
#include <algorithm>
#include "abs_species_tags.h"
#include "constants.h"
#include "jacobian.h"

/** Return the deriveative type based on string input 
 * 
 * @param[in] var Variable in AllLineShapeVars()
 * @param[in] coeff Coefficient in AllLineShapeCoeffs()
 * 
 * @return Derivative
 */
JacPropMatType select_derivativeLineShape(const String& var,
                                          const String& coeff);

/** All available line shape coefficients */
ArrayOfString AllLineShapeCoeffs();

/** All available line shape variables */
ArrayOfString AllLineShapeVars();

/** Computations of line shape derived parameters
 * 
 * Defines many classes and IO routines for line 
 * shape parameters to comply with everything 
 * from no line mixing Doppler to coefficient-based
 * line mixing Hartman-Tran profiles
 */
namespace LineShape {

/** Temperature models
 * 
 * Each input here should correspond to a
 * method of how to compute the variable
 * given the coefficients and Interpolation
 * data available to SingleSpeciesModel
 */
enum class TemperatureModel : Index {
  None,   // 0
  T0,     // Constant, X0
  T1,     // Standard, X0 * (T0/T) ^ X1
  T2,     // X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0));
  T3,     // X0 + X1 * (T - T0)
  T4,     // (X0 + X1 * (T0/T - 1)) * (T0/T)^X2;
  T5,     // X0 * (T0/T)^(0.25 + 1.5*X1)
  LM_AER, // X(200) = X0; X(250) = X1; X(298) = X2; X(340) = X3;  Linear interpolation in between
  DPL     // X0 * (T0/T) ^ X1 + X2 * (T0/T) ^ X3
  // ALWAYS ADD NEW AT THE END
};

/** Turns selected TemperatureModel type into a string
 * 
 * This function takes the input TemperatureModel
 * and returns it as a string
 *  
 * @param[in] type The temperature model type
 * 
 * @return String of temperature model type
 */
inline String temperaturemodel2string(TemperatureModel type) noexcept {
  switch (type) {
    case TemperatureModel::None:
      return "#";
    case TemperatureModel::T0:
      return "T0";
    case TemperatureModel::T1:
      return "T1";
    case TemperatureModel::T2:
      return "T2";
    case TemperatureModel::T3:
      return "T3";
    case TemperatureModel::T4:
      return "T4";
    case TemperatureModel::T5:
      return "T5";
    case TemperatureModel::LM_AER:
      return "LM_AER";
    case TemperatureModel::DPL:
      return "DPL";
  }
  std::terminate();  // Not allowed to reach, fix higher level code
}

/** Turns predefined strings into a TemperatureModel type
 * 
 * This function either acts as the inverse of temperaturemodel2string
 * or it throws a runtime error
 * 
 * @param[in] type The TemperatureModel string
 * 
 * @return The actual TemperatureModel type
 */
inline TemperatureModel string2temperaturemodel(const String& type) {
  if (type == "#")
    return TemperatureModel::None;
  else if (type == String("T0"))
    return TemperatureModel::T0;
  else if (type == String("T1"))
    return TemperatureModel::T1;
  else if (type == String("T2"))
    return TemperatureModel::T2;
  else if (type == String("T3"))
    return TemperatureModel::T3;
  else if (type == String("T4"))
    return TemperatureModel::T4;
  else if (type == String("T5"))
    return TemperatureModel::T5;
  else if (type == String("LM_AER"))
    return TemperatureModel::LM_AER;
  else if (type == String("DPL"))
    return TemperatureModel::DPL;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}

/** List of possible shape variables
 * 
 * Should correspond to strings in AllLineShapeVars()
 */
enum class Variable {
  G0 = 0,   // Pressure broadening speed-independent
  D0 = 1,   // Pressure f-shifting speed-dependent
  G2 = 2,   // Pressure broadening speed-dependent
  D2 = 3,   // Pressure f-shifting speed-independent
  FVC = 4,  // Frequency of velocity-changing collisions
  ETA = 5,  // Correlation
  Y = 6,    // First order line mixing coefficient
  G = 7,    // Second order line mixing coefficient
  DV = 8    // Second order line mixing f-shifting
  // ALWAYS ADD NEW AT THE END
};

/** Output operator for Variable to be human-readable */
inline std::ostream& operator<<(std::ostream& os, Variable v) {
  switch (v) {
    case Variable::G0:
      os << "G0";
      break;
    case Variable::D0:
      os << "D0";
      break;
    case Variable::G2:
      os << "G2";
      break;
    case Variable::D2:
      os << "D2";
      break;
    case Variable::FVC:
      os << "FVC";
      break;
    case Variable::ETA:
      os << "ETA";
      break;
    case Variable::Y:
      os << "Y";
      break;
    case Variable::G:
      os << "G";
      break;
    case Variable::DV:
      os << "DV";
      break;
  }
  return os;
}

/** Turns selected Variable type into a string
 * 
 * This function takes the input Variable
 * and returns it as a string
 *  
 * @param[in] type The Variable type
 * 
 * @return String of Variable type
 */
inline String variable2string(Variable type) noexcept {
#define VARIABLE2STRINGDEF(X) \
  case Variable::X:           \
    return #X
  switch (type) {
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

/** Turns predefined strings into a Variable type
 * 
 * This function either acts as the inverse of variable2string
 * or it throws a runtime error
 * 
 * @param[in] type The Variable string
 * 
 * @return The actual Variable type
 */
inline Variable string2variable(const String& type) {
#define STRING2VARIABLEDEF(X) \
  if (type == #X) return Variable::X
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

/** Coefficients and temperature model for SingleSpeciesModel 
 * 
 * NOTE: Developer should always add new coefficients at the end
 */
struct ModelParameters {
  TemperatureModel type;
  Numeric X0;
  Numeric X1;
  Numeric X2;
  Numeric X3;
};

String modelparameters2metadata(const ModelParameters mp, const Numeric T0);

/** Get a coefficient from ModelParameters by name
 * 
 * Will throw a runtime_error if type is bad
 * 
 * @param[in] mp The model parameters
 * @param[in] type The coefficient by name
 * 
 * @return a reference to the coefficient
 */
inline Numeric& SingleModelParameter(ModelParameters& mp, const String& type) {
  if (type == "X0")
    return mp.X0;
  else if (type == "X1")
    return mp.X1;
  else if (type == "X2")
    return mp.X2;
  else if (type == "X3")
    return mp.X3;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
  std::terminate();
}

/** Output operator for ModelParameters */
inline std::ostream& operator<<(std::ostream& os, const ModelParameters& mp) {
  os << temperaturemodel2string(mp.type) << ' ' << mp.X0 << ' ' << mp.X1 << ' '
  << mp.X2 << ' ' << mp.X3 << ' ';
  return os;
}

/** Input operator for ModelParameters */
inline std::istream& operator>>(std::istream& is, ModelParameters& mp) {
  String tmp;
  is >> tmp >> mp.X0 >> mp.X1 >> mp.X2 >> mp.X3;
  mp.type = string2temperaturemodel(tmp);
  return is;
}

/** Current max number of coefficients */
constexpr Index nmaxTempModelParams = 4;

/** Current max number of line shape variables */
constexpr Index nVars = 9;

/** Compute the line shape parameters for a single broadening species */
class SingleSpeciesModel {
 private:
  std::array<ModelParameters, nVars> X;

  /** Line mixing as done by AER data in ARTS
   * 
   * Uses piece-wise linear interpolation and extrapolates at the edges
   * 
   * var must be G or Y
   * 
   * @param[in] T The temperature
   * @param[in] var The variable
   * 
   * @return The broadening parameter at temperature
   */
  Numeric special_linemixing_aer(Numeric T, ModelParameters mp) const noexcept {
    if (T < 250)
      return mp.X0 + (T - 200) * (mp.X1 - mp.X0) / (250 - 200);
    else if (T > 296)
      return mp.X2 + (T - 296) * (mp.X3 - mp.X2) / (340 - 296);
    else
      return mp.X1 + (T - 250) * (mp.X2 - mp.X1) / (296 - 250);
  }
  
  /** The temperature derivative of special_linemixing_aer
   * 
   * @param[in] T The temperature
   * @param[in] var The variable
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  Numeric special_linemixing_aer_dT(Numeric T, ModelParameters mp) const noexcept {
    if (T < 250)
      return (mp.X1 - mp.X0) / (250 - 200);
    else if (T > 296)
      return (mp.X3 - mp.X2) / (340 - 296);
    else
      return (mp.X2 - mp.X1) / (296 - 250);
  }

 public:
  /** Default initialization */
  constexpr SingleSpeciesModel(
      ModelParameters G0 = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters D0 = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters G2 = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters D2 = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters FVC = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters ETA = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters Y = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters G = {TemperatureModel::None, 0, 0, 0, 0},
      ModelParameters DV = {TemperatureModel::None, 0, 0, 0, 0})
      : X({G0, D0, G2, D2, FVC, ETA, Y, G, DV}) {}

#define x0 X[Index(var)].X0
#define x1 X[Index(var)].X1
#define x2 X[Index(var)].X2
#define x3 X[Index(var)].X3

/** Compute the broadening parameter at the input
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return The broadening parameter at temperature
 */
Numeric compute(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return x0;
    case TemperatureModel::T1:
      return x0 * pow(T0 / T, x1);
    case TemperatureModel::T2:
      return x0 * pow(T0 / T, x1) * (1 + x2 * log(T / T0));
    case TemperatureModel::T3:
      return x0 + x1 * (T - T0);
    case TemperatureModel::T4:
      return (x0 + x1 * (T0 / T - 1.)) * pow(T0 / T, x2);
    case TemperatureModel::T5:
      return x0 * pow(T0 / T, 0.25 + 1.5 * x1);
    case TemperatureModel::LM_AER:
      return special_linemixing_aer(T, X[Index(var)]);
    case TemperatureModel::DPL:
      return x0 * pow(T0 / T, x1) + x2 * pow(T0 / T, x3);
  }
  std::terminate();
}

/** Derivative of compute(...) wrt x0
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return Derivative of compute(...) wrt x0
 */
Numeric compute_dX0(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 1;
    case TemperatureModel::T1:
      return pow(T0 / T, x1);
    case TemperatureModel::T2:
      return pow(T0 / T, x1) * (1 + x2 * log(T / T0));
    case TemperatureModel::T3:
      return 1;
    case TemperatureModel::T4:
      return pow(T0 / T, x2);
    case TemperatureModel::T5:
      return pow(T0 / T, 1.5 * x1 + 0.25);
    case TemperatureModel::LM_AER:
      return 0;
    case TemperatureModel::DPL:
      return pow(T0 / T, x1);
  }
  std::terminate();
}

/** Derivative of compute(...) wrt x1
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return Derivative of compute(...) wrt x1
 */
Numeric compute_dX1(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 0;
    case TemperatureModel::T1:
      return x0 * pow(T0 / T, x1) * log(T0 / T);
    case TemperatureModel::T2:
      return x0 * pow(T0 / T, x1) * (x2 * log(T / T0) + 1.) * log(T0 / T);
    case TemperatureModel::T3:
      return (T - T0);
    case TemperatureModel::T4:
      return pow(T0 / T, x2) * (T0 / T - 1.);
    case TemperatureModel::T5:
      return 1.5 * x0 * pow(T0 / T, 1.5 * x1 + 0.25) * log(T0 / T);
    case TemperatureModel::LM_AER:
      return 0;
    case TemperatureModel::DPL:
      return x0 * pow(T0 / T, x1) * log(T0 / T);
  }
  std::terminate();
}

/** Derivative of compute(...) wrt x2
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return Derivative of compute(...) wrt x2
 */
Numeric compute_dX2(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 0;
    case TemperatureModel::T1:
      return 0;
    case TemperatureModel::T2:
      return x0 * pow(T0 / T, x1) * log(T / T0);
    case TemperatureModel::T3:
      return 0;
    case TemperatureModel::T4:
      return pow(T0 / T, x2) * (x0 + x1 * (T0 / T - 1)) * log(T0 / T);
    case TemperatureModel::T5:
      return 0;
    case TemperatureModel::LM_AER:
      return 0;
    case TemperatureModel::DPL:
      return pow(T0 / T, x3);
  }
  std::terminate();
}

/** Derivative of compute(...) wrt x3
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return Derivative of compute(...) wrt x3
 */
Numeric compute_dX3(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 0;
    case TemperatureModel::T1:
      return 0;
    case TemperatureModel::T2:
      return 0;
    case TemperatureModel::T3:
      return 0;
    case TemperatureModel::T4:
      return 0;
    case TemperatureModel::T5:
      return 0;
    case TemperatureModel::LM_AER:
      return 0;
    case TemperatureModel::DPL:
      return x2 * pow(T0 / T, x3) * log(T0 / T);
  }
  std::terminate();
}

/** Derivative of compute(...) wrt T
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return Derivative of compute(...) wrt T
 */
Numeric compute_dT(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 0;
    case TemperatureModel::T1:
      return -x0 * x1 * pow(T0 / T, x1) / T;
    case TemperatureModel::T2:
      return -x0 * x1 * pow(T0 / T, x1) * (x2 * log(T / T0) + 1.) / T +
      x0 * x2 * pow(T0 / T, x1) / T;
    case TemperatureModel::T3:
      return x1;
    case TemperatureModel::T4:
      return -x2 * pow(T0 / T, x2) * (x0 + x1 * (T0 / T - 1.)) / T -
      T0 * x1 * pow(T0 / T, x2) / pow(T, 2);
    case TemperatureModel::T5:
      return -x0 * pow(T0 / T, 1.5 * x1 + 0.25) * (1.5 * x1 + 0.25) / T;
    case TemperatureModel::LM_AER:
      return special_linemixing_aer_dT(T, X[Index(var)]);
    case TemperatureModel::DPL:
      return -x0 * x1 * pow(T0 / T, x1) / T + -x2 * x3 * pow(T0 / T, x3) / T;
  }
  std::terminate();
}

/** Derivative of compute(...) wrt T0
 * 
 * @param[in] T The temperature
 * @param[in] T0 The temperature used to derive the coefficients
 * @param[in] var The variable
 * 
 * @return Derivative of compute(...) wrt T0
 */
Numeric compute_dT0(Numeric T, Numeric T0, Variable var) const noexcept {
  using std::log;
  using std::pow;
  
  switch (X[Index(var)].type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 0;
    case TemperatureModel::T1:
      return x0 * x1 * pow(T0 / T, x1) / T0;
    case TemperatureModel::T2:
      return x0 * x1 * pow(T0 / T, x1) * (x2 * log(T / T0) + 1.) / T0 -
      x0 * x2 * pow(T0 / T, x1) / T0;
    case TemperatureModel::T3:
      return -x1;
    case TemperatureModel::T4:
      return x2 * pow(T0 / T, x2) * (x0 + x1 * (T0 / T - 1.)) / T0 +
      x1 * pow(T0 / T, x2) / T;
    case TemperatureModel::T5:
      return x0 * pow(T0 / T, 1.5 * x1 + 0.25) * (1.5 * x1 + 0.25) / T0;
    case TemperatureModel::LM_AER:
      return 0;
    case TemperatureModel::DPL:
      return x0 * x1 * pow(T0 / T, x1) / T0 + x2 * x3 * pow(T0 / T, x3) / T0;
  }
  std::terminate();
}

#undef x0
#undef x1
#undef x2
#undef x3

#define ACCESS_INTERNAL(VARPOS)                                             \
  ModelParameters& VARPOS() noexcept { return X[Index(Variable::VARPOS)]; } \
  constexpr ModelParameters VARPOS() const noexcept { return X[Index(Variable::VARPOS)]; }
  ACCESS_INTERNAL(G0);
  ACCESS_INTERNAL(D0);
  ACCESS_INTERNAL(G2);
  ACCESS_INTERNAL(D2);
  ACCESS_INTERNAL(FVC);
  ACCESS_INTERNAL(ETA);
  ACCESS_INTERNAL(Y);
  ACCESS_INTERNAL(G);
  ACCESS_INTERNAL(DV);
#undef ACCESS_INTERNAL

  /** Get internal Data reference */
  std::array<ModelParameters, nVars>& Data() noexcept { return X; }
  
  /** Get const internal Data reference */
  const std::array<ModelParameters, nVars>& Data() const noexcept { return X; }

  /** Set variable to a different ModelParameters
   * 
   * @param[in] var The variable
   * @param[in] x The new ModelParameters for var
   */
  void Set(Variable var, const ModelParameters& x) noexcept {
#define MODELPARAMCASESETTER(X) \
  case Variable::X:             \
    X() = x;                    \
    break
    switch (var) {
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

  /** Get variable by type
   * 
   * @param[in] var The variable type
   * 
   * @return ModelParameters copy
   */
  ModelParameters Get(Variable var) const noexcept {
#define MODELPARAMCASEGETTER(X) \
  case Variable::X:             \
    return X();
    switch (var) {
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

  /** Binary read for SingleSpeciesModel */
  bifstream& read(bifstream& bif) {
    for (auto& data: X) {
      Index x;
      bif >> x >> data.X0 >> data.X1 >> data.X2 >> data.X3;
      data.type = TemperatureModel(x);
    }
    return bif;
  }

  /** Binary write for SingleSpeciesModel */
  bofstream& write(bofstream& bof) const {
    for (auto& data: X) {
      bof << Index(data.type) << data.X0 << data.X1 << data.X2 << data.X3;
    }
    return bof;
  }
  
  bool MatchTypes(const SingleSpeciesModel& other) const noexcept {
    return std::equal (X.cbegin(), X.cend(), other.X.cbegin(), other.X.cend(), [](auto& a, auto& b){return a.type == b.type;});
  }
};

/** Output operator for SingleSpeciesModel */
inline std::ostream& operator<<(std::ostream& os,
                                const SingleSpeciesModel& ssm) {
  for (const auto& mp : ssm.Data())
    if (mp.type not_eq TemperatureModel::None)
      os << mp.X0 << ' ' << mp.X1 << ' ' << mp.X2 << ' ' << mp.X3 << ' ';
  return os;
}

/** Input operator for SingleSpeciesModel */
inline std::istream& operator>>(std::istream& is, SingleSpeciesModel& ssm) {
  for (auto& mp : ssm.Data())
    if(mp.type not_eq TemperatureModel::None)
      is >> mp.X0 >> mp.X1 >> mp.X2 >> mp.X3;
  return is;
}

/** Type of line shape to compute */
enum class Type {
  DP,    // Doppler
  LP,    // Lorentz
  VP,    // Voigt
  SDVP,  // Speed-dependent Voigt
  HTP,   // Hartmann-Tran
};

/** Turns selected Type into a string
 * 
 * This function takes the input Type
 * and returns it as a string
 *  
 * @param[in] type The Type
 * 
 * @return String of Type
 */
inline String shapetype2string(Type type) noexcept {
  switch (type) {
    case Type::DP:
      return "DP";
    case Type::LP:
      return "LP";
    case Type::VP:
      return "VP";
    case Type::SDVP:
      return "SDVP";
    case Type::HTP:
      return "HTP";
  }
  std::terminate();  // Not allowed to reach, fix higher level code
}

/** Turns selected Type into a human readable string
 * 
 * This function takes the input Type
 * and returns it as a string
 *  
 * @param[in] type The Type
 * 
 * @return String of Type
 */
inline String shapetype2metadatastring(Type type) noexcept {
  switch (type) {
    case Type::DP:
      return "The line shape type is the Doppler profile\n";
    case Type::LP:
      return "The line shape type is the Lorentz profile.\n";
    case Type::VP:
      return "The line shape type is the Voigt profile.\n";
    case Type::SDVP:
      return "The line shape type is the speed-dependent Voigt profile.\n";
    case Type::HTP:
      return "The line shape type is the Hartmann-Tran profile.\n";
  }
  std::terminate();  // Not allowed to reach, fix higher level code
}

/** Turns predefined strings into a Type
 * 
 * This function either acts as the inverse of shapetype2string
 * or it throws a runtime error
 * 
 * @param[in] type The Type string
 * 
 * @return The actual Type
 */
inline Type string2shapetype(const String& type) {
  if (type == "DP")
    return Type::DP;
  else if (type == String("LP"))
    return Type::LP;
  else if (type == String("VP"))
    return Type::VP;
  else if (type == String("SDVP"))
    return Type::SDVP;
  else if (type == String("HTP"))
    return Type::HTP;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}

/** Main output of Model */
struct Output {
  Numeric G0,  // Pressure broadening speed-independent
      D0,      // Pressure f-shifting speed-independent
      G2,      // Pressure broadening speed-dependent
      D2,      // Pressure f-shifting speed-dependent
      FVC,     // Frequency of velocity-changing collisions
      ETA,     // Correlation
      Y,       // First order line mixing coefficient
      G,       // Second order line mixing coefficient
      DV;      // Second order line mixing f-shifting
};

/** Output operator for LineShape::Output */
inline std::ostream& operator<<(std::ostream& os, Output x) {
  return os << "G0: " << x.G0 << " D0: " << x.D0 << " G2: " << x.G2
            << " D2: " << x.D2 << " FVC: " << x.FVC << " ETA: " << x.ETA
            << " Y: " << x.Y << " G: " << x.G << " DV: " << x.DV;
}

/** Output to be used by mirroring calls */
constexpr Output mirroredOutput(Output x) noexcept {
  return {x.G0, -x.D0, x.G2, -x.D2, x.FVC, x.ETA, x.Y, x.G, -x.DV};
}

/** Output turned negative */
constexpr Output negativeOutput(Output x) noexcept {
  return {-x.G0, -x.D0, -x.G2, -x.D2, -x.FVC, -x.ETA, -x.Y, -x.G, -x.DV};
}

/** Output turned from SI to CGS units */
constexpr Output si2cgs(Output x) noexcept {
  using Conversion::freq2kaycm;
  return {freq2kaycm(x.G0),
          freq2kaycm(x.D0),
          freq2kaycm(x.D2),
          freq2kaycm(x.G2),
          freq2kaycm(x.FVC),
          x.ETA,
          x.Y,
          x.G,
          freq2kaycm(x.DV)};
}

/** Diff of two output */
constexpr Output differenceOutput(Output y, Output x) noexcept {
  return {y.G0 - x.G0,
          y.D0 - x.D0,
          y.G2 - x.G2,
          y.D2 - x.D2,
          y.FVC - x.FVC,
          y.ETA - x.ETA,
          y.Y - x.Y,
          y.G - x.G,
          y.DV - x.DV};
}



/** Returns a VMR vector for this model's main calculations
 * 
 * Sets a vector that matches the mdata size of VMRs based
 * on atmospheric species and VMRs
 * 
 * Only checks the first species in inner atmosphere
 * 
 * Renormalizes the values to unity.  If this renormalization
 * is impossible then it throws an error
 * 
 * Returns 0s if type is Doppler line shape
 * 
 * @param[in] atmospheric_vmrs VMRS in atmosphere
 * @param[in] atmospheric_species Species in atmosphere
 * @param[in] self An ID of whichever species is self
 * @param[in] lineshape_species Species affecting lineshape
 * @param[in] self_in_list Affects lineshape by itself?
 * @param[in] bath_in_list Affected lineshape by environment?
 * @param[in] type The type of line shape
 */
Vector vmrs(const ConstVectorView& atmospheric_vmrs,
            const ArrayOfArrayOfSpeciesTag& atmospheric_species,
            const QuantumIdentifier& self,
            const ArrayOfSpeciesTag& lineshape_species,
            bool self_in_list,
            bool bath_in_list,
            Type type);

/** Name for bath broadening in printing and reading user input */
static constexpr const char* const bath_broadening = "AIR";

/** Name for self broadening in printing and reading user input */
static constexpr const char* const self_broadening = "SELF";


/** Main line shape model class
 * 
 * Computes the line shape parameters for a given atmosphere
 */
class Model {
 private:
   std::vector<SingleSpeciesModel> mdata;

 public:
  /** Default init just sets the size */
  Model(Index n=0) noexcept : mdata(n) {}
  
  /** Init from copying a vector */
  Model(const std::vector<SingleSpeciesModel>& assm) noexcept : mdata(assm) {}
  
  /** Init from copying itself */
  Model(const Model& m) noexcept : Model(m.mdata) {}
  
  /** Init from moving a vector */
  Model(std::vector<SingleSpeciesModel>&& assm) noexcept : mdata(std::move(assm)) {}
  
  /** Init from moving a itself */
  Model(Model&& m) noexcept : Model(std::move(m.mdata)) {}
  
  /** Copy and equals */
  Model& operator=(const Model& m) {mdata=m.mdata; return *this;}
  
  /** Move and equals */
  Model& operator=(Model&& m) {mdata=std::move(m.mdata); return *this;}
  
  /** Standard HITRAN init
   * 
   * @param[in] sgam Self pressure broadening coefficient
   * @param[in] nself Self pressure broadening exponent
   * @param[in] agam Air pressure broadening coefficient
   * @param[in] nair Air pressure broadening exponent
   * @param[in] psf All pressure shifting coefficient
   * @param[in] interp The interpolation variable for AER type line mixing
   */
  Model(Numeric sgam,
        Numeric nself,
        Numeric agam,
        Numeric nair,
        Numeric psf,
        std::array<Numeric, 12> aer_interp =
             {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) noexcept : mdata(2) {
    mdata.front().G0() = {TemperatureModel::T1, sgam, nself, 0, 0};
    mdata.front().D0() = {TemperatureModel::T5, psf, nair, 0, 0};

    mdata.back().G0() = {TemperatureModel::T1, agam, nair, 0, 0};
    mdata.back().D0() = {TemperatureModel::T5, psf, nair, 0, 0};
    
    if (std::any_of(aer_interp.cbegin(), aer_interp.cend(), [](auto x){return x not_eq 0;})) {
      mdata.front().Y().type = TemperatureModel::LM_AER;
      mdata.front().Y().X0 = aer_interp[4];
      mdata.front().Y().X1 = aer_interp[5];
      mdata.front().Y().X2 = aer_interp[6];
      mdata.front().Y().X3 = aer_interp[7];
      mdata.front().G().type = TemperatureModel::LM_AER;
      mdata.front().G().X0 = aer_interp[8];
      mdata.front().G().X1 = aer_interp[9];
      mdata.front().G().X2 = aer_interp[10];
      mdata.front().G().X3 = aer_interp[11];
      
      mdata.back().Y().type = TemperatureModel::LM_AER;
      mdata.back().Y().X0 = aer_interp[4];
      mdata.back().Y().X1 = aer_interp[5];
      mdata.back().Y().X2 = aer_interp[6];
      mdata.back().Y().X3 = aer_interp[7];
      mdata.back().G().type = TemperatureModel::LM_AER;
      mdata.back().G().X0 = aer_interp[8];
      mdata.back().G().X1 = aer_interp[9];
      mdata.back().G().X2 = aer_interp[10];
      mdata.back().G().X3 = aer_interp[11];
    }
  }
  
  /** The Model is good to use
   * 
   * @return true/false
   */
  bool OK(Type type, bool self, bool bath, const std::vector<SpeciesTag>& species) const noexcept {
    Index n = mdata.size();
    Index k = species.size();
    Index m = Index(self) + Index(bath);
    bool needs_any = type not_eq Type::DP;
    if (n not_eq k or m > n or (needs_any and not n))
      return false;
    else
      return true;
  }

#define LSPC(XVAR, PVAR)                                                     \
  Numeric XVAR(                                                              \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Vector& vmrs) \
      const noexcept {                                                       \
    return PVAR *                                                            \
           std::inner_product(                                               \
               mdata.cbegin(),                                               \
               mdata.cend(),                                                 \
               vmrs.begin(),                                                 \
               0.0,                                                          \
               std::plus<Numeric>(),                                         \
               [=](auto& x, auto vmr) -> Numeric {    \
                 return vmr * x.compute(T, T0, Variable::XVAR);              \
               });                                                           \
  }
  LSPC(G0, P)
  LSPC(D0, P)
  LSPC(G2, P)
  LSPC(D2, P)
  LSPC(FVC, P)
  LSPC(ETA, 1)
  LSPC(Y, P)
  LSPC(G, P* P)
  LSPC(DV, P* P)
#undef LSPC

#define LSPCV(XVAR, PVAR)                                            \
  Numeric d##XVAR##_dVMR(Numeric T,                                  \
                         Numeric T0,                                 \
                         Numeric P [[maybe_unused]],                 \
                         const Index deriv_pos) const noexcept {     \
    if (deriv_pos not_eq -1)                                         \
      return PVAR * mdata[deriv_pos].compute(T, T0, Variable::XVAR); \
    else                                                             \
      return 0;                                                      \
  }
  LSPCV(G0, P)
  LSPCV(D0, P)
  LSPCV(G2, P)
  LSPCV(D2, P)
  LSPCV(FVC, P)
  LSPCV(ETA, 1)
  LSPCV(Y, P)
  LSPCV(G, P* P)
  LSPCV(DV, P* P)
#undef LSPCV

#define LSPCT(XVAR, PVAR)                                                    \
  Numeric d##XVAR##_dT(                                                      \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], const Vector& vmrs) \
      const noexcept {                                                       \
    return PVAR *                                                            \
           std::inner_product(                                               \
               mdata.cbegin(),                                               \
               mdata.cend(),                                                 \
               vmrs.begin(),                                                 \
               0.0,                                                          \
               std::plus<Numeric>(),                                         \
               [=](auto& x, auto vmr) -> Numeric {                           \
                 return vmr * x.compute_dT(T, T0, Variable::XVAR);           \
               });                                                           \
  }
  LSPCT(G0, P)
  LSPCT(D0, P)
  LSPCT(G2, P)
  LSPCT(D2, P)
  LSPCT(FVC, P)
  LSPCT(ETA, 1)
  LSPCT(Y, P)
  LSPCT(G, P* P)
  LSPCT(DV, P* P)
#undef LSPCT

// All shape model derivatives
#define LSPDC(XVAR, DERIV, PVAR)                                     \
  Numeric d##XVAR##DERIV(Numeric T,                                  \
                         Numeric T0,                                 \
                         Numeric P [[maybe_unused]],                 \
                         Index deriv_pos,                            \
                         const Vector& vmrs) const noexcept {        \
    if (deriv_pos not_eq -1)                                         \
      return vmrs[deriv_pos] * PVAR *                                \
             mdata[deriv_pos].compute##DERIV(T, T0, Variable::XVAR); \
    else                                                             \
      return 0;                                                      \
  }
  LSPDC(G0, _dT0, P)
  LSPDC(G0, _dX0, P)
  LSPDC(G0, _dX1, P)
  LSPDC(G0, _dX2, P)
  LSPDC(G0, _dX3, P)
  LSPDC(D0, _dT0, P)
  LSPDC(D0, _dX0, P)
  LSPDC(D0, _dX1, P)
  LSPDC(D0, _dX2, P)
  LSPDC(D0, _dX3, P)
  LSPDC(G2, _dT0, P)
  LSPDC(G2, _dX0, P)
  LSPDC(G2, _dX1, P)
  LSPDC(G2, _dX2, P)
  LSPDC(G2, _dX3, P)
  LSPDC(D2, _dT0, P)
  LSPDC(D2, _dX0, P)
  LSPDC(D2, _dX1, P)
  LSPDC(D2, _dX2, P)
  LSPDC(D2, _dX3, P)
  LSPDC(FVC, _dT0, P)
  LSPDC(FVC, _dX0, P)
  LSPDC(FVC, _dX1, P)
  LSPDC(FVC, _dX2, P)
  LSPDC(FVC, _dX3, P)
  LSPDC(ETA, _dT0, 1)
  LSPDC(ETA, _dX0, 1)
  LSPDC(ETA, _dX1, 1)
  LSPDC(ETA, _dX2, 1)
  LSPDC(ETA, _dX3, 1)
  LSPDC(Y, _dT0, P)
  LSPDC(Y, _dX0, P)
  LSPDC(Y, _dX1, P)
  LSPDC(Y, _dX2, P)
  LSPDC(Y, _dX3, P)
  LSPDC(G, _dT0, P* P)
  LSPDC(G, _dX0, P* P)
  LSPDC(G, _dX1, P* P)
  LSPDC(G, _dX2, P* P)
  LSPDC(G, _dX3, P* P)
  LSPDC(DV, _dT0, P* P)
  LSPDC(DV, _dX0, P* P)
  LSPDC(DV, _dX1, P* P)
  LSPDC(DV, _dX2, P* P)
  LSPDC(DV, _dX3, P* P)
#undef LSPDC

  /** Compute all shape parameters
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] vmrs The VMR vector as derived by this.vmrs()
   * 
   * @return Shape parameters
   */
                                                              Output
      GetParams(Numeric T, Numeric T0, Numeric P, const Vector& vmrs) const
      noexcept {
    return {G0(T, T0, P, vmrs),
            D0(T, T0, P, vmrs),
            G2(T, T0, P, vmrs),
            D2(T, T0, P, vmrs),
            FVC(T, T0, P, vmrs),
            ETA(T, T0, P, vmrs),
            Y(T, T0, P, vmrs),
            G(T, T0, P, vmrs),
            DV(T, T0, P, vmrs)};
  }

  /** Derivative of GetParams(...) wrt T
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] vmrs The VMR vector as derived by this.vmrs()
   * 
   * @return Derivative of GetParams(...) wrt T
   */
  Output GetTemperatureDerivs(Numeric T,
                              Numeric T0,
                              Numeric P,
                              const Vector& vmrs) const noexcept {
    return {dG0_dT(T, T0, P, vmrs),
            dD0_dT(T, T0, P, vmrs),
            dG2_dT(T, T0, P, vmrs),
            dD2_dT(T, T0, P, vmrs),
            dFVC_dT(T, T0, P, vmrs),
            dETA_dT(T, T0, P, vmrs),
            dY_dT(T, T0, P, vmrs),
            dG_dT(T, T0, P, vmrs),
            dDV_dT(T, T0, P, vmrs)};
  }

  /** Derivative of GetParams(...) wrt VMR
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] pos Position of species in mspecies
   * 
   * @return Derivative of GetParams(...) wrt VMR
   */
  Output GetVMRDerivs(Numeric T, Numeric T0, Numeric P, const Index pos) const
      noexcept {
    return {dG0_dVMR(T, T0, P, pos),
            dD0_dVMR(T, T0, P, pos),
            dG2_dVMR(T, T0, P, pos),
            dD2_dVMR(T, T0, P, pos),
            dFVC_dVMR(T, T0, P, pos),
            dETA_dVMR(T, T0, P, pos),
            dY_dVMR(T, T0, P, pos),
            dG_dVMR(T, T0, P, pos),
            dDV_dVMR(T, T0, P, pos)};
  }
  
  /** Derivative of GetParams(...) wrt Coefficient
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] pos Position of species in mspecies
   * @param[in] vmrs The VMR vector as derived by this.vmrs()
   * @param[in] deriv The derivative
   * 
   * @return Derivative of GetParams(...) wrt Coefficient
   */
  Numeric GetInternalDeriv(Numeric T,
                           Numeric T0,
                           Numeric P,
                           Index pos,
                           const Vector& vmrs,
                           JacPropMatType deriv) const noexcept {
    if (pos < 0) return 0;

#define RETURNINTERNALDERIVATIVE(TYPE)         \
  case JacPropMatType::LineShape##TYPE##X0:    \
    return d##TYPE##_dX0(T, T0, P, pos, vmrs); \
  case JacPropMatType::LineShape##TYPE##X1:    \
    return d##TYPE##_dX1(T, T0, P, pos, vmrs); \
  case JacPropMatType::LineShape##TYPE##X2:    \
    return d##TYPE##_dX2(T, T0, P, pos, vmrs); \
  case JacPropMatType::LineShape##TYPE##X3:    \
    return d##TYPE##_dX3(T, T0, P, pos, vmrs)
    switch (deriv) {
      RETURNINTERNALDERIVATIVE(G0);
      RETURNINTERNALDERIVATIVE(D0);
      RETURNINTERNALDERIVATIVE(G2);
      RETURNINTERNALDERIVATIVE(D2);
      RETURNINTERNALDERIVATIVE(FVC);
      RETURNINTERNALDERIVATIVE(ETA);
      RETURNINTERNALDERIVATIVE(Y);
      RETURNINTERNALDERIVATIVE(G);
      RETURNINTERNALDERIVATIVE(DV);
      default:
        return 0;
    }
#undef RETURNINTERNALDERIVATIVE
  }

  /** Number of species in Model */
  Index nelem() const { return Index(mdata.size()); }
  
  /** Resize function for Model 
   * 
   * Just resizes, does nothing with the new data
   * 
   * @param[in] n New size of mspecies and mdata
   */
  void resize(Index n) {mdata.resize(n);}
  
  /** Reserve function for Model 
   * 
   * Just reserves, does nothing with the new data
   * 
   * @param[in] n New reserves of mspecies and mdata
   */
  void reserve(Index n) {mdata.reserve(n);}
  
  
  /** The line shape model data */
  const std::vector<SingleSpeciesModel>& Data() const noexcept { return mdata; }
  
  /** The line shape model data reference */
  std::vector<SingleSpeciesModel>& Data() noexcept { return mdata; }

  /** Remove species and data at position
   * 
   * Uses standard algorithm "erase" to remove,
   * see it for behavior when an error occurs
   * 
   * @param[in] i Index of position to remove
   */
  void Remove(Index i, ArrayOfSpeciesTag& specs) {
    mdata.erase(mdata.begin() + i);
    specs.erase(specs.begin() + i);
  }

  /** Sets the same line mixing model to all species
   * 
   * Looks at x and sets it Y, G and DV values to all the
   * Y, G, and DV values in this
   * 
   * If LM_AER is used, Interp data is written over as well
   * otherwise it remains untouched
   * 
   * @param[in] x Model containing new line mixing data
   */
  void SetLineMixingModel(SingleSpeciesModel x) {
    for (auto& ssm : mdata) {
      ssm.Y() = x.Y();
      ssm.G() = x.G();
      ssm.DV() = x.DV();
    }
  }
  
  bool Match(const Model& other) const noexcept {
    return std::equal (mdata.cbegin(), mdata.cend(), other.mdata.cbegin(), other.mdata.cend(), [](auto& a, auto& b) {return a.MatchTypes(b);});
  }
  
  friend
  std::istream& from_linefunctiondata(std::istream& data,
                                      Type& type,
                                      bool& self,
                                      bool& bath,
                                      Model& m,
                                      ArrayOfSpeciesTag& species);
  
  friend
  std::istream& from_artscat4(std::istream& is,
                              Type& type,
                              bool& self,
                              bool& bath,
                              Model& m,
                              ArrayOfSpeciesTag& species,
                              const QuantumIdentifier& qid);
  
  
  /** Binary read for Model */
  bifstream& read(bifstream& bif) {
    for (auto& data: mdata)
      data.read(bif);
    return bif;
  }
  
  /** Binary write for Model */
  bofstream& write(bofstream& bof) const {
    for (auto& data: mdata)
      data.write(bof);
    return bof;
  }
};  // Model;

std::ostream& operator<<(std::ostream&, const Model&);
std::istream& operator>>(std::istream&, Model&);

String ModelShape2MetaData(const Model& m);
Model MetaData2ModelShape(const String& s);

ArrayOfString ModelMetaDataArray(const Model& m, const bool self, const bool bath, const ArrayOfSpeciesTag& sts, const Numeric T0);

std::istream& from_artscat4(std::istream& is,
                            Type& type,
                            bool& self,
                            bool& bath,
                            Model& m,
                            ArrayOfSpeciesTag& species,
                            const QuantumIdentifier& qid);

std::istream& from_linefunctiondata(std::istream& data,
                                    Type& type,
                                    bool& self,
                                    bool& bath,
                                    Model& m,
                                    ArrayOfSpeciesTag& species);

/** Legacy reading of old deprecated LineMixingData class */
std::istream& from_linemixingdata(std::istream& data, Model& lsc);

/** Legacy reading of old deprecated PressureBroadeningData class */
std::istream& from_pressurebroadeningdata(std::istream& data,
                                          LineShape::Type& type,
                                          bool& self,
                                          bool& bath,
                                          Model& m,
                                          ArrayOfSpeciesTag& species,
                                          const QuantumIdentifier& qid);

/** Legacy dealing with reading old LineFunctionData */
namespace LegacyLineFunctionData {
/** Length per variable for temperature model */
inline Index temperaturemodel2legacynelem(TemperatureModel type) noexcept {
  switch (type) {
    case TemperatureModel::None:
      return 0;
    case TemperatureModel::T0:
      return 1;
    case TemperatureModel::T1:
      return 2;
    case TemperatureModel::T2:
      return 3;
    case TemperatureModel::T3:
      return 2;
    case TemperatureModel::T4:
      return 3;
    case TemperatureModel::T5:
      return 2;
    case TemperatureModel::LM_AER:
      return 12;
    case TemperatureModel::DPL:
      return 4;
  }
  std::terminate();  // Not allowed to reach, fix higher level code
}

/** Line shape models from string */
inline std::vector<Variable> lineshapetag2variablesvector(String type) {
  if (type == String("DP"))
    return {};
  else if (type == String("LP"))
    return {Variable::G0, Variable::D0};
  else if (type == String("VP"))
    return {Variable::G0, Variable::D0};
  else if (type == String("SDVP"))
    return {Variable::G0, Variable::D0, Variable::G2, Variable::D2};
  else if (type == String("HTP"))
    return {Variable::G0,
            Variable::D0,
            Variable::G2,
            Variable::D2,
            Variable::FVC,
            Variable::ETA};
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}

/** Line mixing models from string */
inline std::vector<Variable> linemixingtag2variablesvector(String type) {
  if (type == "#")
    return {};
  else if (type == "LM1")
    return {Variable::Y};
  else if (type == "LM2")
    return {Variable::Y, Variable::G, Variable::DV};
  else if (type == "INT")
    return {};
  else if (type == "ConstG")
    return {Variable::G};
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}
};  // namespace LegacyLineFunctionData

/** Legacy dealing with reading old LineMixingData */
namespace LegacyLineMixingData {
/** Line mixing types that used to exist */
enum class TypeLM {
  LM_NONE,                  // Reserved for no line mixing
  LM_LBLRTM,                // Reserved for LBLRTM line mixing
  LM_LBLRTM_O2NonResonant,  // Reserved for the non-resonant O2 line in LBLRTM
  LM_1STORDER,  // Reserved for Tretyakov et al. 2005 1st order of line mixing
  LM_2NDORDER,  // Reserved for Makarov et al. 2011 second order of line mixing
  LM_BYBAND  // Reserved for Paris data of relaxation matrix line mixing for band
};

/** Line mixing types from string */
inline LegacyLineMixingData::TypeLM string2typelm(String type) {
  if (type == "NA")  // The standard case
    return TypeLM::LM_NONE;
  else if (type == "LL")  // The LBLRTM case
    return TypeLM::LM_LBLRTM;
  else if (type == "NR")  // The LBLRTM O2 non-resonant case
    return TypeLM::LM_LBLRTM_O2NonResonant;
  else if (type == "L2")  // The 2nd order case
    return TypeLM::LM_2NDORDER;
  else if (type == "L1")  // The 2nd order case
    return TypeLM::LM_1STORDER;
  else if (type == "BB")  // The band class
    return TypeLM::LM_BYBAND;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}

/** Line mixing types to number */
inline Index typelm2nelem(LegacyLineMixingData::TypeLM type) {
  switch (type) {
    case TypeLM::LM_NONE:  // The standard case
      return 0;
    case TypeLM::LM_LBLRTM:  // The LBLRTM case
      return 12;
    case TypeLM::LM_LBLRTM_O2NonResonant:  // Nonresonant is just a tag
      return 1;
    case TypeLM::LM_2NDORDER:  // The 2nd order case
      return 10;
    case TypeLM::LM_1STORDER:  // The 2nd order case
      return 3;
    case TypeLM::LM_BYBAND:  // The band class
      return 1;
  }
  std::terminate();
}

/** LineShape::Model from legacy input vector */
Model vector2modellm(Vector x, LegacyLineMixingData::TypeLM type);
};  // namespace LegacyLineMixingData

/** Legacy dealing with reading old PressureBroadeningData */
namespace LegacyPressureBroadeningData {
/** Pressure broadening types that used to exist */
enum class TypePB {
  PB_NONE,                      // No pressure broadening
  PB_AIR_BROADENING,            // Air broadening and self broadening only
  PB_AIR_AND_WATER_BROADENING,  // Air, water, and self broadening
  PB_PLANETARY_BROADENING,  // Gas broadening as done for solar system planets
  // PB_SD_AIR_VOLUME,                 // HTP in air for SD limit NOT SUPPORTED
  // PB_HTP_AIR_VOLUME,                // HTP in air NOT SUPPORTED
  // PB_VOIGT_TEST_WATER,              // Voigt parameters for testing NOT SUPPORTED
  // PB_SD_TEST_WATER,                 // SD parameters for testing NOT SUPPORTED
  // PB_PURELY_FOR_TESTING             // Testing tag for new input structures --- can be changed by anyone... NOT SUPPORTED
};

/** Pressure broadening types from string */
inline LegacyPressureBroadeningData::TypePB string2typepb(String type) {
  if (type == "NA")  // The none case
    return TypePB::PB_NONE;
  else if (type == "N2")  // Air Broadening is N2 broadening mostly...
    return TypePB::PB_AIR_BROADENING;
  else if (type == "WA")  // Water and Air Broadening
    return TypePB::PB_AIR_AND_WATER_BROADENING;
  else if (type == "AP")  // Planetary broadening
    return TypePB::PB_PLANETARY_BROADENING;
  else {
    std::ostringstream os;
    os << "Type: " << type << ", is not accepted.  "
       << "See documentation for accepted types\n";
    throw std::runtime_error(os.str());
  }
}

/** Pressure broadening if self exist */
inline Index self_listed(const QuantumIdentifier& qid,
                         LegacyPressureBroadeningData::TypePB t) {
  if (t == TypePB::PB_PLANETARY_BROADENING and
      (qid.Species() == SpeciesTag(String("N2")).Species() or
       qid.Species() == SpeciesTag(String("O2")).Species() or
       qid.Species() == SpeciesTag(String("H2O")).Species() or
       qid.Species() == SpeciesTag(String("CO2")).Species() or
       qid.Species() == SpeciesTag(String("H2")).Species() or
       qid.Species() == SpeciesTag(String("He")).Species()))
    return true;
  else if (t == TypePB::PB_AIR_AND_WATER_BROADENING and
           qid.Species() == SpeciesTag(String("H2O")).Species())
    return true;
  else
    return false;
}

/** Pressure broadening types to number of elements */
inline Index typepb2nelem(LegacyPressureBroadeningData::TypePB type) {
  switch (type) {
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

/** LineShape::Model from legacy input vector */
void vector2modelpb(LineShape::Type& mtype,
                    bool& self,
                    bool& bath,
                    Model& m,
                    ArrayOfSpeciesTag& species,
                    Vector x,
                    LegacyPressureBroadeningData::TypePB type,
                    bool self_in_list);
};  // namespace LegacyPressureBroadeningData

};  // namespace LineShape

#endif  // lineshapemodel_h

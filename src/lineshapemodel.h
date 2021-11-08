/* Copyright (C) 2018
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

#include "constants.h"
#include "enums.h"
#include "file.h"
#include "jacobian.h"
#include "species_tags.h"
#include <algorithm>
#include <numeric>
#include <utility>

/** Return the derivative type based on string input 
 * 
 * @param[in] var Variable good with LineShape::toVariable()
 * @param[in] coeff Coefficient good with Options::toLineShapeCoeff()
 * 
 * @return Derivative type
 */
Jacobian::Line select_derivativeLineShape(const String& var,
                                          const String& coeff);

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
 * 
 * FIXME:  The python API breaks if this is a char type even though it can be????
 */
ENUMCLASS(TemperatureModel, char,
  None,   // 0
  T0,     // Constant, X0
  T1,     // Standard, X0 * (T0/T) ^ X1
  T2,     // X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0));
  T3,     // X0 + X1 * (T - T0)
  T4,     // (X0 + X1 * (T0/T - 1)) * (T0/T)^X2;
  T5,     // X0 * (T0/T)^(0.25 + 1.5*X1)
  LM_AER, // X(200) = X0; X(250) = X1; X(298) = X2; X(340) = X3;  Linear interpolation in between
  DPL,    // X0 * (T0/T) ^ X1 + X2 * (T0/T) ^ X3
  POLY    // X0 + X1 * T + X2 * T ^ 2 + X3 * T ^ 3
)

/** List of possible shape variables
 * 
 * Should correspond to strings in AllLineShapeVars()
 */
ENUMCLASS(Variable, char,
  G0,   // Pressure broadening speed-independent
  D0,   // Pressure f-shifting speed-dependent
  G2,   // Pressure broadening speed-dependent
  D2,   // Pressure f-shifting speed-independent
  FVC,  // Frequency of velocity-changing collisions
  ETA,  // Correlation
  Y,    // First order line mixing coefficient
  G,    // Second order line mixing coefficient
  DV    // Second order line mixing f-shifting
)

/** Coefficients and temperature model for SingleSpeciesModel 
 * 
 * NOTE: Developer should always add new coefficients at the end
 */
struct ModelParameters {
  static constexpr Index N = 4;
  static_assert(Index(Options::LineShapeCoeff::FINAL) == N, "Must update either LineShapeCoeff options or ModelParameters");
  TemperatureModel type;
  Numeric X0;
  Numeric X1;
  Numeric X2;
  Numeric X3;
  
  constexpr ModelParameters(TemperatureModel intype=TemperatureModel::None,
                            Numeric inX0=std::numeric_limits<Numeric>::quiet_NaN(),
                            Numeric inX1=std::numeric_limits<Numeric>::quiet_NaN(),
                            Numeric inX2=std::numeric_limits<Numeric>::quiet_NaN(),
                            Numeric inX3=std::numeric_limits<Numeric>::quiet_NaN())
  noexcept : type(intype), X0(inX0), X1(inX1), X2(inX2), X3(inX3) {}
  
  template <typename VectorType> constexpr
  ModelParameters(TemperatureModel intype, VectorType&& v) ARTS_NOEXCEPT :
  ModelParameters(intype) {
    const auto n = std::size(v);
    ARTS_ASSERT(n <= N, "Must have at most ", N, " inputs, got: ", n)
    switch (n) {
      case 4: X3 = v[3]; [[fallthrough]];
      case 3: X2 = v[2]; [[fallthrough]];
      case 2: X1 = v[1]; [[fallthrough]];
      case 1: X0 = v[0];
    }
  }
  
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
  [[nodiscard]] constexpr Numeric special_linemixing_aer(Numeric T) const noexcept {
    if (T < 250.0)
      return X0 + (T - 200.0) * (X1 - X0) / (250.0 - 200.0);
    if (T > 296.0)
      return X2 + (T - 296.0) * (X3 - X2) / (340.0 - 296.0);
    return X1 + (T - 250.0) * (X2 - X1) / (296.0 - 250.0);
  }
  
  /** The temperature derivative of special_linemixing_aer
   * 
   * @param[in] T The temperature
   * @param[in] var The variable
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  [[nodiscard]] constexpr Numeric special_linemixing_aer_dT(Numeric T) const noexcept {
    if (T < 250.0)
      return (X1 - X0) / (250.0 - 200.0);
    if (T > 296.0)
      return (X3 - X2) / (340.0 - 296.0);
    return (X2 - X1) / (296.0 - 250.0);
  }
  
  /** The derivative of special_linemixing_aer wrt X0
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX0(Numeric T) noexcept {
    if (T < 250.0)
      return 1 - (T - 200.0) / (250.0 - 200.0);
    return 0;
  }
  
  /** The derivative of special_linemixing_aer wrt X1
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX1(Numeric T) noexcept {
    if (T < 250.0)
      return (T - 200.0) / (250.0 - 200.0);
    if (T > 296.0)
      return 0;
    return 1 - (T - 250.0) / (296.0 - 250.0);
  }
  
  /** The derivative of special_linemixing_aer wrt X2
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX2(Numeric T) noexcept {
    if (T < 250.0)
      return 0;
    if (T > 296.0)
      return 1 - (T - 296.0)  / (340.0 - 296.0);
    return (T - 250.0) / (296.0 - 250.0);
  }
  
  /** The derivative of special_linemixing_aer wrt X3
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX3(Numeric T) noexcept {
    if (T > 296.0)
      return (T - 296.0) / (340.0 - 296.0);
    return 0;
  }
  
  [[nodiscard]] Numeric at(Numeric T, Numeric T0) const noexcept;
  
  [[nodiscard]] [[nodiscard]] Numeric dX0(Numeric T, Numeric T0) const noexcept;
  
  [[nodiscard]] Numeric dX1(Numeric T, Numeric T0) const noexcept;
  
  [[nodiscard]] Numeric dX2(Numeric T, Numeric T0) const noexcept;
  
  [[nodiscard]] Numeric dX3(Numeric T, Numeric T0) const noexcept;
  
  [[nodiscard]] Numeric dT(Numeric T, Numeric T0) const noexcept;
  
  [[nodiscard]] Numeric dT0(Numeric T, Numeric T0) const noexcept;
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
Numeric& SingleModelParameter(ModelParameters& mp, const String& type);

constexpr bool modelparameterEmpty(const ModelParameters mp) noexcept {
  switch(mp.type) {
    case TemperatureModel::None:   // 0
      return true;
    case TemperatureModel::T0:     // Constant, X0
      return (mp.X0 == 0);
    case TemperatureModel::T1:     // Standard, X0 * (T0/T) ^ X1
      return (mp.X0 == 0);
    case TemperatureModel::T2:     // X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0));
      return (mp.X0 == 0);
    case TemperatureModel::T3:     // X0 + X1 * (T - T0)
      return (mp.X0 == 0 and mp.X1 == 0);
    case TemperatureModel::T4:     // (X0 + X1 * (T0/T - 1)) * (T0/T)^X2;
      return (mp.X0 == 0 and mp.X1 == 0);
    case TemperatureModel::T5:     // X0 * (T0/T)^(0.25 + 1.5*X1)
      return (mp.X0 == 0);
    case TemperatureModel::LM_AER: // X(200) = X0; X(250) = X1; X(298) = X2; X(340) = X3;  Linear interpolation in between
      return (mp.X0 == 0 and mp.X1 == 0 and mp.X2 == 0 and mp.X3 == 0);
    case TemperatureModel::DPL:    // X0 * (T0/T) ^ X1 + X2 * (T0/T) ^ X3
      return (mp.X0 == 0 and mp.X2 == 0);
    case TemperatureModel::POLY:
      return (mp.X0 == 0 and mp.X1 == 0 and mp.X2 == 0 and mp.X3 == 0);
    case TemperatureModel::FINAL:
      return true;
  }
  return true;
}

constexpr Numeric modelparameterFirstExponent(const ModelParameters mp) noexcept {
  switch(mp.type) {
    case TemperatureModel::None:   // 0
      return 0;
    case TemperatureModel::T0:     // Constant, X0
      return 0;
    case TemperatureModel::T1:     // Standard, X0 * (T0/T) ^ X1
      return mp.X1;
    case TemperatureModel::T2:     // X0 * (T0/T) ^ X1 * (1 + X2 * log(T/T0));
      return mp.X1;
    case TemperatureModel::T3:     // X0 + X1 * (T - T0)
      return 0;
    case TemperatureModel::T4:     // (X0 + X1 * (T0/T - 1)) * (T0/T)^X2;
      return mp.X2;
    case TemperatureModel::T5:     // X0 * (T0/T)^(0.25 + 1.5*X1)
      return (0.25 + 1.5*mp.X1);
    case TemperatureModel::LM_AER: // X(200) = X0; X(250) = X1; X(298) = X2; X(340) = X3;  Linear interpolation in between
      return 0;
    case TemperatureModel::DPL:    // X0 * (T0/T) ^ X1 + X2 * (T0/T) ^ X3
      return mp.X1;
    case TemperatureModel::POLY:
      return 0;
    case TemperatureModel::FINAL:
      return std::numeric_limits<Numeric>::quiet_NaN();
  }
  return std::numeric_limits<Numeric>::quiet_NaN();
}

/** Output operator for ModelParameters */
std::ostream& operator<<(std::ostream& os, const ModelParameters& mp);

/** Input operator for ModelParameters */
std::istream& operator>>(std::istream& is, ModelParameters& mp);

/** Current max number of line shape variables */
constexpr Index nVars = Index(Variable::FINAL);

/** Compute the line shape parameters for a single broadening species */
class SingleSpeciesModel {
 private:
  std::array<ModelParameters, nVars> X;

 public:
  /** Default initialization */
  constexpr SingleSpeciesModel(
    ModelParameters G0=ModelParameters{},
    ModelParameters D0=ModelParameters{},
    ModelParameters G2=ModelParameters{},
    ModelParameters D2=ModelParameters{},
    ModelParameters FVC=ModelParameters{},
    ModelParameters ETA=ModelParameters{},
    ModelParameters Y=ModelParameters{},
    ModelParameters G=ModelParameters{},
    ModelParameters DV=ModelParameters{})
      : X({G0, D0, G2, D2, FVC, ETA, Y, G, DV}) {}

#define ACCESS_INTERNAL(VARPOS)                                             \
  constexpr ModelParameters& VARPOS() noexcept { return std::get<Index(Variable::VARPOS)>(X); } \
  constexpr ModelParameters VARPOS() const noexcept { return std::get<Index(Variable::VARPOS)>(X); }
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
  constexpr std::array<ModelParameters, nVars>& Data() noexcept { return X; }
  
  /** Get const internal Data reference */
  [[nodiscard]] constexpr const std::array<ModelParameters, nVars>& Data() const noexcept { return X; }

  /** Set variable to a different ModelParameters
   * 
   * @param[in] var The variable
   * @param[in] x The new ModelParameters for var
   */
  constexpr void Set(Variable var, const ModelParameters& x) noexcept {
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
      case Variable::FINAL: break;
    }
#undef MODELPARAMCASESETTER
  }

  /** Get variable by type
   * 
   * @param[in] var The variable type
   * 
   * @return ModelParameters copy
   */
  [[nodiscard]] constexpr ModelParameters Get(Variable var) const noexcept {
  #define MODELPARAMCASEGETTER(X) case Variable::X: out = X(); break;  
  ModelParameters out{};
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
    case Variable::FINAL: break;
  }
  return out;
  #undef MODELPARAMCASEGETTER
  }

  /** Binary read for SingleSpeciesModel */
  bifstream& read(bifstream& bif);

  /** Binary write for SingleSpeciesModel */
  bofstream& write(bofstream& bof) const;
  
  [[nodiscard]] bool MatchTypes(const SingleSpeciesModel& other) const noexcept;
};

/** Output operator for SingleSpeciesModel */
std::ostream& operator<<(std::ostream& os, const SingleSpeciesModel& ssm);

/** Input operator for SingleSpeciesModel */
std::istream& operator>>(std::istream& is, SingleSpeciesModel& ssm);

/** Type of line shape to compute */
ENUMCLASS(Type, char,
  DP,    // Doppler
  LP,    // Lorentz
  VP,    // Voigt
  SDVP,  // Speed-dependent Voigt
  HTP    // Hartmann-Tran
)

/** Turns selected Type into a human readable string
 * 
 * This function takes the input Type
 * and returns it as a string
 *  
 * @param[in] type The Type
 * 
 * @return std::string_view of Type
 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr std::string_view shapetype2metadatastring(Type type) noexcept {
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
    case Type::FINAL:
      return "There's an error.\n";
  }
}
#pragma GCC diagnostic pop

/** Main output of Model */
struct Output {
  Numeric G0{0},  // Pressure broadening speed-independent
      D0{0},      // Pressure f-shifting speed-independent
      G2{0},      // Pressure broadening speed-dependent
      D2{0},      // Pressure f-shifting speed-dependent
      FVC{0},     // Frequency of velocity-changing collisions
      ETA{0},     // Correlation
      Y{0},       // First order line mixing coefficient
      G{0},       // Second order line mixing coefficient
      DV{0};      // Second order line mixing f-shifting
  constexpr Output() noexcept = default;
  constexpr Output(Numeric g0, Numeric d0, Numeric g2,
                   Numeric d2, Numeric fvc, Numeric eta,
                   Numeric y, Numeric g, Numeric dv) noexcept :
    G0(g0), D0(d0), G2(g2), D2(d2), FVC(fvc), ETA(eta), Y(y), G(g), DV(dv) {}
  constexpr Output(const Output&) noexcept = default;
  constexpr Output(Output&&) noexcept = default;
  constexpr Output& operator=(const Output&) noexcept = default;
  constexpr Output& operator=(Output&&) noexcept = default;
};

/** Output operator for LineShape::Output */
std::ostream& operator<<(std::ostream& os, Output x);

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
 * @param[in] lineshape_species Species affecting lineshape
 */
Vector vmrs(const ConstVectorView& atmospheric_vmrs,
            const ArrayOfArrayOfSpeciesTag& atmospheric_species,
            const ArrayOfSpecies& lineshape_species) ARTS_NOEXCEPT;

/** Returns a mass vector for this model's main calculations
 * 
 * Sets a vector that matches the mdata size of VMRs based
 * on atmospheric species and VMRs
 * 
 * Only checks the first species in inner atmosphere
 * 
 * Renormalizes the values to unity.  If this renormalization
 * is impossible then it throws an error
 * 
 * @param[in] atmospheric_vmrs VMRS in atmosphere
 * @param[in] atmospheric_species Species in atmosphere
 * @param[in] lineshape_species Species affecting lineshape
 */
Vector mass(const ConstVectorView& atmospheric_vmrs,
            const ArrayOfArrayOfSpeciesTag& atmospheric_species,
            const ArrayOfSpecies& lineshape_species,
            const SpeciesIsotopologueRatios& ir) ARTS_NOEXCEPT;

/** Name for bath broadening in printing and reading user input */
static constexpr std::string_view bath_broadening = "AIR";

/** Name for self broadening in printing and reading user input */
static constexpr std::string_view self_broadening = "SELF";

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
  explicit Model(const std::vector<SingleSpeciesModel>& assm) noexcept : mdata(assm) {}
  
  /** Init from copying itself */
  Model(const Model& m) noexcept : Model(m.mdata) {}
  
  /** Init from moving a vector */
  explicit Model(std::vector<SingleSpeciesModel>&& assm) noexcept : mdata(std::move(assm)) {}
  
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
        std::array<Numeric, 12> aer_interp = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) noexcept;
  
  /** The Model is good to use
   * 
   * @return true/false
   */
  [[nodiscard]] bool OK(Type type, bool self, bool bath,
          const std::size_t nspecies) const noexcept;

#define LSPC(XVAR, PVAR)                                                       \
  Numeric XVAR(                                                                \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], const ConstVectorView& vmrs) \
      const noexcept;
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
                         const Index deriv_pos) const noexcept;
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

#define LSPCT(XVAR, PVAR)                                                       \
  Numeric d##XVAR##_dT(                                                         \
      Numeric T, Numeric T0, Numeric P [[maybe_unused]], const ConstVectorView& vmrs)  \
      const noexcept;
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
                         const ConstVectorView& vmrs) const noexcept;
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
  [[nodiscard]] Output GetParams(Numeric T,
                   Numeric T0,
                   Numeric P,
                   const ConstVectorView& vmrs) const noexcept;

  /** Compute all shape parameters
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] k The position of the single species shape parameters
   * 
   * @return Shape parameters
   */
  [[nodiscard]] Output GetParams(Numeric T,
                   Numeric T0,
                   Numeric P,
                   size_t k) const noexcept;

  /** Derivative of GetParams(...) wrt T
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] vmrs The VMR vector as derived by this.vmrs()
   * 
   * @return Derivative of GetParams(...) wrt T
   */
  [[nodiscard]] Output GetTemperatureDerivs(Numeric T,
                              Numeric T0,
                              Numeric P,
                              ConstVectorView vmrs) const noexcept;

  /** Derivative of GetParams(...) wrt VMR
   * 
   * @param[in] T The temperature
   * @param[in] T0 The temperature used to derive the coefficients
   * @param[in] P The pressure
   * @param[in] pos Position of species in mspecies
   * 
   * @return Derivative of GetParams(...) wrt VMR
   */
  [[nodiscard]] Output GetVMRDerivs(Numeric T, Numeric T0, Numeric P, const Index pos) const noexcept;
  
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
  [[nodiscard]] Numeric GetInternalDeriv(Numeric T,
                           Numeric T0,
                           Numeric P,
                           Index pos,
                           const Vector& vmrs,
                           Jacobian::Line deriv) const noexcept;

  /** Number of species in Model */
  [[nodiscard]] Index nelem() const { return Index(mdata.size()); }
  
  /** Number of species in Model */
  [[nodiscard]] Index size() const { return Index(mdata.size()); }
  
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
  
  /** Get a SingleSpeciesModel
   * 
   * @param[in] i Position in mdata
   * @return reference to SingleSpeciesModel
   */
  SingleSpeciesModel& operator[](Index i) {return mdata[i];}
  
  /** Get a SingleSpeciesModel
   * 
   * @param[in] i Position in mdata
   * @return reference to SingleSpeciesModel
   */
  const SingleSpeciesModel& operator[](Index i) const {return mdata[i];}
  
  
  /** The line shape model data */
  [[nodiscard]] const std::vector<SingleSpeciesModel>& Data() const noexcept { return mdata; }
  
  /** The line shape model data reference */
  std::vector<SingleSpeciesModel>& Data() noexcept { return mdata; }

  /** Remove species and data at position
   * 
   * Uses standard algorithm "erase" to remove,
   * see it for behavior when an error occurs
   * 
   * @param[in] i Index of position to remove
   */
  void Remove(Index i, ArrayOfSpeciesTag& specs);
  void Remove(Index i, ArrayOfSpecies& specs);

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
  void SetLineMixingModel(SingleSpeciesModel x);
  
  [[nodiscard]] bool Match(const Model& other) const noexcept;
  
  friend
  std::istream& from_linefunctiondata(std::istream& data,
                                      Type& type,
                                      bool& self,
                                      bool& bath,
                                      Model& m,
                                      ArrayOfSpecies& species);
  
  friend
  std::istream& from_artscat4(std::istream& is,
                              Type& type,
                              bool& self,
                              bool& bath,
                              Model& m,
                              ArrayOfSpecies& species,
                              const QuantumIdentifier& qid);
  
  
  /** Binary read for Model */
  bifstream& read(bifstream& bif);
  
  /** Binary write for Model */
  bofstream& write(bofstream& bof) const;
};  // Model;


Model hitran_model(Numeric sgam,
                   Numeric nself,
                   Numeric agam,
                   Numeric nair,
                   Numeric psf);

Model lblrtm_model(Numeric sgam,
                   Numeric nself,
                   Numeric agam,
                   Numeric nair,
                   Numeric psf,
                   std::array<Numeric, 12> aer_interp = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

std::ostream& operator<<(std::ostream&, const Model&);

std::istream& operator>>(std::istream&, Model&);

String ModelShape2MetaData(const Model& m);

Model MetaData2ModelShape(const String& s);

ArrayOfString ModelMetaDataArray(const Model& m, const bool self, const ArrayOfSpecies& sts, const Numeric T0);

std::istream& from_artscat4(std::istream& is,
                            Type& type,
                            bool& self,
                            bool& bath,
                            Model& m,
                            ArrayOfSpecies& species,
                            const QuantumIdentifier& qid);

std::istream& from_linefunctiondata(std::istream& data,
                                    Type& type,
                                    bool& self,
                                    bool& bath,
                                    Model& m,
                                    ArrayOfSpecies& species);

/** Legacy reading of old deprecated LineMixingData class */
std::istream& from_linemixingdata(std::istream& data, Model& lsc);

/** Legacy reading of old deprecated PressureBroadeningData class */
std::istream& from_pressurebroadeningdata(std::istream& data,
                                          LineShape::Type& type,
                                          bool& self,
                                          bool& bath,
                                          Model& m,
                                          ArrayOfSpecies& species,
                                          const QuantumIdentifier& qid);

/** Legacy dealing with reading old LineFunctionData */
namespace LegacyLineFunctionData {
/** Length per variable for temperature model */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr Index temperaturemodel2legacynelem(TemperatureModel type) noexcept {
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
    case TemperatureModel::POLY:
      return 4;
    case TemperatureModel::FINAL: break;
  }
}
#pragma GCC diagnostic pop

/** Line shape models from string */
std::vector<Variable> lineshapetag2variablesvector(String type);

/** Line mixing models from string */
std::vector<Variable> linemixingtag2variablesvector(String type);
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
LegacyLineMixingData::TypeLM string2typelm(String type);

/** Line mixing types to number */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr Index typelm2nelem(LegacyLineMixingData::TypeLM type) {
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
}
#pragma GCC diagnostic pop

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
LegacyPressureBroadeningData::TypePB string2typepb(String type);

/** Pressure broadening if self exist */
Index self_listed(const QuantumIdentifier& qid,
                  LegacyPressureBroadeningData::TypePB t);

/** Pressure broadening types to number of elements */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
constexpr Index typepb2nelem(LegacyPressureBroadeningData::TypePB type)  {
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
}
#pragma GCC diagnostic pop

/** LineShape::Model from legacy input vector */
void vector2modelpb(LineShape::Type& mtype,
                    bool& self,
                    bool& bath,
                    Model& m,
                    ArrayOfSpecies& species,
                    Vector x,
                    LegacyPressureBroadeningData::TypePB type,
                    bool self_in_list,
                    Species::Species self_spec);
};  // namespace LegacyPressureBroadeningData
};  // namespace LineShape

using LineShapeModelParameters = LineShape::ModelParameters;
using LineShapeModel = LineShape::Model;
using LineShapeSingleSpeciesModel = LineShape::SingleSpeciesModel;

using LineShapeType = LineShape::Type;
using LineShapeVariable = LineShape::Variable;
using LineShapeTemperatureModel = LineShape::TemperatureModel;

#endif  // lineshapemodel_h


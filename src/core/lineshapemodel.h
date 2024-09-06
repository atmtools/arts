#pragma once

#include <arts_conversions.h>
#include <atm.h>
#include <enumsLineShapeTemperatureModelOld.h>
#include <enumsLineShapeTypeOld.h>
#include <enumsLineShapeVariableOld.h>
#include <format_tags.h>
#include <species_tags.h>

#include <algorithm>
#include <string_view>
#include <utility>

/** Is the Type independent per broadener?
 * 
 * If it is, the line shape is computed as F_1(...) + F_2(...) + /// + F_N(...)
 * If it is not, then it is computed as F(avg(..._1 + ..._2 + /// + ..._N))
 *
 * Where F is the lineshape, the subindex is the line parameter (such as pressure broadening)
 * and avg(///) is an averaging function
 *  
 * @param[in] type The LineShape::Type
 * 
 * @return true If we compute the line shape per broadener and then sum it up
 * @return false If we average line parameters for all broadeners and then compute the line shape
 */
constexpr bool independent_per_broadener(LineShapeTypeOld in) {
  using enum LineShapeTypeOld;
  constexpr std::array data{SplitLP, SplitVP, SplitSDVP, SplitHTP};
  static_assert(std::is_sorted(data.begin(), data.end()), "Not sorted");
  return std::binary_search(data.begin(), data.end(), in);
}

/** Computations of line shape derived parameters
 * 
 * Defines many classes and IO routines for line 
 * shape parameters to comply with everything 
 * from no line mixing Doppler to coefficient-based
 * line mixing Hartman-Tran profiles
 */
namespace LineShape {

/** Coefficients and temperature model for SingleSpeciesModel 
 * 
 * NOTE: Developer should always add new coefficients at the end
 */
struct ModelParameters {
  static constexpr Index N = 4;
  static_assert(4 == N,
                "Must update either LineShapeCoeff options or ModelParameters");
  LineShapeTemperatureModelOld type;
  Numeric X0;
  Numeric X1;
  Numeric X2;
  Numeric X3;

  constexpr ModelParameters(
      LineShapeTemperatureModelOld intype = LineShapeTemperatureModelOld::None,
      Numeric inX0 = std::numeric_limits<Numeric>::quiet_NaN(),
      Numeric inX1 = std::numeric_limits<Numeric>::quiet_NaN(),
      Numeric inX2 = std::numeric_limits<Numeric>::quiet_NaN(),
      Numeric inX3 = std::numeric_limits<Numeric>::quiet_NaN()) noexcept
      : type(intype), X0(inX0), X1(inX1), X2(inX2), X3(inX3) {}

  template <typename VectorType>
  constexpr ModelParameters(LineShapeTemperatureModelOld intype,
                            VectorType&& v) 
      : ModelParameters(intype) {
    const auto n = std::size(v);
    ARTS_ASSERT(n <= N, "Must have at most ", N, " inputs, got: ", n)
    switch (n) {
      case 4:
        X3 = v[3];
        [[fallthrough]];
      case 3:
        X2 = v[2];
        [[fallthrough]];
      case 2:
        X1 = v[1];
        [[fallthrough]];
      case 1:
        X0 = v[0];
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
  [[nodiscard]] constexpr Numeric special_linemixing_aer(
      Numeric T) const noexcept {
    if (T < 250.0) return X0 + (T - 200.0) * (X1 - X0) / (250.0 - 200.0);
    if (T > 296.0) return X2 + (T - 296.0) * (X3 - X2) / (340.0 - 296.0);
    return X1 + (T - 250.0) * (X2 - X1) / (296.0 - 250.0);
  }

  /** The temperature derivative of special_linemixing_aer
   * 
   * @param[in] T The temperature
   * @param[in] var The variable
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  [[nodiscard]] constexpr Numeric special_linemixing_aer_dT(
      Numeric T) const noexcept {
    if (T < 250.0) return (X1 - X0) / (250.0 - 200.0);
    if (T > 296.0) return (X3 - X2) / (340.0 - 296.0);
    return (X2 - X1) / (296.0 - 250.0);
  }

  /** The derivative of special_linemixing_aer wrt X0
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX0(Numeric T) noexcept {
    if (T < 250.0) return 1 - (T - 200.0) / (250.0 - 200.0);
    return 0;
  }

  /** The derivative of special_linemixing_aer wrt X1
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX1(Numeric T) noexcept {
    if (T < 250.0) return (T - 200.0) / (250.0 - 200.0);
    if (T > 296.0) return 0;
    return 1 - (T - 250.0) / (296.0 - 250.0);
  }

  /** The derivative of special_linemixing_aer wrt X2
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX2(Numeric T) noexcept {
    if (T < 250.0) return 0;
    if (T > 296.0) return 1 - (T - 296.0) / (340.0 - 296.0);
    return (T - 250.0) / (296.0 - 250.0);
  }

  /** The derivative of special_linemixing_aer wrt X3
   * 
   * @param[in] T The temperature
   * 
   * @return The temperature derivative of the broadening parameter at temperature
   */
  static constexpr Numeric special_linemixing_aer_dX3(Numeric T) noexcept {
    if (T > 296.0) return (T - 296.0) / (340.0 - 296.0);
    return 0;
  }

  [[nodiscard]] Numeric at(Numeric T, Numeric T0) const noexcept;

  [[nodiscard]] Numeric dX0(Numeric T, Numeric T0) const noexcept;

  [[nodiscard]] Numeric dX1(Numeric T, Numeric T0) const noexcept;

  [[nodiscard]] Numeric dX2(Numeric T, Numeric T0) const noexcept;

  [[nodiscard]] Numeric dX3(Numeric T, Numeric T0) const noexcept;

  [[nodiscard]] Numeric dT(Numeric T, Numeric T0) const noexcept;

  [[nodiscard]] Numeric dT0(Numeric T, Numeric T0) const noexcept;

  friend std::ostream& operator<<(std::ostream& os, const ModelParameters& mp);

  friend std::istream& operator>>(std::istream& is, ModelParameters& mp);
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

bool modelparameterEmpty(const ModelParameters mp) noexcept;

/** Current max number of line shape variables */
constexpr Index nVars = 9;

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
  constexpr Output(Numeric g0,
                   Numeric d0,
                   Numeric g2,
                   Numeric d2,
                   Numeric fvc,
                   Numeric eta,
                   Numeric y,
                   Numeric g,
                   Numeric dv) noexcept
      : G0(g0),
        D0(d0),
        G2(g2),
        D2(d2),
        FVC(fvc),
        ETA(eta),
        Y(y),
        G(g),
        DV(dv) {}
  constexpr Output(const Output&) noexcept            = default;
  constexpr Output(Output&&) noexcept                 = default;
  constexpr Output& operator=(const Output&) noexcept = default;
  constexpr Output& operator=(Output&&) noexcept      = default;
  constexpr Output& operator-=(Output&& other) noexcept {
    static_assert(nVars == 9, "Must update");
    G0  -= other.G0;
    D0  -= other.D0;
    G2  -= other.G2;
    D2  -= other.D2;
    FVC -= other.FVC;
    ETA -= other.ETA;
    Y   -= other.Y;
    G   -= other.G;
    DV  -= other.DV;
    return *this;
  }

  //! Turns of line mixing if true.  Return *this
  constexpr Output& no_linemixing(bool do_no_linemixing) {
    if (do_no_linemixing) {
      Y = G = DV = 0;
    }
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, Output x);
};

/** Compute the line shape parameters for a single broadening species */
class SingleSpeciesModel {
 private:
  std::array<ModelParameters, nVars> X;

 public:
  /** Default initialization */
  constexpr SingleSpeciesModel(ModelParameters G0  = ModelParameters{},
                               ModelParameters D0  = ModelParameters{},
                               ModelParameters G2  = ModelParameters{},
                               ModelParameters D2  = ModelParameters{},
                               ModelParameters FVC = ModelParameters{},
                               ModelParameters ETA = ModelParameters{},
                               ModelParameters Y   = ModelParameters{},
                               ModelParameters G   = ModelParameters{},
                               ModelParameters DV  = ModelParameters{})
      : X({G0, D0, G2, D2, FVC, ETA, Y, G, DV}) {}

#define ACCESS_INTERNAL(VARPOS)                              \
  constexpr ModelParameters& VARPOS() noexcept {             \
    return std::get<Index(LineShapeVariableOld::VARPOS)>(X); \
  }                                                          \
  constexpr ModelParameters VARPOS() const noexcept {        \
    return std::get<Index(LineShapeVariableOld::VARPOS)>(X); \
  }
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
  [[nodiscard]] constexpr const std::array<ModelParameters, nVars>& Data()
      const noexcept {
    return X;
  }

  /** Set variable to a different ModelParameters
   * 
   * @param[in] var The variable
   * @param[in] x The new ModelParameters for var
   */
  constexpr void Set(LineShapeVariableOld var,
                     const ModelParameters& x) noexcept {
#define MODELPARAMCASESETTER(X) \
  case LineShapeVariableOld::X: \
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
  [[nodiscard]] constexpr ModelParameters Get(
      LineShapeVariableOld var) const noexcept {
#define MODELPARAMCASEGETTER(X) \
  case LineShapeVariableOld::X: \
    out = X();                  \
    break;
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
    }
    return out;
#undef MODELPARAMCASEGETTER
  }

  /** Binary read for SingleSpeciesModel */
  bifstream& read(bifstream& bif);

  /** Binary write for SingleSpeciesModel */
  bofstream& write(bofstream& bof) const;

  [[nodiscard]] std::pair<bool, bool> MatchTypes(
      const SingleSpeciesModel& other) const noexcept;

  friend std::ostream& operator<<(std::ostream& os,
                                  const SingleSpeciesModel& ssm);

  friend std::istream& operator>>(std::istream& is, SingleSpeciesModel& ssm);

  [[nodiscard]] Output at(Numeric T, Numeric T0, Numeric P) const noexcept;

  [[nodiscard]] Output dT(Numeric T, Numeric T0, Numeric P) const noexcept;

  [[nodiscard]] Output dT0(Numeric T, Numeric T0, Numeric P) const noexcept;
};

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
 * @param[in] atmospheric_vmrs VMRS in atmosphere
 * @param[in] atmospheric_species Species in atmosphere
 * @param[in] lineshape_species Species affecting lineshape
 */
Vector vmrs(const ConstVectorView& atmospheric_vmrs,
            const ArrayOfArrayOfSpeciesTag& atmospheric_species,
            const ArrayOfSpeciesEnum& lineshape_species) ;

/** Returns a VMR vector for a list of Species from a point in the atmosphere
 * 
 * Renormalizes the values to unity if possible or return zeroes
 * 
 * @param[in] atm_point As WSV
 * @param[in] lineshape_species Species affecting lineshape
 */
Vector vmrs(const AtmPoint& atm_point,
            const ArrayOfSpeciesEnum& lineshape_species) ;

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
            const ArrayOfSpeciesEnum& lineshape_species,
            const SpeciesIsotopologueRatios& ir) ;

/** Name for bath broadening in printing and reading user input */
inline constexpr std::string_view bath_broadening = "AIR";

/** Name for self broadening in printing and reading user input */
inline constexpr std::string_view self_broadening = "SELF";

/** Main line shape model class
 * 
 * Computes the line shape parameters for a given atmosphere
 */
class Model {
 private:
  std::vector<SingleSpeciesModel> mdata;

 public:
  /** Default init just sets the size */
  Model(Index n = 0) noexcept : mdata(n) {}

  /** Init from copying a vector */
  explicit Model(const std::vector<SingleSpeciesModel>& assm) noexcept
      : mdata(assm) {}

  /** Init from copying itself */
  Model(const Model& m) noexcept : Model(m.mdata) {}

  /** Init from moving a vector */
  explicit Model(std::vector<SingleSpeciesModel>&& assm) noexcept
      : mdata(std::move(assm)) {}

  /** Init from moving a itself */
  Model(Model&& m) noexcept : Model(std::move(m.mdata)) {}

  /** Copy and equals */
  Model& operator=(const Model& m) = default;

  /** Move and equals */
  Model& operator=(Model&& m) = default;

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
        std::array<Numeric, 12> aer_interp = {
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}) noexcept;

  /** The Model is good to use
   * 
   * @return true/false
   */
  [[nodiscard]] bool OK(LineShapeTypeOld type,
                        bool self,
                        bool bath,
                        const std::size_t nspecies) const noexcept;

  /** Number of species in Model */
  [[nodiscard]] Size size() const { return mdata.size(); }

  /** Resize function for Model 
   * 
   * Just resizes, does nothing with the new data
   * 
   * @param[in] n New size of mspecies and mdata
   */
  void resize(Size n) { mdata.resize(n); }

  /** Reserve function for Model 
   * 
   * Just reserves, does nothing with the new data
   * 
   * @param[in] n New reserves of mspecies and mdata
   */
  void reserve(Index n) { mdata.reserve(n); }

  /** Get a SingleSpeciesModel
   * 
   * @param[in] i Position in mdata
   * @return reference to SingleSpeciesModel
   */
  SingleSpeciesModel& operator[](Index i) { return mdata[i]; }

  /** Get a SingleSpeciesModel
   * 
   * @param[in] i Position in mdata
   * @return reference to SingleSpeciesModel
   */
  const SingleSpeciesModel& operator[](Index i) const { return mdata[i]; }

  [[nodiscard]] auto begin() { return mdata.begin(); }
  [[nodiscard]] auto end() { return mdata.end(); }

  [[nodiscard]] auto begin() const { return mdata.begin(); }
  [[nodiscard]] auto end() const { return mdata.end(); }

  [[nodiscard]] auto cbegin() const { return mdata.cbegin(); }
  [[nodiscard]] auto cend() const { return mdata.cend(); }

  /** The line shape model data */
  [[nodiscard]] const std::vector<SingleSpeciesModel>& Data() const noexcept {
    return mdata;
  }

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
  void Remove(Index i, ArrayOfSpeciesEnum& specs);

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

  [[nodiscard]] std::pair<bool, bool> Match(const Model& other) const noexcept;

  friend std::istream& from_linefunctiondata(std::istream& data,
                                             LineShapeTypeOld& type,
                                             bool& self,
                                             bool& bath,
                                             Model& m,
                                             ArrayOfSpeciesEnum& species);

  friend std::istream& from_artscat4(std::istream& is,
                                     LineShapeTypeOld& type,
                                     bool& self,
                                     bool& bath,
                                     Model& m,
                                     ArrayOfSpeciesEnum& species,
                                     const QuantumIdentifier& qid);

  friend std::ostream& operator<<(std::ostream&, const Model&);

  friend std::istream& operator>>(std::istream&, Model&);

  /** Binary read for Model */
  bifstream& read(bifstream& bif);

  /** Binary write for Model */
  bofstream& write(bofstream& bof) const;

  [[nodiscard]] Numeric G0(Numeric T,
                           Numeric T0,
                           Numeric P,
                           const Vector& vmrs) const ;
  [[nodiscard]] Numeric D0(Numeric T,
                           Numeric T0,
                           Numeric P,
                           const Vector& vmrs) const ;
  [[nodiscard]] Numeric G2(Numeric T,
                           Numeric T0,
                           Numeric P,
                           const Vector& vmrs) const ;
  [[nodiscard]] Numeric D2(Numeric T,
                           Numeric T0,
                           Numeric P,
                           const Vector& vmrs) const ;
  [[nodiscard]] Numeric ETA(Numeric T,
                            Numeric T0,
                            Numeric P,
                            const Vector& vmrs) const ;
  [[nodiscard]] Numeric FVC(Numeric T,
                            Numeric T0,
                            Numeric P,
                            const Vector& vmrs) const ;
  [[nodiscard]] Numeric Y(Numeric T,
                          Numeric T0,
                          Numeric P,
                          const Vector& vmrs) const ;
  [[nodiscard]] Numeric G(Numeric T,
                          Numeric T0,
                          Numeric P,
                          const Vector& vmrs) const ;
  [[nodiscard]] Numeric DV(Numeric T,
                           Numeric T0,
                           Numeric P,
                           const Vector& vmrs) const ;

  [[nodiscard]] Numeric dG0dT(Numeric T,
                              Numeric T0,
                              Numeric P,
                              const Vector& vmrs) const ;
  [[nodiscard]] Numeric dD0dT(Numeric T,
                              Numeric T0,
                              Numeric P,
                              const Vector& vmrs) const ;
  [[nodiscard]] Numeric dG2dT(Numeric T,
                              Numeric T0,
                              Numeric P,
                              const Vector& vmrs) const ;
  [[nodiscard]] Numeric dD2dT(Numeric T,
                              Numeric T0,
                              Numeric P,
                              const Vector& vmrs) const ;
  [[nodiscard]] Numeric dETAdT(Numeric T,
                               Numeric T0,
                               Numeric P,
                               const Vector& vmrs) const ;
  [[nodiscard]] Numeric dFVCdT(Numeric T,
                               Numeric T0,
                               Numeric P,
                               const Vector& vmrs) const ;
  [[nodiscard]] Numeric dYdT(Numeric T,
                             Numeric T0,
                             Numeric P,
                             const Vector& vmrs) const ;
  [[nodiscard]] Numeric dGdT(Numeric T,
                             Numeric T0,
                             Numeric P,
                             const Vector& vmrs) const ;
  [[nodiscard]] Numeric dDVdT(Numeric T,
                              Numeric T0,
                              Numeric P,
                              const Vector& vmrs) const ;
};  // Model;

Model hitran_model(
    Numeric sgam, Numeric nself, Numeric agam, Numeric nair, Numeric psf);

Model lblrtm_model(Numeric sgam,
                   Numeric nself,
                   Numeric agam,
                   Numeric nair,
                   Numeric psf,
                   std::array<Numeric, 12> aer_interp = {
                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

String ModelShape2MetaData(const Model& m);

Model MetaData2ModelShape(const String& s);

ArrayOfString ModelMetaDataArray(const Model& m,
                                 const bool self,
                                 const ArrayOfSpeciesEnum& sts,
                                 const Numeric T0);

std::istream& from_artscat4(std::istream& is,
                            LineShapeTypeOld& type,
                            bool& self,
                            bool& bath,
                            Model& m,
                            ArrayOfSpeciesEnum& species,
                            const QuantumIdentifier& qid);

std::istream& from_linefunctiondata(std::istream& data,
                                    LineShapeTypeOld& type,
                                    bool& self,
                                    bool& bath,
                                    Model& m,
                                    ArrayOfSpeciesEnum& species);

/** Legacy reading of old deprecated LineMixingData class */
std::istream& from_linemixingdata(std::istream& data, Model& lsc);

/** Legacy reading of old deprecated PressureBroadeningData class */
std::istream& from_pressurebroadeningdata(std::istream& data,
                                          LineShapeTypeOld& type,
                                          bool& self,
                                          bool& bath,
                                          Model& m,
                                          ArrayOfSpeciesEnum& species,
                                          const QuantumIdentifier& qid);

/** Legacy dealing with reading old LineFunctionData */
namespace LegacyLineFunctionData {
/** Line shape models from string */
std::vector<LineShapeVariableOld> lineshapetag2variablesvector(String type);

/** Line mixing models from string */
std::vector<LineShapeVariableOld> linemixingtag2variablesvector(String type);
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

/** LineShape::Model from legacy input vector */
void vector2modelpb(LineShapeTypeOld& mtype,
                    bool& self,
                    bool& bath,
                    Model& m,
                    ArrayOfSpeciesEnum& species,
                    Vector x,
                    LegacyPressureBroadeningData::TypePB type,
                    bool self_in_list,
                    SpeciesEnum self_spec);
};  // namespace LegacyPressureBroadeningData
};  // namespace LineShape

using LineShapeModelParameters    = LineShape::ModelParameters;
using LineShapeModel              = LineShape::Model;
using LineShapeSingleSpeciesModel = LineShape::SingleSpeciesModel;

template <>
struct std::formatter<LineShapeModelParameters> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const LineShapeModelParameters& v,
                              FmtContext& ctx) const {
    const std::string_view sep = tags.sep();
    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, v.type, sep, v.X0, sep, v.X1, sep, v.X2, sep, v.X3);
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<LineShapeSingleSpeciesModel> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const LineShapeSingleSpeciesModel& v,
                              FmtContext& ctx) const {
    return tags.format(ctx, v.Data());
  }
};

template <>
struct std::formatter<LineShapeModel> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const LineShapeModel& v, FmtContext& ctx) const {
    return tags.format(ctx, v.Data());
  }
};

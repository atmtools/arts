/** 
  @file   jacobian.h
  @author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  @date   2004-09-14

  @brief  Routines for setting up the jacobian.
 */

#ifndef jacobian_h
#define jacobian_h

#include "arts_conversions.h"
#include "arts_options.h"
#include "array.h"
#include "bifstream.h"
#include "enums.h"
#include "interpolation.h"
#include "logic.h"
#include "matpack_arrays.h"
#include "matpack_data.h"
#include "matpack_sparse.h"
#include "methods.h"
#include "mystring.h"
#include "ppath_struct.h"
#include "quantum_numbers.h"  
#include "species_tags.h"
#include <iostream>
#include <map>
#include <stdexcept>
#include <utility>

namespace Jacobian {

/** Holds the type of the target quantity */
ENUMCLASS(Type, char,
          Atm,
          Line,
          Sensor,
          Special
          )

/** Holds the Atmosphere-related targets */
ENUMCLASS(Atm, char,
          Temperature,
          WindMagnitude, WindU, WindV, WindW,
          MagneticMagnitude, MagneticU, MagneticV, MagneticW,
          Electrons,
          Particulates
          )

/** Holds the Line-related targets */
ENUMCLASS(Line, char,
          VMR,
          Strength,
          Center,
          ShapeG0X0, ShapeG0X1, ShapeG0X2, ShapeG0X3,
          ShapeD0X0, ShapeD0X1, ShapeD0X2, ShapeD0X3,
          ShapeG2X0, ShapeG2X1, ShapeG2X2, ShapeG2X3,
          ShapeD2X0, ShapeD2X1, ShapeD2X2, ShapeD2X3,
          ShapeFVCX0, ShapeFVCX1, ShapeFVCX2, ShapeFVCX3,
          ShapeETAX0, ShapeETAX1, ShapeETAX2, ShapeETAX3,
          ShapeYX0, ShapeYX1, ShapeYX2, ShapeYX3,
          ShapeGX0, ShapeGX1, ShapeGX2, ShapeGX3,
          ShapeDVX0, ShapeDVX1, ShapeDVX2, ShapeDVX3,
          ECS_SCALINGX0, ECS_SCALINGX1, ECS_SCALINGX2, ECS_SCALINGX3,
          ECS_BETAX0, ECS_BETAX1, ECS_BETAX2, ECS_BETAX3,
          ECS_LAMBDAX0, ECS_LAMBDAX1, ECS_LAMBDAX2, ECS_LAMBDAX3,
          ECS_DCX0, ECS_DCX1, ECS_DCX2, ECS_DCX3,
          NLTE
          )
static_assert(
    Index(Line::FINAL) == 4 + 13 * Index(Options::LineShapeCoeff::FINAL),
    "Either you have added some \"Line\" parameter(s) or changed the temperature model");

/** Holds the Sensor-related targets */
ENUMCLASS(Sensor, char,
          FrequencyShift,
          FrequencyStretch,
          Polyfit,
          Sinefit,
          PointingZenithInterp,
          PointingZenithRecalc
          )

/** Holds special targets that require careful extra manipulation */
ENUMCLASS(Special, char,
          ArrayOfSpeciesTagVMR,
          ScatteringString,
          SurfaceString
          )

/** Holds all information required for individual partial derivatives */
struct Target {
  /**! Type of quantity, never set manually */
  Type type{Type::FINAL};

  /** Type of atm quantity */
  Atm atm{Atm::FINAL};

  /** Type of line quantity */
  Line line{Line::FINAL};

  /** Type of sensor quantity */
  Sensor sensor{Sensor::FINAL};

  /** Type of special quantity */
  Special special{Special::FINAL};

  /** Perturbations for methods where theoretical computations are impossible or plain slow */
  Numeric perturbation{std::numeric_limits<Numeric>::quiet_NaN()};

  /** ID for the Line types of partial derivatives */
  QuantumIdentifier qid{};

  /** ID for some of the Special types of partial derivatives */
  ArrayOfSpeciesTag species_array_id{0};

  /** ID for some of the Special types of partial derivatives */
  String string_id{};

  /** Species ID for line parameters */
  Species::Species species_id{Species::Species::FINAL};

  /** Atmospheric type */
  explicit Target(Atm atype) : type(Type::Atm), atm(atype) {
    ARTS_ASSERT(good_enum(atype))
  }

  /** Line type */
  explicit Target(Line ltype,
                  const QuantumIdentifier& iqid,
                  Species::Species specid)
      : type(Type::Line), line(ltype), qid(iqid), species_id(specid) {
    ARTS_ASSERT(good_enum(ltype))
  }

  /** Sensor type */
  explicit Target(Sensor stype) : type(Type::Sensor), sensor(stype) {
    ARTS_ASSERT(good_enum(stype))
  }

  /** Special type */
  explicit Target(Special stype, ArrayOfSpeciesTag aostid)
      : type(Type::Special), special(stype), species_array_id(std::move(aostid)) {
    ARTS_ASSERT(stype == Special::ArrayOfSpeciesTagVMR,
                "Only for Special::ArrayOfSpeciesTagVMR, but you fed: ",
                special)
  }

  /** Special type */
  explicit Target(Special stype, String sid)
      : type(Type::Special), special(stype), string_id(std::move(sid)) {
    ARTS_ASSERT(
        stype == Special::SurfaceString or stype == Special::ScatteringString,
        "Only for Special::ScatteringString or Special::SurfaceString, but you fed: ",
        special)
  }

  /** A defaultable none-type */
  explicit Target() = default;

  /** Checks if the type of jacobian is the input atmospheric parameter */
  bool operator==(Atm other) const noexcept { return other == atm; }

  /** Checks if the type of jacobian is the input line parameter */
  bool operator==(Line other) const noexcept { return other == line; }

  /** Checks if the type of jacobian is the input sensor parameter */
  bool operator==(Sensor other) const noexcept { return other == sensor; }

  /** Checks if the type of jacobian is the input sensor parameter */
  bool operator==(Special other) const noexcept { return other == special; }

  /** Checks if the type is correct */
  bool operator==(Type other) const noexcept { return other == type; }

  [[nodiscard]] bool sameTargetType(const Target& other) const noexcept {
    return type == other.type and atm == other.atm and
           line == other.line and sensor == other.sensor and
           special == other.special;
  }

  /** Return type as string */
  [[nodiscard]] std::string_view TargetType() const noexcept {
    return toString(type);
  }

  /** Sets target based on a string */
  void TargetType(const std::string_view& s) noexcept { type = toType(s); }

  /** Sets sub target based on a string */
  void TargetSubType(const std::string_view& s) noexcept {
    atm = Atm::FINAL;
    line = Line::FINAL;
    sensor = Sensor::FINAL;
    special = Special::FINAL;

    switch (type) {
      case Type::Atm:
        atm = toAtm(s);
        break;
      case Type::Line:
        line = toLine(s);
        break;
      case Type::Sensor:
        sensor = toSensor(s);
        break;
      case Type::Special:
        special = toSpecial(s);
        break;
      case Type::FINAL: { /* leave last, don't use default */
      }
    }
  }

  [[nodiscard]] std::string_view TargetSubType() const noexcept {
    switch (type) {
      case Type::Atm:
        return toString(atm);
      case Type::Line:
        return toString(line);
      case Type::Sensor:
        return toString(sensor);
      case Type::Special:
        return toString(special);
      case Type::FINAL: { /* leave last, don't use default */
      }
    }
    return "BAD SUBTYPE";
  }

  /** Are we good? */
  [[nodiscard]] bool TargetTypeOK() const noexcept { return good_enum(type); }

  /** Are we good? */
  [[nodiscard]] bool TargetSubTypeOK() const noexcept {
    // We can only hold one valid enum at a time, and it must be of the correct type
    if (1 == (good_enum(special) + good_enum(sensor) + good_enum(line) +
              good_enum(atm))) {
      switch (type) {
        case Type::Special:
          return good_enum(special);
        case Type::Sensor:
          return good_enum(sensor);
        case Type::Line:
          return good_enum(line);
        case Type::Atm:
          return good_enum(atm);
        case Type::FINAL: { /* leave last, don't use default */
        }
      }
    }

    return false;
  }

  /** Special species case */
  [[nodiscard]] bool isSpeciesVMR() const noexcept {
    return line == Line::VMR or special == Special::ArrayOfSpeciesTagVMR;
  }

  /** Special wind case */
  [[nodiscard]] bool isWind() const noexcept {
    return atm == Atm::WindMagnitude or atm == Atm::WindU or
           atm == Atm::WindV or atm == Atm::WindW;
  }

  /** Special magnetic field case */
  [[nodiscard]] bool isMagnetic() const noexcept {
    return atm == Atm::MagneticMagnitude or atm == Atm::MagneticU or
           atm == Atm::MagneticV or atm == Atm::MagneticW;
  }

  /** Special frequency case */
  [[nodiscard]] bool isFrequency() const noexcept {
    return sensor == Sensor::FrequencyStretch or
           sensor == Sensor::FrequencyShift;
  }

  /** Special pointing case */
  [[nodiscard]] bool isPointing() const noexcept {
    return sensor == Sensor::PointingZenithInterp or
           sensor == Sensor::PointingZenithRecalc;
  }

  /** Does this type need the QuantumIdentifier? */
  [[nodiscard]] bool needQuantumIdentity() const noexcept {
    return type == Type::Line;
  }

  /** Does this type need the ArrayOfSpeciesTag? */
  [[nodiscard]] bool needArrayOfSpeciesTag() const noexcept {
    return special == Special::ArrayOfSpeciesTagVMR;
  }

  /** Does this type need the String? */
  [[nodiscard]] bool needString() const noexcept {
    return special == Special::ScatteringString or
           special == Special::SurfaceString;
  }

  friend std::ostream& operator<<(std::ostream& os, const Target& x);
};  // Target
} // namespace Jacobian
using JacobianTarget = Jacobian::Target;
using ArrayOfJacobianTarget = Array<Jacobian::Target>;


/** Deals with internal derivatives, Jacobian definition, and OEM calculations */
class RetrievalQuantity {
 public:
  /** Default constructor. Needed by make_array. */
  RetrievalQuantity()
      : msubtag(),
        msubsubtag(),
        mmode(),
        mgrids(),
        mjac() { /* Nothing to do here. */
  }

  /** Constructor that sets the values
   * 
   * @param[in] subtag The sub-derivative
   * @param[in] subsubtag The sub-sub-derivative
   * @param[in] mode The mode of the derivative
   * @param[in] analytical Boolean analytical tag
   * @param[in] perturbation The size of the perturbation required
   * @param[in] grids The retrieval grid
   */
  RetrievalQuantity(Jacobian::Target  target,
                    String  subtag,
                    String  subsubtag,
                    String  mode,
                    const Numeric& perturbation,
                    ArrayOfVector  grids)
      : msubtag(std::move(subtag)),
        msubsubtag(std::move(subsubtag)),
        mmode(std::move(mode)),
        mgrids(std::move(grids)),
        mjac(std::move(target)) {
    mjac.perturbation = perturbation;
  }
  
  /** Returns the sub-tag
   * 
   * Subtag. Eg. for gas species: O3, ClO
   * 
   * @return A representation of the sub-tag
   */
  [[nodiscard]] const String& Subtag() const { return msubtag; }
  
  /** Sets the sub-tag
   * 
   * @param[in] st A sub-tag
   */
  void Subtag(const String& st) { msubtag = st; }
  
  /** Returns the sub-sub-tag
   * 
   * SubSubtag. Eg. for scat species fields: mass_density, mass_flux, ...
   * 
   * @return A representation of the sub-sub-tag
   */
  [[nodiscard]] const String& SubSubtag() const { return msubsubtag; }
  
  /** Sets the sub-sub-tag
   * 
   * @param[in] sst A sub-sub-tag
   */
  void SubSubtag(const String& sst) { msubsubtag = sst; }
  
  /** Returns the mode
   * 
   * Calculation mode. Eg. "abs", "rel", "vmr", "nd", "From propagation matrix". 
   * Note that the latter of these only supports "vmr" for abs species.
   * 
   * @return A representation of the mode
   */
  [[nodiscard]] const String& Mode() const { return mmode; }
  
  /** Sets the mode
   * 
   * @param[in] m A mode
   */
  void Mode(const String& m) { mmode = m; }
  
  /** Returns the grids of the retrieval
   * 
   * Grids. Definition grids for the jacobian, eg. p, lat and lon.
   * 
   * @return The grids of the retrieval
   */
  [[nodiscard]] const ArrayOfVector& Grids() const { return mgrids; }
  
  /** Sets the grids of the retrieval
   * 
   * @param[in] g The grids of the retrieval
   */
  void Grids(const ArrayOfVector& g) { mgrids = g; }

  /** Number of elements in the grids 
   * 
   * The multiplicative accumulation of grid elements
   * 
   * @return The number of elements if each grid represents a dimension
   */
  [[nodiscard]] Index nelem() const {
    Index i = 1;
    for (Index j = 0; j < mgrids.nelem(); ++j) {
      i *= mgrids[j].nelem();
    }
    return i;
  }
  
  /** Get the Jacobian Target */
  Jacobian::Target& Target() {return mjac;}
  
  /** Get the Jacobian Target */
  [[nodiscard]] const Jacobian::Target& Target() const {return mjac;}
  
  /** Set the Jacobian Target */
  void Target(const Jacobian::Target& jac) {mjac=jac;}
  
  /** Returns the identity of this Jacobian
   * 
   * QuantumIdentifier as necessary for matching line specific parameters to Jacobian grid
   * 
   * @return The identity of this Jacobian
   */
  [[nodiscard]] const QuantumIdentifier& QuantumIdentity() const {
    return mjac.qid;
  }
  
  /** Return line type */
  [[nodiscard]] Jacobian::Line LineType() const noexcept {return mjac.line;}
  
  /** Return atm type equality */
  bool operator==(Jacobian::Atm other) const noexcept {return mjac==other;}
  
  /** Return line type equality */
  bool operator==(Jacobian::Line other) const noexcept {return mjac==other;}
  
  /** Return sensor type equality */
  bool operator==(Jacobian::Sensor other) const noexcept {return mjac==other;}
  
  /** Return special type equality */
  bool operator==(Jacobian::Special other) const noexcept {return mjac==other;}
  
  /** Return special type equality */
  bool operator==(Jacobian::Type other) const noexcept {return mjac==other;}
  
  /** Return special type equality */
  bool operator==(const ArrayOfSpeciesTag& st) const noexcept {return mjac==Jacobian::Special::ArrayOfSpeciesTagVMR and mjac.species_array_id == st;}
  
  /** Sets the identity of this Jacobian
   * 
   * @param[in] qi The identity of this Jacobian
   */
  void QuantumIdentity(const QuantumIdentifier& qi) { mjac.qid = qi; }
  
  /** Returns if this is a propagation matrix type */
  [[nodiscard]] bool propmattype() const noexcept {
    return mjac == Jacobian::Type::Line or 
           mjac == Jacobian::Type::Atm or 
           mjac == Jacobian::Special::ArrayOfSpeciesTagVMR;
  }
  
  [[nodiscard]] bool is_wind() const noexcept {return mjac.isWind();}
  
  [[nodiscard]] bool is_mag() const noexcept {return mjac.isMagnetic();}
  
  /** Transformation
   * 
   * FIXMEDOC@Simon The transformations are yours to fix and document
   * FIXMEDOC@Patrick The transformations are yours to fix and document
   */
  void SetTransformationFunc(const String& s) { transformation_func = s; }
  void SetTFuncParameters(const Vector& p) { tfunc_parameters = p; }
  void SetTransformationMatrix(const Matrix& A) { transformation_matrix = A; }
  void SetOffsetVector(const Vector& b) { offset_vector = b; }
  [[nodiscard]] bool HasAffine() const { return !transformation_matrix.empty(); }
  [[nodiscard]] const String& TransformationFunc() const { return transformation_func; }
  [[nodiscard]] const Vector& TFuncParameters() const { return tfunc_parameters; }
  [[nodiscard]] const Matrix& TransformationMatrix() const { return transformation_matrix; }
  [[nodiscard]] const Vector& OffsetVector() const { return offset_vector; }

  /** Checks that all the internal variables of *this match with those of the input
   * 
   * @param[in] a Another retrieval quantity object
   * @return true if the input performs the exact same Jacobian calculations
   * @return false otherwise
   */
  [[nodiscard]] bool HasSameInternalsAs(const RetrievalQuantity& a) const {
    return a.msubtag == msubtag and
           a.msubsubtag == msubsubtag and a.mmode == mmode and
           a.mjac.sameTargetType(mjac);
  }
  
  String& SubTag() {return msubtag;}
  String& SubSubTag() {return msubsubtag;}
  String& Mode() {return mmode;}
  ArrayOfVector& Grids() {return mgrids;}
  String& TransformationFunc() {return transformation_func;}
  Vector& TFuncParameters() {return tfunc_parameters;}
  Matrix& Transformation() {return transformation_matrix;}
  Vector& Offset() {return offset_vector;}

  [[nodiscard]] const String& SubTag() const { return msubtag; }
  [[nodiscard]] const String& SubSubTag() const { return msubsubtag; }
  [[nodiscard]] const Matrix& Transformation() const { return transformation_matrix; }
  [[nodiscard]] const Vector& Offset() const { return offset_vector; }

  friend ostream& operator<<(ostream& os, const RetrievalQuantity& ot);

 private:
  String msubtag;
  String msubsubtag;
  String mmode;
  ArrayOfVector mgrids;
  Jacobian::Target mjac;

  String transformation_func;
  Vector tfunc_parameters;

  Matrix transformation_matrix;
  Vector offset_vector;
};

using ArrayOfRetrievalQuantity = Array<RetrievalQuantity>;

// A macro to loop analytical jacobian quantities
#define FOR_ANALYTICAL_JACOBIANS_DO(what_to_do)                             \
  for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {              \
    if (not(jacobian_quantities[iq] == Jacobian::Type::Sensor) and          \
        not(jacobian_quantities[iq] == Jacobian::Special::SurfaceString)) { \
      what_to_do                                                            \
    }                                                                       \
  }
// A macro to loop analytical jacobian quantities
#define FOR_ANALYTICAL_JACOBIANS_DO2(what_to_do)                  \
  for (Index iq = 0; iq < jacobian_quantities.nelem(); iq++) {    \
    if (not(jacobian_quantities[iq] == Jacobian::Type::Sensor)) { \
      what_to_do                                                  \
    }                                                             \
  }

//======================================================================
//             Index ranges and transformation functions
//======================================================================

/** Determines the index range inside x and the Jacobian for each retrieval quantity

    The ranges are given as an ArrayOfArrayOfIndex, where outermost dimension
    corresponds to retrieval quantity. The inner dimension has throughout size
    2, where element 0 is the first index and element 1 is the last index of
    the range.

    @param[out]   jis            Indices, as described above
    @param[in]    any_affine     True if at least one  quantity has affine
                                 transformation. 
    @param[in]   jqs             The WSV jacobian_quantities
    @param[in]   before_affine   Set to true to get indices without any affine
                                 transformation. Default is false.

    @author Simon Pfreundschuh and Patrick Eriksson 
    @date   2017-12-30
 */
void jac_ranges_indices(ArrayOfArrayOfIndex& jis,
                        bool& any_affine,
                        const ArrayOfRetrievalQuantity& jqs,
                        const bool& before_affine = false);

/** Applies both functional and affine transformations.
 *
 *  @param[in] jacobian As the WSV jacobian
 *  @param[in] jqs As the WSV jacobian_quantities
 *
 *  @author Simon Pfreundschuh and Patrick Eriksson 
 *  @date   2017-12-30
 */
void transform_jacobian(Matrix& jacobian,
                        const Vector x,
                        const ArrayOfRetrievalQuantity& jqs);
 
/** Handles transformations of the state vector
 *
 * Applies both functional and affine transformations.
 *
 *  @param[in,out] x    As the WSV x
 *  @param[in]     jqs  As the WSV jacobian_quantities
 *
 *  @author Simon Pfreundschuh and Patrick Eriksson 
 *  @date   2017-12-30
 */
void transform_x(Vector& x, const ArrayOfRetrievalQuantity& jqs);
 
/** Handles back-transformations of the state vector
 *
 * Applies both functional and affine transformations.
 *
 *  @param[in] x As the WSV x
 *  @param[in] jqs As the WSV jacobian_quantities
 *
 *  @author Simon Pfreundschuh and Patrick Eriksson 
 *  @date   2017-12-30
 */
void transform_x_back(Vector& x_t,
                      const ArrayOfRetrievalQuantity& jqs,
                      bool revert_functional_transforms = true);

//======================================================================
//             Functions related to calculation of Jacobian
//======================================================================

/** Check that the retrieval grids are defined for each atmosphere dim

   Use this version for atmospheric fields.

   This function checks for the given atmosphere dimension that 
     I)  the retrieval grids are defined 
     II) and that they are covered by the corresponding atmospheric grid. 
   If not the return is false and an output string is created to print 
   the error to the user. If the grids are ok they are stored in an array 
   and true is  returned.
   
   @param[out] grids        The array of retrieval grids.
   @param[in,out] os        The output string stream.
   @param[in] p_grid        The atmospheric pressure grid
   @param[in] lat_grid      The atmospheric latitude grid
   @param[in] lon_grid      The atmospheric longitude grid
   @param[in] p_retr        The pressure retrieval grid.
   @param[in] lat_retr      The latitude retrieval grid.
   @param[in] lon_retr      The longitude retrieval grid.
   @param[in] p_retr_name   The control file name used for the pressure retrieval grid.
   @param[in] lat_retr_name The control file name for the latitude retrieval grid.
   @param[in] lon_retr_name The control file name for the longitude retrieval grid.
   @param[in] dim           The atmosphere dimension

   @return              Boolean for check.
   
   @author Mattias Ekstrom
   @date   2005-05-11
 */
bool check_retrieval_grids(ArrayOfVector& grids,
                           ostringstream& os,
                           const Vector& p_grid,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Vector& p_retr,
                           const Vector& lat_retr,
                           const Vector& lon_retr,
                           const String& p_retr_name,
                           const String& lat_retr_name,
                           const String& lon_retr_name,
                           const Index& dim);


/** Check that the retrieval grids are defined for each atmosphere dim

   Use this version for surface variables

   This function checks for the given atmosphere dimension that 
     I)  the retrieval grids are defined 
     II) and that they are covered by the corresponding atmospheric grid. 
   If not the return is false and an output string is created to print 
   the error to the user. If the grids are ok they are stored in an array 
   and true is  returned.
   
   @param[out] grids        The array of retrieval grids.
   @param[in,out] os        The output string stream.
   @param[in] lat_grid      The atmospheric latitude grid
   @param[in] lon_grid      The atmospheric longitude grid
   @param[in] lat_retr      The latitude retrieval grid.
   @param[in] lon_retr      The longitude retrieval grid.
   @param[in] lat_retr_name The control file name for the latitude retrieval grid.
   @param[in] lon_retr_name The control file name for the longitude retrieval grid.
   @param[in] dim           The atmosphere dimension
   @return              Boolean for check.
   
   @author Mattias Ekstrom
   @date   2005-05-11
 */
bool check_retrieval_grids(ArrayOfVector& grids,
                           ostringstream& os,
                           const Vector& lat_grid,
                           const Vector& lon_grid,
                           const Vector& lat_retr,
                           const Vector& lon_retr,
                           const String& lat_retr_name,
                           const String& lon_retr_name,
                           const Index& dim);

/** Maps jacobian data for points along the propagation path, to
    jacobian retrieval grid data.

    @param[out]  diy_dx              One element of the WSV *diy_dx*.
    @param[in]   jacobian_quantity   One element of of the WSV *jacobian_quantities*.
    @param[in]   diy_dpath           Jacobians along the propagation path.
    @param[in]   atmosphere_dim      As the WSV.
    @param[in]   ppath               As the WSV.
    @param[in]   ppath_p             The pressure at each ppath point.

    @author Patrick Eriksson 
    @date   2009-10-08
 */
void diy_from_path_to_rgrids(Tensor3View diy_dx,
                             const RetrievalQuantity& jacobian_quantity,
                             ConstTensor3View diy_dpath,
                             const Index& atmosphere_dim,
                             const Ppath& ppath,
                             ConstVectorView ppath_p);

/** diy_from_pos_to_rgrids

    Maps jacobian data for a surface position, to jacobian retrieval grid data.

    @param[out]   diy_dx             One element of the WSV *diy_dx*.
    @param[in]   jacobian_quantity   One element of of the WSV *jacobian_quantities*.
    @param[in]   diy_dpos            Jacobian for the position itself.
    @param[in]   atmosphere_dim      As the WSV.
    @param[in]   rtp_pos             As the WSV.

    @author Patrick Eriksson 
    @date   2018-04-10
 */
void diy_from_pos_to_rgrids(Tensor3View diy_dx,
                            const RetrievalQuantity& jacobian_quantity,
                            ConstMatrixView diy_dpos,
                            const Index& atmosphere_dim,
                            ConstVectorView rtp_pos);

/** Help function for analytical jacobian calculations
 * 
 *  The size of the by-path computations of the Jacobian
 *  is set and zeroed
 * 
 *  @param[in]    jacobian_quantities   As the WSV.
 *  @param[in]    np                    The path grid count
 *  @param[in]    nf                    The frequency grid count
 *  @param[in]    ns                    The Stokes grid count
 *  @param[in]    active                If true, the middle dimension is np times nf, otherwise it is nf
 * 
 *  @return       diy_dpath             As expected by the methods using diy_dpath (size: jacobian_quantities.nelem())
 * 
 *  @author Richard Larsson 
 *  @date   2020-11-12
 */
ArrayOfTensor3 get_standard_diy_dpath(const ArrayOfRetrievalQuantity& jacobian_quantities, Index np, Index nf, Index ns, bool active);

/** Help function for analytical jacobian calculations
 * 
 *  The size of the computations of the Jacobian is set and zeroed
 * 
 *  @param[in]    jacobian_quantities   As the WSV.
 *  @param[in]    np                    The path grid count
 *  @param[in]    nf                    The frequency grid count
 *  @param[in]    ns                    The Stokes grid count
 *  @param[in]    active                If true, the middle dimension is np times nf, otherwise it is nf
 * 
 *  @return       diy_dx                As expected by the methods using diy_dpath (size: jacobian_quantities.nelem())
 * 
 *  @author Richard Larsson 
 *  @date   2020-11-12
 */
ArrayOfTensor3 get_standard_starting_diy_dx(const ArrayOfRetrievalQuantity& jacobian_quantities, Index np, Index nf, Index ns, bool active);
 
/** Help function for analytical jacobian calculations

    The function determines which terms in jacobian_quantities that are 
    analytical absorption species. 

    The output Array has the position in abs_species for
    Jacobian::Line::VMR and Jacobian::Special::ArrayOfSpeciesTagVMR
    quantities.  It has -9999 for Jacobian::Atm::Particulates and
    for Jacobian::Atm::Electrons.  It has -1 for all other types of
    partial derivatives

    @param[in]    jacobian_quantities   As the WSV.
    @param[in]    abs_species           As the WSV.
    
    @return       ArrayOfIndex          With information as above (size: jacobian_quantities.nelem())

    @author Richard Larsson 
    @date   2020-11-12
*/
ArrayOfIndex get_pointers_for_analytical_species(const ArrayOfRetrievalQuantity& jacobian_quantities,
                                                 const ArrayOfArrayOfSpeciesTag& abs_species);

/** Help function for analytical jacobian calculations

    The function determines which terms in jacobian_quantities that are 
    analytical absorption species. 

    If the scat species is there, it will have an index pointing to it,
    otherwise the index will be -1
    
    cloudbox_on must be true or all output is -1

    @param[in]    jacobian_quantities   As the WSV.
    @param[in]    scat_species          As the WSV.
    @param[in]    cloudbox_on           As the WSV.
    
    @return       ArrayOfIndex          With information as above (size: jacobian_quantities.nelem())

    @author Richard Larsson 
    @date   2020-11-12
 */
ArrayOfIndex get_pointers_for_scat_species(const ArrayOfRetrievalQuantity& jacobian_quantities,
                                           const ArrayOfString& scat_species,
                                           const bool cloudbox_on);

/** Checks if analytical calculations are needed at all
    
    The template argument is either 1 or 2 for checks with the macros
    FOR_ANALYTICAL_JACOBIANS_DO or FOR_ANALYTICAL_JACOBIANS_DO2,
    respectively
    
    @param[in]    jacobian_quantities   As the WSV.
    
    @return       true if we need analytical calculations
    @return       false if we do not need analytical calculations

    @author Richard Larsson 
    @date   2020-11-12
*/
template <std::size_t N>
Index do_analytical_jacobian(const ArrayOfRetrievalQuantity& jacobian_quantities) {
  static_assert(N == 1 or N == 2, "FOR_ANALYTICAL_JACOBIANS_DO or FOR_ANALYTICAL_JACOBIANS_DO2");
  if constexpr (N == 1) FOR_ANALYTICAL_JACOBIANS_DO(return 1;)
  else if constexpr (N == 2) FOR_ANALYTICAL_JACOBIANS_DO2(return 1;)
  return 0;
}

/** Adopts grid positions to extrapolation used for jacobians

  The standard interpolation scheme applies a linear extrapolation, while for
  the jacobians the extrapolation can be seen as a "closest" interpolation.
  That is, for points outisde the covered grid, the value at closest end point
  is taken. And here extrapolation to +-Inf is allowed.

  This function modifies grid positions to the jacobian extrapolation approach.
  For efficiency, the input grid positions are not asserted at all, and
  "extrapolation points" are identified simply  by a fd outside [0,1].

  @param[in/out] gp   Array of grid positions.

  @author Patrick Eriksson 
  @date   2015-09-10
 */
void jacobian_type_extrapol(ArrayOfGridPos& gp);

/** Calculates polynomial basis functions

   The basis function is b(x) = 1 for poly_coeff = 0. For higher
   coefficients, x^poly_coeff - m, where first the range covered by
   *x* is normalised to [-1,1] and m is selected in such way that
   sum(b) = 0.
   
   @param[out] b            Calculated basis function.
   @param[in] x            The grid over which the fit shall be performed.
   @param[in] poly_coeff   Polynomial coefficient.
   
   @author Patrick Eriksson
   @date   2008-11-07
 */
void polynomial_basis_func(Vector& b, const Vector& x, const Index& poly_coeff);

/** Calculate baseline fit
 *
 * Computes the baseline fit from a given state vector.
 *
 * Given a retrieval quantitiy which is either a polynomial or a sine baseline fit
 * this function computes the baseline offset in y_baseline.
 *
 *  @param[out] y_baseline The computed baseline offset. Computed baseline offset are
 *                         accumulated into this vector, so it must be
 *                         initialized externally!
 *  @param[in] x State vector consisten with given retrieval quantity
 *  @param[in] mblock_index The index of the measurement block.
 *  @param[in] sensor_response Must be consistent with size of y_baseline.
 *  @param[in] sensor_response_pol_grid Must be consistent with size of y_baseline.
 *  @param[in] sensor_response_f_grid Must be consistent with size of y_baseline.
 *  @param[in] sensor_dlos_grid Must be consistent with size of y_baseline.
 *  @param[in] rq The poly- or sinefit retrieval quantity
 *  @param[in] rq_index The index of the retrieval quantity
 *  @param[in] jacobian_indices
 */
void calcBaselineFit(Vector& y_baseline,
                     const Vector& x,
                     const Index& mblock_index,
                     const Sparse& sensor_response,
                     const ArrayOfIndex& sensor_response_pol_grid,
                     const Vector& sensor_response_f_grid,
                     const Matrix& sensor_response_dlos_grid,
                     const RetrievalQuantity& rq,
                     const Index rq_index,
                     const ArrayOfArrayOfIndex& jacobian_indices);

/** Scale factor for conversion between gas species units.

    The function finds the factor with which the total absorption of a
    gas species shall be multiplicated to match the selected
    (jacobian) unit. 

    @param[out]   x     Scale factor
    @param[in]   unit   Unit selected.
    @param[in]   vmr    VMR value.
    @param[in]   p      Pressure
    @param[in]   t      Temperature.

    @author Patrick Eriksson 
    @date   2009-10-08
 */
void vmrunitscf(Numeric& x,
                const String& unit,
                const Numeric& vmr,
                const Numeric& p,
                const Numeric& t);

/** Scale factor for conversion of derivatives with respect to VMR.

    The function finds the factor with which a partial derivative with respect
    to gas species VMR shall be multiplicated to match the selected (jacobian)
    unit. The function was implemented for scaling of *dpropmat_clearsky_dx*
    but it could also be of use in other contexts.

    @param[out]   x      Scale factor
    @param[in]   unit   Unit selected.
    @param[in]   vmr    VMR value.
    @param[in]   p      Pressure
    @param[in]   t      Temperature.

    @author Patrick Eriksson 
    @date   2015-12-11
 */
void dxdvmrscf(Numeric& x,
               const String& unit,
               const Numeric& vmr,
               const Numeric& p,
               const Numeric& t);

//======================================================================
//             Propmat partials descriptions
//======================================================================

/** Returns the temperature perturbation if it exists
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return Temperature perturbation
 */
Numeric temperature_perturbation(const ArrayOfRetrievalQuantity& js) noexcept;

/** Returns the frequency perturbation if it exists
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return Frequency perturbation
 */
Numeric frequency_perturbation(const ArrayOfRetrievalQuantity& js) noexcept;

/** Returns the magnetic field perturbation if it exists
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return Magnetic field perturbation
 */
Numeric magnetic_field_perturbation(
    const ArrayOfRetrievalQuantity& js) noexcept;

/** Returns a string of the retrieval quantity propagation matrix type
 * 
 * Only use for debugging purpose
 * 
 * @param[in] rq A retrieval quantity
 * @return Text representation
 */
String propmattype_string(const RetrievalQuantity& rq);

//======================================================================
//             Propmat partials boolean functions
//======================================================================

/** Returns if the Retrieval quantity is a wind parameter in propagation matrix calculations
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_wind_parameter(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a frequency parameter in propagation matrix calculations
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_frequency_parameter(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a derived magnetic parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_derived_magnetic_parameter(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a NLTE parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_nlte_parameter(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a G0 derivative
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_G0(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a D0 derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_D0(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a G0 derivative
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_G2(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a D2 derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_D2(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a FVC derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_FVC(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a ETA derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_ETA(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a Y derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_Y(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a G derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_G(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a DV derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_pressure_broadening_DV(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a X0 derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_lineshape_parameter_X0(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a X1 derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_lineshape_parameter_X1(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a X2 derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_lineshape_parameter_X2(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a G0, D0, G2, D2, ETA, FVC derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_lineshape_parameter_bar_linemixing(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is a G0, D0, G2, D2, ETA, FVC, Y, G, DV derivative
 * 
 * See LineShape::Model for details on parameter
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_lineshape_parameter(const RetrievalQuantity& t) noexcept;

/** Returns if the Retrieval quantity is related to the absorption line
 * 
 * @param[in] t A retrieval quantity
 * @return true if it is
 * @return false if it is not
 */
bool is_line_parameter(const RetrievalQuantity& t) noexcept;

/** Returns if the array supports CIA derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_CIA(const ArrayOfRetrievalQuantity& js);

/** Returns if the array supports HITRAN cross-section derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_hitran_xsec(const ArrayOfRetrievalQuantity& js);


/** Returns if the array supports continuum derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_continuum(const ArrayOfRetrievalQuantity& js);

/** Returns if the array supports relaxation matrix derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_relaxation_matrix(const ArrayOfRetrievalQuantity& js);

/** Returns if the array supports lookup table derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_lookup(const ArrayOfRetrievalQuantity& js);

/** Returns if the array supports Zeeman derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_zeeman(const ArrayOfRetrievalQuantity& js);

/** Returns if the array supports Faraday derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_faraday(const ArrayOfRetrievalQuantity& js);

/** Returns if the array supports propagation matrix derivatives
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool supports_propmat_clearsky(const ArrayOfRetrievalQuantity& js);

/** Returns if the Retrieval quantity is VMR derivative for all the species in the species tags
 * 
 * Very slow compared to index input
 * 
 * @param[in] rq A retrieval quantity
 * @param[in] ast A list of species tags
 * @return true if all species in the tags list matches the species in the retrieval quantity
 * @return false otherwise
 */
bool species_match(const RetrievalQuantity& rq, const ArrayOfSpeciesTag& ast);

/** Returns if the Retrieval quantity is VMR derivative for all the species in the species tags
 * 
 * @param[in] rq A retrieval quantity
 * @param[in] species An index-mapped species
 * @return true the species match the species in the retrieval quantity
 * @return false otherwise
 */
bool species_match(const RetrievalQuantity& rq, const Species::Species species);

/** Returns if the Retrieval quantity is VMR derivative for all the species in the species tags
 * 
 * @param[in] rq A retrieval quantity
 * @param[in] species An index-mapped species
 * @param[in] iso An index-mapped isotopologue
 * @return true the species and isotopologue match the species and isotopologue in the retrieval quantity
 * @return false otherwise
 */
bool species_iso_match(const RetrievalQuantity& rq,
                       const Species::IsotopeRecord& ir);

/** Returns if the array wants the temperature derivative
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool do_temperature_jacobian(const ArrayOfRetrievalQuantity& js) noexcept;

/** Deals with whether or not we should do a VMR derivative */
struct jacobianVMRcheck {
  bool test;
  const QuantumIdentifier& qid;
};

/** Returns the required info for VMR Jacobian
 * 
 * FIXME: The entire existence of this function is a logical error of the programmer...
 * 
 * @param[in] js As jacobian_quantities WSV
 * @param[in] line_qid A line identifier
 * @return true and retrieval quantity quantum identity pointer if available
 * @return false and line_qid pointer if not available
 */
jacobianVMRcheck do_vmr_jacobian(const ArrayOfRetrievalQuantity& js,
                                 const QuantumIdentifier& line_qid) noexcept;

/** Returns if the array wants a line center derivative
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool do_line_center_jacobian(const ArrayOfRetrievalQuantity& js) noexcept;

/** Returns if the array wants a wind-based frequency derivative derivative
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool do_wind_jacobian(const ArrayOfRetrievalQuantity& js) noexcept;

/** Returns if the array wants a frequency derivative
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool do_frequency_jacobian(const ArrayOfRetrievalQuantity& js) noexcept;

/** Returns if the array wants a magnetic derivative
 * 
 * @param[in] js As jacobian_quantities WSV
 * @return true if it does
 * @return false if it does not
 */
bool do_magnetic_jacobian(const ArrayOfRetrievalQuantity& js) noexcept;

#endif  // jacobian_h

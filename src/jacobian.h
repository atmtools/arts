/* Copyright (C) 2004-2012 Mattias Ekstrom <ekstrom@rss.chalmers.se>

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

/** 
  @file   jacobian.h
  @author Mattias Ekstrom <ekstrom@rss.chalmers.se>
  @date   2004-09-14

  @brief  Routines for setting up the jacobian.
 */

#ifndef jacobian_h
#define jacobian_h

#include <iostream>
#include <map>
#include <stdexcept>
#include "species_tags.h"
#include "agenda_class.h"
#include "array.h"
#include "bifstream.h"
#include "enums.h"
#include "interpolation.h"
#include "logic.h"
#include "matpackI.h"
#include "methods.h"
#include "mystring.h"
#include "ppath.h"
#include "quantum.h"  

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
          NLTE,
          SpecialParameter1,
          SpecialParameter2,
          SpecialParameter3
          )

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
class Target {
private:
  /**! Type of quantity, never set manually */
  Type mtype;
  
  /** Type of atm quantity */
  Atm matm;
  
  /** Type of line quantity */
  Line mline;
  
  /** Type of sensor quantity */
  Sensor msensor;
  
  /** Type of special quantity */
  Special mspecial;
  
  /** Perturbations for methods where theoretical computations are impossible or plain slow */
  Numeric mperturbation;
  
  /** ID for the Line types of partial derivatives */
  QuantumIdentifier mqid;
  
  /** ID for some of the Special types of partial derivatives */
  ArrayOfSpeciesTag maostid;
  
  /** ID for some of the Special types of partial derivatives */
  String msid;
  
  /** Species ID for line parameters */
  Species::Species mspecid;
public:
  /** Atmospheric type */
  explicit Target (Atm type) :
  mtype(Type::Atm), matm(type), mline(Line::FINAL),
  msensor(Sensor::FINAL), mspecial(Special::FINAL),
  mperturbation(std::numeric_limits<Numeric>::quiet_NaN()),
  mqid(), maostid(0), msid(), mspecid(Species::Species::FINAL) {
    ARTS_ASSERT(good_enum(type))
  }
  
  /** Line type */
  explicit Target (Line type, const QuantumIdentifier& qid, Species::Species specid) :
  mtype(Type::Line), matm(Atm::FINAL), mline(type),
  msensor(Sensor::FINAL), mspecial(Special::FINAL),
  mperturbation(std::numeric_limits<Numeric>::quiet_NaN()),
  mqid(qid), maostid(0), msid(), mspecid(specid) {
    ARTS_ASSERT(good_enum(type))
  }
  
  /** Sensor type */
  explicit Target (Sensor type) :
  mtype(Type::Sensor), matm(Atm::FINAL), mline(Line::FINAL),
  msensor(type), mspecial(Special::FINAL),
  mperturbation(std::numeric_limits<Numeric>::quiet_NaN()),
  mqid(), maostid(0), msid(), mspecid(Species::Species::FINAL) {
    ARTS_ASSERT(good_enum(type))
  }
  
  /** Special type */
  explicit Target (Special type, const ArrayOfSpeciesTag& aostid) :
  mtype(Type::Special), matm(Atm::FINAL), mline(Line::FINAL),
  msensor(Sensor::FINAL), mspecial(type),
  mperturbation(std::numeric_limits<Numeric>::quiet_NaN()),
  mqid(), maostid(aostid), msid(), mspecid(Species::Species::FINAL) {
    ARTS_ASSERT (type == Special::ArrayOfSpeciesTagVMR,
      "Only for Special::ArrayOfSpeciesTagVMR, but you fed: ", mspecial)
  }
  
  /** Special type */
  explicit Target (Special type, const String& sid) :
  mtype(Type::Special), matm(Atm::FINAL), mline(Line::FINAL),
  msensor(Sensor::FINAL), mspecial(type),
  mperturbation(std::numeric_limits<Numeric>::quiet_NaN()),
  mqid(), maostid(0), msid(sid), mspecid(Species::Species::FINAL) {
    ARTS_ASSERT (type == Special::SurfaceString or type == Special::ScatteringString,
      "Only for Special::ScatteringString or Special::SurfaceString, but you fed: ", mspecial)
  }
  
  /** A defaultable none-type */
  explicit Target () : 
  mtype(Type::FINAL), matm(Atm::FINAL), mline(Line::FINAL),
  msensor(Sensor::FINAL), mspecial(Special::FINAL),
  mperturbation(std::numeric_limits<Numeric>::quiet_NaN()),
  mqid(), maostid(0), msid(), mspecid(Species::Species::FINAL) {}
  
  /** Perturbation */
  void Perturbation(Numeric x) noexcept {mperturbation=x;}
  Numeric& Perturbation() noexcept {return mperturbation;}
  Numeric Perturbation() const noexcept {return mperturbation;}
  
  /** ID */
  void QuantumIdentity(const QuantumIdentifier& x) noexcept {mqid=x;}
  QuantumIdentifier& QuantumIdentity() noexcept {return mqid;}
  const QuantumIdentifier& QuantumIdentity() const noexcept {return mqid;}
  
  /** ID */
  void SpeciesList(const ArrayOfSpeciesTag& x) noexcept {maostid = x;}
  ArrayOfSpeciesTag& SpeciesList() noexcept {return maostid;}
  const ArrayOfSpeciesTag& SpeciesList() const noexcept {return maostid;}
  
  /** ID */
  void StringKey(const String& x) noexcept {msid = x;}
  String& StringKey() noexcept {return msid;}
  const String& StringKey() const noexcept {return msid;}
  
  /** Return the atm type */
  Atm AtmType() const noexcept {return matm;}
  
  /** Return the line type */
  Line LineType() const noexcept {return mline;}
  
  /** Return the sensor type */
  Sensor SensorType() const noexcept {return msensor;}
  
  /** Return the special type */
  Special SpecialType() const noexcept {return mspecial;}
  
  /** Checks if the type of jacobian is the input atmospheric parameter */
  bool operator==(Atm other) const noexcept {return other == matm;}
  
  /** Checks if the type of jacobian is the input line parameter */
  bool operator==(Line other) const noexcept {return other == mline;}
  
  /** Checks if the type of jacobian is the input sensor parameter */
  bool operator==(Sensor other) const noexcept {return other == msensor;}
  
  /** Checks if the type of jacobian is the input sensor parameter */
  bool operator==(Special other) const noexcept {return other == mspecial;}
  
  /** Checks if the type is correct */
  bool operator==(Type other) const noexcept {return other == mtype;}
  
  bool sameTargetType(const Target& other) const noexcept {
    return mtype    == other.mtype    and
           matm     == other.matm     and
           mline    == other.mline    and
           msensor  == other.msensor  and
           mspecial == other.mspecial;
  }
  
  /** Return type as string */
  std::string_view TargetType() const noexcept {return toString(mtype);}
  
  /** Sets target based on a string */
  void TargetType(const std::string_view& s) noexcept {mtype = toType(s);}
  
  /** Sets sub target based on a string */
  void TargetSubType(const std::string_view& s) noexcept {
    matm = Atm::FINAL;
    mline = Line::FINAL;
    msensor = Sensor::FINAL;
    mspecial = Special::FINAL;
    
    switch (mtype) {
      case Type::Atm: matm = toAtm(s); break;
      case Type::Line: mline = toLine(s); break;
      case Type::Sensor: msensor = toSensor(s); break;
      case Type::Special: mspecial = toSpecial(s); break;
      case Type::FINAL: { /* leave last, don't use default */ }
    }
  }
  
  std::string_view TargetSubType() const noexcept {
    switch (mtype) {
      case Type::Atm: return toString(matm);
      case Type::Line: return toString(mline);
      case Type::Sensor: return toString(msensor);
      case Type::Special: return toString(mspecial);
      case Type::FINAL: { /* leave last, don't use default */ }
    }
    return "BAD SUBTYPE";
  }
  
  /** Are we good? */
  bool TargetTypeOK() const noexcept {return good_enum(mtype);}

  /** Are we good? */
  bool TargetSubTypeOK() const noexcept {
    // We can only hold one valid enum at a time, and it must be of the correct type
    if (1 == (good_enum(mspecial) + good_enum(msensor) + good_enum(mline) + good_enum(matm))) {
      switch (mtype) {
        case Type::Special:
          return good_enum(mspecial);
        case Type::Sensor:
          return good_enum(msensor);
        case Type::Line:
          return good_enum(mline);
        case Type::Atm:
          return good_enum(matm);
        case Type::FINAL: { /* leave last, don't use default */ }
      }
    }
    
    return false;
  }
  
  /** Special species case */
  bool isSpeciesVMR() const noexcept {
    return mline == Line::VMR or mspecial == Special::ArrayOfSpeciesTagVMR;
  }
  
  /** Special wind case */
  bool isWind() const noexcept {
    return matm == Atm::WindMagnitude or 
           matm == Atm::WindU or
           matm == Atm::WindV or
           matm == Atm::WindW;
  }
  
  /** Special magnetic field case */
  bool isMagnetic() const noexcept {
    return matm == Atm::MagneticMagnitude or 
           matm == Atm::MagneticU or
           matm == Atm::MagneticV or
           matm == Atm::MagneticW;
  }
  
  /** Special frequency case */
  bool isFrequency() const noexcept {
    return msensor == Sensor::FrequencyStretch or 
           msensor == Sensor::FrequencyShift;
  }
  
  /** Special pointing case */
  bool isPointing() const noexcept {
    return msensor == Sensor::PointingZenithInterp or 
           msensor == Sensor::PointingZenithRecalc;
  }
  
  /** Does this type need the QuantumIdentifier? */
  bool needQuantumIdentity() const noexcept {
    return mtype == Type::Line;
  }
  
  /** Does this type need the ArrayOfSpeciesTag? */
  bool needArrayOfSpeciesTag() const noexcept {
    return mspecial == Special::ArrayOfSpeciesTagVMR;
  }
  
  /** Does this type need the String? */
  bool needString() const noexcept {
    return mspecial == Special::ScatteringString or
           mspecial == Special::SurfaceString;
  }
  
  const Species::Species& LineSpecies() const noexcept {return mspecid;}
  Species::Species& LineSpecies() noexcept {return mspecid;}
  void LineSpecies(Species::Species x) noexcept {mspecid = x;}
};  // Target

/** Output operator 
 *
 * @param[in,out] os A stream
 * @param[in] x A Target
 * @return A modified stream
 */
std::ostream& operator<<(std::ostream& os, const Target& x);
};  // Jacobian
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
  RetrievalQuantity(const Jacobian::Target& target,
                    const String& subtag,
                    const String& subsubtag,
                    const String& mode,
                    const Numeric& perturbation,
                    const ArrayOfVector& grids)
      : msubtag(subtag),
        msubsubtag(subsubtag),
        mmode(mode),
        mgrids(grids),
        mjac(target) {
    mjac.Perturbation() = perturbation;
  }
  
  /** Returns the sub-tag
   * 
   * Subtag. Eg. for gas species: O3, ClO
   * 
   * @return A representation of the sub-tag
   */
  const String& Subtag() const { return msubtag; }
  
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
  const String& SubSubtag() const { return msubsubtag; }
  
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
  const String& Mode() const { return mmode; }
  
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
  const ArrayOfVector& Grids() const { return mgrids; }
  
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
  Index nelem() const {
    Index i = 1;
    for (Index j = 0; j < mgrids.nelem(); ++j) {
      i *= mgrids[j].nelem();
    }
    return i;
  }
  
  /** Get the Jacobian Target */
  Jacobian::Target& Target() {return mjac;}
  
  /** Get the Jacobian Target */
  const Jacobian::Target& Target() const {return mjac;}
  
  /** Set the Jacobian Target */
  void Target(const Jacobian::Target& jac) {mjac=jac;}
  
  /** Returns the identity of this Jacobian
   * 
   * QuantumIdentifier as necessary for matching line specific parameters to Jacobian grid
   * 
   * @return The identity of this Jacobian
   */
  const QuantumIdentifier& QuantumIdentity() const {
    return mjac.QuantumIdentity();
  }
  
  /** Return line type */
  Jacobian::Line LineType() const noexcept {return mjac.LineType();}
  
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
  bool operator==(const ArrayOfSpeciesTag& st) const noexcept {return mjac==Jacobian::Special::ArrayOfSpeciesTagVMR and mjac.SpeciesList() == st;}
  
  /** Sets the identity of this Jacobian
   * 
   * @param[in] qi The identity of this Jacobian
   */
  void QuantumIdentity(const QuantumIdentifier& qi) { mjac.QuantumIdentity() = qi; }
  
  /** Returns if this is a propagation matrix type */
  bool propmattype() const noexcept {
    return mjac == Jacobian::Type::Line or 
           mjac == Jacobian::Type::Atm or 
           mjac == Jacobian::Special::ArrayOfSpeciesTagVMR;
  }
  
  bool is_wind() const noexcept {return mjac.isWind();}
  
  bool is_mag() const noexcept {return mjac.isMagnetic();}
  
  /** Transformation
   * 
   * FIXMEDOC@Simon The transformations are yours to fix and document
   * FIXMEDOC@Patrick The transformations are yours to fix and document
   */
  void SetTransformationFunc(const String& s) { transformation_func = s; }
  void SetTFuncParameters(const Vector& p) { tfunc_parameters = p; }
  void SetTransformationMatrix(const Matrix& A) { transformation_matrix = A; }
  void SetOffsetVector(const Vector& b) { offset_vector = b; }
  bool HasAffine() const { return !transformation_matrix.empty(); }
  const String& TransformationFunc() const { return transformation_func; }
  const Vector& TFuncParameters() const { return tfunc_parameters; }
  const Matrix& TransformationMatrix() const { return transformation_matrix; }
  const Vector& OffsetVector() const { return offset_vector; }

  /** Checks that all the internal variables of *this match with those of the input
   * 
   * @param[in] a Another retrieval quantity object
   * @return true if the input performs the exact same Jacobian calculations
   * @return false otherwise
   */
  bool HasSameInternalsAs(const RetrievalQuantity& a) const {
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

/** Output operator for RetrievalQuantity */
ostream& operator<<(ostream& os, const RetrievalQuantity& ot);

typedef Array<RetrievalQuantity> ArrayOfRetrievalQuantity;

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
 *  @return       diy_dpath             As expected by the methods using diy_dpath (size: jacobain_quantities.nelem())
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
 *  @return       diy_dx                As expected by the methods using diy_dpath (size: jacobain_quantities.nelem())
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
    
    @return       ArrayOfIndex          With information as above (size: jacobain_quantities.nelem())

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
    
    @return       ArrayOfIndex          With information as above (size: jacobain_quantities.nelem())

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

  This function modifies grid positions to jacobaina extrapolation approach.
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

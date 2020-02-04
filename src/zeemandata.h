/* Copyright (C) 2018 Richard Larsson

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
 * @file   zeemandata.h
 * @author Richard Larsson <larsson (at) mps.mpg.de>
 * @date   2018-04-06
 * 
 * @brief Headers and class definition of Zeeman modeling
 * 
 * This file serves to describe the Zeeman splitting
 * implementations using various up-to-speed methods
 */

#ifndef zeemandata_h
#define zeemandata_h

#include <limits>
#include "constants.h"
#include "file.h"
#include "mystring.h"
#include "quantum.h"
#include "wigner_functions.h"

/** Implements Zeeman modeling */
namespace Zeeman {
  
/** Zeeman polarization selection */
enum class Polarization { SigmaMinus, Pi, SigmaPlus };

/** Gives the change of M given a polarization type
 * 
 * @param[in] type The polarization type
 * 
 * @return The change in M
 */
constexpr Index dM(Polarization type) noexcept {
  switch (type) {
    case Polarization::SigmaMinus:
      return -1;
    case Polarization::Pi:
      return 0;
    case Polarization::SigmaPlus:
      return 1;
  }
  return std::numeric_limits<Index>::max();
}

/** Gives the lowest M for a polarization type of this transition
 * 
 * Since the polarization determines the change in M, this
 * function gives the first M of interest in the range of M
 * possible for a given transition
 * 
 * The user has to ensure that Ju and Jl is a valid transition
 * 
 * @param[in] Ju J of the upper state
 * @param[in] Jl J of the upper state
 * @param[in] type The polarization type
 * 
 * @return The lowest M value
 */
constexpr Rational start(Rational Ju, Rational Jl, Polarization type) noexcept {
  switch (type) {
    case Polarization::SigmaMinus:
      if (Ju < Jl)
        return -Ju;
      else if (Ju == Jl)
        return -Ju + 1;
      else
        return -Ju + 2;
    case Polarization::Pi:
      return -std::min(Ju, Jl);
    case Polarization::SigmaPlus:
      return -Ju;
  }
  return Rational(std::numeric_limits<Index>::max());
}

/** Gives the largest M for a polarization type of this transition
 * 
 * Since the polarization determines the change in M, this
 * function gives the last M of interest in the range of M
 * possible for a given transition
 * 
 * The user has to ensure that Ju and Jl is a valid transition
 * 
 * @param[in] Ju J of the upper state
 * @param[in] Jl J of the upper state
 * @param[in] type The polarization type
 * 
 * @return The largest M value
 */
constexpr Rational end(Rational Ju, Rational Jl, Polarization type) noexcept {
  switch (type) {
    case Polarization::SigmaMinus:
      return Ju + 1;
    case Polarization::Pi:
      return std::min(Ju, Jl);
    case Polarization::SigmaPlus:
      if (Ju < Jl)
        return Ju + 1;
      else if (Ju == Jl)
        return Ju;
      else
        return Jl;
  }
  return Rational(std::numeric_limits<Index>::max());
}

/** Gives the number of elements of the polarization type of this transition
 * 
 * The user has to ensure that Ju and Jl is a valid transition
 * 
 * @param[in] Ju J of the upper state
 * @param[in] Jl J of the upper state
 * @param[in] type The polarization type
 * 
 * @return The number of elements
 */
constexpr Index nelem(Rational Ju, Rational Jl, Polarization type) noexcept {
  return (end(Ju, Jl, type) - start(Ju, Jl, type)).toIndex() + 1;
}

/** Gives the upper state M value at an index
 * 
 * The user has to ensure that Ju and Jl is a valid transition
 * 
 * The user has to ensure n is less than the number of elements
 * 
 * @param[in] Ju J of the upper state
 * @param[in] Jl J of the upper state
 * @param[in] type The polarization type
 * @param[in] n The position
 * 
 * @return The upper state M
 */
constexpr Rational Mu(Rational Ju,
                      Rational Jl,
                      Polarization type,
                      Index n) noexcept {
  return start(Ju, Jl, type) + n;
}


/** Gives the lower state M value at an index
 * 
 * The user has to ensure that Ju and Jl is a valid transition
 * 
 * The user has to ensure n is less than the number of elements
 * 
 * @param[in] Ju J of the upper state
 * @param[in] Jl J of the upper state
 * @param[in] type The polarization type
 * @param[in] n The position
 * 
 * @return The lower state M
 */
constexpr Rational Ml(Rational Ju,
                      Rational Jl,
                      Polarization type,
                      Index n) noexcept {
  return Mu(Ju, Jl, type, n) + dM(type);
}

/** The renormalization factor of a polarization type
 * 
 * The polarization comes from some geometry.  This function
 * returns the factor we need to compute that geometry and to
 * turn it into something that normalizes every possible M
 * for this type into some strength that sums to unity
 *  
 * @param[in] type The polarization type
 * 
 * @return Rescale factor
 */
constexpr Numeric PolarizationFactor(Polarization type) noexcept {
  switch (type) {
    case Polarization::SigmaMinus:
      return .75;
    case Polarization::Pi:
      return 1.5;
    case Polarization::SigmaPlus:
      return .75;
  }
  return std::numeric_limits<Numeric>::max();
}

/** Checks if the quantum numbers are good for this transition
 * 
 * Given some Hund state, various quantum numbers must
 * be defined to allow the Zeeman calculations to work
 *  
 * @param[in] qns Quantum numbers of a level
 * 
 * @return If the numbers can be used to compute simple Zeeman effect
 */
constexpr bool GoodHundData(const QuantumNumbers& qns) noexcept {
  if (qns[QuantumNumberType::Hund].isUndefined()) return false;
  switch (Hund(qns[QuantumNumberType::Hund].toIndex())) {
    case Hund::CaseA:
      if (qns[QuantumNumberType::Omega].isUndefined() or
          qns[QuantumNumberType::J].isUndefined() or
          qns[QuantumNumberType::Lambda].isUndefined() or
          qns[QuantumNumberType::S].isUndefined())
        return false;
      break;
    case Hund::CaseB:
      if (qns[QuantumNumberType::N].isUndefined() or
          qns[QuantumNumberType::J].isUndefined() or
          qns[QuantumNumberType::Lambda].isUndefined() or
          qns[QuantumNumberType::S].isUndefined())
        return false;
      break;
    default:
      return false;
  }
  return true;
}

/** Computes the Zeeman splitting coefficient
 * 
 * The level should be Hund case b type and all
 * the values have to be defined
 *  
 * @param[in] N The N quantum number of the level
 * @param[in] J The J quantum number of the level
 * @param[in] Lambda The Lambda quantum number of the level
 * @param[in] S The S quantum number of the level
 * @param[in] GS The spin Landé coefficient of the molecule
 * @param[in] GL The Landé coefficient of the molecule
 * 
 * @return Zeeman splitting coefficient of the level
 */
constexpr Numeric SimpleGCaseB(Rational N,
                               Rational J,
                               Rational Lambda,
                               Rational S,
                               Numeric GS,
                               Numeric GL) noexcept {
  auto JJ = J * (J + 1);
  auto NN = N * (N + 1);
  auto SS = S * (S + 1);
  auto LL = Lambda * Lambda;

  if (JJ == 0)
    return 0.0;
  else if (NN not_eq 0) {
    auto T1 = ((JJ + SS - NN) / JJ / 2).toNumeric();
    auto T2 = ((JJ - SS + NN) * LL / NN / JJ / 2).toNumeric();
    return GS * T1 + GL * T2;
  } else {
    auto T1 = ((JJ + SS - NN) / JJ / 2).toNumeric();
    return GS * T1;
  }
}

/** Computes the Zeeman splitting coefficient
 * 
 * The level should be Hund case a type and all
 * the values have to be defined
 *  
 * @param[in] Omega The Omega quantum number of the level
 * @param[in] J The J quantum number of the level
 * @param[in] Lambda The Lambda quantum number of the level
 * @param[in] Sigma The Sigma quantum number of the level
 * @param[in] GS The spin Landé coefficient of the molecule
 * @param[in] GL The Landé coefficient of the molecule
 * 
 * @return Zeeman splitting coefficient of the level
 */
constexpr Numeric SimpleGCaseA(Rational Omega,
                               Rational J,
                               Rational Lambda,
                               Rational Sigma,
                               Numeric GS,
                               Numeric GL) noexcept {
  auto JJ = J * (J + 1);

  if (JJ == 0)
    return 0.0;
  else {
    auto DIV = Omega / JJ;
    auto T1 = (Sigma * DIV).toNumeric();
    auto T2 = (Lambda * DIV).toNumeric();
    return GS * T1 + GL * T2;
  }
}

/** Computes the Zeeman splitting coefficient
 * 
 * The level should be Hund case a or b type and all
 * the quantum numbers have to be defined
 *  
 * @param[in] qns Quantum numbers of a level
 * @param[in] GS The spin Landé coefficient of the molecule
 * @param[in] GS The Landé coefficient of the molecule
 * 
 * @return If the numbers can be used to compute simple Zeeman effect
 */
constexpr Numeric SimpleG(const QuantumNumbers& qns,
                          const Numeric& GS,
                          const Numeric& GL) noexcept{
  if (not GoodHundData(qns))
    return 0;

  switch (Hund(qns[QuantumNumberType::Hund].toIndex())) {
    case Hund::CaseA:
      return SimpleGCaseA(qns[QuantumNumberType::Omega],
                          qns[QuantumNumberType::J],
                          qns[QuantumNumberType::Lambda],
                          qns[QuantumNumberType::S],
                          GS,
                          GL);
    case Hund::CaseB:
      return SimpleGCaseB(qns[QuantumNumberType::N],
                          qns[QuantumNumberType::J],
                          qns[QuantumNumberType::Lambda],
                          qns[QuantumNumberType::S],
                          GS,
                          GL);
  }

  return 0;
}

/** Main storage for Zeeman splitting coefficients
 * 
 * The splitting data has an upper (gu) and lower (gl)
 * component and this stores both of them to not confuse
 * them elsewhere
 */
struct SplittingData {
  Numeric gu, gl;
};

/** Main Zeeman Model
 * 
 * This model contains the splitting coefficients
 * of an energy level.  Various detailed and simplified
 * initialization routines are defined.  Is also used
 * as the interface for all Zeeman computations
 */
class Model {
 private:
  SplittingData mdata;

 public:
   /** Default copy/init of Model from its only private variable */
  constexpr Model(SplittingData gs = {0, 0}) noexcept : mdata(gs) {}
  
  /** Attempts to compute Zeeman input if available
   * 
   * Will first attempt advanced initialization from
   * specialized functions for special species.  If
   * this fails, will attempt simple initialization
   * from pure Hund-cases.  If this fails, will throw
   * a runtime_error.
   * 
   * @param[in] qid Transition type quantum id
   */
  Model(const QuantumIdentifier& qid) noexcept;

  /** Returns true if the Model represents no Zeeman effect */
  constexpr bool empty() const noexcept {
    return mdata.gu == mdata.gl and mdata.gu == 0;
  }

  /** Returns the upper state g */
  Numeric& gu() noexcept { return mdata.gu; }
  
  /** Returns the lower state g */
  Numeric& gl() noexcept { return mdata.gl; }
  
  /** Returns the upper state g */
  constexpr Numeric gu() const noexcept { return mdata.gu; }
  
  /** Returns the lower state g */
  constexpr Numeric gl() const noexcept { return mdata.gl; }

  /** Gives the strength of one subline of a given polarization
   * 
   * The user has to ensure that Ju and Jl is a valid transition
   * 
   * The user has to ensure n is less than the number of elements
   * 
   * @param[in] Ju J of the upper state
   * @param[in] Jl J of the upper state
   * @param[in] type The polarization type
   * @param[in] n The position
   * 
   * @return The relative strength of the Zeeman subline
   */
  Numeric Strength(Rational Ju, Rational Jl, Polarization type, Index n) const {
    using Constant::pow2;

    auto ml = Ml(Ju, Jl, type, n);
    auto mu = Mu(Ju, Jl, type, n);
    auto dm = dM(type);
    return PolarizationFactor(type) * pow2(wigner3j(Jl, 1, Ju, ml, -dm, -mu));
  }
  
  /** Gives the splitting of one subline of a given polarization
   * 
   * The user has to ensure that Ju and Jl is a valid transition
   * 
   * The user has to ensure n is less than the number of elements
   * 
   * @param[in] Ju J of the upper state
   * @param[in] Jl J of the upper state
   * @param[in] type The polarization type
   * @param[in] n The position
   * 
   * @return The splitting of the Zeeman subline
   */
  constexpr Numeric Splitting(Rational Ju, Rational Jl, Polarization type, Index n) const
      noexcept {
    using Constant::bohr_magneton;
    using Constant::h;
    constexpr Numeric C = bohr_magneton / h;

    return C * (Ml(Ju, Jl, type, n).toNumeric() * gl() -
                Mu(Ju, Jl, type, n).toNumeric() * gu());
  }

  /** Output operator for Zeeman::Model */
  friend inline std::ostream& operator<<(std::ostream& os, const Model& m);
  
  /** Input operator for Zeeman::Model */
  friend inline std::istream& operator>>(std::istream& is, Model& m);
  
  /** Output operator for Zeeman::Model */
  friend inline std::ostream& operator<<(bofstream& bof, const Model& m);
  
  /** Input operator for Zeeman::Model */
  friend inline std::istream& operator>>(bifstream& bif, Model& m);
};  // Model;

/** Returns a simple Zeeman model 
 * 
 * Will use the simple Hund case provided
 * by input.  Throws if the input is bad
 * 
 * @param[in] qid Transition type quantum id
 * 
 * @return Zeeman model
 */
Model GetSimpleModel(const QuantumIdentifier& qid) noexcept;

/** Returns an advanced Zeeman model 
 * 
 * Will look at available Quantum numbers
 * and use best approximation for the model
 * to use.  If no good approximation is available,
 * it returns Model({0, 0}).
 * 
 * @param[in] qid Transition type quantum id
 * 
 * @return Zeeman model
 */
Model GetAdvancedModel(const QuantumIdentifier& qid) noexcept;

inline std::ostream& operator<<(std::ostream& os, const Model& m) {
  os << m.mdata.gu << ' ' << m.mdata.gl;
  return os;
}

inline std::istream& operator>>(std::istream& is, Model& m) {
  is >> double_imanip() >> m.mdata.gu >> m.mdata.gl;
  return is;
}

inline std::ostream& operator<<(bofstream& bof, const Model& m) {
  bof << m.mdata.gu << m.mdata.gl;
  return bof;
}

inline std::istream& operator>>(bifstream& bif, Model& m) {
  bif >> m.mdata.gu >> m.mdata.gl;
  return bif;
}

/** Polarization vector for Zeeman Propagation Matrix
 * 
 * Meant to contain the polarization state in two vectors
 * representing [a,b,c,d] and [u,v,w] of PropagationMatrix
 * class
 */
class PolarizationVector {
 private:
  Eigen::RowVector4d att;  // attenuation vector
  Eigen::RowVector3d dis;  // dispersion vector

 public:
  /** Default init of class */
  PolarizationVector(Numeric a = 1,
                     Numeric b = 0,
                     Numeric c = 0,
                     Numeric d = 0,
                     Numeric u = 0,
                     Numeric v = 0,
                     Numeric w = 0) noexcept
      : att(a, b, c, d), dis(u, v, w){};

  /** Returns the attenuation vector */
  const Eigen::RowVector4d& attenuation() const noexcept { return att; }
  
  /** Returns the dispersion vector */
  const Eigen::RowVector3d& dispersion() const noexcept { return dis; }
  
  /** Returns the attenuation vector */
  Eigen::RowVector4d& attenuation() noexcept { return att; }
  
  /** Returns the dispersion vector */
  Eigen::RowVector3d& dispersion() noexcept { return dis; }

  /** Returns the true propagation matrix
   * 
   * Use only for debug printing if possible
   */
  Eigen::Matrix4d matrix() const noexcept {
    return (Eigen::Matrix4d() << att[0],
            att[1],
            att[2],
            att[3],
            att[1],
            att[0],
            dis[0],
            dis[1],
            att[2],
            -dis[0],
            att[0],
            dis[2],
            att[3],
            -dis[1],
            -dis[2],
            att[0])
        .finished();
  }
};

/** PolarizationVector for each Polarization
 * 
 * Contains the polarization vectors for each
 * possible polarization
 */
struct AllPolarizationVectors {
  PolarizationVector sm, pi, sp;
};

/** Computes the polarization of each polarization type
 * 
 * @param[in] theta The angle along the magnetic field
 * @param[in] eta The angle counter-clockwise in the magnetic field plane
 * 
 * @return The polarization vectors of all Zeeman polarization
 */
inline AllPolarizationVectors AllPolarization(Numeric theta,
                                              Numeric eta) noexcept {
  const Numeric ST = std::sin(theta), CT = std::cos(theta), ST2 = ST * ST,
                CT2 = CT * CT, ST2C2E = ST2 * std::cos(2 * eta),
                ST2S2E = ST2 * std::sin(2 * eta);

  AllPolarizationVectors pv;
  pv.sm = PolarizationVector(
      1 + CT2, ST2C2E, ST2S2E, 2 * CT, 4 * CT, 2 * ST2S2E, -2 * ST2C2E);
  pv.pi =
      PolarizationVector(ST2, -ST2C2E, -ST2S2E, 0, 0, -2 * ST2S2E, 2 * ST2C2E);
  pv.sp = PolarizationVector(
      1 + CT2, ST2C2E, ST2S2E, -2 * CT, -4 * CT, 2 * ST2S2E, -2 * ST2C2E);
  return pv;
}

/** The derivative of AllPolarization wrt theta
 * 
 * @param[in] theta The angle along the magnetic field
 * @param[in] eta The angle counter-clockwise in the magnetic field plane
 * 
 * @return The derivative of AllPolarization wrt theta
 */
inline AllPolarizationVectors AllPolarization_dtheta(
    Numeric theta, const Numeric eta) noexcept {
  const Numeric ST = std::sin(theta), CT = std::cos(theta),
                C2E = std::cos(2 * eta), S2E = std::sin(2 * eta), dST = CT,
                dST2 = 2 * ST * dST, dCT = -ST, dST2C2E = dST2 * C2E,
                dST2S2E = dST2 * S2E, dCT2 = 2 * CT * dCT;

  AllPolarizationVectors pv;
  pv.sm = PolarizationVector(
      dCT2, dST2C2E, dST2S2E, 2 * dCT, 4 * dCT, 2 * dST2S2E, -2 * dST2C2E);
  pv.pi = PolarizationVector(
      dST2, -dST2C2E, -dST2S2E, 0, 0, -2 * dST2S2E, 2 * dST2C2E);
  pv.sp = PolarizationVector(
      dCT2, dST2C2E, dST2S2E, -2 * dCT, -4 * dCT, 2 * dST2S2E, -2 * dST2C2E);
  return pv;
}

/** The derivative of AllPolarization wrt eta
 * 
 * @param[in] theta The angle along the magnetic field
 * @param[in] eta The angle counter-clockwise in the magnetic field plane
 * 
 * @return The derivative of AllPolarization wrt eta
 */
inline AllPolarizationVectors AllPolarization_deta(Numeric theta,
                                                   Numeric eta) noexcept {
  const Numeric ST = std::sin(theta), ST2 = ST * ST, C2E = std::cos(2 * eta),
                S2E = std::sin(2 * eta), dST2C2E = -2 * ST2 * S2E,
                dST2S2E = 2 * ST2 * C2E;

  AllPolarizationVectors pv;
  pv.sm =
      PolarizationVector(0, dST2C2E, dST2S2E, 0, 0, 2 * dST2S2E, -2 * dST2C2E);
  pv.pi = PolarizationVector(
      0, -dST2C2E, -dST2S2E, 0, 0, -2 * dST2S2E, 2 * dST2C2E);
  pv.sp =
      PolarizationVector(0, dST2C2E, dST2S2E, 0, 0, 2 * dST2S2E, -2 * dST2C2E);
  return pv;
}

/** Selects the polarization vector depending on polarization type 
 * 
 * @param[in] data The pre-computed polarization vectors
 * @param[in] type The type of polarization to select
 */
inline const PolarizationVector& SelectPolarization(
    const AllPolarizationVectors& data, Polarization type) noexcept {
      switch (type) {
    case Polarization::SigmaMinus:
      return data.sm;
    case Polarization::Pi:
      return data.pi;
    case Polarization::SigmaPlus:
      return data.sp;
  }
  std::terminate();
}

/** Contains derived values useful for Zeeman calculations
 * 
 * These are called derived since ARTS operates on a 3D grid
 * and uses the Zenith and Azimuth angles to describe LOS,
 * whereas the Zeeman effect becomes a lot clearer to work
 * with when in a defined LOS-plane.  In short, these are
 * derived 'coordinates'.
 */
struct Derived {
  Numeric H, theta, eta, dH_du, dH_dv, dH_dw, dtheta_du, dtheta_dv, dtheta_dw,
      deta_du, deta_dv, deta_dw;
};

/** Computes the derived plane from ARTS grids
 * 
 * @param[in] u Magnetic field u-parameter
 * @param[in] v Magnetic field b-parameter
 * @param[in] w Magnetic field w-parameter
 * @param[in] z Zenith angle
 * @param[in] a Azimuth angle
 * 
 * @return The derived plane
 */
inline Derived FromGrids(
    Numeric u, Numeric v, Numeric w, Numeric z, Numeric a) noexcept {
  // Constants evaluated once for both eta and theta
  auto cz = std::cos(z), ca = std::cos(a), sz = std::sin(z), sa = std::sin(a);

  /*
   *   H = ||{u, v, w}||
   */
  auto H = std::hypot(std::hypot(u, v), w);
  auto dH_du = H > 0 ? u / H : 0;
  auto dH_dv = H > 0 ? v / H : 0;
  auto dH_dw = H > 0 ? w / H : 0;

  /*
   *              ( / cos(a) sin(z) \   / u \   // || / u \ || )
   *  theta = acos( | sin(a) sin(z) | o | v |  //  || | v | || )
   *              ( \        cos(z) /   \ w / //   || \ w / || )
   */
  auto x = u * sz * ca + v * sa * sz + w * cz,
       d = std::sqrt(1 - std::pow(x / H, 2)) * H * H * H;

  auto theta = H > 0 ? std::acos(x / H) : std::acos(0);
  auto dtheta_du = d not_eq 0 ? (u * x - H * H * sz * ca) / d : 0;
  auto dtheta_dv = d not_eq 0 ? (v * x - H * H * sa * sz) / d : 0;
  auto dtheta_dw = d not_eq 0 ? (w * x - H * H * cz) / d : 0;

  /*
   *              ( / -sin(a) \   [ / u \   / cos(a) sin(z) \   / u \   / u \ ]   / cos(a) sin(z) \   // / -sin(a) \   [ / u \   / cos(a) sin(z) \   / u \   / u \ ] )
   *    eta = atan( |  cos(a) | x [ | v | - | sin(a) sin(z) | o | v | * | v | ] o | sin(a) sin(z) |  //  |  cos(a) | o [ | v | - | sin(a) sin(z) | o | v | * | v | ] )
   *              ( \    0    /   [ \ w /   \        cos(z) /   \ w /   \ w / ]   \        cos(z) / //   \    0    /   [ \ w /   \        cos(z) /   \ w /   \ w / ] )
   */
  auto p = std::pow(u * sa - v * ca, 2) +
           std::pow(u * ca * cz + v * sa * cz - w * sz, 2);

  auto eta = std::atan2(u * ca * cz + v * sa * cz - w * sz, u * sa - v * ca);
  auto deta_du = p not_eq 0 ? (-v * cz + w * sa * sz) / p : 0;
  auto deta_dv = p not_eq 0 ? (u * cz - w * sz * ca) / p : 0;
  auto deta_dw = p not_eq 0 ? -(u * sa - v * ca) * sz / p : 0;

  return {H,
          theta,
          eta,
          dH_du,
          dH_dv,
          dH_dw,
          dtheta_du,
          dtheta_dv,
          dtheta_dw,
          deta_du,
          deta_dv,
          deta_dw};
}

/** Sets Derived from predefined Derived parameters
 * 
 * @param[in] H Derived magnetic field strength
 * @param[in] theta Derived magnetic field theta angle
 * @param[in] eta Derived magnetic field eta angle
 * 
 * @return The pre-derived plane
 */
inline Derived FromPreDerived(Numeric H,
                              Numeric theta,
                              Numeric eta) noexcept {
  return {H, theta, eta, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}
};  // namespace Zeeman

#endif /* zeemandata_h */

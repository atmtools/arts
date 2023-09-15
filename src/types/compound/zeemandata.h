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

#include "arts_conversions.h"
#include "file.h"
#include "mystring.h"
#include <matpack.h>
#include <rtepack.h>
#include "quantum_numbers.h"

#include <limits>

/** Implements Zeeman modeling */
namespace Zeeman {
  
/** Zeeman polarization selection */
enum class Polarization : char { SigmaMinus, Pi, SigmaPlus, None };

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
    case Polarization::None:
      return 0;
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
      return -min(Ju, Jl);
    case Polarization::SigmaPlus:
      return -Ju;
    case Polarization::None:
      return 0;
  }
  return std::numeric_limits<Index>::max();
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
      return min(Ju, Jl);
    case Polarization::SigmaPlus:
      if (Ju < Jl)
        return Ju + 1;
      else if (Ju == Jl)
        return Ju;
      else
        return Jl;
    case Polarization::None:
      return 0;
  }
  return std::numeric_limits<Index>::max();
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
constexpr Index size(Rational Ju, Rational Jl, Polarization type) noexcept {
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
    case Polarization::None:
      return 1.0;
  }
  return std::numeric_limits<Numeric>::max();
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
  if (NN not_eq 0) {
    auto T1 = ((JJ + SS - NN) / JJ / 2).toNumeric();
    auto T2 = ((JJ - SS + NN) * LL / NN / JJ / 2).toNumeric();
    return GS * T1 + GL * T2;
  }      
  auto T1 = ((JJ + SS - NN) / JJ / 2).toNumeric();
  return GS * T1;
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

  if (JJ == Rational(0))
    return 0.0;
  auto DIV = Omega / JJ;
  auto T1 = (Sigma * DIV).toNumeric();
  auto T2 = (Lambda * DIV).toNumeric();
  return GS * T1 + GL * T2;
 
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
   constexpr Model(SplittingData gs = {NAN, NAN}) noexcept : mdata(gs) {}
   
   /** Default copy/init of Model from its only private variable */
   constexpr Model(Numeric gu, Numeric gl) noexcept : Model(SplittingData{gu, gl}) {}
  
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
  explicit Model(const QuantumIdentifier& qid) noexcept;

  /** Returns true if the Model represents no Zeeman effect */
  [[nodiscard]] /* constexpr */ bool empty() const noexcept {
    return std::isnan(mdata.gu) and std::isnan(mdata.gl);
  }

  /** Returns the upper state g */
  constexpr Numeric& gu() noexcept { return mdata.gu; }
  
  /** Returns the lower state g */
  constexpr Numeric& gl() noexcept { return mdata.gl; }

  /** Sets the upper state g */
  constexpr void gu(Numeric x) noexcept { mdata.gu = x; }
  
  /** Sets the lower state g */
  constexpr void gl(Numeric x) noexcept { mdata.gl = x; }
  
  /** Returns the upper state g */
  [[nodiscard]] constexpr Numeric gu() const noexcept { return mdata.gu; }
  
  /** Returns the lower state g */
  [[nodiscard]] constexpr Numeric gl() const noexcept { return mdata.gl; }

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
  [[nodiscard]] Numeric Strength(Rational Ju, Rational Jl, Polarization type, Index n) const ARTS_NOEXCEPT;
  
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
  [[nodiscard]] constexpr Numeric Splitting(Rational Ju, Rational Jl, Polarization type, Index n) const
      noexcept {
    using Constant::bohr_magneton;
    using Constant::h;
    constexpr Numeric C = bohr_magneton / h;

    return C * (Ml(Ju, Jl, type, n) * gl() - Mu(Ju, Jl, type, n) * gu());
  }

  /** Output operator for Zeeman::Model */
  friend std::ostream& operator<<(std::ostream& os, const Model& m);
  
  /** Input operator for Zeeman::Model */
  friend std::istream& operator>>(std::istream& is, Model& m);
  
  /** Output operator for Zeeman::Model */
  friend std::ostream& operator<<(bofstream& bof, const Model& m);
  
  /** Input operator for Zeeman::Model */
  friend std::istream& operator>>(bifstream& bif, Model& m);
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
Model GetSimpleModel(const QuantumIdentifier& qid) ARTS_NOEXCEPT;

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
Model GetAdvancedModel(const QuantumIdentifier& qid) ARTS_NOEXCEPT;

/** Polarization vector for Zeeman Propagation Matrix
 * 
 * Meant to contain the polarization state in two vectors
 * representing [a,b,c,d] and [u,v,w] of PropagationMatrix
 * class
 */
struct PolarizationVector {
  Vector4 att{0, 0, 0, 0};  // attenuation vector
  Vector3 dis{0, 0, 0};     // dispersion vector

  /** Default init of class */
  PolarizationVector(Numeric a = 1,
                     Numeric b = 0,
                     Numeric c = 0,
                     Numeric d = 0,
                     Numeric u = 0,
                     Numeric v = 0,
                     Numeric w = 0) noexcept
      : att({a, b, c, d}), dis({u, v, w}) {};
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
AllPolarizationVectors AllPolarization(Numeric theta,
                                       Numeric eta) noexcept;

/** The derivative of AllPolarization wrt theta
 * 
 * @param[in] theta The angle along the magnetic field
 * @param[in] eta The angle counter-clockwise in the magnetic field plane
 * 
 * @return The derivative of AllPolarization wrt theta
 */
AllPolarizationVectors AllPolarization_dtheta(
    Numeric theta, const Numeric eta) noexcept;

/** The derivative of AllPolarization wrt eta
 * 
 * @param[in] theta The angle along the magnetic field
 * @param[in] eta The angle counter-clockwise in the magnetic field plane
 * 
 * @return The derivative of AllPolarization wrt eta
 */
AllPolarizationVectors AllPolarization_deta(Numeric theta,
                                            Numeric eta) noexcept;

/** Selects the polarization vector depending on polarization type 
 * 
 * @param[in] data The pre-computed polarization vectors
 * @param[in] type The type of polarization to select
 */
const PolarizationVector& SelectPolarization(
    const AllPolarizationVectors& data, Polarization type) noexcept;

/** Sums the Zeeman components into a propagation matrix
 *
 * @param[in,out] pm The propagation matrix
 * @param[in] abs The complex absorption vector
 * @param[in] polvec The polarization vector
 */
void sum_propmat(PropmatVectorView pm, const ConstComplexVectorView &abs,
                 const PolarizationVector &polvec);

/** Sums the Zeeman components into a source vector
 *
 * @param[in,out] sv The source vector
 * @param[in] abs The complex absorption vector
 * @param[in] polvec The polarization vector
 */
void sum_stokvec(StokvecVectorView sv, const ConstComplexVectorView &abs,
                 const PolarizationVector &polvec);

/** Sums the Zeeman components derivatives into a propagation matrix
 *
 * @param[in,out] pm The propagation matrix derivative
 * @param[in] abs The complex absorption vector
 * @param[in] dabs The complex absorption vector derivative w.r.t. H
 * @param[in] polvec The polarization vector
 * @param[in] dpolvec_dtheta The polarization vector derivative w.r.t. theta
 * @param[in] dpolvec_deta The polarization vector derivative w.r.t. eta
 * @param[in] dH The derivative w.r.t. H
 * @param[in] dtheta The derivative w.r.t. theta
 * @param[in] deta The derivative w.r.t. eta
 */
void dsum_propmat(PropmatVectorView dpm, const ConstComplexVectorView &abs,
                  const ConstComplexVectorView &dabs,
                  const PolarizationVector &polvec,
                  const PolarizationVector &dpolvec_dtheta,
                  const PolarizationVector &dpolvec_deta, const Numeric dH,
                  const Numeric dtheta, const Numeric deta);

/** Sums the Zeeman components derivatives into a source vector
 *
 * @param[in,out] dsv The source vector derivative
 * @param[in] abs The complex absorption vector
 * @param[in] dabs The complex absorption vector derivative w.r.t. H
 * @param[in] polvec The polarization vector
 * @param[in] dpolvec_dtheta The polarization vector derivative w.r.t. theta
 * @param[in] dpolvec_deta The polarization vector derivative w.r.t. eta
 * @param[in] dH The derivative w.r.t. H
 * @param[in] dtheta The derivative w.r.t. theta
 * @param[in] deta The derivative w.r.t. eta
 */
void dsum_stokvec(StokvecVectorView dsv, const ConstComplexVectorView &abs,
                  const ConstComplexVectorView &dabs,
                  const PolarizationVector &polvec,
                  const PolarizationVector &dpolvec_dtheta,
                  const PolarizationVector &dpolvec_deta, const Numeric dH,
                  const Numeric dtheta, const Numeric deta);

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
 * When done and if everything is well-defined:
 * \f[ H = \sqrt{u^2 + v^2 + w^2}, \f]
 * \f[ \theta = \arccos \left( \vec{n} \cdot \vec{n}_H \right), \f]
 * \f[ \eta = \arctan\left( \frac{y}{x} \right), \f]
 * \f[ \frac{\partial H}{\partial \vec{H}} = \vec{n}_H, \f]
 * \f[ \frac{\partial \theta}{\partial \vec{H}} = \frac{\vec{n}_H \cos{\theta} - \vec{n}}{H\sin\theta}, \f]
 * \f[ \frac{\partial \eta}{\partial \vec{H}} = \frac{\vec{n}\times\vec{n}_H}{H\left(x^2 + y^2\right)} \f]
 * 
 * With these helpers (some defined, others not):
 * \f[ \vec{H} = \left[\begin{array}{l} v \\ u \\ w \end{array}\right], \f]
 * \f[ \vec{n}_H = \frac{\vec{H}}{H} , \f]
 * \f[ \vec{n} = \left[\begin{array}{r} \cos a \sin z \\ \sin a\sin z \\ \cos z \end{array}\right], \f]
 * \f[ \vec{e}_v = \left[\begin{array}{r} \cos a \cos z \\ \sin a\cos z \\ -\sin z \end{array}\right], \f]
 * \f[ y = \left\{\vec{e}_v \times \left[\vec{n}_H - \left(\vec{n}_H\cdot\vec{n}\right)\vec{n}\right]\right\} \cdot \vec{n}, \f]
 * \f[ x = \vec{e}_v \cdot \left[\vec{n}_H - \left(\vec{n}_H\cdot\vec{n}\right)\vec{n}\right] \f]
 * 
 * Note that all other values are zero if \f$ H \f$ is zero, that \f$ \frac{\partial \theta}{\partial \vec{H}} \f$
 * is zero if \f$ \sin{\theta} \f$ is zero, that \f$ \frac{\partial \eta}{\partial \vec{H}} \f$ is zero
 * if \f$ x \f$ and \f$ y \f$ are zero, and that the atan2(y, x) function is used for \f$ \eta \f$ to
 * compensate for when \f$ x \f$ is zero.
 * 
 * @param[in] u Magnetic field u-parameter
 * @param[in] v Magnetic field b-parameter
 * @param[in] w Magnetic field w-parameter
 * @param[in] z Zenith angle
 * @param[in] a Azimuth angle
 * 
 * @return The derived plane
 */
Derived FromGrids(Numeric u, Numeric v, Numeric w, Numeric z, Numeric a) noexcept;

/** Sets Derived from predefined Derived parameters
 * 
 * @param[in] H Derived magnetic field strength
 * @param[in] theta Derived magnetic field theta angle
 * @param[in] eta Derived magnetic field eta angle
 * 
 * @return The pre-derived plane
 */
constexpr Derived FromPreDerived(Numeric H,
                                 Numeric theta,
                                 Numeric eta) noexcept {
  return {H, theta, eta, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}
};  // namespace Zeeman

// Typedef to make it easier to use
using ZeemanModel = Zeeman::Model;

#endif /* zeemandata_h */

#pragma once

#include <arts_constants.h>
#include <matpack.h>
#include <quantum.h>
#include <rtepack.h>

#include <limits>

/** Implements Zeeman modeling */
namespace lbl::zeeman {
/** Zeeman polarization selection */
enum class pol : char { sm, pi, sp, no };

/** Gives the change of M given a polarization type
 * 
 * @param[in] type The polarization type
 * 
 * @return The change in M
 */
constexpr Index dM(pol type) noexcept {
  switch (type) {
    case pol::sm: return -1;
    case pol::pi: return 0;
    case pol::sp: return 1;
    case pol::no: return 0;
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
constexpr Rational start(Rational Ju, Rational Jl, pol type) noexcept {
  switch (type) {
    case pol::sm:
      if (Ju < Jl)
        return -Ju;
      else if (Ju == Jl)
        return -Ju + 1;
      else
        return -Ju + 2;
    case pol::pi: return -minr(Ju, Jl);
    case pol::sp: return -Ju;
    case pol::no: return 0;
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
constexpr Rational end(Rational Ju, Rational Jl, pol type) noexcept {
  switch (type) {
    case pol::sm: return Ju + 1;
    case pol::pi: return minr(Ju, Jl);
    case pol::sp:
      if (Ju < Jl)
        return Ju + 1;
      else if (Ju == Jl)
        return Ju;
      else
        return Jl;
    case pol::no: return 0;
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
constexpr Index size(Rational Ju, Rational Jl, pol type) noexcept {
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
constexpr Rational Mu(Rational Ju, Rational Jl, pol type, Index n) noexcept {
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
constexpr Rational Ml(Rational Ju, Rational Jl, pol type, Index n) noexcept {
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
constexpr Numeric polarization_factor(pol type) noexcept {
  switch (type) {
    case pol::sm: return .75;
    case pol::pi: return 1.5;
    case pol::sp: return .75;
    case pol::no: return 1.0;
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

  if (JJ == 0) return 0.0;
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

  if (JJ == Rational(0)) return 0.0;
  auto DIV = Omega / JJ;
  auto T1  = (Sigma * DIV).toNumeric();
  auto T2  = (Lambda * DIV).toNumeric();
  return GS * T1 + GL * T2;
}

/** Main storage for Zeeman splitting coefficients
 * 
 * The splitting data has an upper (gu) and lower (gl)
 * component and this stores both of them to not confuse
 * them elsewhere
 */
struct data {
  Numeric gu{0}, gl{0};
  constexpr auto operator<=>(const data &) const noexcept = default;
};

/** Main Zeeman Model
 * 
 * This model contains the splitting coefficients
 * of an energy level.  Various detailed and simplified
 * initialization routines are defined.  Is also used
 * as the interface for all Zeeman computations
 */
struct model {
  data mdata{};
  bool on{false};

  /** Default init */
  constexpr model() noexcept                               = default;
  constexpr model(model &&) noexcept                       = default;
  constexpr model(const model &) noexcept                  = default;
  constexpr model &operator=(model &&) noexcept            = default;
  constexpr model &operator=(const model &) noexcept       = default;
  constexpr auto operator<=>(const model &) const noexcept = default;
  constexpr model(data d) noexcept : mdata(d), on(not empty()) {};

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
  explicit model(const QuantumIdentifier &qid) noexcept;

  /** Returns true if the Model represents no Zeeman effect */
  [[nodiscard]] constexpr bool empty() const noexcept {
    return 0.0 == mdata.gu and mdata.gl == 0.0;
  }

  /** Returns the upper state g */
  constexpr Numeric &gu() noexcept { return mdata.gu; }

  /** Returns the lower state g */
  constexpr Numeric &gl() noexcept { return mdata.gl; }

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
  [[nodiscard]] Numeric Strength(Rational Ju,
                                 Rational Jl,
                                 pol type,
                                 Index n) const;

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
  [[nodiscard]] Numeric Strength(const QuantumState &qn,
                                 pol type,
                                 Index n) const;

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
  [[nodiscard]] constexpr Numeric Splitting(Rational Ju,
                                            Rational Jl,
                                            pol type,
                                            Index n) const noexcept {
    using Constant::bohr_magneton;
    using Constant::h;
    constexpr Numeric C = bohr_magneton / h;

    return (type == pol::no)
               ? 0.0
               : C * (Ml(Ju, Jl, type, n) * gl() - Mu(Ju, Jl, type, n) * gu());
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
  [[nodiscard]] Numeric Splitting(const QuantumState &qn,
                                  pol type,
                                  Index n) const noexcept;

  /** Gives the number of lines for a given polarization type
   * 
   * @param[in] qn The quantum numbers of the local state
   * @param[in] type The polarization type
   * 
   * @return The splitting of the Zeeman subline
   */
  [[nodiscard]] Index size(const QuantumState &qn, pol type) const noexcept;

  /** Input operator for Zeeman::Model */
  friend std::istream &operator>>(std::istream &is, model &m);
};  // Model;

/** Returns a simple Zeeman model 
 * 
 * Will use the simple Hund case provided
 * by input.  Throws if the input is bad
 * 
 * @param[in] qid Transition type quantum id
 * 
 * @return Zeeman model data
 */
data GetSimpleModel(const QuantumIdentifier &qid);

/** Returns an advanced Zeeman model 
 * 
 * Will look at available Quantum numbers
 * and use best approximation for the model
 * to use.  If no good approximation is available,
 * it returns Model({0, 0}).
 * 
 * @param[in] qid Transition type quantum id
 * 
 * @return Zeeman model data
 */
data GetAdvancedModel(const QuantumIdentifier &qid);

struct magnetic_angles {
  Numeric u, v, w, sa, ca, sz, cz, H, uct, duct;

  magnetic_angles(const Vector3 mag = {0, 0, 0}, const Vector2 los = {0, 0});

  [[nodiscard]] Numeric theta() const;
  [[nodiscard]] Numeric dtheta_du() const;
  [[nodiscard]] Numeric dtheta_dv() const;
  [[nodiscard]] Numeric dtheta_dw() const;
  [[nodiscard]] Numeric eta() const;
  [[nodiscard]] Numeric deta_du() const;
  [[nodiscard]] Numeric deta_dv() const;
  [[nodiscard]] Numeric deta_dw() const;
};

Propmat norm_view(pol p, Vector3 mag, Vector2 los);

Propmat dnorm_view_du(pol p, Vector3 mag, Vector2 los);

Propmat dnorm_view_dv(pol p, Vector3 mag, Vector2 los);

Propmat dnorm_view_dw(pol p, Vector3 mag, Vector2 los);

constexpr Propmat scale(const Propmat &a, const Complex F) noexcept {
  return {a.A() * F.real(),
          a.B() * F.real(),
          a.C() * F.real(),
          a.D() * F.real(),
          a.U() * F.imag(),
          a.V() * F.imag(),
          a.W() * F.imag()};
}

constexpr Propmat scale(const Propmat &a,
                        const Propmat &da,
                        const Complex F,
                        const Complex dF) noexcept {
  return {da.A() * F.real() + a.A() * dF.real(),
          da.B() * F.real() + a.B() * dF.real(),
          da.C() * F.real() + a.C() * dF.real(),
          da.D() * F.real() + a.D() * dF.real(),
          da.U() * F.imag() + a.U() * dF.imag(),
          da.V() * F.imag() + a.V() * dF.imag(),
          da.W() * F.imag() + a.W() * dF.imag()};
}
};  // namespace lbl::zeeman

template <>
struct std::formatter<lbl::zeeman::model> {
  format_tags tags;

  [[nodiscard]] constexpr auto &inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto &inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context &ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const lbl::zeeman::model &v,
                              FmtContext &ctx) const {
    if (tags.help) {
      if (v.on) {
        tags.format(
            ctx, "<on>; Upper state: "sv, v.gu(), "; Lower state: "sv, v.gl());
      } else {
        tags.format(ctx, "<off>"sv);
      }
    } else if (tags.io) {
      tags.format(ctx, Index{v.on}, ' ', v.gu(), ' ', v.gl());
    } else {
      const auto sep = tags.sep();
      tags.add_if_bracket(ctx, '[');
      tags.format(ctx, v.on, sep, v.gu(), sep, v.gl());
      tags.add_if_bracket(ctx, ']');
    }

    return ctx.out();
  }
};

using ZeemanLineModel = lbl::zeeman::model;

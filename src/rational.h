/* Copyright (C) 2012 
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

/**
 * @file   rational.h
 * @author Richard Larsson
 * @date   2012-10-31
 * 
 * @brief  Contains the rational class definition
 **/

#ifndef rational_h
#define rational_h

#include <cassert>
#include <ostream>
#include "array.h"
#include "bifstream.h"
#include "bofstream.h"
#include "math_funcs.h"
#include "matpack.h"

/** Implements rational numbers to work with other ARTS types */
class Rational {
 public:
  /** Initialization call
   * 
   * @param[in] nom Nominator
   * @param[in] denom Denominator
   */
  constexpr Rational(const Index nom = 0, const Index denom = 1)
      : mnom(nom), mdenom(denom) {}
  
  Rational(const String& s);

  /** Nominator */
  constexpr Index Nom() const { return mnom; }

  /** Denominator */
  constexpr Index Denom() const { return mdenom; }

  /** Is the object not defined
   * 
   * @return true If Denom() is 0
   * @return false Otherwise
   */
  constexpr bool isUndefined() const { return (mdenom == 0); }

  /** Is the object defined
   * 
   * @return true If Denom() is not 0
   * @return false Otherwise
   */
  constexpr bool isDefined() const { return mdenom not_eq 0; }

  /** Is the object an Index
   * 
   * @return true If Nom() % Denom() is 0
   * @return false Otherwise
   */
  constexpr bool isIndex() const {
    return isDefined() and not bool(mnom % mdenom);
  }

  /** Converts the value to index by direct division
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @return constexpr Index Of *this
   */
  constexpr Index toIndex() const {
    return isIndex() ? mnom / mdenom
                     : throw std::logic_error(
                           "Rational is not representative of an integer");
  }

  /** Converts this to a Numeric
   * 
   * @return constexpr Numeric of *this
   */
  constexpr Numeric toNumeric() const {
    return Numeric(mnom) / Numeric(mdenom);
  }

  /** Same as int(toIndex()) */
  constexpr int toInt() const { return int(toIndex()); }

  /** Simplify by reducing to smallest possible denominator */
  Rational& Simplify();

  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  Rational& operator+=(const Rational& a) {
    mnom = mnom * a.Denom() + a.Nom() * mdenom;
    mdenom *= a.Denom();
    return *this;
  }

  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  Rational& operator+=(const Index& a) {
    mnom += mdenom * a;
    return *this;
  }

  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  Rational& operator-=(const Rational& a) {
    mnom = mnom * a.Denom() - a.Nom() * mdenom;
    mdenom *= a.Denom();
    return *this;
  }

  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  Rational& operator-=(const Index& a) {
    mnom -= mdenom * a;
    return *this;
  }

  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  Rational& operator/=(const Rational& a) {
    mnom *= a.Denom();
    mdenom *= a.Nom();
    return *this;
  }

  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  Rational& operator/=(const Index& a) {
    mdenom *= a;
    return *this;
  }

  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  Rational& operator*=(const Rational& a) {
    mnom *= a.Nom();
    mdenom *= a.Denom();
    return *this;
  }

  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  Rational& operator*=(const Index& a) {
    mnom *= a;
    return *this;
  }

 /** Add one if possible */
  constexpr Rational operator++(int) const {
    return isDefined()
               ? Rational(mnom + mdenom, mdenom)
               : throw std::logic_error("Undefined Rational in operator++");
  }

  /** Remove one if possible */
  constexpr Rational operator--(int) const {
    return isDefined()
               ? Rational(mnom - mdenom, mdenom)
               : throw std::logic_error("Undefined Rational in operator--");
  }

 /** Add one if possible */
  Rational& operator++() {
    mnom += mdenom;
    return *this;
  }

  /** Remove one if possible */
  Rational& operator--() {
    mnom -= mdenom;
    return *this;
  }

  /** Cast to bool */
  explicit constexpr operator bool() const {
    return isDefined() and bool(mnom);
  }

  /** Cast to Numeric */
  explicit constexpr operator Numeric() const { return toNumeric(); }

  /** Cast to Index */
  explicit constexpr operator Index() const { return toIndex(); }

  /** Cast to int */
  explicit constexpr operator int() const { return toInt(); }
  
  /** Binary read for Rational */
  bifstream& read(bifstream& bif) {
    bif >> mnom >> mdenom;
    return bif;
  }
  
  /** Binary write for Rational */
  bofstream& write(bofstream& bof) const {
    bof << mnom << mdenom;
    return bof;
  }

 private:
  // Rational is supposed to be used rationally ( mnom / mdenom )
  Index mnom;
  Index mdenom;

  /** Makes the sign of mdenom positive */
  Rational& fixSign() {
    if (mdenom < 0) {
      mnom = -mnom;
      mdenom = -mdenom;
    }
    return *this;
  }
};  // Rational;

// An undefined rational to be used everywhere for a undefined rationals
#define RATIONAL_UNDEFINED Rational(0, 0)

/** Negative
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational Negative a
 */
constexpr Rational operator-(Rational a) {
  return Rational(-a.Nom(), a.Denom());
}

/** Positive
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational a
 */
constexpr Rational operator+(Rational a) { return a; }

/** Addition
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(Rational a, Rational b) {
  return (a.Denom() == b.Denom())
             ? Rational(a.Nom() + b.Nom(), a.Denom())
             : Rational(a.Nom() * b.Denom() + b.Nom() * a.Denom(),
                        a.Denom() * b.Denom());
}

/** Addition
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(Rational a, Index b) {
  return Rational(a.Nom() + b * a.Denom(), a.Denom());
}

/** Addition
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(Index b, Rational a) { return operator+(a, b); }

/** Subtraction
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a - b
 */
constexpr Rational operator-(Rational a, Rational b) {
  return (a.Denom() == b.Denom())
             ? Rational(a.Nom() - b.Nom(), a.Denom())
             : Rational(a.Nom() * b.Denom() - b.Nom() * a.Denom(),
                        a.Denom() * b.Denom());
}

/** Subtraction
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a - b
 */
constexpr Rational operator-(Rational a, Index b) {
  return Rational(a.Nom() - b * a.Denom(), a.Denom());
}

/** Subtraction
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b - a
 */
constexpr Rational operator-(Index b, Rational a) {
  return Rational(-a.Nom() + b * a.Denom(), a.Denom());
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(Rational a, Rational b) {
  return Rational(a.Nom() * b.Denom(), a.Denom() * b.Nom());
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(Rational a, Index b) {
  return Rational(a.Nom(), a.Denom() * b);
}

/** Division
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b / a
 */
constexpr Rational operator/(Index b, Rational a) {
  return Rational(a.Denom() * b, a.Nom());
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(Rational a, Rational b) {
  return Rational(a.Nom() * b.Nom(), a.Denom() * b.Denom());
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(Rational a, Index b) {
  return Rational(a.Nom() * b, a.Denom());
}

/** Multiplication
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(Index b, Rational a) { return operator*(a, b); }

/** Remainder
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a % b
 */
constexpr Rational operator%(Rational a, Rational b) {
  return (a.Denom() == b.Denom())
             ? Rational(a.Nom() % b.Nom(), a.Denom())
             : Rational((a.Nom() * b.Denom()) % (a.Denom() * b.Nom()),
                        a.Denom() * b.Denom());
}

/** Remainder
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a % b
 */
constexpr Rational operator%(Rational a, Index b) {
  return Rational(a.Nom() % b, a.Denom());
}

/** Remainder
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b % a
 */
constexpr Rational operator%(Index b, Rational a) {
  return Rational((b * a.Denom()) % a.Nom(), a.Denom());
}

/** Equality
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If equal
 * @return false Otherwise
 */
constexpr bool operator==(Rational a, Rational b) {
  return a.isDefined() and b.isDefined() and
         a.Nom() * b.Denom() == a.Denom() * b.Nom();
}

/** Inequality
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If not equal
 * @return false Otherwise
 */
constexpr bool operator!=(Rational a, Rational b) {
  return not operator==(a, b);
}

/** Less than
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a < b
 * @return false Otherwise
 */
constexpr bool operator<(Rational a, Rational b) {
  return a.isDefined() and b.isDefined() and
         a.Nom() * b.Denom() < a.Denom() * b.Nom();
}

/** More than
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a > b
 * @return false Otherwise
 */
constexpr bool operator>(Rational a, Rational b) { return operator<(b, a); }

/** Less than or equal to
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a <= b
 * @return false Otherwise
 */
constexpr bool operator<=(Rational a, Rational b) {
  return not operator>(a, b);
}

/** More than or equal to
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a >= b
 * @return false Otherwise
 */
constexpr bool operator>=(Rational a, Rational b) {
  return not operator<(a, b);
}

/** Not
 * 
 * @param[in] a Any Rational
 * @return true If a.Nom() and a.isDefined()
 * @return false Otherwise
 */
constexpr bool operator!(Rational a) { return a.Nom() and a.isDefined(); }

/** Factorial
 * 
 * @param[in] r Any Rational
 * @return Numeric Factorial of the Rational
 */
inline Numeric fac(Rational r) { return (::fac(r.toIndex())); }

/** Square root
 * 
 * @param[in] r Any Rational
 * @return Numeric Square root of the Rational
 */
inline Numeric sqrt(Rational r) { return std::sqrt(r.toNumeric()); }

/** Power of
 * 
 * @param[in] base Any Rational
 * @param[in] exp Any Rational
 * @return Numeric base to the power of exp
 */
inline Numeric pow(Rational base, Rational exp) {
  return std::pow(Numeric(base), Numeric(exp));
}

/** Output operator */
std::ostream& operator<<(std::ostream& os, const Rational& a);

/** Input operator */
std::istream& operator>>(std::istream& os, Rational& a);

/** Absolute
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational Absolute value of the Rational
 */
constexpr Rational abs(Rational a) {
  return a < 0 ? -a : a;
}

/** Maximum
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational Largest of a and b
 */
constexpr Rational& max(Rational& a, Rational& b) {
  return a < b ? b : a;
}  // Let other operators find out if this is allowed instead

typedef Array<Rational> ArrayOfRational;

#endif  // rational_h

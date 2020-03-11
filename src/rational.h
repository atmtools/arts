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


/** Returns the greatest common denominator of two numbers
 * 
 * @param[in] a number a
 * @param[in] b number b
 * @return num such that Rational(a/num, b/num) is the same as Rational(a, b)
 */
constexpr Index gcd(Index a, Index b) {
  if (b == 0)
    return a;
  else
    return gcd(b, a%b);
}


/** Implements rational numbers to work with other ARTS types */
class Rational {
 public:
  /** Initialization call
   * 
   * @param[in] nom Nominator
   * @param[in] denom Denominator
   */
  constexpr Rational(const Index nom = 0, const Index denom = 1)
  : mnom(denom ? nom : 0), mdenom(denom) {
    const auto div = gcd(nom, denom);
    if (div) {
      mnom /= div;
      mdenom /= div;
    }
  }
  
  /** Initialization call
   * 
   * Sets the rational from the string. Formats accepted are
   * 
   * Numeric:  1.234567890
   * Fraction: 12345/67890
   * Index:    1234567890
   * 
   * Note that overflow is possible and we do not care to capture it
   * 
   * @param[in] s String of the value
   */
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

  /** Is the object a n-scaled Index
   * 
   * @param[in]  n Scale
   * @return true If n*Nom() % Denom() is 0
   * @return false Otherwise
   */
  constexpr bool isIndex(int n=1) const {
    return isDefined() and not bool((n*mnom) % mdenom);
  }

  /** Converts the value to index by n-scaled division
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @param[in]  n Scale to *this
   * @return constexpr Index Of *this
   */
  constexpr Index toIndex(int n=1) const {
    return isIndex(n) ? (n*mnom) / mdenom
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
  
  /** Converts the value to int by n-scaled division in Index form
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @param[in]  n Scale to *this
   * @return constexpr int Of *this
   */
  constexpr int toInt(int n=1) const { return int(toIndex(n)); }
  
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
  
  /** Makes the sign of mdenom positive */
  constexpr Rational& fixSign() {
    if (mdenom < 0) {
      mnom = -mnom;
      mdenom = -mdenom;
    }
    return *this;
  }

 private:
  // Rational is supposed to be used rationally ( mnom / mdenom )
  Index mnom;
  Index mdenom;
};  // Rational;

/** Returns the rational reduced by the greates 
 * 
 * @param[in] a Any Rational
 * @return a / gcd(a)
 */
constexpr Rational reduce_by_gcd(Rational a) {
  const Index div = gcd(a.Nom(), a.Denom());
  if (div)
    return Rational(a.Nom() / div, a.Denom());
  else
    return a;
}

/** Rational from Numeric
 * 
 * Performs basic rounding
 * 
 * @param[in] x Numeric value
 * @param[in] maxdec Maximum number of decimals
 */
constexpr Rational numeric2rational(Numeric x, size_t maxdec=4) {
  Index nom=0, denom=1;
  
  // Keep track of sign independently
  const bool signchange = x < 0;
  x = x < 0 ? -x : x;
  
  // Add numbers by keeping the floor
  size_t i=0;
  do {
    const Index xi=Index(x);
    nom += xi;
    x = 10 * (x - Numeric(xi));
    nom *= 10;
    denom *= 10;
    i++;
  } while (i<=maxdec);
  
  // Fix possible rounding error
  if (x >= 5)
    nom += 10;
  
  // Change sign or not
  if (signchange)
    return Rational(-nom, denom);
  else
    return Rational(nom, denom);
}


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
  return Rational(a.Nom() % (a.Denom() * b), a.Denom());
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
 * @param[in] exp Any Numeric
 * @return Numeric base to the power of exp
 */
inline Numeric pow(Rational base, Numeric exp) {
  return std::pow(base.toNumeric(), exp);
}

/** Power of
 * 
 * @param[in] base Any Numeric
 * @param[in] exp Any Rational
 * @return Numeric base to the power of exp
 */
inline Numeric pow(Numeric base, Rational exp) {
  return std::pow(base, exp.toNumeric());
}

/** Power of
 * 
 * @param[in] base Any Rational
 * @param[in] exp Any Rational
 * @return Numeric base to the power of exp
 */
inline Numeric pow(Rational base, Rational exp) {
  return pow(base, exp.toNumeric());
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

/** Returns common operator n/2
 * 
 * @param[in] n Any positive integer
 * @return Rational(n, 2)
 */
constexpr Rational operator ""_2(unsigned long long int n) {
  return Rational(n, 2);
};

/** Returns true if even integer
 * 
 * @param[in]  r Any rational
 * @return  true if r is even, otherwise false
 */
constexpr bool even(Rational r) {
  if (r % 2)
    return false;
  else
    return true;
}

#endif  // rational_h

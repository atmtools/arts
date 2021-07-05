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
constexpr Index gcd(Index a, Index b) noexcept {
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
  constexpr Rational(const Index nom = 0, const Index denom = 1) noexcept
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
  explicit Rational(const String& s);

  /** Nominator */
  constexpr Index Nom() const noexcept { return mnom; }

  /** Denominator */
  constexpr Index Denom() const noexcept { return mdenom; }
  
  /** Nominator */
  constexpr Index& Nom() noexcept { return mnom; }
  
  /** Denominator */
  constexpr Index& Denom() noexcept { return mdenom; }
  
  /** Nominator */
  constexpr void Nom(Index x) noexcept { mnom = x; }
  
  /** Denominator */
  constexpr void Denom(Index x) noexcept { mdenom = x; }
  
  /** Simplify by reducing the values locally */
  void simplify_in_place() noexcept;

  /** Is the object not defined
   * 
   * @return true If Denom() is 0
   * @return false Otherwise
   */
  constexpr bool isUndefined() const noexcept { return (mdenom == 0); }

  /** Is the object defined
   * 
   * @return true If Denom() is not 0
   * @return false Otherwise
   */
  constexpr bool isDefined() const noexcept { return not isUndefined(); }

  /** Is the object a n-scaled Index
   * 
   * @param[in]  n Scale
   * @return true If n*Nom() % Denom() is 0
   * @return false Otherwise
   */
  constexpr bool isIndex(int n=1) const noexcept {
    return isDefined() and not bool((n*mnom) % mdenom);
  }

  /** Converts the value to index by n-scaled division
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @param[in]  n Scale to *this
   * @return constexpr Index Of *this
   */
  constexpr Index toIndex(int n=1) const noexcept {
    return (n*mnom) / mdenom;;
  }

  /** Converts this to a Numeric
   * 
   * @return constexpr Numeric of *this
   */
  constexpr Numeric toNumeric() const noexcept {
    return Numeric(mnom) / Numeric(mdenom);
  }
  
  /** Converts the value to int by n-scaled division in Index form
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @param[in]  n Scale to *this
   * @return constexpr int Of *this
   */
  constexpr int toInt(int n=1) const noexcept { return int(toIndex(n)); }
  
  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  constexpr Rational& operator+=(const Rational& a) noexcept {
    mnom = mnom * a.Denom() + a.Nom() * mdenom;
    mdenom *= a.Denom();
    return *this;
  }

  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  constexpr Rational& operator+=(const Index& a) noexcept {
    mnom += mdenom * a;
    return *this;
  }
  
  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  constexpr Rational& operator+=(const int& a) noexcept {
    mnom += mdenom * a;
    return *this;
  }

  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  constexpr Rational& operator-=(const Rational& a) noexcept {
    mnom = mnom * a.Denom() - a.Nom() * mdenom;
    mdenom *= a.Denom();
    return *this;
  }

  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  constexpr Rational& operator-=(const Index& a) noexcept {
    mnom -= mdenom * a;
    return *this;
  }
  
  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  constexpr Rational& operator-=(const int& a) noexcept {
    mnom -= mdenom * a;
    return *this;
  }

  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  constexpr Rational& operator/=(const Rational& a) noexcept {
    mnom *= a.Denom();
    mdenom *= a.Nom();
    return *this;
  }

  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  constexpr Rational& operator/=(const Index& a) noexcept {
    mdenom *= a;
    return *this;
  }
  
  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  constexpr Rational& operator/=(const int& a) noexcept {
    mdenom *= a;
    return *this;
  }

  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  constexpr Rational& operator*=(const Rational& a) noexcept {
    mnom *= a.Nom();
    mdenom *= a.Denom();
    return *this;
  }

  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  constexpr Rational& operator*=(const Index& a) noexcept {
    mnom *= a;
    return *this;
  }
  
  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  constexpr Rational& operator*=(const int& a) noexcept {
    mnom *= a;
    return *this;
  }

 /** Add one if possible */
  constexpr Rational operator++(int) const noexcept{
    return Rational(mnom + mdenom, mdenom);
  }

  /** Remove one if possible */
  constexpr Rational operator--(int) const noexcept {
    return Rational(mnom - mdenom, mdenom);
  }

 /** Add one if possible */
  constexpr Rational& operator++() noexcept {
    mnom += mdenom;
    return *this;
  }

  /** Remove one if possible */
  constexpr Rational& operator--() noexcept {
    mnom -= mdenom;
    return *this;
  }

  /** Cast to bool */
  explicit constexpr operator bool() const noexcept {
    return isDefined() and bool(mnom);
  }

  /** Cast to Numeric */
  explicit constexpr operator Numeric() const noexcept { return toNumeric(); }

  /** Cast to Index */
  explicit constexpr operator Index() const noexcept { return toIndex(); }

  /** Cast to int */
  explicit constexpr operator int() const noexcept { return toInt(); }
  
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
  constexpr Rational& fixSign() noexcept {
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
constexpr Rational reduce_by_gcd(const Rational a) noexcept {
  const Index div = gcd(a.Nom(), a.Denom());
  if (div)
    return Rational(a.Nom() / div, a.Denom() / div);
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
constexpr Rational numeric2rational(Numeric x, size_t maxdec=4) noexcept {
  Index nom=0, denom=1;
  
  // Keep track of sign independently
  const bool signchange = x < 0;
  x = signchange ? -x : x;
  
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
constexpr Rational operator-(const Rational a) noexcept {
  return Rational(-a.Nom(), a.Denom());
}

/** Positive
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational a
 */
constexpr Rational operator+(const Rational a) noexcept { return a; }

/** Addition
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(const Rational a, const Rational b) noexcept {
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
constexpr Rational operator+(const Rational a, Index b) noexcept {
  return Rational(a.Nom() + b * a.Denom(), a.Denom());
}

/** Addition
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(const Rational a, int b) noexcept {
  return Rational(a.Nom() + b * a.Denom(), a.Denom());
}

/** Addition
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(Index b, const Rational a) noexcept { return operator+(a, b); }

/** Addition
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(int b, const Rational a) noexcept { return operator+(a, b); }

/** Subtraction
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a - b
 */
constexpr Rational operator-(const Rational a, const Rational b) noexcept {
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
constexpr Rational operator-(const Rational a, Index b) noexcept {
  return Rational(a.Nom() - b * a.Denom(), a.Denom());
}

/** Subtraction
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a - b
 */
constexpr Rational operator-(const Rational a, int b) noexcept {
  return Rational(a.Nom() - b * a.Denom(), a.Denom());
}

/** Subtraction
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b - a
 */
constexpr Rational operator-(Index b, const Rational a) noexcept {
  return Rational(-a.Nom() + b * a.Denom(), a.Denom());
}

/** Subtraction
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b - a
 */
constexpr Rational operator-(int b, const Rational a) noexcept {
  return Rational(-a.Nom() + b * a.Denom(), a.Denom());
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(const Rational a, const Rational b) noexcept {
  return Rational(a.Nom() * b.Denom(), a.Denom() * b.Nom());
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(const Rational a, Index b) noexcept {
  return Rational(a.Nom(), a.Denom() * b);
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(const Rational a, int b) noexcept {
  return Rational(a.Nom(), a.Denom() * b);
}

/** Division
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b / a
 */
constexpr Rational operator/(Index b, const Rational a) noexcept {
  return Rational(a.Denom() * b, a.Nom());
}

/** Division
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b / a
 */
constexpr Rational operator/(int b, const Rational a) noexcept {
  return Rational(a.Denom() * b, a.Nom());
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(const Rational a, const Rational b) noexcept {
  return Rational(a.Nom() * b.Nom(), a.Denom() * b.Denom());
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(const Rational a, Index b) noexcept {
  return Rational(a.Nom() * b, a.Denom());
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(const Rational a, int b) noexcept {
  return Rational(a.Nom() * b, a.Denom());
}

/** Multiplication
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(Index b, const Rational a) noexcept { return operator*(a, b); }

/** Multiplication
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(int b, const Rational a) noexcept { return operator*(a, b); }

/** Remainder
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a % b
 */
constexpr Rational operator%(const Rational a, const Rational b) noexcept {
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
constexpr Rational operator%(const Rational a, Index b) noexcept {
  return Rational(a.Nom() % (a.Denom() * b), a.Denom());
}

/** Remainder
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a % b
 */
constexpr Rational operator%(const Rational a, int b) noexcept {
  return Rational(a.Nom() % (a.Denom() * b), a.Denom());
}

/** Remainder
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b % a
 */
constexpr Rational operator%(Index b, const Rational a) noexcept {
  return Rational((b * a.Denom()) % a.Nom(), a.Denom());
}

/** Remainder
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b % a
 */
constexpr Rational operator%(int b, const Rational a) noexcept {
  return Rational((b * a.Denom()) % a.Nom(), a.Denom());
}

/** Equality
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If equal
 * @return false Otherwise
 */
constexpr bool operator==(const Rational a, const Rational b) noexcept {
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
constexpr bool operator!=(const Rational a, const Rational b) noexcept {
  return not (a.isDefined() and b.isDefined() and operator==(a, b));
}

/** Less than
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a < b
 * @return false Otherwise
 */
constexpr bool operator<(const Rational a, const Rational b) noexcept {
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
constexpr bool operator>(const Rational a, const Rational b) noexcept { return operator<(b, a) and a.isDefined() and b.isDefined(); }

/** Less than or equal to
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a <= b
 * @return false Otherwise
 */
constexpr bool operator<=(const Rational a, const Rational b) noexcept {
  return not operator>(a, b) and a.isDefined() and b.isDefined();
}

/** More than or equal to
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return true If a >= b
 * @return false Otherwise
 */
constexpr bool operator>=(const Rational a, const Rational b) noexcept {
  return not operator<(a, b) and a.isDefined() and b.isDefined();
}

/** Not
 * 
 * @param[in] a Any Rational
 * @return true If a.Nom() and a.isDefined()
 * @return false Otherwise
 */
constexpr bool operator!(const Rational a) noexcept { return a.Nom() and a.isDefined(); }

/** Square root
 * 
 * @param[in] r Any Rational
 * @return Numeric Square root of the Rational
 */
inline Numeric sqrt(const Rational r) { return std::sqrt(r.toNumeric()); }

/** Power of
 * 
 * @param[in] base Any Rational
 * @param[in] exp Any Numeric
 * @return Numeric base to the power of exp
 */
inline Numeric pow(const Rational base, Numeric exp) {
  return std::pow(base.toNumeric(), exp);
}

/** Power of
 * 
 * @param[in] base Any Numeric
 * @param[in] exp Any Rational
 * @return Numeric base to the power of exp
 */
inline Numeric pow(Numeric base, const Rational exp) {
  return std::pow(base, exp.toNumeric());
}

/** Power of
 * 
 * @param[in] base Any Rational
 * @param[in] exp Any Rational
 * @return Numeric base to the power of exp
 */
inline Numeric pow(const Rational base, const Rational exp) {
  return pow(base, exp.toNumeric());
}

/** Output operator */
std::ostream& operator<<(std::ostream& os, const Rational& a);

/** Input operator */
std::istream& operator>>(std::istream& is, Rational& a);

/** less
 * 
 * @param[in] a Any Index
 * @param[in] b Any Rational
 * @return True if a is less than b
 */
constexpr bool operator<(const Index a, const Rational b) noexcept {
  return Rational(a, 1) < b and b.isDefined();
}

/** less
 * 
 * @param[in] a Any Index
 * @param[in] b Any Rational
 * @return True if a is less than b
 */
constexpr bool operator<(const int a, const Rational b) noexcept {
  return Rational(a, 1) < b and b.isDefined();
}

/** less
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is less than b
 */
constexpr bool operator<(const Rational a, const Index b) noexcept {
  return a < Rational(b, 1) and a.isDefined();
}

/** less
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is less than b
 */
constexpr bool operator<(const Rational a, const int b) noexcept {
  return a < Rational(b, 1) and a.isDefined();
}

/** more
 * 
 * @param[in] a Any Index
 * @param[in] b Any Rational
 * @return True if a is more than b
 */
constexpr bool operator>(const Index a, const Rational b) noexcept {
  return Rational(a, 1) > b and b.isDefined();
}

/** more
 * 
 * @param[in] a Any Index
 * @param[in] b Any Rational
 * @return True if a is more than b
 */
constexpr bool operator>(const int a, const Rational b) noexcept {
  return Rational(a, 1) > b and b.isDefined();
}

/** more
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is more than b
 */
constexpr bool operator>(const Rational a, const Index b) noexcept {
  return a > Rational(b, 1) and a.isDefined();
}

/** more
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is more than b
 */
constexpr bool operator>(const Rational a, const int b) noexcept {
  return a > Rational(b, 1) and a.isDefined();
}

/** equal
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is equal to b
 */
constexpr bool operator==(const Rational a, const Index b) noexcept {
  return a == Rational(b, 1) and a.isDefined();
}

/** equal
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is equal to b
 */
constexpr bool operator==(const Rational a, const int b) noexcept {
  return a == Rational(b, 1) and a.isDefined();
}

/** not equal
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is not equal to b
 */
constexpr bool operator!=(const Rational a, const Index b) noexcept {
  return not(a == b) and a.isDefined();
}

/** not equal
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return True if a is not equal to b
 */
constexpr bool operator!=(const Rational a, const int b) noexcept {
  return not(a == b) and a.isDefined();
}

/** Absolute
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational Absolute value of the Rational
 */
constexpr Rational abs(const Rational a) noexcept {
  return a < 0 ? -a : a;
}

/** Maximum
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational Largest of a and b
 */
constexpr Rational max(const Rational a, const Rational b) noexcept {
  return a < b ? b : a;
}  // Let other operators find out if this is allowed instead

/** Minimum
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational Smallest of a and b
 */
constexpr Rational min(const Rational a, const Rational b) noexcept {
  return a < b ? a : b;
}  // Let other operators find out if this is allowed instead

typedef Array<Rational> ArrayOfRational;

/** Returns common operator n/2
 * 
 * @param[in] n Any positive integer
 * @return Rational(n, 2)
 */
constexpr Rational operator ""_2(unsigned long long int n) noexcept {
  return Rational(n, 2);
};

/** Returns true if even integer
 * 
 * @param[in]  r Any rational
 * @return  true if r is even, otherwise false
 */
constexpr bool iseven(const Rational r) noexcept {
  if (r % 2)
    return false;
  else
    return true;
}

#endif  // rational_h


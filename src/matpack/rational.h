/**
 * @file   rational.h
 * @author Richard Larsson
 * @date   2012-10-31
 * 
 * @brief  Contains the rational class definition
 **/

#ifndef rational_h
#define rational_h

#include "array.h"
#include "bifstream.h"
#include "bofstream.h"
#include "math_funcs.h"
#include "matpack_concepts.h"

#include <numeric>
#include <ostream>


using std::gcd;

/** Implements rational numbers to work with other ARTS types */
struct Rational {
  Index numer;
  Index denom;

  /** Initialization call
   * 
   * @param[in] nom Nominator
   * @param[in] denom Denominator
   */
  constexpr Rational(const Index n = 0, const Index d = 1) noexcept
  : numer(d ? n / gcd(n, d) : 0), denom(d / gcd(n, d)) {}
  
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
  
  /** Simplify by reducing the values locally */
  void simplify_in_place() noexcept;

  /** Is the object not defined
   * 
   * @return true If Denom() is 0
   * @return false Otherwise
   */
  [[nodiscard]] constexpr bool isUndefined() const noexcept { return (denom == 0); }

  /** Is the object defined
   * 
   * @return true If Denom() is not 0
   * @return false Otherwise
   */
  [[nodiscard]] constexpr bool isDefined() const noexcept { return not isUndefined(); }

  /** Is the object a n-scaled Index
   * 
   * @param[in]  n Scale
   * @return true If n*Nom() % Denom() is 0
   * @return false Otherwise
   */
  [[nodiscard]] constexpr bool isIndex(int n=1) const noexcept {
    return isDefined() and not bool((n*numer) % denom);
  }

  /** Converts the value to index by n-scaled division
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @param[in]  n Scale to *this
   * @return constexpr Index Of *this
   */
  [[nodiscard]] constexpr Index toIndex(int n=1) const noexcept {
    return (n*numer) / denom;;
  }

  /** Converts this to a Numeric
   * 
   * @return constexpr Numeric of *this
   */
  [[nodiscard]] constexpr Numeric toNumeric() const noexcept {
    return Numeric(numer) / Numeric(denom);
  }
  
  /** Converts the value to int by n-scaled division in Index form
   * 
   * Throws a logic error if *this is not an Index
   * 
   * @param[in]  n Scale to *this
   * @return constexpr int Of *this
   */
  [[nodiscard]] constexpr int toInt(int n=1) const noexcept { return int(toIndex(n)); }
  
  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  constexpr Rational& operator+=(const Rational& a) noexcept {
    numer = numer * a.denom + a.numer * denom;
    denom *= a.denom;
    return *this;
  }

  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  constexpr Rational& operator+=(const Index& a) noexcept {
    numer += denom * a;
    return *this;
  }
  
  /** Add to this
   * 
   * @param[in] a To add
   * @return Rational& *this
   */
  constexpr Rational& operator+=(const int& a) noexcept {
    numer += denom * a;
    return *this;
  }

  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  constexpr Rational& operator-=(const Rational& a) noexcept {
    numer = numer * a.denom - a.numer * denom;
    denom *= a.denom;
    return *this;
  }

  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  constexpr Rational& operator-=(const Index& a) noexcept {
    numer -= denom * a;
    return *this;
  }
  
  /** Remove from this
   * 
   * @param[in] a To remove
   * @return Rational& *this
   */
  constexpr Rational& operator-=(const int& a) noexcept {
    numer -= denom * a;
    return *this;
  }

  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  constexpr Rational& operator/=(const Rational& a) noexcept {
    numer *= a.denom;
    denom *= a.numer;
    return *this;
  }

  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  constexpr Rational& operator/=(const Index& a) noexcept {
    denom *= a;
    return *this;
  }
  
  /** Divide by this
   * 
   * @param[in] a To divide by
   * @return Rational& *this
   */
  constexpr Rational& operator/=(const int& a) noexcept {
    denom *= a;
    return *this;
  }

  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  constexpr Rational& operator*=(const Rational& a) noexcept {
    numer *= a.numer;
    denom *= a.denom;
    return *this;
  }

  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  constexpr Rational& operator*=(const Index& a) noexcept {
    numer *= a;
    return *this;
  }
  
  /** Multiply by this
   * 
   * @param[in] a To multiply by
   * @return Rational& *this
   */
  constexpr Rational& operator*=(const int& a) noexcept {
    numer *= a;
    return *this;
  }

 /** Add one if possible */
  constexpr Rational operator++(int) const noexcept{
    return {numer + denom, denom};
  }

  /** Remove one if possible */
  constexpr Rational operator--(int) const noexcept {
    return {numer - denom, denom};
  }

 /** Add one if possible */
  constexpr Rational& operator++() noexcept {
    numer += denom;
    return *this;
  }

  /** Remove one if possible */
  constexpr Rational& operator--() noexcept {
    numer -= denom;
    return *this;
  }

  /** Cast to bool */
  explicit constexpr operator bool() const noexcept {
    return isDefined() and bool(numer);
  }

  /** Cast to Numeric */
  explicit constexpr operator Numeric() const noexcept { return toNumeric(); }

  /** Cast to Index */
  explicit constexpr operator Index() const noexcept { return toIndex(); }

  /** Cast to int */
  explicit constexpr operator int() const noexcept { return toInt(); }
  
  /** Binary read for Rational */
  bifstream& read(bifstream& bif) {
    bif >> numer >> denom;
    return bif;
  }
  
  /** Binary write for Rational */
  bofstream& write(bofstream& bof) const {
    bof << numer << denom;
    return bof;
  }
  
  /** Makes the sign of denom positive */
  constexpr Rational& fixSign() noexcept {
    if (denom < 0) {
      numer = -numer;
      denom = -denom;
    }
    return *this;
  }

  friend std::ostream& operator<<(std::ostream& os, const Rational& a);
  
  friend std::istream& operator>>(std::istream& is, Rational& a);
};  // Rational;

/** Returns the rational reduced by the greates 
 * 
 * @param[in] a Any Rational
 * @return a / gcd(a)
 */
constexpr Rational reduce_by_gcd(Rational a) noexcept {
  const Index div = gcd(a.numer, a.denom);
  if (div)
    return {a.numer / div, a.denom / div};
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
    const auto xi=Index(x);
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
  return signchange ? Rational{-nom, denom} : Rational{nom, denom};
}


// An undefined rational to be used everywhere for a undefined rationals
#define RATIONAL_UNDEFINED Rational(0, 0)

/** Negative
 * 
 * @param[in] a Any Rational
 * @return constexpr Rational Negative a
 */
constexpr Rational operator-(const Rational a) noexcept {
  return {-a.numer, a.denom};
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
  return (a.denom == b.denom)
             ? Rational(a.numer + b.numer, a.denom)
             : Rational(a.numer * b.denom + b.numer * a.denom,
                        a.denom * b.denom);
}

/** Addition
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(const Rational a, Index b) noexcept {
  return {a.numer + b * a.denom, a.denom};
}

/** Addition
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a + b
 */
constexpr Rational operator+(const Rational a, int b) noexcept {
  return {a.numer + b * a.denom, a.denom};
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
  return (a.denom == b.denom)
             ? Rational(a.numer - b.numer, a.denom)
             : Rational(a.numer * b.denom - b.numer * a.denom,
                        a.denom * b.denom);
}

/** Subtraction
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a - b
 */
constexpr Rational operator-(const Rational a, Index b) noexcept {
  return {a.numer - b * a.denom, a.denom};
}

/** Subtraction
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a - b
 */
constexpr Rational operator-(const Rational a, int b) noexcept {
  return {a.numer - b * a.denom, a.denom};
}

/** Subtraction
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b - a
 */
constexpr Rational operator-(Index b, const Rational a) noexcept {
  return {-a.numer + b * a.denom, a.denom};
}

/** Subtraction
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b - a
 */
constexpr Rational operator-(int b, const Rational a) noexcept {
  return {-a.numer + b * a.denom, a.denom};
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(const Rational a, const Rational b) noexcept {
  return {a.numer * b.denom, a.denom * b.numer};
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(const Rational a, Index b) noexcept {
  return {a.numer, a.denom * b};
}

/** Division
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a / b
 */
constexpr Rational operator/(const Rational a, int b) noexcept {
  return {a.numer, a.denom * b};
}

/** Division
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b / a
 */
constexpr Rational operator/(Index b, const Rational a) noexcept {
  return {a.denom * b, a.numer};
}

/** Division
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b / a
 */
constexpr Rational operator/(int b, const Rational a) noexcept {
  return {a.denom * b, a.numer};
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Rational
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(const Rational a, const Rational b) noexcept {
  return {a.numer * b.numer, a.denom * b.denom};
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(const Rational a, Index b) noexcept {
  return {a.numer * b, a.denom};
}

/** Multiplication
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a * b
 */
constexpr Rational operator*(const Rational a, int b) noexcept {
  return {a.numer * b, a.denom};
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
  return (a.denom == b.denom)
             ? Rational(a.numer % b.numer, a.denom)
             : Rational((a.numer * b.denom) % (a.denom * b.numer),
                        a.denom * b.denom);
}

/** Remainder
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a % b
 */
constexpr Rational operator%(const Rational a, Index b) noexcept {
  return {a.numer % (a.denom * b), a.denom};
}

/** Remainder
 * 
 * @param[in] a Any Rational
 * @param[in] b Any Index
 * @return constexpr Rational a % b
 */
constexpr Rational operator%(const Rational a, int b) noexcept {
  return {a.numer % (a.denom * b), a.denom};
}

/** Remainder
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b % a
 */
constexpr Rational operator%(Index b, const Rational a) noexcept {
  return {(b * a.denom) % a.numer, a.denom};
}

/** Remainder
 * 
 * @param[in] b Any Index
 * @param[in] a Any Rational
 * @return constexpr Rational b % a
 */
constexpr Rational operator%(int b, const Rational a) noexcept {
  return {(b * a.denom) % a.numer, a.denom};
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
         a.numer * b.denom == a.denom * b.numer;
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
         a.numer * b.denom < a.denom * b.numer;
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
 * @return true If a.numer and a.isDefined()
 * @return false Otherwise
 */
constexpr bool operator!(const Rational a) noexcept { return a.numer and a.isDefined(); }

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

using ArrayOfRational = Array<Rational>;

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
  return 0 == (r % 2);
}

/** Multiplication with numeric */
constexpr Numeric operator*(Rational y, Numeric x) noexcept {return y.toNumeric() * x;}
constexpr Numeric operator*(Numeric x, Rational y) noexcept {return x * y.toNumeric();}

/** Division with numeric */
constexpr Numeric operator/(Rational y, Numeric x) noexcept {return y.toNumeric() / x;}
constexpr Numeric operator/(Numeric x, Rational y) noexcept {return x / y.toNumeric();}

/** Addition with numeric */
constexpr Numeric operator+(Rational y, Numeric x) noexcept {return y.toNumeric() + x;}
constexpr Numeric operator+(Numeric x, Rational y) noexcept {return x + y.toNumeric();}

/** Subtraction with numeric */
constexpr Numeric operator-(Rational y, Numeric x) noexcept {return y.toNumeric() - x;}
constexpr Numeric operator-(Numeric x, Rational y) noexcept {return x - y.toNumeric();}

#endif  // rational_h


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

/** Contains the rational class definition
 * \file   rational.h
 * 
 * \author Richard Larsson
 * \date   2012-10-31
 **/

#ifndef rational_h
#define rational_h

#include <cassert>
#include <ostream>
#include "matpack.h"
#include "array.h"
#include "math_funcs.h"

class Rational
{
public:
    
    // Defining an object
    constexpr Rational() : mnom(0), mdenom(1) {}
    constexpr Rational(const Index n1) : mnom(n1), mdenom(1) {}
    constexpr Rational(const Index n1, const Index n2) : mnom(n1), mdenom(n2) {}
    
    // Reading values of object
    constexpr Index Nom() const {return mnom;}
    constexpr Index Denom() const {return mdenom;}

    // Status of object
    constexpr bool isUndefined() const { return (mdenom == 0); }
    constexpr bool isDefined() const { return mdenom not_eq 0; }
    constexpr bool isIndex() const {return isDefined() and not bool(mnom%mdenom);}

    // Converting object
    constexpr Index toIndex() const {return isIndex() ? mnom/mdenom : throw std::logic_error("Rational is not also Index");}
    constexpr Numeric toNumeric() const {return Numeric(mnom)/Numeric(mdenom); }
    constexpr int toInt() const {return int(toIndex());}

    // Useful. Keep this around and/or improve.
    void Simplify();
    
    // Assignment operators
    Rational& operator+=(const Rational& a) {mnom = mnom*a.Denom() + a.Nom()*mdenom;mdenom *= a.Denom(); return *this;}
    Rational& operator+=(const Index& a)    {mnom+= mdenom*a; return *this;}
    Rational& operator-=(const Rational& a) {mnom = mnom*a.Denom() - a.Nom()*mdenom;mdenom *= a.Denom(); return *this;}
    Rational& operator-=(const Index& a)    {mnom-= mdenom*a; return *this;}
    Rational& operator/=(const Rational& a) {mnom*=a.Denom();mdenom*=a.Nom(); return *this;}
    Rational& operator/=(const Index& a)    {mdenom*=a; return *this;}
    Rational& operator*=(const Rational& a) {mnom*=a.Nom();mdenom*=a.Denom(); return *this;}
    Rational& operator*=(const Index& a)    {mnom*=a; return *this;}

    // Iterative operators
    constexpr Rational operator++(int) const {return isDefined() ? Rational(mnom+mdenom, mdenom) : throw std::logic_error("Undefined Rational in operator++");}
    constexpr Rational operator--(int) const {return isDefined() ? Rational(mnom-mdenom, mdenom) : throw std::logic_error("Undefined Rational in operator--");}
    Rational& operator++() {mnom += mdenom; return *this;}
    Rational& operator--() {mnom -= mdenom; return *this;}
    
    // Boolean
    explicit constexpr operator bool() const {return bool(mnom) and bool(mdenom);}
    
private:
    // Rational is supposed to be used rationally ( mnom / mdenom )
    Index mnom;
    Index mdenom;
};

#define RATIONAL_UNDEFINED Rational(0, 0)

// Arithmetic operations
constexpr Rational operator-(Rational a)
{
  return a.isDefined() ? 
         Rational(-a.Nom(), a.Denom()) :
  throw std::logic_error("Undefined Rational in operator-");
}

constexpr Rational operator+(Rational a) 
{
  return a.isDefined() ? 
         a :
  throw std::logic_error("Undefined Rational in operator+");
}

constexpr Rational operator+(Rational a, Rational b)
{
  return (a.isDefined() and b.isDefined()) ?
         Rational(a.Nom()*b.Denom() + b.Nom()*a.Denom(), a.Denom() * b.Denom()) :
  throw std::logic_error("Undefined Rational in operator+");
}

constexpr Rational operator+(Rational a, Index b) 
{
  return a.isDefined() ? 
         Rational(a.Nom()+b*a.Denom(),a.Denom()) :
  throw std::logic_error("Undefined Rational in operator+");
}

constexpr Rational operator+(Index b, Rational a) {return operator+(a, b);}

constexpr Rational operator-(Rational a, Rational b)
{
  return (a.isDefined() and b.isDefined()) ? 
         Rational(a.Nom()*b.Denom() - b.Nom()*a.Denom(), a.Denom() * b.Denom()) :
  throw std::logic_error("Undefined Rational in operator-");
}

constexpr Rational operator-(Rational a, Index b)
{
  return a.isDefined() ?
         Rational( a.Nom()-b*a.Denom(), a.Denom()) :
  throw std::logic_error("Undefined Rational in operator-");
}

constexpr Rational operator-(Index b, Rational a)
{
  return a.isDefined() ?
         Rational(-a.Nom()+b*a.Denom(), a.Denom()) :
  throw std::logic_error("Undefined Rational in operator-");
}

constexpr Rational operator/(Rational a, Rational b)
{
  return (a.isDefined() and b.isDefined()) ?
         Rational(a.Nom()*b.Denom(),a.Denom() * b.Nom()) :
  (throw std::logic_error("Undefined Rational in operator/"));
}

constexpr Rational operator/(Rational a, Index b)
{
  return a.isDefined() ?
         Rational(a.Nom(),a.Denom()*b) :
  throw std::logic_error("Undefined Rational in operator/");
}

constexpr Rational operator/(Index b, Rational a)
{
  return a.isDefined() ?
         Rational(a.Denom()*b,a.Nom()) :
  throw std::logic_error("Undefined Rational in operator/");
}

constexpr Rational operator*(Rational a, Rational b)
{
  return (a.isDefined() and b.isDefined()) ?
         Rational(a.Nom()*b.Nom(),a.Denom() * b.Denom()) :
  throw std::logic_error("Undefined Rational in operator*");
}

constexpr Rational operator*(Rational a, Index b)
{
  return a.isDefined() ?
         Rational(a.Nom()*b,a.Denom()) :
  throw std::logic_error("Undefined Rational in operator*");
}

constexpr Rational operator*(Index b, Rational a) {return operator*(a, b);}

constexpr Rational operator%(Rational a, Rational b)
{
  return (a.isDefined() and b.isDefined()) ?
         Rational((a.Nom()*b.Denom())%(a.Denom()*b.Nom()), a.Denom()*b.Denom()) :
  throw std::logic_error("Undefined Rational in operator%");
}

constexpr Rational operator%(Rational a, Index b)
{
  return a.isDefined() ?
         Rational(a.Nom()%b,a.Denom()) :
  throw std::logic_error("Undefined Rational in operator%");
}

constexpr Rational operator%(Index b, Rational a)
{
  return a.isDefined() ? 
         Rational((b*a.Denom())%a.Nom(), a.Denom()) :
  throw std::logic_error("Undefined Rational in operator%");
}

// Boolean operations
constexpr bool operator==(Rational a, Rational b) {return a.isDefined() and b.isDefined() and a.Nom()*b.Denom()==a.Denom()*b.Nom();}
constexpr bool operator!=(Rational a, Rational b) {return not operator==(a, b);}
constexpr bool operator<(Rational a, Rational b)  {return a.isDefined() and b.isDefined() and a.Nom()*b.Denom()<a.Denom()*b.Nom();}
constexpr bool operator>(Rational a, Rational b)  {return operator<(b, a);}
constexpr bool operator<=(Rational a, Rational b) {return not operator>(a, b);}
constexpr bool operator>=(Rational a, Rational b) {return not operator<(a, b);}
constexpr bool operator!(Rational a) {return a.Nom() and a.isDefined();}

inline Numeric fac(const Rational& r) { return (::fac(r.toIndex())); }

inline Numeric sqrt(const Rational& r) {return std::sqrt(r.toNumeric());}

std::ostream& operator<<(std::ostream& os, const Rational& a);

std::istream& operator>>(std::istream& os, Rational& a);

Rational abs(const Rational& a);

typedef Array<Rational> ArrayOfRational;

#endif // rational_h

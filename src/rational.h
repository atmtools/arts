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
    Rational() : mnom(0), mdenom(1) {}
    Rational(const Index& n1) : mnom(n1), mdenom(1) {}
    Rational(const Index& n1, const Index& n2) : mnom(n1), mdenom(n2) {}
    
    // Reading values of object
    Index Nom() const {return mnom;}
    Index Denom() const {return mdenom;}

    bool isUndefined() const { return (mnom == 0 && mdenom == 0); }

    // Converting object
    Index toIndex() const;
    Numeric toNumeric() const { return (Numeric)mnom/(Numeric)mdenom; }

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
    Rational  operator++(int) {mnom += mdenom; return *this;}
    Rational  operator--(int) {mnom -= mdenom; return *this;}
    Rational& operator++()    {mnom += mdenom; return *this;}
    Rational& operator--()    {mnom -= mdenom; return *this;}
    
    // Boolean
    explicit operator bool() const {return bool(mnom);}
    
private:
    // Rational is supposed to be used rationally ( mnom / mdenom )
    Index mnom;
    Index mdenom;
};

#define RATIONAL_UNDEFINED Rational(0, 0)

// Arithmetic operations
inline Rational operator-(const Rational& a) {return Rational(-a.Nom(), a.Denom());}
inline Rational operator+(const Rational& a) {return a;}
inline Rational operator+(const Rational& a, const Rational& b) {return Rational(a.Nom()*b.Denom() + b.Nom()*a.Denom(),a.Denom() * b.Denom());}
inline Rational operator+(const Rational& a, const Index& b) {return Rational(a.Nom()+b*a.Denom(),a.Denom());}
inline Rational operator+(const Index& b, const Rational& a) {return Rational(a.Nom()+b*a.Denom(),a.Denom());}
inline Rational operator-(const Rational& a, const Rational& b) {return Rational(a.Nom()*b.Denom() - b.Nom()*a.Denom(),a.Denom() * b.Denom());}
inline Rational operator-(const Rational& a, const Index& b) {return Rational(a.Nom()-b*a.Denom(),a.Denom());}
inline Rational operator-(const Index& b, const Rational& a) {return Rational(-a.Nom()+b*a.Denom(),a.Denom());}
inline Rational operator/(const Rational& a, const Rational& b) {return Rational(a.Nom()*b.Denom(),a.Denom() * b.Nom());}
inline Rational operator/(const Rational& a, const Index& b) {return Rational(a.Nom(),a.Denom()*b);}
inline Rational operator/(const Index& b, const Rational& a) {return Rational(a.Denom()*b,a.Nom());}
inline Rational operator*(const Rational& a, const Rational& b) {return Rational(a.Nom()*b.Nom(),a.Denom() * b.Denom());}
inline Rational operator*(const Rational& a, const Index& b) {return Rational(a.Nom()*b,a.Denom());}
inline Rational operator*(const Index& b, const Rational& a) {return Rational(a.Nom()*b,a.Denom());}
inline Rational operator%(const Rational& a, const Rational& b) {return Rational((a.Nom()*b.Denom())%(a.Denom()*b.Nom()),a.Denom()*b.Denom());}
inline Rational operator%(const Rational& a, const Index& b) {return Rational(a.Nom()%b,a.Denom());}

// Boolean operations
inline bool operator==(const Rational& a, const Rational& b) {return a.Denom()!=0 && b.Denom()!=0 && a.Nom()*b.Denom()==a.Denom()*b.Nom();}
inline bool operator!=(const Rational& a, const Rational& b) {return !operator==(a, b);}
inline bool operator<(const Rational& a, const Rational& b)  {return a.Denom()!=0 && b.Denom()!=0 && a.Nom()*b.Denom()<a.Denom()*b.Nom();}
inline bool operator>(const Rational& a, const Rational& b)  {return operator<(b, a);}
inline bool operator<=(const Rational& a, const Rational& b) {return !operator>(a, b);}
inline bool operator>=(const Rational& a, const Rational& b) {return !operator<(a, b);}
inline bool operator!(const Rational& a) {return a.Nom()==0 && a.Denom()!=0;}

inline Numeric fac(const Rational& r) { return (::fac(r.toIndex())); }

std::ostream& operator<<(std::ostream& os, const Rational& a);

std::istream& operator>>(std::istream& os, Rational& a);

Rational abs(const Rational& a);

typedef Array<Rational> ArrayOfRational;

#endif // rational_h

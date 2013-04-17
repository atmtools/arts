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

class Rational
{
public:
    
    // Defining an object
    Rational(){mnom=0;mdenom=1;}
    Rational(const Index& n1){mnom=n1;mdenom=1;}
    Rational(const Index& n1, const Index& n2){mnom=n1;mdenom=n2;Simplify();}
    
    // Reading values of object
    Index Nom() const {return mnom;}
    Index Denom() const {return mdenom;}

    bool isUndefined() const { return (mnom == 0 && mdenom == 0); }

    // Converting object
    operator Index() const{Rational b(*this); b.Simplify();assert(b.Denom()==1);return b.Nom();}
    operator Numeric() const{return (Numeric)mnom/(Numeric)mdenom;}
    operator bool() const{return mnom==0?false:true;}
    
    // Useful. Keep this around and/or improve.
    void Simplify();
    
    // Assigning values to object.
    void operator=(const Rational& a){mnom = a.Nom();mdenom = a.Denom();Simplify();}
    void operator=(const Index& a){mnom = a;mdenom = 1;}
    void operator=(const int& a){mnom = a;mdenom = 1;}
    
    // Iterative operators
    Rational  operator++(int){mnom += mdenom; return *this;}
    Rational  operator--(int){mnom -= mdenom; return *this;}
    Rational& operator++()   {mnom += mdenom; return *this;}
    Rational& operator--()   {mnom -= mdenom; return *this;}
    
    // Changing of signs
    Rational operator-() const{return Rational(-mnom,mdenom);}
    Rational operator+() const{return *this;}
    
    // Changing the object via Index
    Rational operator+=(const Index& a){mnom+=a*mdenom;Simplify(); return *this;}
    Rational operator-=(const Index& a){mnom-=a*mdenom;Simplify(); return *this;}
    Rational operator*=(const Index& a){mnom*=a;Simplify(); return *this;}
    Rational operator/=(const Index& a){mdenom*=a;Simplify(); return *this;}
    
    // Changing the object via other Rational
    Rational operator+=(const Rational& a){mnom = mnom*a.Denom() + a.Nom()*mdenom;mdenom *= a.Denom();Simplify(); return *this;}
    Rational operator-=(const Rational& a){mnom = mnom*a.Denom() - a.Nom()*mdenom;mdenom *= a.Denom();Simplify(); return *this;}
    Rational operator/=(const Rational& a){mnom*=a.Denom();mdenom*=a.Nom();Simplify(); return *this;}
    Rational operator*=(const Rational& a){mnom*=a.Nom();mdenom*=a.Denom();Simplify(); return *this;}
    
    // Changing the object via int
    Rational operator+=(const int& a){mnom+=a*mdenom;Simplify(); return *this;}
    Rational operator-=(const int& a){mnom-=a*mdenom;Simplify(); return *this;}
    Rational operator*=(const int& a){mnom*=a;Simplify(); return *this;}
    Rational operator/=(const int& a){mdenom*=a;Simplify(); return *this;}
    
    // Arithmetic operations involving Index
    Rational operator+(const Index& a) const{return Rational(mnom+a*mdenom,mdenom);}
    Rational operator-(const Index& a) const{return Rational(mnom-a*mdenom,mdenom);}
    Rational operator*(const Index& a) const{return Rational(mnom*a,mdenom);}
    Rational operator/(const Index& a) const{return Rational(mnom,mdenom*a);}
    Rational operator%(const Index& a) const{return Rational(mnom%(a*mdenom),mdenom);}
    
    // Arithmetic operations involving int
    Rational operator+(const int& a) const{return Rational(mnom+a*mdenom,mdenom);}
    Rational operator-(const int& a) const{return Rational(mnom-a*mdenom,mdenom);}
    Rational operator*(const int& a) const{return Rational(mnom*a,mdenom);}
    Rational operator/(const int& a) const{return Rational(mnom,mdenom*a);}
    Rational operator%(const int& a) const{return Rational(mnom%(a*mdenom),mdenom);}
    
    // Arithmetic operations involving other Rational
    Rational operator+(const Rational& a) const{return Rational(mnom*a.Denom() + a.Nom()*mdenom,mdenom * a.Denom());}
    Rational operator-(const Rational& a) const{return Rational(mnom*a.Denom() - a.Nom()*mdenom,mdenom * a.Denom());}
    Rational operator/(const Rational& a) const{return Rational(mnom*a.Denom(),mdenom * a.Nom());}
    Rational operator*(const Rational& a) const{return Rational(mnom*a.Nom(),mdenom * a.Denom());}
    Rational operator%(const Rational& a) const{return Rational((mnom*a.Denom())%(mdenom*a.Nom()),mdenom*a.Denom());}
        
    // Boolean operations involving Index
    bool  operator<(const Index& a) const{return mnom<a*mdenom;}
    bool  operator>(const Index& a) const{return mnom>a*mdenom;}
    bool operator<=(const Index& a) const{return mnom<=a*mdenom;}
    bool operator>=(const Index& a) const{return mnom>=a*mdenom;}
    bool operator==(const Index& a) const{Rational b(*this); b.Simplify(); return b.Denom()!=1?false:b.Nom()==a;}
    bool operator!=(const Index& a) const{Rational b(*this); b.Simplify(); return b.Denom()!=1?false:b.Nom()!=a;}
    
    // Boolean operations involving int
    bool operator<(const int& a)  const{return mnom<a*mdenom;}
    bool operator>(const int& a)  const{return mnom>a*mdenom;}
    bool operator<=(const int& a) const{return mnom<=a*mdenom;}
    bool operator>=(const int& a) const{return mnom>=a*mdenom;}
    bool operator==(const int& a) const{Rational b(*this); b.Simplify(); return b.Denom()!=1?false:b.Nom()==a;}
    bool operator!=(const int& a) const{Rational b(*this); b.Simplify(); return b.Denom()!=1?false:b.Nom()!=a;}
    
    // Boolean operations involving other Rational
    bool operator<(const Rational& a)  const{return mnom*a.Denom()<mdenom*a.Nom();}
    bool operator>(const Rational& a)  const{return mnom*a.Denom()>mdenom*a.Nom();}
    bool operator<=(const Rational& a) const{return mnom*a.Denom()<=mdenom*a.Nom();}
    bool operator>=(const Rational& a) const{return mnom*a.Denom()>=mdenom*a.Nom();}
    bool operator==(const Rational& a) const{return mnom*a.Denom()==mdenom*a.Nom();}
    bool operator!=(const Rational& a) const{return mnom*a.Denom()!=mdenom*a.Nom();}

private:
    // Rational is supposed to be used rationally ( mnom / mdenom )
    Index mnom;
    Index mdenom;
};

#define RATIONAL_UNDEFINED Rational(0, 0)

ostream& operator<<(ostream& os, const Rational& a);

Rational abs(const Rational& a);

typedef Array<Rational> ArrayOfRational;

#endif // rational_h

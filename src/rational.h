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

#include "matpack.h"
#include <cassert>
#include <ostream>

class Rational
{
public:
    
    // Defining an object
    Rational(){top=0;low=1;}
    Rational(const Index& n1){top=n1;low=1;}
    Rational(const Index& n1, const Index& n2){assert(n2>0);top=n1;low=n2;Simplify();}
    
    // Reading values of object
    Index Top() const{return top;}
    Index Low() const {return low;}
    
    // Converting object
    operator Index() const{Rational b(*this); b.Simplify();assert(b.Low()==1);return b.Top();}
    operator Numeric() const{return (Numeric)top/(Numeric)low;}
    operator bool() const{return top==0?false:true;}
    
    // Useful. Keep this around and/or improve.
    void Simplify()
    {
        Index ii = low;
        
        while( ii != 0 )
        {
            if( top % ii == 0 && low % ii == 0 )
            {
                top /= ii;
                low /= ii;
                
                if( top == 1 || low == 1 || top == 0 || top == -1 )
                    break;
            }
            if( low < ii )
                ii = low;
            else
                ii--;
        }
    }
    
    // Assigning values to object.
    void operator=(const Rational& a){top = a.Top();low = a.Low();Simplify();}
    void operator=(const Index& a){top = a;low = 1;}
    void operator=(const int& a){top = a;low = 1;}
    
    // Iterative operators
    Rational  operator++(int){top += low; return *this;}
    Rational  operator--(int){top -= low; return *this;}
    Rational& operator++()   {top += low; return *this;}
    Rational& operator--()   {top -= low; return *this;}
    
    // Changing of signs
    Rational operator-() const{return Rational(-top,low);}
    Rational operator+() const{return *this;}
    
    // Changing the object via Index
    Rational operator+=(const Index& a){top+=a*low;Simplify(); return *this;}
    Rational operator-=(const Index& a){top-=a*low;Simplify(); return *this;}
    Rational operator*=(const Index& a){top*=a;Simplify(); return *this;}
    Rational operator/=(const Index& a){low*=a;Simplify(); return *this;}
    
    // Changing the object via other Rational
    Rational operator+=(const Rational& a){top = top*a.Low() + a.Top()*low;low *= a.Low();Simplify(); return *this;}
    Rational operator-=(const Rational& a){top = top*a.Low() - a.Top()*low;low *= a.Low();Simplify(); return *this;}
    Rational operator/=(const Rational& a){top*=a.Low();low*=a.Top();Simplify(); return *this;}
    Rational operator*=(const Rational& a){top*=a.Top();low*=a.Low();Simplify(); return *this;}
    
    // Changing the object via int
    Rational operator+=(const int& a){top+=a*low;Simplify(); return *this;}
    Rational operator-=(const int& a){top-=a*low;Simplify(); return *this;}
    Rational operator*=(const int& a){top*=a;Simplify(); return *this;}
    Rational operator/=(const int& a){low*=a;Simplify(); return *this;}
    
    // Arithmetic operations involving Index
    Rational operator+(const Index& a) const{return Rational(top+a*low,low);}
    Rational operator-(const Index& a) const{return Rational(top-a*low,low);}
    Rational operator*(const Index& a) const{return Rational(top*a,low);}
    Rational operator/(const Index& a) const{return Rational(top,low*a);}
    Rational operator%(const Index& a) const{return Rational(top%(a*low),low);}
    
    // Arithmetic operations involving int
    Rational operator+(const int& a) const{return Rational(top+a*low,low);}
    Rational operator-(const int& a) const{return Rational(top-a*low,low);}
    Rational operator*(const int& a) const{return Rational(top*a,low);}
    Rational operator/(const int& a) const{return Rational(top,low*a);}
    Rational operator%(const int& a) const{return Rational(top%(a*low),low);}
    
    // Arithmetic operations involving other Rational
    Rational operator+(const Rational& a) const{return Rational(top*a.Low() + a.Top()*low,low * a.Low());}
    Rational operator-(const Rational& a) const{return Rational(top*a.Low() - a.Top()*low,low * a.Low());}
    Rational operator/(const Rational& a) const{return Rational(top*a.Low(),low * a.Top());}
    Rational operator*(const Rational& a) const{return Rational(top*a.Top(),low * a.Low());}
    Rational operator%(const Rational& a) const{return Rational((top*a.Low())%(low*a.Top()),low*a.Low());}
        
    // Boolean operations involving Index
    bool  operator<(const Index& a) const{return top<a*low;}
    bool  operator>(const Index& a) const{return top>a*low;}
    bool operator<=(const Index& a) const{return top<=a*low;}
    bool operator>=(const Index& a) const{return top>=a*low;}
    bool operator==(const Index& a) const{Rational b(*this); b.Simplify(); return b.Low()!=1?false:b.Top()==a;}
    bool operator!=(const Index& a) const{Rational b(*this); b.Simplify(); return b.Low()!=1?false:b.Top()!=a;}
    
    // Boolean operations involving int
    bool operator<(const int& a)  const{return top<a*low;}
    bool operator>(const int& a)  const{return top>a*low;}
    bool operator<=(const int& a) const{return top<=a*low;}
    bool operator>=(const int& a) const{return top>=a*low;}
    bool operator==(const int& a) const{Rational b(*this); b.Simplify(); return b.Low()!=1?false:b.Top()==a;}
    bool operator!=(const int& a) const{Rational b(*this); b.Simplify(); return b.Low()!=1?false:b.Top()!=a;}
    
    // Boolean operations involving other Rational
    bool operator<(const Rational& a)  const{return top*a.Low()<low*a.Top();}
    bool operator>(const Rational& a)  const{return top*a.Low()>low*a.Top();}
    bool operator<=(const Rational& a) const{return top*a.Low()<=low*a.Top();}
    bool operator>=(const Rational& a) const{return top*a.Low()>=low*a.Top();}
    bool operator==(const Rational& a) const{return top*a.Low()==low*a.Top();}
    bool operator!=(const Rational& a) const{return top*a.Low()!=low*a.Top();}
    
private:
    // Rational is supposed to be used rationally ( top / low )
    Index top;
    Index low;
};

ostream& operator<<(ostream& os, const Rational& a);

Rational abs(const Rational& a);

#endif // rational_h
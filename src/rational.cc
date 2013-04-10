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

/** Contains some additional functionality of the rational class
   \file   rational.cc
   
   \author Richard Larsson
   \date   2012-10-31
**/

#include "rational.h"

void Rational::Simplify()
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

ostream& operator<<(ostream& os, const Rational& a)
{
    os << a.Top() << "/" << a.Low();
    return os;
}

Rational abs(const Rational& a)
{
    return a<0?-a:a;
}

/* Copyright (C) 2002 Stefan Buehler <sbuehler@uni-bremen.de>

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

#include "array.h"
#include "matpackVII.h"
#include "interpolation.h"

/** Define the global joker objekt. */
Joker joker;

void test01()
{
  cout << "Simple 1D interpolation case\n"
       << "----------------------------\n";
  //  Vector og(5,5,-1);		// 5,4,3,2,1
  Vector og(1,5,+1);		// 1, 2, 3, 4, 5
  Vector ng(1,5,0.25);		// 1.0, 1,25, 1.5, 1.75, 2.0

  cout << "og:\n" << og << "\n";
  cout << "ng:\n" << ng << "\n";

  // To store the grid positions:
  ArrayOfGridPos gp(ng.nelem());

  gridpos(gp,og,ng);
  cout << "gp:\n" << gp << "\n";

  // To store interpolation weights:
  Matrix itw(gp.nelem(),2);
  interpweights(itw,gp);
  cout << "itw:\n" << itw << "\n";
}

void test02(Index n)
{
  cout << "Test whether for loop or iterator are quicker\n"
       << "a) for loop\n"
       << "---------------------------------------------\n";

  Vector a(n);
  for (Index i=0; i<a.nelem(); ++i)
    a[i] = i;
}

void test03(Index n)
{
  cout << "Test whether for loop or iterator are quicker\n"
       << "b) iterator\n"
       << "---------------------------------------------\n";

  Vector a(n);
  Iterator1D ai=a.begin();
  const Iterator1D ae=a.end();
  Index i=0;
  for ( ; ai!=ae; ++ai, ++i )
    *ai = i;
}

// Result: Both are almost equally fast, with a slight advantage of
// the for loop if compiler optimization is enabled.

int main()
{
  test01();
}

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

#include <iostream>
#include "array.h"
#include "matpackVII.h"
#include "interpolation.h"

/** Define the global joker objekt. */
Joker joker;

void test01()
{
  cout << "Simple interpolation cases\n"
       << "--------------------------\n";
  //  Vector og(5,5,-1);		// 5,4,3,2,1
  Vector og(1,5,+1);		// 1, 2, 3, 4, 5
  Vector ng(2,5,0.25);		// 2.0, 2,25, 2.5, 2.75, 3.0

  cout << "og:\n" << og << "\n";
  cout << "ng:\n" << ng << "\n";

  // To store the grid positions:
  ArrayOfGridPos gp(ng.nelem());

  gridpos(gp,og,ng);
  cout << "gp:\n" << gp << "\n";

  cout << "1D:\n"
       << "---\n";
  {
    // To store interpolation weights:
    Matrix itw(gp.nelem(),2);
    interpweights(itw,gp);
    
    cout << "itw:\n" << itw << "\n";

    // Original field:
    Vector of(og.nelem(),0);
    of[2] = 10;			// 0, 0, 10, 0, 0

    cout << "of:\n" << of << "\n";

    // Interpolated field:
    Vector nf(ng.nelem());

    interp(nf, itw, of, gp);

    cout << "nf:\n" << nf << "\n";
  }

  cout << "Blue 2D:\n"
       << "--------\n";
  {
    // To store interpolation weights:
    Matrix itw(gp.nelem(),4);
    interpweights(itw,gp,gp);
    
    cout << "itw:\n" << itw << "\n";

    // Original field:
    Matrix of(og.nelem(),og.nelem(),0);
    of(2,2) = 10;			// 0 Matrix with 10 in the middle

    cout << "of:\n" << of << "\n";

    // Interpolated field:
    Vector nf(ng.nelem());

    interp(nf, itw, of, gp, gp);

    cout << "nf:\n" << nf << "\n";
  }

  cout << "Blue 6D:\n"
       << "--------\n";
  {
    // To store interpolation weights:
    Matrix itw(gp.nelem(),64);
    interpweights(itw,gp,gp,gp,gp,gp,gp);
    
    //    cout << "itw:\n" << itw << "\n";

    // Original field:
    Index n = og.nelem();
    Tensor6 of(n,n,n,n,n,n,0);
    of(2,2,2,2,2,2) = 10;			// 0 Tensor with 10 in the middle

    //    cout << "of:\n" << of << "\n";

    // Interpolated field:
    Vector nf(ng.nelem());

    interp(nf, itw, of, gp, gp, gp, gp, gp, gp);

    cout << "nf:\n" << nf << "\n";
  }

  cout << "Green 2D:\n"
       << "---------\n";
  {
    // To store interpolation weights:
    Tensor3 itw(gp.nelem(),gp.nelem(),4);
    interpweights(itw,gp,gp);
    
    for ( Index i=0; i<itw.ncols(); ++i )
      cout << "itw " << i << ":\n" << itw(Range(joker),Range(joker),i) << "\n";

    // Original field:
    Matrix of(og.nelem(),og.nelem(),0);
    of(2,2) = 10;			// 0 Matrix with 10 in the middle

    cout << "of:\n" << of << "\n";

    // Interpolated field:
    Matrix nf(ng.nelem(),ng.nelem());

    interp(nf, itw, of, gp, gp);

    cout << "nf:\n" << nf << "\n";
  }
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

/* Copyright (C) 2003 Mattias Ekström <ekstrom@rss.chalmers.se>

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

#include "sensor.h"
//#include <cmath>
//#include <iostream>
#include "matpackI.h"
//#include "array.h"
//#include "make_array.h"
//#include "mystring.h"
//#include "make_vector.h"
//#include "math_funcs.h"
//#include "describe.h"
//#include "matpackII.cc"

void test1()
{
  //Test antenna_transfer_matrix with sparse matrix
  Sparse H(3,15);
  Vector m_za(-8,5,4);
  Vector a_grid(-10,21,1);
  Vector a(a_grid.nelem());
  Vector f_grid(2,3,2);

  antenna_diagram_gaussian(a, a_grid, 2);

  antenna_transfer_matrix( H, m_za, a, a_grid, f_grid );

  cout << "H:\n" << H << "\n";
  cout << "m_za:\n" << m_za << "\n";
  cout << "a:\n" << a << "\n";
  cout << "a_grid:\n" << a_grid << "\n";
  cout << "f_grid:\n" << f_grid << "\n";
}

int main()
{
  test1();

  return 0;
}

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
#include <iostream>
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
  cout << "\nTest 1:\n";

  Sparse H(3,15);
  Vector m_za(-8,5,4);
  Matrix srm(21,2);
  Vector a_grid(-10,21,1);
  srm(Range(joker),0) = a_grid;
  Vector f_grid(2,3,2);

  antenna_diagram_gaussian(srm, 2);

  antenna_transfer_matrix( H, m_za, srm, f_grid );

  cout << "H:\n" << H << "\n";
  cout << "m_za:\n" << m_za << "\n";
  cout << "srm:\n" << srm << "\n";
  cout << "f_grid:\n" << f_grid << "\n";
}

void test2()
{
  //Test spectrometer_transfer_matrix
  cout << "\nTest 2:\n";

  Sparse H(4,5);
  Matrix srm(7,2);
  Vector x_srm(-3,7,1);
//  Vector values_srm(x_srm.nelem());
//  antenna_diagram_gaussian(values_srm, x_srm, 2);
  srm(Range(joker),0) = x_srm;
  antenna_diagram_gaussian(srm, 2);
//  srm(Range(joker),2) = values_srm;
  Vector f_grid(1,5,2);
  Vector cf_grid(2,4,2);

  spectrometer_transfer_matrix( H, srm, cf_grid, f_grid);

  cout << "H:\n" << H << "\n";
  cout << "srm:\n" << srm << "\n";
  cout << "cf_grid:\n" << cf_grid << "\n";
  cout << "f_grid:\n" << f_grid << "\n";
}

int main()
{
  test1();
  test2();

  return 0;
}

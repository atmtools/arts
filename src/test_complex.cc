/* Copyright (C) 2002-2012 Stefan Buehler <sbuehler@ltu.se>

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

/*!
  \file   test_complex.cc
  \author  <sbuehler@ltu.se>
  \date   Sat Dec 14 19:42:25 2002
  
  \brief  Test the complex numbers.
*/

#include <iostream>
#include "complex.h"

using std::cout;


void test01()
{
  Complex a;
  Complex b (3., 0.);

  //cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";

  a = b;
  cout << "a = " << a << "\n";

  
  a = Complex (1., 1.);
  cout << "a = " << a << "\n";
  cout << "a.abs() = " << abs (a) << "\n";
  cout << "a.arg()= "<< arg (a)*57.2957914331333 <<" degree"<< "\n";
  

  Complex c;
  c = a + b;
  cout << "c = " << c << "\n";

  cout << "a += b: " << (a += b) << "\n";

  cout << "c + 3 = " << (c + 3.) << "\n";

  cout << "c + 3i = c + Complex (0,3) = " << (c + Complex (0., 3.)) << "\n";

   Complex d;
   a = Complex (1., 1.);
   b = Complex (-1., -1.);
   d = a * b;
   cout << "d = " << d << "\n";
   Complex e;
   e = a / b;
   cout << "e = " << e << "\n";
   Complex f;
   f = pow(a,0)+pow(b,0);
   cout << "f = " << f << "\n";
   Complex g;
   g = pow(a,1)+pow(b,1);
   cout << "g = " << g << "\n";
   Complex h;
   h=pow(a,2)+pow(b,2)+pow(a,3)+pow(b,3);
   cout << "h = " << h << "\n";
   
   ComplexMatrix A(4, 4);
   ComplexMatrix B(4, 4);
   ComplexMatrix X1(4, 1);
   ComplexVector X2(4);
   ComplexMatrix C(4, 4);
   for(Index i = 0; i < 4; i++)
   {
     X1(i, 0) = Complex((Numeric)i * 5.0 + 2.0, (Numeric)i + 1.0);
     X2[i] = Complex((Numeric)i * 5.0 + 2.0, (Numeric)i + 1.0);
     for(Index j = 0; j < 4; j++)
     {
       A(i, j) = Complex((Numeric)i + (Numeric)j, 0.0) + Complex(1.0, 1.0);
       B(i, j) = Complex(2.0 * (Numeric)i + 4.0 * (Numeric)j, 0.0) + Complex(3.0, 3.0);
     }
   }
     
   mult(C, A, B);
   std::cout<<C<<"\n";
   C = 0;
   mult(C(ComplexRange(2,2), ComplexRange(2,2)), A(ComplexRange(0,2), ComplexRange(0,2)), B(ComplexRange(1,2), ComplexRange(1,2)));
   std::cout<<C<<"\n";
   C = 0;
   mult(C(ComplexRange(0,4), ComplexRange(2,1)), A, X1);
   mult(C(joker, 1), A, X2);
   std::cout<<C<<"\n";
}

int main()
{
  test01();
  return 0;
}

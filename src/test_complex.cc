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

/*!
  \file   test_complex.cc
  \author  <sbuehler@uni-bremen.de>
  \date   Sat Dec 14 19:42:25 2002
  
  \brief  Test the complex numbers.
*/

#include "complex.h"

void test01()
{
  Complex a;
  Complex b(3,4);

  //cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";

  a = b;
  cout << "a = " << a << "\n";

  a.Re() = 2;
  a.Im() = 5;
  cout << "a = " << a << "\n";

  // cout << "b.Abs() = " << b.Abs() << "\n";
  // cout << "b.Pha() = " << b.Pha() << "\n";
  Complex c;
  c = a + b;
  cout << "c = " << c << "\n";

  cout << "a+=b: " << (a+=b) << "\n";

  cout << "c+3 = " << (c+3) << "\n";

  cout << "c+3i = c+Complex(0,3) = " << (c+Complex(0,3)) << "\n";
   Complex d;
   c = a * b;
   cout << "d = " << d << "\n";
   Complex e;
   e = a * b;
   cout << "e = " << e << "\n";
}

int main()
{
  test01();
  return 0;
}

/* Copyright (C) 2001 Stefan Buehler <sbuehler@uni-bremen.de>

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

#include "matpackI.h"
#include "matpackIII.h"
#include "matpackIV.h"

void test1()
{
  cout << "Test Tensor4:\n";

  Tensor4 a(2, 3, 4, 5);

  Index fill = 0;

  for (Index i = 0; i < a.nbooks(); ++i)
    for (Index j = 0; j < a.npages(); ++j)
      for (Index k = 0; k < a.nrows(); ++k)
        for (Index l = 0; l < a.ncols(); ++l)
          a(i, j, k, l) = ++fill;

  cout << "a =" << endl << a << endl;
}

int main()
{
  test1();
  return 0;
}

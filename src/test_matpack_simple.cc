/* Copyright (C) 2019 Oliver Lemke  <oliver.lemke@uni-hamburg.de>

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

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include "matpack_data.h"

// Simple element access operator benchmark
void test1() {
  Matrix m(500, 500);
  for (Index k = 0; k < 1000; k++)
    for (Index i = 0; i < m.nrows(); i++)
      for (Index j = 0; j < m.ncols(); j++) m(i, j) = 2.;
}

int main() {
  test1();

  return 0;
}

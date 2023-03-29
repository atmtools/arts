/* Copyright (C) 2002-2012
   Stefan Buehler <sbuehler@ltu.se>
   Wolfram-Andre Haas <wolhaas@hermes.fho-emden.de>

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

#include <cmath>
#include <iostream>

#include "array.h"
#include "matpack_arrays.h"
#include "matpack_data.h"
#include "matpack_math.h"

using std::cout;
using std::setprecision;

/* --------------------------------------------------------------
     Helper functions for filling tensors with running numbers,
     starting with start = 1 (default).
   -------------------------------------------------------------- */

Tensor4 fill_tensor4(Index b, Index p, Index r, Index c, Index start = 1) {
  Tensor4 t(b, p, r, c);

  Index fill = start;

  for (Index i = 0; i < b; i++)
    for (Index j = 0; j < p; j++)
      for (Index k = 0; k < r; k++)
        for (Index l = 0; l < c; l++) t(i, j, k, l) = (Numeric)(fill++);

  return t;
}

Tensor5 fill_tensor5(
    Index s, Index b, Index p, Index r, Index c, Index start = 1) {
  Tensor5 t(s, b, p, r, c);

  Index fill = start;

  for (Index i = 0; i < s; i++)
    for (Index j = 0; j < b; j++)
      for (Index k = 0; k < p; k++)
        for (Index l = 0; l < r; l++)
          for (Index m = 0; m < c; m++) t(i, j, k, l, m) = (Numeric)(fill++);

  return t;
}

/*********************
 * Tests for Tensor4 *
 *********************/

void test1() {
  cout << "Test Tensor4:\n\n";

  Tensor4 a = fill_tensor4(2, 3, 3, 4);  // 2 books, 3 pages, 3 rows, 4 columns

  cout << "Dimensions of tensor a:\n"
       << "book = " << a.nbooks() << ", page = " << a.npages()
       << ", row = " << a.nrows() << ", column = " << a.ncols() << '\n';

  cout << "\nmin(a) = " << min(a) << ", max(a) = " << max(a) << "\n\na=\n"
       << a << "\n\n";

  cout << "Second book:\n"
       << a(1, Range(joker), Range(joker), Range(joker)) << "\n\n";

  cout << "First letter of every page of every book:\n"
       << a(Range(joker), Range(joker), 0, 0) << "\n\n";

  cout << "First rows of second page of every book\n"
       << a(Range(joker), 1, 0, Range(joker)) << "\n\n";

  cout << "Third column of third page of first book\n"
       << a(0, 2, Range(joker), 2) << "\n\n";

  cout << "Substract 10 from each element of the second book:\n"
       << (a(1, Range(joker), Range(joker), Range(joker)) -= 10) << "\n\n";

  // Create a sub-tensor of Tensor4
  Tensor4View b = a(Range(joker), Range(0, 2), Range(joker), Range(1, 3));

  cout << "b =\n" << b << '\n';
}

void test2() {
  Tensor4 a = fill_tensor4(2, 2, 4, 5, 4);  // Fill a with running numbers,
                                            // starting with 4

  cout << "a =\n" << a << "\n\n";

  transform(a, sqrt, a);

  cout << "After taking the square-root:\n" << setprecision(3) << a << '\n';
}

void test3() {
  ArrayOfTensor4 a(2);
  Tensor4 b(2, 2, 2, 3, 1.5);

  a[0].resize(2, 1, 4, 5);
  a[0] = 3;
  a[1] = b;

  cout << "a =\n" << a << '\n';
}

/*********************
 * Tests for Tensor5 *
 *********************/

void test4() {
  cout << "Test Tensor5:\n\n";

  Tensor5 a = fill_tensor5(2, 2, 3, 3, 4);  // 2 shelves,
                                            // 2 books, 3 pages,
                                            // 3 rows, 4 columns

  cout << "Dimensions of tensor a:\n"
       << "shelf = " << a.nshelves() << ", book = " << a.nbooks()
       << ", page = " << a.npages() << ", row = " << a.nrows()
       << ", column = " << a.ncols() << '\n';

  cout << "\nmin(a) = " << min(a) << ", max(a) = " << max(a) << "\n\na=\n"
       << a << "\n\n";

  cout << "Second shelf:\n"
       << a(1, Range(joker), Range(joker), Range(joker), Range(joker))
       << "\n\n";

  cout << "First rows of every second page of every book and every shelf:\n"
       << a(Range(joker), Range(joker), Range(0, joker, 2), 0, Range(joker))
       << "\n\n";

  cout << "First and second letter of third column "
       << "of every page and every book of first shelf:\n"
       << a(0, Range(joker), Range(joker), Range(0, 2), 2) << "\n\n";

  cout << "Last two letters of last row of second and third page "
       << "of every book and every shelf:\n"
       << a(Range(joker),
            Range(joker),
            Range(1, 2),
            a.nrows() - 1,
            Range(a.ncols() - 2, 2))
       << "\n\n";

  // Assignment between Tensor5 and Tensor4

  Tensor4 b(3, 2, a.nrows(), a.ncols(), 4);

  // Copy the first two pages of a shelf 1, book 1 to b book 2, page 1-2
  b(1, Range(joker), Range(joker), Range(joker)) =
      a(0, 0, Range(0, 2), Range(joker), Range(joker));

  cout << "b =\n" << b << '\n';
}

void test5() {
  Tensor5 u = fill_tensor5(2, 3, 3, 2, 3);
  Tensor5 v = fill_tensor5(2, 3, 3, 2, 3, 5);  // Fill v with running numbers,
                                               // starting with 5

  cout << "u =\n" << u << "\n\nv =\n" << v << "\n\n";

  u += v;

  cout << "After element-vise addition with tensor v:\n"
       << "u =\n"
       << u << '\n';
}

/*********************
 * Tests for Tensor6 *
 *********************/

void fill_tensor6(
    Tensor6& x, Index v, Index s, Index b, Index p, Index r, Index c) {
  // Lets fill the tensor with special values, so that we can
  // immediately see the vitrine, shelf, etc.
  // 345123 shall mean vitrine 3, shelf 4, book 5, and so on.
  //
  // Will work only if all dimensions are smaller 10.

  x.resize(v, s, b, p, r, c);

  for (Index is = 0; is < s; is++)
    for (Index iv = 0; iv < v; iv++)
      for (Index ib = 0; ib < b; ib++)
        for (Index ip = 0; ip < p; ip++)
          for (Index ir = 0; ir < r; ir++)
            for (Index ic = 0; ic < c; ic++)
              x(iv, is, ib, ip, ir, ic) =
                  (Numeric)(ic + ir * 10 + ip * 100 + ib * 1000 + is * 10000 +
                            iv * 100000);
}

void test6() {
  cout << "Test Tensor6:\n\n";

  Tensor6 a;
  fill_tensor6(a, 3, 2, 2, 3, 3, 4);  // 3 vitrines, 2 shelves,
                                      // 2 books, 3 pages,
                                      // 3 rows, 4 columns

  cout << "Dimensions of tensor a:\n"
       << "vitrine = " << a.nvitrines() << "\n"
       << "shelf   = " << a.nshelves() << "\n"
       << "book    = " << a.nbooks() << "\n"
       << "page    = " << a.npages() << "\n"
       << "row     = " << a.nrows() << "\n"
       << "column  = " << a.ncols() << "\n\n";

  cout << "a(1,1,1,1,1,1) = " << a(1, 1, 1, 1, 1, 1) << "\n\n";

  cout << "a(1,1,1,1,Range(joker),1) = " << a(1, 1, 1, 1, Range(joker), 1)
       << "\n\n";
}

/*********************
 * Tests for Tensor7 *
 *********************/

void fill_tensor7(
    Tensor7& x, Index l, Index v, Index s, Index b, Index p, Index r, Index c) {
  // Lets fill the tensor with special values, so that we can
  // immediately see the vitrine, shelf, etc.
  // 345123 shall mean vitrine 3, shelf 4, book 5, and so on.
  //
  // Will work only if all dimensions are smaller 10.

  x.resize(l, v, s, b, p, r, c);

  for (Index il = 0; il < l; il++)
    for (Index is = 0; is < s; is++)
      for (Index iv = 0; iv < v; iv++)
        for (Index ib = 0; ib < b; ib++)
          for (Index ip = 0; ip < p; ip++)
            for (Index ir = 0; ir < r; ir++)
              for (Index ic = 0; ic < c; ic++)
                x(il, iv, is, ib, ip, ir, ic) =
                    (Numeric)(ic + ir * 10 + ip * 100 + ib * 1000 + is * 10000 +
                              iv * 100000 + il * 1000000);
}

void test7() {
  cout << "Test Tensor7:\n\n";

  Tensor7 a;
  fill_tensor7(a, 2, 3, 2, 2, 3, 3, 4);  // 2 libraries,
                                         // 3 vitrines, 2 shelves,
                                         // 2 books, 3 pages,
                                         // 3 rows, 4 columns

  cout << "Dimensions of tensor a:\n"
       << "libraries = " << a.nlibraries() << "\n"
       << "vitrine   = " << a.nvitrines() << "\n"
       << "shelf     = " << a.nshelves() << "\n"
       << "book      = " << a.nbooks() << "\n"
       << "page      = " << a.npages() << "\n"
       << "row       = " << a.nrows() << "\n"
       << "column    = " << a.ncols() << "\n\n";

  cout << "a(1,1,1,1,1,1,1) = " << setprecision(10) << a(1, 1, 1, 1, 1, 1, 1)
       << "\n\n";

  cout << "a(1,1,1,1,1,Range(joker),1) = " << setprecision(10)
       << a(1, 1, 1, 1, 1, Range(joker), 1) << "\n\n";
}

void test8() {
  cout << "Test Tensor7:\n"
       << "The output of this test should be all ones!\n\n";

  Tensor7 a;
  fill_tensor7(a, 2, 3, 2, 2, 3, 3, 4);  // 2 libraries,
                                         // 3 vitrines, 2 shelves,
                                         // 2 books, 3 pages,
                                         // 3 rows, 4 columns

  Index I(1);     // Select element 1.
  Range R(1, 1);  // Select one element starting at index 1.

  cout << setprecision(10) << a(I, I, I, I, I, I, I) << "\n";
  cout << setprecision(10) << a(I, I, I, I, I, I, R) << "\n";
  cout << setprecision(10) << a(I, I, I, I, I, R, I) << "\n";
  cout << setprecision(10) << a(I, I, I, I, I, R, R) << "\n";
  cout << setprecision(10) << a(I, I, I, I, R, I, I) << "\n";
  cout << setprecision(10) << a(I, I, I, I, R, I, R) << "\n";
  cout << setprecision(10) << a(I, I, I, I, R, R, I) << "\n";
  cout << setprecision(10) << a(I, I, I, I, R, R, R) << "\n";
  cout << setprecision(10) << a(I, I, I, R, I, I, I) << "\n";
  cout << setprecision(10) << a(I, I, I, R, I, I, R) << "\n";
  cout << setprecision(10) << a(I, I, I, R, I, R, I) << "\n";
  cout << setprecision(10) << a(I, I, I, R, I, R, R) << "\n";
  cout << setprecision(10) << a(I, I, I, R, R, I, I) << "\n";
  cout << setprecision(10) << a(I, I, I, R, R, I, R) << "\n";
  cout << setprecision(10) << a(I, I, I, R, R, R, I) << "\n";
  cout << setprecision(10) << a(I, I, I, R, R, R, R) << "\n";
  cout << setprecision(10) << a(I, I, R, I, I, I, I) << "\n";
  cout << setprecision(10) << a(I, I, R, I, I, I, R) << "\n";
  cout << setprecision(10) << a(I, I, R, I, I, R, I) << "\n";
  cout << setprecision(10) << a(I, I, R, I, I, R, R) << "\n";
  cout << setprecision(10) << a(I, I, R, I, R, I, I) << "\n";
  cout << setprecision(10) << a(I, I, R, I, R, I, R) << "\n";
  cout << setprecision(10) << a(I, I, R, I, R, R, I) << "\n";
  cout << setprecision(10) << a(I, I, R, I, R, R, R) << "\n";
  cout << setprecision(10) << a(I, I, R, R, I, I, I) << "\n";
  cout << setprecision(10) << a(I, I, R, R, I, I, R) << "\n";
  cout << setprecision(10) << a(I, I, R, R, I, R, I) << "\n";
  cout << setprecision(10) << a(I, I, R, R, I, R, R) << "\n";
  cout << setprecision(10) << a(I, I, R, R, R, I, I) << "\n";
  cout << setprecision(10) << a(I, I, R, R, R, I, R) << "\n";
  cout << setprecision(10) << a(I, I, R, R, R, R, I) << "\n";
  cout << setprecision(10) << a(I, I, R, R, R, R, R) << "\n";
  cout << setprecision(10) << a(I, R, I, I, I, I, I) << "\n";
  cout << setprecision(10) << a(I, R, I, I, I, I, R) << "\n";
  cout << setprecision(10) << a(I, R, I, I, I, R, I) << "\n";
  cout << setprecision(10) << a(I, R, I, I, I, R, R) << "\n";
  cout << setprecision(10) << a(I, R, I, I, R, I, I) << "\n";
  cout << setprecision(10) << a(I, R, I, I, R, I, R) << "\n";
  cout << setprecision(10) << a(I, R, I, I, R, R, I) << "\n";
  cout << setprecision(10) << a(I, R, I, I, R, R, R) << "\n";
  cout << setprecision(10) << a(I, R, I, R, I, I, I) << "\n";
  cout << setprecision(10) << a(I, R, I, R, I, I, R) << "\n";
  cout << setprecision(10) << a(I, R, I, R, I, R, I) << "\n";
  cout << setprecision(10) << a(I, R, I, R, I, R, R) << "\n";
  cout << setprecision(10) << a(I, R, I, R, R, I, I) << "\n";
  cout << setprecision(10) << a(I, R, I, R, R, I, R) << "\n";
  cout << setprecision(10) << a(I, R, I, R, R, R, I) << "\n";
  cout << setprecision(10) << a(I, R, I, R, R, R, R) << "\n";
  cout << setprecision(10) << a(I, R, R, I, I, I, I) << "\n";
  cout << setprecision(10) << a(I, R, R, I, I, I, R) << "\n";
  cout << setprecision(10) << a(I, R, R, I, I, R, I) << "\n";
  cout << setprecision(10) << a(I, R, R, I, I, R, R) << "\n";
  cout << setprecision(10) << a(I, R, R, I, R, I, I) << "\n";
  cout << setprecision(10) << a(I, R, R, I, R, I, R) << "\n";
  cout << setprecision(10) << a(I, R, R, I, R, R, I) << "\n";
  cout << setprecision(10) << a(I, R, R, I, R, R, R) << "\n";
  cout << setprecision(10) << a(I, R, R, R, I, I, I) << "\n";
  cout << setprecision(10) << a(I, R, R, R, I, I, R) << "\n";
  cout << setprecision(10) << a(I, R, R, R, I, R, I) << "\n";
  cout << setprecision(10) << a(I, R, R, R, I, R, R) << "\n";
  cout << setprecision(10) << a(I, R, R, R, R, I, I) << "\n";
  cout << setprecision(10) << a(I, R, R, R, R, I, R) << "\n";
  cout << setprecision(10) << a(I, R, R, R, R, R, I) << "\n";
  cout << setprecision(10) << a(I, R, R, R, R, R, R) << "\n";
  cout << setprecision(10) << a(R, I, I, I, I, I, I) << "\n";
  cout << setprecision(10) << a(R, I, I, I, I, I, R) << "\n";
  cout << setprecision(10) << a(R, I, I, I, I, R, I) << "\n";
  cout << setprecision(10) << a(R, I, I, I, I, R, R) << "\n";
  cout << setprecision(10) << a(R, I, I, I, R, I, I) << "\n";
  cout << setprecision(10) << a(R, I, I, I, R, I, R) << "\n";
  cout << setprecision(10) << a(R, I, I, I, R, R, I) << "\n";
  cout << setprecision(10) << a(R, I, I, I, R, R, R) << "\n";
  cout << setprecision(10) << a(R, I, I, R, I, I, I) << "\n";
  cout << setprecision(10) << a(R, I, I, R, I, I, R) << "\n";
  cout << setprecision(10) << a(R, I, I, R, I, R, I) << "\n";
  cout << setprecision(10) << a(R, I, I, R, I, R, R) << "\n";
  cout << setprecision(10) << a(R, I, I, R, R, I, I) << "\n";
  cout << setprecision(10) << a(R, I, I, R, R, I, R) << "\n";
  cout << setprecision(10) << a(R, I, I, R, R, R, I) << "\n";
  cout << setprecision(10) << a(R, I, I, R, R, R, R) << "\n";
  cout << setprecision(10) << a(R, I, R, I, I, I, I) << "\n";
  cout << setprecision(10) << a(R, I, R, I, I, I, R) << "\n";
  cout << setprecision(10) << a(R, I, R, I, I, R, I) << "\n";
  cout << setprecision(10) << a(R, I, R, I, I, R, R) << "\n";
  cout << setprecision(10) << a(R, I, R, I, R, I, I) << "\n";
  cout << setprecision(10) << a(R, I, R, I, R, I, R) << "\n";
  cout << setprecision(10) << a(R, I, R, I, R, R, I) << "\n";
  cout << setprecision(10) << a(R, I, R, I, R, R, R) << "\n";
  cout << setprecision(10) << a(R, I, R, R, I, I, I) << "\n";
  cout << setprecision(10) << a(R, I, R, R, I, I, R) << "\n";
  cout << setprecision(10) << a(R, I, R, R, I, R, I) << "\n";
  cout << setprecision(10) << a(R, I, R, R, I, R, R) << "\n";
  cout << setprecision(10) << a(R, I, R, R, R, I, I) << "\n";
  cout << setprecision(10) << a(R, I, R, R, R, I, R) << "\n";
  cout << setprecision(10) << a(R, I, R, R, R, R, I) << "\n";
  cout << setprecision(10) << a(R, I, R, R, R, R, R) << "\n";
  cout << setprecision(10) << a(R, R, I, I, I, I, I) << "\n";
  cout << setprecision(10) << a(R, R, I, I, I, I, R) << "\n";
  cout << setprecision(10) << a(R, R, I, I, I, R, I) << "\n";
  cout << setprecision(10) << a(R, R, I, I, I, R, R) << "\n";
  cout << setprecision(10) << a(R, R, I, I, R, I, I) << "\n";
  cout << setprecision(10) << a(R, R, I, I, R, I, R) << "\n";
  cout << setprecision(10) << a(R, R, I, I, R, R, I) << "\n";
  cout << setprecision(10) << a(R, R, I, I, R, R, R) << "\n";
  cout << setprecision(10) << a(R, R, I, R, I, I, I) << "\n";
  cout << setprecision(10) << a(R, R, I, R, I, I, R) << "\n";
  cout << setprecision(10) << a(R, R, I, R, I, R, I) << "\n";
  cout << setprecision(10) << a(R, R, I, R, I, R, R) << "\n";
  cout << setprecision(10) << a(R, R, I, R, R, I, I) << "\n";
  cout << setprecision(10) << a(R, R, I, R, R, I, R) << "\n";
  cout << setprecision(10) << a(R, R, I, R, R, R, I) << "\n";
  cout << setprecision(10) << a(R, R, I, R, R, R, R) << "\n";
  cout << setprecision(10) << a(R, R, R, I, I, I, I) << "\n";
  cout << setprecision(10) << a(R, R, R, I, I, I, R) << "\n";
  cout << setprecision(10) << a(R, R, R, I, I, R, I) << "\n";
  cout << setprecision(10) << a(R, R, R, I, I, R, R) << "\n";
  cout << setprecision(10) << a(R, R, R, I, R, I, I) << "\n";
  cout << setprecision(10) << a(R, R, R, I, R, I, R) << "\n";
  cout << setprecision(10) << a(R, R, R, I, R, R, I) << "\n";
  cout << setprecision(10) << a(R, R, R, I, R, R, R) << "\n";
  cout << setprecision(10) << a(R, R, R, R, I, I, I) << "\n";
  cout << setprecision(10) << a(R, R, R, R, I, I, R) << "\n";
  cout << setprecision(10) << a(R, R, R, R, I, R, I) << "\n";
  cout << setprecision(10) << a(R, R, R, R, I, R, R) << "\n";
  cout << setprecision(10) << a(R, R, R, R, R, I, I) << "\n";
  cout << setprecision(10) << a(R, R, R, R, R, I, R) << "\n";
  cout << setprecision(10) << a(R, R, R, R, R, R, I) << "\n";
  cout << setprecision(10) << a(R, R, R, R, R, R, R) << "\n";
}

void test9() {
  cout << "Test Tensor7:\n"
       << "The output of this test should be 128\n\n";

  Tensor7 a(2, 3, 2, 2, 3, 3, 4, 0.0);  // 2 libraries,
                                        // 3 vitrines, 2 shelves,
                                        // 2 books, 3 pages,
                                        // 3 rows, 4 columns
                                        // Fill with zeroes.

  Index I(1);     // Select element 1.
  Range R(1, 1);  // Select one element starting at index 1.

  a(I, I, I, I, I, I, I) += 1;
  a(I, I, I, I, I, I, R) += 1;
  a(I, I, I, I, I, R, I) += 1;
  a(I, I, I, I, I, R, R) += 1;
  a(I, I, I, I, R, I, I) += 1;
  a(I, I, I, I, R, I, R) += 1;
  a(I, I, I, I, R, R, I) += 1;
  a(I, I, I, I, R, R, R) += 1;
  a(I, I, I, R, I, I, I) += 1;
  a(I, I, I, R, I, I, R) += 1;
  a(I, I, I, R, I, R, I) += 1;
  a(I, I, I, R, I, R, R) += 1;
  a(I, I, I, R, R, I, I) += 1;
  a(I, I, I, R, R, I, R) += 1;
  a(I, I, I, R, R, R, I) += 1;
  a(I, I, I, R, R, R, R) += 1;
  a(I, I, R, I, I, I, I) += 1;
  a(I, I, R, I, I, I, R) += 1;
  a(I, I, R, I, I, R, I) += 1;
  a(I, I, R, I, I, R, R) += 1;
  a(I, I, R, I, R, I, I) += 1;
  a(I, I, R, I, R, I, R) += 1;
  a(I, I, R, I, R, R, I) += 1;
  a(I, I, R, I, R, R, R) += 1;
  a(I, I, R, R, I, I, I) += 1;
  a(I, I, R, R, I, I, R) += 1;
  a(I, I, R, R, I, R, I) += 1;
  a(I, I, R, R, I, R, R) += 1;
  a(I, I, R, R, R, I, I) += 1;
  a(I, I, R, R, R, I, R) += 1;
  a(I, I, R, R, R, R, I) += 1;
  a(I, I, R, R, R, R, R) += 1;
  a(I, R, I, I, I, I, I) += 1;
  a(I, R, I, I, I, I, R) += 1;
  a(I, R, I, I, I, R, I) += 1;
  a(I, R, I, I, I, R, R) += 1;
  a(I, R, I, I, R, I, I) += 1;
  a(I, R, I, I, R, I, R) += 1;
  a(I, R, I, I, R, R, I) += 1;
  a(I, R, I, I, R, R, R) += 1;
  a(I, R, I, R, I, I, I) += 1;
  a(I, R, I, R, I, I, R) += 1;
  a(I, R, I, R, I, R, I) += 1;
  a(I, R, I, R, I, R, R) += 1;
  a(I, R, I, R, R, I, I) += 1;
  a(I, R, I, R, R, I, R) += 1;
  a(I, R, I, R, R, R, I) += 1;
  a(I, R, I, R, R, R, R) += 1;
  a(I, R, R, I, I, I, I) += 1;
  a(I, R, R, I, I, I, R) += 1;
  a(I, R, R, I, I, R, I) += 1;
  a(I, R, R, I, I, R, R) += 1;
  a(I, R, R, I, R, I, I) += 1;
  a(I, R, R, I, R, I, R) += 1;
  a(I, R, R, I, R, R, I) += 1;
  a(I, R, R, I, R, R, R) += 1;
  a(I, R, R, R, I, I, I) += 1;
  a(I, R, R, R, I, I, R) += 1;
  a(I, R, R, R, I, R, I) += 1;
  a(I, R, R, R, I, R, R) += 1;
  a(I, R, R, R, R, I, I) += 1;
  a(I, R, R, R, R, I, R) += 1;
  a(I, R, R, R, R, R, I) += 1;
  a(I, R, R, R, R, R, R) += 1;
  a(R, I, I, I, I, I, I) += 1;
  a(R, I, I, I, I, I, R) += 1;
  a(R, I, I, I, I, R, I) += 1;
  a(R, I, I, I, I, R, R) += 1;
  a(R, I, I, I, R, I, I) += 1;
  a(R, I, I, I, R, I, R) += 1;
  a(R, I, I, I, R, R, I) += 1;
  a(R, I, I, I, R, R, R) += 1;
  a(R, I, I, R, I, I, I) += 1;
  a(R, I, I, R, I, I, R) += 1;
  a(R, I, I, R, I, R, I) += 1;
  a(R, I, I, R, I, R, R) += 1;
  a(R, I, I, R, R, I, I) += 1;
  a(R, I, I, R, R, I, R) += 1;
  a(R, I, I, R, R, R, I) += 1;
  a(R, I, I, R, R, R, R) += 1;
  a(R, I, R, I, I, I, I) += 1;
  a(R, I, R, I, I, I, R) += 1;
  a(R, I, R, I, I, R, I) += 1;
  a(R, I, R, I, I, R, R) += 1;
  a(R, I, R, I, R, I, I) += 1;
  a(R, I, R, I, R, I, R) += 1;
  a(R, I, R, I, R, R, I) += 1;
  a(R, I, R, I, R, R, R) += 1;
  a(R, I, R, R, I, I, I) += 1;
  a(R, I, R, R, I, I, R) += 1;
  a(R, I, R, R, I, R, I) += 1;
  a(R, I, R, R, I, R, R) += 1;
  a(R, I, R, R, R, I, I) += 1;
  a(R, I, R, R, R, I, R) += 1;
  a(R, I, R, R, R, R, I) += 1;
  a(R, I, R, R, R, R, R) += 1;
  a(R, R, I, I, I, I, I) += 1;
  a(R, R, I, I, I, I, R) += 1;
  a(R, R, I, I, I, R, I) += 1;
  a(R, R, I, I, I, R, R) += 1;
  a(R, R, I, I, R, I, I) += 1;
  a(R, R, I, I, R, I, R) += 1;
  a(R, R, I, I, R, R, I) += 1;
  a(R, R, I, I, R, R, R) += 1;
  a(R, R, I, R, I, I, I) += 1;
  a(R, R, I, R, I, I, R) += 1;
  a(R, R, I, R, I, R, I) += 1;
  a(R, R, I, R, I, R, R) += 1;
  a(R, R, I, R, R, I, I) += 1;
  a(R, R, I, R, R, I, R) += 1;
  a(R, R, I, R, R, R, I) += 1;
  a(R, R, I, R, R, R, R) += 1;
  a(R, R, R, I, I, I, I) += 1;
  a(R, R, R, I, I, I, R) += 1;
  a(R, R, R, I, I, R, I) += 1;
  a(R, R, R, I, I, R, R) += 1;
  a(R, R, R, I, R, I, I) += 1;
  a(R, R, R, I, R, I, R) += 1;
  a(R, R, R, I, R, R, I) += 1;
  a(R, R, R, I, R, R, R) += 1;
  a(R, R, R, R, I, I, I) += 1;
  a(R, R, R, R, I, I, R) += 1;
  a(R, R, R, R, I, R, I) += 1;
  a(R, R, R, R, I, R, R) += 1;
  a(R, R, R, R, R, I, I) += 1;
  a(R, R, R, R, R, I, R) += 1;
  a(R, R, R, R, R, R, I) += 1;
  a(R, R, R, R, R, R, R) += 1;

  cout << "The sum is: " << a(I, I, I, I, I, I, I) << "\n";
}

int main() {
  test9();
  return 0;
}

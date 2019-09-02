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

#include <cmath>
#include <iostream>
#include "array.h"
#include "interpolation.h"
#include "interpolation_poly.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "xml_io.h"

void test01() {
  cout << "Simple interpolation cases\n"
       << "--------------------------\n";
  //  Vector og(5,5,-1);                // 5,4,3,2,1
  Vector og(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng(2, 5, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0

  cout << "og:\n" << og << "\n";
  cout << "ng:\n" << ng << "\n";

  // To store the grid positions:
  ArrayOfGridPos gp(ng.nelem());

  gridpos(gp, og, ng);
  cout << "gp:\n" << gp << "\n";

  cout << "1D:\n"
       << "---\n";
  {
    // To store interpolation weights:
    Matrix itw(gp.nelem(), 2);
    interpweights(itw, gp);

    cout << "itw:\n" << itw << "\n";

    // Original field:
    Vector of(og.nelem(), 0);
    of[2] = 10;  // 0, 0, 10, 0, 0

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
    Matrix itw(gp.nelem(), 4);
    interpweights(itw, gp, gp);

    cout << "itw:\n" << itw << "\n";

    // Original field:
    Matrix of(og.nelem(), og.nelem(), 0);
    of(2, 2) = 10;  // 0 Matrix with 10 in the middle

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
    Matrix itw(gp.nelem(), 64);
    interpweights(itw, gp, gp, gp, gp, gp, gp);

    //    cout << "itw:\n" << itw << "\n";

    // Original field:
    Index n = og.nelem();
    Tensor6 of(n, n, n, n, n, n, 0);
    of(2, 2, 2, 2, 2, 2) = 10;  // 0 Tensor with 10 in the middle

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
    Tensor3 itw(gp.nelem(), gp.nelem(), 4);
    interpweights(itw, gp, gp);

    for (Index i = 0; i < itw.ncols(); ++i)
      cout << "itw " << i << ":\n"
           << itw(Range(joker), Range(joker), i) << "\n";

    // Original field:
    Matrix of(og.nelem(), og.nelem(), 0);
    of(2, 2) = 10;  // 0 Matrix with 10 in the middle

    cout << "of:\n" << of << "\n";

    // Interpolated field:
    Matrix nf(ng.nelem(), ng.nelem());

    interp(nf, itw, of, gp, gp);

    cout << "nf:\n" << nf << "\n";
  }

  cout << "Green 6D:\n"
       << "---------\n";
  {
    // To store interpolation weights:
    Tensor7 itw(gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                64);
    interpweights(itw, gp, gp, gp, gp, gp, gp);

    // Original field:
    Tensor6 of(og.nelem(),
               og.nelem(),
               og.nelem(),
               og.nelem(),
               og.nelem(),
               og.nelem(),
               0);
    of(2, 2, 2, 2, 2, 2) = 10;  // 0 Tensor with 10 in the middle

    cout << "Middle slice of of:\n"
         << of(2, 2, 2, 2, Range(joker), Range(joker)) << "\n";

    // Interpolated field:
    Tensor6 nf(
        ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem());

    interp(nf, itw, of, gp, gp, gp, gp, gp, gp);

    cout << "Last slice of nf:\n"
         << nf(4, 4, 4, 4, Range(joker), Range(joker)) << "\n";
  }
}

void test02(Index n) {
  cout << "Test whether for loop or iterator are quicker\n"
       << "a) for loop\n"
       << "---------------------------------------------\n";

  Vector a(n);
  for (Index i = 0; i < a.nelem(); ++i) a[i] = (Numeric)i;
}

void test03(Index n) {
  cout << "Test whether for loop or iterator are quicker\n"
       << "b) iterator\n"
       << "---------------------------------------------\n";

  Vector a(n);
  Iterator1D ai = a.begin();
  const Iterator1D ae = a.end();
  Index i = 0;
  for (; ai != ae; ++ai, ++i) *ai = (Numeric)i;
}

// Result: Both are almost equally fast, with a slight advantage of
// the for loop if compiler optimization is enabled.

void test04() {
  cout << "Green type interpolation of all pages of a Tensor3\n";

  // The original Tensor is called a, the new one n.

  // 10 pages, 20 rows, 30 columns, all grids are: 1,2,3
  Vector a_pgrid(1, 3, 1), a_rgrid(1, 3, 1), a_cgrid(1, 3, 1);
  Tensor3 a(a_pgrid.nelem(), a_rgrid.nelem(), a_cgrid.nelem());

  a = 0;
  // Put some simple numbers in the middle of each page:
  a(0, 1, 1) = 10;
  a(1, 1, 1) = 20;
  a(2, 1, 1) = 30;

  // New row and column grids:
  // 1, 1.5, 2, 2.5, 3
  Vector n_rgrid(1, 5, .5), n_cgrid(1, 5, .5);
  Tensor3 n(a_pgrid.nelem(), n_rgrid.nelem(), n_cgrid.nelem());

  // So, n has the same number of pages as a, but more rows and columns.

  // Get the grid position arrays:
  ArrayOfGridPos n_rgp(n_rgrid.nelem());  // For row grid positions.
  ArrayOfGridPos n_cgp(n_cgrid.nelem());  // For column grid positions.

  gridpos(n_rgp, a_rgrid, n_rgrid);
  gridpos(n_cgp, a_cgrid, n_cgrid);

  // Get the interpolation weights:
  Tensor3 itw(n_rgrid.nelem(), n_cgrid.nelem(), 4);
  interpweights(itw, n_rgp, n_cgp);

  // Do a "green" interpolation for all pages of the Tensor a:

  for (Index i = 0; i < a.npages(); ++i) {
    // Select the current page of both a and n:
    ConstMatrixView ap = a(i, Range(joker), Range(joker));
    MatrixView np = n(i, Range(joker), Range(joker));

    // Do the interpolation:
    interp(np, itw, ap, n_rgp, n_cgp);

    // Note that this is efficient, because interpolation weights and
    // grid positions are re-used.
  }

  cout << "Original field:\n";
  for (Index i = 0; i < a.npages(); ++i)
    cout << "page " << i << ":\n" << a(i, Range(joker), Range(joker)) << "\n";

  cout << "Interpolated field:\n";
  for (Index i = 0; i < n.npages(); ++i)
    cout << "page " << i << ":\n" << n(i, Range(joker), Range(joker)) << "\n";
}

void test05() {
  cout << "Very simple interpolation case\n";

  Vector og(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng(2, 5, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0

  cout << "Original grid:\n" << og << "\n";
  cout << "New grid:\n" << ng << "\n";

  // To store the grid positions:
  ArrayOfGridPos gp(ng.nelem());

  gridpos(gp, og, ng);
  cout << "Grid positions:\n" << gp;

  // To store interpolation weights:
  Matrix itw(gp.nelem(), 2);
  interpweights(itw, gp);

  cout << "Interpolation weights:\n" << itw << "\n";

  // Original field:
  Vector of(og.nelem(), 0);
  of[2] = 10;  // 0, 0, 10, 0, 0

  cout << "Original field:\n" << of << "\n";

  // Interpolated field:
  Vector nf(ng.nelem());

  interp(nf, itw, of, gp);

  cout << "New field:\n" << nf << "\n";
}

void test06() {
  cout << "Simple extrapolation cases\n"
       << "--------------------------\n";
  //  Vector og(5,5,-1);                // 5,4,3,2,1
  Vector og(1, 5, +1);               // 1, 2, 3, 4, 5
  Vector ng{0.9, 1.5, 3, 4.5, 5.1};  // 0.9, 1.5, 3, 4.5, 5.1

  cout << "og:\n" << og << "\n";
  cout << "ng:\n" << ng << "\n";

  // To store the grid positions:
  ArrayOfGridPos gp(ng.nelem());

  gridpos(gp, og, ng);
  cout << "gp:\n" << gp << "\n";

  cout << "1D:\n"
       << "---\n";
  {
    // To store interpolation weights:
    Matrix itw(gp.nelem(), 2);
    interpweights(itw, gp);

    cout << "itw:\n" << itw << "\n";

    // Original field:
    Vector of(og.nelem(), 0);
    for (Index i = 0; i < og.nelem(); ++i)
      of[i] = (Numeric)(10 * (i + 1));  // 10, 20, 30, 40, 50

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
    Matrix itw(gp.nelem(), 4);
    interpweights(itw, gp, gp);

    cout << "itw:\n" << itw << "\n";

    // Original field:
    Matrix of(og.nelem(), og.nelem(), 0);
    of(2, 2) = 10;  // 0 Matrix with 10 in the middle

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
    Matrix itw(gp.nelem(), 64);
    interpweights(itw, gp, gp, gp, gp, gp, gp);

    //    cout << "itw:\n" << itw << "\n";

    // Original field:
    Index n = og.nelem();
    Tensor6 of(n, n, n, n, n, n, 0);
    of(2, 2, 2, 2, 2, 2) = 10;  // 0 Tensor with 10 in the middle

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
    Tensor3 itw(gp.nelem(), gp.nelem(), 4);
    interpweights(itw, gp, gp);

    for (Index i = 0; i < itw.ncols(); ++i)
      cout << "itw " << i << ":\n"
           << itw(Range(joker), Range(joker), i) << "\n";

    // Original field:
    Matrix of(og.nelem(), og.nelem(), 0);
    of(2, 2) = 10;  // 0 Matrix with 10 in the middle

    cout << "of:\n" << of << "\n";

    // Interpolated field:
    Matrix nf(ng.nelem(), ng.nelem());

    interp(nf, itw, of, gp, gp);

    cout << "nf:\n" << nf << "\n";
  }

  cout << "Green 6D:\n"
       << "---------\n";
  {
    // To store interpolation weights:
    Tensor7 itw(gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                gp.nelem(),
                64);
    interpweights(itw, gp, gp, gp, gp, gp, gp);

    // Original field:
    Tensor6 of(og.nelem(),
               og.nelem(),
               og.nelem(),
               og.nelem(),
               og.nelem(),
               og.nelem(),
               0);
    of(2, 2, 2, 2, 2, 2) = 10;  // 0 Tensor with 10 in the middle

    cout << "Middle slice of of:\n"
         << of(2, 2, 2, 2, Range(joker), Range(joker)) << "\n";

    // Interpolated field:
    Tensor6 nf(
        ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem());

    interp(nf, itw, of, gp, gp, gp, gp, gp, gp);

    cout << "Last slice of nf:\n"
         << nf(4, 4, 4, 4, Range(joker), Range(joker)) << "\n";
  }
}

// Test polynomial interpolation (included by CE)
void test07() {
  // FileType ftype = FILE_TYPE_ASCII;
  Vector new_x(0, 21, +0.25);
  Vector x(0, 10, +1);

  ArrayOfGridPos gp(new_x.nelem());
  gridpos(gp, x, new_x);

  Vector y1(x.nelem());
  Vector y2(x.nelem());
  Vector y3(x.nelem());

  for (Index i = 0; i < x.nelem(); i++) {
    // linear function
    y1[i] = 3 * x[i];
    // cubic function
    y2[i] = pow(x[i], 3) + 2;
    // trigonometric function
    y3[i] = sin(x[i]);
  }

  // Linear interpolation:
  Matrix itw(gp.nelem(), 2);
  interpweights(itw, gp);

  Vector y1_lin(new_x.nelem());
  Vector y2_lin(new_x.nelem());
  Vector y3_lin(new_x.nelem());

  interp(y1_lin, itw, y1, gp);
  interp(y2_lin, itw, y2, gp);
  interp(y3_lin, itw, y3, gp);

  cout << "y1_lin = [" << y1_lin << "];\n";
  cout << "y2_lin = [" << y2_lin << "];\n";
  cout << "y3_lin = [" << y3_lin << "];\n";

  // Cubic interpolation:
  Vector y1_cub(new_x.nelem());
  Vector y2_cub(new_x.nelem());
  Vector y3_cub(new_x.nelem());

  for (Index i = 0; i < new_x.nelem(); i++) {
    y1_cub[i] = interp_poly(x, y1, new_x[i], gp[i]);
    y2_cub[i] = interp_poly(x, y2, new_x[i], gp[i]);
    y3_cub[i] = interp_poly(x, y3, new_x[i], gp[i]);
  }

  cout << "y1_cub = [" << y1_cub << "];\n";
  cout << "y2_cub = [" << y2_cub << "];\n";
  cout << "y3_cub = [" << y3_cub << "];\n";

  // Stefan's new polynomial interpolation routines:
  Index order = 2;

  Vector y1_new(new_x.nelem());
  Vector y2_new(new_x.nelem());
  Vector y3_new(new_x.nelem());

  ArrayOfGridPosPoly gpp(new_x.nelem());
  gridpos_poly(gpp, x, new_x, order);
  Matrix itwp(new_x.nelem(), order + 1);
  interpweights(itwp, gpp);

  interp(y1_new, itwp, y1, gpp);
  interp(y2_new, itwp, y2, gpp);
  interp(y3_new, itwp, y3, gpp);

  cout << "y1_new = [" << y1_new << "];\n";
  cout << "y2_new = [" << y2_new << "];\n";
  cout << "y3_new = [" << y3_new << "];\n";
}

void test08() {
  cout << "Very simple interpolation case for the "
       << "new higher order polynomials.\n";

  Vector og(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng(2, 9, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0 ... 4.0

  cout << "Original grid:\n" << og << "\n";
  cout << "New grid:\n" << ng << "\n";

  // To store the grid positions:
  ArrayOfGridPosPoly gp(ng.nelem());

  Index order = 0;  // Interpolation order.

  gridpos_poly(gp, og, ng, order);
  cout << "Grid positions:\n" << gp;

  // To store interpolation weights:
  Matrix itw(gp.nelem(), order + 1);
  interpweights(itw, gp);

  cout << "Interpolation weights:\n" << itw << "\n";

  // Original field:
  Vector of(og.nelem(), 0);
  of[2] = 10;  // 0, 0, 10, 0, 0

  cout << "Original field:\n" << of << "\n";

  // Interpolated field:
  Vector nf(ng.nelem());

  interp(nf, itw, of, gp);

  cout << "New field (order=" << order << "):\n" << nf << "\n";

  cout << "All orders systematically:\n";
  for (order = 0; order < 5; ++order) {
    gridpos_poly(gp, og, ng, order);
    itw.resize(gp.nelem(), order + 1);
    interpweights(itw, gp);
    interp(nf, itw, of, gp);

    cout << "order " << order << ": ";
    for (Index i = 0; i < nf.nelem(); ++i) cout << setw(8) << nf[i] << " ";
    cout << "\n";
  }
}

int main() { test08(); }

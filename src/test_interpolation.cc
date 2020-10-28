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
#include "interpolation_lagrange.h"
#include "interpolation_poly.h"
#include "math_funcs.h"
#include "matpackVII.h"
#include "xml_io.h"
#include "auto_md.h"

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

void test09() {
  Vector og(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng(2, 9, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0 ... 4.0
  Vector yi{5, -2, 50, 2, 1};
  std::cout << yi << '\n';
  
  for (auto x: ng) {
    GridPos gp;
    gridpos(gp, og, x);
    const Interpolation::Lagrange lag(x, og);
    std::cout << "gp " << gp << "lag: " << lag << '\n';
    
    Vector iwgp(2);
    interpweights(iwgp, gp);
    const Vector iwlag = interpweights(lag);
    std::cout << "gp " << iwgp << "\nlag: " << iwlag << '\n';
    
    
    std::cout << "gp " << interp(iwgp, yi, gp) << "\nlag: " << interp(yi, iwlag, lag) << '\n' << '\n';
  }
}

void test10() {
  Vector xi(1, 100, +1);    // 1, 2, 3, 4, 5 ... 100
  Vector xn(2, 900, +0.1);
  Vector yi(100);
  for (Index i=0; i<100; i++) yi[i] = std::exp(-xi[i]/3.14);
  
  Index order = 10;
  {
    for (auto x: xn) {
      GridPosPoly gp;
      gridpos_poly(gp, xi, x, order);
      const Interpolation::Lagrange lag(x, xi, order);
      
      Vector iwgp(order+1);
      interpweights(iwgp, gp);
      const Vector iwlag = interpweights(lag);
      
      std::cout << order << ' ' << x << ' ' << interp(iwgp, yi, gp) << " " << interp(yi, iwlag, lag) << '\n';
    }
  }
}

void test11() {
  constexpr int N=20;
  Vector x0(1, N, +1);    // 1, 2, 3, 4, 5 ... 10
  Vector x1(1, N, +2);    // 1, 2, 3, 4, 5 ... 10
  Vector x0n(100);
  Vector x1n(100);
  nlinspace(x0n, x0[0], x0[N-1], 100);
  nlinspace(x1n, x1[0], x1[N-1], 100);
  Matrix yi(N, N);
  for (Index i=0; i<N; i++) for (Index j=0; j<N; j++) yi(i, j) = std::exp(-x0[i]/3.14 + x1[j]/3.14);
  
  std::cerr << x0 << '\n';
  std::cerr << x1 << '\n';
  std::cerr << yi << '\n';
  
  constexpr Index order = 5;
  {
      const auto lag0 = Interpolation::LagrangeVector(x0n, x0, order, 0.5);
      const auto lag1 = Interpolation::LagrangeVector(x1n, x1, order, 0.5);
      const auto iwlag = interpweights(lag0, lag1);
      std::cout << x0n << '\n' << x1n << '\n' << reinterp(yi, iwlag, lag0, lag1) << '\n';
  }
}

constexpr bool is_around(Numeric x, Numeric x0, Numeric e=1e-12) {
  return x - x0 < e and x - x0 > -e;
}

void test12() {
  constexpr int N=10;
  constexpr std::array<Numeric, N> y{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  constexpr std::array<Numeric, N> y2{1, 2*2, 3*3, 4*4, 5*5, 6*6, 7*7, 8*8, 9*9, 10*10};
  constexpr std::array<Numeric, N> y3{1, 2*2*2, 3*3*3, 4*4*4, 5*5*5, 6*6*6, 7*7*7, 8*8*8, 9*9*9, 10*10*10};
  constexpr Numeric x = 1.5;
  
  // Set up the interpolation Lagranges
  constexpr Interpolation::FixedLagrange<1> lin(x, 1, 2);
  constexpr Interpolation::FixedLagrange<2> sqr(x, 1, 2, 3);
  constexpr Interpolation::FixedLagrange<3> cub(x, 1, 2, 3, 4);
  
  // Set up the interpolation weights
  constexpr auto lin_iw = interpweights(lin);
  constexpr auto sqr_iw = interpweights(sqr);
  constexpr auto cub_iw = interpweights(cub);
  
  // Get interpolation value
  constexpr auto lin_lin = interp(y, lin_iw, lin);
  constexpr auto lin_sqr = interp(y2, lin_iw, lin);
  constexpr auto lin_cub = interp(y3, lin_iw, lin);
  constexpr auto sqr_lin = interp(y, sqr_iw, sqr);
  constexpr auto sqr_sqr = interp(y2, sqr_iw, sqr);
  constexpr auto sqr_cub = interp(y3, sqr_iw, sqr);
  constexpr auto cub_lin = interp(y, cub_iw, cub);
  constexpr auto cub_sqr = interp(y2, cub_iw, cub);
  constexpr auto cub_cub = interp(y3, cub_iw, cub);
  
  // Compile-time check that these values are good
  static_assert(is_around(lin_lin, x));
  static_assert(is_around(sqr_sqr, x*x));
  static_assert(is_around(cub_cub, x*x*x));
  
  // Should output 1.5 2.5 4.5 1.5 2.25 3 1.5 2.25 3.375
  std::cout << lin_lin << ' ' << lin_sqr << ' ' << lin_cub << ' '
            << sqr_lin << ' ' << sqr_sqr << ' ' << sqr_cub << ' '
            << cub_lin << ' ' << cub_sqr << ' ' << cub_cub << '\n';
}

void test13() {
  constexpr int N=10;
  constexpr std::array<Numeric, N> y{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  constexpr std::array<Numeric, N> y2{1, 2*2, 3*3, 4*4, 5*5, 6*6, 7*7, 8*8, 9*9, 10*10};
  constexpr std::array<Numeric, N> y3{1, 2*2*2, 3*3*3, 4*4*4, 5*5*5, 6*6*6, 7*7*7, 8*8*8, 9*9*9, 10*10*10};
  constexpr Numeric x = 1.5;
  constexpr Numeric dx = 1e-2;
  
  // Set up the interpolation Lagranges
  constexpr Interpolation::FixedLagrange<1> lin(x, 1, 2);
  constexpr Interpolation::FixedLagrange<2> sqr(x, 1, 2, 3);
  constexpr Interpolation::FixedLagrange<3> cub(x, 1, 2, 3, 4);
  constexpr Interpolation::FixedLagrange<1> dlin(x+dx, 1, 2);
  constexpr Interpolation::FixedLagrange<2> dsqr(x+dx, 1, 2, 3);
  constexpr Interpolation::FixedLagrange<3> dcub(x+dx, 1, 2, 3, 4);
  
  // Set up the interpolation weights
  constexpr auto lin_iw = interpweights(lin);
  constexpr auto sqr_iw = interpweights(sqr);
  constexpr auto cub_iw = interpweights(cub);
  constexpr auto dlin_iw = dinterpweights(lin);
  constexpr auto dsqr_iw = dinterpweights(sqr);
  constexpr auto dcub_iw = dinterpweights(cub);
  constexpr auto alt_lin_iw = interpweights(dlin);
  constexpr auto alt_sqr_iw = interpweights(dsqr);
  constexpr auto alt_cub_iw = interpweights(dcub);
  
  // Get interpolation value
  constexpr auto lin_lin = interp(y, lin_iw, lin);
  constexpr auto sqr_sqr = interp(y2, sqr_iw, sqr);
  constexpr auto cub_cub = interp(y3, cub_iw, cub);
  constexpr auto dlin_lin = interp(y, dlin_iw, dlin);
  constexpr auto dsqr_sqr = interp(y2, dsqr_iw, dsqr);
  constexpr auto dcub_cub = interp(y3, dcub_iw, dcub);
  constexpr auto alt_lin_lin = interp(y, alt_lin_iw, dlin);
  constexpr auto alt_sqr_sqr = interp(y2, alt_sqr_iw, dsqr);
  constexpr auto alt_cub_cub = interp(y3, alt_cub_iw, dcub);
  
  // Compile-time check that these values are good
  static_assert(is_around(lin_lin, x));
  static_assert(is_around(sqr_sqr, x*x));
  static_assert(is_around(cub_cub, x*x*x));
  static_assert(is_around(dlin_lin, 1));
  static_assert(is_around(dsqr_sqr, 2*x));
  static_assert(is_around(dcub_cub, 3*x*x));
  
  // Should be
  // 1.5  2.25   3.375
  // 1    3      6.75
  // 1.51 2.2801 3.44295
  // ~1   ~3     ~6.75 (ofc, worse for less linear cases)
  std::cout << lin_lin << ' ' << ' ' << sqr_sqr << ' ' << ' ' << ' ' << cub_cub << '\n' 
            << dlin_lin << ' ' << ' ' << ' ' << ' ' << dsqr_sqr << ' ' << ' ' << ' ' << ' ' << ' ' << ' '  << dcub_cub << '\n'
            << alt_lin_lin << ' ' << alt_sqr_sqr << ' ' << alt_cub_cub << '\n'
            << (alt_lin_lin - lin_lin)/dx << ' ' << ' ' << ' ' << ' ' << (alt_sqr_sqr - sqr_sqr)/dx << ' ' << ' ' << ' ' << (alt_cub_cub - cub_cub)/dx << '\n';
}

void test14() {
  Vector y{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  Vector y2{1, 2*2, 3*3, 4*4, 5*5, 6*6, 7*7, 8*8, 9*9, 10*10};
  Vector y3{1, 2*2*2, 3*3*3, 4*4*4, 5*5*5, 6*6*6, 7*7*7, 8*8*8, 9*9*9, 10*10*10};
  
  Vector x{0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};
  
  // Set up the interpolation Lagranges
  auto lin = Interpolation::FixedLagrangeVector<1>(x, y, 1.0, Interpolation::LagrangeType::Linear);
  auto sqr = Interpolation::FixedLagrangeVector<2>(x, y, 1.0, Interpolation::LagrangeType::Linear);
  auto cub = Interpolation::FixedLagrangeVector<3>(x, y, 1.0, Interpolation::LagrangeType::Linear);
  
  // Set up the interpolation weights
  auto lin_iw = interpweights(lin);
  auto sqr_iw = interpweights(sqr);
  auto cub_iw = interpweights(cub);
  auto dlin_iw = dinterpweights(lin);
  auto dsqr_iw = dinterpweights(sqr);
  auto dcub_iw = dinterpweights(cub);
  
  // Get interpolation value
  auto lin_lin = reinterp(y, lin_iw, lin);
  auto sqr_sqr = reinterp(y2, sqr_iw, sqr);
  auto cub_cub = reinterp(y3, cub_iw, cub);
  auto dlin_lin = reinterp(y, dlin_iw, lin);
  auto dsqr_sqr = reinterp(y2, dsqr_iw, sqr);
  auto dcub_cub = reinterp(y3, dcub_iw, cub);
  
  // Print the values and derivatives
  std::cout << lin_lin << '\n' << sqr_sqr << '\n' << cub_cub << '\n';
  std::cout << dlin_lin << '\n' << dsqr_sqr << '\n' << dcub_cub << '\n';
}

void test15() {
  Verbosity verbosity;
  const Index old_size = 50000;
  const Index new_size = 100000;
  Vector xold;
  VectorNLinSpace(xold, old_size, 1, 100, verbosity);
  Vector xnew;
  VectorNLinSpace(xnew, new_size, 1, 100, verbosity);
  Vector yold;
  VectorNLinSpace(yold, old_size, 100, 200, verbosity);

  const Index f_order = 3;
  chk_interpolation_grids(
      "Frequency interpolation for cross sections", xold, xnew, f_order);

  constexpr Index N = 500;
  std::vector<TimeStep> new_time(N);
  std::vector<TimeStep> old_time(N);
  
  std::cout << "========== New interp ==========" << std::endl;
  for (Index i=0; i<N; i++) {
    Time now_new;
    const auto lag_interp =
    Interpolation::LagrangeVector(xnew, xold, f_order, 0.5);
    const auto iw_interp = interpweights(lag_interp);
    const auto ynew = reinterp(yold, iw_interp, lag_interp);
    new_time[i] = Time() - now_new;
  }
  std::sort(new_time.begin(), new_time.end());
  std::cerr << "Median time: " << new_time[N/2] << '\n';
  
  std::cout << "========== Old interp ==========" << std::endl;
  ArrayOfGridPosPoly f_gp(xnew.nelem()), T_gp(1);
  Matrix itw(f_gp.nelem(), f_order + 1);
  Vector ynew_oldinterp(new_size);
  for (Index i=0; i<N; i++) {
    Time now_old;
    gridpos_poly(f_gp, xold, xnew, f_order);
    interpweights(itw, f_gp);
    interp(ynew_oldinterp, itw, yold, f_gp);
    old_time[i] = Time() - now_old;
  }
  std::sort(old_time.begin(), old_time.end());
  std::cerr << "Median time: " << old_time[N/2] << '\n';
}

int main() { test15(); }

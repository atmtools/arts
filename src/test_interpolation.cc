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
#include "arts_conversions.h"
#include "auto_md.h"
#include "gridded_fields.h"
#include "interpolation.h"
#include "interp.h"
#include "math_funcs.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_math.h"
#include "nonstd.h"
#include "xml_io.h"

void test01() {
  cout << "Simple interpolation cases\n"
       << "--------------------------\n";
  //  Vector og(5,5,-1);                // 5,4,3,2,1
  Vector og=uniform_grid(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng=uniform_grid(2, 5, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0

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
    Tensor7 itw(gp.nelem(), gp.nelem(), gp.nelem(), gp.nelem(), gp.nelem(),
                gp.nelem(), 64);
    interpweights(itw, gp, gp, gp, gp, gp, gp);

    // Original field:
    Tensor6 of(og.nelem(), og.nelem(), og.nelem(), og.nelem(), og.nelem(),
               og.nelem(), 0);
    of(2, 2, 2, 2, 2, 2) = 10;  // 0 Tensor with 10 in the middle

    cout << "Middle slice of of:\n"
         << of(2, 2, 2, 2, Range(joker), Range(joker)) << "\n";

    // Interpolated field:
    Tensor6 nf(ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(),
               ng.nelem());

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
  auto ai = a.begin();
  const auto ae = a.end();
  Index i = 0;
  for (; ai != ae; ++ai, ++i) *ai = (Numeric)i;
}

// Result: Both are almost equally fast, with a slight advantage of
// the for loop if compiler optimization is enabled.

void test04() {
  cout << "Green type interpolation of all pages of a Tensor3\n";

  // The original Tensor is called a, the new one n.

  // 10 pages, 20 rows, 30 columns, all grids are: 1,2,3
  Vector a_pgrid=uniform_grid(1, 3, 1), a_rgrid=uniform_grid(1, 3, 1), a_cgrid=uniform_grid(1, 3, 1);
  Tensor3 a(a_pgrid.nelem(), a_rgrid.nelem(), a_cgrid.nelem());

  a = 0;
  // Put some simple numbers in the middle of each page:
  a(0, 1, 1) = 10;
  a(1, 1, 1) = 20;
  a(2, 1, 1) = 30;

  // New row and column grids:
  // 1, 1.5, 2, 2.5, 3
  Vector n_rgrid=uniform_grid(1, 5, .5), n_cgrid=uniform_grid(1, 5, .5);
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

  Vector og=uniform_grid(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng=uniform_grid(2, 5, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0

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
  Vector og=uniform_grid(1, 5, +1);               // 1, 2, 3, 4, 5
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
    Tensor7 itw(gp.nelem(), gp.nelem(), gp.nelem(), gp.nelem(), gp.nelem(),
                gp.nelem(), 64);
    interpweights(itw, gp, gp, gp, gp, gp, gp);

    // Original field:
    Tensor6 of(og.nelem(), og.nelem(), og.nelem(), og.nelem(), og.nelem(),
               og.nelem(), 0);
    of(2, 2, 2, 2, 2, 2) = 10;  // 0 Tensor with 10 in the middle

    cout << "Middle slice of of:\n"
         << of(2, 2, 2, 2, Range(joker), Range(joker)) << "\n";

    // Interpolated field:
    Tensor6 nf(ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(), ng.nelem(),
               ng.nelem());

    interp(nf, itw, of, gp, gp, gp, gp, gp, gp);

    cout << "Last slice of nf:\n"
         << nf(4, 4, 4, 4, Range(joker), Range(joker)) << "\n";
  }
}

void test09() {
  // Original and new grid
  Vector og=uniform_grid(1, 5, +1);    // 1, 2, 3, 4, 5
  Vector ng=uniform_grid(2, 9, 0.25);  // 2.0, 2,25, 2.5, 2.75, 3.0 ... 4.0
  
  // Original data
  Vector yi{5, -2, 50, 2, 1};
  
  std::cout << og << '\n' << yi << '\n' << '\n';

  // Do the fit point-by-point:
  for (auto x : ng) {
    // Simple linear interpolation:
    GridPos gp;
    gridpos(gp, og, x);
    
    // General Lagrange interpolation, special case interpolation order 1:
    const LagrangeInterpolation lag(0, x, og, 1);
    
    std::cout << "gp point:  " << gp << "lag point:\n" << lag;

    // Linear interpolation weights:
    Vector iwgp(2);
    interpweights(iwgp, gp);
    
    // General Lagrange interpolation weights:
    const Vector iwlag = interpweights(lag);
    
    std::cout << "gp iw:  " << iwgp << "\nlag iw: " << iwlag << '\n';

    std::cout << "gp res:  " << interp(iwgp, yi, gp)
              << "\nlag res: " << interp(yi, iwlag, lag) << '\n'
              << '\n';
  }
  
  // Linear interpolation of all points
  ArrayOfGridPos gp(ng.nelem());
  gridpos(gp, og, ng);
  Matrix gp_iw(gp.nelem(), 2);
  interpweights(gp_iw, gp);
  Vector gp_y(ng.nelem());
  interp(gp_y, gp_iw, yi, gp);
  std::cout << "gp:  " << gp_y << '\n';
  
  // General Lagrange interpolation, special case interpolation order 1:
  const auto lag = my_interp::lagrange_interpolation_list<LagrangeInterpolation>(ng, og, 1);
  std::cout << "lag: " << reinterp(yi, interpweights(lag), lag) << '\n';
}

void test11() {
  constexpr int N = 20;
  Vector x0=uniform_grid(1, N, +1);  // 1, 2, 3, 4, 5 ... 10
  Vector x1=uniform_grid(1, N, +2);  // 1, 2, 3, 4, 5 ... 10
  Vector x0n(100);
  Vector x1n(100);
  nlinspace(x0n, x0[0], x0[N - 1], 100);
  nlinspace(x1n, x1[0], x1[N - 1], 100);
  Matrix yi(N, N);
  for (Index i = 0; i < N; i++)
    for (Index j = 0; j < N; j++)
      yi(i, j) = std::exp(-x0[i] / 3.14 + x1[j] / 3.14);

  std::cerr << x0 << '\n';
  std::cerr << x1 << '\n';
  std::cerr << yi << '\n';

  constexpr Index order = 5;
  {
    const auto lag0 =my_interp::lagrange_interpolation_list<LagrangeInterpolation>(
        x0n, x0, order, 0.5);
    const auto lag1 = my_interp::lagrange_interpolation_list<LagrangeInterpolation>(
        x1n, x1, order, 0.5);
    const auto iwlag = interpweights(lag0, lag1);
    std::cout << x0n << '\n'
              << x1n << '\n'
              << reinterp(yi, iwlag, lag0, lag1) << '\n';
  }
}

constexpr bool is_around(Numeric x, Numeric x0, Numeric e = 1e-12) {
  return x - x0 < e and x - x0 > -e;
}

void test12() {
  constexpr int N = 10;
  constexpr matpack::matpack_constant_data<Numeric, N> y{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  constexpr matpack::matpack_constant_data<Numeric, N> y2{1,     2 * 2, 3 * 3, 4 * 4, 5 * 5,
                                      6 * 6, 7 * 7, 8 * 8, 9 * 9, 10 * 10};
  constexpr matpack::matpack_constant_data<Numeric, N> y3{
      1,         2 * 2 * 2, 3 * 3 * 3, 4 * 4 * 4, 5 * 5 * 5,
      6 * 6 * 6, 7 * 7 * 7, 8 * 8 * 8, 9 * 9 * 9, 10 * 10 * 10};
  constexpr Numeric x = 1.5;

  // Set up the interpolation Lagranges
  constexpr FixedLagrangeInterpolation<1> lin(0, x, matpack::matpack_constant_data<Numeric, 2>{1, 2});
  constexpr FixedLagrangeInterpolation<2> sqr(0, x, matpack::matpack_constant_data<Numeric, 3>{1, 2, 3});
  constexpr FixedLagrangeInterpolation<3> cub(0, x, matpack::matpack_constant_data<Numeric, 4>{1, 2, 3, 4});

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
  static_assert(is_around(sqr_sqr, x * x));
  static_assert(is_around(cub_cub, x * x * x));

  // Should output 1.5 2.5 4.5 1.5 2.25 3 1.5 2.25 3.375
  std::cout << lin_lin << ' ' << lin_sqr << ' ' << lin_cub << ' ' << sqr_lin
            << ' ' << sqr_sqr << ' ' << sqr_cub << ' ' << cub_lin << ' '
            << cub_sqr << ' ' << cub_cub << '\n';
}

void test13() {
  constexpr int N = 10;
  constexpr matpack::matpack_constant_data<Numeric, N> y{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  constexpr matpack::matpack_constant_data<Numeric, N> y2{1,     2 * 2, 3 * 3, 4 * 4, 5 * 5,
                                      6 * 6, 7 * 7, 8 * 8, 9 * 9, 10 * 10};
  constexpr matpack::matpack_constant_data<Numeric, N> y3{
      1,         2 * 2 * 2, 3 * 3 * 3, 4 * 4 * 4, 5 * 5 * 5,
      6 * 6 * 6, 7 * 7 * 7, 8 * 8 * 8, 9 * 9 * 9, 10 * 10 * 10};
  constexpr Numeric x = 2;
  constexpr Numeric dx = 1e-2;

  // Set up the interpolation Lagranges
  constexpr FixedLagrangeInterpolation<1, true> lin(
      0, x, matpack::matpack_constant_data<Numeric, 2>{1, 2});
  constexpr FixedLagrangeInterpolation<2, true> sqr(
      0, x, matpack::matpack_constant_data<Numeric, 3>{1, 2, 3});
  constexpr FixedLagrangeInterpolation<3, true> cub(
      0, x, matpack::matpack_constant_data<Numeric, 4>{1, 2, 3, 4});
  constexpr FixedLagrangeInterpolation<1, true> dlin(
      0, x + dx, matpack::matpack_constant_data<Numeric, 2>{1, 2});
  constexpr FixedLagrangeInterpolation<2, true> dsqr(
      0, x + dx, matpack::matpack_constant_data<Numeric, 3>{1, 2, 3});
  constexpr FixedLagrangeInterpolation<3, true> dcub(
      0, x + dx, matpack::matpack_constant_data<Numeric, 4>{1, 2, 3, 4});

  // Set up the interpolation weights
  constexpr auto lin_iw = interpweights(lin);
  constexpr auto sqr_iw = interpweights(sqr);
  constexpr auto cub_iw = interpweights(cub);
  constexpr auto dlin_iw = dinterpweights<0>(lin);
  constexpr auto dsqr_iw = dinterpweights<0>(sqr);
  constexpr auto dcub_iw = dinterpweights<0>(cub);
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
  static_assert(is_around(sqr_sqr, x * x));
  static_assert(is_around(cub_cub, x * x * x));
  static_assert(is_around(dlin_lin, 1));
  static_assert(is_around(dsqr_sqr, 2 * x));
  static_assert(is_around(dcub_cub, 3 * x * x));

  std::cout << lin_lin << ' ' << ' ' << sqr_sqr << ' ' << ' ' << ' ' << cub_cub
            << '\n'
            << dlin_lin << ' ' << ' ' << ' ' << ' ' << dsqr_sqr << ' ' << ' '
            << ' ' << ' ' << ' ' << ' ' << dcub_cub << '\n'
            << alt_lin_lin << ' ' << alt_sqr_sqr << ' ' << alt_cub_cub << '\n'
            << (alt_lin_lin - lin_lin) / dx << ' ' << ' ' << ' ' << ' '
            << (alt_sqr_sqr - sqr_sqr) / dx << ' ' << ' ' << ' '
            << (alt_cub_cub - cub_cub) / dx << '\n';
}

void test14() {
  Vector y{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  Vector y2{1, 2 * 2, 3 * 3, 4 * 4, 5 * 5, 6 * 6, 7 * 7, 8 * 8, 9 * 9, 10 * 10};
  Vector y3{1,         2 * 2 * 2, 3 * 3 * 3, 4 * 4 * 4, 5 * 5 * 5,
            6 * 6 * 6, 7 * 7 * 7, 8 * 8 * 8, 9 * 9 * 9, 10 * 10 * 10};

  Vector x{0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5};

  // Set up the interpolation Lagranges
  auto lin = my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<1, true>>(
      x, y);
  auto sqr = my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<2, true>>(
      x, y);
  auto cub = my_interp::lagrange_interpolation_list<FixedLagrangeInterpolation<3, true>>(
      x, y);

  // Set up the interpolation weights
  auto lin_iw = interpweights(lin);
  auto sqr_iw = interpweights(sqr);
  auto cub_iw = interpweights(cub);
  auto dlin_iw = dinterpweights<0>(lin);
  auto dsqr_iw = dinterpweights<0>(sqr);
  auto dcub_iw = dinterpweights<0>(cub);

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

void test16() {
  const Vector xasc{2, 3, 4};
  const Vector xdes{4, 3, 2};
  Vector x{2.25, 3.25, 2.35};

  auto asc = my_interp::lagrange_interpolation_list<my_interp::Lagrange<-1, true>>(x, xasc, 1, 1.0);
  auto des = my_interp::lagrange_interpolation_list<my_interp::Lagrange<-1, true>>(x, xdes, 1, 1.0);
  auto fasc = my_interp::lagrange_interpolation_list<my_interp::Lagrange<1, true>>(x, xasc, 1.0);
  auto fdes = my_interp::lagrange_interpolation_list<my_interp::Lagrange<1, true>>(x, xdes, 1.0);
  for (Index i = 0; i < x.size(); i++) {
    std::cout << x[i] << ": " << asc[i] << ' ' << '-' << ' ' << des[i] << ' '
              << " xstart [asc des]: " << xasc[asc[i].pos] << ' '
              << xdes[des[i].pos] << '\n';
  }
  std::cout << '\n';
  for (Index i = 0; i < x.size(); i++) {
    std::cout << x[i] << ": " << fasc[i] << ' ' << '-' << ' ' << fdes[i] << ' '
              << " xstart [asc des]: " << xasc[fasc[i].pos] << ' '
              << xdes[fdes[i].pos] << '\n';
  }
}

void test17() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, 0, Constant::two_pi, verbosity);
  Vector y = x;
  for (auto& f : y) f = std::sin(f);
  for (Numeric n = -Constant::two_pi; n <= 2 * Constant::two_pi; n += 0.1) {
    auto lag = my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, my_interp::cycle_0_p2pi>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, my_interp::cycle_0_p2pi>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' ' << std::sin(n) << ' '
              << interp(y, dlag_iw, lag) << ' ' << interp(y, dflag_iw, flag)
              << ' ' << std::cos(n) << '\n';
  }
}

void test18() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, -180, 180, verbosity);
  Vector y = x;
  for (auto& f : y) f = Conversion::sind(f);
  for (Numeric n = -3 * 180; n <= 3 * 180; n += 0.1) {
    auto lag = my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, my_interp::cycle_m180_p180>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, my_interp::cycle_m180_p180>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' ' << Conversion::sind(n) << ' '
              << interp(y, dlag_iw, lag) << ' ' << interp(y, dflag_iw, flag)
              << ' ' << Conversion::cosd(n) << '\n';
  }
}

template <my_interp::cycle_limit cl>
struct zero_to_half {
  static constexpr Numeric bound = cl == my_interp::cycle_limit::upper ? 0.5 : 0.0;
};

void test19() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, 0, 0.5, verbosity);
  Vector y = x;
  for (auto& f : y) f = Conversion::sind(720 * f);
  for (Numeric n = -0.5; n <= 1.5; n += 0.01) {
    auto lag =
        my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, zero_to_half>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, zero_to_half>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' ' << Conversion::sind(720 * n)
              << ' ' << interp(y, dlag_iw, lag) << ' '
              << interp(y, dflag_iw, flag) << ' ' << Conversion::cosd(720 * n)
              << '\n';
  }
}

template <my_interp::cycle_limit cl>
struct zero_dot_123_to_zero_dot_456 {
  static constexpr Numeric bound = cl == my_interp::cycle_limit::upper ? 0.456 : 0.123;
};

void test20() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, -0.123, 0.456, verbosity);
  Vector y = x;
  for (auto& f : y) f = Conversion::sind(360 / (0.456 + 0.123) * f);
  for (Numeric n = -0.5; n <= 1.5; n += 0.01) {
    auto lag = my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, zero_dot_123_to_zero_dot_456>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, zero_dot_123_to_zero_dot_456>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' '
              << Conversion::sind(360 / (0.456 + 0.123) * n) << ' '
              << interp(y, dlag_iw, lag) << ' ' << interp(y, dflag_iw, flag)
              << ' ' << Conversion::cosd(360 / (0.456 + 0.123) * n) << '\n';
  }
}

void test21() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, 0.05, 0.45, verbosity);
  Vector y = x;
  for (auto& f : y) f = Conversion::sind(720 * f);
  for (Numeric n = -0.5; n <= 1.5; n += 0.01) {
    auto lag =
        my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, zero_to_half>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, zero_to_half>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' ' << Conversion::sind(720 * n)
              << ' ' << interp(y, dlag_iw, lag) << ' '
              << interp(y, dflag_iw, flag) << ' ' << Conversion::cosd(720 * n)
              << '\n';
  }
}

void test22() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, Constant::two_pi, 0, verbosity);
  Vector y = x;
  for (auto& f : y) f = std::sin(f);
  for (Numeric n = -Constant::two_pi; n <= 2 * Constant::two_pi; n += 0.1) {
    auto lag = my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, my_interp::cycle_0_p2pi>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, my_interp::cycle_0_p2pi>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' ' << std::sin(n) << ' '
              << interp(y, dlag_iw, lag) << ' ' << interp(y, dflag_iw, flag)
              << ' ' << std::cos(n) << '\n';
  }
}

void test23() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, 0.45, 0.05, verbosity);
  Vector y = x;
  for (auto& f : y) f = Conversion::sind(720 * f);
  for (Numeric n = -0.5; n <= 1.5; n += 0.01) {
    auto lag =
        my_interp::Lagrange<-1, true, my_interp::GridType::Cyclic, zero_to_half>(0, n, x, 1);
    auto flag = my_interp::Lagrange<1, true, my_interp::GridType::Cyclic, zero_to_half>(
        0, n, x);
    auto lag_iw = interpweights(lag);
    auto flag_iw = interpweights(flag);
    auto dlag_iw = dinterpweights<0>(lag);
    auto dflag_iw = dinterpweights<0>(flag);
    std::cout << n << ' ' << interp(y, lag_iw, lag) << ' '
              << interp(y, flag_iw, flag) << ' ' << Conversion::sind(720 * n)
              << ' ' << interp(y, dlag_iw, lag) << ' '
              << interp(y, dflag_iw, flag) << ' ' << Conversion::cosd(720 * n)
              << '\n';
  }
}

void test25() {
  const Index N = 500;
  Vector x(N);
  Verbosity verbosity;
  VectorNLinSpace(x, N, 30, 150, verbosity);
  Vector y = x;
  for (auto& f : y) f = 15*f*f + f*f*f;
  for (Numeric n = 0; n <= 180; n += 0.01) {
    const auto lag =
    my_interp::Lagrange<-1, true, my_interp::GridType::CosDeg>(0, n, x, 5);
    std::cout << n << ' ' << interp(y, interpweights(lag), lag) << ' ' << interp(y, dinterpweights<0>(lag), lag) << '\n';
  }
}

template <my_interp::cycle_limit cl>
struct m2_to_p2 {
  static constexpr Numeric bound = cl == my_interp::cycle_limit::upper ? 2 : -2;
};

void test26() {
  auto f = [](Numeric x){return Math::pow2(my_interp::cyclic_clamp<m2_to_p2>(x)) - 4;};
  auto df = [](Numeric x){return 2*my_interp::cyclic_clamp<m2_to_p2>(x);};
  
  constexpr Index N = 9;
  constexpr matpack::matpack_constant_data<Numeric, N> xi{-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2};
  constexpr matpack::matpack_constant_data<Numeric, N> yi{f(-2), f(-1.5), f(-1), f(-0.5), f(0), f(0.5), f(1), f(1.5), f(2)};
  constexpr Index O1 = 3;
  
  // Test for a few values of interpolation
  {
    constexpr Numeric x = -1.75;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
    constexpr FixedLagrangeInterpolation<O1, true> lin(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(lin), lin));
    static_assert(df(x) == interp(yi, dinterpweights<0>(lin), lin));
  }
  {
    constexpr Numeric x = -1.25;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
    constexpr FixedLagrangeInterpolation<O1, true> lin(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(lin), lin));
    static_assert(df(x) == interp(yi, dinterpweights<0>(lin), lin));
  }
  {
    constexpr Numeric x = -0.25;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
    constexpr FixedLagrangeInterpolation<O1, true> lin(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(lin), lin));
    static_assert(df(x) == interp(yi, dinterpweights<0>(lin), lin));
  }
  {
    constexpr Numeric x = 1;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
    constexpr FixedLagrangeInterpolation<O1, true> lin(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(lin), lin));
    static_assert(is_around(df(x), interp(yi, dinterpweights<0>(lin), lin)));
  }
  {
    constexpr Numeric x = -2;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
    constexpr FixedLagrangeInterpolation<O1, true> lin(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(lin), lin));
    static_assert(df(x) == interp(yi, dinterpweights<0>(lin), lin));
  }
  {
    constexpr Numeric x = 0;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
    constexpr FixedLagrangeInterpolation<O1, true> lin(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(lin), lin));
    static_assert(df(x) == interp(yi, dinterpweights<0>(lin), lin));
  }
  {
    constexpr Numeric x = -4;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
  }
  {
    constexpr Numeric x = 4;
    constexpr my_interp::Lagrange<O1, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, x, xi);
    static_assert(f(x) == interp(yi, interpweights(cyc), cyc));
    static_assert(df(x) == interp(yi, dinterpweights<0>(cyc), cyc));
  }

  constexpr Index O2 = 3;
  std::cout << "x f(x) interp(x) df(x) dinterp(x)\n";
  for (Numeric X=-3; X<3; X+=0.025) {
    const my_interp::Lagrange<O2, true, my_interp::GridType::Cyclic, m2_to_p2> cyc(0, X, xi);
    const FixedLagrangeInterpolation<O2, true> lin(0, X, xi);
    std::cout << X << ' ' << f(X) << ' ' << interp(yi, interpweights(cyc), cyc) << ' ' << df(X) << ' ' << interp(yi, dinterpweights<0>(cyc), cyc) 
    << ' ' << interp(yi, interpweights(lin), lin) << ' ' << interp(yi, dinterpweights<0>(lin), lin) << '\n';
  }
}

void test27() {
  for (auto a : my_interp::enumstrs::GridTypeNames)
    std::cout << a.size() << '\n';
  for (auto a : my_interp::enumstrs::GridTypeNames)
    std::cout << a << '\n';
  for (auto a : my_interp::enumstrs::GridTypeNames)
    std::cout << a << '\n';
}

void test28() {
  using Math::pow2;
  using Math::pow3;
  using Conversion::sind;
  using Conversion::cosd;
  
  // Old Grid of pressure, latitude, and longitude
  Vector pre; VectorNLogSpace(pre, 10, 1e5, 1e-1, Verbosity()); 
  Vector lat; VectorNLinSpace(lat, 5, -80, 80, Verbosity());
  Vector lon; VectorNLinSpace(lon, 4, -170, 170, Verbosity());
  
  // New Grids (pressure will be reduced)
  Vector newpre(1, pre[pre.nelem()/2]); 
  Vector newlat; VectorNLinSpace(newlat, 4, - 90,  90, Verbosity());
  Vector newlon; VectorNLinSpace(newlon, 3, -180, 180, Verbosity());
  
  // Old Data given some values
  Tensor3 data(pre.nelem(), lat.nelem(), lon.nelem());
  for (Index i=0; i<pre.nelem(); i++) {
    for (Index j=0; j<lat.nelem(); j++) {
      for (Index k=0; k<lon.nelem(); k++) {
        data(i, j, k) =
          std::log(pre[i]) * sind(lat[j]) * pow3(cosd(lon[k]));
      }
    }
  }
  
  // Create Lagrange interpolation values
  const ArrayOfLagrangeLogInterpolation lag_pre =
    my_interp::lagrange_interpolation_list<LagrangeLogInterpolation>(newpre, pre, 2, 1e9);
  const auto lag_lat =
    my_interp::lagrange_interpolation_list<my_interp::Lagrange<-1, false, my_interp::GridType::SinDeg>>(newlat, lat, 3, 1e9);
  const auto lag_lon =
    my_interp::lagrange_interpolation_list<my_interp::Lagrange<-1, false, my_interp::GridType::Cyclic, my_interp::cycle_m180_p180>>(newlon, lon, 1, 1e9);
  
  // Create the interpolation weights
  const auto lag_iw = interpweights(lag_pre, lag_lat, lag_lon);
  
  // Create and reduce the new data to dim 1 and 2
  const Matrix newdata = reinterp(data, lag_iw, lag_pre, lag_lat, lag_lon).reduce_rank<1, 2>();
  
  // Print the data
  std::cout << data(pre.nelem()/2, joker, joker) 
            << '\n' << '\n' << newdata << '\n';
}

int main() { test09(); }

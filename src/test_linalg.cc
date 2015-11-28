/* Copyright (C) 2002-2012 Claudia Emde <claudia.emde@dlr.de>

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
   USA.

*/

/*!
  \file   test_linalg.cc
  \author Claudia Emde <claudia.emde@dlr.de>
  \date   Thu May  2 14:37:57 2002

  \brief  Test for the linear algebra functions.


*/

#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include "lin_alg.h"
#include "make_vector.h"
#include "array.h"
#include "test_utils.h"

using std::cout;
using std::endl;
using std::abs;
using std::setw;

//!
/*! The function tests the LU-decompusition method for solving a
  1D linear equation
  system. It uses the functions 'ludcmp' and  'lubacksub'.
*/
void test_lusolve1D(void)
{
  Matrix a(1,1);
  ArrayOfIndex indx(1);
  Matrix orig(1,1);
  Matrix b(1,1);

 /* Assign test-matrix element. */
  a(0,0) = 3;

 /* ------------------------------------------------------------------------
     Test the function ludcmp.
     ----------------------------------------------------------------------- */
  cout << "\n LU decomposition test \n";
  cout << "initial matrix: \n";

  cout << " " << a(0,0) << endl;

   /* input: Test-matrix a,
      output: Decomposed matrix b (includes upper and lower triangle, cp.
      Numerical Recipies)
      and index which gives information about pivoting. */
  ludcmp(b, indx, a);

  cout << "\n after decomposition: ";
  cout << b(0,0) << endl;

  /* Seperate b into the two triangular matrices. */
  Matrix l(1,1,0.0);
  Matrix u(1,1,0.0);
  Matrix lu(1,1,0.0);

  l(0,0) = 1.0;
  u(0,0) = b(0,0);

  /*-------------------------------------------------------------------
    end of ludcmp test
     ------------------------------------------------------------------*/
  /*--------------------------------------------------------------------
    test backsubstitution routine lubacksub
    -------------------------------------------------------------------*/

   Vector c(1);
   c[0] = 6;
   cout << indx[0] << "  " << c[0] << endl;

   Vector x(1);
   lubacksub(x, b, c, indx);

   cout << "\n solution vector x: ";
   cout << x[0] << endl;

}


//!
/*! The function tests the LU-decompusition method for solving a
  linear equation
  system. It uses the functions 'ludcmp' and  'lubacksub'.
*/
void test_lusolve4D(void)
{
  Matrix a(4,4);
  ArrayOfIndex indx(4);
  Matrix orig(4,4);
  Matrix b(4,4);

  /* Assign test-matrix elements. */

  a(0,0) = 1;
  a(0,1) = 3;
  a(0,2) = 5;
  a(0,3) = 6;
  a(1,0) = 2;
  a(1,1) = 3;
  a(1,2) = 4;
  a(1,3) = 4;
  a(2,0) = 1;
  a(2,1) = 2;
  a(2,2) = 5;
  a(2,3) = 1;
  a(3,0) = 7;
  a(3,1) = 2;
  a(3,2) = 4;
  a(3,3) = 3;


  /* ------------------------------------------------------------------------
     Test the function ludcmp.
     ----------------------------------------------------------------------- */


  cout << "\n LU decomposition test \n";
  cout << "initial matrix: \n";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << a(i,j);
    }
   cout << "\n";

   /* input: Test-matrix a,
      output: Decomposed matrix b (includes upper and lower triangle, cp.
      Numerical Recipies)
      and index which gives information about pivoting. */
   ludcmp(b, indx, a);

  cout << "\n after decomposition";
  for( Index i = 0; i<4; i++)
    {      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << b(i,j);
    }
  cout << "\n";

  /* Seperate b into the two triangular matrices. */
  Matrix l(4,4,0.0);
  Matrix u(4,4,0.0);
  Matrix lu(4,4,0.0);


  for(Index i = 0; i<4; i++) l(i,i) = 1.0;
  l(1,0) = b(1,0);
  l(2,Range(0,2)) = b(2, Range(0,2));
  l(3,Range(0,3)) = b(3, Range(0,3));

  cout << "\n Matrix L";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << l(i,j);
    }
  cout << "\n";


  u(0,Range(0,4)) = b(0,Range(0,4));
  u(1,Range(1,3)) = b(1,Range(1,3));
  u(2,Range(2,2)) = b(2,Range(2,2));
  u(3,Range(3,1)) = b(3,Range(3,1));


   cout << "\n Matrix U";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << u(i,j);
    }
  cout << "\n";


  /* Test, if LU = a. */
  mult(lu,l,u);

  cout << "\n product L*U";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << lu(i,j);
    }
   cout << "\n";

   /*-------------------------------------------------------------------
     end of ludcmp test
     ------------------------------------------------------------------*/


   /*--------------------------------------------------------------------
     test backsubstitution routine lubacksub
     -------------------------------------------------------------------*/

   Vector c(4);
   c[0] = 2;
   c[1] = 5;
   c[2] = 6;
   c[3] = 7;

   cout << "\n  vector indx";
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << indx[i] << "  " << c[i];
     }
   cout << endl;

   Vector x(4);
   lubacksub(x, b, c, indx);

   cout << "\n solution vector x" << endl;
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << x[i];
     }
   cout << endl;

   cout << "\n test solution LU*x";
     Vector y(4);
     mult(y,a,x);
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << y[i];
     }
   cout <<"\n";
}

//! Test ludcmp and lubacksub by solving a linear system of equations.
/*!
  Generates a random, square (n,n)-matrix A and a length-n vector x0 and
  solves A*x = A*x0. The maximum relative, component-wise error in abs(x0 - x)
  is written to standard out. The numbers of tests performed is controlled by
  the parameter ntests. If verbose == true, also A, x0 and x are written to
  standard out.

  \param[in] ntests Number of tests to be performed.
  \param[in] dim    Dimensionality of the equation system.
  \param[in] verbose Controls verbosity of output. If true, for each test the
                    matrix A and the vectors x0 and x are written to standard out.
                    Otherwise only the maximum relative error in each component of
                    x is written out.
  \return void
*/
void test_solve_linear_system( Index ntests,
                               Index dim,
                               bool verbose )
{
    Matrix A(dim,dim);
    Matrix LU(dim,dim);
    Vector x0(dim);
    Vector x(dim);
    Vector b(dim);
    ArrayOfIndex indx(dim);

    // initialize random seed
    srand((unsigned int) time(0));

    cout << endl << endl << "Testing linear system solution: n = " << dim;
    cout << ", ntests = " << ntests << endl;
    cout << endl << setw(10) << "Test no." << setw(20) << "lubacksub(...)";
    cout << setw(20) << "solve(...)" << endl << endl;

    for (Index i = 0; i < ntests; i++)
    {

        // Generate linear system, make sure the determinant
        // is non-zero.
        random_fill_matrix_pos_def( A, 10, false );
        random_fill_vector( x0, 10, false );
        mult( b, A, x0);

        // Test ludcmp/lubacksub.
        ludcmp(LU, indx, A);
        lubacksub(x, LU, b, indx);

        Numeric err = 0.0;
        err = get_maximum_error( x, x0, true );

        cout  << setw(10) << i << setw(20) << err;

        // Test solve.
        solve(x, A, b);

        err = get_maximum_error( x, x0, true );

        cout  << setw(20) << err << endl;

        if (verbose)
        {
            cout << endl;
            cout << "A:" << endl << A << endl << endl;
            cout << "x0:" << endl << x0 << endl << endl;
            cout << "x:" << endl << x << endl << endl;
            cout << "Permutation Vector:" << endl << indx << endl;
            cout << endl;
        }
    }
}

//! Test matrix inversion.
/*!
  Generates a random, square (n,n)-matrix A and computes its inverse Ainv and
  I = Ainv*A. The maximum absolute error in I with respect to the identity matrix
  is written to standard out. The number of tests performed is controlled by the
  parameter ntests. If verbose == true, also A, Ainv and A*Ainv are written to
  standard out.

  \param[in] ntests Number of tests to be performed.
  \param[in] dim    Size of matrix A.
  \param[in] verbose Controls verbosity of output. If true, for each test the
                    matrices A, Ainv and I = Ainv*A are written to standard out.
                    Otherwise only the maximum absolute error in I = Ainv*A w.r.t.
                    the identity matrix is written out.
  \return void
*/

void test_inv(Index ntests, Index dim, bool verbose=false)
{
    Matrix A(dim,dim);
    Matrix Ainv(dim,dim);
    Matrix I0(dim,dim);
    id_mat(I0);
    Matrix I(dim,dim);

    // initialize random seed
    srand((unsigned int) time(0));

    cout << endl << endl << "Testing matrix inversion: n = " << dim;
    cout << ", ntests = " << ntests << endl << endl;
    cout << setw(10) << "Test no." << setw(20) << "Max. rel. error" << endl << endl;

    for (Index i = 0; i < ntests; i++)
    {
        // Generate random matrix, make sure the determinant
        // is non-zero.
        random_fill_matrix_pos_def(A, 10, false);

        inv(Ainv,A);
        mult(I, Ainv, A);

        Numeric err = get_maximum_error( I, I0, false );
        // Print results.
        cout  << setw(10) << i << setw(20) << err << endl;

        if (verbose)
        {
            cout << endl;
            cout << "A:" << endl << A << endl << endl;
            cout << "Ainv:" << endl << Ainv << endl << endl;
            cout << "A*Ainv:" << endl << I << endl << endl;
        }

    }
}

//! Test for the matrix exponential function (4D matrix)
/*!

 */
void test_matrix_exp4D(void)
{
  Matrix A(4,4);
  Matrix F(4,4);
  A(0,0) = 1;
  A(0,1) = 3;
  A(0,2) = 5;
  A(0,3) = 6;
  A(1,0) = 2;
  A(1,1) = 3;
  A(1,2) = 4;
  A(1,3) = 4;
  A(2,0) = 1;
  A(2,1) = 2;
  A(2,2) = 5;
  A(2,3) = 1;
  A(3,0) = 7;
  A(3,1) = 2;
  A(3,2) = 4;
  A(3,3) = 3;

  /* set parameter for accuracy */
  Index q = 8;

  /*execute matrix exponential function*/
  matrix_exp(F,A,q);


  cout << "\n Exponential of Matrix K";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << F(i,j);
    }
  cout << "\n";
}


//! Test for the matrix exponential function (3D matrix)
/*!

 */
void test_matrix_exp1D(void)
{
  Matrix A(1,1);
  Matrix F(1,1);
  A(0,0) = 5;

  /* set parameter for accuracy */
  Index q = 8;

  /*execute matrix exponential function*/
  matrix_exp(F,A,q);

   cout << "\n Exponential of Matrix A:\n";
   cout << F(0,0);
   cout <<"\n";
}

//! Test for the matrix exponential function (3D matrix)
/*!

 */
void test_matrix_exp3D(void)
{
  Matrix A(3,3);
  Matrix F(3,3);
  A(0,0) = 1;
  A(0,1) = 3;
  A(0,2) = 5;
  A(1,0) = 2;
  A(1,1) = 3;
  A(1,2) = 4;
  A(2,0) = 1;
  A(2,1) = 2;
  A(2,2) = 5;


  /* set parameter for accuracy */
  Index q = 8;

  /*execute matrix exponential function*/
  matrix_exp(F,A,q);


  cout << "\n Exponential of Matrix A";
  for( Index i = 0; i<3; i++)
    {
      cout << "\n";
      for (Index j = 0; j<3; j++)
        cout << " " << F(i,j);
    }
   cout << "\n";
}

int main(void)
{
  //test_lusolve4D();
  test_inv( 20, 1000 );
  test_solve_linear_system( 20, 1000, false );
  // test_matrix_exp1D();
  return(0);
}

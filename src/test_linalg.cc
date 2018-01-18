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
#include "array.h"
#include "test_utils.h"
#include <random>

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

void test_real_diagonalize(Index ntests, Index dim)
{
    Matrix A(dim,dim),F1(dim,dim),F2(dim,dim),tmp1(dim,dim),tmp2(dim,dim),P(dim,dim);
    Vector Wr(dim),Wi(dim);
    
    const Matrix ZEROES(dim,dim,0);
    
    // initialize random seed
    srand((unsigned int) time(0));
    
    cout << endl << endl << "Testing diagonalize: n = " << dim;
    cout << ", ntests = " << ntests << endl;
    cout << setw(10) << "Test no." << setw(25) << "Max. rel. expm error";
    cout << setw(25) << "Max. abs. P^-1*A*P-W" << endl<< endl;
    
    for (Index i = 0; i < ntests; i++)
    {
        
        // Generate a matrix that does not have complex answers...
        random_fill_matrix_symmetric( A, 10, false );
        
        // Use the two expm methods to test that diagonalize works
        matrix_exp(F1,A,10);
        matrix_exp2(F2,A);
        
        Numeric err1 = 0.0,err2 = 0.0;
        err1 = get_maximum_error( F1, F2, true );
        
        
        // diagonalize directly to test that P^-1*A*P gives a diagonal matrix
        diagonalize(P,Wr,Wi,A);
        
        // P^-1*A*P
        inv(tmp1,P);
        mult(tmp2,tmp1,A);
        mult(tmp1,tmp2,P);
        
        // Minus W as diagonal matrix
        for(Index j=0;j<dim;j++)
        {
            tmp1(j,j)-=Wr[j];
        }
        
        err2 = get_maximum_error( ZEROES, tmp1, false );
        
        cout  << setw(10) << i << setw(25) << err1<< setw(25) << err2 << endl;
        
    }
    
}

void test_complex_diagonalize(Index ntests, Index dim)
{
    ComplexMatrix A(dim,dim),tmp1(dim,dim),tmp2(dim,dim),P(dim,dim);
    ComplexVector W(dim);
    
    const ComplexMatrix ZEROES(dim,dim,0);
    
    // initialize random seed
    srand((unsigned int) time(0));
    
    cout << endl << endl << "Testing diagonalize: n = " << dim;
    cout << ", ntests = " << ntests << endl;
    cout << setw(10) << "Test no.";
    cout << setw(25) << "Max. abs. P^-1*A*P-W" << endl<< endl;
    
    for (Index i = 0; i < ntests; i++)
    {
        
        // Generate a matrix that does not have complex answers...
        random_fill_matrix_symmetric( A, 10, false );
        
        Numeric err = 0.0;
        
        // diagonalize directly to test that P^-1*A*P gives a diagonal matrix
        diagonalize(P,W,A);
        
        // P^-1*A*P
        inv(tmp1,P);
        mult(tmp2,tmp1,A);
        mult(tmp1,tmp2,P);
        
        // Minus W as diagonal matrix
        for(Index j=0;j<dim;j++)
        {
            tmp1(j,j)-=W[j];
        }
        
        err = get_maximum_error( ZEROES, tmp1, false );
        
        cout  << setw(10) << i << setw(25) << err << endl;
        
    }
}

void test_matrix_exp_propmat(Index nruns, Index ndiffs)
{
  for(Index i = 0; i < nruns; i++)
  {
    Numeric a = ((Numeric)rand() / RAND_MAX), 
            b = ((Numeric)rand() / RAND_MAX), 
            c = ((Numeric)rand() / RAND_MAX), 
            d = ((Numeric)rand() / RAND_MAX), 
            u = ((Numeric)rand() / RAND_MAX), 
            v = ((Numeric)rand() / RAND_MAX), 
            w = ((Numeric)rand() / RAND_MAX);
  
    Matrix A(4, 4), F0(4, 4), F1(4, 4), F2(4, 4);
    
    A(0, 0) = A(1, 1) = A(2, 2) = A(3, 3) = a;
    A(0, 1) = A(1, 0) = b;
    A(0, 2) = A(2, 0) = c;
    A(0, 3) = A(3, 0) = d;
    A(1, 2) =  u;  
    A(2, 1) = -A(1, 2);
    A(1, 3) = v;
    A(3, 1) = -A(1, 3);
    A(2, 3) =  w;
    A(3, 2) = -A(2, 3);
    
    A *= -1.0;
    
    Tensor3 dF1(ndiffs, 4, 4), dF2(ndiffs, 4, 4), dF3(ndiffs, 4, 4), dF4(ndiffs, 4, 4), dA1(ndiffs, 4, 4), dA2(ndiffs, 4, 4);
    
    Numeric da1 = ((Numeric)rand() / RAND_MAX), 
    db1 = ((Numeric)rand() / RAND_MAX), 
    dc1 = ((Numeric)rand() / RAND_MAX), 
    dd1 = ((Numeric)rand() / RAND_MAX), 
    du1 = ((Numeric)rand() / RAND_MAX), 
    dv1 = ((Numeric)rand() / RAND_MAX), 
    dw1 = ((Numeric)rand() / RAND_MAX);
    
    Numeric da2 = ((Numeric)rand() / RAND_MAX), 
    db2 = ((Numeric)rand() / RAND_MAX), 
    dc2 = ((Numeric)rand() / RAND_MAX), 
    dd2 = ((Numeric)rand() / RAND_MAX), 
    du2 = ((Numeric)rand() / RAND_MAX), 
    dv2 = ((Numeric)rand() / RAND_MAX), 
    dw2 = ((Numeric)rand() / RAND_MAX);
    
    dA1(joker, 0, 0) = dA1(joker, 1, 1) = dA1(joker, 2, 2) = dA1(joker, 3, 3) = da1;
    dA1(joker, 0, 1) = dA1(joker, 1, 0) = db1;
    dA1(joker, 0, 2) = dA1(joker, 2, 0) = dc1;
    dA1(joker, 0, 3) = dA1(joker, 3, 0) = dd1;
    dA1(joker, 1, 2) =  du1;
    dA1(joker, 2, 1) = -du1;
    dA1(joker, 1, 3) =  dv1;
    dA1(joker, 3, 1) = -dv1;
    dA1(joker, 2, 3) =  dw1;
    dA1(joker, 3, 2) = -dw1;
    
    dA1 *= -1.0;
    
    dA2(joker, 0, 0) = dA2(joker, 1, 1) = dA2(joker, 2, 2) = dA2(joker, 3, 3) = da2;
    dA2(joker, 0, 1) = dA2(joker, 1, 0) = db2;
    dA2(joker, 0, 2) = dA2(joker, 2, 0) = dc2;
    dA2(joker, 0, 3) = dA2(joker, 3, 0) = dd2;
    dA2(joker, 1, 2) =  du2;
    dA2(joker, 2, 1) = -du2;
    dA2(joker, 1, 3) =  dv2;
    dA2(joker, 3, 1) = -dv2;
    dA2(joker, 2, 3) =  dw2;
    dA2(joker, 3, 2) = -dw2;
    
    dA2 *= -1.0;
    
    special_matrix_exp_and_dmatrix_exp_dx_for_rt(F1, dF1, dF2, A, dA1, dA2, 10);
    
    cayley_hamilton_fitted_method_4x4_propmat_to_transmat__explicit(F1, dF1, dF2, A, dA1, dA2);
    cayley_hamilton_fitted_method_4x4_propmat_to_transmat__eigen(F2, dF3, dF4, A, dA1, dA2);
    
    matrix_exp(F0, A);
    cayley_hamilton_fitted_method_4x4_propmat_to_transmat__eigen(F1, A);
    cayley_hamilton_fitted_method_4x4_propmat_to_transmat__explicit(F2, A);

    std::cout<<MapToEigen(F1)<<"\n\n";
    std::cout<<MapToEigen(F2)<<"\n\n";
    std::cout<<MapToEigen(F1)-MapToEigen(F2)<<"\n\n";
    std::cout<<MapToEigen(dF1(0, joker, joker))<<"\n\n";
    std::cout<<MapToEigen(dF3(0, joker, joker))<<"\n\n";
    std::cout<<MapToEigen(dF3(0, joker, joker))-MapToEigen(dF1(0, joker, joker))<<"\n\n";
    std::cout<<MapToEigen(dF2(0, joker, joker))<<"\n\n";
    std::cout<<MapToEigen(dF4(0, joker, joker))<<"\n\n";
    std::cout<<MapToEigen(dF4(0, joker, joker))-MapToEigen(dF2(0, joker, joker))<<"\n\n";
  }
}

int main(void)
{
  // test_lusolve4D();
  // test_inv( 20, 1000 );
  // test_solve_linear_system( 20, 1000, false );
  // test_matrix_exp1D();
  //  test_real_diagonalize(20,100);
  //  test_complex_diagonalize(20,100);
  test_matrix_exp_propmat(100000, 10);
  return(0);
}

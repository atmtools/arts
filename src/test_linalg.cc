/* Copyright (C) 2002-2008 Claudia Emde <claudia.emde@dlr.de>
                      
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
#include "lin_alg.h"
#include "make_vector.h"
#include "array.h"


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

  cout << "\n after decomposition";
  cout << " " << b(0,0) << endl;
  
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
   
   cout << "\n solution vector x";
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

   Vector x(4);
   lubacksub(x, b, c, indx);

   cout << "\n solution vector x" << endl;
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << x[i];
     }
  
   cout << "\n test solution LU*x";
     Vector y(4);
   mult(y,lu,x);
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << y[i];
     }
   cout <<"\n";   
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
  test_lusolve1D();
  // test_matrix_exp1D();
  return(0);
}

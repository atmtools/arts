/* Copyright (C) 2002 Claudia Emde <claudia@sat.physik.uni-bremen.de>
                      
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
  \author Claudia Emde <claudia@sat.physik.uni-bremen.de>
  \date   Thu May  2 14:37:57 2002
  
  \brief  Test for the linear algebra functions.
  
  
*/

#include "lin_alg.h"

int main(void)
{
  Matrix a(4,4);
  ArrayOfIndex indx(4);
  Numeric d;
  Matrix orig(4,4);

  // assign matrix elements
  
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

  orig = a;

  cout << "\n LU decomposition test \n";
  cout << "initial matrix: \n";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << a(i,j);
    }
   cout << "\n";

  ludcmp(a, indx, d);

  cout << "\n after decomposition";
  for( Index i = 0; i<4; i++)
    {      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << a(i,j);
    }
  cout << "\n";

  Matrix l(4,4,0.0);
  Matrix u(4,4,0.0);
  Matrix lu(4,4,0.0);
  
  
  for(Index i = 0; i<4; i++) l(i,i) = 1.0;
  l(1,0) = a(1,0);
  l(2,Range(0,2)) = a(2, Range(0,2));
  l(3,Range(0,3)) = a(3, Range(0,3));
  
  cout << "\n Matrix L";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << l(i,j);
    }
  cout << "\n";


  u(0,Range(0,4)) = a(0,Range(0,4));
  u(1,Range(1,3)) = a(1,Range(1,3));
  u(2,Range(2,2)) = a(2,Range(2,2));
  u(3,Range(3,1)) = a(3,Range(3,1));


   cout << "\n Matrix U";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << u(i,j);
    }
  cout << "\n";

  mult(lu,l,u);

  cout << "\n product L*U";
  for( Index i = 0; i<4; i++)
    {
      cout << "\n";
      for (Index j = 0; j<4; j++)
        cout << " " << lu(i,j);
    }
   cout << "\n";


   //test backsubstitution routine

   Vector b(4);
   b[0] = 2;
   b[1] = 5;
   b[2] = 6;
   b[3] = 7;

   cout << "\n  vector indx";
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << indx[i] << "  " << b[i];
     }


   lubacksub(a,b,indx);

   cout << "\n solution vector x";
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << b[i];
     }
  
   cout << "\n test solution LU*x";
     Vector y(4);
   mult(y,lu,b);
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << y[i];
     }
   
   
   cout << "\n";
 
   /*============================================================
     Now test the function lusolve, which combines LU-
     decomposition and backsubstitution to solve linear equation
     systems.
     ============================================================*/

   // assign matrix elements
  
   Matrix K(4,4);
   
   K(0,0) = 1;
   K(0,1) = 3;
   K(0,2) = 5;
   K(0,3) = 6;
   K(1,0) = 2;
   K(1,1) = 3;
   K(1,2) = 4;
   K(1,3) = 4;
   K(2,0) = 1;
   K(2,1) = 2;
   K(2,2) = 5;
   K(2,3) = 1;
   K(3,0) = 7;
   K(3,1) = 2;
   K(3,2) = 4;
   K(3,3) = 3;

  

   // assign vector elements
   Vector y2(4);
   
   y2[0] = 2;
   y2[1] = 5;
   y2[2] = 6;
   y2[3] = 7;
   
   // solution vector
   Vector x(4);

   lusolve(x,K,y2);
   cout << "\n Test lusolve function:";
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << x[i];
     }
    cout << "\n";
   
   // test solution:
   Vector y_test(4);
   
   mult(y_test,K,x);

   cout << "\n Test lusolve function:";
   for (Index i=0; i<4; i++)
     {
       cout << "\n";
       cout << y_test[i];
     }
    cout << "\n";
   
     
   return 0;
}

/**
   This file is at the moment obsolete. The reason why I have not
   removed it is that maybe I can translate the examples and add them to
   demonstrate_vecmat.cc. 
   
   See demonstrate_vecmat.cc for usage examples of matrix/vector
   functions and mathematical functions.
   
   \date   2001-01-08
   \author Stefan Buehler
*/


#include "arts.h"
#include "vecmat.h"
#include "math_funcs.h"          

void test_math()
{
  Numeric     a = 1;
  VECTOR      x(4),            // Gives a vector of length 4
              y(3,5.0),        // Gives the vector [5,5,5]
              z(3,"1 2 3");    // Gives the vector [1,2,3]
  MATRIX      A(9,7),          // Gives a matrix with 9 rows and 7 columns
              B(2,3,2.3),      // 2 rows, 3 columns, all elements 2.3
              C(2, 3,"1 2 3"   // Gives the matrix |1,2,3|
                     "4 5 6"); //                  |4,5,6|

  // Reallocation
  x.resize(9);            // Gives a vector of length 9
  A = MATRIX(4,5);        // Gives a matrix with 4 rows and 5 columns
  A = transpose(B);       // Transpose of matrix

  // To get dimensions
  x.size();               // Gives the length of the vector
  A.nrows();              // Gives the number of rows of the matrix
  A.ncols();              // Gives the number of columns of the matrix

  // ??
  x = 1.1;                // Sets all elements of x to 1.1
  x = y;                  // ?? 
  A = 1.1;                // Sets all elements of A to 1.1
  A = B;                  // ??
  
  // Sum and difference
  x = y + 3.3;            // Elementwise vector-scalar operations
  x = 3.3 + y;
  x = y - 3.3;    
  x = 3.3 - y;
  A = B + 1.1;            // Elementwise matrix-scalar operations
  A = 1.1 + B;
  A = B - 3.3;    
  A = 3.3 - B;    
  x = y + z;              // Elementwise sum of two vectors
  x = y - z;              // Elementwise difference between two vectors
  A = B + C;              // Elementwise sum of two matrices
  A = B - C;              // Elementwise difference between two matrices

  // Multiplication and divison
  x = y * 3.3;            // Elementwise vector-scalar operations
  x = 3.3 * y;
  x = y / 3.3;    
  x = 3.3 / y;
  A = B * 1.1;            // Elementwise matrix-scalar operations
  A = 1.1 * B;
  A = B / 3.3;    
  A = 3.3 / B;    
  x = emult(y,z);         // Elementwise multiplication of two vectors
  x = ediv(y,z);          // Elementwise division between two vectors
  A = emult(B,C);         // Elementwise multiplication of two matrices
  A = ediv(B,C);          // Elementwise division between two matrices
  x = B*y;                // Matrix-vector multiplication
  A = B*transpose(C);     // Matrix-matrx multiplication
  a = dot_prod(y,z);      // Scalar product of two vectors.

  // Conversion between vectors and matrices
  x = row(2,B);        // Makes a vector of row 2 of B.
  row(x,2,B);
  x = col(2,B);        // Makes a vector of column 2 of B.
  col(x,2,B);
  A = to_matrix(y);    // Makes a column matrix of a vector
  A = to_matrix(y);
  x = to_vector(A);    // Converts a matrix to a vector. The matrix can
  to_vector(x,A);      // either be a column (n x 1) or row (1 x n) vector.
 
  // Regions
  // ??
}

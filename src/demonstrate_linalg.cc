#include <limits.h>
#include "arts.h"
#include "linalg.h"

int
main()
{
  INDEX n = 5;
  VEC x(n);
  ARRAY<string> a(n);
  ARRAY<bool>   b(n);

  cout << "Int limits: " << INT_MIN << ", " << INT_MAX << "\n";
  cout << "Long limits: " << LONG_MIN << ", " << LONG_MAX << "\n";

  for (INDEX i = 0; i < n; ++i)
    {
      x[i] = i;
      a[i] = "aa";
      b[i] = true;
    }
  a[3] = "b";

  cout << "Size of x: " << x.size() << "\n";
  cout << "x: ";
  print_vector(cout,x);

  cout << "One-norm of x: ";
  cout << one_norm(x) << "\n";

  x.resize(n+1);
  cout << "New size of x: " << x.size() << "\n";
  cout << "x: ";
  print_vector(cout,x);

  erase_vector_element(x,3);
  cout << "After erasing element 3: " << x.size() << "\n";
  cout << "x: ";
  print_vector(cout,x);

  cout << "a: ";
  print_vector(cout,a);

  cout << "b: ";
  print_vector(cout,b);
  //  cout << one_norm(a) << "\n";

  INDEX i = max_index(x);
  cout << "max_index(x): " << i << "\n";

  // This proves that you can use ARRAYs like any STL container:
  cout << "Index of a==b: "
	    << std::find(a.begin(), a.end(), "b") - a.begin()
	    << "\n";

  // Subrange:
  cout << "x(2:5): ";
  print_vector(cout,x(2,5));

  // Same with explicit storage of the subrange:
  VEC::subrange_type xs(x(2,5));
  cout << "Same: ";
  print_vector(cout,xs);



  // Demonstrate matrix features:

  MAT A(n, 2*n);

  for (INDEX i = 0; i < A.nrows(); ++i)
    for (INDEX j = 0; j < A.ncols(); ++j)
      A[i][j] = i+j;

  cout << "A: ";
  print_all_matrix(cout,A);


  // Demonstrate row selection:
  // A[3] is of type MAT::OneD, therefore ar = A[3] does not
  // work. However, all vector algorithms work for A[3], for example
  // copy() as in the example below.
  // Fazit, you can use A[3] like a vector, but you cannot assign it
  // to a vector directly, you have to use copy.
  VEC ar(A.ncols());
  copy(A[3],ar);
  cout << "ar: ";
  print_vector(cout,ar);

  // Demonstrate columns selection:
  // We can use the helper function columns to aparently transpose the
  // matrix. Also this demonstrate that the vector algorithm
  // print_vector works on the thus selected beast:
  cout << "A(*,4): ";
  print_vector(cout,columns(A)[4]);

  // Demonstrate submatrix selection:
  cout << "A(0:3,0:2): ";
  print_all_matrix(cout,A.sub_matrix(0,3,0,2));

  // How to copy a submatrix to a matrix:
  MAT::submatrix_type S1 = A.sub_matrix(0,3,0,2);

  // Test [][] indexing for submatrix:
  cout << "S1[1][1]: " << S1[1][1] << "\n";

  MAT S2(S1.nrows(),S1.ncols());
  copy(S1,S2);
  cout << "The same copied: ";
  print_all_matrix(cout,S2);


  // Sparse Matrices
  //===================

  SPARSE_MAT S(n,2*n);

  // How many nonzero elements?
  cout << "S.nnz(): " << S.nnz() << "\n";

  S[1][2]  = 1;
  S[3][5]  = 2;
  
  cout << "S.nnz(): " << S.nnz() << "\n";
  cout << "S: ";
  print_all_matrix(S);

  MAT At(2*n,n);
  copy(trans(A),At);
  cout << "At: ";
  print_all_matrix(At);

  // Calculate S*transpose(A):
  MAT Z(n,n);
  cout << "Capacity of Z: " << Z.capacity() << "\n";

  mult(S,At,Z);
  //  copy(trans(S),Z);

  cout << "Z: ";
  print_all_matrix(Z);
}


#include <limits.h>
#include "arts.h"
#include "vecmat.h"
#include "math_funcs.h"
#include "make_vector.h"

int
main()
{
  INDEX n = 5;
  ARRAY<string> a(n);
  ARRAY<bool>   b(n);

  cout << "Int limits: " << INT_MIN << ", " << INT_MAX << "\n";
  cout << "Long limits: " << LONG_MIN << ", " << LONG_MAX << "\n";

  // Demonstrate vector features:
  cout << "\n"
       << "----------------\n";
  cout << "Vector features:\n";
  cout << "----------------\n";

  // Define and initialize:
  VECTOR x(n,3);
  cout << "x: " << x << "\n";

  for (INDEX i = 0; i < n; ++i)
    {
      x[i] = i;
      a[i] = "aa";
      b[i] = true;
    }
  a[3] = "Stefan";

  cout << "Size of x: " << x.size() << "\n";
  cout << "x: " << x << "\n";

  // Subrange:
  cout << "x(2:5): " << x(2,5) << "\n";

  cout << "One-norm of x: ";
  cout << one_norm(x) << "\n";

  resize(x,n+1);
  cout << "New size of x: " << x.size() << "\n";
  cout << "x: " << x << "\n";

  for (INDEX i = 0; i < x.size(); ++i)
    x[i] = i;
  
  cout << "Re-initialize, x: " << x << "\n";

  erase_vector_element(x,3);
  cout << "Size after erasing element 3: " << x.size() << "\n";
  cout << "x: " << x << "\n";

  cout << "a: " << a << "\n";

  // Copy a:
  ARRAY<string> aa(n);
  copy(a,aa);
  cout << "aa: " << aa << "\n";
  
  cout << "b: " << b << "\n";
  //  cout << one_norm(a) << "\n";

  INDEX i = max_index(x);
  cout << "max_index(x): " << i << "\n";

  // This proves that you can use ARRAYs like any STL container:
  cout << "Index of a==Stefan: "
	    << std::find(a.begin(), a.end(), "Stefan") - a.begin()
	    << "\n";

  // Subrange:
  cout << "x(2:5): " << x(2,5) << "\n";

  // Same with explicit storage of the subrange:
  VECTOR::subrange_type xs(x(2,5));
  cout << "Same: " << xs << "\n";

  x(2,5)[1] = 9999;

  cout << "After changing middle element of sub range, x: " << x << "\n";

  // Demonstrate matrix features:
  cout << "\n"
       << "----------------\n";
  cout << "Matrix features:\n";
  cout << "----------------\n";

  MATRIX A(n, 2*n);

  {
    Numeric s = 0;
    for (INDEX i = 0; i < A.nrows(); ++i)
      for (INDEX j = 0; j < A.ncols(); ++j)
	{
	  A[i][j] = s;
	  ++s;
	}
  }

  cout << "A:\n" << A << "\n";


  // Demonstrate row selection:
  // A[3] is of type MATRIX::OneD, therefore assigning it to a vector
  // with `=' will not work. However, you can use all vector
  // algorithms on it, for example copy(). See also next example
  // below. Here, we use the print_vector() algorithm directly:
  cout << "A(3,*): ";
  print_vector(cout,A[3]);
  cout << "\n";
  
  // Demonstrate columns selection:
  // We can use the helper function columns to aparently transpose the
  // matrix. 
  // As an alternative to the example above, here I use copy() to copy
  // the selected column to a vector, then I can use the << operator
  // on the vector.
  // Fazit, you can use A[3] and columns(A)[4] like vectors, but you
  // cannot assign them to a vector directly, you have to use copy.
  VECTOR ac(A.nrows());
  copy(columns(A)[4],ac);
  cout << "A(*,4): " << ac << "\n";

  // Demonstrate submatrix selection:
  cout << "A(0:3,0:2):\n" << A.sub_matrix(0,3,0,2) << "\n";

  // How to copy a submatrix to a matrix:
  MATRIX::submatrix_type S1 = A.sub_matrix(0,3,0,2);

  // Test [][] indexing for submatrix:
  cout << "S1[1][1]:\n" << S1[1][1] << "\n";

  MATRIX S2(S1.nrows(),S1.ncols());
  copy(S1,S2);
  cout << "A copy of S1:\n" << S2 << "\n";

  // The trans() adapter:
  cout << "Trans(A): " << trans(A).nrows() << "x" << trans(A).ncols() << "\n";
  cout << "After setting (2,0) to 1000.\n";
  trans(A)(2,0) = 1000;
  print_all_matrix(cout,trans(A));
  cout << "\n";
  // We have to use the (,) form for indexing here, since [][] does
  // not work with the trans() adapter!


  // Sparse Matrices
  //===================
  cout << "\n"
       << "----------------\n";
  cout << "Sparse Matrices:\n";
  cout << "----------------\n";

  SPARSE S(n,2*n);

  // How many nonzero elements?
  cout << "S.nnz():\n" << S.nnz() << "\n";

  S[1][2]  = 1;
  S[3][5]  = 2;
  
  cout << "S.nnz():\n" << S.nnz() << "\n";
  cout << "S:\n" << S << "\n";

  MATRIX At(2*n,n);
  copy(trans(A),At);
  cout << "At:\n" << At << "\n";

  // Calculate S*transpose(A):
  MATRIX Z(n,n);
  cout << "Capacity of Z:\n" << Z.capacity() << "\n";

  mult(S,At,Z);
  //  copy(trans(S),Z);

  cout << "Z = S*At:\n" << Z << "\n";


  // Demonstrate basic math functions:
  cout << "\n"
       << "---------------------\n";
  cout << "Basic math functions:\n";
  cout << "---------------------\n";

  x[0] = .1;
  x[4] = 5;
  cout << "x: " << x << "\n";
  
  // Non-return version of math functions:
  VECTOR y(x.size());
  transf(x,sqrt,y);
  cout << "Sqrt(x): " << y << "\n";

  // Return version of another math functions:
  cout << "Exp(x): " << transf(x,exp) << "\n";

  // Can I call the non-return version with the same variable as
  // return argument?
  copy(x,y);
  transf(y,log,y);
  cout << "Log(x): " << y << "\n";

  copy(x,y);
  ele_mult(y,y,y);
  cout << "x^2: " << y << "\n";

  // Demonstrate nlinspace (return version):
  cout << "Nlinspace(1,2,5): " << nlinspace(1,2,5) << "\n";

  // Linear interpolation:
  x = nlinspace(1,2,5);
  y = make_vector(10,20,30,40,50);
  VECTOR xi = nlinspace(1,2,10);
  VECTOR yi(10);
  interp_lin_vector(yi,x,y,xi);
  cout << "y interpolated: " << yi << "\n";

  resize(A,5,5);
  rand_matrix_uniform(A,1,2);
  cout << "Random matrix A(1-2):\n" << A << "\n";

  // Transform applied to Matrix:
  resize(A,2,2);
  A[0][0] = 1;
  A[0][1] = 2;
  A[1][0] = 3;
  A[1][1] = 4;

  cout << "A:\n"       << A              << "\n";
  cout << "Sqrt(A):\n" << transf(A,sqrt) << "\n";
  
  // Demonstrate, how std::bind1st and std::bind2nd can be used to use
  // transf even with functions taking another (scalar) argument:
  //  cout << "pow(A,3):\n" << transf(A,std::binder2nd(sqrt,3)) << "\n";
  // FIXME: Dunno how to make this work.


  // Demonstrate SYMMETRIC
  cout << "\nSymmetric Matrices:\n";
  SYMMETRIC As(3,3);

  As[0][0] = 1;
  As[0][1] = 2;
  cout << "As:\n"       << As              << "\n";


  resize(As,4,4);
  mtl::set(As,1);
  cout << "As:\n"       << As              << "\n";




  // Demonstrate Choleski factorization
  cout << "\nCholesky Factorization:\n";
  // Copying a symmetric positive definite matrix    
  SYMMETRIC C(4,4);
  
  C[0][0] = 3.33;
  C[0][1] = 1.33;
  C[0][2] = -0.67;
  C[0][3] = -2.67;
 
  C[1][1] = 0.67;
  C[1][2] = 0;
  C[1][3] = -0.67;

  C[2][2] = 0.67;
  C[2][3] = 1.33;

  C[3][3] = 3.33;
  cout << "Symmetric matrix:\n";
  cout << "C:\n"       << C              << "\n";

  // Getting the Choleski factor of C
  MATRIX R(4,4);
  chol(R,C);
  cout << "R:\n"       << R              << "\n";
  
  // Demonstrate BANDEDSYM
  cout << "\nBanded Symmetric Matrices:\n";
  SYMMETRIC Abs(4,1);

  // Try resize:
  resize(Abs,5,1);

  mtl::set(Abs,1); 

  // This does not seem to work either... :-(
//    {
//      SYMMETRIC::iterator ri = Abs.begin();
//      while (ri != Abs.end()) {
//        SYMMETRIC::Row::iterator i = (*ri).begin();
//        while (i != (*ri).end()) {
//  	*i = 1;
//  	++i;
//        }
//        ++ri;
//      }
//    }

  cout << "Abs (all elements set to 1):\n"       << Abs              << "\n";
  cout << "Only the diagonal +- 1 should be set.\n"
       << "Dunno why this does not work.\n";


}


#include "matpackI.h"
//#include "matpackII.h"

void fill_with_junk(VectorView x)
{
  x = 999;
}

void fill_with_junk(MatrixView x)
{
  x = 888;
}

int test1()
{
  Vector v(20);

  cout << "v.nelem() = " << v.nelem() << "\n";

  for (Index i=0; i<v.nelem(); ++i )
    v[i] = i;

  cout << "v.begin() = " << *v.begin() << "\n";

  cout << "v = \n" << v << "\n";

  fill_with_junk(v[Range(1,8,2)][Range(2,joker)]);
  //  fill_with_junk(v);

  Vector v2 = v[Range(2,4)];

  cout << "v2 = \n" << v2 << "\n";

  for (Index i=0 ; i<1000; ++i)
    {
      Vector v3(1000);
      v3 = i;
    }

  v2[Range(joker)] = 88;

  v2[Range(0,2)] = 77;

  cout << "v = \n" << v << "\n";
  cout << "v2 = \n" << v2 << "\n";
  cout << "v2.nelem() = \n" << v2.nelem() << "\n";

  Vector v3;
  v3.resize(v2.nelem());
  v3 = v2;

  cout << "\nv3 = \n" << v3 << "\n";
  fill_with_junk(v2);
  cout << "\nv3 after junking v2 = \n" << v3 << "\n";
  v3 *= 2;
  cout << "\nv3 after *2 = \n" << v3 << "\n";

  Matrix M(10,15);
  {
    Numeric n=0;
    for (Index i=0; i<M.nrows(); ++i)
      for (Index j=0; j<M.ncols(); ++j)
	M(i,j) = ++n;
  }

  cout << "\nM =\n" << M << "\n";

  cout << "\nM(Range(2,4),Range(2,4)) =\n" << M(Range(2,4),Range(2,4)) << "\n";

  cout << "\nM(Range(2,4),Range(2,4))(Range(1,2),Range(1,2)) =\n"
       << M(Range(2,4),Range(2,4))(Range(1,2),Range(1,2)) << "\n";

  cout << "\nM(1,Range(joker)) =\n" << M(1,Range(joker)) << "\n";

  cout << "\nFilling M(1,Range(1,2)) with junk.\n";
  fill_with_junk(M(1,Range(1,2)));
    
  cout << "\nM(Range(0,4),Range(0,4)) =\n" << M(Range(0,4),Range(0,4)) << "\n";

  cout << "\nFilling M(Range(4,2,2),Range(6,3)) with junk.\n";

  MatrixView s = M(Range(4,2,2),Range(6,3));
  fill_with_junk(s);

  cout << "\nM =\n" << M << "\n";

  const Matrix C = M;

  cout << "\nC(Range(3,4,2),Range(2,3,3)) =\n"
       << C(Range(3,4,2),Range(2,3,3)) << "\n";

  cout << "\nC(Range(3,4,2),Range(2,3,3)).transpose() =\n"
       << transpose(C(Range(3,4,2),Range(2,3,3))) << "\n";

  return 0;
}

void test2()
{
  Vector v(50000000);

  cout << "v.nelem() = " << v.nelem() << "\n";

  cout << "Filling\n";
//   for (Index i=0; i<v.nelem(); ++i )
//     v[i] = sqrt(i);
  v = 1.;
  cout << "Done\n";

}


// void test3()
// {
//   SparseMatrix M(10,15);

//   cout << "M.nrows(), M.ncols() = "
//        << M.nrows() << ", " << M.ncols() << "\n";

//   for (Index i=0; i<10; ++i)
//     M(i,i) = i+1;

//   cout << "\nM = \n" << M;

//   const SparseMatrix S(M);

//   cout << "\nS(2,0) = " << S(2,0) << "\n";

//   cout << "\nS = \n" << S;

// }

void test4()
{
  Vector a(10);
  Vector b(a.nelem());
  
  for ( Index i=0; i<a.nelem(); ++i )
    {
      a[i] = i+1;
      b[i] = a.nelem()-i;
    }

  cout << "a = \n" << a << "\n";
  cout << "b = \n" << b << "\n";
  cout << "a*b \n= " << a*b << "\n";

  Matrix A(11,6);
  Matrix B(10,20);
  Matrix C(20,5);

  B = 2;
  C = 3;
  mult(A(Range(1,joker),Range(1,joker)),B,C);

  //  cout << "\nB =\n" << B << "\n";
  //  cout << "\nC =\n" << C << "\n";
  cout << "\nB*C =\n" << A << "\n";
  
}

void test5()
{
  Vector a(10);
  Vector b(20);
  Matrix M(10,20);

  // Fill b and M with a constant number:
  b = 1;
  M = 2;

  cout << "b = \n" << b << "\n";
  cout << "M =\n" << M << "\n";

  mult(a,M,b);    // a = M*b
  cout << "\na = M*b = \n" << a << "\n";

  mult(transpose(b),transpose(a),M);    // b^t = a^t * M
  cout << "\nb^t = a^t * M = \n" <<  transpose(b) << "\n";
  
}

void test6()
{
  Index n = 5000;
  Vector x(1,n,1), y(n);
  Matrix M(n,n);
  M = 1;
  //  cout << "x = \n" << x << "\n";

  cout << "Transforming.\n";
  //  transform(x,sin,x);
  // transform(transpose(y),sin,transpose(x));
  //  cout << "sin(x) =\n" << y << "\n";
  for (Index i=0; i<1000; ++i)
    {
      //      mult(y,M,x);
      transform(y,sin,static_cast<MatrixView>(x));
      x+=1;
    }
  //  cout << "y =\n" << y << "\n";
  
  cout << "Done.\n";
}

void test7()
{
  Vector x(1,20000000,1);
  Vector y(x.nelem());
  transform(y,sin,x);
  cout << "min(sin(x)), max(sin(x)) = " << min(y) << ", " << max(y) << "\n";
}

int main()
{
  test7();
  return 0;
}

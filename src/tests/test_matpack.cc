#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <iostream>
#include "array.h"
#include "exceptions.h"
#include "math_funcs.h"
#include "matpack.h"
#include "matpack_mdspan_helpers_eigen.h"
#include "mystring.h"
#include "test_utils.h"
#include "time.h"
#include "wigner_functions.h"

using std::cout;
using std::runtime_error;

Numeric by_reference(const Numeric& x) { return x + 1; }

Numeric by_value(Numeric x) { return x + 1; }

void fill_with_junk(StridedVectorView x) { x = 999; }

void fill_with_junk(StridedMatrixView x) { x = 888; }

int test1() {
  Vector v(20);

  cout << "v.size() = " << v.size() << "\n";

  for (Size i = 0; i < v.size(); ++i) v[i] = (Numeric)i;

  cout << "v.begin() = " << *v.begin() << "\n";

  cout << "v = \n" << std::format("{}", v) << "\n";

  fill_with_junk(v[StridedRange(1, 8, 2)][Range(2, v.size()-2)]);
  //  fill_with_junk(v);

  Vector v2{v[Range(2, 4)]};

  cout << "v2 = \n" << std::format("{}", v2) << "\n";

  for (Index i = 0; i < 1000; ++i) {
    Vector v3(1000);
    v3 = (Numeric)i;
  }

  v2[joker] = 88;

  v2[Range(0, 2)] = 77;

  cout << "v = \n" << std::format("{}", v) << "\n";
  cout << "v2 = \n" << std::format("{}", v2) << "\n";
  cout << "v2.size() = \n" << v2.size() << "\n";

  Vector v3;
  v3.resize(v2.size());
  v3 = v2;

  cout << "\nv3 = \n" << std::format("{}", v3) << "\n";
  fill_with_junk((VectorView)v2);
  cout << "\nv3 after junking v2 = \n" << std::format("{}", v3) << "\n";
  v3 *= 2;
  cout << "\nv3 after *2 = \n" << std::format("{}", v3) << "\n";

  Matrix M(10, 15);
  {
    Numeric n = 0;
    for (Index i = 0; i < M.nrows(); ++i)
      for (Index j = 0; j < M.ncols(); ++j) M[i, j] = ++n;
  }

  cout << "\nM =\n" << std::format("{}", M) << "\n";

  cout << "\nM(Range(2,4),Range(2,4)) =\n"
       << std::format("{}", M[Range(2, 4), Range(2, 4)]) << "\n";

  cout << "\nM(Range(2,4),Range(2,4))(Range(1,2),Range(1,2)) =\n"
       << std::format("{}", M[Range(2, 4), Range(2, 4)][Range(1, 2), Range(1, 2)]) << "\n";

  cout << "\nM(1,joker) =\n" << std::format("{}", M[1, joker]) << "\n";

  cout << "\nFilling M(1,Range(1,2)) with junk.\n";
  fill_with_junk(M[1, Range(1, 2)]);

  cout << "\nM(Range(0,4),Range(0,4)) =\n"
       <<std::format("{}",  M[Range(0, 4), Range(0, 4)]) << "\n";

  cout << "\nFilling M(Range(4,2,2),Range(6,3)) with junk.\n";

  StridedMatrixView s = M[StridedRange(4, 2, 2), Range(6, 3)];
  fill_with_junk(s);

  cout << "\nM =\n" << std::format("{}", M) << "\n";

  const Matrix C = M;

  cout << "\nC(Range(3,4,2),Range(2,3,3)) =\n"
       << std::format("{}", C[StridedRange(3, 4, 2), StridedRange(2, 3, 3)]) << "\n";

  cout << "\nC(Range(3,4,2),Range(2,3,3)).transpose() =\n"
       << std::format("{}", transpose(C[StridedRange(3, 4, 2), StridedRange(2, 3, 3)])) << "\n";

  return 0;
}

void test2() {
  Vector v(50000000);

  cout << "v.size() = " << v.size() << "\n";

  cout << "Filling\n";
  //   for (Index i=0; i<v.size(); ++i )
  //     v[i] = sqrt(i);
  v = 1.;
  cout << "Done\n";
}

void test5() {
  Vector a(10);
  Vector b(20);
  Matrix M(10, 20);

  // Fill b and M with a constant number:
  b = 1;
  M = 2;

  cout << "b = \n" << std::format("{}", b) << "\n";
  cout << "M =\n" << std::format("{}", M) << "\n";

  mult(a, M, b);  // a = M*b
  cout << "\na = M*b = \n" << std::format("{}", a) << "\n";

  mult(transpose(b.view_as(1, b.size())), transpose(a.view_as(1, a.size())), M);  // b^t = a^t * M
  cout << "\nb^t = a^t * M = \n" << std::format("{}", transpose(b.view_as(1, b.size()))) << "\n";
}

void test6() {
  Index n = 5000;
  Vector x=matpack::uniform_grid(1, n, 1), y(n);
  Matrix M(n, n);
  M = 1;
  //  cout << "x = \n" << x << "\n";

  cout << "Transforming.\n";
  //  transform(x,sin,x);
  // transform(transpose(y),sin,transpose(x));
  //  cout << "sin(x) =\n" << y << "\n";
  for (Index i = 0; i < 1000; ++i) {
    //      mult(y,M,x);
 y.unary_transform(x, [](auto z){return std::sin(z);});
    x += 1;
  }
  //  cout << "y =\n" << y << "\n";

  cout << "Done.\n";
}

void test7() {
  Vector x=matpack::uniform_grid(1, 20000000, 1);
  Vector y(x.size());
 y.unary_transform(x, [](auto z){return std::sin(z);});
  cout << "min(sin(x)), max(sin(x)) = " << min(y) << ", " << max(y) << "\n";
}

void test8() {
  Vector x(80000000);
  for (Size i = 0; i < x.size(); ++i) x[i] = (Numeric)i;
  cout << "Done."
       << "\n";
}

void test9() {
  // Initialization of Matrix with view of other Matrix:
  Matrix A(4, 8);
  Matrix B(A[joker, Range(0, 3)]);
  cout << "B = " << std::format("{}", B) << "\n";
}

void test10() {
  // Initialization of Matrix with a vector (giving a 1 column Matrix).

  // At the moment doing this with a non-const Vector will result in a
  // warning message.
  Vector v=matpack::uniform_grid(1, 8, 1);
  Matrix M(v.view_as(1, v.size()));
  cout << "M = " << std::format("{}", M) << "\n";
}

void test11() {
  // Assignment between Vector and Matrix:

  // At the moment doing this with a non-const Vector will result in a
  // warning message.
  Vector v=matpack::uniform_grid(1, 8, 1);
  Matrix M(v.size(), 1);
  M = v.view_as(1, v.size());
  cout << "M = " << std::format("{}", M) << "\n";
}

void test12() {
  // Copying of Arrays

  Array<String> sa(3);
  sa[0] = "It's ";
  sa[1] = "a ";
  sa[2] = "test.";

  Array<String> sb(sa), sc(sa.size());

  cout << "sb = \n" << sb << "\n";

  sc = sa;

  cout << "sc = \n" << sc << "\n";
}

void test13() {
  // Mix vector and one-column matrix in += operator.
  const Vector v=matpack::uniform_grid(1, 8, 1);  // The const is necessary here to
                            // avoid compiler warnings about
                            // different conversion paths.
  Matrix M(v.view_as(1, v.size()));
  M += v.view_as(1, v.size());
  cout << "M = \n" << std::format("{}", M) << "\n";
}

void test14() {
  // Test explicit Array constructors.
  Array<String> a{"Test"};
  Array<Index> b{1, 2};
  Array<Numeric> c{1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0};
  cout << "a = \n" << a << "\n";
  cout << "b = \n" << b << "\n";
  cout << "c = \n" << c << "\n";
}

void test15() {
  // Test String.
  String a = "Nur ein Test.";
  cout << "a = " << a << "\n";
  String b(a, 5, -1);
  cout << "b = " << b << "\n";
}

/*void test16()
{
  // Test interaction between Array<Numeric> and Vector.
  Vector a;
  Array<Numeric> b;
  b.push_back(1);
  b.push_back(2);
  b.push_back(3);
  a.resize(b.size());
  a = b;
  cout << "b =\n" << b << "\n";
  cout << "a =\n" << a << "\n";
}*/

void test17() {
  // Test Sum.
  Vector a=matpack::uniform_grid(1, 10, 1);
  cout << "a.sum() = " << sum(a) << "\n";  // FIXME: sum() is no longer a member function
}

void test18() {
  // Test elementvise square of a vector:
  Vector a=matpack::uniform_grid(1, 10, 1);
  a *= a;
  cout << "a *= a =\n" << std::format("{}", a) << "\n";
}

void test19() {
  // There exists no explicit filling constructor of the form
  // Vector a(3,1.7).
  // But you can use the more general filling constructor with 3 arguments.

  Vector a=matpack::uniform_grid(1, 10, 1);
  Vector b=matpack::uniform_grid(5.3, 10, 0);
  cout << "a =\n" << std::format("{}", a) << "\n";
  cout << "b =\n" << std::format("{}", b) << "\n";
}

void test20() {
  // Test initialization list constructor:
  Vector a{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  cout << "a =\n" << std::format("{}", a) << "\n";
}

void test21() {
  Numeric s = 0;
  // Test speed of call by reference:
  cout << "By reference:\n";
  for (Index i = 0; i < (Index)1e8; ++i) {
    s += by_reference(s);
    s -= by_reference(s);
  }
  cout << "s = " << s << "\n";
}

void test22() {
  Numeric s = 0;
  // Test speed of call by value:
  cout << "By value:\n";
  for (Index i = 0; i < (Index)1e8; ++i) {
    s += by_value(s);
    s -= by_value(s);
  }
  cout << "s = " << s << "\n";
}

void test23() {
  // Test constructors that fills with constant:
  Vector a(10, 3.5);
  cout << "a =\n" << std::format("{}", a) << "\n";
  Matrix b(10, 10, 4.5);
  cout << "b =\n" << std::format("{}", b) << "\n";
}

void test24() {
  // Try element-vise multiplication of Matrix and Vector:
  Matrix a(5, 1, 2.5);
  Vector b=matpack::uniform_grid(1, 5, 1);
  a *= Matrix{b.view_as(1, b.size())};
  cout << "a*=b =\n" << std::format("{}", a) << "\n";
  a /= MatrixView{b.view_as(1, b.size())};
  cout << "a/=b =\n" << std::format("{}", a) << "\n";
  a += MatrixView{b.view_as(1, b.size())};
  cout << "a+=b =\n" << std::format("{}", a) << "\n";
  a -= ConstMatrixView{b.view_as(1, b.size())};
  cout << "a-=b =\n" << std::format("{}", a) << "\n";
}

void test25() {
  // Test min and max for Array:
  Array<Index> a{1, 2, 3, 4, 5, 6, 5, 4, 3, 2, 1};
  cout << "min/max of a = " << min(a) << "/" << max(a) << "\n";
}

void test26() {
  cout << "Test filling constructor for Array:\n";
  Array<String> a(4, "Hello");
  cout << "a =\n" << a << "\n";
}

void test27() {
  cout << "Test Arrays of Vectors:\n";
  Array<Vector> a;
  a.push_back({1.0, 2.0});
  a.push_back(matpack::uniform_grid(1.0, 10, 1.0));
  cout << "a =\n" << std::format("{}", a) << "\n";
}

void test28() {
  cout << "Test default constructor for Matrix:\n";
  Matrix a;
  Matrix b(a);
  cout << "b =\n" << std::format("{}", b) << "\n";
}

void test29() {
  cout << "Test Arrays of Matrix:\n";
  ArrayOfMatrix a;
  Matrix b;

  b.resize(2, 2);
  b[0, 0] = 1;
  b[0, 1] = 2;
  b[1, 0] = 3;
  b[1, 1] = 4;
  a.push_back(b);
  b *= 2;
  a.push_back(b);

  a[0].resize(2, 3);
  a[0] = 4;

  a.resize(3);
  a[2].resize(4, 5);
  a[2] = 5;

  cout << "a =\n" << std::format("{}", a) << "\n";
}

void test30() {
  cout << "Test Matrices of size 0:\n";
  Matrix a(0, 0);
  //  cout << "a(0,0) =\n" << a(0,0) << "\n";
  a.resize(2, 2);
  a = 1;
  cout << "a =\n" << std::format("{}", a) << "\n";

  Matrix b(3, 0);
  //  cout << "b(0,0) =\n" << b(0,0) << "\n";
  b.resize(b.nrows(), b.ncols() + 3);
  b = 2;
  cout << "b =\n" << std::format("{}", b) << "\n";

  Matrix c(0, 3);
  //  cout << "c(0,0) =\n" << c(0,0) << "\n";
  c.resize(c.nrows() + 3, c.ncols());
  c = 3;
  cout << "c =\n" << std::format("{}", c) << "\n";
}

void test31() {
  cout << "Test Tensor3:\n";

  Tensor3 a(2, 3, 4, 1.0);

  Index fill = 0;

  // Fill with some numbers
  for (Index i = 0; i < a.npages(); ++i)
    for (Index j = 0; j < a.nrows(); ++j)
      for (Index k = 0; k < a.ncols(); ++k) a[i, j, k] = (Numeric)(++fill);

  cout << "a =\n" << std::format("{}", a) << "\n";

  cout << "Taking out first row of first page:\n"
       << std::format("{}", a[0, 0, joker]) << "\n";

  cout << "Taking out last column of second page:\n"
       << std::format("{}", a[1, joker, a.ncols() - 1]) << "\n";

  cout << "Taking out the first letter on every page:\n"
       << std::format("{}", a[joker, 0, 0]) << "\n";

  cout << "Taking out first page:\n"
       << std::format("{}", a[0, joker, joker]) << "\n";

  cout << "Taking out last row of all pages:\n"
       << std::format("{}", a[joker, a.nrows() - 1, joker]) << "\n";

  cout << "Taking out second column of all pages:\n"
       << std::format("{}", a[joker, joker, 1]) << "\n";

  a *= 2;

  cout << "After element-vise multiplication with 2:\n" << std::format("{}", a) << "\n";

  a.unary_transform(a, [](auto z){return std::sqrt(z);});

  cout << "After taking the square-root:\n" << std::format("{}", a) << "\n";

  Index s = 200;
  cout << "Let's allocate a large tensor, "
       << (Numeric)(s * s * s * 8) / 1024. / 1024. << " MB...\n";

  a.resize(s, s, s);

  cout << "Set it to 1...\n";

  a = 1;

  cout << "a(900,900,900) = " << a[90, 90, 90] << "\n";

  fill = 0;

  cout << "Fill with running numbers, using for loops...\n";
  for (Index i = 0; i < a.npages(); ++i)
    for (Index j = 0; j < a.nrows(); ++j)
      for (Index k = 0; k < a.ncols(); ++k) a[i, j, k] = (Numeric)(++fill);

  cout << "Max(a) = ...\n";

  cout << max(a) << "\n";
}

void test32() {
  cout << "Test von X = A*X:\n";
  Matrix X(3, 3), A(3, 3), B(3, 3);

  for (Index j = 0; j < A.nrows(); ++j)
    for (Index k = 0; k < A.ncols(); ++k) {
      X[j, k] = 1;
      A[j, k] = (Numeric)(j + k);
    }
  cout << "A:\n" << std::format("{}", A) << "\n";
  cout << "X:\n" << std::format("{}", X) << "\n";

  mult(B, A, X);
  cout << "B = A*X:\n" << std::format("{}", B) << "\n";
  mult(X, A, X);
  cout << "X = A*X:\n" << std::format("{}", X) << "\n";

  cout << "This is not the same, and should not be, because you\n"
       << "are not allowed to use the same matrix as input and output\n"
       << "for mult!\n";
}

void test33() {
  cout << "Making things look bigger than they are...\n";

  {
    cout << "1. Make a scalar look like a vector:\n";
    Numeric a = 3.1415;  // Just any number here.
    VectorView av(a);
    cout << "a, viewed as a vector: " << std::format("{}", av) << "\n";
    cout << "Describe a: " << "DEPR" << "\n";  // FIXME: Removed describe for Numeric
    av[0] += 1;
    cout << "a, after the first element\n"
         << "of the vector has been increased by 1: " << a << "\n";
  }

  {
    cout
        << "\n2. Make a vector look like a matrix:\n"
        << "This is an exception, because the new dimension is added at the end.\n";
    Vector a{1, 2, 3, 4, 5};
    MatrixView am{a.view_as(1, a.size())};
    cout << "a, viewed as a matrix:\n" << std::format("{}", am) << "\n";
    cout << "Trasnpose view:\n" << std::format("{}", transpose(am)) << "\n";

    cout << "\n3. Make a matrix look like a tensor3:\n";
    Tensor3View at3 = Tensor3View{am};
    cout << "at3 = \n" << std::format("{}", at3) << "\n";
    at3[0, 2, 0] += 1;
    cout << "a after Increasing element at3(0,2,0) by 1: \n" << std::format("{}", a) << "\n\n";

    Tensor4View at4 = Tensor4View{at3};
    cout << "at4 = \n" << std::format("{}", at4) << "\n";

    Tensor5View at5 = Tensor5View{at4};
    cout << "at5 = \n" << std::format("{}", at5) << "\n";

    Tensor6View at6 = Tensor6View{at5};
    cout << "at6 = \n" << std::format("{}", at6) << "\n";

    Tensor7View at7 = Tensor7View{at6};
    cout << "at7 = \n" << std::format("{}", at7) << "\n";

    at7[0, 0, 0, 0, 0, 2, 0] -= 1;

    cout << "After subtracting 1 from at7(0,0,0,0,0,2,0)\n"
         << "a = " << std::format("{}", a) << "\n";

    cout << "\nAll in one go:\n";
    Numeric b = 3.1415;  // Just any number here.
    Tensor7View bt7 = Tensor7View{Tensor6View(
        Tensor5View(Tensor4View(Tensor3View(MatrixView(VectorView(b))))))};
    cout << "bt7:\n" << std::format("{}", bt7) << "\n";
  }
}

void test35() {
  cout << "Test the new copy semantics.\n";
  Vector a=matpack::uniform_grid(1, 4, 1);
  Vector b;

  b = a;
  cout << "b = " << std::format("{}", b) << "\n";

  Vector aa=matpack::uniform_grid(1, 5, 1);
  ConstVectorView c = aa;
  b = c;
  cout << "b = " << std::format("{}", b) << "\n";

  Vector aaa=matpack::uniform_grid(1, 6, 1);
  VectorView d = aaa;
  b = d;
  cout << "b = " << std::format("{}", b) << "\n";
}

void test36() {
  cout << "Test using naked joker on Vector.\n";
  Vector a=matpack::uniform_grid(1, 4, 1);
  VectorView b = a[joker];
  cout << "a = " << std::format("{}", a) << "\n";
  cout << "b = " << std::format("{}", b) << "\n";
}

void test37(const Index& i) {
  Vector v1=matpack::uniform_grid(5e-15, 10, 0.42e-15 / 11);
  Vector v2 = v1;
  //  Index i = 10000000;

  v1 /= (Numeric)i;
  v2 /= (Numeric)i;
  cout.precision(12);
  //  cout.setf(ios_base::scientific,ios_base::floatfield);
  v1 *= v1;
  v2 *= v2;
  cout << std::format("{}", v1) << '\n';
  cout << std::format("{}", v2) << '\n';
}

void test38() {
  Vector v(5, 0.);
  Numeric* const a = v.data_handle();

  a[4] = 5.;

  cout << std::format("{}", v) << '\n';
  cout << '\n' << "========================" << '\n' << '\n';

  Matrix m(5, 5, 0.);
  Numeric* const b = m.data_handle();

  b[4] = 5.;

  cout << std::format("{}", m) << '\n';
  cout << '\n' << "========================" << '\n' << '\n';

  Tensor3 t3(5, 6, 7, 0.);
  Numeric* const c = t3.data_handle();

  c[6] = 5.;

  cout << std::format("{}", t3) << '\n';
}

void test39() {
  Vector v1=matpack::uniform_grid(1, 5, 1), v2(5);

  v2 = v1 * 2;
  // Unfortunately, this thing compiles, but at least it gives an
  // assertion failure at runtime. I think what happens is that it
  // implicitly creates a one element vector out of the "2", then
  // tries to do a scalar product with v1.

  cout << std::format("{}", v2) << '\n';
}

void test40() {
  Vector v(100);

  v = 5;

  cout << std::format("{}", v) << '\n';
}

void test41() {
  const Vector v1(10, 1);

  ConstVectorView vv = v1[Range(0, 5)];

  cout << "Vector:     " << std::format("{}", v1) << '\n';
  cout << "VectorView: " << std::format("{}", vv) << '\n';

  try {
    //vv = 2;  FIXME: No more runtime errors for simple types
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << '\n';
    exit(EXIT_FAILURE);
  }

  cout << "VectorView: " << std::format("{}", vv) << '\n';
  cout << "Vector:     " << std::format("{}", v1) << '\n';
}

// Test behaviour of VectorView::operator MatrixView, for which I have fixed
// a bug today.
// SAB 2013-01-18
void test42() {
  cout << "test42\n";
  Vector x{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  cout << "x: " << std::format("{}", x) << '\n';

  StridedVectorView y = x[StridedRange(2, 4, 2)];
  cout << "y: " << std::format("{}", y) << '\n';
}


// Test function for internal use, must return 2 as last n2r
constexpr Rational test_numeric2rational(const Index i, const Index maxi, const Rational r=0, const Rational n2r=0) {
  if (i > maxi)
    return n2r;
  else {
    return 
    (r == n2r) ?
    test_numeric2rational(i+1, maxi, Rational(maxi + i, maxi), numeric2rational(1 + Numeric(i)/Numeric(maxi), 12)) :
    throw std::logic_error("Fail to initialize");;
  }
}


void test43() {
  // Simple construction compile-time test
  constexpr Rational r=3_2;  // should be 3/2
  static_assert(r.numer == 3, "Rational fail to initialize properly");
  static_assert(r.denom == 2, "Rational fail to initialize properly");

  // Simple expression compile-time test
  constexpr Rational r2 = 1 + r;  // should be 5/2
  static_assert(r2.numer == 5, "Rational fail to initialize properly");
  static_assert(r2.denom == 2, "Rational fail to initialize properly");

  // Div-zero by index generates an undefined rational:  0/0
  constexpr Rational r3 = r / 0;  // should be undefined
  static_assert(r3.numer == 0, "Rational fail to initialize properly");
  static_assert(r3.denom == 0, "Rational fail to initialize properly");
  static_assert(r3 not_eq r3, "Undefined rational does not equal self");
  static_assert(r3.isUndefined(), "Rational fail to initialize properly");
  
  // Mul-zero gives 0/1
  constexpr Rational r4 = r * 0;  // should be undefined
  static_assert(r4.numer == 0, "Rational fail to initialize properly");
  static_assert(r4.denom == 1, "Rational fail to initialize properly");

  // Complicated expression compile-time test
  constexpr Rational r5 = (Rational(1, 9) * r) % 3;  // should be 1/6
  static_assert(r5.numer == 1, "Rational fail to initialize properly");
  static_assert(r5.denom == 6, "Rational fail to initialize properly");

  // The simplify operation simplifies expressions
  constexpr Rational r6 = 10 % r;  // should be 1
  static_assert(r6.numer == 1, "Rational fail to initialize properly");
  static_assert(r6.denom == 1, "Rational fail to initialize properly");
  static_assert(r6.toInt() == 1, "Rational fail to initialize properly");
  static_assert(r6.toIndex() == 1, "Rational fail to initialize properly");
  static_assert(r6.toNumeric() == 1e0, "Rational fail to initialize properly");
  
  constexpr Index rattest=1<<8;
  constexpr Rational r7 = test_numeric2rational(0, rattest);
  static_assert(r7 == 2, "Rational fail to initialize properly");
}

void test46() {
  Vector v(5, 0.);
  nlinspace(v, 1, 10, 10);
  VectorView v1 = v;
  auto compare_func = [](Numeric n) { return n != 0; };
  cout << std::any_of(v1.begin(), v1.end(), compare_func) << '\n';
  v1 = 0.;
  cout << std::any_of(v1.begin(), v1.end(), compare_func) << '\n';
}

//! Perform matrix multiplication test.
/*!
  Creates to random m times k matrix A and a random k times n matrix
  B and computes C = A * B, C = A^T * B, C = A * B^T, C = A^T * B^T unsing
  the mult function from matpackI.cc. The result is verified using the
  mult_general functions. In addition to the operations above, ntests numbers
  of random submatrices of A and B are generated and their multiplication is
  also tested. Returns the maximum element-wise, relative error that occurred in
  the tests.

  \param[in] k Dimension over which the multiplication is performed.
  \param[in] m Number of rows of A.
  \param[in] n Number of coluns of B.
  \param[in] ntests Number of standard multiplications to be performed.
  \param[in] nsubtests Number of submatrix multiplications to be performed.
  \param[in] bool If verbose  is true the results of each test are printed to
  standard out.

  \return The maximum element-wise, relative error that occured in the tests.
*/
Numeric matrix_mult(
    Index k, Index m, Index n, Index ntests, Index nsubtests, bool verbose) {
  Numeric max_err = 0;
  Matrix A(m, k), B(k, n), C(m, n), C_ref(m, n);

  for (Index i = 0; i < ntests; i++) {
    random_fill_matrix(A, 1000, false);
    random_fill_matrix(B, 1000, false);

    if (verbose) {
      cout.precision(15);
      cout << "MATRIX MULTIPLICATION: m = " << m << ", k = " << k;
      cout << ", n = " << n << '\n';
    }

    // A * B

    random_fill_matrix(C, 10, false);
    random_fill_matrix(C_ref, 10, false);

    mult(C, A, B);
    mult(C_ref, A, B);

    Numeric err_mul = get_maximum_error(C, C_ref, true);
    if (err_mul > max_err) max_err = err_mul;

    if (verbose) {
      cout << "\t"
           << "A * B: max. rel. error = " << err_mul << '\n';
    }

    // A^T * B

    Matrix AT(k, m);
    random_fill_matrix(AT, 100.0, false);

    mult(C, transpose(AT), B);
    mult(C_ref, transpose(AT), B);
    Numeric err_trans1 = get_maximum_error(C, C_ref, true);
    if (err_trans1 > max_err) max_err = err_trans1;

    if (verbose) {
      cout << "\t"
           << "A^T * B: max. rel. err = " << err_trans1 << '\n';
    }

    // A * B^T

    Matrix BT(n, k);
    random_fill_matrix(BT, 100.0, false);

    mult(C, A, transpose(BT));
    mult(C_ref, A, transpose(BT));
    Numeric err_trans2 = get_maximum_error(C, C_ref, true);
    if (err_trans2 > max_err) max_err = err_trans2;

    if (verbose) {
      cout << "\t"
           << "A * B^T: max. rel. err = " << err_trans2 << '\n';
    }

    // A^T * B^T

    mult(C, transpose(AT), transpose(BT));
    mult(C_ref, transpose(AT), transpose(BT));
    Numeric err_trans3 = get_maximum_error(C, C_ref, true);
    if (err_trans3 > max_err) max_err = err_trans3;

    if (verbose) {
      cout << "\t"
           << "A^T * B^T: max. rel. err = " << err_trans3 << '\n';
      cout << '\n';
    }
  }

  // Multiplication of submatrices.
  Rand<Index> k_rand(1, k - 1), m_rand(1, m - 1), n_rand(1, n - 1);

  for (Index i = 0; i < nsubtests - 1; i++) {
    Index k1(k_rand()), m1(m_rand()), n1(n_rand());

    Range r1(random_range(m1));
    Range r2(random_range(k1));
    Range r3(random_range(n1));

    StridedMatrixView C_sub(C[r1, r3]);
    StridedMatrixView C_sub_ref(C_ref[r1, r3]);

    StridedConstMatrixView A_sub(A[r1, r2]);
    StridedConstMatrixView B_sub(B[r2, r3]);

    mult(C_sub, A_sub, B_sub);
    mult(C_sub_ref, A_sub, B_sub);

    Numeric err = get_maximum_error(C_sub, C_sub_ref, true);
    if (err > max_err) max_err = err;

    if (verbose) {
      cout << "\t"
           << "Submatrix multiplication: max. rel. err = " << err;
      cout << '\n';
    }
  }

  if (verbose) {
    cout << '\n';
  }

  return max_err;
}

//! Perform matrix multiplication tests.
/*!
  Use matrix_mult function to perform matrix multiplication tests,
  for different combinations of dimensions including:
      k = m = n = 0
      k = m = n = 1
      k = 10, m = n = 1
      k = 10, m = 100, n = 10
      k = 10, m = 100, n = 100
      k = m = n = 100

 */
Numeric test_matrix_multiplication(bool verbose) {
  Numeric max_err, err;

  // Trivial case k = m = n = 0.
  max_err = matrix_mult(0, 0, 0, 10, 0, verbose);

  // k = 1, m = 1, n = 1.
  err = matrix_mult(1, 1, 1, 10, 0, verbose);
  if (err > max_err) max_err = err;

  // k = 10, m = 1, n = 10.
  err = matrix_mult(10, 1, 10, 20, 20, verbose);
  if (err > max_err) max_err = err;

  // k = 10, m = 1, n = 10.
  err = matrix_mult(10, 1, 10, 20, 20, verbose);
  if (err > max_err) max_err = err;

  // k = 10, m = 1, n = 1.
  err = matrix_mult(10, 1, 1, 20, 20, verbose);
  if (err > max_err) max_err = err;

  // k = 100, m = 100, n = 100.
  err = matrix_mult(200, 200, 200, 20, 20, verbose);
  if (err > max_err) max_err = err;

  // k = 10, m = 100, n = 10.
  err = matrix_mult(10, 100, 10, 20, 20, verbose);
  if (err > max_err) max_err = err;

  // k = 10, m = 100, n = 100.
  err = matrix_mult(10, 100, 100, 20, 20, verbose);
  if (err > max_err) max_err = err;

  // k = 100, m = 100, n = 100, 100 submatrix multiplications.
  err = matrix_mult(100, 100, 100, 0, 100, verbose);
  if (err > max_err) max_err = err;

  return max_err;
}

//! Check if the empty function is working correctly
void test_empty() {
  {
    cout << "Array" << '\n';
    ArrayOfIndex v;
    v.resize(0);
    cout << v.empty() << '\n';
    v.resize(1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Vector" << '\n';
    Vector v;
    v.resize(0);
    cout << v.empty() << '\n';
    v.resize(1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Matrix" << '\n';
    Matrix v;
    v.resize(0, 1);
    cout << v.empty() << '\n';
    v.resize(1, 0);
    cout << v.empty() << '\n';
    v.resize(1, 1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Tensor3" << '\n';
    Tensor3 v;
    v.resize(0, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 0, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 0);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Tensor4" << '\n';
    Tensor4 v;
    v.resize(0, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 0, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 0, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 0);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Tensor5" << '\n';
    Tensor5 v;
    v.resize(0, 1, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 0, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 0, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 0, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 0);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Tensor6" << '\n';
    Tensor6 v;
    v.resize(0, 1, 1, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 0, 1, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 0, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 0, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 0, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 1, 0);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 1, 1);
    cout << v.empty() << '\n';
  }
  {
    cout << "Tensor7" << '\n';
    Tensor7 v;
    v.resize(0, 1, 1, 1, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 0, 1, 1, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 0, 1, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 0, 1, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 0, 1, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 1, 0, 1);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 1, 1, 0);
    cout << v.empty() << '\n';
    v.resize(1, 1, 1, 1, 1, 1, 1);
    cout << v.empty() << '\n';
  }
}

void test47() {
  // Selecting empty matpack dimensions with Joker shouldn't fail
  Vector v;
  Matrix m;

  VectorView vv = v[joker];
  MatrixView mv = m[joker, joker];

  std::cout << "vv.size: " << vv.size() << '\n';
  std::cout << "mv.nrows: " << mv.nrows() << " ";
  std::cout << "mv.ncols: " << mv.ncols() << '\n';
}

void test_wigner_error() {
  try {
    wigner3j(1, 0, 1, 0, 0, 0);
  } catch(std::exception& e) {
    std::cerr << e.what() << '\n';
  }

  make_wigner_ready(250, 20000000, 3);
  wigner3j(1, 0, 1, 0, 0, 0);
}

void test_pow_negative_one() {
  std::vector<Index> x(30);
  std::iota(x.begin(), x.end(), -15);
  for (auto& i: x) std::cout << "-1^" << i << '=' << pow_negative_one(i) << '\n';
}

Matrix build_test_matrix(Index rows, Index cols) {
  Matrix a(rows, cols);
  for (Index i = 0; i < rows; i++) {
    for (Index j = 0; j < cols; j++) {
      a[i, j] = static_cast<Numeric>(10 * i + j + 1);
    }
  }
  return a;
}

void test_mult() {
  std::cout << "# MULT TEST ######################################\n";
  {
    std::cout << "TEST 1 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A (3, 3, 1);
    mult(y, A, x);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "x (in):\n" << std::format("{}", x) << '\n'
              << "y (out):\n" << std::format("{}", y) << '\n';
  }
  {
    std::cout << "TEST 2 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A = build_test_matrix(3, 3);
    mult(y, A, x);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "x (in):\n" << std::format("{}", x) << '\n'
              << "y (out):\n" << std::format("{}", y) << '\n';
  }
  {
    std::cout << "TEST 3 (y = A * x; NOT STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2});
    Matrix A = build_test_matrix(3, 2);
    mult(y, A, x);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "x (in):\n" << std::format("{}", x) << '\n'
              << "y (out):\n" << std::format("{}", y) << '\n';
  }
    {
    std::cout << "TEST 4 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A (3, 3, 1);
    mult(y, A, x);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "x (in):\n" << std::format("{}", x) << '\n'
              << "y (out):\n" << std::format("{}", y) << '\n';
  }
  {
    std::cout << "TEST 5 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2, 3});
    Matrix A = build_test_matrix(3, 3);
    mult(y, A, x);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "x (in):\n" << std::format("{}", x) << '\n'
              << "y (out):\n" << std::format("{}", y) << '\n';
  }
  {
    std::cout << "TEST 6 (y = A * x; STRIDED):\n";
    Vector y(3);
    Vector x(std::vector<Numeric>{1, 2});
    Matrix A = build_test_matrix(3, 2);
    mult(y, A, x);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "x (in):\n" << std::format("{}", x) << '\n'
              << "y (out):\n" << std::format("{}", y) << '\n';
  }
  {
    std::cout << "TEST 7 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 3);
    Matrix B(3, 3, 1);
    Matrix A (3, 3, 1);
    mult(C, A, B);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "B (in):\n" << std::format("{}", B) << '\n'
              << "C (out):\n" << std::format("{}", C) << '\n';
  }
  {
    std::cout << "TEST 8 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 4);
    Matrix A (3, 2, 1);
    Matrix B(2, 4, 1);
    mult(C, A, B);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "B (in):\n" << std::format("{}", B) << '\n'
              << "C (out):\n" << std::format("{}", C) << '\n';
  }
  {
    std::cout << "TEST 9 (C = A * B; NON STRIDED):\n";
    Matrix C(3, 4);
    Matrix A = build_test_matrix(3, 2);
    Matrix B = build_test_matrix(2, 4);
    mult(C, A, B);
    std::cout << "A (in):\n" << std::format("{}", A) << '\n'
              << "B (in):\n" << std::format("{}", B) << '\n'
              << "C (out):\n" << std::format("{}", C) << '\n';
  }
  std::cout << "#/MULT TEST ######################################\n";
}

int main() {
  //   test1();
  //   test2();
  //   test3();
  //   test4();
  //   test5();
  //   test6();
  //   test7();
  //   test8();
  //   test9();
  //   test10();
  //   test11();
  //   test12();
  //   test13();
  //   test14();
  //   test15();
  //   test16();
  //   test17();
  //   test18();
  //   test19();
  //   test20();
  //   test21();
  //   test22();
  //   test23();
  //   test24();
  //   test25();
  //   test26();
  //   test27();
  //   test28();
  //   test29();
  //   test30();
  //   test31();
  //   test32();
  //   test33();
  //   test34();
  //   test35();
  //   test36();
  //  Index i = 10000000;
  //  test37(i);
  //  test38();
  //  test39();
  //  test40();
  //  test41();
  //    test42();
  //   test43();
  //    test44();
  //    test45();
  //    test46();
  //  test47();
  //test48();
  //test_wigner_error();
  //test_pow_negative_one();
  //test_concepts();
  //test_mult();

  //    const double tolerance = 1e-9;
  //    double error;
  //
  //    // Matrix Vector Multiplication.
  //    error = test_matrix_vector_multiplication(false);
  //    cout << "Matrix Vector Multiplication: ";
  //    if (error > tolerance)
  //        cout << "FAILED, maximum error: " << error << '\n';
  //    else
  //        cout << "PASSED." << '\n';
  //
  //    // Matrix Matrix Multiplication.
  //    error = test_matrix_multiplication(false);
  //
  //    cout << "Matrix Matrix Multiplication: ";
  //    if (error > tolerance)
  //        cout << "FAILED, maximum error: " << error << '\n';
  //    else
  //        cout << "PASSED." << '\n';

  //test_diagonal( 100 );
  //test_empty();

  return 1;
}

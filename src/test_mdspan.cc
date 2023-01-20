#include <algorithm>
#include <complex>
#include <iostream>
#include <numeric>
#include <vector>

#include "debug.h"
#include "logic.h"
#include "matpack_concepts.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "matpack_iter.h"
#include "matpack_math.h"
#include "matpack_view.h"
#include "matpack_complex.h"

using std::cout;
using std::endl;

Numeric by_reference(const Numeric &x) { return x + 1; }

Numeric by_value(Numeric x) { return x + 1; }

void fill_with_junk(VectorView x) { x = 999; }

void fill_with_junk(MatrixView x) { x = 888; }

int test1() {
  Vector v(20);

  cout << "v.nelem() = " << v.nelem() << "\n";

  for (Index i = 0; i < v.nelem(); ++i)
    v[i] = (Numeric)i;

  cout << "v.begin() = " << *v.begin() << "\n";

  cout << "v = \n" << v << "\n";
  
  fill_with_junk(v[Range(1, 8, 2)][Range(2, joker)]);
  //  fill_with_junk(v);

  Vector v2 = v[Range(2, 4)];
  
  cout << "v2 = \n" << v2 << "\n";

  for (Index i = 0; i < 1000; ++i) {
    Vector v3(1000);
    v3 = (Numeric)i;
  }

  v2[Range(joker)] = 88;

  v2[Range(0, 2)] = 77;

  cout << "v = \n" << v << "\n";
  cout << "v2 = \n" << v2 << "\n";
  cout << "v2.nelem() = \n" << v2.nelem() << "\n";

  Vector v3;
  v3.resize(v2.nelem());
  v3 = v2;

  cout << "\nv3 = \n" << v3 << "\n";
  fill_with_junk((VectorView)v2);
  cout << "\nv3 after junking v2 = \n" << v3 << "\n";
  v3 *= 2;
  cout << "\nv3 after *2 = \n" << v3 << "\n";

  Matrix M(10, 15);
  {
    Numeric n = 0;
    for (Index i = 0; i < M.nrows(); ++i)
      for (Index j = 0; j < M.ncols(); ++j)
        M(i, j) = ++n;
  }

  cout << "\nM =\n" << M << "\n";

  cout << "\nM(Range(2,4),Range(2,4)) =\n"
       << M(Range(2, 4), Range(2, 4)) << "\n";

  cout << "\nM(Range(2,4),Range(2,4))(Range(1,2),Range(1,2)) =\n"
       << M(Range(2, 4), Range(2, 4))(Range(1, 2), Range(1, 2)) << "\n";

  cout << "\nM(1,Range(joker)) =\n" << M(1, Range(joker)) << "\n";

  cout << "\nFilling M(1,Range(1,2)) with junk.\n";
  fill_with_junk(M(1, Range(1, 2)));

  cout << "\nM(Range(0,4),Range(0,4)) =\n"
       << M(Range(0, 4), Range(0, 4)) << "\n";

  cout << "\nFilling M(Range(4,2,2),Range(6,3)) with junk.\n";

  MatrixView s = M(Range(4, 2, 2), Range(6, 3));
  fill_with_junk(s);

  cout << "\nM =\n" << M << "\n";

  const Matrix C = M;

  cout << "\nC(Range(3,4,2),Range(2,3,3)) =\n"
       << C(Range(3, 4, 2), Range(2, 3, 3)) << "\n";

  cout << "\nC(Range(3,4,2),Range(2,3,3)).transpose() =\n"
       << transpose(C(Range(3, 4, 2), Range(2, 3, 3))) << "\n";

  return 0;
}

void test2() {
  Vector v(50000000);

  cout << "v.nelem() = " << v.nelem() << "\n";

  cout << "Filling\n";
  //   for (Index i=0; i<v.nelem(); ++i )
  //     v[i] = sqrt(i);
  v = 1.;
  cout << "Done\n";
}

void test4() {
  Vector a(10);
  Vector b(a.nelem());

  for (Index i = 0; i < a.nelem(); ++i) {
    a[i] = (Numeric)(i + 1);
    b[i] = (Numeric)(a.nelem() - i);
  }

  cout << "a = \n" << a << "\n";
  cout << "b = \n" << b << "\n";
  cout << "a*b \n= " << a * b << "\n";

  Matrix A(11, 6);
  Matrix B(10, 20);
  Matrix C(20, 5);

  B = 2;
  C = 3;
  mult(A(Range(1, joker), Range(1, joker)), B, C);

  //  cout << "\nB =\n" << B << "\n";
  //  cout << "\nC =\n" << C << "\n";
  cout << "\nB*C =\n" << A << "\n";
}

void test5() {
  Vector a(10);
  Vector b(20);
  Matrix M(10, 20);

  // Fill b and M with a constant number:
  b = 1;
  M = 2;

  cout << "b = \n" << b << "\n";
  cout << "M =\n" << M << "\n";

  mult(a, M, b);  // a = M*b
  cout << "\na = M*b = \n" << a << "\n";

  mult(transpose((MatrixView)b), transpose((MatrixView)a), M);  // b^t = a^t * M
  cout << "\nb^t = a^t * M = \n" << transpose((MatrixView)b) << "\n";
}

void test6() {
  Index n = 5000;
  Vector x=uniform_grid(1, n, 1), y(n);  // FIXME: Function call rather than constructor OK?
  Matrix M(n, n);
  M = 1;
  //  cout << "x = \n" << x << "\n";

  cout << "Transforming.\n";
  //  transform(x,sin,x);
  // transform(transpose(y),sin,transpose(x));
  //  cout << "sin(x) =\n" << y << "\n";
  for (Index i = 0; i < 1000; ++i) {
    //      mult(y,M,x);
    transform(y, sin, static_cast<MatrixView>(x));
    x += 1;
  }
  //  cout << "y =\n" << y << "\n";

  cout << "Done.\n";
}

void test7() {
  Vector x=uniform_grid(1, 20000000, 1);  // FIXME: Function call rather than constructor OK?
  Vector y(x.nelem());
  transform(y, sin, x);
  cout << "min(sin(x)), max(sin(x)) = " << min(y) << ", " << max(y) << "\n";
  const auto [mn, mx] = minmax(y);
  cout << "minmax(sin(x)) = " << mn << ", " << mx << "\n";
}

void test8() {
  Vector x(80000000);
  for (Index i = 0; i < x.nelem(); ++i) x[i] = (Numeric)i;
  cout << "Done."
       << "\n";
}

void test9() {
  // Initialization of Matrix with view of other Matrix:
  Matrix A(4, 8);
  Matrix B(A(Range(joker), Range(0, 3)));
  cout << "B = " << B << "\n";
}

void test10() {
  // Initialization of Matrix with a vector (giving a 1 column Matrix).

  // At the moment doing this with a non-const Vector will result in a
  // warning message.
  Vector v=uniform_grid(1, 8, 1);  // FIXME: Function call rather than constructor OK?
  Matrix M(v);
  cout << "M = " << M << "\n";
}

void test11() {
  // Assignment between Vector and Matrix:

  // At the moment doing this with a non-const Vector will result in a
  // warning message.
  Vector v=uniform_grid(1, 8, 1);  // FIXME: Function call rather than constructor OK?
  Matrix M(v.nelem(), 1);
  M = MatrixView{v};  // FIXME:  Implicit cast to larger view is not allowed anymore OK?
  cout << "M = " << M << "\n";
}

void test13() {
  // Mix vector and one-column matrix in += operator.
  const Vector v=uniform_grid(1, 8, 1);  // The const is necessary here to
                            // avoid compiler warnings about
                            // different conversion paths.  // FIXME: Function call rather than constructor OK?
  Matrix M(v);
  M += MatrixView{v};  // FIXME:  Implicit cast to larger view is not allowed anymore OK?
  cout << "M = \n" << M << "\n";
}

void test17() {
  // Test Sum.
  Vector a=uniform_grid(1, 10, 1);  // FIXME: Function call rather than constructor OK?
  cout << "sum(a) = " << sum(a) << "\n";  // FIXME: Function call rather than member call OK?
}

void test18() {
  // Test elementvise square of a vector:
  Vector a=uniform_grid(1, 10, 1);  // FIXME: Function call rather than constructor OK?
  a *= a;
  cout << "a *= a =\n" << a << "\n";
}

void test19() {
  // There exists no explicit filling constructor of the form
  // Vector a(3,1.7).
  // But you can use the more general filling constructor with 3 arguments.

  Vector a=uniform_grid(1, 10, 1);  // FIXME: Function call rather than constructor OK?
  Vector b=uniform_grid(5.3, 10, 0);  // FIXME: Function call rather than constructor OK?
  cout << "a =\n" << a << "\n";
  cout << "b =\n" << b << "\n";
}

void test20() {
  // Test initialization list constructor:
  Vector a{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  cout << "a =\n" << a << "\n";
}

void test23() {
  // Test constructors that fills with constant:
  Vector a(10, 3.5);
  cout << "a =\n" << a << "\n";
  Matrix b(10, 10, 4.5);
  cout << "b =\n" << b << "\n";
}

void test24() {
  // Try element-vise multiplication of Matrix and Vector:
  Matrix a(5, 1, 2.5);
  Vector b=uniform_grid(1, 5, 1);  // FIXME: Function call rather than constructor OK?
  a *= MatrixView{b};  // FIXME:  Implicit cast to larger view is not allowed anymore OK?
  cout << "a*=b =\n" << a << "\n";
  a /= ConstMatrixView{b};  // FIXME:  Implicit cast to larger view is not allowed anymore OK?
  cout << "a/=b =\n" << a << "\n";
  a += ExhaustiveMatrixView{b};  // FIXME:  Implicit cast to larger view is not allowed anymore OK?
  cout << "a+=b =\n" << a << "\n";
  a -= ExhaustiveConstMatrixView{b};  // FIXME:  Implicit cast to larger view is not allowed anymore OK?
  cout << "a-=b =\n" << a << "\n";
}

void test28() {
  cout << "Test default constructor for Matrix:\n";
  Matrix a;
  Matrix b(a);
  cout << "b =\n" << b << "\n";
}

void test30() {
  cout << "Test Matrices of size 0:\n";
  Matrix a(0, 0);
  //  cout << "a(0,0) =\n" << a(0,0) << "\n";
  a.resize(2, 2);
  a = 1;
  cout << "a =\n" << a << "\n";

  Matrix b(3, 0);
  //  cout << "b(0,0) =\n" << b(0,0) << "\n";
  b.resize(b.nrows(), b.ncols() + 3);
  b = 2;
  cout << "b =\n" << b << "\n";

  Matrix c(0, 3);
  //  cout << "c(0,0) =\n" << c(0,0) << "\n";
  c.resize(c.nrows() + 3, c.ncols());
  c = 3;
  cout << "c =\n" << c << "\n";
}

void test31() {
  cout << "Test Tensor3:\n";

  Tensor3 a(2, 3, 4, 1.0);

  Index fill = 0;

  // Fill with some numbers
  for (Index i = 0; i < a.npages(); ++i)
    for (Index j = 0; j < a.nrows(); ++j)
      for (Index k = 0; k < a.ncols(); ++k) a(i, j, k) = (Numeric)(++fill);

  cout << "a =\n" << a << "\n";

  cout << "Taking out first row of first page:\n"
       << a(0, 0, Range(joker)) << "\n";

  cout << "Taking out last column of second page:\n"
       << a(1, Range(joker), a.ncols() - 1) << "\n";

  cout << "Taking out the first letter on every page:\n"
       << a(Range(joker), 0, 0) << "\n";

  cout << "Taking out first page:\n"
       << a(0, Range(joker), Range(joker)) << "\n";

  cout << "Taking out last row of all pages:\n"
       << a(Range(joker), a.nrows() - 1, Range(joker)) << "\n";

  cout << "Taking out second column of all pages:\n"
       << a(Range(joker), Range(joker), 1) << "\n";

  a *= 2;

  cout << "After element-vise multiplication with 2:\n" << a << "\n";

  transform(a, sqrt, a);

  cout << "After taking the square-root:\n" << a << "\n";

  Index s = 200;
  cout << "Let's allocate a large tensor, "
       << (Numeric)(s * s * s * 8) / 1024. / 1024. << " MB...\n";

  a.resize(s, s, s);

  cout << "Set it to 1...\n";

  a = 1;

  cout << "a(900,900,900) = " << a(90, 90, 90) << "\n";

  fill = 0;

  cout << "Fill with running numbers, using for loops...\n";
  for (Index i = 0; i < a.npages(); ++i)
    for (Index j = 0; j < a.nrows(); ++j)
      for (Index k = 0; k < a.ncols(); ++k) a(i, j, k) = (Numeric)(++fill);

  cout << "Max(a) = ...\n";

  cout << max(a) << "\n";
}

void junk4(Tensor4View&& a) { cout << "Describe a: " << describe(a) << "\n"; }

void junk2(ConstVectorView&& a) { cout << "Describe a: " << describe(a) << "\n"; }

void test34() {
  cout << "Test, if dimension expansion works implicitly.\n";

  Tensor3 t3(2, 3, 4);
  junk4(Tensor4View{t3});  // FIXME:  Implicit cast to larger view is not allowed anymore OK?

  Numeric x;
  junk2(ConstVectorView(x));
}

void test35() {
  cout << "Test the new copy semantics.\n";
  Vector a=uniform_grid(1, 4, 1);  // FIXME: Function call rather than constructor OK?
  Vector b;

  b = a;
  cout << "b = " << b << "\n";

  Vector aa=uniform_grid(1, 5, 1);  // FIXME: Function call rather than constructor OK?
  ConstVectorView c = aa;
  b = c;
  cout << "b = " << b << "\n";

  Vector aaa=uniform_grid(1, 6, 1);  // FIXME: Function call rather than constructor OK?
  VectorView d = aaa;
  b = d;
  cout << "b = " << b << "\n";
}

void test36() {
  cout << "Test using naked joker on Vector.\n";
  Vector a=uniform_grid(1, 4, 1);  // FIXME: Function call rather than constructor OK?
  VectorView b = a[joker];
  cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";
}

void test37() {
  Vector v1=uniform_grid(5e-15, 10, 0.42e-15 / 11);  // FIXME: Function call rather than constructor OK?
  Vector v2 = v1;
  Index i = 10'000'000;

  v1 /= (Numeric)i;
  v2 /= (Numeric)i;
  cout.precision(12);
  //  cout.setf(ios_base::scientific,ios_base::floatfield);
  v1 *= v1;
  v2 *= v2;
  cout << v1 << endl;
  cout << v2 << endl;
}

void test38() {
  Vector v(5, 0.);
  Numeric* const a = v.data_handle();  // FIXME:  Use data_handle instead of get_c_array OK?

  a[4] = 5.;

  cout << v << endl;
  cout << endl << "========================" << endl << endl;

  Matrix m(5, 5, 0.);
  Numeric* const b = m.data_handle();

  b[4] = 5.;

  cout << m << endl;
  cout << endl << "========================" << endl << endl;

  Tensor3 t3(5, 6, 7, 0.);
  Numeric* const c = t3.data_handle();

  c[6] = 5.;

  cout << t3 << endl;
}

void test39() {
  Vector v1=uniform_grid(1, 5, 1), v2(5);  // FIXME: Function call rather than constructor OK?

  v2 = v1 * 2;  // FIXME: This both compiles and works properly now OK?
  // Unfortunately, this thing compiles, but at least it gives an
  // assertion failure at runtime. I think what happens is that it
  // implicitly creates a one element vector out of the "2", then
  // tries to do a scalar product with v1.

  cout << v2 << endl;
}

void test40() {
  Vector v(100);

  v = 5;

  cout << v << endl;
}

void test41() {
  const Vector v1(10, 1);

  ConstVectorView vv = v1[Range(0, 5)];

  cout << "Vector:     " << v1 << endl;
  cout << "VectorView: " << vv << endl;

  try {
    // vv = 2;    // FIXME: This doesn't work anymore OK?
  } catch (const std::runtime_error& e) {
    std::cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  cout << "VectorView: " << vv << endl;
  cout << "Vector:     " << v1 << endl;
}

// Test behaviour of VectorView::operator MatrixView, for which I have fixed
// a bug today.
// SAB 2013-01-18
void test42() {
  cout << "test42\n";
  Vector x{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  cout << "x: " << x << endl;

  VectorView y = x[Range(2, 4, 2)];
  cout << "y: " << y << endl;

  ConstMatrixView z((MatrixView)y);
  cout << "z:\n" << z << endl;

  cout << "Every other z:\n" << z(Range(1, 2, 2), joker) << endl;

  // Try write access also:
  MatrixView zz(y);
  zz(Range(1, 2, 2), joker) = 0;
  cout << "zz:\n" << zz << endl;
  cout << "New x: " << x << endl;
}

void test44() {
#define docheck(fn, val, expect)                                           \
  cout << #fn << "(" << val << ") = " << fn(x) << " (expected " << #expect \
       << ")" << std::endl;

  Vector x{1, 2, 3};
  docheck(is_increasing, x, 1) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 1)

          x = {3, 2, 1};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 1)
      docheck(is_sorted, x, 0)

          x = {1, 2, 2};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 1)

          x = {2, 2, 1};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

          x = {1, NAN, 2};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

          x = {2, NAN, 1};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

          x = {NAN, NAN, NAN};
  docheck(is_increasing, x, 0) docheck(is_decreasing, x, 0)
      docheck(is_sorted, x, 0)

#undef docheck
}

void test_view() {
  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    ConstTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }
    
    Numeric d=1.0;
    for (auto a: y) for (auto b: a) for (auto c: b) {ARTS_ASSERT(c == d, c, ' ', d); d+=1.0;}
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    Tensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto& g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }
    
    Numeric d=1.0;
    for (const auto& a: y) for (auto b: a) for (auto c: b) {ARTS_ASSERT(c == d, c, ' ', d); d+=1.0;}
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    ExhaustiveConstTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric>);
    }
    
    Numeric d=1.0;
    for (auto a: y) for (auto b: a) for (auto c: b) {ARTS_ASSERT(c == d, c, ' ', d); d+=1.0;}
  }

  //! Simply test that some standard operators work
  {
    std::vector<Numeric> x{1,2,3,4,5, 6, 7, 8};
    ExhaustiveTensor3View y{matpack::exhaustive_mdspan<Numeric, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    for (Index i=0; i<2; i++) {
      [[maybe_unused]] auto& g = y(i, joker, matpack::matpack_strided_access{joker})(i, i);
      [[maybe_unused]] auto& h = y[i][i][i];
      static_assert(std::same_as<decltype(g), decltype(h)>);
      static_assert(std::same_as<decltype(g), Numeric&>);
    }
    
    Numeric d=1.0;
    for (const auto& a: y) for (auto b: a) for (auto c: b) {ARTS_ASSERT(c == d, c, ' ', d); d+=1.0;}
  }

  //! Allow assignment between views
  {
    std::vector<std::complex<Numeric>> x{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ExhaustiveComplexTensor3View y{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yr{matpack::exhaustive_mdspan<Numeric, 3>{xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yi{matpack::exhaustive_mdspan<Numeric, 3>{xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() == yi)
    ARTS_ASSERT(y.imag() != yr)

    y.imag() = y.real();
    
    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() == y.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() != yi)
    ARTS_ASSERT(y.imag() == yr)
  }
  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    std::vector<Numeric> xr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> xi{0, 0, 0, 0, 0, 0, 0, 0};
    ExhaustiveComplexTensor3View y{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yr{matpack::exhaustive_mdspan<Numeric, 3>{xr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View yi{matpack::exhaustive_mdspan<Numeric, 3>{xi.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    auto z = y(joker, matpack::matpack_strided_access(0,-1,2), joker);

    std::vector<Numeric> reszr{1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<Numeric> reszi{1, 2, 0, 0, 5, 6, 0, 0};
    const ExhaustiveTensor3View zr{matpack::exhaustive_mdspan<Numeric, 3>{reszr.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    const ExhaustiveTensor3View zi{matpack::exhaustive_mdspan<Numeric, 3>{reszi.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(z.real() != z.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() == yi)
    ARTS_ASSERT(y.real() == zr)
    ARTS_ASSERT(y.imag() != zi)

    z.imag() = z.real();

    ARTS_ASSERT(yr != yi)
    ARTS_ASSERT(y.real() != y.imag())
    ARTS_ASSERT(z.real() == z.imag())
    ARTS_ASSERT(y.real() == yr)
    ARTS_ASSERT(y.imag() != yi)
    ARTS_ASSERT(y.real() == zr)
    ARTS_ASSERT(y.imag() == zi)
  }

  //! Test basic +=, -=, *=, /=, = for arithmetic change or assignment
  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ExhaustiveComplexTensor3View y{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{x.data(), std::array{Index{2}, Index{2}, Index{2}}}};
    
    std::vector<std::complex<Numeric>> copy = x;
    const ExhaustiveComplexTensor3View test{matpack::exhaustive_mdspan<std::complex<Numeric>, 3>{copy.data(), std::array{Index{2}, Index{2}, Index{2}}}};

    ARTS_ASSERT(y == test)

    y += 3;
    for (auto& a: copy) a += 3;
    ARTS_ASSERT(y == test)

    y /= 3;
    for (auto& a: copy) a /= 3;
    ARTS_ASSERT(y == test)

    y.imag() -= 2;
    for (auto& a: copy) a -= std::complex<Numeric>(0, 2);
    ARTS_ASSERT(y == test)

    y -= 5;
    for (auto& a: copy) a -= 5;
    ARTS_ASSERT(y == test)

    y *= 5;
    for (auto& a: copy) a *= 5;
    ARTS_ASSERT(y == test)

    y = 42;
    for (auto& a: copy) a = 42;
    ARTS_ASSERT(y == test)
  }

  //! Test external assignment
  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexVectorView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{x.data(), 8}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexVectorView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{copy.data(), 8}};

    for(auto& a: copy) a+=std::complex<Numeric>(3, 2);

    ARTS_ASSERT(y != test)
    
    y = copy;

    ARTS_ASSERT(y == test)
  }
}

void test_eigen() {
  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y);

    ARTS_ASSERT(y != test)
    for (auto& a: copy) a *= 2;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexVectorView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 1>{x.data(), 8}};

    const std::complex<Numeric> z = y * y;
    ARTS_ASSERT(z == std::transform_reduce(x.begin(), x.end(), x.begin(), std::complex<Numeric>{0}))
  }

  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y) + y;

    ARTS_ASSERT(y != test)
    for (auto& a: copy) a = 2.0 * a + a;
    ARTS_ASSERT(y == test)
  }

  {
    std::vector<std::complex<Numeric>> x{1,2,3,4,5, 6, 7, 8};
    ComplexMatrixView y{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{x.data(), std::array{Index{2}, Index{4}}}};
    std::vector<std::complex<Numeric>> copy = x;
    const ComplexMatrixView test{matpack::exhaustive_mdspan<std::complex<Numeric>, 2>{copy.data(), std::array{Index{2}, Index{4}}}};

    ARTS_ASSERT(y == test)

    y = 2 * matpack::eigen::as_eigen(y) - y;
    
    ARTS_ASSERT(y == test)
  }
}

void test_data() {
  matpack::matpack_data<Numeric, 3> x(4, 2, 3, 2.1);
  std::cout << x << '\n' << '\n';
  std::cout << x[0]  << '\n' << '\n';
  std::cout << x[0][0]  << '\n' << '\n';
  std::cout << x[0][0][0]  << '\n' << '\n';
  x[0][0][0] += 1;
  std::cout << x << '\n' << '\n';

  auto y = std::move(x).flatten();
  std::cout << y << '\n';
  std::cout << x << '\n';

  matpack::matpack_data<Numeric, 1> z;
  z.swap(y);
}

void test_complex() {
  {
    Complex x{0, 0};
    const Complex y{0, 0};
    ARTS_ASSERT(x == y)
    ARTS_ASSERT(real_val(x) == real_val(y))
    ARTS_ASSERT(imag_val(x) == imag_val(y))
  }
  
  {
    Complex x{1, 1};
    x = 2. + x;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = 2. - x;
    ARTS_ASSERT(x == (Complex{-1, -1}))
    x = 2. * x;
    ARTS_ASSERT(x == (Complex{-2, -2}))
    x = 2. / x;
    ARTS_ASSERT(x == (Complex{-0.5, 0.5}))
  }

  {
    Complex x{1, 1};
    x = x + 2.;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = x - 2.;
    ARTS_ASSERT(x == (Complex{1, 1}))
    x = x * 2.;
    ARTS_ASSERT(x == (Complex{2, 2}))
    x = x / 2.;
    ARTS_ASSERT(x == (Complex{1, 1}))
  }

  {
    Complex x{1, 1};
    x = 2 + x;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = 2 - x;
    ARTS_ASSERT(x == (Complex{-1, -1}))
    x = 2 * x;
    ARTS_ASSERT(x == (Complex{-2, -2}))
    x = 2 / x;
    ARTS_ASSERT(x == (Complex{-0.5, 0.5}))
  }

  {
    Complex x{1, 1};
    x = x + 2;
    ARTS_ASSERT(x == (Complex{3, 1}))
    x = x - 2;
    ARTS_ASSERT(x == (Complex{1, 1}))
    x = x * 2;
    ARTS_ASSERT(x == (Complex{2, 2}))
    x = x / 2;
    ARTS_ASSERT(x == (Complex{1, 1}))
  }

  {
    Complex x{3, 3};
    const std::complex<int> y{3, 3};
    x = x + y;
    ARTS_ASSERT(x == (Complex{6, 6}))
    x = x - y;
    ARTS_ASSERT(x == (Complex{3, 3}))
    x = x * y;
    ARTS_ASSERT(x == (Complex{0, 18}))
    x = x / y;
    ARTS_ASSERT(x == (Complex{3, 3}))
  }
}

void test_math() {
  {
    Vector x{0, 1, 2, 3, 4};
    ARTS_ASSERT(min(x) == 0)
    ARTS_ASSERT(min(ExhaustiveVectorView{x}) == 0)
    ARTS_ASSERT(min(ExhaustiveConstVectorView{x}) == 0)
    ARTS_ASSERT(min(VectorView{x}) == 0)
    ARTS_ASSERT(min(ConstVectorView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.begin())>);

    static_assert(std::random_access_iterator<decltype(x.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.elem_end())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.elem_end())>);

    static_assert(std::random_access_iterator<decltype(x.end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.end())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.end())>);
  }

/*
  {
    Matrix x(3, 3, 0);
    ARTS_ASSERT(min(x) == 0)
    ARTS_ASSERT(min(ExhaustiveVectorView{x}) == 0)
    ARTS_ASSERT(min(ExhaustiveConstVectorView{x}) == 0)

    static_assert(std::random_access_iterator<decltype(x.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.elem_begin())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.elem_begin())>);

    static_assert(std::random_access_iterator<decltype(x.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ExhaustiveConstVectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(VectorView{x}.begin())>);
    static_assert(std::random_access_iterator<decltype(ConstVectorView{x}.begin())>);
  } */
}

void test_mult() {
  {
    Matrix A(4,4,1);
    for (Index i=0; i<4; i++) for (Index j=0; j<4; j++) A(i, j) = static_cast<Numeric>(i);
    Vector y({1,2,3,4}), x(4);

    Vector eig_x=A*y;
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, eig_x, " vs ", x)
  }
  {
    ComplexMatrix A(4,4,1);
    for (Index i=0; i<4; i++) for (Index j=0; j<4; j++) A(i, j) = {static_cast<Numeric>(i), 2*static_cast<Numeric>(i)};
    ComplexVector y({1,2,3,4}), x(4);

    ComplexVector eig_x=A*y;
    mult(x, A, y);
    ARTS_USER_ERROR_IF(eig_x not_eq x, eig_x, " vs ", x)
  }
}

#define EXECUTE_TEST(X) \
std::cout << "#########################################################\n";\
std::cout << "Executing test: " #X << '\n'; \
std::cout << "#########################################################\n";\
X();\
std::cout << "#########################################################\n";

int main() {
  EXECUTE_TEST(test1)
  EXECUTE_TEST(test2)
  EXECUTE_TEST(test4)
  EXECUTE_TEST(test5)
  EXECUTE_TEST(test6)
  EXECUTE_TEST(test7)
  EXECUTE_TEST(test8)
  EXECUTE_TEST(test9)
  EXECUTE_TEST(test10)
  EXECUTE_TEST(test11)
  EXECUTE_TEST(test13)
  EXECUTE_TEST(test17)
  EXECUTE_TEST(test18)
  EXECUTE_TEST(test19)
  EXECUTE_TEST(test20)
  EXECUTE_TEST(test23)
  EXECUTE_TEST(test24)
  EXECUTE_TEST(test28)
  EXECUTE_TEST(test30)
  EXECUTE_TEST(test31)
  EXECUTE_TEST(test34)
  EXECUTE_TEST(test35)
  EXECUTE_TEST(test36)
  EXECUTE_TEST(test37)
  EXECUTE_TEST(test38)
  EXECUTE_TEST(test39)
  EXECUTE_TEST(test40)
  EXECUTE_TEST(test41)
  EXECUTE_TEST(test42)
  EXECUTE_TEST(test44)

  EXECUTE_TEST(test_view)
  EXECUTE_TEST(test_eigen)
  EXECUTE_TEST(test_data)
  EXECUTE_TEST(test_complex)
  EXECUTE_TEST(test_math)
  EXECUTE_TEST(test_mult)

  matpack::flat_shape_pos<4> pos{{3,4,5,6}};

  for (Index i=0; i<100; i++) {
    std::cout << pos << " " << pos.to_index() << '\n';
    pos++;
  }
}

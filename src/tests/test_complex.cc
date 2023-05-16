/*!
  \file   test_complex.cc
  \author  <sbuehler@ltu.se>
  \date   Sat Dec 14 19:42:25 2002
  
  \brief  Test the complex numbers.
*/

#include "matpack_complex.h"
#include "matpack_data.h"
#include "matpack_eigen.h"
#include "matpack_math.h"
#include <iostream>

using std::cout;

void test01() {
  Complex a;
  Complex b(3., 0.);

  //cout << "a = " << a << "\n";
  cout << "b = " << b << "\n";

  a = b;
  cout << "a = " << a << "\n";

  a = Complex(1., 1.);
  cout << "a = " << a << "\n";
  cout << "a.abs() = " << abs(a) << "\n";
  cout << "a.arg()= " << arg(a) * 57.2957914331333 << " degree"
       << "\n";

  Complex c;
  c = a + b;
  cout << "c = " << c << "\n";

  cout << "a += b: " << (a += b) << "\n";

  cout << "c + 3 = " << (c + 3.) << "\n";

  cout << "c + 3i = c + Complex (0,3) = " << (c + Complex(0., 3.)) << "\n";

  Complex d;
  a = Complex(1., 1.);
  b = Complex(-1., -1.);
  d = a * b;
  cout << "d = " << d << "\n";
  Complex e;
  e = a / b;
  cout << "e = " << e << "\n";
  Complex f;
  f = pow(a, 0) + pow(b, 0);
  cout << "f = " << f << "\n";
  Complex g;
  g = pow(a, 1) + pow(b, 1);
  cout << "g = " << g << "\n";
  Complex h;
  h = pow(a, 2) + pow(b, 2) + pow(a, 3) + pow(b, 3);
  cout << "h = " << h << "\n";

  ComplexMatrix A(4, 4);
  ComplexMatrix B(4, 4);
  ComplexMatrix X1(4, 1);
  ComplexVector X2(4);
  ComplexMatrix C(4, 4);
  for (Index i = 0; i < 4; i++) {
    X1(i, 0) = Complex((Numeric)i * 5.0 + 2.0, (Numeric)i + 1.0);
    X2[i] = Complex((Numeric)i * 5.0 + 2.0, (Numeric)i + 1.0);
    for (Index j = 0; j < 4; j++) {
      A(i, j) = Complex((Numeric)i + (Numeric)j, 0.0) + Complex(1.0, 1.0);
      B(i, j) =
          Complex(2.0 * (Numeric)i + 4.0 * (Numeric)j, 0.0) + Complex(3.0, 3.0);
    }
  }

  std::cout << "\n";
  mult(C, A, B);
  std::cout << matpack::eigen::mat(A) << "\nx\n"
            << matpack::eigen::mat(B) << "\n=\n"
            << matpack::eigen::mat(C) << "\n\n";
  mult(C, C, C);
  mult(C, C, C);
  std::cout
      << "Same matrix can be both input and output, here is the above to the power of 4:\n"
      << matpack::eigen::mat(C) << "\n\n";
  mult(C(joker, 1), A, X2);
  std::cout << matpack::eigen::mat(A) << "\nx\n"
            << matpack::eigen::row_vec(X2) << "\n=\n"
            << matpack::eigen::row_vec(C(joker, 1)) << "\n\n";
  mult(C(Range(1, 3), Range(0, 3)),
       A(Range(1, 3), Range(1, 2)),
       B(Range(1, 2), Range(1, 3)));
  std::cout << matpack::eigen::mat(A(Range(1, 3), Range(1, 2))) << "\nx\n"
            << matpack::eigen::mat(B(Range(1, 2), Range(1, 3))) << "\n=\n"
            << matpack::eigen::mat(C(Range(1, 3), Range(0, 3))) << "\n\n";
}

void test02() {
  constexpr Numeric nul = 0.f;
  constexpr Numeric one = 1.f;
  constexpr Numeric two = 2.f;
  constexpr Numeric tre = 3.f;
  constexpr Numeric fyr = 4.f;

  constexpr Complex r(one, nul);
  constexpr Complex i(nul, one);

  constexpr Complex a(one, two);
  constexpr Complex b(tre, two);
  constexpr Complex c(tre, fyr);

  // Test that every operator is treated as Numeric regardless of basic types
  static_assert(tre + a == 3 + a, "Bad operator+ int Complex");
  static_assert(tre - a == 3 - a, "Bad operator- int Complex");
  static_assert(tre * a == 3 * a, "Bad operator* int Complex");
  static_assert(tre / a == 3 / a, "Bad operator/ int Complex");
  static_assert(a + fyr == a + 4, "Bad operator+ Complex int");
  static_assert(a - fyr == a - 4, "Bad operator- Complex int");
  static_assert(a * fyr == a * 4, "Bad operator* Complex int");
  static_assert(a / fyr == a / 4, "Bad operator/ Complex int");
  static_assert(tre + a == 3.f + a, "Bad operator+ float Complex");
  static_assert(tre - a == 3.f - a, "Bad operator- float Complex");
  static_assert(tre * a == 3.f * a, "Bad operator* float Complex");
  static_assert(tre / a == 3.f / a, "Bad operator/ float Complex");
  static_assert(a + fyr == a + 4.f, "Bad operator+ Complex float");
  static_assert(a - fyr == a - 4.f, "Bad operator- Complex float");
  static_assert(a * fyr == a * 4.f, "Bad operator* Complex float");
  static_assert(a / fyr == a / 4.f, "Bad operator/ Complex float");

  // Test the most basic test of (1, 0)*z and (0, 1)*z
  static_assert(r * a == Complex(one * one, one * two),
                "Bad Complex(1, 0) Complex ");
  static_assert(i * a == Complex(-one * two, one * one),
                "Bad Complex(0, 1) Complex ");

  // Test key helper function
  static_assert(abs2(c) == tre * tre + fyr * fyr, "Bad Numeric abs2(Complex)");

  // Test Complex-Complex operators
  static_assert(a + b == Complex(one + tre, two + two),
                "Bad operator+ Complex Complex");
  static_assert(a - b == Complex(one - tre, two - two),
                "Bad operator- Complex Complex");
  static_assert(a * b == Complex(one * tre - two * two, one * two + two * tre),
                "Bad operator* Complex Complex");
  static_assert(
      a / b == Complex(one * tre + two * two, two * tre - one * two) / abs2(b),
      "Bad operator/ Complex Complex");
}

void test03()
{
  constexpr Index n=4;
  constexpr Index col=2;
  constexpr Index row=2;
  static_assert(n>=4, "Req for size of val to follow");
  static_assert(col<n, "Req for size of val to follow");
  static_assert(row<n, "Req for size of val to follow");
  
  constexpr Range r0 = Range(joker);
  constexpr Range r1 = Range(1, 2);
  constexpr Range r2 = Range(1, 2, 2);
  constexpr Range r3 = Range(n-1, 2, -1);
  constexpr Range r4 = Range(n-1, 2, -2);
  constexpr Range r5 = Range(n-1, n, -1);
  const auto ranges={r0, r1, r2, r3, r4, r5};
  
  ComplexVector V(n, static_cast<Numeric>(n));
  for (Index i=0; i<n; i++)
    V[i] = Complex(Numeric(i), Numeric(i+2*n));
  const ComplexVector cV=V;
  
  ComplexMatrix M(n, n);
  for (Index i=0; i<n; i++)
    for (Index j=0; j<n; j++)
      M(i, j) = Complex(Numeric(i), Numeric(j+2*n));
  const ComplexMatrix cM=M;
  
  std::cout<<"Vector: ";
  std::cout<<matpack::eigen::row_vec(V).transpose()<<'\n';
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  
  std::cout<<"Matrix:\n";
  std::cout<<matpack::eigen::mat(M)<<"\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"Vector\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"Vector["<<r<<"]: "<<matpack::eigen::row_vec(V[r]).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(V[r].real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(V[r].imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"const Vector\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"const Vector["<<r<<"]: "<<matpack::eigen::row_vec(cV[r]).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(V[r].real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(V[r].imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"RowView Matrix\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"Matrix("<<row<<", "<<r<<"): "<<matpack::eigen::row_vec(M(row,r)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(M(row,r).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(M(row,r).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"ColView Matrix\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"Matrix("<<r<<", "<<col<<"): "<<matpack::eigen::row_vec(M(r, col)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(M(r, col).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(M(r, col).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"RowView Matrix-transpose\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"transpose(Matrix)("<<row<<", "<<r<<"): "<<matpack::eigen::row_vec(transpose(M)(row,r)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(transpose(M)(row,r).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(transpose(M)(row,r).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"ColView Matrix transpose\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"transpose(Matrix)("<<r<<", "<<col<<"): "<<matpack::eigen::row_vec(transpose(M)(r, col)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(transpose(M)(r, col).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(transpose(M)(r, col).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"Sub-Matrix\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto rr: ranges) {
    for (auto rc: ranges) {
      std::cout<<"Matrix("<<rr<<", "<<rc<<"):\n"<<matpack::eigen::mat(M(rr,rc)).transpose()<<'\n';
      std::cout<<"Real:\n"<<matpack::eigen::mat(M(rr,rc).real()).transpose()<<'\n';
      std::cout<<"Imag:\n"<<matpack::eigen::mat(M(rr,rc).imag()).transpose()<<'\n';
      std::cout<<"---------------------------------------------------\n";
    }
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"Sub-Matrix transpose\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto rr: ranges) {
    for (auto rc: ranges) {
      std::cout<<"transpose(Matrix)("<<rr<<", "<<rc<<"):\n"<<matpack::eigen::mat(transpose(M)(rr,rc)).transpose()<<'\n';
      std::cout<<"Real:\n"<<matpack::eigen::mat(transpose(M)(rr,rc).real()).transpose()<<'\n';
      std::cout<<"Imag:\n"<<matpack::eigen::mat(transpose(M)(rr,rc).imag()).transpose()<<'\n';
      std::cout<<"---------------------------------------------------\n";
    }
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"RowView Matrix const\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"const Matrix("<<row<<", "<<r<<"): "<<matpack::eigen::row_vec(cM(row,r)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(cM(row,r).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(cM(row,r).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"ColView Matrix const\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"const Matrix("<<r<<", "<<col<<"): "<<matpack::eigen::row_vec(cM(r, col)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(cM(r, col).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(cM(r, col).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"RowView Matrix-transpose const\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"transpose(const Matrix)("<<row<<", "<<r<<"): "<<matpack::eigen::row_vec(transpose(cM)(row,r)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(transpose(cM)(row,r).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(transpose(cM)(row,r).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"ColView Matrix-transpose const\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto r: ranges) {
    std::cout<<"transpose(const Matrix)("<<r<<", "<<col<<"): "<<matpack::eigen::row_vec(transpose(cM)(r, col)).transpose()<<'\n';
    std::cout<<"Real: "<<matpack::eigen::row_vec(transpose(cM)(r, col).real()).transpose()<<'\n';
    std::cout<<"Imag: "<<matpack::eigen::row_vec(transpose(cM)(r, col).imag()).transpose()<<'\n';
    std::cout<<"---------------------------------------------------\n";
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"Sub-Matrix const\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto rr: ranges) {
    for (auto rc: ranges) {
      std::cout<<"Matrix("<<rr<<", "<<rc<<"):\n"<<matpack::eigen::mat(M(rr,rc)).transpose()<<'\n';
      std::cout<<"Real:\n"<<matpack::eigen::mat(M(rr,rc).real()).transpose()<<'\n';
      std::cout<<"Imag:\n"<<matpack::eigen::mat(M(rr,rc).imag()).transpose()<<'\n';
      std::cout<<"---------------------------------------------------\n";
    }
  }
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  std::cout<<"Sub-Matrix transpose const\n";
  std::cout<<"----------------------------------------------------------------------------------------------------\n";
  for (auto rr: ranges) {
    for (auto rc: ranges) {
      std::cout<<"transpose(const Matrix)("<<rr<<", "<<rc<<"):\n"<<matpack::eigen::mat(transpose(cM)(rr,rc)).transpose()<<'\n';
      std::cout<<"Real:\n"<<matpack::eigen::mat(transpose(cM)(rr,rc).real()).transpose()<<'\n';
      std::cout<<"Imag:\n"<<matpack::eigen::mat(transpose(cM)(rr,rc).imag()).transpose()<<'\n';
      std::cout<<"---------------------------------------------------\n";
    }
  }
}

int main() {
//   test01();
  test03();

  Complex a(0, 0);
  real_val(a) += 1;
  imag_val(a) -= 2;
  imag_val(a) *= 2;
  std::cout << a << '\n';
  return 0;
}

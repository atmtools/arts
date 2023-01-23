#include "artstime.h"
#include "matpack_arrays.h"
#include "matpack_eigen.h"
#include "matpack_math.h"
#include "matpack_view.h"

#include <cstdlib>
#include <ostream>
#include <stdexcept>

struct Timing {
  std::string_view name;
  Timing(const char * c) : name(c) {}
  TimeStep dt{};
  template <typename Function> void operator()(Function&& f) {
    Time start{};
    f();
    Time end{};
    dt = end - start;
  }
};

std::ostream& operator<<(std::ostream& os, const std::vector<Timing>& vt) {
  for (auto& t: vt) if (t.name not_eq "dummy") os << t.name  << " : " << t.dt << '\n';
  return os;
}

std::vector<Timing> test_sum(Index N) {
  Numeric X;
  Vector a(N, 1);
  auto b = static_cast<ExhaustiveVectorView>(a);
  VectorView c = a;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  {
    X = sum(a) + sum(a) + sum(a) + sum(a) + sum(a);
  };
  xvec.push_back(X);

  out.emplace_back("sum(vec)")([&]() {
    X = sum(a);
  });
  xvec.push_back(X);

  out.emplace_back("sum(exview)")([&]() {
    X = sum(b);
  });
  xvec.push_back(X);

  out.emplace_back("sum(view)")([&]() {
    X = sum(c);
  });
  xvec.push_back(X);

  out.emplace_back("mean(vec)")([&]() {
    X = mean(a);
  });
  xvec.push_back(X);

  out.emplace_back("mean(exview)")([&]() {
    X = mean(b);
  });
  xvec.push_back(X);

  out.emplace_back("mean(view)")([&]() {
    X = mean(c);
  });
  xvec.push_back(X);

  out.emplace_back("nanmean(vec)")([&]() {
    X = nanmean(a);
  });
  xvec.push_back(X);

  out.emplace_back("nanmean(exview)")([&]() {
    X = nanmean(b);
  });
  xvec.push_back(X);

  out.emplace_back("nanmean(view)")([&]() {
    X = nanmean(c);
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing>  test_dot(Index N) {
  Vector a(N, 1);
  Vector b(N, 1);

  auto ae=static_cast<ExhaustiveVectorView>(a);
  auto be=static_cast<ExhaustiveVectorView>(b);

  VectorView av=a;
  VectorView bv=b;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    X = a * b + a * b + a * b + a * b + a * b + a * b;
  };
  xvec.push_back(X);

  out.emplace_back("vec * vec")([&]() {
    X = a * b;
  });
  xvec.push_back(X);

  out.emplace_back("exview * exview")([&]() {
    X = ae * be;
  });
  xvec.push_back(X);

  out.emplace_back("view * view")([&]() {
    X = av * bv;
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing>  test_vec_mult(Index N) {
  Vector x(N, 1);
  Matrix A(N, N, 1);
  Vector y(N, 1);

  auto xe=static_cast<ExhaustiveVectorView>(x);
  MatrixView Ae=static_cast<ExhaustiveMatrixView>(A);
  auto ye=static_cast<ExhaustiveVectorView>(y);

  VectorView xv=x;
  MatrixView Av=A;
  VectorView yv=y;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    mult(x, A, y);
    X = x[0];
    mult(x, A, y);
    X += x[0];
    mult(x, A, y);
    X += x[0];
    mult(x, A, y);
    X += x[0];
    mult(x, A, y);
    X += x[0];
  }
  xvec.push_back(X);

  out.emplace_back("vec = mat * vec")([&]() {
    mult(x, A, y);
    X = x[0] + x[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview = exmview * exview")([&]() {
    mult(xe, Ae, ye);
    X = xe[0] + xe[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view = mview * view")([&]() {
    mult(xv, Av, yv);
    X = xv[0] + xv[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing> test_mat_multiply(Index N) {
  Matrix A(N, N), C(N, N, 1), B(N, N, 1);
  auto Ae=static_cast<ExhaustiveMatrixView>(A), Be=static_cast<ExhaustiveMatrixView>(B), Ce=static_cast<ExhaustiveMatrixView>(C);
  MatrixView Av=A, Bv=B, Cv=C;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    mult(A, B, C);
    X = A(0, 0);
  }
  xvec.push_back(X);

  out.emplace_back("mat = mat * mat")([&]() {
    mult(A, B, C);
    X = A(0, 0) + A(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview = exmview * exmview")([&]() {
    mult(Ae, Be, Ce);
    X = Ae(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview = mview * mview")([&]() {
    mult(Av, Bv, Cv);
    X = Av(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing> test_elementary_ops_vec(Index N) {
  Vector A(N, 1);
  auto Ae = static_cast<ExhaustiveVectorView>(A);
  VectorView Av = A;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    A = 1;
    A += 1;
    A *= 2;
    A /= 2;
    A -= 1;
    X = A[0] + A[N-1];
  }
  xvec.push_back(X);

  out.emplace_back("vec = 1")([&]() {
    A = 1;
    X = A[0] + A[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview = 1")([&]() {
    Ae = 1;
    X = Ae[0] + Ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view = 1")([&]() {
    Av = 1;
    X = Av[0] + Av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec += 1.33")([&]() {
    A += 1.33;
    X = A[0] + A[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview += 1.33")([&]() {
    Ae += 1.33;
    X = Ae[0] + Ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view += 1.33")([&]() {
    Av += 1.33;
    X = Av[0] + Av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec -= 1.5")([&]() {
    A -= 1.5;
    X = A[0] + A[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview -= 1.5")([&]() {
    Ae -= 1.5;
    X = Ae[0] + Ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view -= 1.5")([&]() {
    Av -= 1.5;
    X = Av[0] + Av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec *= 3.5")([&]() {
    A *= 3.5;
    X = A[0] + A[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview *= 3.5")([&]() {
    Ae *= 3.5;
    X = Ae[0] + Ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view *= 3.5")([&]() {
    Av *= 3.5;
    X = Av[0] + Av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec /= 7.77")([&]() {
    A /= 7.77;
    X = A[0] + A[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview /= 7.77")([&]() {
    Ae /= 7.77;
    X = Ae[0] + Ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view /= 7.77")([&]() {
    Av /= 7.77;
    X = Av[0] + Av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing> test_elementary_ops_mat(Index N) {
  Matrix A(N, N, 1);
  auto Ae = static_cast<ExhaustiveMatrixView>(A);
  MatrixView Av = A;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    A = 1;
    A += 1;
    A *= 2;
    A /= 2;
    A -= 1;
    X = A(0, 0) + A(N-1, N-1);
  }
  xvec.push_back(X);

  out.emplace_back("mat = 1")([&]() {
    A = 1;
    X = A(0, 0) + A(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview = 1")([&]() {
    Ae = 1;
    X = Ae(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview = 1")([&]() {
    Av = 1;
    X = Av(0, 0) + Av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat += 1.33")([&]() {
    A += 1.33;
    X = A(0, 0) + A(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview += 1.33")([&]() {
    Ae += 1.33;
    X = Ae(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview += 1.33")([&]() {
    Av += 1.33;
    X = Av(0, 0) + Av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat -= 1.5")([&]() {
    A -= 1.5;
    X = A(0, 0) + A(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview -= 1.5")([&]() {
    Ae -= 1.5;
    X = Ae(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview -= 1.5")([&]() {
    Av -= 1.5;
    X = Av(0, 0) + Av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat *= 3.5")([&]() {
    A *= 3.5;
    X = A(0, 0) + A(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview *= 3.5")([&]() {
    Ae *= 3.5;
    X = Ae(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview *= 3.5")([&]() {
    Av *= 3.5;
    X = Av(0, 0) + Av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat /= 7.77")([&]() {
    A /= 7.77;
    X = A(0, 0) + A(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview /= 7.77")([&]() {
    Ae /= 7.77;
    X = Ae(0, 0) + Ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview /= 7.77")([&]() {
    Av /= 7.77;
    X = Av(0, 0) + Av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing> test_ops_vec(Index N) {
  Vector a(N, 1.1);
  Vector b(N, 2.2);
  auto ae=static_cast<ExhaustiveVectorView>(a);
  auto be=static_cast<ExhaustiveVectorView>(b);
  VectorView av=a;
  VectorView bv=b;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    a += b;
    a -= b;
    a *= b;
    a /= b;
    X = a[0] + a[N-1];
  }
  xvec.push_back(X);

  out.emplace_back("vec += vec")([&]() {
    a += b;
    X = a[0] + a[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview += exview")([&]() {
    ae += be;
    X = ae[0] + ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view += view")([&]() {
    av += bv;
    X = av[0] + av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec -= vec")([&]() {
    a -= b;
    X = a[0] + a[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview -= exview")([&]() {
    ae -= be;
    X = ae[0] + ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view -= view")([&]() {
    av -= bv;
    X = av[0] + av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec *= vec")([&]() {
    a *= b;
    X = a[0] + a[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview *= exview")([&]() {
    ae *= be;
    X = ae[0] + ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view *= exview")([&]() {
    av *= bv;
    X = av[0] + av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("vec /= vec")([&]() {
    a /= b;
    X = a[0] + a[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("exview /= exview")([&]() {
    ae /= be;
    X = ae[0] + ae[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("view /= view")([&]() {
    av /= bv;
    X = av[0] + av[N-1];
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

std::vector<Timing> test_ops_mat(Index N) {
  Matrix a(N, N, 1.1);
  Matrix b(N, N, 2.2);
  auto ae=static_cast<ExhaustiveMatrixView>(a);
  auto be=static_cast<ExhaustiveMatrixView>(b);
  MatrixView av=a;
  MatrixView bv=b;

  std::vector<Numeric> xvec;
  std::vector<Timing> out;

  Numeric X;
  {
    a += b;
    a -= b;
    a *= b;
    a /= b;
    X = a(0, 0) + a(N-1, N-1);
  }
  xvec.push_back(X);

  out.emplace_back("mat += mat")([&]() {
    a += b;
    X = a(0, 0) + a(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview += exmview")([&]() {
    ae += be;
    X = ae(0, 0) + ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview += mview")([&]() {
    av += bv;
    X = av(0, 0) + av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat -= mat")([&]() {
    a -= b;
    X = a(0, 0) + a(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview -= exmview")([&]() {
    ae -= be;
    X = ae(0, 0) + ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview -= mview")([&]() {
    av -= bv;
    X = av(0, 0) + av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat *= mat")([&]() {
    a *= b;
    X = a(0, 0) + a(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview *= exmview")([&]() {
    ae *= be;
    X = ae(0, 0) + ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview *= mview")([&]() {
    av *= bv;
    X = av(0, 0) + av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mat /= mat")([&]() {
    a /= b;
    X = a(0, 0) + a(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("exmview /= exmview")([&]() {
    ae /= be;
    X = ae(0, 0) + ae(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("mview /= mview")([&]() {
    av /= bv;
    X = av(0, 0) + av(N-1, N-1);
  });
  xvec.push_back(X);

  out.emplace_back("dummy")([xvec](){return xvec[xvec.size()-1];});

  return out;
}

int main(int argc, char** c) {
  std::array <Index, 8> N;
  if (static_cast<std::size_t>(argc) < 1 + 1 + N.size()) {
    std::cerr << "Expects PROGNAME NREPEAT NSIZE..., wehere NSIZE is " << N.size() << " indices\n";
    return EXIT_FAILURE;
  }

  const auto n = static_cast<Index>(std::atoll(c[1]));
  for (std::size_t i=0; i<N.size(); i++)  N[i] = static_cast<Index>(std::atoll(c[2 + i]));

  for (Index i=0; i<n; i++) {
    std::cout << N[0] << " input test_sum\n" << test_sum(N[0]) << '\n';
    std::cout << N[1] << " input test_dot\n"  << test_dot(N[1]) << '\n';
    std::cout << N[2] << " input test_vec_mult\n"  << test_vec_mult(N[2]) << '\n';
    std::cout << N[3] << " input test_mat_multiply\n"  << test_mat_multiply(N[3]) << '\n';
    std::cout << N[4] << " input test_elementary_ops_vec\n"  << test_elementary_ops_vec(N[4]) << '\n';
    std::cout << N[5] << " input test_elementary_ops_mat\n"  << test_elementary_ops_mat(N[5]) << '\n';
    std::cout << N[6] << " input test_ops_vec\n"  << test_ops_vec(N[6]) << '\n';
    std::cout << N[7] << " input test_ops_mat\n"  << test_ops_mat(N[7]) << '\n';
  }
}

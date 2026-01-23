#include <rng.h>
#include <rtepack.h>

#include <iomanip>
#include <iostream>

#include "artstime.h"
#include "configtypes.h"
#include "lin_alg.h"
#include "rtepack_transmission.h"

void test_expm() {
  const Numeric A = 0.01;
  for (Size i = 0; i < 100; i++) {
    auto rng = RandomNumberGenerator{}.get(0.0, A);
    auto rng2 = RandomNumberGenerator{}.get(-A, A);

    const Propmat k{rng(), rng2(), rng2(), rng2(), rng2(), rng2(), rng2()};
    std::println("{:B,}", k);

    Matrix k_expm = to_matrix(k);
    k_expm *= -1;
    Matrix t_expm(4, 4);

    matrix_exp(t_expm, k_expm, 80);
    auto m = rtepack::tran(k, k, 1)();
    auto mmat = Matrix(m);
    mmat -= t_expm;

    std::println("{:B,}\n{:B,}\n{:B,}\n", t_expm, Matrix(m), mmat);
  }
}

void test_dexpm() {
  constexpr Numeric A = 0.1;
  auto rng = RandomNumberGenerator{}.get(0.0, A);
  auto rng2 = RandomNumberGenerator{}.get(-A, A);
  auto drng = RandomNumberGenerator{}.get(0.0, A);
  auto drng2 = RandomNumberGenerator{}.get(-A, A);

  const auto gen = [a = rng(),
                    da = drng(),
                    b = rng2(),
                    db = drng2(),
                    c = rng2(),
                    dc = drng2(),
                    d = rng2(),
                    dd = drng2(),
                    u = rng2(),
                    du = drng2(),
                    v = rng2(),
                    dv = drng2(),
                    w = rng2(),
                    dw = drng2()](Numeric x, bool deriv) -> Propmat {
    return deriv ? Propmat{da, db, dc, dd, du, dv, dw}
                 : Propmat{a + da * x,
                           b + db * x,
                           c + dc * x,
                           d + dd * x,
                           u + du * x,
                           v + dv * x,
                           w + dw * x};
  };

  const Propmat k = gen(0.0, false);
  const PropmatVector dk(1, gen(0.0, true));
  const Vector dr{0};
  Muelmat t{};
  MuelmatVector dt(1, Muelmat{});

  const Numeric x = 1e-9;
  const Propmat k2 = gen(x, false);
  const PropmatVector dk2{};
  const Vector dr2{};
  Muelmat t2{};
  MuelmatVector dt2{};

  Muelmat m = (t2 - t);
  for (Size i = 0; i < 4; i++) {
    for (Size j = 0; j < 4; j++) {
      m[i][j] /= dt[0][i][j] * x * 2.0;
    }
  }
  std::cout << std::format("{}", dt) << '\n';
  std::cout << std::format("{}", m) << '\n';
}

void test_inv() {
  constexpr Numeric A = 0.1;
  auto rng = RandomNumberGenerator{}.get(0.0, A);
  auto rng2 = RandomNumberGenerator{}.get(-A, A);
  Muelmat m;
  const Propmat k{rng(), rng2(), rng2(), rng2(), rng2(), rng2(), rng2()};

  Matrix k_inv = Vector{k.A(),
                        k.B(),
                        k.C(),
                        k.D(),
                        k.B(),
                        k.A(),
                        k.U(),
                        k.V(),
                        k.C(),
                        -k.U(),
                        k.A(),
                        k.W(),
                        k.D(),
                        -k.V(),
                        -k.W(),
                        k.A()}
                     .reshape(4, 4);
  Matrix inv_k(4, 4);

  constexpr Size N = 1000000;
  Numeric sumup{};
  {
    DebugTime t("old inv");
    for (Size i = 0; i < N; i++) {
      inv(inv_k, k_inv);
      sumup += inv_k[0][0];
    }
  }
  {
    DebugTime t("inv");
    for (Size i = 0; i < N; i++) {
      m = inv(k);
      sumup += m[0][0];
    }
  }
  std::cout << sumup << '\n';

  for (Size i = 0; i < 4; i++) {
    for (Size j = 0; j < 4; j++) {
      inv_k[i][j] /= m[i][j];
    }
  }
  std::print(std::cout, "{}\n", inv_k);
}

int main() {
  test_expm();
  test_dexpm();
  test_inv();
  return 0;
}

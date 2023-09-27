#include <rng.h>
#include <rtepack.h>
#include <iomanip>
#include <iostream>

#include "artstime.h"
#include "lin_alg.h"

int main() {
  const Numeric A = 0.1;
  for (Size i=0; i<100; i++) {
  auto rng = RandomNumberGenerator{}.get(0.0, A);
  auto rng2 = RandomNumberGenerator{}.get(-A, A);

  Muelmat t_test;
  MuelmatVector dt;

  const Propmat k{rng(), rng2(), rng2(), rng2(), rng2(), rng2(), rng2()};
  for (auto&a:k.data) std::cout << -a << ',';
  std::cout << '\n';
  const PropmatVector dk;

  const Vector dr;

  Matrix k_expm = Vector{k.A(),
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
  k_expm *= -1;
  Matrix t_expm(4, 4);
  Size N = 100'000'000;
  N = 1;

  {
    //DebugTime expm{"expm"};
    for (Size i=0; i<N; i++)
      matrix_exp(t_expm, k_expm, 80);
  }
  {
    //DebugTime expm{"cheig"};
    for (Size i=0; i<N; i++)
      rtepack::two_level_exp_test(t_test, dt, dt, k, k, dk, dk, 1.0, dr, dr);
  }

  Matrix t_test_diff = Vector{t_test.data}
                      .reshape(4, 4);

                    t_test_diff -= t_expm;
                    t_test_diff /= t_expm;
                    std::cout << t_expm[0][0] << ' ';
                    std::cout << t_test[0][0] << '\n';
                    std::cout  << '\n';

  }
}
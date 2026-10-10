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
    auto rng  = RandomNumberGenerator{}.get(0.0, A);
    auto rng2 = RandomNumberGenerator{}.get(-A, A);

    const Propmat k{rng(), rng2(), rng2(), rng2(), rng2(), rng2(), rng2()};
    std::println("{:B,}", k);

    Matrix k_expm  = to_matrix(k);
    k_expm        *= -1;
    Matrix t_expm(4, 4);

    matrix_exp(t_expm, k_expm, 80);
    auto m     = rtepack::tran(k, k, 1)();
    auto mmat  = Matrix(m);
    mmat      -= t_expm;

    std::println("{:B,}\n{:B,}\n{:B,}\n", t_expm, Matrix(m), mmat);
  }
}

void test_dexpm() {
  constexpr Numeric A     = 0.1;
  auto              rng   = RandomNumberGenerator{}.get(0.0, A);
  auto              rng2  = RandomNumberGenerator{}.get(-A, A);
  auto              drng  = RandomNumberGenerator{}.get(0.0, A);
  auto              drng2 = RandomNumberGenerator{}.get(-A, A);

  const auto gen = [a  = rng(),
                    da = drng(),
                    b  = rng2(),
                    db = drng2(),
                    c  = rng2(),
                    dc = drng2(),
                    d  = rng2(),
                    dd = drng2(),
                    u  = rng2(),
                    du = drng2(),
                    v  = rng2(),
                    dv = drng2(),
                    w  = rng2(),
                    dw = drng2()](Numeric x, bool deriv) -> Propmat {
    return deriv ? Propmat{da, db, dc, dd, du, dv, dw}
                 : Propmat{a + da * x, b + db * x, c + dc * x, d + dd * x, u + du * x, v + dv * x, w + dw * x};
  };

  const PropmatVector dk(1, gen(0.0, true));
  const Vector        dr{0};
  Muelmat             t{};
  MuelmatVector       dt(1, Muelmat{});

  const Numeric       x = 1e-9;
  const PropmatVector dk2{};
  const Vector        dr2{};
  Muelmat             t2{};
  MuelmatVector       dt2{};

  Muelmat m = (t2 - t);
  for (Size i = 0; i < 4; i++) {
    for (Size j = 0; j < 4; j++) { m[i][j] /= dt[0][i][j] * x * 2.0; }
  }
  std::cout << std::format("{}", dt) << '\n';
  std::cout << std::format("{}", m) << '\n';
}

void test_inv() {
  constexpr Numeric A    = 0.1;
  auto              rng  = RandomNumberGenerator{}.get(0.0, A);
  auto              rng2 = RandomNumberGenerator{}.get(-A, A);
  Muelmat           m;
  const Propmat     k{rng(), rng2(), rng2(), rng2(), rng2(), rng2(), rng2()};

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
  Numeric        sumup{};
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
      m      = inv(k);
      sumup += m[0][0];
    }
  }
  std::cout << sumup << '\n';

  for (Size i = 0; i < 4; i++) {
    for (Size j = 0; j < 4; j++) { inv_k[i][j] /= m[i][j]; }
  }
  std::print(std::cout, "{}\n", inv_k);
}

//! The compact scattering-plane Mueller matrix: its layout, and products that are those of the full matrix
void test_compact_planar_muelmat() {
  const rtepack::compact_planar_muelmat f{1.0, -0.3, 0.8, 0.6, 0.2, 0.5};
  const rtepack::muelmat                F = f.expand();
  const rtepack::muelmat expected{1.0, -0.3, 0.0, 0.0, -0.3, 0.8, 0.0, 0.0, 0.0, 0.0, 0.6, 0.2, 0.0, 0.0, -0.2, 0.5};
  const rtepack::muelmat L  = rtepack::stokes_rotation(std::cos(0.8), std::sin(0.8));
  const rtepack::muelmat LF = L * f, FL = f * L, FF = f * f, LF_ = L * F, FL_ = F * L, FF_ = F * F;
  for (Index i = 0; i < 4; i++) {
    for (Index j = 0; j < 4; j++) {
      const bool layout = F[i, j] == expected[i, j];
      // Equal up to floating-point contraction of the inlined products
      const auto close = [](Numeric a, Numeric b) {
        return std::abs(a - b) <= 4 * std::numeric_limits<Numeric>::epsilon();
      };
      const bool products = close(LF[i, j], LF_[i, j]) and close(FL[i, j], FL_[i, j]) and close(FF[i, j], FF_[i, j]);
      ARTS_USER_ERROR_IF(not layout, "compact_planar_muelmat expands wrongly at [{}, {}]", i, j)
      ARTS_USER_ERROR_IF(
          not products, "compact_planar_muelmat products must be those of the full matrix at [{}, {}]", i, j)
    }
  }
  ARTS_USER_ERROR_IF(
      f.F11() != 1.0 or f.F12() != -0.3 or f.F22() != 0.8 or f.F33() != 0.6 or f.F34() != 0.2 or f.F44() != 0.5,
      "compact_planar_muelmat accessors are in the wrong order")
  std::cout << "compact_planar_muelmat: layout, accessors and products OK\n";
}

//! rotated is the product of the rotations and the matrix, for a full and a compact one alike, and mirror is
//! D A D with D = diag(1, 1, -1, -1)
void test_rotated_and_mirror() {
  const rtepack::compact_planar_muelmat f{1.0, -0.3, 0.8, 0.6, 0.2, 0.5};
  const rtepack::muelmat A{0.9, -0.4, 0.3, 0.1, 0.2, 0.7, -0.6, 0.5, -0.8, 0.35, 0.45, -0.25, 0.15, -0.55, 0.65, 0.95};
  const Numeric          c1 = std::cos(0.8), s1 = std::sin(0.8), c2 = std::cos(-2.1), s2 = std::sin(-2.1);
  const rtepack::muelmat Z  = rtepack::rotated(f, c1, s1, c2, s2);
  const rtepack::muelmat Z_ = rtepack::stokes_rotation(c2, s2) * f * rtepack::stokes_rotation(c1, s1);
  const rtepack::muelmat Zf = rtepack::rotated(f.expand(), c1, s1, c2, s2);
  const rtepack::muelmat R  = rtepack::rotated(A, c1, s1, c2, s2);
  const rtepack::muelmat R_ = rtepack::stokes_rotation(c2, s2) * A * rtepack::stokes_rotation(c1, s1);
  const rtepack::muelmat D{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, -1.0};
  const rtepack::muelmat M = rtepack::mirror(Z), M_ = D * Z * D;
  for (Index i = 0; i < 4; i++) {
    for (Index j = 0; j < 4; j++) {
      // Equal up to the rounding of the products
      const auto close = [](Numeric a, Numeric b) {
        return std::abs(a - b) <= 4 * std::numeric_limits<Numeric>::epsilon();
      };
      ARTS_USER_ERROR_IF(not close(Z[i, j], Z_[i, j]) or not close(Zf[i, j], Z_[i, j]),
                         "rotated must be the product of the rotations and F at [{}, {}]",
                         i,
                         j)
      ARTS_USER_ERROR_IF(
          not close(R[i, j], R_[i, j]), "rotated must be the product of the rotations and A at [{}, {}]", i, j)
      ARTS_USER_ERROR_IF((M[i, j] != M_[i, j]), "mirror must be D A D at [{}, {}]", i, j)
    }
  }
  std::cout << "rotated and mirror: the products OK\n";
}

int main() {
  test_expm();
  test_dexpm();
  test_inv();
  test_compact_planar_muelmat();
  test_rotated_and_mirror();
  return 0;
}

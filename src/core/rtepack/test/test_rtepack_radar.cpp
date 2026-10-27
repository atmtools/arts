#include <rtepack.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
using namespace rtepack;

constexpr Size np = 4;
constexpr Size nf = 1;
constexpr Size nq = 3;

//! Polarized levels whose neighbouring propagation matrices do not commute
struct Levels {
  std::vector<propmat_vector> K_out, K_ret;
  std::vector<propmat_matrix> dK_out, dK_ret;
  std::vector<muelmat_vector> Z;
  std::vector<muelmat_matrix> dZ;
  Vector                      r{300.0, 250.0, 400.0};

  Levels()
      : K_out(np, propmat_vector(nf)),
        K_ret(np, propmat_vector(nf)),
        dK_out(np, propmat_matrix(nq, nf)),
        dK_ret(np, propmat_matrix(nq, nf)),
        Z(np, muelmat_vector(nf)),
        dZ(np, muelmat_matrix(nq, nf)) {
    for (Size i = 0; i < np; i++) {
      const Numeric s = static_cast<Numeric>(i + 1);
      K_out[i][0]     = propmat{1e-3 * s, 2e-4, -1e-4 * s, 5e-5, 1e-4, -2e-4 / s, 3e-5 * s};
      K_ret[i][0]     = propmat{1.1e-3 * s, -1e-4 * s, 2e-4, -5e-5, 2e-4 / s, 1e-4, -3e-5};
      Z[i][0] =
          1e-5 *
          muelmat{1.0 + s, 0.2, 0.1 * s, 0.0, 0.2, 0.9, -0.1, 0.05 * s, -0.1, 0.1, 0.8 * s, 0.2, 0.0, -0.05, 0.1, 0.7};

      // The three derivative directions perturb K_out, K_ret and Z in turn
      dK_out[i][0, 0] = propmat{1e-4, -3e-5 * s, 2e-5, 1e-5, -2e-5, 4e-5, 1e-5};
      dK_ret[i][1, 0] = propmat{-2e-4, 1e-5, -4e-5 * s, 2e-5, 3e-5, -1e-5, 2e-5};
      dZ[i][2, 0]     = 1e-6 * muelmat{1.0, 0.3, 0.0, -0.2, 0.1, 1.0, 0.2, 0.0, 0.0, -0.1, 1.0, 0.1, 0.2, 0.0, 0.1, s};
    }
  }

  void perturb(Size level, Size q, Numeric h) {
    switch (q) {
      case 0:  K_out[level][0] += h * dK_out[level][0, 0]; break;
      case 1:  K_ret[level][0] += h * dK_ret[level][1, 0]; break;
      default: Z[level][0] += h * dZ[level][2, 0]; break;
    }
  }

  void solve(stokvec_matrix&                  I,
             stokvec_tensor4&                 dI,
             const std::vector<radar_sample>& samples,
             Size                             ngate,
             TransmittanceOption              option) const {
    const stokvec transmitted{1.0, 0.3, -0.2, 0.1};
    I.resize(ngate, nf);
    dI.resize(ngate, np, nq, nf);
    radar_single_scattering(I, dI, transmitted, K_out, dK_out, K_ret, dK_ret, Z, dZ, r, samples, option);
  }
};

const char* name(TransmittanceOption option) { return option == TransmittanceOption::magop ? "magop" : "constant"; }

bool test_derivatives(TransmittanceOption option) {
  const std::vector<radar_sample> samples{{.layer = 1, .x = 0.3, .weight = 0.4, .gate = 0},
                                          {.layer = 2, .x = 0.5, .weight = 0.6, .gate = 0},
                                          {.layer = 2, .x = 0.9, .weight = 0.2, .gate = 1},
                                          {.layer = 3, .x = 0.25, .weight = 0.8, .gate = 1}};
  constexpr Size                  ngate = 2;

  const Levels    levels;
  stokvec_matrix  I;
  stokvec_tensor4 dI;
  levels.solve(I, dI, samples, ngate, option);

  // The derivatives are linear in the direction, so central differences are
  // accurate to the second order in h.
  constexpr Numeric h = 1e-3;
  stokvec_matrix    Ip, Im;
  stokvec_tensor4   dummy;
  bool              ok = true;
  for (Size level = 0; level < np; level++) {
    for (Size q = 0; q < nq; q++) {
      Levels plus = levels, minus = levels;
      plus.perturb(level, q, h);
      minus.perturb(level, q, -h);
      plus.solve(Ip, dummy, samples, ngate, option);
      minus.solve(Im, dummy, samples, ngate, option);

      for (Size gate = 0; gate < ngate; gate++) {
        for (Size s = 0; s < 4; s++) {
          const Numeric numeric  = (Ip[gate, 0][s] - Im[gate, 0][s]) / (2.0 * h);
          const Numeric analytic = dI[gate, level, q, 0][s];
          const Numeric scale    = std::max(std::abs(I[gate, 0][s]), 1e-12);
          if (not std::isfinite(analytic) or std::abs(analytic - numeric) > 1e-7 * scale) {
            std::cerr << "Radar derivative for " << name(option) << ", level " << level << ", direction " << q
                      << ", gate " << gate << ", Stokes " << s << ": analytic " << analytic << ", numeric " << numeric
                      << '\n';
            ok = false;
          }
        }
      }
    }
  }
  return ok;
}

bool test_layer_continuity(TransmittanceOption option) {
  // The far end of a layer and the near end of the next layer are the same
  // point, so the partial-layer transmissions must agree with the full ones.
  const std::vector<radar_sample> samples{{.layer = 1, .x = 1.0, .weight = 1.0, .gate = 0},
                                          {.layer = 2, .x = 0.0, .weight = 1.0, .gate = 1},
                                          {.layer = 2, .x = 1.0, .weight = 1.0, .gate = 2},
                                          {.layer = 3, .x = 0.0, .weight = 1.0, .gate = 3}};

  const Levels    levels;
  stokvec_matrix  I;
  stokvec_tensor4 dI;
  levels.solve(I, dI, samples, 4, option);

  bool ok = true;
  for (Size gate = 0; gate < 4; gate += 2) {
    for (Size s = 0; s < 4; s++) {
      const Numeric far  = I[gate, 0][s];
      const Numeric near = I[gate + 1, 0][s];
      if (std::abs(far - near) > 1e-12 * std::max(std::abs(far), 1e-12)) {
        std::cerr << "Radar return for " << name(option) << " is discontinuous at level " << gate / 2 + 1 << ", Stokes "
                  << s << ": " << far << " vs " << near << '\n';
        ok = false;
      }
    }
  }

  for (Size gate = 0; gate < 4; gate++)
    for (Size level = 0; level < np; level++)
      for (Size q = 0; q < nq; q++)
        for (Size s = 0; s < 4; s++)
          if (not std::isfinite(dI[gate, level, q, 0][s])) {
            std::cerr << "Radar derivative for " << name(option) << " is not finite at a level\n";
            return false;
          }
  return ok;
}

bool test_magnus_is_used() {
  const std::vector<radar_sample> samples{{.layer = 3, .x = 0.5, .weight = 1.0, .gate = 0}};

  const Levels    levels;
  stokvec_matrix  constant, magnus;
  stokvec_tensor4 dI;
  levels.solve(constant, dI, samples, 1, TransmittanceOption::constant);
  levels.solve(magnus, dI, samples, 1, TransmittanceOption::magop);

  if (constant[0, 0] == magnus[0, 0]) {
    std::cerr << "The Magnus operator does not change the return of non-commuting layers\n";
    return false;
  }
  return true;
}
}  // namespace

int main() {
  bool ok = test_magnus_is_used();
  for (const auto option : {TransmittanceOption::constant, TransmittanceOption::magop}) {
    ok = test_derivatives(option) and ok;
    ok = test_layer_continuity(option) and ok;
  }
  return ok ? 0 : 1;
}

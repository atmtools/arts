#include <arts_constexpr_math.h>
#include <lin_alg.h>

#include <array>
#include <cmath>

#include "lbl_zeeman.h"

namespace lbl::zeeman {
namespace {
using Math::pow2;
using std::sqrt;

/*** Hyperfine Hamiltonian, matrix elements from Gerry et al. (1989) ***/

Numeric get_x(Numeric j) { return pow2(j + 0.5); }

Numeric get_R(Numeric f, Numeric i, Numeric j) { return f * (f + 1) - i * (i + 1) - j * (j + 1); }

Numeric get_W(Numeric R, Numeric i, Numeric j) {
  if (i == 0.5 or j == 0.5) return 0.0;
  const Numeric lj = j * (j + 1);
  return (3 * R * (R + 1) - 4 * i * (i + 1) * lj) / (8 * i * (2 * i - 1) * lj * (2 * j - 1) * (2 * j + 3));
}

Numeric get_P(Numeric f, Numeric i, Numeric j) {
  return sqrt((f + i + j + 1) * (i + j - f) * (f + j - i) * (f + i - j + 1)) / (4 * j);
}

Numeric get_T(Numeric R, Numeric i, Numeric j) {
  if (i == 0.5) return 0.0;
  return (R + j + 1) / (4 * i * (2 * i - 1) * (j - 1) * (j + 1));
}

Numeric get_Z(Numeric P1, Numeric P2, Numeric i, Numeric j) { return P1 * P2 / (i * (2 * i - 1) * sqrt(2 * j - 1)); }

/** Hyperfine Hamiltonian elements of one F
 *
 * The element <J - dJ, o1| H |J, o2>, with o = 0 for |Omega| = 1/2 and
 * o = 1 for |Omega| = 3/2, and e the e/f label of the J - dJ level.
 */
struct hyperfine_hamiltonian {
  const hund_a_hyperfine& c;
  Numeric                 f, i;

  [[nodiscard]] Numeric h_f() const { return c.a - 0.5 * (c.b_F + c.c); }
  [[nodiscard]] Numeric h_f_alt() const { return c.a + 0.5 * (c.b_F + c.c); }

  [[nodiscard]] Numeric diag_om1(Numeric j, Numeric ef) const {
    const Numeric x  = get_x(j);
    const Numeric R  = get_R(f, i, j);
    const Numeric lj = j * (j + 1);
    const Numeric W  = get_W(R, i, j);
    const Numeric pq = c.p + 2 * c.q;

    return -0.5 * c.A + (c.B - 0.5 * c.A_D) * x - c.D * (x * x + x - 1) + h_f() * R / (4 * lj) + 0.5 * c.c_I * R +
           ef * (j + 0.5) * (0.5 * pq + c.d * R / (4 * lj)) + c.eQq0 * W * (0.75 - lj);
  }

  [[nodiscard]] Numeric diag_om2(Numeric j) const {
    const Numeric x  = get_x(j);
    const Numeric R  = get_R(f, i, j);
    const Numeric lj = j * (j + 1);
    const Numeric W  = get_W(R, i, j);

    return 0.5 * c.A + (c.B + 0.5 * c.A_D) * (x - 2) - c.D * (x * x - 3 * x + 3) + 3 * h_f_alt() * R / (4 * lj) +
           0.5 * c.c_I * R + c.eQq0 * W * (27.0 / 4.0 - lj);
  }

  [[nodiscard]] Numeric diag_om21(Numeric j, Numeric ef) const {
    const Numeric x  = get_x(j);
    const Numeric R  = get_R(f, i, j);
    const Numeric lj = j * (j + 1);
    const Numeric W  = get_W(R, i, j);

    return -(c.B - 2 * c.D * (x - 1) - 0.5 * c.gamma - c.b_F * R / (4 * lj) + ef * 0.5 * (j + 0.5) * (c.q + c.eQq2 * W)) *
           sqrt(x - 1);
  }

  [[nodiscard]] Numeric minoneJ_om1(Numeric j, Numeric ef) const {
    const Numeric R = get_R(f, i, j);
    const Numeric T = get_T(R, i, j);
    const Numeric P = get_P(f, i, j);

    return -(h_f() + ef * 0.5 * c.d + 1.5 * c.eQq0 * T) * P;
  }

  [[nodiscard]] Numeric minoneJ_om2(Numeric j) const {
    const Numeric P = get_P(f, i, j);

    return -P * (h_f_alt() + 4.5 * c.eQq0) * sqrt((4 * j * j - 9) / (4 * j * j - 1));
  }

  [[nodiscard]] Numeric minoneJ_om12(Numeric j, Numeric ef) const {
    const Numeric R = get_R(f, i, j);
    const Numeric T = get_T(R, i, j);
    const Numeric P = get_P(f, i, j);

    return (c.b_F + ef * c.eQq2 * T * (j - 0.5)) * 0.5 * P * sqrt((2 * j - 3) * (2 * j + 1)) * (2 * j + 3) / (2 * j - 1);
  }

  [[nodiscard]] Numeric minoneJ_om21(Numeric j, Numeric ef) const {
    const Numeric R = get_R(f, i, j);
    const Numeric T = get_T(R, i, j);
    const Numeric P = get_P(f, i, j);

    return -(c.b_F - ef * c.eQq2 * T * (j + 0.5)) * 0.5 * P * sqrt((2 * j - 3) * (2 * j + 1));
  }

  [[nodiscard]] Numeric Z(Numeric j) const { return get_Z(get_P(f, i, j), get_P(f, i, j - 1), i, j); }

  [[nodiscard]] Numeric mintwoJ_om(Numeric j, Numeric om) const {
    return 3 * c.eQq0 * Z(j) * (j * j - om * om) * (pow2(j - 1) - om * om) /
           sqrt((2 * j - 3) * (2 * j - 1) * (2 * j + 1));
  }

  [[nodiscard]] Numeric mintwoJ_om12(Numeric j, Numeric ef) const { return -ef * c.eQq2 * Z(j) * sqrt(2 * j + 3) / 8; }

  [[nodiscard]] Numeric mintwoJ_om21(Numeric j, Numeric ef) const { return -ef * c.eQq2 * Z(j) * sqrt(2 * j - 5) / 8; }

  [[nodiscard]] Numeric operator()(Numeric j, Index dj, Index o1, Index o2, bool e) const {
    const Numeric ef = e ? -1.0 : 1.0;

    switch (dj) {
      case 0:
        if (o1 != o2) return diag_om21(j, ef);
        return o1 == 0 ? diag_om1(j, ef) : diag_om2(j);
      case 1:
        if (o1 == 0) return o2 == 0 ? minoneJ_om1(j, ef) : minoneJ_om12(j, ef);
        return o2 == 0 ? minoneJ_om21(j, ef) : minoneJ_om2(j);
      case 2:
        if (o1 == 0) return o2 == 0 ? mintwoJ_om(j, 0.5) : mintwoJ_om12(j, ef);
        return o2 == 0 ? mintwoJ_om21(j, ef) : mintwoJ_om(j, 1.5);
      default: return 0.0;
    }
  }
};

/*** Zeeman Hamiltonian, matrix elements from Carrington et al. (1967) and Brown et al. (2003) ***/

/** The signed Hund's case (a) states of one J: Sigma, Lambda, and Omega = Lambda + Sigma */
constexpr std::array<Numeric, 4> signed_sigma{-0.5, 0.5, -0.5, 0.5};
constexpr std::array<Numeric, 4> signed_lambda{-1.0, -1.0, 1.0, 1.0};

Numeric get_Pzee(Numeric f, Numeric i, Numeric j) { return sqrt((f - i + j) * (f + j + i + 1)); }

Numeric get_Qzee(Numeric f, Numeric i, Numeric j) { return sqrt((j + i - f) * (f - j + i + 1)); }

Numeric get_Vzee(Numeric sig, Index pm) {
  const Numeric d_s   = std::abs(sig) * (std::abs(sig) + 1);
  const Numeric d_sig = pm == 1 ? sig * (sig + 1) : sig * (sig - 1);
  return sqrt(d_s - d_sig);
}

/** Zeeman Hamiltonian elements of one F in the signed basis, in units of the Bohr magneton
 *
 * The orbital and spin g-factors include the rotational g-factor.
 */
struct zeeman_hamiltonian {
  Numeric f, i;
  Numeric gl, ge, gi, gr;

  //! <J, o| H |J, o>
  [[nodiscard]] Numeric diag(Numeric j, Index o) const {
    const Numeric sig = signed_sigma[o], lam = signed_lambda[o], om = lam + sig;
    const Numeric lj = j * (j + 1);
    const Numeric lf = f * (f + 1);
    const Numeric li = i * (i + 1);
    const Numeric U  = lf - li + lj;
    const Numeric T  = lf + li - lj;

    return ((gl * lam + ge * sig) * om * U - lj * T * gi - lj * U * gr) / (2 * lj * lf);
  }

  //! <J, o| H |J + 1, o>
  [[nodiscard]] Numeric plusJ(Numeric j, Index o) const {
    const Numeric sig = signed_sigma[o], lam = signed_lambda[o], om = lam + sig;
    const Numeric P   = get_Pzee(f + 1, i, j);
    const Numeric Q   = get_Qzee(f - 1, i, j);
    const Numeric fac = sqrt((j + om + 1) * (j - om + 1));
    const Numeric lf  = f * (f + 1);

    return -(gl * lam + ge * sig) * fac * P * Q / (2 * lf * (j + 1) * sqrt((2 * j + 1) * (2 * j + 3)));
  }

  //! <J, o| H |J, o + 1>, Sigma of o is -1/2
  [[nodiscard]] Numeric off_om(Numeric j, Index o) const {
    const Numeric sig = signed_sigma[o], lam = signed_lambda[o], om = lam + sig;
    const Numeric V   = get_Vzee(sig, 1);
    const Numeric lf  = f * (f + 1);
    const Numeric li  = i * (i + 1);
    const Numeric lj  = j * (j + 1);
    const Numeric U   = lf - li + lj;

    return ge * sqrt((j - om) * (j + om + 1)) * V * U / (4 * lf * lj);
  }

  //! <J, o| H |J - 1, o + 1> if Sigma of o is -1/2, else <J, o| H |J - 1, o - 1>
  [[nodiscard]] Numeric minJ_off_om(Numeric j, Index o) const {
    const Numeric sig  = signed_sigma[o], lam = signed_lambda[o], om = lam + sig;
    const Index   pm   = sig < 0 ? 1 : -1;
    const Numeric V    = get_Vzee(sig, pm);
    const Numeric P    = get_Pzee(f, i, j);
    const Numeric Q    = get_Qzee(f, i, j);
    const Numeric lf   = f * (f + 1);
    const Numeric noem = 4 * lf * j * sqrt(4 * j * j - 1);

    if (pm == 1) return -ge * sqrt((j - om) * (j - om - 1)) * V * P * Q / noem;
    return ge * sqrt((j + om) * (j + om - 1)) * V * P * Q / noem;
  }

  //! The signed-basis block <J| H |J>
  [[nodiscard]] Matrix block(Numeric j) const {
    Matrix G(4, 4, 0.0);
    for (Index o = 0; o < 4; o++) G[o, o] = diag(j, o);
    for (Index o : {0, 2}) G[o, o + 1] = G[o + 1, o] = off_om(j, o);
    return G;
  }

  //! The signed-basis block <J| H |J + 1>
  [[nodiscard]] Matrix block_plus(Numeric j) const {
    Matrix G(4, 4, 0.0);
    for (Index o = 0; o < 4; o++) G[o, o] = plusJ(j, o);
    for (Index o : {0, 2}) G[o + 1, o] = minJ_off_om(j + 1, o);
    for (Index o : {1, 3}) G[o - 1, o] = minJ_off_om(j + 1, o);
    return G;
  }
};

/** Transform a signed-basis block of the Zeeman Hamiltonian to the parity basis
 *
 * The rows of the block belong to J.  Returns the |Omega| = 1/2, 3/2 block
 * of the e levels if e, else of the f levels.
 */
Matrix parity_block(const Matrix& G, Numeric j, bool e) {
  const Numeric s   = 1.0 / sqrt(2.0);
  const Numeric fas = std::pow(-1.0, j - 0.5);

  Matrix U(4, 4, 0.0);
  U[1, 0] = s;
  U[2, 0] = s * fas;
  U[0, 1] = s;
  U[3, 1] = s * fas;
  U[1, 2] = s;
  U[2, 2] = -s * fas;
  U[0, 3] = s;
  U[3, 3] = -s * fas;

  Matrix GU(4, 4), H(4, 4);
  mult(GU, G, U);
  mult(H, transpose(U), GU);

  const Index k = e ? 0 : 2;
  Matrix      out(2, 2);
  for (Index a = 0; a < 2; a++)
    for (Index b = 0; b < 2; b++) out[a, b] = H[k + a, k + b];
  return out;
}
}  // namespace

std::optional<Numeric> hund_a_hyperfine::g(Rational F, Rational J, Rational Omega, bool e) const {
  const Rational J0 = abs(F - I);

  if (Omega != Rational(1, 2) and Omega != Rational(3, 2)) return std::nullopt;
  if (J < Omega) return std::nullopt;
  if (J < J0 or J > F + I or not iseven(Rational{2} * (J - J0))) return std::nullopt;

  const Numeric f  = F.toNumeric();
  const Numeric i  = I.toNumeric();
  const Numeric j0 = J0.toNumeric();
  const Index   nj = (Rational{2} * minr(F, I) + Rational{1}).toIndex();
  const Index   n  = 2 * nj;

  // The e/f labels alternate with J within a parity block
  const bool e0   = e == ((J - J0).toIndex() % 2 == 0);
  const auto is_e = [e0](Index jq) { return e0 == (jq % 2 == 0); };

  // Hyperfine Hamiltonian in the basis |J0 + jq, o> at index 2 jq + o
  const hyperfine_hamiltonian hyp{.c = *this, .f = f, .i = i};
  Matrix                      H(n, n, 0.0);
  for (Index jq1 = 0; jq1 < nj; jq1++) {
    for (Index o1 = 0; o1 < 2; o1++) {
      for (Index jq2 = jq1; jq2 < std::min(jq1 + 3, nj); jq2++) {
        for (Index o2 = jq2 == jq1 ? o1 : 0; o2 < 2; o2++) {
          H[2 * jq1 + o1, 2 * jq2 + o2] = H[2 * jq2 + o2, 2 * jq1 + o1] =
              hyp(j0 + static_cast<Numeric>(jq2), jq2 - jq1, o1, o2, is_e(jq1));
        }
      }
    }
  }

  Matrix U(n, n);
  Vector E(n);
  diagonalize_symmetric(U, E, H);

  // Zeeman Hamiltonian in the same basis: the diagonal blocks <J|H|J>
  // and the element <J, 3/2| H |J + 1, 1/2> of the blocks <J|H|J + 1>.
  const zeeman_hamiltonian zee{.f = f, .i = i, .gl = gL + gR, .ge = gS + gR, .gi = gI, .gr = gR};
  std::vector<Matrix>      Gd(nj);
  Vector                   Go(nj - 1);
  for (Index jq = 0; jq < nj; jq++) {
    const Numeric j = j0 + static_cast<Numeric>(jq);
    Gd[jq]          = parity_block(zee.block(j), j, is_e(jq));
    if (jq + 1 < nj) Go[jq] = parity_block(zee.block_plus(j), j, is_e(jq))[1, 0];
  }

  Index   nfound = 0;
  Numeric g      = 0.0;
  for (Index k = 0; k < n; k++) {
    Numeric lj_exp = 0.0, om_exp = 0.0;
    for (Index jq = 0; jq < nj; jq++) {
      const Numeric j = j0 + static_cast<Numeric>(jq);
      for (Index o = 0; o < 2; o++) {
        const Numeric u2  = pow2(U[2 * jq + o, k]);
        lj_exp           += u2 * j * (j + 1);
        om_exp           += u2 * (0.5 + static_cast<Numeric>(o));
      }
    }
    const Numeric j_exp = -0.5 + sqrt(0.25 + lj_exp);

    if (std::round(2 * j_exp) != 2 * J.toNumeric() or (om_exp > 1) != (Omega == Rational(3, 2))) continue;

    nfound++;
    g = 0.0;
    for (Index jq = 0; jq < nj; jq++) {
      for (Index a = 0; a < 2; a++)
        for (Index b = 0; b < 2; b++) g += U[2 * jq + a, k] * Gd[jq][a, b] * U[2 * jq + b, k];
      if (jq + 1 < nj) g += U[2 * jq + 1, k] * Go[jq] * U[2 * jq + 2, k];
    }
  }

  if (nfound != 1) return std::nullopt;
  return g;
}
}  // namespace lbl::zeeman

/* Polarized delta-M with the Nakajima-Tanaka corrections (TMS and IMS)
   against VDISORT itself without delta-M.

   Where the expected answer comes from.  One layer with a forward-peaked,
   polarizing phase matrix (a Legendre series of scattering-plane matrices,
   F12 and F34 included) under an unpolarized beam
   is solved by VDISORT without delta-M on enough streams to carry the
   whole series: the reference.  The same layer is then solved with
   classical delta-M on few streams, at user directions near and away from
   the forward peak, without and with the cached polarized TMS and IMS
   corrections.  The reference shares the transport solver, which the RT3
   and RT4 comparisons validate, but not the delta-M scaling, the removed
   peak or the corrections, which this tests in all four Stokes components.

   The delta-M split is DISORT's classical one (removed moments r_l = 1)
   applied to the whole matrix.  With N streams, f = chi_N of F11 and
   c_l = (2 l + 1) chi_l the series of the layer,
     c_l = (1 - f) c'_l + f r_l,
   with the transport series c'_l = (c_l - f (2 l + 1) id) / (1 - f) for
   l < N (zero beyond) and the removed peak r_l = (2 l + 1) id for l < N and
   c_l / f beyond (id the identity), so the split is exact for every l.
   The phase matrices, on the streams and into the user directions, are
   the Fourier modes of the laboratory-frame matrices of these series by
   vector geometry (lab-frame.h).

   The phase matrix must be regular at forward and backward scattering
   (F12 = F34 = 0 at both, F22 = F33 forward and F22 = -F33 backward),
   where the scattering plane is undefined: otherwise its laboratory-frame Q
   and U elements jump at the forward peak, and no discrete-ordinate
   solution converges there.  The layer's matrix is therefore built in angle
   space, with x = cos(Theta) and h the Henyey-Greenstein F11,
     F11 = F22 = h,  F33 = F44 = x h,  F12 = -0.3 (1 - x^2) h,
     F34 = 0.2 x (1 - x^2) h,
   and its series follows exactly from multiplying Legendre series by x.

   TMS subtracts the beam's single scattering as the solver carries it,
   from the transport beam phase matrices into the user directions
   (NFourier modes).  For F11 that is the whole transport matrix; the
   polarized elements also contain the rotation angles, which are not
   band-limited, and subtracting the full matrix instead made Q 2 to 4
   times worse than delta-M without corrections. */

#include <arts_constants.h>
#include <legendre.h>
#include <vdisort.h>

#include <algorithm>
#include <cmath>
#include <format>
#include <iostream>
#include <stdexcept>
#include <string>

#include "lab-frame.h"

namespace {
using vdisort_test::fourier_modes;
using vdisort_test::lab_frame;
using vdisort_test::legendre_matrix;

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

constexpr Numeric g = 0.75, omega = 0.9, mu0 = 0.5, phi0 = 0.0;
constexpr Index   degree = 160;  // g^160 = 1e-20

//! The Legendre coefficients of x f(x) from those of f: x P_l = ((l + 1) P_l+1 + l P_l-1) / (2 l + 1)
Vector times_x(const Vector& a) {
  const auto n = static_cast<Index>(a.size());
  Vector     b(n + 1, 0.0);
  for (Index l = 0; l < n; l++) {
    const auto k  = static_cast<Numeric>(2 * l + 1);
    b[l + 1]     += a[l] * static_cast<Numeric>(l + 1) / k;
    if (l > 0) b[l - 1] += a[l] * static_cast<Numeric>(l) / k;
  }
  return b;
}

//! The layer's series, {F11, F12, F22, F33, F34, F44} per degree, regular at forward and backward scattering
CompactPlanarMuelmatVector layer_series() {
  Vector h(degree + 1);
  for (Index l = 0; l <= degree; l++) h[l] = static_cast<Numeric>(2 * l + 1) * std::pow(g, l);
  const Vector xh = times_x(h), xxh = times_x(xh), xxxh = times_x(xxh);
  const auto   at = [](const Vector& v, Index l) { return l < static_cast<Index>(v.size()) ? v[l] : 0.0; };
  CompactPlanarMuelmatVector c(degree + 1);
  for (Index l = 0; l <= degree; l++) {
    const Numeric f11 = h[l], f33 = at(xh, l);
    const Numeric f12 = -0.3 * (h[l] - at(xxh, l)), f34 = 0.2 * (at(xh, l) - at(xxxh, l));
    c[l] = {f11, f12, f11, f33, f34, f33};
  }
  return c;
}

//! The identity in compact form times 2 l + 1
CompactPlanarMuelmat unit(Index l) { return static_cast<Numeric>(2 * l + 1) * CompactPlanarMuelmat::id(); }

struct split {
  Numeric                    f;
  CompactPlanarMuelmatVector transport, removed;
};

split delta_m_split(const CompactPlanarMuelmatVector& c, Index N) {
  const Numeric              f = c[N].F11() / static_cast<Numeric>(2 * N + 1);
  CompactPlanarMuelmatVector transport(N), removed(degree + 1);
  for (Index l = 0; l <= degree; l++) {
    removed[l] = l < N ? unit(l) : c[l] / f;
    if (l < N) transport[l] = (c[l] - f * unit(l)) / (1 - f);
  }
  return {.f = f, .transport = std::move(transport), .removed = std::move(removed)};
}

//! VDISORT's signed streams: N positive double-Gauss nodes, then their negatives
Vector streams(Index nquad) {
  const Index n = nquad / 2;
  Vector      mu(nquad), w(n);
  Legendre::PositiveDoubleGaussLegendre(mu[Range{0, n}], w);
  for (Index i = 0; i < n; i++) mu[n + i] = -mu[i];
  return mu;
}

//! The combined phase matrices of a series into the directions out, from the streams and from the beam
struct phases {
  vdisort::phase_matrix_data      diffuse;
  vdisort::beam_phase_matrix_data beam;
};

phases combined(const CompactPlanarMuelmatVector& series, const Vector& out, const Vector& mu, Index nfourier) {
  constexpr Index          nphi = 160;
  const auto               F    = legendre_matrix(series);
  const auto               d    = fourier_modes(F, out, mu, nfourier, nphi);
  const auto               b    = fourier_modes(F, out, Vector{-mu0}, nfourier, nphi);
  const auto               no = static_cast<Index>(out.size()), ni = static_cast<Index>(mu.size());
  rtepack::muelmat_tensor4 Cd(nfourier, 1, no, ni), Sd(nfourier, 1, no, ni);
  rtepack::muelmat_tensor3 Cb(nfourier, 1, no), Sb(nfourier, 1, no);
  for (Index m = 0; m < nfourier; m++)
    for (Index i = 0; i < no; i++) {
      for (Index j = 0; j < ni; j++) {
        Cd[m, 0, i, j] = d.C[m, i, j];
        Sd[m, 0, i, j] = d.S[m, i, j];
      }
      Cb[m, 0, i] = b.C[m, i, 0];
      Sb[m, 0, i] = b.S[m, i, 0];
    }
  return {.diffuse = vdisort::combine_phase_matrices(Cd, Sd), .beam = vdisort::combine_beam_phase_matrices(Cb, Sb)};
}

const rtepack::stokvec beam{1.0, 0.0, 0.0, 0.0};

vdisort::main_data solver(Index nquad, Numeric layer_tau, Numeric layer_omega, const CompactPlanarMuelmatVector& s) {
  const Vector mu    = streams(nquad);
  const Index  NF    = nquad;
  auto         phase = combined(s, mu, mu, NF);
  return vdisort::main_data(nquad,
                            NF,
                            AscendingGrid{Vector{layer_tau}},
                            Vector{layer_omega},
                            std::move(phase.diffuse),
                            rtepack::stokvec_tensor3(2, NF, nquad / 2, rtepack::stokvec{}),
                            rtepack::stokvec_tensor3(2, NF, nquad / 2, rtepack::stokvec{}),
                            rtepack::stokvec_matrix(1, 0),
                            {},
                            mu0,
                            beam,
                            phi0,
                            std::move(phase.beam));
}

//! The user directions: upward at the top, downward at the bottom (around the beam at -mu0)
const Vector up_mu{0.25, 0.6, 0.95}, down_mu{0.4, 0.5, 0.55, 0.85};
const Vector user_phi{0.0, 0.15, 0.5, Constant::pi / 2, Constant::pi};

Vector user_directions() {
  Vector u(up_mu.size() + down_mu.size());
  for (Size i = 0; i < up_mu.size(); i++) u[i] = up_mu[i];
  for (Size i = 0; i < down_mu.size(); i++) u[up_mu.size() + i] = -down_mu[i];
  return u;
}

//! The radiances [phi, user] that matter: upward ones at the top, downward ones at the bottom
using field = rtepack::stokvec_matrix;

field reference(Index nquad, Numeric tau) {
  const auto                 c = layer_series();
  CompactPlanarMuelmatVector truncated(nquad);
  for (Index l = 0; l < nquad; l++) truncated[l] = c[l];
  const auto   v    = solver(nquad, tau, omega, truncated);
  const Vector user = user_directions();
  const auto   p    = combined(truncated, user, streams(nquad), nquad);

  // One evaluation at both levels and every azimuth, with the user modes formed once
  rtepack::stokvec_tensor3 u(2, user_phi.size(), user.size());
  v.ungridded_u_user(u, AscendingGrid{Vector{0.0, tau}}, user_phi, user, p.diffuse, p.beam);
  field out(user_phi.size(), user.size());
  for (Size k = 0; k < user_phi.size(); k++)
    for (Size i = 0; i < user.size(); i++) out[k, i] = u[user[i] > 0 ? 0 : 1, k, i];
  return out;
}

struct delta_m_fields {
  field uncorrected, corrected, tms, ims;
};

delta_m_fields delta_m(Index nquad, Numeric tau) {
  const auto    c     = layer_series();
  const auto    s     = delta_m_split(c, nquad);
  const Numeric scale = 1 - omega * s.f;
  const auto    v     = solver(nquad, scale * tau, omega * (1 - s.f) / scale, s.transport);
  const Vector  user  = user_directions();
  const auto    p     = combined(s.transport, user, streams(nquad), nquad);

  const auto phase = [](const CompactPlanarMuelmatVector& series) {
    return [F = legendre_matrix(series)](Index, Numeric out_mu, Numeric out_phi, Numeric in_mu, Numeric in_phi) {
      return lab_frame(F, in_mu, in_phi, out_mu, out_phi);
    };
  };
  const vdisort::delta_m_correction_cache correction(AscendingGrid{Vector{tau}},
                                                     Vector{omega},
                                                     Vector{s.f},
                                                     mu0,
                                                     phi0,
                                                     beam,
                                                     user,
                                                     user_phi,
                                                     phase(c),
                                                     phase(s.removed),
                                                     48,
                                                     96);

  delta_m_fields out{field(user_phi.size(), user.size()),
                     field(user_phi.size(), user.size()),
                     field(user_phi.size(), user.size()),
                     field(user_phi.size(), user.size())};
  for (Size k = 0; k < user_phi.size(); k++) {
    for (const Numeric level : {0.0, tau}) {
      vdisort::user_u_data    u, uc;
      rtepack::stokvec_vector ims, tms;
      v.u_user(u, scale * level, user_phi[k], user, p.diffuse, p.beam);
      v.u_user_corr(uc, ims, tms, level, static_cast<Index>(k), correction, p.diffuse, p.beam);
      for (Size i = 0; i < user.size(); i++)
        if ((user[i] > 0) == (level == 0.0)) {
          out.uncorrected[k, i] = u.intensities[i];
          out.corrected[k, i]   = uc.intensities[i];
          out.tms[k, i]         = tms[i];
          out.ims[k, i]         = ims[i];
        }
    }
  }
  return out;
}

//! max |a - b| per Stokes component / max |I of b|
Vector4 deviation(const field& a, const field& b) {
  Numeric scale = 0.0;
  for (const auto& x : b | by_elem) scale = std::max(scale, std::abs(x.I()));
  Vector4 d{};
  for (Index k = 0; k < a.nrows(); k++)
    for (Index i = 0; i < a.ncols(); i++)
      for (Index s = 0; s < 4; s++) d[s] = std::max(d[s], std::abs(a[k, i][s] - b[k, i][s]) / scale);
  return d;
}

std::string format(const Vector4& d) {
  return std::format("I {:9.3e}  Q {:9.3e}  U {:9.3e}  V {:9.3e}", d[0], d[1], d[2], d[3]);
}

Numeric max_I(const field& f) {
  Numeric scale = 0.0;
  for (const auto& x : f | by_elem) scale = std::max(scale, std::abs(x.I()));
  return scale;
}

/* A layer of tau = 1: the reference converges, and at every NQuad the
   corrected radiances are closer to it than delta-M alone in I (by orders
   of magnitude), Q and U, and no farther in V, whose error comes from
   multiple scattering through the small F34.  The corrected error falls
   with NQuad. */
void test_corrections() {
  const Numeric tau   = 1.0;
  const auto    ref   = reference(64, tau);
  const auto    check = reference(48, tau);
  const auto    self  = deviation(check, ref);
  std::cout << std::format("tau 1, reference 48 vs 64 streams:          {}\n", format(self));
  // 48 streams truncate the series at 1e-4 of the phase function, 64 at 1e-6
  require(self[0] < 3e-5 and std::max({self[1], self[2], self[3]}) < 1e-6,
          "The reference must converge in the number of streams");

  Vector4 polarization{};
  for (const auto& x : ref | by_elem)
    for (Index s = 1; s < 4; s++) polarization[s] = std::max(polarization[s], std::abs(x[s]) / max_I(ref));
  std::cout << std::format("tau 1, reference max |Q|, |U|, |V| / max I:  {:.2e}, {:.2e}, {:.2e}\n",
                           polarization[1],
                           polarization[2],
                           polarization[3]);
  require(polarization[1] > 1e-2 and polarization[2] > 1e-2 and polarization[3] > 1e-5,
          "Every Stokes component must be exercised");

  Vector4 previous{1, 1, 1, 1};
  for (const Index nquad : {8, 12, 16, 24}) {
    const auto x = delta_m(nquad, tau);
    const auto u = deviation(x.uncorrected, ref), c = deviation(x.corrected, ref);
    std::cout << std::format(
        "NQuad {:2d}, f = {:.2e}: delta-M          {}\n", nquad, delta_m_split(layer_series(), nquad).f, format(u));
    std::cout << std::format("                     with TMS and IMS {}\n", format(c));
    require(c[0] < u[0] / 100 and c[1] < u[1] / 1.5 and c[2] < u[2] / 1.5 and c[3] <= 1.01 * u[3],
            std::format("NQuad {}: the corrections must improve I, Q and U and not worsen V", nquad));
    require(c[0] < previous[0] and c[1] < previous[1] and c[2] < 2 * previous[2],
            std::format("NQuad {}: the corrected error must fall with NQuad", nquad));
    previous = c;
  }
}

/* A thin layer (tau = 0.02), where single scattering dominates: TMS alone
   must remove most of delta-M's error in every component it affects,
   as it replaces the truncated single scattering by the exact one.  The
   reference has 48 streams: they agree with 64 to 1e-7 in Q, U and V at
   tau = 1, far below the errors compared here. */
void test_thin_layer() {
  const Numeric tau = 0.02;
  const auto    ref = reference(48, tau);
  for (const Index nquad : {8, 16}) {
    const auto x        = delta_m(nquad, tau);
    field      tms_only = x.corrected;
    for (Index k = 0; k < tms_only.nrows(); k++)
      for (Index i = 0; i < tms_only.ncols(); i++) tms_only[k, i] -= x.ims[k, i];
    const auto u = deviation(x.uncorrected, ref), t = deviation(tms_only, ref);
    std::cout << std::format("tau 0.02, NQuad {:2d}: delta-M          {}\n", nquad, format(u));
    std::cout << std::format("                    with TMS         {}\n", format(t));
    require(t[0] < u[0] / 100 and t[1] < u[1] / 10 and t[2] < u[2] / 10,
            std::format("NQuad {}: TMS must remove most of the single-scattering error in I, Q and U", nquad));
  }
}
}  // namespace

int main() try {
  test_corrections();
  test_thin_layer();
  std::cout << "vdisort delta-M test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

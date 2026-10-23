/* Azimuthally random orientation (ARO) of the T-matrix

   1. Spheres against Mie: ARTS's TRO-to-ARO conversion (scattering::tro_lab_frame) of Mie's scattering matrix at
      the exact scattering angles, in all 16 elements and for any tilt, and the extinction and absorption against
      Mie's cross sections, to the requested accuracy.
   2. A uniform tilt distribution is random orientation: tilts on Gauss-Legendre nodes in cos(beta) against
      tro_lab_frame of tmatrix::random's analytic orientation average, to the requested accuracy, for a spheroid, an
      odd Chebyshev particle and a generalized Chebyshev particle, whose random orientation has no Fortran reference.
   3. The averages over alpha and the absorption integral are exact: against fixed_batch with twice their nodes, to
      rounding, also for a vertical axis, whose single node is the shortcut of invariance.  Half the nodes in alpha or
      phi', or one fewer in cos(theta'), fail; the highest harmonics within the band limits are below rounding.
   4. The extinction matrix is the coherent transmission's: K = (lambda / 2) [M(I - i S) - M(I + i S)] of the
      forward amplitude S, with M ARTS's Mueller matrix of a Jones matrix, which fixes K34's sign in ARTS's
      convention; to rounding.
   5. Mirror symmetry, Z(-delta_aa) = mirror(Z(delta_aa)), to rounding, and a lossless particle absorbs nothing, to
      the accuracy of its T-matrix.
*/
#include <arts_constants.h>
#include <legendre.h>
#include <mie.h>
#include <phase_matrix.h>
#include <tmatrix.h>
#include <wigner_functions.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <print>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using tmatrix::Geometry;

constexpr Numeric pi      = Constant::pi;
constexpr Numeric lambda  = 2.0 * pi;
constexpr Numeric epsilon = std::numeric_limits<Numeric>::epsilon();

int failures = 0;

void require(bool ok, const std::string& what) {
  if (not ok) {
    std::println("FAIL {}", what);
    ++failures;
  }
}

const Vector za{0, 30, 60, 90, 120, 150, 180}, daa{-180, -90, -30, 0, 30, 90, 180};

//! The ARO phase matrices of a scattering matrix F(theta [rad]) per particle, by ARTS's TRO-to-ARO conversion
MuelmatTensor3 tro_reference(auto&& F) {
  const auto one   = std::make_shared<const Vector>(Vector{1.0});
  const auto phase = scattering::tro_lab_frame(
      one,
      one,
      std::make_shared<const Vector>(za),
      std::make_shared<const Vector>(daa),
      std::make_shared<const scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(za)),
      [&F](Numeric theta, rtepack::compact_planar_muelmat_matrix& tro) { tro[0, 0] = F(theta); });
  MuelmatTensor3 out(za.size(), daa.size(), za.size());
  for (Size i = 0; i < za.size(); ++i)
    for (Size a = 0; a < daa.size(); ++a)
      for (Size s = 0; s < za.size(); ++s) out[i, a, s] = phase[0, 0, i, a, s];
  return out;
}

Numeric largest(const MuelmatTensor3& Z) {
  Numeric out = 0.0;
  for (Size i = 0; i < Z.size(); ++i)
    for (Numeric e : Z.elem_at(i).data) out = std::max(out, std::abs(e));
  return out;
}

Numeric distance(const MuelmatTensor3& Z, const MuelmatTensor3& R) {
  Numeric out = 0.0;
  for (Size i = 0; i < Z.size(); ++i)
    for (Size e = 0; e < 16; ++e) out = std::max(out, std::abs(Z.elem_at(i).data[e] - R.elem_at(i).data[e]));
  return out;
}

//! 1. A sphere of any tilt against Mie
void sphere(Numeric r, Complex m, Numeric tilt) {
  const Numeric accuracy = 1e-6;
  const auto    aro      = tmatrix::azimuthally_random(
      r, lambda, tmatrix::Shape{}, m.real(), m.imag(), za, daa, za, {tilt}, {1.0}, accuracy);
  const auto            reference = tro_reference([&](Numeric theta) {
    return scattering::MieSphere(lambda, r, m, Vector{Conversion::rad2deg(theta)}).get_scattering_matrix_compact()[0];
  });
  scattering::MieSphere mie(lambda, r, m, Vector{0.0});
  const Numeric         phase = distance(aro.phase, reference) / largest(reference);
  Numeric               cross = 0.0;
  for (Size i = 0; i < za.size(); ++i) {
    const auto K = aro.extinction[i];
    const auto a = aro.absorption[i];
    cross        = std::max({cross,
                             std::abs(K.A() / mie.get_extinction_coeff() - 1.0),
                             std::abs(a.I() / mie.get_absorption_coeff() - 1.0),
                             std::abs(K.B()) / K.A(),
                             std::abs(K.W()) / K.A(),
                             std::abs(a.Q()) / K.A()});
  }
  std::println("sphere r={} m={} tilt {}: phase {:.1e}, cross sections {:.1e}", r, m.real(), tilt, phase, cross);
  require(phase <= accuracy and cross <= accuracy, std::format("sphere r={} tilt {} against Mie", r, tilt));
}

//! 2. A uniform tilt distribution against random orientation
void uniform(Numeric r, const tmatrix::Shape& shape) {
  const Numeric accuracy = 1e-6;
  const Complex m{1.5, 0.02};
  const auto    tro = tmatrix::random(r, lambda, shape, m.real(), m.imag(), 0, accuracy, EquivalentRadius::volume);
  const auto    fixed =
      tmatrix::fixed(r, lambda, shape, m.real(), m.imag(), 0, 0, 0, 0, 0, 0, accuracy, EquivalentRadius::volume);
  // Gauss-Legendre in cos(beta), with nodes to spare for products of four Wigner functions of degree nmax
  const Index nodes = 4 * fixed.order + 4;
  Vector      x(nodes), w(nodes), tilts(nodes);
  Legendre::GaussLegendre(x, w);
  for (Index k = 0; k < nodes; ++k) tilts[k] = Conversion::rad2deg(std::acos(x[k]));
  const auto aro = tmatrix::azimuthally_random(
      r, lambda, shape, m.real(), m.imag(), za, daa, za, tilts, w, accuracy, EquivalentRadius::volume);

  // The Legendre series of the differential cross section, sum_l a_l sqrt((2 l + 1) / 4 pi) P_l(cos(theta))
  const auto    reference = tro_reference([&](Numeric theta) {
    rtepack::compact_planar_muelmat F{};
    Numeric                         p_previous = 0.0, p = 1.0;
    const Numeric                   u = std::cos(theta);
    for (Size l = 0; l < tro.legendre.size(); ++l) {
      if (l > 0) {
        const Numeric ln   = static_cast<Numeric>(l);
        const Numeric next = ((2 * ln - 1) * u * p - (ln - 1) * p_previous) / ln;
        p_previous         = std::exchange(p, next);
      }
      F += (std::sqrt((2.0 * static_cast<Numeric>(l) + 1.0) / (4.0 * pi)) * p * tro.scattering / (4.0 * pi)) *
           tro.legendre[l];
    }
    return F;
  });
  const Numeric phase     = distance(aro.phase, reference) / largest(reference);
  Numeric       cross     = 0.0;
  for (Size i = 0; i < za.size(); ++i) {
    const auto K = aro.extinction[i];
    const auto a = aro.absorption[i];
    cross        = std::max({cross,
                             std::abs(K.A() / tro.extinction - 1.0),
                             std::abs(a.I() - (tro.extinction - tro.scattering)) / tro.extinction,
                             std::abs(K.B()) / K.A(),
                             std::abs(K.W()) / K.A(),
                             std::abs(a.Q()) / K.A()});
  }
  const auto name = std::format("uniform tilt {} {} {}", shape.kind(), shape.degree(), shape.parameter());
  std::println("{}: phase {:.1e}, cross sections {:.1e}", name, phase, cross);
  require(phase <= accuracy and cross <= accuracy, name);
}

//! The extinction (A, B, W) of a forward amplitude by the coherent transmission, (lambda / 2) [M(I - iS) - M(I + iS)]
rtepack::compact_planar_propmat transmission_extinction(const Jones& S) {
  const Jones   iS{Complex{0.0, 1.0} * S};
  const Muelmat K = 0.5 * lambda * (rtepack::to_muelmat(Jones{} - iS) - rtepack::to_muelmat(Jones{} + iS));
  return {K[0, 0], K[0, 1], K[2, 3]};
}

//! 3 and 4. The averages over alpha, the absorption integral and the extinction against fixed_batch
void exact(Numeric tilt, Numeric r, Numeric accuracy) {
  const auto    shape = tmatrix::Shape::spheroid(1.7);
  const Complex m{1.5, 0.05};
  const Vector  za_inc{20, 90, 160}, delta{-120, 0, 75}, za_scat{0, 70, 180};
  const auto    aro = tmatrix::azimuthally_random(
      r, lambda, shape, m.real(), m.imag(), za_inc, delta, za_scat, {tilt}, {1.0}, accuracy);
  const Index nmax = aro.order;
  // Twice the nodes of each average and integral
  const Index M = 8 * nmax + 2, n_theta = 2 * nmax + 2, n_phi = 4 * nmax + 2;
  Vector      x(n_theta), wx(n_theta);
  Legendre::GaussLegendre(x, wx);
  const auto alpha   = [M](Index j) { return 360.0 * static_cast<Numeric>(j) / static_cast<Numeric>(M); };
  const auto azimuth = [](Numeric phi) {
    const Numeric out = std::fmod(phi, 360.0);
    return out < 0 ? out + 360.0 : out;
  };

  std::vector<Geometry> phase_geometries, forward_geometries, node_geometries;
  for (Index j = 0; j < M; ++j)
    for (Numeric zi : za_inc) {
      forward_geometries.push_back({.theta_incident  = zi,
                                    .theta_scattered = zi,
                                    .phi_incident    = 90,
                                    .phi_scattered   = 90,
                                    .alpha           = alpha(j),
                                    .beta            = tilt});
      for (Numeric a : delta)
        for (Numeric zs : za_scat)
          phase_geometries.push_back({.theta_incident  = zi,
                                      .theta_scattered = zs,
                                      .phi_incident    = 90,
                                      .phi_scattered   = azimuth(90 - a),
                                      .alpha           = alpha(j),
                                      .beta            = tilt});
      for (Index q = 0; q < n_theta; ++q)
        for (Index p = 0; p < n_phi; ++p)
          node_geometries.push_back({.theta_incident  = Conversion::rad2deg(std::acos(x[q])),
                                     .theta_scattered = zi,
                                     .phi_incident    = 360.0 * static_cast<Numeric>(p) / static_cast<Numeric>(n_phi),
                                     .phi_scattered   = 90,
                                     .alpha           = alpha(j),
                                     .beta            = tilt});
    }
  const auto phase   = tmatrix::fixed_batch(r, lambda, shape, m.real(), m.imag(), phase_geometries, accuracy);
  const auto forward = tmatrix::fixed_batch(r, lambda, shape, m.real(), m.imag(), forward_geometries, accuracy);
  const auto nodes   = tmatrix::fixed_batch(r, lambda, shape, m.real(), m.imag(), node_geometries, accuracy);

  // The references, and the moduli of their terms, the scale of their rounding
  const Size     ni = za_inc.size(), na = delta.size(), ns = za_scat.size(), per = n_theta * n_phi;
  MuelmatTensor3 Z(ni, na, ns, Muelmat{0.0});
  std::vector<rtepack::compact_planar_propmat> K(ni);
  std::vector<rtepack::compact_planar_stokvec> a(ni);
  Numeric                                      scale_Z = 0.0, transmission = 0.0;
  std::vector<Numeric>                         scale_a(ni, 0.0);
  for (Index j = 0; j < M; ++j)
    for (Size i = 0; i < ni; ++i) {
      const Size                            row = j * ni + i;
      const auto                            k   = transmission_extinction(forward[row].amplitude);
      const auto                            S   = forward[row].amplitude;
      const rtepack::compact_planar_propmat explicit_k{lambda * (S[0, 0] + S[1, 1]).imag(),
                                                       lambda * (S[0, 0] - S[1, 1]).imag(),
                                                       lambda * (S[0, 0] - S[1, 1]).real()};
      for (Size e = 0; e < 3; ++e)
        transmission = std::max(transmission, std::abs(k.data[e] - explicit_k.data[e]) / k.A());
      K[i] += (1.0 / static_cast<Numeric>(M)) * k;
      for (Size aa = 0; aa < na; ++aa)
        for (Size s = 0; s < ns; ++s) {
          const Muelmat& z  = phase[(row * na + aa) * ns + s].phase;
          Z[i, aa, s]      += (1.0 / static_cast<Numeric>(M)) * z;
          for (Numeric e : z.data) scale_Z = std::max(scale_Z, std::abs(e));
        }
      Numeric I = 0.0, Q = 0.0;
      for (Size n = 0; n < per; ++n) {
        const Muelmat& z  = nodes[row * per + n].phase;
        const Numeric  w  = wx[n / n_phi] * 2.0 * pi / static_cast<Numeric>(n_phi);
        I                += w * z[0, 0];
        Q                += w * z[1, 0];
        scale_a[i]       += (1.0 / static_cast<Numeric>(M)) * w * std::abs(z[0, 0]);
      }
      a[i] += (1.0 / static_cast<Numeric>(M)) * rtepack::compact_planar_stokvec{k.A() - I, k.B() - Q};
    }

  Numeric worst_Z = distance(aro.phase, Z) / scale_Z, worst_K = 0.0, worst_a = 0.0;
  for (Size i = 0; i < ni; ++i) {
    for (Size e = 0; e < 3; ++e)
      worst_K = std::max(worst_K, std::abs(aro.extinction[i].data[e] - K[i].data[e]) / K[i].A());
    for (Size e = 0; e < 2; ++e)
      worst_a = std::max(worst_a, std::abs(aro.absorption[i].data[e] - a[i].data[e]) / (K[i].A() + scale_a[i]));
  }
  // Rounding: a few ulp of the largest term, grown as the square root of the number of terms summed, M amplitudes
  // for the phase and extinction, M n_theta n_phi phase matrices for the absorption
  const Numeric average  = 16 * epsilon * std::sqrt(static_cast<Numeric>(M));
  const Numeric integral = 16 * epsilon * std::sqrt(static_cast<Numeric>(M * n_theta * n_phi));
  std::println(
      "exact tilt {} nmax {}: phase {:.1e}, extinction {:.1e} (bound {:.1e}), absorption {:.1e} (bound {:.1e}), "
      "transmission {:.1e}",
      tilt,
      nmax,
      worst_Z,
      worst_K,
      average,
      worst_a,
      integral,
      transmission);
  require(worst_Z <= average and worst_K <= average and worst_a <= integral and transmission <= 16 * epsilon,
          std::format("exact quadratures at tilt {}", tilt));
}

//! 5. Mirror symmetry of the azimuthal average, and no absorption without loss
void symmetry() {
  const Numeric r = 1.6, eps = 0.6, accuracy = 1e-6;
  const auto    aro = tmatrix::azimuthally_random(
      r, lambda, tmatrix::Shape::spheroid(eps), 1.5, 0.0, za, daa, za, {35.0, 80.0}, {0.3, 0.7}, accuracy);
  Numeric mirror = 0.0, absorption = 0.0;
  for (Size i = 0; i < za.size(); ++i) {
    for (Size a = 0; a < daa.size(); ++a)
      for (Size s = 0; s < za.size(); ++s) {
        const Muelmat& Z  = aro.phase[i, a, s];
        const Muelmat  Zm = rtepack::mirror(aro.phase[i, daa.size() - 1 - a, s]);
        for (Size e = 0; e < 16; ++e) mirror = std::max(mirror, std::abs(Z.data[e] - Zm.data[e]));
      }
    absorption = std::max(absorption, std::abs(aro.absorption[i].I()) / aro.extinction[i].A());
  }
  mirror /= largest(aro.phase);
  std::println("mirror symmetry {:.1e}, lossless absorption {:.1e}", mirror, absorption);
  require(mirror <= 64 * epsilon, "mirror symmetry of the azimuthal average");
  require(absorption <= accuracy, "a lossless particle absorbs");
}
}  // namespace

int main() try {
  make_wigner_ready(250, 20000000, 6);
  for (Numeric tilt : {0.0, 30.0, 90.0}) sphere(1.2, {1.33, 0.01}, tilt);
  sphere(2.5, {1.78, 0.05}, 45.0);
  uniform(1.5, tmatrix::Shape::spheroid(1.6));
  uniform(1.5, tmatrix::Shape::chebyshev(3, 0.12));
  // The coefficients of the original's DROP routine, as test data of a generalized Chebyshev particle
  uniform(1.5,
          tmatrix::Shape::generalized_chebyshev(
              Vector{-0.0481, 0.0359, -0.1263, 0.0244, 0.0091, -0.0099, 0.0015, 0.0025, -0.0016, -0.0002, 0.0010}));
  for (Numeric tilt : {0.0, 90.0}) exact(tilt, 1.8, 1e-3);
  exact(40.0, 3.0, 1e-2);
  symmetry();
  if (failures) throw std::runtime_error(std::format("{} comparisons failed", failures));
} catch (const std::exception& e) {
  std::println(stderr, "{}", e.what());
  return 1;
}

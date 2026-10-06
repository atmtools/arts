/* Tests of vdisort::scattering_optics and vdisort::main_data_from_path
   (vdisort_arts.h), and of the Stokes and azimuth conventions that relate
   ARTS's scattering data to VDISORT, RT3 and RT4.

   Where the expected answers come from:
   - V1: closed forms of the Fourier coefficients of Rayleigh scattering for
     m = 0, 1, 2, derived below from the dipole (Jones-matrix) picture in the
     meridional basis, independently of any rotation formula.  ARTS's
     GasScatterer gives the laboratory-frame phase matrix;
     vdisort::scattering_optics takes its Fourier coefficients.
   - V2: ARTS's own laboratory-frame phase matrix (to_lab_frame in
     scattering/phase_matrix.h, Mishchenko-style spherical trigonometry)
     against the vector-geometry construction of lab-frame.h (this
     directory) for a polarizing Mie particle (F12, F34 != 0), at generic
     directions.  The two formalisms agree element by element when ARTS's
     propagation zenith angle za and azimuth aa map to mu = cos(za) and
     phi = -aa, with the same I, Q, U, V and the same F, including F34.
     The other azimuth sense and a negated F34 must both fail.  Exact and
     near-forward pairs must keep the diagonal of F (phase_matrix.h used
     Z22 = -F22 there).
   - V3: vdisort::scattering_optics for the same particle as a Legendre
     series (exact from its Gauss-Legendre nodes), which takes the Fourier
     modes of ARTS's laboratory-frame phase matrix, against the Fourier
     coefficients of the vector-geometry phase matrix of lab-frame.h (V2's
     mapping) of the same series, diffuse and beam column.
   - V4: main_data_from_path against the path data it is given.
   - V5: azimuthally randomly oriented particle data: Rayleigh's own
     laboratory-frame phase matrix, gridded on the streams and over all
     scattering zenith angles, must give the GasScatterer's coefficients to
     the accuracy of its 1 deg grids; data without the phase integral, and
     optics that depend on the direction, must be errors. */
#include <arts_constants.h>
#include <arts_conversions.h>
#include <legendre.h>
#include <physics_funcs.h>
#include <vdisort_arts.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include "lab-frame.h"

namespace {
constexpr Numeric pi = Constant::pi;

void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

void require_error(const std::function<void()>& f, const std::string& what) {
  bool threw = false;
  try {
    f();
  } catch (const std::exception&) { threw = true; }
  require(threw, std::format("Expected an error: {}", what));
}

AtmPoint air(Numeric p, Numeric t) {
  AtmPoint a;
  a.pressure    = p;
  a.temperature = t;
  return a;
}

constexpr Numeric cross_section = 1e-30;  // m^2

ArrayOfScatteringSpecies rayleigh_species() {
  ArrayOfScatteringSpecies s;
  s.add(GasScatterer{scattering::ConstantGasScattering{cross_section}, scattering::RayleighGasScattering{0.0}});
  return s;
}

/* Rayleigh scattering from (mu_i, phi) into (mu_o, 0), signed cosines, in the
   meridional basis e_v = (mu cos, mu sin, -s), e_h = (-sin, cos, 0).  The
   dipole radiates the part of E normal to k_out, so the Jones matrix is
   J = [[a, b], [d, e]] with a = e_vo . e_vi = A + B cos(phi), b = e_vo . e_hi =
   -mu_o sin(phi), d = e_ho . e_vi = mu_i sin(phi), e = e_ho . e_hi = cos(phi),
   A = s_o s_i, B = mu_o mu_i.  For a real J, with Q = |E_v|^2 - |E_h|^2 and
   U = 2 Re(E_v E_h*),
     M11 = (a2 + b2 + d2 + e2) / 2,  M12 = (a2 - b2 + d2 - e2) / 2,  M13 = ab + de,
     M21 = (a2 + b2 - d2 - e2) / 2,  M22 = (a2 - b2 - d2 + e2) / 2,  M23 = ab - de,
     M31 = ad + be,                  M32 = ad - be,                  M33 = ae + bd,
     M44 = ae - bd,
   and Z = 3/2 M (so that Z11 = 3/4 (1 + cos^2 Theta)).  Expanding in cos m phi
   and sin m phi gives, with a = mu_o^2 and b = mu_i^2 below and without the
   factor 2 - delta_m0 (C^m, S^m = (1 / 2 pi) int Z {cos, sin}(m phi) dphi):
     m = 0: C11 = 3/8 (3 - a - b + 3ab), C12 = 3/8 (1 - 3a)(1 - b),
            C21 = 3/8 (1 - a)(1 - 3b), C22 = 9/8 (1 - a)(1 - b), C44 = 3/2 mu_o mu_i;
     m = 1: C11 = C12 = C21 = C22 = 3/4 A B, C33 = C44 = 3/4 A,
            S13 = S23 = -3/4 mu_o A, S31 = S32 = 3/4 mu_i A;
     m = 2: C11 = 3/16 (1 - a)(1 - b), C12 = -3/16 (1 - a)(1 + b),
            C21 = -3/16 (1 + a)(1 - b), C22 = 3/16 (1 + a)(1 + b), C33 = 3/4 mu_o mu_i,
            S13 = 3/8 mu_i (1 - a), S23 = -3/8 mu_i (1 + a),
            S31 = -3/8 mu_o (1 - b), S32 = 3/8 mu_o (1 + b);
   every other element is 0, and every m >= 3 coefficient is 0. */
std::pair<rtepack::muelmat, rtepack::muelmat> rayleigh_modes(Index m, Numeric mo, Numeric mi) {
  const Numeric    a = mo * mo, b = mi * mi;
  const Numeric    A = std::sqrt((1 - a) * (1 - b)), B = mo * mi;
  rtepack::muelmat C{0.0}, S{0.0};
  if (m == 0) {
    C[0, 0] = 3.0 / 8.0 * (3 - a - b + 3 * a * b);
    C[0, 1] = 3.0 / 8.0 * (1 - 3 * a) * (1 - b);
    C[1, 0] = 3.0 / 8.0 * (1 - a) * (1 - 3 * b);
    C[1, 1] = 9.0 / 8.0 * (1 - a) * (1 - b);
    C[3, 3] = 1.5 * mo * mi;
  } else if (m == 1) {
    C[0, 0] = C[0, 1] = C[1, 0] = C[1, 1] = 0.75 * A * B;
    C[2, 2] = C[3, 3] = 0.75 * A;
    S[0, 2] = S[1, 2] = -0.75 * mo * A;
    S[2, 0] = S[2, 1] = 0.75 * mi * A;
  } else if (m == 2) {
    C[0, 0] = 3.0 / 16.0 * (1 - a) * (1 - b);
    C[0, 1] = -3.0 / 16.0 * (1 - a) * (1 + b);
    C[1, 0] = -3.0 / 16.0 * (1 + a) * (1 - b);
    C[1, 1] = 3.0 / 16.0 * (1 + a) * (1 + b);
    C[2, 2] = 0.75 * mo * mi;
    S[0, 2] = 3.0 / 8.0 * mi * (1 - a);
    S[1, 2] = -3.0 / 8.0 * mi * (1 + a);
    S[2, 0] = -3.0 / 8.0 * mo * (1 - b);
    S[2, 1] = 3.0 / 8.0 * mo * (1 + b);
  }
  return {C, S};
}

//! VDISORT's 16 streams (upward first), and the extra cosines +-0.35 and +-1
Vector test_cosines() {
  Vector mu(16), inv(16), w(8);
  disort_common::initialize_streams(mu, inv, w);
  Vector out(20);
  for (Index i = 0; i < 16; i++) out[i] = mu[i];
  out[16] = 0.35;
  out[17] = 1.0;
  out[18] = -0.35;
  out[19] = -1.0;
  return out;
}

//! V1: Rayleigh against the dipole closed forms, m = 0 .. 3, diffuse and beam column
void test_rayleigh() {
  const auto   species = rayleigh_species();
  const auto   atm     = air(9e4, 270.0);
  const auto   sigma   = cross_section * number_density(atm.pressure, atm.temperature);
  const Vector mu      = test_cosines();
  Vector       mu_in(mu.size() + 1);
  std::ranges::copy(mu, mu_in.begin());
  mu_in[mu.size()] = -0.6;  // a beam

  {
    const auto f = vdisort::scattering_optics(species, atm, 50e9, mu, mu_in, 4, 1e-12);
    require(std::abs(f.extinction - sigma) < 1e-14 * sigma and std::abs(f.scattering - sigma) < 1e-14 * sigma,
            "V1: Rayleigh extinction and scattering must be sigma");
    Numeric d = 0.0;
    for (Index m = 0; m < 4; m++) {
      for (Index o = 0; o < static_cast<Index>(mu.size()); o++) {
        for (Index i = 0; i < static_cast<Index>(mu_in.size()); i++) {
          const auto [C, S] = rayleigh_modes(m, mu[o], mu_in[i]);
          for (Index a = 0; a < 4; a++)
            for (Index b = 0; b < 4; b++)
              d = std::max({d, std::abs(f.cosine[m, o, i][a, b] - C[a, b]), std::abs(f.sine[m, o, i][a, b] - S[a, b])});
        }
      }
    }
    std::cout << std::format(
        "V1 Rayleigh GasScatterer, m = 0..3, 16 streams, +-0.35, +-1 and the beam -0.6: max |C^m, S^m - closed "
        "form| {:.1e} (tolerance 1e-13)\n",
        d);
    require(d <= 1e-13,
            std::format("V1: the Fourier coefficients of ARTS's Rayleigh GasScatterer must equal the dipole closed "
                        "forms to 1e-13, got {:.2e}",
                        d));
  }

  require_error([&] { (void)vdisort::scattering_optics(species, atm, 50e9, Vector{0.0}, mu, 2, 1e-3); }, "mu = 0");
  require_error([&] { (void)vdisort::scattering_optics(species, atm, 50e9, mu, mu, 0, 1e-3); }, "nfourier 0");
  const auto none = vdisort::scattering_optics(ArrayOfScatteringSpecies{}, atm, 50e9, mu, mu, 2, 1e-3);
  require(none.extinction == 0.0 and
              stdr::all_of(none.cosine | by_elem,
                           [](const auto& z) { return stdr::all_of(z | by_elem, [](Numeric x) { return x == 0.0; }); }),
          "V1: an empty species array must give zero coefficients");
}

/* A polarizing Mie particle: radius 0.25 um, n = 1.44 + 0.01i, at 0.951 um
   (size parameter 1.65), on the nodes of a 2000-point Gauss-Legendre rule
   in cos(Theta).  Its F12 and F34 are both clearly non-zero. */
constexpr Numeric mie_frequency = Constant::speed_of_light / 0.951e-6;
constexpr Index   mie_nodes     = 2000;

using TroSeries = scattering::PhaseMatrixData<Numeric, scattering::Format::TRO, scattering::Representation::Spectral>;

//! The degree of the particle's Legendre series, converged to rounding for size parameter 1.65
constexpr Index mie_degree = 40;

struct mie_case {
  //! The particle as gridded data on the Gauss-Legendre nodes
  ArrayOfScatteringSpecies species;
  //! The particle as its Legendre series to mie_degree
  ArrayOfScatteringSpecies spectral_species;
  TroSeries                series;
  AtmPoint                 atm;
};

mie_case mie() {
  Vector x(mie_nodes), w(mie_nodes), angles(mie_nodes);
  Legendre::GaussLegendre(x, w);
  stdr::reverse(x);  // ascending scattering angles: ARTS interpolates TRO data on an ascending grid
  for (Index i = 0; i < mie_nodes; i++) angles[i] = Conversion::rad2deg(std::acos(x[i]));
  ComplexMatrix index(1, 1);
  index[0, 0] = Complex{1.44, 0.01};
  Vector     t_grid{280.0}, f_grid{mie_frequency}, d{0.5e-6};
  const auto habit = ParticleHabit::sphere(
      t_grid, f_grid, d, scattering::ZenithAngleGrid{scattering::IrregularZenithAngleGrid(angles)}, index, 1000.0);
  const auto prop = ScatteringSpeciesProperty{"mie", ParticulateProperty::NumberDensity};
  mie_case   c{.species = {}, .spectral_species = {}, .series = {}, .atm = air(1e5, 280.0)};
  c.species.add(ScatteringHabit{habit, scattering::PSD{scattering::MonodispersePSD{prop}}, 1.0, 3.0});
  c.atm[prop] = 1e9;

  // a_l = 2 pi sqrt((2 l + 1) / 4 pi) sum_i w_i F(x_i) P_l(x_i), exact for data on the nodes
  using Gridded  = scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Gridded>;
  using Spectral = scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Spectral>;
  const auto& gridded = std::get<Gridded>(habit[0]);
  c.series            = TroSeries(gridded.phase_matrix->get_t_grid(), gridded.phase_matrix->get_f_grid(), mie_degree);
  Vector p(mie_degree + 1);
  for (Index i = 0; i < mie_nodes; i++) {
    Legendre::legendre_polynomials(p, x[i]);
    for (Index l = 0; l <= mie_degree; l++) {
      const Numeric f = 2 * pi * std::sqrt((2.0 * static_cast<Numeric>(l) + 1.0) / (4 * pi)) * w[i] * p[l];
      for (Index k = 0; k < 6; k++) c.series[0, 0, l, k] += f * (*gridded.phase_matrix)[0, 0, i, k];
    }
  }
  const Spectral spectral(gridded.properties,
                          c.series,
                          gridded.extinction_matrix.to_spectral(),
                          gridded.absorption_vector.to_spectral(),
                          gridded.backscatter_matrix,
                          gridded.forwardscatter_matrix);
  c.spectral_species.add(ScatteringHabit{
      ParticleHabit{std::vector<Spectral>{spectral}}, scattering::PSD{scattering::MonodispersePSD{prop}}, 1.0, 3.0});
  return c;
}

//! The scattering matrix of the particle's Legendre series at cos(Theta), times its number density
vdisort_test::tro_matrix series_of(const mie_case& c) {
  return [&c](Numeric cos_theta) {
    const Matrix  f = scattering::tro_legendre::evaluate(c.series.coefficients(0, 0),
                                                        Vector{Conversion::rad2deg(std::acos(cos_theta))});
    const Numeric n = 1e9;
    return vdisort_test::tro_elements{
        .F11 = n * f[0, 0], .F12 = n * f[0, 1], .F22 = n * f[0, 2], .F33 = n * f[0, 3], .F34 = n * f[0, 4], .F44 = n * f[0, 5]};
  };
}

//! ARTS's laboratory-frame Z for the incident propagation zenith za_in, za_out and aa_out - aa_in [deg]
rtepack::muelmat arts_lab_frame(const mie_case& c, Numeric za_in, Numeric delta_aa, Numeric za_out) {
  const auto bulk = c.species.get_bulk_scattering_properties_aro_gridded(
      c.atm,
      Vector{mie_frequency},
      Vector{za_in},
      Vector{delta_aa},
      std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(Vector{za_out})));
  rtepack::muelmat Z{0.0};
  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 4; j++) Z[i, j] = (*bulk.phase_matrix)[0, 0, 0, 0, 0, 4 * i + j];
  return Z;
}

//! The scattering matrix of the species at cos(Theta), via its TRO data at that exact angle
vdisort_test::tro_matrix tro_of(const mie_case& c, Numeric sign34 = 1.0) {
  return [&c, sign34](Numeric cos_theta) {
    const auto bulk = c.species.get_bulk_scattering_properties_tro_gridded(
        c.atm,
        Vector{mie_frequency},
        std::make_shared<scattering::ZenithAngleGrid>(
            scattering::IrregularZenithAngleGrid(Vector{Conversion::rad2deg(std::acos(cos_theta))})));
    const auto& p = *bulk.phase_matrix;
    return vdisort_test::tro_elements{.F11 = p[0, 0, 0, 0],
                                      .F12 = p[0, 0, 0, 1],
                                      .F22 = p[0, 0, 0, 2],
                                      .F33 = p[0, 0, 0, 3],
                                      .F34 = sign34 * p[0, 0, 0, 4],
                                      .F44 = p[0, 0, 0, 5]};
  };
}

Numeric max_abs(const rtepack::muelmat& z) {
  Numeric x = 0.0;
  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 4; j++) x = std::max(x, std::abs(z[i, j]));
  return x;
}

Numeric max_diff(const rtepack::muelmat& a, const rtepack::muelmat& b) {
  Numeric x = 0.0;
  for (Index i = 0; i < 4; i++)
    for (Index j = 0; j < 4; j++) x = std::max(x, std::abs(a[i, j] - b[i, j]));
  return x;
}

//! V2: ARTS's laboratory frame against vector geometry
void test_lab_frame(const mie_case& c) {
  const auto                             F      = tro_of(c);
  const auto                             F_flip = tro_of(c, -1.0);
  const auto                             ref    = F(std::cos(1.2));
  const auto                             rad    = [](Numeric deg) { return Conversion::deg2rad(deg); };
  std::mt19937_64                        rng(42);
  std::uniform_real_distribution<double> za(20.0, 160.0), aa(5.0, 175.0), coin(0.0, 1.0);

  Numeric right = 0.0, sense = 0.0, flipped = 0.0;
  for (Index n = 0; n < 300; n++) {
    const Numeric zi = za(rng), zo = za(rng), daa = coin(rng) < 0.5 ? aa(rng) : 360.0 - aa(rng);
    const auto    Za    = arts_lab_frame(c, zi, daa, zo);
    const Numeric scale = max_abs(Za);
    const Numeric ci = std::cos(rad(zi)), co = std::cos(rad(zo));
    right   = std::max(right, max_diff(Za, vdisort_test::lab_frame(F, ci, 0.0, co, -rad(daa))) / scale);
    sense   = std::max(sense, max_diff(Za, vdisort_test::lab_frame(F, ci, 0.0, co, rad(daa))) / scale);
    flipped = std::max(flipped, max_diff(Za, vdisort_test::lab_frame(F_flip, ci, 0.0, co, -rad(daa))) / scale);
  }
  std::cout << std::format(
      "V2 ARTS laboratory frame vs vector geometry, Mie F12 / F11 = {:.2f}, F34 / F11 = {:.3f} at 69 deg, 300 "
      "generic directions:\n"
      "    mu = cos(za), phi = -aa, same F:   max |dZ| / max |Z| {:.1e} (tolerance 1e-10)\n"
      "    phi = +aa (other azimuth sense):   max |dZ| / max |Z| {:.1e}\n"
      "    F34 negated:                       max |dZ| / max |Z| {:.1e}\n",
      ref.F12 / ref.F11,
      ref.F34 / ref.F11,
      right,
      sense,
      flipped);
  require(right <= 1e-10,
          std::format("V2: ARTS's laboratory-frame phase matrix must equal the vector-geometry phase matrix of the "
                      "same F with mu = cos(za) and phi = -aa to 1e-10, got {:.2e}",
                      right));
  require(sense > 1e-2 and flipped > 1e-2,
          "V2: the other azimuth sense and a negated F34 must both miss ARTS's phase matrix by more than 1e-2");

  /* Exact forward (rays equal) and a near-forward pair (Theta = 8.7e-4 rad).
     The diagonal must be that of F: for forward-scattered Q, Z22 = F22
     (phase_matrix.h had Z22 = -F22, Z33 = -F33 there).  The difference of
     all elements is printed, not asserted. */
  for (auto [z0, daa] : {std::pair{60.0, 0.0}, std::pair{85.0, 0.05}}) {
    const auto Za       = arts_lab_frame(c, z0, daa, z0);
    const auto Zv       = vdisort_test::lab_frame(F, std::cos(rad(z0)), 0.0, std::cos(rad(z0)), -rad(daa));
    Numeric    diagonal = 0.0;
    for (Index i = 1; i < 4; i++) diagonal = std::max(diagonal, std::abs(Za[i, i] / Za[0, 0] - Zv[i, i] / Zv[0, 0]));
    std::cout << std::format(
        "V2 forward pair za {} deg, daa {} deg: Z22 / Z11 = {:+.6f}, Z33 / Z11 = {:+.6f} (vector geometry {:+.6f}, "
        "{:+.6f}); diagonal ratios differ by {:.1e} (tolerance 1e-5); all elements by {:.1e} of max |Z|\n",
        z0,
        daa,
        Za[1, 1] / Za[0, 0],
        Za[2, 2] / Za[0, 0],
        Zv[1, 1] / Zv[0, 0],
        Zv[2, 2] / Zv[0, 0],
        diagonal,
        max_diff(Za, Zv) / max_abs(Zv));
    require(diagonal <= 1e-5,
            std::format("V2: ARTS's phase matrix must keep the diagonal of F for (near) forward scattering (za {} deg, "
                        "daa {} deg); the diagonal ratios differ by {:.2e}",
                        z0,
                        daa,
                        diagonal));
  }
}

//! V3: vdisort::scattering_optics against the Fourier coefficients of the vector-geometry phase matrix
void test_fourier_against_vector_geometry(const mie_case& c) {
  Vector mu(8), inv(8), w(4);
  disort_common::initialize_streams(mu, inv, w);
  Vector mu_in(9);
  mu_in[Range(0, 8)] = mu;
  mu_in[8]           = -0.6;

  constexpr Index NF = 4, N = 512;
  const auto      f = vdisort::scattering_optics(c.spectral_species, c.atm, mie_frequency, mu, mu_in, NF, 1e-9);

  // sigma = 2 pi int F11 dx = sqrt(4 pi) a_0 of the series
  const Numeric sigma = 1e9 * std::sqrt(4 * pi) * c.series[0, 0, 0, 0].real();

  // Incidence at phi_k and scattering at 0, the convention of scattering_optics; 512 midpoints resolve the
  // modes of the series to rounding
  const auto F = series_of(c);
  Numeric    d = 0.0, scale = 0.0;
  for (Index o = 0; o < static_cast<Index>(mu.size()); o++) {
    for (Index i = 0; i < static_cast<Index>(mu_in.size()); i++) {
      rtepack::muelmat_vector C(NF, rtepack::muelmat{0.0}), S(NF, rtepack::muelmat{0.0});
      for (Index k = 0; k < N; k++) {
        const Numeric          phi = 2 * pi * (static_cast<Numeric>(k) + 0.5) / N;
        const rtepack::muelmat Z   = 4 * pi / sigma / N * vdisort_test::lab_frame(F, mu_in[i], phi, mu[o], 0.0);
        for (Index m = 0; m < NF; m++) {
          C[m] += std::cos(static_cast<Numeric>(m) * phi) * Z;
          S[m] += std::sin(static_cast<Numeric>(m) * phi) * Z;
        }
      }
      for (Index m = 0; m < NF; m++) {
        scale = std::max(scale, max_abs(C[m]));
        d     = std::max({d, max_diff(C[m], f.cosine[m, o, i]), max_diff(S[m], f.sine[m, o, i])});
      }
    }
  }
  std::cout << std::format(
      "V3 vdisort::scattering_optics (ARTS's laboratory-frame Fourier modes) vs the Fourier coefficients of the "
      "vector-geometry phase matrix, Mie Legendre series of degree 40, m = 0..3, 8 streams and the beam -0.6: max "
      "|dC|, |dS| / max |C| {:.1e} (tolerance 1e-10)\n",
      d / scale);
  require(d <= 1e-10 * scale,
          std::format("V3: VDISORT's Fourier coefficients must be those of the vector-geometry phase matrix to 1e-10, "
                      "got {:.2e}",
                      d / scale));
}

//! V4: main_data_from_path against its inputs
void test_path() {
  const auto                  species = rayleigh_species();
  ArrayOfPropagationPathPoint ray_path;
  ArrayOfAtmPoint             atm_path;
  ArrayOfPropmatVector        propmat;
  const Vector                alt{4000.0, 1500.0, 0.0}, temp{235.0, 260.0, 285.0}, pres{6e4, 8.5e4, 1e5};
  for (Index l = 0; l < 3; l++) {
    PropagationPathPoint pp;
    pp.pos = {alt[l], 0.0, 0.0};
    pp.los = {180.0, 0.0};
    ray_path.push_back(pp);
    atm_path.push_back(air(pres[l], temp[l]));
    PropmatVector pm(1);
    pm[0].A() = 2e-5 * static_cast<Numeric>(1 + l);
    propmat.push_back(pm);
  }
  const AscendingGrid freq_grid{Vector{89e9}};

  const vdisort::path_settings s{.nquad = 8, .nfourier = 1};
  const auto                   v = vdisort::main_data_from_path(
      ray_path, atm_path, propmat, freq_grid, 0, species, s, vdisort::lambertian_surface{.albedo = 0.3}, 290.0, 2.7);

  Numeric tau = 0.0, dev = 0.0;
  for (Index l = 0; l < 2; l++) {
    const Numeric gas = 2e-5 * (1.5 + static_cast<Numeric>(l));
    const Numeric sca =
        0.5 * cross_section * (number_density(pres[l], temp[l]) + number_density(pres[l + 1], temp[l + 1]));
    tau += (gas + sca) * (alt[l] - alt[l + 1]);
    dev  = std::max({dev, std::abs(v.tau()[l] - tau) / tau, std::abs(v.omega()[l] - sca / (gas + sca))});
  }
  require(dev < 1e-14, std::format("V4: tau and omega from the path data, deviation {:.1e}", dev));

  // Gas-only layer and an empty species array: no scattering, the thermal source alone
  const auto g = vdisort::main_data_from_path(ray_path,
                                              atm_path,
                                              propmat,
                                              freq_grid,
                                              0,
                                              ArrayOfScatteringSpecies{},
                                              s,
                                              vdisort::fresnel_surface{.refractive_index = Complex{2.0, 0.1}},
                                              290.0,
                                              2.7);
  require(stdr::all_of(g.omega(), [](Numeric x) { return x == 0.0; }), "V4: no species must give omega = 0");

  auto bad      = propmat;
  bad[1][0].V() = 1e-7;
  require_error(
      [&] {
        (void)vdisort::main_data_from_path(
            ray_path, atm_path, bad, freq_grid, 0, species, s, vdisort::lambertian_surface{}, 290.0, 2.7);
      },
      "a polarized gas propagation matrix");
  auto zero = propmat;
  for (auto& p : zero) p[0].A() = 0.0;
  require_error(
      [&] {
        (void)vdisort::main_data_from_path(ray_path,
                                           atm_path,
                                           zero,
                                           freq_grid,
                                           0,
                                           ArrayOfScatteringSpecies{},
                                           s,
                                           vdisort::lambertian_surface{},
                                           290.0,
                                           2.7);
      },
      "a layer without optical thickness");
  std::cout << std::format("V4 path builder: tau and omega to {:.1e}, gas-only, 2 error paths\n", dev);
}
//! V5: ARO particle data in VDISORT
void test_aro_data() {
  using ARO = scattering::SingleScatteringData<Numeric, scattering::Format::ARO, scattering::Representation::Gridded>;
  const auto gas = rayleigh_species();
  const auto atm = air(9e4, 270.0);
  Vector     mu(4), inv(4), w(2);
  disort_common::initialize_streams(mu, inv, w);

  // The streams' propagation zenith angles, and 1 deg grids that hold them (or, without the ends, do not span)
  Vector streams(4);
  for (Index i = 0; i < 4; i++) streams[i] = Conversion::rad2deg(std::acos(mu[i]));
  stdr::sort(streams);
  const auto with_streams = [&](Numeric first, Numeric last) {
    std::vector<Numeric> za;
    for (Numeric x = first; x <= last + 1e-9; x += 1.0) za.push_back(x);
    for (Numeric x : streams) za.push_back(x);
    stdr::sort(za);
    return Vector(za);
  };
  const Vector delta = nlinspace(-180.0, 180.0, 361);

  const auto aro_species = [&](const Vector& za_scat, Numeric extinction_change) {
    auto bulk = gas.get_bulk_scattering_properties_aro_gridded(
        atm, Vector{50e9}, streams, delta, std::make_shared<scattering::ZenithAngleGrid>(scattering::IrregularZenithAngleGrid(za_scat)));
    bulk.extinction_matrix[0, 0, 1, 0] *= 1.0 + extinction_change;
    const auto& pm = *bulk.phase_matrix;
    ARO         ssd(scattering::ParticleProperties{.name = "rayleigh", .mass = 1e-12, .d_veq = 1e-6, .d_max = 1e-6},
            pm,
            bulk.extinction_matrix,
            bulk.absorption_vector,
            scattering::BackscatterMatrixData<Numeric, scattering::Format::ARO>(pm.get_t_grid(), pm.get_f_grid(), pm.get_za_inc_grid()),
            scattering::ForwardscatterMatrixData<Numeric, scattering::Format::ARO>(pm.get_t_grid(), pm.get_f_grid(), pm.get_za_inc_grid()));
    const auto prop = ScatteringSpeciesProperty{"aro", ParticulateProperty::NumberDensity};
    ArrayOfScatteringSpecies species;
    species.add(ScatteringHabit{ParticleHabit{std::vector<ARO>{ssd}}, scattering::PSD{scattering::MonodispersePSD{prop}}, 1.0, 3.0});
    auto a  = atm;
    a[prop] = 1.0;
    return std::pair{species, a};
  };

  const auto [species, a] = aro_species(with_streams(0.0, 180.0), 0.0);
  const auto f            = vdisort::scattering_optics(species, a, 50e9, mu, mu, 2, 1e-3);
  const auto g            = vdisort::scattering_optics(gas, atm, 50e9, mu, mu, 2, 1e-3);
  Numeric    d = 0.0, scale = 0.0;
  for (Index m = 0; m < 2; m++) {
    for (Index o = 0; o < 4; o++) {
      for (Index i = 0; i < 4; i++) {
        scale = std::max(scale, max_abs(g.cosine[m, o, i]));
        d     = std::max({d, max_diff(f.cosine[m, o, i], g.cosine[m, o, i]), max_diff(f.sine[m, o, i], g.sine[m, o, i])});
      }
    }
  }
  std::cout << std::format(
      "V5 ARO data (Rayleigh's laboratory frame on 1 deg grids) vs GasScatterer, m = 0, 1, 4 streams: max |dC|, |dS| / "
      "max |C| {:.1e} (tolerance 1e-3); extinction {:.1e}\n",
      d / scale,
      std::abs(f.extinction - g.extinction) / g.extinction);
  require(d <= 1e-3 * scale and std::abs(f.extinction - g.extinction) <= 1e-14 * g.extinction,
          "V5: ARO data must give VDISORT the coefficients of the medium they tabulate");

  require_error(
      [&] {
        const auto [s2, a2] = aro_species(with_streams(1.0, 179.0), 0.0);
        (void)vdisort::scattering_optics(s2, a2, 50e9, mu, mu, 2, 1e-3);
      },
      "ARO data without all scattering zenith angles have no phase integral");
  require_error(
      [&] {
        const auto [s2, a2] = aro_species(with_streams(0.0, 180.0), 0.1);
        (void)vdisort::scattering_optics(s2, a2, 50e9, mu, mu, 2, 1e-3);
      },
      "an extinction that depends on the direction");
}
}  // namespace

int main() try {
  test_rayleigh();
  const auto c = mie();
  test_lab_frame(c);
  test_fourier_against_vector_geometry(c);
  test_path();
  test_aro_data();
  std::cout << "vdisort-arts test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

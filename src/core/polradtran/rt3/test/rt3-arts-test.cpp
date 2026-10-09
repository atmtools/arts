/* Tests of rt3::scattering_optics and rt3::problem_from_path (rt3_arts.h).

   Where the expected answers come from:
   - B1: the Legendre series of Rayleigh scattering, F11 = 3/4 (1 + x^2),
     F12 = -3/4 (1 - x^2), F22 = F11, F33 = F44 = 3/2 x, in RT3's columns
     (F11, F12, F33, F34, F22, F44): [[1, -1/2, 0, 0, 1, 0],
     [0, 0, 3/2, 0, 0, 3/2], [1/2, 1/2, 0, 0, 1/2, 0]] (rayleigh.sca of
     Evans' runtesta).  ARTS's GasScatterer gives its series in closed form,
     so it must be exact to round-off.
   - B2: Evans' runmietest series (Table 3 of Evans and Stephens 1991,
     3rdparty/polradtran/runmietest): Mie scattering at 0.951 um by a
     gamma distribution of spheres with effective radius 0.2 um, effective
     variance 0.07 and refractive index 1.44 (the L = 13 problem of Garcia
     and Siewert 1989), printed with 8 decimals.  ARTS's Mie code
     (ParticleHabit::sphere) for that population, integrated over the
     Hansen and Travis (1974) gamma distribution
     n(r) ~ r^((1 - 3 b) / b) exp(-r / (a b)), must reproduce F11, F12 and
     F33 (and F22 = F11, F44 = F33) to half a unit in the 8th decimal.  F34
     comes out with the opposite sign: for the same physical sphere ARTS's
     F34 is -1 times Evans'.  The Mie data are on the nodes of a 96-point
     Gauss-Legendre rule in cos(Theta), so their Legendre series follows
     from that rule exactly; the particles hold those series.  This pins the relative sign of F34 against a
     reference external to ARTS; see rt3_arts.h for what it implies for V.
   - B3: the path builder must reproduce the path data it is given and the
     scattering-weighted layer mean of the level series. */
#include <arts_constants.h>
#include <arts_conversions.h>
#include <legendre.h>
#include <physics_funcs.h>
#include <rt3_arts.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

namespace rt3 = polradtran::rt3;

namespace {
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

constexpr Numeric cross_section = 1e-30;  // m^2

ArrayOfScatteringSpecies rayleigh_species() {
  ArrayOfScatteringSpecies s;
  s.add(GasScatterer{scattering::ConstantGasScattering{cross_section}, scattering::RayleighGasScattering{0.0}});
  return s;
}

AtmPoint air(Numeric p, Numeric t) {
  AtmPoint a;
  a.pressure    = p;
  a.temperature = t;
  return a;
}

//! B1: Rayleigh through scattering_optics
void test_rayleigh() {
  const auto species = rayleigh_species();
  const auto atm     = air(7e4, 250.0);
  const auto sigma   = cross_section * number_density(atm.pressure, atm.temperature);

  const std::array<std::array<Numeric, 6>, 3> ref{
      {{1.0, -0.5, 0.0, 0.0, 1.0, 0.0}, {0.0, 0.0, 1.5, 0.0, 0.0, 1.5}, {0.5, 0.5, 0.0, 0.0, 0.5, 0.0}}};
  for (const Index degree : {2, 6}) {
    const auto s = rt3::scattering_optics(species, atm, 89e9, degree, 1e-12);
    require(s.legendre.shape() == (std::array<Index, 2>{degree + 1, 6}), "B1: legendre shape");
    Numeric d = 0.0;
    for (Index l = 0; l <= degree; l++)
      for (Index k = 0; k < 6; k++) d = std::max(d, std::abs(s.legendre[l, k] - (l < 3 ? ref[l][k] : 0.0)));
    const Numeric dk = std::max(std::abs(s.extinction - sigma), std::abs(s.scattering - sigma)) / sigma;
    std::cout << std::format(
        "B1 Rayleigh GasScatterer, degree {}: max |legendre - closed form| {:.1e}; extinction and scattering "
        "{:.1e}\n",
        degree,
        d,
        dk);
    require(d < 1e-14 and dk < 1e-14,
            std::format("B1: the Legendre series of ARTS's Rayleigh GasScatterer must be RT3's rayleigh.sca to "
                        "round-off (degree {}), got {:.1e} and {:.1e}",
                        degree,
                        d,
                        dk));
  }

  const auto z = rt3::scattering_optics(ArrayOfScatteringSpecies{}, atm, 89e9, 4, 1e-3);
  require(z.extinction == 0.0 and z.scattering == 0.0 and z.legendre[0, 0] == 1.0 and
              stdr::count(z.legendre | by_elem, 0.0) == 5 * 6 - 1,
          "B1: an empty species array must give no scattering and the isotropic series");
  require_error([&] { (void)rt3::scattering_optics(species, atm, 89e9, -1, 1e-3); }, "negative degree");
}

//! Table 3 of Evans and Stephens (1991): P1 = F11, P2 = F12, P3 = F33, P4 = F34
constexpr std::array<std::array<Numeric, 4>, 12> evans_table3{{
    {1.00000000, -.32071711, .71206342, -.01882245},
    {1.45529318, -.20350675, 1.76014119, -.04725108},
    {1.05402631, .24638948, 1.06682431, .00894436},
    {.39758994, .18605748, .39651104, .04505815},
    {.11659302, .07124848, .09576412, .00958275},
    {.02387477, .01700757, .01765088, .00215761},
    {.00395010, .00302534, .00261549, .00029195},
    {.00053888, .00043592, .00032713, .00003502},
    {.00006372, .00005326, .00003583, .00000337},
    {.00000667, .00000572, .00000351, .00000029},
    {.00000063, .00000055, .00000031, .00000002},
    {.00000006, .00000005, .00000003, .00000000},
}};

/* B2: the Garcia and Siewert particles as one ScatteringHabit per radius of
   an nr-point Gauss-Legendre rule on [0.002, 1.6] um, each with a
   monodisperse PSD whose number density is the quadrature weight times
   n(r).  The phase matrices are on the nodes of a Gauss-Legendre rule, whose
   projection gives each particle's Legendre series exactly. */
void test_mietest() {
  using Gridded =
      scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Gridded>;
  using Spectral =
      scattering::SingleScatteringData<Numeric, scattering::Format::TRO, scattering::Representation::Spectral>;
  constexpr Numeric wavelength = 0.951e-6, reff = 0.2, veff = 0.07;
  constexpr Index   nr = 120, nang = 96, degree = 13;
  const Numeric     frequency = Constant::speed_of_light / wavelength;

  Vector x(nr), w(nr);
  {
    const scattering::GaussLegendreQuadrature q(nr);
    for (Index i = 0; i < nr; i++) {
      x[i] = q.get_nodes()[i];
      w[i] = q.get_weights()[i];
    }
  }
  constexpr Numeric rmin = 0.002, rmax = 1.6;
  // The nodes of the projection rule, as scattering angles
  Vector xa(nang), wa(nang), angles(nang);
  Legendre::GaussLegendre(xa, wa);
  stdr::reverse(xa);  // ascending scattering angles, the same nodes as the projection
  for (Index i = 0; i < nang; i++) angles[i] = Conversion::rad2deg(std::acos(xa[i]));
  const auto za = scattering::ZenithAngleGrid{scattering::IrregularZenithAngleGrid(angles)};

  ArrayOfScatteringSpecies species;
  AtmPoint                 atm = air(1e5, 280.0);
  ComplexMatrix            index(1, 1);
  index[0, 0] = Complex{1.44, 0.0};
  // n(r) relative to its maximum at the mode r = (1 - 3 b) a, in log form
  const Numeric mode  = (1.0 - 3.0 * veff) * reff;
  const auto    log_n = [&](Numeric r) {
    return (1.0 - 3.0 * veff) / veff * std::log(r / mode) - (r - mode) / (reff * veff);
  };
  for (Index i = 0; i < nr; i++) {
    const Numeric r    = 0.5 * (rmax - rmin) * x[i] + 0.5 * (rmax + rmin);
    const auto    prop = ScatteringSpeciesProperty{std::format("mie{}", i), ParticulateProperty::NumberDensity};
    Vector        t_grid{280.0}, f_grid{frequency}, diameter{2e-6 * r};
    const auto    gridded = std::get<Gridded>(ParticleHabit::sphere(t_grid, f_grid, diameter, za, index, 1000.0)[0]);

    // a_l = 2 pi sqrt((2 l + 1) / 4 pi) sum_i w_i F(x_i) P_l(x_i), the rule exact for these data
    scattering::PhaseMatrixData<Numeric, scattering::Format::TRO, scattering::Representation::Spectral> series(
        gridded.phase_matrix->get_t_grid(), gridded.phase_matrix->get_f_grid(), degree);
    Vector p(degree + 1);
    for (Index j = 0; j < nang; j++) {
      Legendre::legendre_polynomials(p, xa[j]);
      for (Index l = 0; l <= degree; l++) {
        const Numeric f =
            2.0 * Constant::pi * std::sqrt((2.0 * static_cast<Numeric>(l) + 1.0) / (4.0 * Constant::pi)) * wa[j] * p[l];
        for (Index k = 0; k < 6; k++) series[0, 0, l, k] += f * (*gridded.phase_matrix)[0, 0, j, k];
      }
    }
    const Spectral spectral(gridded.properties,
                            series,
                            gridded.extinction_matrix.to_spectral(),
                            gridded.absorption_vector.to_spectral(),
                            gridded.backscatter_matrix,
                            gridded.forwardscatter_matrix);
    species.add(ScatteringHabit{
        ParticleHabit{std::vector<Spectral>{spectral}}, scattering::PSD{scattering::MonodispersePSD{prop}}, 1.0, 3.0});
    atm[prop] = 0.5 * (rmax - rmin) * w[i] * std::exp(log_n(r));
  }

  const auto s = rt3::scattering_optics(species, atm, frequency, degree, 1e-10);

  // RT3 columns: 0 F11, 1 F12, 2 F33, 3 F34, 4 F22, 5 F44
  Numeric same = 0.0, f34_opposite = 0.0, f34_same = 0.0, spheres = 0.0, tail = 0.0;
  for (Index l = 0; l <= degree; l++) {
    const bool tabulated = l < static_cast<Index>(evans_table3.size());
    for (Index k = 0; k < 3; k++)
      same = std::max(same, std::abs(s.legendre[l, k] - (tabulated ? evans_table3[l][k] : 0.0)));
    const Numeric p4 = tabulated ? evans_table3[l][3] : 0.0;
    f34_opposite     = std::max(f34_opposite, std::abs(s.legendre[l, 3] + p4));
    f34_same         = std::max(f34_same, std::abs(s.legendre[l, 3] - p4));
    spheres          = std::max(
        {spheres, std::abs(s.legendre[l, 4] - s.legendre[l, 0]), std::abs(s.legendre[l, 5] - s.legendre[l, 2])});
    if (not tabulated)
      for (Index k = 0; k < 6; k++) tail = std::max(tail, std::abs(s.legendre[l, k]));
  }
  constexpr Numeric tol = 6e-9;  // half a unit in the 8th decimal, plus 1e-9 for the radius quadrature
  std::cout << std::format(
      "B2 runmietest (Evans and Stephens 1991, Table 3) from ARTS's Mie code: max |F11, F12, F33 - Evans| {:.2e}, "
      "max |F34 + Evans| {:.2e} (tolerance {:.0e}); max |F34 - Evans| {:.2e}; F22 = F11, F44 = F33 to {:.0e}; "
      "degrees 12-13 {:.0e}; albedo {:.15f}\n",
      same,
      f34_opposite,
      tol,
      f34_same,
      spheres,
      tail,
      s.scattering / s.extinction);
  require(same <= tol and f34_opposite <= tol and spheres == 0.0,
          std::format("B2: ARTS's Mie series of the runmietest particles must equal Evans' Table 3 in F11, F12 and "
                      "F33, and minus Table 3 in F34, to {:.0e}; got {:.2e} and {:.2e}",
                      tol,
                      same,
                      f34_opposite));
  require(f34_same > 1e4 * tol, "B2: the F34 comparison must be able to tell the two signs apart");
}

struct path_data {
  ArrayOfPropagationPathPoint ray_path;
  ArrayOfAtmPoint             atm_path;
  ArrayOfPropmatVector        propmat;
  AscendingGrid               freq_grid{Vector{89e9}};
};

path_data make_path(const Vector& altitude, const Vector& temperature, const Vector& pressure) {
  path_data d;
  for (Index l = 0; l < static_cast<Index>(altitude.size()); l++) {
    PropagationPathPoint pp;
    pp.pos = {altitude[l], 0.0, 0.0};
    pp.los = {180.0, 0.0};
    d.ray_path.push_back(pp);
    d.atm_path.push_back(air(pressure[l], temperature[l]));
    PropmatVector pm(1);
    pm[0].A() = 1e-4 * static_cast<Numeric>(1 + l);
    d.propmat.push_back(pm);
  }
  return d;
}

/* B3: Rayleigh scattering plus a Henyey-Greenstein scatterer (g = 0.5)
   whose extinction and albedo are atmospheric fields, so that the mix, and
   with it the normalised series, changes from level to level and the layer
   series is a non-trivial scattering-weighted mean. */
void test_path() {
  auto       species = rayleigh_species();
  const auto hg_ext  = ScatteringSpeciesProperty{"hg", ParticulateProperty::Extinction};
  const auto hg_ssa  = ScatteringSpeciesProperty{"hg", ParticulateProperty::SingleScatteringAlbedo};
  species.add(HenyeyGreensteinScatterer{hg_ext, hg_ssa, 0.5});
  auto d = make_path(Vector{2000.0, 1000.0, 0.0}, Vector{240.0, 260.0, 280.0}, Vector{6e4, 8e4, 1e5});
  for (Index l = 0; l < 3; l++) {
    d.atm_path[l][hg_ext] = 1e-5 * static_cast<Numeric>(l * l);
    d.atm_path[l][hg_ssa] = 0.9 - 0.2 * static_cast<Numeric>(l);
  }

  const rt3::path_settings s{.nstokes = 4, .nmu = 6, .quad = polradtran::quadrature_type::double_gauss};
  const auto               p      = rt3::problem_from_path(d.ray_path,
                                                           d.atm_path,
                                                           d.propmat,
                                                           d.freq_grid,
                                                           0,
                                                           species,
                                                           s,
                                                           polradtran::lambertian_surface{.albedo = 0.1},
                                                           285.0,
                                                           3.0);
  const Index              degree = rt3::max_legendre_degree(6, polradtran::quadrature_type::double_gauss);
  require(p.scattering_sets.size() == 2 and p.layer_scattering_index == ArrayOfIndex({0, 1}) and
              p.scattering_sets[0].legendre.nrows() == degree + 1 and p.thermal and p.direct_flux == 0.0,
          "B3: one scattering set per layer, of RT3's maximum degree, thermal and no beam");

  Numeric dev = 0.0;
  for (Index l = 0; l < 3; l++)
    dev = std::max({dev,
                    std::abs(p.height[l] - d.ray_path[l].altitude()),
                    std::abs(p.temperature[l] - d.atm_path[l].temperature)});
  for (Index l = 0; l < 2; l++)
    dev = std::max(dev, std::abs(p.gas_extinction[l] - 1e-4 * (1.5 + static_cast<Numeric>(l))) / 1e-4);

  Numeric mix = 0.0;
  for (Index l = 0; l < 2; l++) {
    const auto a = rt3::scattering_optics(species, d.atm_path[l], 89e9, degree, 1e-3);
    const auto b = rt3::scattering_optics(species, d.atm_path[l + 1], 89e9, degree, 1e-3);
    if (l == 0) mix = std::max(mix, std::abs(a.legendre[1, 0] - b.legendre[1, 0]));
    const auto& t = p.scattering_sets[l];
    dev           = std::max({dev,
                              std::abs(t.extinction - 0.5 * (a.extinction + b.extinction)) / t.extinction,
                              std::abs(t.scattering - 0.5 * (a.scattering + b.scattering)) / t.extinction});
    for (Index i = 0; i <= degree; i++)
      for (Index k = 0; k < 6; k++)
        dev = std::max(dev,
                       std::abs(t.legendre[i, k] - (a.scattering * a.legendre[i, k] + b.scattering * b.legendre[i, k]) /
                                                       (a.scattering + b.scattering)));
  }
  require(dev < 1e-14, std::format("B3: path data and scattering-weighted layer means, deviation {:.1e}", dev));
  require(mix > 0.1, "B3: the level series must differ for the layer mean to be a test");

  // With delta_m and gauss the automatic degree reaches 2 nmu_total
  const rt3::path_settings sd{.nmu = 6, .quad = polradtran::quadrature_type::gauss, .delta_m = true};
  const auto               pd = rt3::problem_from_path(
      d.ray_path, d.atm_path, d.propmat, d.freq_grid, 0, species, sd, polradtran::lambertian_surface{}, 285.0, 3.0);
  require(pd.scattering_sets[0].legendre.nrows() ==
              std::max(rt3::max_legendre_degree(6, polradtran::quadrature_type::gauss), Index{12}) + 1,
          "B3: delta-M degree");

  const auto build = [&](const path_data& x) {
    (void)rt3::problem_from_path(
        x.ray_path, x.atm_path, x.propmat, x.freq_grid, 0, species, s, polradtran::lambertian_surface{}, 285.0, 3.0);
  };
  auto e              = d;
  e.propmat[1][0].B() = 1e-6;
  require_error([&] { build(e); }, "a polarized gas propagation matrix");
  e                    = d;
  e.ray_path[2].pos[0] = 1000.0;
  require_error([&] { build(e); }, "equal altitudes");

  // A species whose phase-function integral is 3/4 of its scattering coefficient
  ArrayOfScatteringSpecies mismatched;
  mismatched.add(
      ScatteringGeneralSpectralTRO{ScatteringGeneralSpectralTROFunc{[](const AtmPoint&, const Vector& f, Index degree) {
        ScatteringTroSpectralVector v{.phase_matrix      = SpecmatMatrix(f.size(), degree + 1, Specmat{0.0}),
                                      .extinction_matrix = PropmatVector(f.size(), Propmat{1e-4}),
                                      .absorption_vector = StokvecVector(f.size(), Stokvec{})};
        for (Size i = 0; i < f.size(); i++) (*v.phase_matrix)[i, 0][0, 0] = 0.75e-4 / std::sqrt(4.0 * Constant::pi);
        return v;
      }}});
  require_error([&] { (void)rt3::scattering_optics(mismatched, d.atm_path[0], 89e9, 0, 1e-6); },
                "a phase-function integral that misses the scattering coefficient");
  std::cout << std::format("B3 path builder: path data and layer means to {:.1e}, delta-M degree, 3 error paths\n",
                           dev);
}
}  // namespace

int main() try {
  test_rayleigh();
  test_mietest();
  test_path();
  std::cout << "rt3-arts test passed\n";
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

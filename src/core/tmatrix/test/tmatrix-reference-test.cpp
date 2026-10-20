/* The C++ T-matrix against outputs of the original Fortran, without building it

   1. The recorded outputs of the double-precision Fortran (tmatrix_fortran_cases.h), with
      the tolerances of tmatrix_comparison.h: spheroids, cylinders, Chebyshev and generalized Chebyshev particles in
      fixed orientation, and spheroids, cylinders and Chebyshev particles of all size distributions in random
      orientation.
   2. The printed outputs of the original programs (tmatrix_ampld.ref, tmatrix_amplq.ref, tmatrix_tmd.ref and
      tmatrix_tmq.ref, of the double- and extended-precision programs), every printed number.  The inputs are those
      of the programs' main routines.  The tolerance of a printed number is half a unit of its last printed digit
      plus the bound of tmatrix_comparison.h, as the extended-precision programs also store the T-matrix in single
      precision.  Not compared: the van der Mee and Hovenier test, which is not ported, and the third footer of
      tmatrix_tmd.ref, which the ARTS2 wrapper printed from the last particle's cross sections.

   Argument: the directory of the printed outputs.
*/
#include <wigner_functions.h>

#include <fstream>
#include <regex>

#include "tmatrix_comparison.h"
#include "tmatrix_fortran_cases.h"

namespace {
using namespace tmatrix;
using namespace tmatrix::comparison;

//! A printed number and half a unit of its last printed digit, of the F (0.5000) or D (0.21504D+02) format
struct Printed {
  Numeric value, half;
};

Printed printed(const std::string& token) {
  const auto d        = token.find_first_of("Dd");
  const auto mantissa = token.substr(0, d);
  const auto dot      = mantissa.find('.');
  const int  decimals = dot == std::string::npos ? 0 : static_cast<int>(mantissa.size() - dot - 1);
  const int  exponent = d == std::string::npos ? 0 : std::stoi(token.substr(d + 1));
  return {std::stod(mantissa) * std::pow(10.0, exponent), 0.5 * std::pow(10.0, exponent - decimals)};
}

std::string slurp(const std::string& path) {
  std::ifstream in(path);
  if (not in) throw std::runtime_error(std::format("cannot read {}", path));
  return {std::istreambuf_iterator<char>(in), {}};
}

std::vector<std::string> lines(const std::string& text) {
  std::vector<std::string> out;
  std::istringstream       is(text);
  for (std::string line; std::getline(is, line);) out.push_back(line);
  return out;
}

//! The printed value after a label, such as "R1="
Printed after(const std::string& line, const std::string& label) {
  const auto at = line.find(label);
  if (at == std::string::npos) throw std::runtime_error(std::format("no {} in \"{}\"", label, line));
  std::istringstream is(line.substr(at + label.size()));
  std::string        token;
  is >> token;
  return printed(token);
}

//! Mishchenko's Stokes convention to ARTS's: D M D with D = diag(1, 1, -1, 1)
Numeric arts_sign(Index a, Index b) { return (a == 2) != (b == 2) ? -1.0 : 1.0; }

/** ampld.lp.f and amplq.lp.f: a prolate spheroid of axis ratio 0.5, LAM = 2 acos(-1), MRR = 1.5, MRI = 0.02,
    DDELT = 0.001, NDGS = 2, of the volume-equivalent radius 10 (RAT = 1, ampld) or of the surface-equivalent one
    (RAT = 0.1, amplq), at THET0 = 56, THET = 65, PHI0 = 114, PHI = 128, ALPHA = 145, BETA = 52, which the programs
    move off the singular directions */
void published_fixed(const std::string& path, EquivalentRadius equivalent) {
  const auto       text = slurp(path);
  const Particle   c{.radius             = 10.0,
                     .wavelength         = 2.0 * Constant::pi,
                     .aspect_ratio       = 0.5,
                     .refractive_index   = {1.5, 0.02},
                     .shape              = {ParticleShape::spheroid},
                     .equivalent         = equivalent,
                     .accuracy           = 0.001,
                     .surface_quadrature = 2};
  const std::array geometry{fortran_geometry({.theta_incident  = 56,
                                              .theta_scattered = 65,
                                              .phi_incident    = 114,
                                              .phi_scattered   = 128,
                                              .alpha           = 145,
                                              .beta            = 52})};
  const auto       port    = fixed_port(c, geometry);
  const auto&      r       = port.results[0];
  const auto       name    = std::format("{} ({})", path.substr(path.find_last_of('/') + 1), c.describe());
  Numeric          worst_s = 0.0, worst_z = 0.0;

  std::smatch match;
  if (std::regex_search(text, match, std::regex(R"(NMAX=\s*(\d+))")))
    require(r.order == std::stoi(match[1]), std::format("{}: n_max {} against {}", name, r.order, match[1].str()));
  if (std::regex_search(text, match, std::regex(R"(RAT=\s*(\S+))"))) {
    // The volume- over the surface-equivalent radius, in closed form
    const auto rat = printed(match[1]);
    require(std::abs(c.factor() - rat.value) <= rat.half + rounding,
            std::format("{}: radius ratio {} against {}", name, c.factor(), rat.value));
  }
  const std::regex amplitude(R"(S([12])([12])=\s*(\S+) \+ i\*\s*(\S+))");
  Size             elements = 0;
  for (auto it = std::sregex_iterator(text.begin(), text.end(), amplitude); it != std::sregex_iterator(); ++it) {
    const Index   a = std::stoi((*it)[1]) - 1, b = std::stoi((*it)[2]) - 1;
    const Printed re = printed((*it)[3]), im = printed((*it)[4]);
    worst_s = std::max({worst_s,
                        std::abs(r.amplitude[a, b].real() - re.value) / (port.amplitude[0] + re.half),
                        std::abs(r.amplitude[a, b].imag() - im.value) / (port.amplitude[0] + im.half)});
    ++elements;
  }
  require(elements == 4, name + ": four amplitudes printed");
  const auto all = lines(text);
  const auto at  = std::ranges::find(all, "PHASE MATRIX");
  require(std::distance(at, all.end()) > 4, name + ": the phase matrix printed");
  for (Index a = 0; a < 4 and std::distance(at, all.end()) > 4; ++a) {
    std::istringstream row(*(at + 1 + a));
    for (Index b = 0; b < 4; ++b) {
      std::string token;
      row >> token;
      const Printed z = printed(token);
      worst_z = std::max(worst_z, std::abs(r.phase[a, b] - arts_sign(a, b) * z.value) / (port.phase[0] + z.half));
    }
  }
  require(worst_s <= 1.0, std::format("{}: amplitude at {:.2f} of its tolerance", name, worst_s));
  require(worst_z <= 1.0, std::format("{}: phase matrix at {:.2f} of its tolerance", name, worst_z));
  std::println("{}: amplitude {:.3f} and phase {:.3f} of their tolerances", name, worst_s, worst_z);
}

/** tmd.lp.f and tmq.lp.f: randomly oriented oblate spheroids of axis ratio 2, LAM = 0.5, DDELT = 0.001, NDGS = 2,
    MRR = 1.53 and MRI = 0.008 (tmd) or MRR = 1.33 and MRI = 0.001 (tmq), in the power law of the
    surface-equivalent effective radii 1 and 0.5 (AXMAX = 1, NPNAX = 2) and the variance B = 0.1, with the printed
    numbers of size quadrature points, at 19 scattering angles */
void published_random(const std::string& path, Complex m) {
  const auto all    = lines(slurp(path));
  const auto file   = path.substr(path.find_last_of('/') + 1);
  Size       blocks = 0;
  for (Size i = 0; i < all.size(); ++i) {
    if (not all[i].starts_with("R1=")) continue;
    ++blocks;
    const Printed r1 = after(all[i], "R1="), r2 = after(all[i], "R2=");
    const Printed reff = after(all[i + 1], "REFF="), veff = after(all[i + 1], "VEFF=");
    const Index   nodes = static_cast<Index>(after(all[i + 2], "AVERAGING =").value);
    Size          k     = i + 3;
    while (k < all.size() and not all[k].starts_with("CEXT=")) ++k;
    if (k == all.size()) throw std::runtime_error(std::format("{}: no cross sections after R1=", file));
    const Printed cext = after(all[k], "CEXT="), csca = after(all[k], "CSCA="), w = after(all[k], "W="),
                  cos = after(all[k], "<COS>=");

    const Particle     c{.radius             = reff.value,
                         .wavelength         = 0.5,
                         .aspect_ratio       = 2.0,
                         .refractive_index   = m,
                         .shape              = {ParticleShape::spheroid},
                         .equivalent         = EquivalentRadius::surface_area,
                         .accuracy           = 0.001,
                         .surface_quadrature = 2};
    const Distribution d{.kind            = TMatrixSizeDistribution::power_law,
                         .b               = veff.value,
                         .gamma           = 0.5,
                         .lower_ratio     = 0.9999999,
                         .upper_ratio     = 1.0000001,
                         .size_quadrature = nodes - 2};
    const auto         port = random_port(c, d, 19);
    const auto&        r    = port.result;
    const auto&        e    = port.mean.expansion;
    const auto         name = std::format("{} reff={}", file, reff.value);

    const auto check = [&name](const char* what, Numeric got, Printed ref, Numeric bound) {
      require(std::abs(got - ref.value) <= bound + ref.half,
              std::format("{}: {} {} against {}, beyond {:.1e}", name, what, got, ref.value, bound + ref.half));
    };
    // The roots of the power law's radii, to the precision of double and of the original's ZEROIN
    check("r_1", port.nodes.lower, r1, rounding * r1.value);
    check("r_2", port.nodes.upper, r2, rounding * r2.value);
    check("effective radius", r.effective_radius, reff, rounding * reff.value);
    check("effective variance", r.effective_variance, veff, variance_bound(port, veff.value));
    check("extinction", r.extinction, cext, port.cross * cext.value);
    check("scattering", r.scattering, csca, port.cross * csca.value);
    check("albedo", r.albedo, w, 2 * port.cross * w.value);
    check("asymmetry", r.asymmetry, cos, port.matrix);

    // The expansion coefficients, Mishchenko's, which the original ends at the first l with |g_1| < 1e-6
    Size l = 0, j = k + 1;
    while (j < all.size() and all[j].find("ALPHA1") == std::string::npos) ++j;
    Numeric worst = 0.0;
    for (++j; j < all.size(); ++j, ++l) {
      std::istringstream row(all[j]);
      Index              degree;
      if (not(row >> degree)) break;
      require(degree == static_cast<Index>(l) and l < e.size(), std::format("{}: coefficients of degree {}", name, l));
      if (l >= e.size()) break;
      for (const Vector* x : {&e.alpha1, &e.alpha2, &e.alpha3, &e.alpha4, &e.beta1, &e.beta2}) {
        std::string token;
        row >> token;
        const Printed p     = printed(token);
        const Numeric bound = single * std::abs((*x)[l]) + port.tail + p.half;
        worst               = std::max(worst, std::abs((*x)[l] - p.value) / bound);
      }
    }
    require(l > 2 and worst <= 1.0,
            std::format("{}: {} expansion coefficients at {:.2f} of their tolerance", name, l, worst));

    // The scattering matrix, F11, F22, F33, F44, F12, F34 as printed, F34 in Mishchenko's convention
    while (j < all.size() and all[j].find("F11") == std::string::npos) ++j;
    Numeric worst_f = 0.0;
    Size    angles  = 0;
    for (++j; j < all.size() and angles < r.phase.size(); ++j, ++angles) {
      std::istringstream row(all[j]);
      std::string        token;
      row >> token;
      const auto& F = r.phase[angles];
      for (const Numeric got : {F[0], F[2], F[3], F[5], F[1], -F[4]}) {
        row >> token;
        const Printed p = printed(token);
        worst_f         = std::max(worst_f, std::abs(got - p.value) / (port.matrix + p.half));
      }
    }
    require(angles == 19 and worst_f <= 1.0,
            std::format("{}: scattering matrix at {} angles at {:.2f} of its tolerance", name, angles, worst_f));
    std::println("{}: {} expansion coefficients at {:.3f} and scattering matrix at {:.3f} of their tolerances",
                 name,
                 l,
                 worst,
                 worst_f);
  }
  require(blocks == 2, std::format("{}: two size distributions printed", file));
}
}  // namespace

int main(int argc, char** argv) try {
  if (argc != 2) throw std::runtime_error("usage: tmatrix-reference-test <directory of the printed outputs>");
  make_wigner_ready(250, 20000000, 6);
  const std::string published = argv[1];

  for (const auto& r : fortran_fixed) {
    const auto geometries = r.geometries();
    compare("fixed " + r.particle.describe(), fixed_port(r.particle, geometries), r.results());
  }
  for (const auto& r : fortran_random)
    compare(std::format("random {} distribution={} nk={} ndgs={}",
                        r.particle.describe(),
                        r.distribution.kind,
                        r.distribution.size_quadrature,
                        r.particle.surface_quadrature),
            random_port(r.particle, r.distribution, static_cast<Index>(r.phase.size())),
            r.result());
  std::println("{} fixed and {} random recorded cases", fortran_fixed.size(), fortran_random.size());
  require(not fortran_fixed.empty() and not fortran_random.empty(), "recorded cases of both orientations");

  published_fixed(published + "/tmatrix_ampld.ref", EquivalentRadius::volume);
  published_fixed(published + "/tmatrix_amplq.ref", EquivalentRadius::surface_area);
  published_random(published + "/tmatrix_tmd.ref", {1.53, 0.008});
  published_random(published + "/tmatrix_tmq.ref", {1.33, 0.001});

  if (failures) throw std::runtime_error(std::format("{} comparisons failed", failures));
} catch (const std::exception& e) {
  std::println(stderr, "{}", e.what());
  return 1;
}

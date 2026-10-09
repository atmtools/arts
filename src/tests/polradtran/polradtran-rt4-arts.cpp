/* ARTS's RT4 (src/core/polradtran/rt4, the library with the ARTS changes) on Evans'
   two RT4 benchmark scripts, against his expected outputs.

   polradtran-rt4-arts <polradtran folder> <scatcnv program> <work directory>

   - runtestc: horizontally oriented ice columns at 340 GHz, a cirrus layer
     in a tropical atmosphere over land, 8 Lobatto streams.  The optics are
     Evans' RT4 scattering file cl340d14.dda (DDA), read as rt4.f reads it
     (GET_SCAT_FILE in radscat4.f.orig).
   - runtestr: a 2 mm/h rain layer of spherical drops at 85 GHz over water
     (Fresnel), 8 Gauss streams.  As in the script, Evans' scatcnv (the
     program scatcnv-evans) converts the Mie Legendre series testr.sca to
     the RT4 scattering file testr.rts, which is read the same way.

   So the library gets exactly the optics that rt4.f gets.  The problems
   (layers, surface, sky, wavelength) are read from the scripts.  Evans'
   tables were made with 5-digit Planck constants, which the ARTS3 RT4
   replaces by exact ones; the temperatures given to the library are those
   at which the exact Planck function equals Evans' (evans-scripts.h).  The
   output is converted as rt4.f's OUTPUT_FILE and CONVERT_OUTPUT do: fluxes
   2 pi sum_j w_j mu_j I_j on the quadrature streams, the V and H
   polarizations (I +- Q) / 2, and the effective blackbody temperature of
   2 V, 2 H (and of the flux / pi) with Evans' constants.  Every brightness
   temperature must agree with the table to one unit in its last printed
   digit (0.01 K). */
#include <arts_constants.h>
#include <rt4.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "evans-scripts.h"

namespace rt4 = polradtran::rt4;

namespace {
namespace fs = std::filesystem;
using evans::require;

//! Exact-SI radiation constants in Evans' units, 2 h c^2 [W m-2 sr-1 um^4] and h c / k [um K]
constexpr Numeric planck_c1 = 2.0 * Constant::h * Constant::c * Constant::c * 1e24;
constexpr Numeric planck_c2 = Constant::h * Constant::c / Constant::k * 1e6;

std::string read_file(const fs::path& file) {
  std::ifstream in(file);
  require(in.good(), std::format("Cannot read {}", file.string()));
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

/* An RT4 scattering file, as GET_SCAT_FILE (radscat4.f.orig) reads it, for
   nstokes <= 4 and the m = 0 mode: after the comment lines, NMU, NAZ and the
   quadrature; then for every incident hemisphere L1 and stream J1, outgoing
   hemisphere L2 and stream J2, a line MU1 MU2 M and the 4 x 4 Mueller matrix
   (row: outgoing Stokes, column: incident); then the extinction matrix and
   the emission vector of every hemisphere and stream.  Hemisphere 1 has
   mu > 0 (propagating down), RT4's L = 1, the library's rt4::down. */
rt4::layer_optics read_rt4_scattering(const std::string& text, Index nmu, char quad, Index ns) {
  std::istringstream       in(text);
  std::vector<std::string> lines;
  for (std::string line; std::getline(in, line);) {
    const auto first = line.find_first_not_of(' ');
    if (first != std::string::npos and line[first] != 'C') lines.push_back(line);  // titles: "C ..." or " C ..."
  }
  require(not lines.empty(), "Empty scattering file");
  // NMU, NAZ and the quoted quadrature name
  Index n = 0, naz = -1;
  std::istringstream(lines.front()) >> n >> naz;
  const auto        quote = lines.front().find('\'');
  const std::string q     = quote == std::string::npos ? "" : lines.front().substr(quote + 1);
  std::string       rest;
  for (std::size_t k = 1; k < lines.size(); k++) rest += lines[k] + '\n';
  std::istringstream tokens(rest);
  require(n == nmu and not q.empty() and q.front() == quad and naz == 0,
          std::format("Scattering file with nmu {}, aziorder {} and quadrature {}; expected nmu {}, aziorder 0 and {}",
                      n,
                      naz,
                      q,
                      nmu,
                      quad));

  rt4::layer_optics o{.extinction = Tensor4(2, nmu, ns, ns, 0.0),
                      .absorption = Tensor3(2, nmu, ns, 0.0),
                      .phase      = Tensor6(2, 2, nmu, nmu, ns, ns, 0.0)};
  const auto        number = [&] {
    Numeric x;
    tokens >> x;
    require(not tokens.fail(), "Scattering file ended early");
    return x;
  };
  for (Index l1 = 0; l1 < 2; l1++) {
    for (Index j1 = 0; j1 < nmu; j1++) {
      for (Index l2 = 0; l2 < 2; l2++) {
        for (Index j2 = 0; j2 < nmu; j2++) {
          number(), number();
          require(number() == 0.0, "Only the m = 0 mode is read");
          for (Index i2 = 0; i2 < 4; i2++) {
            for (Index i1 = 0; i1 < 4; i1++) {
              const Numeric x = number();
              if (i2 < ns and i1 < ns) o.phase[l2, l1, j2, j1, i2, i1] = x;
            }
          }
        }
      }
    }
  }
  for (Index l = 0; l < 2; l++) {
    for (Index j = 0; j < nmu; j++) {
      number();
      for (Index i2 = 0; i2 < 4; i2++) {
        for (Index i1 = 0; i1 < 4; i1++) {
          const Numeric x = number();
          if (i2 < ns and i1 < ns) o.extinction[l, j, i2, i1] = x;
        }
      }
    }
  }
  for (Index l = 0; l < 2; l++) {
    for (Index j = 0; j < nmu; j++) {
      number();
      for (Index i = 0; i < 4; i++) {
        const Numeric x = number();
        if (i < ns) o.absorption[l, j, i] = x;
      }
    }
  }
  return o;
}

//! Evans' CONVERT_OUTPUT for UNITS 'T': [I, Q] per micrometre to the effective blackbody temperatures of V and H
std::array<Numeric, 2> brightness_vh(Numeric i, Numeric q, Numeric lambda, bool flux) {
  std::array<Numeric, 2> t{};
  for (Index k = 0; k < 2; k++) {
    Numeric rad = 2.0 * 0.5 * (i + (k == 0 ? q : -q));
    if (flux) rad /= Constant::pi;
    const Numeric s = rad < 0 ? -1.0 : 1.0;
    t[k] = rad == 0 ? 0.0 : s * 1.4388e4 / (lambda * std::log(1.0 + 1.1911e8 / (s * rad * std::pow(lambda, 5))));
  }
  return t;
}

//! The library on a script, and its maximum deviation from the table in units of the last printed digit
Numeric run(const fs::path& folder, const std::string& name, const fs::path& scatcnv, const fs::path& work) {
  const auto s  = evans::read_script(folder / name);
  const auto st = evans::read_rt4_settings(s);
  require(st.units == 'T' and st.polarization == "VH" and st.nstokes == 2,
          std::format("{}: only EBB temperatures in V and H with nstokes 2 are handled", name));
  const auto quad = st.quad == 'L'   ? polradtran::quadrature_type::lobatto
                    : st.quad == 'G' ? polradtran::quadrature_type::gauss
                                     : polradtran::quadrature_type::double_gauss;

  // The scattering files: Evans' data files, or scatcnv run on the script's input as the script does
  fs::create_directories(work);
  for (const auto& [file, body] : s.files) std::ofstream(work / file) << body;
  for (const auto& r : s.runs) {
    if (r.program != "scatcnv") continue;
    std::ofstream answers(work / "scatcnv.in");
    for (const auto& a : r.answers) answers << a << '\n';
    answers.close();
    const auto command =
        std::format("cd \"{}\" && \"{}\" < scatcnv.in > scatcnv.log 2>&1", work.string(), scatcnv.string());
    require(std::system(command.c_str()) == 0,
            std::format("{} failed; see {}", scatcnv.string(), (work / "scatcnv.log").string()));
  }
  const auto optics_of = [&](const std::string& file) {
    const auto path = s.files.contains(file) or fs::exists(work / file) ? work / file : folder / file;
    return read_rt4_scattering(read_file(path), st.nmu, st.quad, st.nstokes);
  };

  const auto  levels = evans::read_layers(s.files.at(st.layer_file));
  const auto  T = [&](Numeric t) { return evans::exact_temperature_of_5digit(st.wavelength, t, planck_c1, planck_c2); };
  const Index nlay        = static_cast<Index>(levels.size()) - 1;
  const Numeric frequency = Constant::c / (st.wavelength * 1e-6);
  const Numeric per_um    = st.wavelength / frequency;

  rt4::problem p;
  p.nstokes             = st.nstokes;
  p.nmu                 = st.nmu;
  p.quad                = quad;
  p.max_delta_tau       = 1e-6;  // rt4.f's MAX_DELTA_TAU
  p.frequency           = frequency;
  p.height              = Vector(nlay + 1);
  p.temperature         = Vector(nlay + 1);
  p.gas_extinction      = Vector(nlay);
  p.layer_optics_index  = ArrayOfIndex(nlay, -1);
  p.sky_temperature     = T(st.sky_temperature);
  p.surface_temperature = T(st.ground_temperature);
  if (st.ground_type == 'F')
    p.ground = polradtran::fresnel_surface{.refractive_index = st.ground_index};
  else
    p.ground = polradtran::lambertian_surface{.albedo = st.albedo};
  std::map<std::string, Index> set_of;
  for (Index l = 0; l <= nlay; l++) {
    p.height[l]      = levels[l].height;
    p.temperature[l] = T(levels[l].temperature);
    if (l == nlay) break;
    p.gas_extinction[l] = levels[l].gas;
    const auto& file    = levels[l].scattering_file;
    if (file.empty()) continue;
    if (not set_of.contains(file)) {
      set_of[file] = static_cast<Index>(p.optics.size());
      p.optics.push_back(optics_of(file));
    }
    p.layer_optics_index[l] = set_of[file];
  }
  const auto r = rt4::solve(p);

  // The table's rows from the solution
  const auto table = evans::read_output(s.files.at(s.check));
  Numeric    worst = 0.0;
  for (const auto& row : table) {
    Index l = -1;
    for (Index k = 0; k <= nlay; k++)
      if (std::abs(p.height[k] - row.z) < 1e-9) l = k;
    require(l >= 0, std::format("{}: no level at Z = {}", name, row.z));
    std::array<Numeric, 2> iq{};
    bool                   flux = std::abs(row.mu) == 2.0;
    if (flux) {
      const auto& rad = row.mu < 0 ? r.up : r.down;
      for (Index j = 0; j < p.nmu; j++)
        for (Index k = 0; k < 2; k++) iq[k] += 2 * Constant::pi * r.weights[j] * r.mu[j] * rad[l, j, k];
    } else {
      Index j = -1;
      for (Index k = 0; k < static_cast<Index>(r.mu.size()); k++)
        if (std::abs(r.mu[k] - std::abs(row.mu)) < 6e-6) j = k;
      require(j >= 0, std::format("{}: no stream at MU = {}", name, row.mu));
      const auto& rad = row.mu < 0 ? r.up : r.down;
      for (Index k = 0; k < 2; k++) iq[k] = rad[l, j, k];
    }
    const auto t = brightness_vh(iq[0] / per_um, iq[1] / per_um, st.wavelength, flux);
    for (Index k = 0; k < 2; k++) worst = std::max(worst, std::abs(t[k] - row.iquv[k]) / row.unit[k]);
  }
  std::cout << std::format(
      "{}: RT4 (ARTS) with Evans' optics, {} {} streams, {} layers: {} brightness temperatures within {:.2f} unit(s) "
      "of the last printed digit (0.01 K) of his table (tolerance 1)\n",
      name,
      st.nmu,
      st.quad == 'L' ? "Lobatto" : "Gauss",
      nlay,
      2 * table.size(),
      worst);
  return worst;
}
}  // namespace

int main(int argc, char** argv) try {
  require(argc == 4, "Usage: polradtran-rt4-arts <polradtran folder> <scatcnv program> <work directory>");
  require(rt4::available(), "This test requires ENABLE_RT4=ON");
  const fs::path folder = fs::absolute(argv[1]), scatcnv = fs::absolute(argv[2]), work = fs::absolute(argv[3]);
  for (const std::string name : {"runtestc", "runtestr"})
    require(run(folder, name, scatcnv, work / name) <= 1.0 + 1e-9,
            std::format("{}: the library must reproduce Evans' table to 0.01 K", name));
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

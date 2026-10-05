#pragma once

/* Evans' PolRadTran benchmark scripts (3rdparty/polradtran: runmietest and
   runtesta for RT3, runtestr and runtestc for RT4, unchanged from the tar)
   as data, so that the tests take the problems and Evans' expected outputs
   from the scripts themselves instead of from copies.

   A script writes its input files with "cat >FILE <<EOF" blocks, runs his
   programs (rt3, rt4, scatcnv) with their answers in "PROGRAM <<EOF" blocks,
   and holds the expected output as the "cat >*.check <<EOF" block. */

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <format>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace evans {
inline void require(bool ok, const std::string& what) {
  if (not ok) throw std::runtime_error(what);
}

//! One program run of a script: the program (rt3, rt4 or scatcnv) and its answers, one per line
struct run {
  std::string              program;
  std::vector<std::string> answers;
};

struct script {
  std::string                        name;   // e.g. runtesta
  std::map<std::string, std::string> files;  // the "cat >FILE <<EOF" blocks
  std::vector<run>                   runs;   // in order
  std::string                        check;  // the name of the expected output among files

  //! The answers of the radiative transfer run (rt3 or rt4), the last run
  const run& solver() const { return runs.back(); }
};

inline script read_script(const std::filesystem::path& path) {
  std::ifstream in(path);
  require(in.good(), std::format("Cannot read {}", path.string()));
  script s;
  s.name = path.filename().string();
  for (std::string line; std::getline(in, line);) {
    const auto here = line.find(" <<EOF");
    if (here == std::string::npos) continue;
    const std::string head = line.substr(0, here);
    std::string       body;
    std::vector<std::string> answers;
    for (std::string l; std::getline(in, l) and l != "EOF";) {
      body += l + '\n';
      answers.push_back(l);
    }
    if (head.starts_with("cat >")) {
      const std::string name = head.substr(5);
      s.files[name]          = body;
      if (name.ends_with(".check")) s.check = name;
    } else {
      s.runs.push_back({.program = head, .answers = std::move(answers)});
    }
  }
  require(not s.runs.empty() and not s.check.empty(),
          std::format("{} has no program run or no .check output", path.string()));
  return s;
}

/* One line of an rt3.f or rt4.f output file: Z, PHI [deg] (rt3.f only; 0
   for rt4.f), MU, then the Stokes values (rt3.f: I, Q, U, V; rt4.f: V, H or
   I, Q).  MU = -2 / +2 are the up / down fluxes; MU < 0 is upwelling.  unit
   is one unit in the last printed digit of each value. */
struct row {
  double                z, phi, mu;
  std::size_t           n;  // number of values
  std::array<double, 4> iquv;
  std::array<double, 4> unit;
};

//! One unit in the last printed digit of a number as written, e.g. 1e-6 for ".173477E+01" and 0.01 for "242.02"
inline double last_digit(const std::string& text) {
  const auto e        = text.find_first_of("Ee");
  const auto mantissa = text.substr(0, e);
  const auto dot      = mantissa.find('.');
  const int  decimals = dot == std::string::npos ? 0 : static_cast<int>(mantissa.size() - dot - 1);
  const int  exponent = e == std::string::npos ? 0 : std::stoi(text.substr(e + 1));
  return std::pow(10.0, exponent - decimals);
}

inline std::vector<row> read_output(const std::string& text) {
  std::istringstream in(text);
  std::vector<row>   rows;
  bool               phi = false;
  for (std::string line; std::getline(in, line);) {
    if (line.empty()) continue;
    if (line.front() == 'C') {
      if (line.find(" MU ") != std::string::npos) phi = line.find(" PHI ") != std::string::npos;
      continue;
    }
    std::istringstream       is(line);
    std::vector<std::string> tokens;
    for (std::string t; is >> t;) tokens.push_back(t);
    const std::size_t first = phi ? 3 : 2;
    require(tokens.size() > first and tokens.size() <= first + 4,
            std::format("Unexpected output line \"{}\"", line));
    row r{};
    r.z   = std::stod(tokens[0]);
    r.phi = phi ? std::stod(tokens[1]) : 0.0;
    r.mu  = std::stod(tokens[first - 1]);
    r.n   = tokens.size() - first;
    for (std::size_t i = 0; i < 4; i++) {
      r.iquv[i] = i < r.n ? std::stod(tokens[first + i]) : 0.0;
      r.unit[i] = i < r.n ? last_digit(tokens[first + i]) : 0.0;
    }
    rows.push_back(r);
  }
  return rows;
}

//! A line of the layer file: the interface height, temperature, and the gas extinction and scattering file of
//! the layer below it (the last line only ends the last layer)
struct level {
  double      height, temperature, gas;
  std::string scattering_file;  // empty for a gas-only layer
};

inline std::vector<level> read_layers(const std::string& text) {
  std::istringstream in(text);
  std::vector<level> levels;
  for (std::string line; std::getline(in, line);) {
    std::istringstream is(line);
    level              l{};
    if (not(is >> l.height >> l.temperature >> l.gas)) continue;
    std::string rest;
    std::getline(is, rest);
    const auto a = rest.find('\''), b = rest.rfind('\'');
    if (a != std::string::npos and b > a) l.scattering_file = rest.substr(a + 1, b - a - 1);
    while (not l.scattering_file.empty() and l.scattering_file.back() == ' ') l.scattering_file.pop_back();
    levels.push_back(l);
  }
  return levels;
}

//! An RT3 scattering file (rt3.f, scatcnv.f): extinction and scattering per unit height, and the Legendre series
//! of (F11, F12, F33, F34, F22, F44), row l.  Leading "C" comment lines are skipped.
struct scattering_file {
  double                             extinction, scattering;
  std::vector<std::array<double, 6>> legendre;
};

inline scattering_file read_scattering(const std::string& text) {
  std::istringstream       in(text);
  std::vector<std::string> lines;
  for (std::string l; std::getline(in, l);)
    if (not l.starts_with("C")) lines.push_back(l);
  require(lines.size() >= 4, "Scattering file too short");
  scattering_file s{};
  double          albedo;
  int             degree;
  std::istringstream(lines[0]) >> s.extinction;
  std::istringstream(lines[1]) >> s.scattering;
  std::istringstream(lines[2]) >> albedo;
  std::istringstream(lines[3]) >> degree;
  require(lines.size() >= static_cast<std::size_t>(degree) + 5, "Scattering file shorter than its degree");
  for (int l = 0; l <= degree; l++) {
    std::istringstream    is(lines[4 + l]);
    int                   index;
    std::array<double, 6> c{};
    is >> index >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5];
    require(not is.fail() and index == l, std::format("Bad scattering file line \"{}\"", lines[4 + l]));
    s.legendre.push_back(c);
  }
  return s;
}

//! Reads the answers of a run in order
class answers {
  const run&  r;
  std::size_t i = 0;

 public:
  explicit answers(const run& x) : r(x) {}
  std::istringstream next() {
    require(i < r.answers.size(), std::format("{}: too few answers", r.program));
    return std::istringstream(r.answers[i++]);
  }
  std::string first() {
    std::string x;
    next() >> x;
    return x;
  }
  template <typename T> T value() {
    T x{};
    next() >> x;
    return x;
  }
};

//! The answers to rt3.f's prompts (USER_INPUT), as the scripts give them
struct rt3_settings {
  int                  nstokes, nmu;
  char                 quad;
  int                  aziorder;
  std::string          layer_file;
  bool                 delta_m;
  int                  src_code;
  double               direct_flux{0.0};  // W m-2 um-1 on the horizontal
  double               direct_mu{1.0};    // as rt3.f computes it from the zenith angle
  double               ground_temperature;
  char                 ground_type;
  double               albedo{0.0};
  std::complex<double> ground_index{1.0, 0.0};
  double               sky_temperature;
  double               wavelength;  // um
  std::string          output;
};

inline rt3_settings read_rt3_settings(const script& s) {
  require(s.solver().program == "rt3", std::format("{} does not run rt3", s.name));
  answers      a(s.solver());
  rt3_settings p{};
  p.nstokes = a.value<int>();
  p.nmu     = a.value<int>();
  p.quad    = a.first().front();
  require(p.quad != 'E', "Extra angles ('E') are not supported");
  p.aziorder   = a.value<int>();
  p.layer_file = a.first();
  p.delta_m    = a.first().front() == 'Y';
  p.src_code   = a.value<int>();
  if (p.src_code == 1 or p.src_code == 3) {
    p.direct_flux = a.value<double>();
    p.direct_mu   = std::abs(std::cos(0.017453292 * a.value<double>()));  // rt3.f's DABS(DCOS(0.017453292D0*(THETA)))
  }
  p.ground_temperature = a.value<double>();
  p.ground_type        = a.first().front();
  if (p.ground_type == 'F')
    p.ground_index = a.value<std::complex<double>>();
  else
    p.albedo = a.value<double>();
  p.sky_temperature = a.value<double>();
  p.wavelength      = a.value<double>();
  a.first();  // units
  a.first();  // output polarization
  a.first();  // number of output levels
  a.first();  // output levels
  a.first();  // number of output azimuths
  p.output = a.first();
  return p;
}

//! The answers to rt4.f's prompts (USER_INPUT), as the scripts give them
struct rt4_settings {
  int                  nstokes, nmu;
  char                 quad;
  std::string          layer_file;
  double               ground_temperature;
  char                 ground_type;
  double               albedo{0.0};
  std::complex<double> ground_index{1.0, 0.0};
  double               sky_temperature;
  double               wavelength;  // um
  char                 units;       // T: EBB brightness temperature
  std::string          polarization;
  std::string          output;
};

inline rt4_settings read_rt4_settings(const script& s) {
  require(s.solver().program == "rt4", std::format("{} does not run rt4", s.name));
  answers      a(s.solver());
  rt4_settings p{};
  p.nstokes = a.value<int>();
  p.nmu     = a.value<int>();
  p.quad    = a.first().front();
  require(p.quad != 'E', "Extra angles ('E') are not supported");
  p.layer_file         = a.first();
  p.ground_temperature = a.value<double>();
  p.ground_type        = a.first().front();
  if (p.ground_type == 'F')
    p.ground_index = a.value<std::complex<double>>();
  else
    p.albedo = a.value<double>();
  p.sky_temperature = a.value<double>();
  p.wavelength      = a.value<double>();
  p.units           = a.first().front();
  p.polarization    = a.first();
  a.first();  // number of output levels
  a.first();  // output levels
  p.output = a.first();
  return p;
}

//! The output file of a script's solver run (its last answer)
inline std::string output_file(const script& s) {
  for (auto it = s.solver().answers.rbegin(); it != s.solver().answers.rend(); ++it) {
    std::istringstream is(*it);
    std::string        x;
    if (is >> x) return x;
  }
  throw std::runtime_error(std::format("{}: no output file", s.name));
}

//! Evans' 5-digit Planck function [W m-2 sr-1 um-1] (PLANCK_FUNCTION as in the tar)
inline double planck_5digit(double lambda_um, double t) {
  return 1.1911e8 / std::pow(lambda_um, 5) / (std::exp(1.4388e4 / (lambda_um * t)) - 1.0);
}

/* The temperature at which the exact Planck function equals planck_5digit(t), given the exact
   2 h c^2 [W m-2 sr-1 um^4] and h c / k [um K].  The ARTS3 RT3 and VDISORT use exact constants, and
   Evans' tables were made with the 5-digit ones; they evaluate the Planck function only at the
   interface, surface and sky temperatures, so these temperatures reproduce his values exactly. */
inline double exact_temperature_of_5digit(double lambda_um, double t, double c1, double c2) {
  if (t == 0.0) return 0.0;
  const double b = planck_5digit(lambda_um, t);
  return c2 / (lambda_um * std::log1p(c1 / (std::pow(lambda_um, 5) * b)));
}
}  // namespace evans

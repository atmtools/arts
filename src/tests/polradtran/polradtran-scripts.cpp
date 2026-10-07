/* Evans' PolRadTran benchmark scripts with his original programs.

   polradtran-scripts <script> <work directory> <program>=<path> ...

   <script> is one of his csh scripts in 3rdparty/polradtran, unchanged:
   runmietest and runtesta (rt3), runtestr (scatcnv, then rt4) and runtestc
   (rt4).  Each <program>=<path> names the executable of a program the
   script runs: rt3-evans, rt4-evans and scatcnv-evans, built from the tar's
   sources (the .orig files and the unchanged rt3.f, rt4.f and scatcnv.f).

   A script writes its input files with "cat >FILE <<EOF" blocks, runs the
   programs with their answers in "PROGRAM <<EOF" blocks, and holds the
   expected output as the "cat >*.check <<EOF" block, which it diffs as
   text.  This runner executes those blocks in <work directory> (no csh
   needed), with the data files they name from the script's folder, as if
   the script ran in the tar's folder, and compares the output with the expected one numerically,
   because compilers print numbers differently (".500000E+01" and "    .00"
   in the tables, "0.500000E+01" and "   0.00" from gfortran):
   - the Z, PHI and MU columns must be equal;
   - every value must agree to one unit in its last printed digit, except
   - in rt3.f's output, the values that are zero by symmetry (U and V at
     PHI = 0 and 180 deg, and of the fluxes, MU = -2 and 2).  rt3.f sums the
     Fourier series in REAL*4 (OUTPUT_FILE), so in the table and here they
     are round-off of about 1e-8 of max I, which differs between compilers.
     They must stay below 1e-7 times max |I| of the table. */
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include "evans-scripts.h"

namespace {
namespace fs = std::filesystem;
using evans::require;

std::string read_file(const fs::path& file) {
  std::ifstream in(file);
  require(in.good(), std::format("Cannot read {}", file.string()));
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}
}  // namespace

int main(int argc, char** argv) try {
  require(argc >= 4, "Usage: polradtran-scripts <script> <work directory> <program>=<path> ...");
  const fs::path                  script = fs::absolute(argv[1]), work = fs::absolute(argv[2]);
  std::map<std::string, fs::path> programs;
  for (int i = 3; i < argc; i++) {
    const std::string a = argv[i];
    const auto        e = a.find('=');
    require(e != std::string::npos, std::format("Expected <program>=<path>, got \"{}\"", a));
    programs[a.substr(0, e)] = fs::absolute(a.substr(e + 1));
  }

  // Execute the script's here-documents
  const auto s      = evans::read_script(script);
  const auto output = evans::output_file(s);
  fs::create_directories(work);
  for (const auto& [name, body] : s.files) std::ofstream(work / name) << body;
  // Evans runs a script in the tar's folder: data files that the written files name in quotes (e.g. the
  // scattering file cl340d14.dda in runtestc's layer file) are taken from the script's folder
  for (const auto& [name, body] : s.files) {
    for (auto a = body.find('\''); a != std::string::npos; a = body.find('\'', a + 1)) {
      const auto b = body.find('\'', a + 1);
      if (b == std::string::npos) break;
      std::string data = body.substr(a + 1, b - a - 1);
      std::erase(data, ' ');
      a = b;
      if (data.empty() or s.files.contains(data) or not fs::exists(script.parent_path() / data)) continue;
      fs::copy_file(script.parent_path() / data, work / data, fs::copy_options::overwrite_existing);
    }
  }
  fs::remove(work / output);
  fs::current_path(work);
  for (std::size_t n = 0; n < s.runs.size(); n++) {
    const auto& r = s.runs[n];
    require(programs.contains(r.program), std::format("{} runs {}, which was not given", s.name, r.program));
    const std::string in = std::format("{}.{}.in", r.program, n), log = std::format("{}.{}.log", r.program, n);
    std::ofstream     answers(work / in);
    for (const auto& a : r.answers) answers << a << '\n';
    answers.close();
    const int status =
        std::system(std::format("\"{}\" < {} > {} 2>&1", programs.at(r.program).string(), in, log).c_str());
    require(
        status == 0,
        std::format("{} failed (status {}); see {}", programs.at(r.program).string(), status, (work / log).string()));
  }

  const auto ref = evans::read_output(s.files.at(s.check));
  const auto got = evans::read_output(read_file(work / output));
  require(ref.size() == got.size(), std::format("{} has {} lines, {} has {}", output, got.size(), s.check, ref.size()));

  const bool rt3   = s.solver().program == "rt3";
  double     max_i = 0.0;
  for (const auto& r : ref) max_i = std::max(max_i, std::abs(r.iquv[0]));

  double      worst = 0.0, worst_zero = 0.0;
  std::size_t values = 0, zeros = 0;
  for (std::size_t n = 0; n < ref.size(); n++) {
    require(ref[n].z == got[n].z and ref[n].phi == got[n].phi and ref[n].mu == got[n].mu and ref[n].n == got[n].n,
            std::format("Line {}: Z, PHI, MU {} {} {} ({} values) in the table, {} {} {} ({} values) in the output",
                        n + 1,
                        ref[n].z,
                        ref[n].phi,
                        ref[n].mu,
                        ref[n].n,
                        got[n].z,
                        got[n].phi,
                        got[n].mu,
                        got[n].n));
    const bool symmetric = rt3 and (std::fmod(ref[n].phi, 180.0) == 0.0 or std::abs(ref[n].mu) == 2.0);
    for (std::size_t c = 0; c < ref[n].n; c++) {
      if (symmetric and c >= 2) {
        worst_zero = std::max({worst_zero, std::abs(ref[n].iquv[c]), std::abs(got[n].iquv[c])});
        zeros++;
      } else {
        worst = std::max(worst, std::abs(got[n].iquv[c] - ref[n].iquv[c]) / ref[n].unit[c]);
        values++;
      }
    }
  }

  std::cout << std::format("{} ({}, {} lines): {} values within {:.2f} unit(s) of the last printed digit (tolerance 1)",
                           s.name,
                           output,
                           ref.size(),
                           values,
                           worst);
  if (rt3)
    std::cout << std::format(
        "; {} zero by symmetry, max |value| / max I {:.1e} (tolerance 1e-7)", zeros, worst_zero / max_i);
  std::cout << '\n';
  require(worst <= 1.0 + 1e-9, "The output must agree with Evans' table to one unit in the last printed digit");
  require(worst_zero <= 1e-7 * max_i, "The values that are zero by symmetry must stay below 1e-7 of max I");
  return 0;
} catch (const std::exception& e) {
  std::cerr << e.what() << '\n';
  return 1;
}

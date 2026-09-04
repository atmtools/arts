#include <arts_constants.h>
#include <disort.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <format>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include "cdisort/cdisort.h"

namespace {
constexpr Index   nquad   = 16;
constexpr Index   nmodes  = 3;
constexpr Index   nlayers = 4;
constexpr Numeric mu0     = 0.5;
constexpr Numeric beam    = Constant::pi;

using clock_type = std::chrono::steady_clock;

Numeric run_cpp_disort() {
  AscendingGrid tau{0.25, 0.5, 0.75, 1.0};
  Vector        omega(nlayers, 0.9);
  Matrix        moments(nlayers, nquad + 1);
  for (Index layer = 0; layer < nlayers; ++layer) {
    moments[layer]    = 0.0;
    moments[layer, 0] = 1.0;
    moments[layer, 2] = 0.1;
  }

  disort::main_data solver(nquad,
                           nquad,
                           nmodes,
                           std::move(tau),
                           std::move(omega),
                           std::move(moments),
                           Matrix(nmodes, nquad / 2, 0.0),
                           Matrix(nmodes, nquad / 2, 0.0),
                           Vector(nlayers, 0.0),
                           Matrix(nlayers, 0),
                           {},
                           mu0,
                           beam,
                           0.0);

  disort::u_data radiance;
  Numeric        checksum = 0.0;
  solver.u(radiance, 0.0, 0.0);
  checksum += radiance.intensities.front();
  for (const Numeric depth : solver.tau()) {
    solver.u(radiance, depth, 0.0);
    checksum += radiance.intensities.front();
  }
  return checksum;
}

Numeric run_cdisort() {
  disort_state  state{};
  disort_output output{};

  state.accur                         = 1e-6;
  state.flag.usrtau                   = FALSE;
  state.flag.usrang                   = FALSE;
  state.flag.spher                    = FALSE;
  state.flag.general_source           = FALSE;
  state.flag.output_uum               = FALSE;
  state.flag.brdf_type                = BRDF_NONE;
  state.flag.ibcnd                    = GENERAL_BC;
  state.flag.planck                   = FALSE;
  state.flag.onlyfl                   = FALSE;
  state.flag.lamber                   = TRUE;
  state.flag.quiet                    = TRUE;
  state.flag.intensity_correction     = FALSE;
  state.flag.old_intensity_correction = FALSE;
  state.nlyr                          = static_cast<int>(nlayers);
  state.nstr                          = static_cast<int>(nquad);
  state.nphase                        = state.nstr;
  state.nmom                          = state.nstr;
  state.numu                          = state.nstr;
  state.nphi                          = 1;

  c_disort_state_alloc(&state);
  c_disort_out_alloc(&state, &output);
  for (Index layer = 0; layer < nlayers; ++layer) {
    state.dtauc[layer]                            = 0.25;
    state.ssalb[layer]                            = 0.9;
    state.pmom[layer * (state.nmom_nstr + 1)]     = 1.0;
    state.pmom[2 + layer * (state.nmom_nstr + 1)] = 0.1;
  }
  state.phi[0]   = 0.0;
  state.bc.umu0  = mu0;
  state.bc.phi0  = 0.0;
  state.bc.fbeam = beam;

  c_disort(&state, &output);
  Numeric checksum = 0.0;
  for (Index level = 0; level <= nlayers; ++level)
    checksum += output.rad[level].rfldir + output.rad[level].rfldn + output.rad[level].flup;

  c_disort_out_free(&state, &output);
  c_disort_state_free(&state);
  return checksum;
}

template <typename Function>
std::vector<double> measure(const Index repetitions, Function&& function, Numeric& checksum) {
  std::vector<double> elapsed(static_cast<std::size_t>(repetitions));
  for (Index repetition = 0; repetition < repetitions; ++repetition) {
    const auto start                               = clock_type::now();
    checksum                                      += function();
    const auto stop                                = clock_type::now();
    elapsed[static_cast<std::size_t>(repetition)]  = std::chrono::duration<double, std::milli>(stop - start).count();
  }
  return elapsed;
}

void print_result(const std::string_view name, std::vector<double> elapsed) {
  stdr::sort(elapsed);
  const double median = elapsed[elapsed.size() / 2];
  const double mean   = std::reduce(elapsed.begin(), elapsed.end()) / static_cast<double>(elapsed.size());
  std::cout << std::format(
      "{:<27} median {:9.3f} ms, mean {:9.3f} ms, minimum {:9.3f} ms\n", name, median, mean, elapsed.front());
}
}  // namespace

int main(const int argc, char** argv) try {
  const Index repetitions = argc > 1 ? static_cast<Index>(std::stoll(argv[1])) : 200;
  ARTS_USER_ERROR_IF(repetitions < 1, "The repetition count must be positive, got {}", repetitions);

  Numeric    checksum = run_cpp_disort() + run_cdisort();
  const auto cpp      = measure(repetitions, run_cpp_disort, checksum);
  const auto c        = measure(repetitions, run_cdisort, checksum);

  std::cout << "DISORT quadrature-output benchmark\n"
               "16 streams, 3 Fourier modes, 4 Rayleigh-scattering layers\n"
               "Includes allocation, setup, solve, and quadrature output; excludes user-angle interpolation.\n";
  print_result("CPP-DISORT", cpp);
  print_result("cDISORT (Fortran-derived)", c);
  std::vector<double> cpp_sorted = cpp;
  std::vector<double> c_sorted   = c;
  stdr::sort(cpp_sorted);
  stdr::sort(c_sorted);
  const double ratio = cpp_sorted[cpp_sorted.size() / 2] / c_sorted[c_sorted.size() / 2];
  std::cout << std::format("CPP-DISORT / cDISORT median ratio: {:.3f} (checksum {:.9g})\n", ratio, checksum);
  return 0;
} catch (const std::exception& error) {
  std::cerr << error.what() << '\n';
  return 1;
}

#include "obsel.h"

#include <algorithm>
#include <exception>
#include <iomanip>
#include <iterator>
#include <numeric>

#include "arts_constants.h"
#include "arts_constexpr_math.h"
#include "arts_omp.h"
#include "configtypes.h"
#include "debug.h"
#include "math_funcs.h"
#include "matpack_view.h"
#include "mystring.h"
#include "rtepack.h"
#include "sorted_grid.h"
#include "sorting.h"

namespace sensor {
Index Obsel::ind(const PosLos& pl) const {
  auto ptr = std::ranges::find(poslos_grid, pl);
  return ptr == poslos_grid.end() ? dont_have
                                  : std::distance(poslos_grid.begin(), ptr);
}

Index Obsel::ind(const Numeric& f) const {
  auto elem = std::ranges::lower_bound(f_grid, f);
  return (elem == f_grid.end() or *elem != f)
             ? dont_have
             : std::distance(f_grid.begin(), elem);
}

//! Question: The sensor should be a Muelmat?  All materials can have different propagation speeds for different polarization states.
Stokvec Obsel::at(Index freq_ind, Index poslos_ind) const {
  return (freq_ind == dont_have or poslos_ind == dont_have)
             ? Stokvec{0., 0., 0., 0.}
             : (f_grid_w[freq_ind] * poslos_grid_w[poslos_ind]) * polarization;
}

Stokvec Obsel::at(const Numeric& f, const PosLos& pl) const {
  return at(ind(f), ind(pl));
}

bool all_ok(const Array<Obsel>& obsels) {
  return std::all_of(obsels.begin(), obsels.end(), [](const Obsel& obsel) {
    return obsel.ok();
  });
}

bool is_exhaustive_like(const Array<Obsel>& obsels) {
  if (not all_ok(obsels)) return false;

  if (obsels.empty()) return true;

  const auto& first         = obsels.front();
  const auto& f_grid        = first.f_grid;
  const auto& poslos_grid   = first.poslos_grid;
  const auto& f_grid_w      = first.f_grid_w;
  const auto& poslos_grid_w = first.poslos_grid_w;

  for (const auto& obsel : obsels | std::views::drop(1)) {
    if (obsel.f_grid.size() != f_grid.size() or
        obsel.poslos_grid.size() != poslos_grid.size() or
        obsel.f_grid_w.size() != f_grid_w.size() or
        obsel.poslos_grid_w.size() != poslos_grid_w.size()) {
      return false;
    }
  }

  return true;
}

AscendingGrid& collect_f_grid(AscendingGrid& f_grid,
                              const Array<Obsel>& obsels,
                              const PosLos& poslos) try {
  ARTS_ASSERT(all_ok(obsels))

  f_grid.clear();

  for (const auto& obsel : obsels) {
    if (obsel.ind(poslos) == obsel.dont_have) continue;

    for (const auto& f : obsel.f_grid) {
      f_grid.unsafe_emplace_back(f);
    }
  }

  std::sort(f_grid.unsafe_begin(), f_grid.unsafe_end());
  f_grid.unsafe_resize(
      std::distance(f_grid.unsafe_begin(),
                    std::unique(f_grid.unsafe_begin(), f_grid.unsafe_end())));

  return f_grid;
}
ARTS_METHOD_ERROR_CATCH

AscendingGrid collect_f_grid(const Array<Obsel>& obsels, const PosLos& poslos) {
  AscendingGrid f_grid;
  collect_f_grid(f_grid, obsels, poslos);
  return f_grid;
}

AscendingGrid& collect_f_grid(AscendingGrid& f_grid,
                              const Array<Obsel>& obsels) try {
  ARTS_ASSERT(all_ok(obsels))

  f_grid.clear();

  for (auto& obsel : obsels) {
    for (const auto& f : obsel.f_grid) {
      f_grid.unsafe_emplace_back(f);
    }
  }

  std::sort(f_grid.unsafe_begin(), f_grid.unsafe_end());
  f_grid.unsafe_resize(
      std::distance(f_grid.unsafe_begin(),
                    std::unique(f_grid.unsafe_begin(), f_grid.unsafe_end())));

  return f_grid;
}
ARTS_METHOD_ERROR_CATCH

AscendingGrid collect_f_grid(const Array<Obsel>& obsels) {
  AscendingGrid f_grid;
  collect_f_grid(f_grid, obsels);
  return f_grid;
}

PosLosVector& collect_poslos(PosLosVector& poslos,
                             const Array<Obsel>& obsels) try {
  ARTS_ASSERT(all_ok(obsels))

  poslos.clear();

  for (const auto& obsel : obsels) {
    for (const auto& pl : obsel.poslos_grid) {
      poslos.emplace_back(pl);
    }
  }

  std::ranges::sort(poslos);
  poslos.resize(
      std::distance(poslos.begin(), std::unique(poslos.begin(), poslos.end())));

  return poslos;
}
ARTS_METHOD_ERROR_CATCH

PosLosVector collect_poslos(const Array<Obsel>& obsels) {
  PosLosVector poslos;
  collect_poslos(poslos, obsels);
  return poslos;
}

Index max_frequency_size(const Array<Obsel>& obsels) {
  return std::transform_reduce(
      obsels.begin(),
      obsels.end(),
      Index{0},
      std::plus{},
      [](const Obsel& obsel) { return obsel.f_grid.size(); });
}

void sumup(Vector& out,
           const StokvecVector& in,
           const AscendingGrid& freqs,
           const Array<Obsel>& obsels,
           const PosLos& poslos) try {
  ARTS_ASSERT(all_ok(obsels))

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Size iel = 0; iel < obsels.size(); iel++) {
      const auto& obsel = obsels[iel];

      const Index iposlos = obsel.ind(poslos);
      if (iposlos == obsel.dont_have) continue;

      const auto first = freqs.begin();
      const auto last  = freqs.end();
      auto cur = std::ranges::lower_bound(first, last, obsel.f_grid.front());
      const auto stop =
          std::ranges::upper_bound(cur, last, obsel.f_grid.back());

      for (; cur != stop; ++cur) {
        const Index ifreq = obsel.ind(*cur);
        if (ifreq == obsel.dont_have) continue;

        const Index iv  = std::distance(first, cur);
        out[iel]       += obsel.at(ifreq, iposlos) * in[iv];
      }
    }
  } else {
    String error{};

#pragma omp parallel for
    for (Size iel = 0; iel < obsels.size(); iel++) {
      try {
        const auto& obsel = obsels[iel];

        const Index iposlos = obsel.ind(poslos);
        if (iposlos == obsel.dont_have) continue;

        const auto first = freqs.begin();
        const auto last  = freqs.end();
        auto cur = std::ranges::lower_bound(first, last, obsel.f_grid.front());
        const auto stop =
            std::ranges::upper_bound(cur, last, obsel.f_grid.back());

        for (; cur != stop; ++cur) {
          const Index ifreq = obsel.ind(*cur);
          if (ifreq == obsel.dont_have) continue;

          const Index iv  = std::distance(first, cur);
          out[iel]       += obsel.at(ifreq, iposlos) * in[iv];
        }
      } catch (std::exception& e) {
#pragma omp critical
        { error += var_string(e.what(), '\n'); }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  }
}
ARTS_METHOD_ERROR_CATCH

void sumup(Vector& out,
           Matrix& out_jac,
           const StokvecVector& in,
           const StokvecMatrix& in_jac,
           const AscendingGrid& freqs,
           const Array<Obsel>& obsels,
           const PosLos& poslos) try {
  ARTS_ASSERT(all_ok(obsels))

  const Index njac = in_jac.nrows();

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Size iel = 0; iel < obsels.size(); iel++) {
      const auto& obsel = obsels[iel];

      const Index iposlos = obsel.ind(poslos);
      if (iposlos == obsel.dont_have) continue;

      const auto first = freqs.begin();
      const auto last  = freqs.end();
      auto cur = std::ranges::lower_bound(first, last, obsel.f_grid.front());
      const auto stop =
          std::ranges::upper_bound(cur, last, obsel.f_grid.back());

      for (; cur != stop; ++cur) {
        const Index ifreq = obsel.ind(*cur);
        if (ifreq == obsel.dont_have) continue;

        const Index iv = std::distance(first, cur);
        const auto scl = obsel.at(ifreq, iposlos);

        out[iel] += scl * in[iv];
        for (Index ijac = 0; ijac < njac; ijac++) {
          out_jac(iel, ijac) += scl * in_jac(ijac, iv);
        }
      }
    }
  } else {
    String error{};

#pragma omp parallel for
    for (Size iel = 0; iel < obsels.size(); iel++) {
      try {
        const auto& obsel = obsels[iel];

        const Index iposlos = obsel.ind(poslos);
        if (iposlos == obsel.dont_have) continue;

        const auto first = freqs.begin();
        const auto last  = freqs.end();
        auto cur = std::ranges::lower_bound(first, last, obsel.f_grid.front());
        const auto stop =
            std::ranges::upper_bound(cur, last, obsel.f_grid.back());

        for (; cur != stop; ++cur) {
          const Index ifreq = obsel.ind(*cur);
          if (ifreq == obsel.dont_have) continue;

          const Index iv = std::distance(first, cur);
          const auto scl = obsel.at(ifreq, iposlos);

          out[iel] += scl * in[iv];
          for (Index ijac = 0; ijac < njac; ijac++) {
            out_jac(iel, ijac) += scl * in_jac(ijac, iv);
          }
        }
      } catch (std::exception& e) {
#pragma omp critical
        { error += var_string(e.what(), '\n'); }
      }
    }
  }
}
ARTS_METHOD_ERROR_CATCH

void exhaustive_sumup(Vector& out,
                      const StokvecVector& in,
                      const Array<Obsel>& obsels,
                      const PosLos& poslos) try {
  ARTS_ASSERT(is_exhaustive_like(obsels))

  const Index nv      = obsels.front().f_grid.size();
  const Index iposlos = obsels.front().ind(poslos);
  if (iposlos == obsels.front().dont_have) return;

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Size iel = 0; iel < obsels.size(); iel++) {
      for (Index iv = 0; iv < nv; iv++) {
        out[iel] += obsels[iel].at(iv, iposlos) * in[iv];
      }
    }
  } else {
    String error{};

#pragma omp parallel for
    for (Size iel = 0; iel < obsels.size(); iel++) {
      for (Index iv = 0; iv < nv; iv++) {
        try {
          out[iel] += obsels[iel].at(iv, iposlos) * in[iv];
        } catch (std::exception& e) {
#pragma omp critical
          { error += var_string(e.what(), '\n'); }
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  }
}
ARTS_METHOD_ERROR_CATCH

void exhaustive_sumup(Vector& out,
                      Matrix& out_jac,
                      const StokvecVector& in,
                      const StokvecMatrix& in_jac,
                      const Array<Obsel>& obsels,
                      const PosLos& poslos) try {
  ARTS_ASSERT(is_exhaustive_like(obsels))

  const Index nv      = obsels.front().f_grid.size();
  const Index njac    = in_jac.nrows();
  const Index iposlos = obsels.front().ind(poslos);
  if (iposlos == obsels.front().dont_have) return;

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Size iel = 0; iel < obsels.size(); iel++) {
      for (Index iv = 0; iv < nv; iv++) {
        const auto scl = obsels[iel].at(iv, iposlos);

        out[iel] += scl * in[iv];
        for (Index ijac = 0; ijac < njac; ijac++) {
          out_jac(iel, ijac) += scl * in_jac(ijac, iv);
        }
      }
    }
  } else {
    String error{};

#pragma omp parallel for
    for (Size iel = 0; iel < obsels.size(); iel++) {
      for (Index iv = 0; iv < nv; iv++) {
        try {
          const auto scl = obsels[iel].at(iv, iposlos);

          out[iel] += scl * in[iv];
          for (Index ijac = 0; ijac < njac; ijac++) {
            out_jac(iel, ijac) += scl * in_jac(ijac, iv);
          }
        } catch (std::exception& e) {
#pragma omp critical
          { error += var_string(e.what(), '\n'); }
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), "{}", error)
  }
}
ARTS_METHOD_ERROR_CATCH

std::ostream& operator<<(std::ostream& os, const PosLos& obsel) {
  return os << "[pos: [" << obsel.pos << "], los: [" << obsel.los << "]]";
}

std::ostream& operator<<(std::ostream& os, const Obsel& obsel) {
  os << "Obsel:\n";
  os << "  frequency grid:                 " << obsel.f_grid << '\n';
  os << "  pos-los grid:                   " << obsel.poslos_grid << '\n';
  os << "  polarization:                   " << obsel.polarization << '\n';
  os << "  frequency grid weights:         " << obsel.f_grid_w << '\n';
  os << "  pos-los grid polarized weigths: " << obsel.poslos_grid_w << '\n';
  return os;
}

std::ostream& operator<<(std::ostream& os, const Array<Obsel>& obsel) {
  for (const auto& o : obsel) {
    os << o << '\n';
  }
  return os;
}

bool Obsel::ok() const {
  return f_grid.size() == f_grid_w.size() and f_grid.size() > 0 and
         poslos_grid.size() == poslos_grid_w.size() and poslos_grid.size() > 0;
}

void Obsel::set_frequency_dirac(const Numeric& f0) {
  f_grid   = Vector{f0};
  f_grid_w = {1.};
}

void boxcar(ExhaustiveVectorView x,
            const Numeric& f0,
            const Numeric& width,
            const Numeric& sl = -0.5,
            const Numeric& su = 0.5) {
  ARTS_USER_ERROR_IF(width <= 0.0, "Width must be positive")
  ARTS_USER_ERROR_IF(su > -sl, "Upper limit must be lower than lower limit")

  const Index N = x.size();
  if (N == 1) {
    x[0] = f0;
  } else {
    const Numeric dx = 1.0 / static_cast<Numeric>(N - 1);
    x.front()        = f0 + sl * width;
    x.back()         = f0 + su * width;
    for (Index i = 1; i < N - 1; i++) {
      x[i] = std::lerp(x.front(), x.back(), static_cast<Numeric>(i) * dx);
    }
  }
}

Vector boxcar(const Numeric& f0,
              const Numeric& width,
              const Index& N,
              const Numeric& sl = -0.5,
              const Numeric& su = 0.5) {
  ARTS_USER_ERROR_IF(N < 1, "N must be greater than 0")

  Vector x(N);
  boxcar(x, f0, width, sl, su);
  return x;
}

void Obsel::set_frequency_boxcar(const Numeric& f0,
                                 const Numeric& width,
                                 const Index& N) {
  f_grid_w.resize(N);
  f_grid   = boxcar(f0, width, N);
  f_grid_w = 1.0 / static_cast<Numeric>(N);
}

void Obsel::set_frequency_boxcar(const Numeric& f0,
                                 const Numeric& width,
                                 const AscendingGrid& in_f_grid,
                                 const bool error_if_empty) {
  const auto fs = std::ranges::lower_bound(in_f_grid, f0 - 0.5 * width);
  const auto fe = std::ranges::upper_bound(in_f_grid, f0 + 0.5 * width);
  std::vector<Numeric> x(fs, fe);

  ARTS_USER_ERROR_IF(error_if_empty and x.empty(),
                     "No frequencies in the boxcar")

  f_grid = std::move(x);
  f_grid_w.resize(f_grid.size());
  f_grid_w = 1.0 / static_cast<Numeric>(f_grid.size());
}

void Obsel::normalize_frequency_weights() { f_grid_w /= sum(f_grid_w); }

Numeric gauss(Numeric f0, Numeric f, Numeric fwhm) {
  return std::exp(-4 * Constant::ln_2 * Math::pow2((f - f0) / fwhm));
}

void Obsel::set_frequency_gaussian(const Numeric& f0,
                                   const Numeric& fwhm,
                                   const Index& Nfwhm,
                                   const Index& Nhwhm) {
  ARTS_USER_ERROR_IF(
      Nfwhm < 1,
      "Nfwhm must be greater than 0 (must have at least one full width half maximum)")
  ARTS_USER_ERROR_IF(
      Nhwhm < 2,
      "dNhwhm must be greater than 1 (must have at least two point per full width half maximum)")

  const Index N      = 2 * Nfwhm * (Nhwhm - 1) + 1;
  const Numeric hwhm = 0.5 * fwhm;
  const Numeric dx   = hwhm / static_cast<Numeric>(Nhwhm - 1);

  Vector x, y;
  x.reserve(N);
  y.reserve(N);
  for (Index i = -Nfwhm; i < Nfwhm; i++) {
    const Numeric fi = f0 + static_cast<Numeric>(i) * hwhm;
    y.emplace_back(gauss(f0, x.emplace_back(fi), fwhm));
    for (Index j = 1; j < Nhwhm - 1; j++) {
      y.emplace_back(
          gauss(f0, x.emplace_back(fi + static_cast<Numeric>(j) * dx), fwhm));
    }
  }

  y.emplace_back(
      gauss(f0, x.emplace_back(static_cast<Numeric>(Nfwhm) * hwhm + f0), fwhm));

  f_grid   = std::move(x);
  f_grid_w = std::move(y);

  normalize_frequency_weights();
}

void Obsel::set_frequency_gaussian(const Numeric& f0,
                                   const Numeric& fwhm,
                                   const AscendingGrid& in_f_grid) {
  f_grid = in_f_grid;

  f_grid_w.resize(f_grid.size());
  std::transform(f_grid.begin(),
                 f_grid.end(),
                 f_grid_w.begin(),
                 [f0, fwhm](const auto& f) { return gauss(f0, f, fwhm); });

  normalize_frequency_weights();
}

Vector lochain_central_frequencies(const DescendingGrid& f0s,
                                   const String& filter) {
  const Index M = f0s.size();
  ARTS_USER_ERROR_IF(M < 1, "Must have at least one central frequency")
  ARTS_USER_ERROR_IF(
      filter.size() != static_cast<Size>(M - 1),
      "Filter must be one smaller than the number of central frequencies: \"{}\" vs {}",
      filter,
      M)

  const Index NF = [M]() {
    Index n = 1;
    for (Index i = 1; i < M; i++) n *= 2;
    return n;
  }();

  std::vector<Numeric> x(NF, f0s[0]);
  for (Index j = 1; j < M; j++) {
    const Index COUNT  = NF / (1 << j);
    Index COUNTER      = 0;
    bool add           = true;
    const Numeric sadd = filter[j - 1] == 'L' ? NAN : 1.0;
    const Numeric smin = filter[j - 1] == 'U' ? NAN : 1.0;
    for (Index iv = 0; iv < NF; iv++) {
      if (add) {
        x[iv] += f0s[j] * sadd;
      } else {
        x[iv] -= f0s[j] * smin;
      }

      COUNTER++;
      if (COUNTER == COUNT) {
        COUNTER = 0;
        add     = not add;
      }
    }
  }

  // Remove NANs
  x.erase(std::remove_if(
              x.begin(), x.end(), [](const auto& v) { return std::isnan(v); }),
          x.end());

  std::ranges::sort(x);

  return x;
}

Vector lochain_central_frequencies(const DescendingGrid& f0s) {
  const Index M = f0s.size();
  ARTS_USER_ERROR_IF(M < 1, "Must have at least one central frequency")

  const Index NF = [M]() {
    Index n = 1;
    for (Index i = 1; i < M; i++) n *= 2;
    return n;
  }();

  Vector x(NF, f0s[0]);
  for (Index j = 1; j < M; j++) {
    const Index COUNT = NF / (1 << j);
    Index COUNTER     = 0;
    bool add          = true;
    for (Index iv = 0; iv < NF; iv++) {
      if (add) {
        x[iv] += f0s[j];
      } else {
        x[iv] -= f0s[j];
      }

      COUNTER++;
      if (COUNTER == COUNT) {
        COUNTER = 0;
        add     = not add;
      }
    }
  }

  std::ranges::sort(x);

  return x;
}

void Obsel::set_frequency_lochain(const DescendingGrid& f0s,
                                  const Numeric& width,
                                  const Index& N,
                                  const String& filter,
                                  const Numeric& lower_width,
                                  const Numeric& upper_width) {
  ARTS_USER_ERROR_IF(N < 1, "Must have at least one frequency per sideband")

  const Vector F0 = filter.empty() ? lochain_central_frequencies(f0s)
                                   : lochain_central_frequencies(f0s, filter);

  Vector x(N * F0.size());
  for (Index i = 0; i < F0.size(); i++) {
    boxcar(x.slice(i * N, N), F0[i], width, lower_width, upper_width);
  }

  std::ranges::sort(x);
  f_grid = std::move(x);

  f_grid_w.resize(f_grid.size());
  f_grid_w = 1.0 / static_cast<Numeric>(f_grid.size());
}

void Obsel::set_frequency_lochain(const DescendingGrid& f0s,
                                  const Numeric& width,
                                  const AscendingGrid& in_f_grid,
                                  const String& filter,
                                  const Numeric& lower_width,
                                  const Numeric& upper_width,
                                  const bool error_if_empty) {
  const Vector F0 = filter.empty() ? lochain_central_frequencies(f0s)
                                   : lochain_central_frequencies(f0s, filter);

  std::vector<Numeric> x;
  for (double f0 : F0) {
    const auto fs =
        std::ranges::lower_bound(in_f_grid, f0 + lower_width * width);
    const auto fe =
        std::ranges::upper_bound(in_f_grid, f0 + upper_width * width);

    ARTS_USER_ERROR_IF(error_if_empty and fs == fe,
                       "No frequencies in the metmm sampler")

    x.insert(x.end(), fs, fe);
  }

  std::ranges::sort(x);
  f_grid = std::move(x);
  f_grid_w.resize(f_grid.size());
  f_grid_w = 1.0 / static_cast<Numeric>(f_grid.size());
}

constexpr Index cutoff_size(const Vector& f_grid_w_new,
                            const Numeric& cutoff,
                            const bool relative) {
  if (relative) {
    const Index N = f_grid_w_new.size();
    Index sz      = N;
    Numeric sum   = 0.0;

    // First simple sum to find the cutoff
    while (--sz >= 0) {
      sum += f_grid_w_new[sz];
      if (sum >= cutoff) break;
    }

    // Then remove all scaled weights below the cutoff
    while (sz > 1) {
      const Numeric v = f_grid_w_new[sz];
      if (v / (1.0 - sum) >= cutoff) break;
      sum += v;
      --sz;
    }

    return std::min(N, sz + 1);
  }

  // Simply keep only the weights above the cutoff by finding the partition point
  return std::distance(
      f_grid_w_new.begin(),
      std::ranges::partition_point(f_grid_w_new, Cmp::ge(cutoff)));
}

bool Obsel::cutoff_frequency_weights(const Numeric& cutoff, bool relative) {
  const Index N = f_grid.size();

  Vector f_grid_new   = f_grid;
  Vector f_grid_w_new = f_grid_w;

  // sorted from large to small weights
  bubble_sort_by(
      [&f_grid_w_new](Index i, Index j) {
        return f_grid_w_new[i] < f_grid_w_new[j];
      },
      f_grid_w_new,
      f_grid_new);

  //! Get and set the new sizes
  const Index sz = cutoff_size(f_grid_w_new, cutoff, relative);
  f_grid_w_new.resize(sz);
  f_grid_new.resize(sz);

  // sorted from small to large frequencies
  bubble_sort_by(
      [&f_grid_new](Index i, Index j) { return f_grid_new[i] > f_grid_new[j]; },
      f_grid_new,
      f_grid_w_new);

  f_grid   = std::move(f_grid_new);
  f_grid_w = std::move(f_grid_w_new);
  normalize_frequency_weights();

  return f_grid.size() != N;
}
}  // namespace sensor

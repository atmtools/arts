#include "obsel.h"

#include <algorithm>
#include <exception>
#include <numeric>
#include <stdexcept>

#include "arts_omp.h"
#include "debug.h"
#include "mystring.h"
#include "rtepack.h"
#include "sorted_grid.h"

namespace sensor {
Index Obsel::ind(const PosLos& pl) const {
  auto ptr = std::ranges::find(poslos_grid, pl);
  return ptr == poslos_grid.end() ? dont_have
                                  : std::distance(poslos_grid.begin(), ptr);
}

Index Obsel::ind(const Numeric& f) const {
  auto d = std::distance(f_grid.begin(), std::ranges::lower_bound(f_grid, f));
  return d >= f_grid.size() or f_grid[d] != f ? dont_have : d;
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

  const auto& first = obsels.front();
  const auto& f_grid = first.f_grid;
  const auto& poslos_grid = first.poslos_grid;
  const auto& f_grid_w = first.f_grid_w;
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
      const auto last = freqs.end();
      auto cur = std::ranges::lower_bound(first, last, obsel.f_grid.front());
      const auto stop =
          std::ranges::upper_bound(cur, last, obsel.f_grid.back());

      for (; cur != stop; ++cur) {
        const Index ifreq = obsel.ind(*cur);
        if (ifreq == obsel.dont_have) continue;

        const Index iv = std::distance(first, cur);
        out[iel] += obsel.at(ifreq, iposlos) * in[iv];
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
        const auto last = freqs.end();
        auto cur = std::ranges::lower_bound(first, last, obsel.f_grid.front());
        const auto stop =
            std::ranges::upper_bound(cur, last, obsel.f_grid.back());

        for (; cur != stop; ++cur) {
          const Index ifreq = obsel.ind(*cur);
          if (ifreq == obsel.dont_have) continue;

          const Index iv = std::distance(first, cur);
          out[iel] += obsel.at(ifreq, iposlos) * in[iv];
        }
      } catch (std::exception& e) {
#pragma omp critical
        { error += var_string(e.what(), '\n'); }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)
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
      const auto last = freqs.end();
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
          out_jac(ijac, iel) += scl * in_jac(ijac, iv);
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
        const auto last = freqs.end();
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
            out_jac(ijac, iel) += scl * in_jac(ijac, iv);
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

  const Index nv = obsels.front().f_grid.size();
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

    ARTS_USER_ERROR_IF(not error.empty(), error)
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

  const Index nv = obsels.front().f_grid.size();
  const Index njac = in_jac.nrows();
  const Index iposlos = obsels.front().ind(poslos);
  if (iposlos == obsels.front().dont_have) return;

  if (arts_omp_in_parallel() or arts_omp_get_max_threads() == 1) {
    for (Size iel = 0; iel < obsels.size(); iel++) {
      for (Index iv = 0; iv < nv; iv++) {
        const auto scl = obsels[iel].at(iv, iposlos);

        out[iel] += scl * in[iv];
        for (Index ijac = 0; ijac < njac; ijac++) {
          out_jac(ijac, iel) += scl * in_jac(ijac, iv);
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
            out_jac(ijac, iel) += scl * in_jac(ijac, iv);
          }
        } catch (std::exception& e) {
#pragma omp critical
          { error += var_string(e.what(), '\n'); }
        }
      }
    }

    ARTS_USER_ERROR_IF(not error.empty(), error)
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
}  // namespace sensor

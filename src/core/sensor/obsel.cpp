#include "obsel.h"

#include <algorithm>
#include <utility>

#include "compare.h"
#include "debug.h"

bool SensorKey::operator==(const SensorKey& other) const {
  return other.elem == elem and other.model == model and other.type == type;
}

namespace sensor {
std::ostream& operator<<(std::ostream& os, const PosLos& obsel) {
  return os << "[pos: [" << obsel.pos << "], los: [" << obsel.los << "]]";
}

std::ostream& operator<<(std::ostream& os, const Obsel& obsel) {
  return os << std::format("{}", obsel);
}

std::ostream& operator<<(std::ostream& os, const Array<Obsel>& obsel) {
  for (const auto& o : obsel) {
    os << o << '\n';
  }
  return os;
}

void Obsel::check() const {
  ARTS_ASSERT(f, "Must exist");

  ARTS_ASSERT(poslos, "Must exist");

  ARTS_USER_ERROR_IF(
      (w.shape() != std::array<Index, 2>{poslos->size(), f->size()}),
      R"(Weight matrix must have the same shape as the poslos times the frequency grids:
  weight_matrix.shape() = {:B,}
  poslos->size()        = {}
  f->size()             = {}
)",
      w.shape(),
      poslos->size(),
      f->size());
}

void Obsel::normalize(Stokvec pol, Numeric new_value) {
  const Numeric sum = std::transform_reduce(
      w.elem_begin(),
      w.elem_end(),
      0.0,
      std::plus<>(),
      [pol](const Stokvec& ws) -> Numeric { return pol * ws; });

  ARTS_USER_ERROR_IF(sum == 0.0, "Cannot normalize, sum is zero");

  w *= new_value / sum;
}

Numeric Obsel::sumup(const StokvecVectorView& i, Index ip) const {
  ARTS_ASSERT(i.size() == f->size(), "Bad size");
  ARTS_ASSERT(ip < poslos->size() and ip >= 0, "Bad index");

  const auto ws = w[ip];

  return std::transform_reduce(i.begin(), i.end(), ws.begin(), 0.0);
}

void Obsel::sumup(VectorView out, const StokvecMatrixView& j, Index ip) const {
  const auto ws = w[ip];

  // j is a matrix of shape JACS x f->size()

  for (Index ij = 0; ij < j.nrows(); ij++) {
    auto jac  = j[ij];
    out[ij]  += std::transform_reduce(jac.begin(), jac.end(), ws.begin(), 0.0);
  }
}

Size Obsel::flat_size(const SensorKeyType& key) const {
  switch (key) {
    using enum SensorKeyType;
    case Frequency:         return f->size();
    case PointingZenith:    return poslos->size();
    case PointingAzimuth:   return poslos->size();
    case PointingAltitude:  return poslos->size();
    case PointingLatitude:  return poslos->size();
    case PointingLongitude: return poslos->size();
    case Weights:           return f->size() * poslos->size() * 4;
  }

  std::unreachable();
}

void Obsel::flat(ExhaustiveVectorView x, const SensorKeyType& key) const {
  switch (key) {
    using enum SensorKeyType;
    case Frequency:
      ARTS_USER_ERROR_IF(x.size() != f_grid().size(),
                         "Bad size. x.size(): {}, f_grid().size(): {}",
                         x.size(),
                         f_grid().size())
      x = f_grid();
      break;
    case PointingZenith:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.los[0]; });
      break;
    case PointingAzimuth:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.los[1]; });
      break;
    case PointingAltitude:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.pos[0]; });
      break;
    case PointingLatitude:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.pos[1]; });
      break;
    case PointingLongitude:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.pos[2]; });
      break;
    case Weights: {
      ARTS_USER_ERROR_IF(w.size() * 4 != x.size(),
                         "Bad size. x.size(): {} for w.size()*4: {}",
                         x.size(),
                         w.size() * 4)
      auto X = x.reshape_as(w.nrows(), w.ncols(), 4);
      for (Index i = 0; i < X.npages(); ++i) {
        for (Index j = 0; j < X.nrows(); ++j) {
          X(i, j, 0) = w(i, j).I();
          X(i, j, 1) = w(i, j).Q();
          X(i, j, 2) = w(i, j).U();
          X(i, j, 3) = w(i, j).V();
        }
      }
    } break;
  }
}

Vector Obsel::flat(const SensorKeyType& key) const {
  Vector out(flat_size(key));
  flat(out, key);
  return out;
}
}  // namespace sensor

SensorSimulations collect_simulations(const ArrayOfSensorObsel& obsels) {
  SensorSimulations out;

  for (const auto& obsel : obsels) {
    out[obsel.f_grid_ptr()].insert(obsel.poslos_grid_ptr());
  }

  return out;
}

void make_exhaustive(ArrayOfSensorObsel& obsels) {
  const SensorSimulations simuls = collect_simulations(obsels);

  // Early return on trivial cases
  if (simuls.size() == 0) return;
  if (simuls.size() == 1 and simuls.begin()->second.size() <= 1) return;

  /*
    NOTE:

    The following code is written for a case where the f-grid and poslos-grid are
    entirerly unknown.  It is more likely that either pos-los or f-grid is known
    and the other is not.  In that case, the code can be optimized to only
    consider the unknown grid.  This is not done here for simplicity.
   */

  AscendingGrid f_grid;
  SensorPosLosVector poslos_grid;

  for (const auto& [f_g, poslos_gs] : simuls) {
    Vector newfs;
    for (auto& f : *f_g) {
      if (not std::ranges::binary_search(f_grid, f)) newfs.push_back(f);
    }
    if (newfs.size() > 0) {
      newfs.reserve(newfs.size() + f_grid.size());
      for (auto f : f_grid) newfs.push_back(f);
      std::ranges::sort(newfs);
      f_grid = newfs;
    }

    for (auto& pl_gs : poslos_gs) {
      for (auto& pl : *pl_gs) {
        if (not std::ranges::contains(poslos_grid, pl))
          poslos_grid.push_back(pl);
      }
    }
  }

  auto f_grid_ptr = std::make_shared<const AscendingGrid>(f_grid);
  auto poslos_grid_ptr =
      std::make_shared<const SensorPosLosVector>(poslos_grid);

  for (auto& obsel : obsels) {
    StokvecMatrix weights(
        f_grid.size(), poslos_grid.size(), Stokvec{0.0, 0.0, 0.0, 0.0});

    for (Index iv = 0; iv < obsel.f_grid().size(); iv++) {
      const Index ivn = std::distance(
          std::ranges::find(f_grid, obsel.f_grid()[iv]), f_grid.end());
      if (ivn == f_grid.size()) continue;

      for (Index ip = 0; ip < obsel.poslos_grid().size(); ip++) {
        const Index ipn = std::distance(
            std::ranges::find(poslos_grid, obsel.poslos_grid()[ip]),
            poslos_grid.end());
        if (ipn == poslos_grid.size()) continue;

        weights(ivn, ipn) = obsel.weight_matrix()(iv, ip);
      }
    }

    obsel = SensorObsel(f_grid_ptr, poslos_grid_ptr, std::move(weights));
  }
}

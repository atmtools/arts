#include "obsel.h"

#include <compare.h>
#include <debug.h>

#include <algorithm>
#include <utility>

bool SensorKey::operator==(const SensorKey& other) const {
  return other.sensor_elem == sensor_elem and
         other.measurement_elem == measurement_elem and other.model == model and
         other.type == type;
}

namespace sensor {
void Obsel::check() const {
  ARTS_ASSERT(f, "Must exist");

  ARTS_ASSERT(poslos, "Must exist");

  ARTS_USER_ERROR_IF(
      (w.shape() != std::array{static_cast<Index>(poslos->size()),
                               static_cast<Index>(f->size())}),
      R"(Weight matrix must have the same shape as the poslos times the frequency grids:
  weight_matrix.shape() = {:B,}
  poslos->size()        = {}
  f->size()             = {}
)",
      w.shape(),
      poslos->size(),
      f->size());
}

void Obsel::normalize(Stokvec pol) {
  const Stokvec x = sum(w);

  if (x.I() != 0.0) pol.I() = std::abs(pol.I() / x.I());
  if (x.is_polarized()) {
    const Numeric hyp = std::hypot(x.Q(), x.U(), x.V());
    pol.Q()           = std::abs(pol.Q() / hyp);
    pol.U()           = std::abs(pol.U() / hyp);
    pol.V()           = std::abs(pol.V() / hyp);
  }

  std::transform(
      w.elem_begin(), w.elem_end(), w.elem_begin(), [pol](auto& e) -> Stokvec {
        return {
            e.I() * pol.I(), e.Q() * pol.Q(), e.U() * pol.U(), e.V() * pol.V()};
      });
}

Numeric Obsel::sumup(const StokvecVectorView& i, Index ip) const {
  ARTS_ASSERT(i.size() == f->size(), "Bad size");
  ARTS_ASSERT(ip < static_cast<Index>(poslos->size()) and ip >= 0, "Bad index");

  const auto ws = w[ip];

  return std::transform_reduce(
      i.begin(), i.end(), ws.begin(), 0.0, std::plus<>(), [](auto& a, auto& b) {
        return dot(a, b);
      });
}

void Obsel::sumup(VectorView out, const StokvecMatrixView& j, Index ip) const {
  const auto ws = w[ip];

  // j is a matrix of shape JACS x f->size()

  for (Index ij = 0; ij < j.nrows(); ij++) {
    auto jac = j[ij];
    out[ij] +=
        std::transform_reduce(jac.begin(),
                              jac.end(),
                              ws.begin(),
                              0.0,
                              std::plus<>{},
                              [](auto& a, auto& b) { return dot(a, b); });
  }
}

Size Obsel::flat_size(const SensorKeyType& key) const {
  switch (key) {
    using enum SensorKeyType;
    case f:   return this->f->size();
    case za:  return poslos->size();
    case aa:  return poslos->size();
    case alt: return poslos->size();
    case lat: return poslos->size();
    case lon: return poslos->size();
  }

  std::unreachable();
}

void Obsel::flat(VectorView x, const SensorKeyType& key) const {
  switch (key) {
    using enum SensorKeyType;
    case f:
      ARTS_USER_ERROR_IF(x.size() != f_grid().size(),
                         "Bad size. x.size(): {}, f_grid().size(): {}",
                         x.size(),
                         f_grid().size())
      x = f_grid();
      break;
    case za:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.los[0]; });
      break;
    case aa:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.los[1]; });
      break;
    case alt:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.pos[0]; });
      break;
    case lat:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.pos[1]; });
      break;
    case lon:
      ARTS_USER_ERROR_IF(x.size() != poslos_grid().size(),
                         "Bad size. x.size(): {}, poslos_grid().size(): {}",
                         x.size(),
                         poslos_grid().size())
      std::transform(poslos_grid().begin(),
                     poslos_grid().end(),
                     x.begin(),
                     [](auto& poslos) { return poslos.pos[2]; });
      break;
  }
}

Vector Obsel::flat(const SensorKeyType& key) const {
  Vector out(flat_size(key));
  flat(out, key);
  return out;
}

Index Obsel::find(const Vector3& pos, const Vector2& los) const {
  const auto&& first       = poslos->begin();
  const auto&& last        = poslos->end();
  const auto same_poslos = [&](const auto& p) {
    return &pos == &p.pos and & los == &p.los;
  };

  const auto&& it = std::find_if(first, last, same_poslos);

  return (it == last) ? dont_have : std::distance(first, it);
}

Index Obsel::find(const AscendingGrid& frequency_grid) const {
  return (f.get() != &frequency_grid) ? dont_have : 0;
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

    for (Size iv = 0; iv < obsel.f_grid().size(); iv++) {
      const Index ivn = std::distance(
          std::ranges::find(f_grid, obsel.f_grid()[iv]), f_grid.end());
      if (ivn == static_cast<Index>(f_grid.size())) continue;

      for (Size ip = 0; ip < obsel.poslos_grid().size(); ip++) {
        const Index ipn = std::distance(
            std::ranges::find(poslos_grid, obsel.poslos_grid()[ip]),
            poslos_grid.end());
        if (ipn == static_cast<Index>(poslos_grid.size())) continue;

        weights[ivn, ipn] = obsel.weight_matrix()[iv, ip];
      }
    }

    obsel = SensorObsel(f_grid_ptr, poslos_grid_ptr, std::move(weights));
  }
}

namespace {
void set_frq(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ConstVectorView x) {
  ARTS_USER_ERROR_IF(x.size() != v.f_grid().size(),
                     "Bad size. x.size(): {}, f_grid().size(): {}",
                     x.size(),
                     v.f_grid().size())

  const auto xs = std::make_shared<const AscendingGrid>(
      x.begin(), x.end(), [](auto& x) { return x; });

  // Must copy, as we may change the shared_ptr later
  const auto fs = v.f_grid_ptr();

  for (auto& elem : sensor) {
    if (elem.f_grid_ptr() == fs) {
      elem.set_f_grid_ptr(xs);  // may change here
    }
  }
}

template <bool pos, Index k>
void set_poslos(const SensorObsel& v,
                ArrayOfSensorObsel& sensor,
                const ConstVectorView x) {
  ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                     "Bad size. x.size(): {}, poslos_grid().size(): {}",
                     x.size(),
                     v.poslos_grid().size())

  SensorPosLosVector xsv = v.poslos_grid();

  std::transform(xsv.begin(),
                 xsv.end(),
                 x.begin(),
                 xsv.begin(),
                 [](auto poslos, Numeric val) {
                   if constexpr (pos) {
                     poslos.pos[k] = val;
                   } else {
                     poslos.los[k] = val;
                   }
                   return poslos;
                 });

  const auto xs = std::make_shared<const SensorPosLosVector>(std::move(xsv));

  // Must copy, as we may change the shared_ptr later
  const auto ps = v.poslos_grid_ptr();

  for (auto& elem : sensor) {
    if (elem.poslos_grid_ptr() == ps) {
      elem.set_poslos_grid_ptr(ps);  // may change here
    }
  }
}

void set_alt(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ConstVectorView x) {
  set_poslos<true, 0>(v, sensor, x);
}

void set_lat(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ConstVectorView x) {
  set_poslos<true, 1>(v, sensor, x);
}

void set_lon(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ConstVectorView x) {
  set_poslos<true, 2>(v, sensor, x);
}

void set_zag(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ConstVectorView x) {
  set_poslos<false, 0>(v, sensor, x);
}

void set_aag(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ConstVectorView x) {
  set_poslos<false, 1>(v, sensor, x);
}
}  // namespace

void unflatten(ArrayOfSensorObsel& sensor,
               const ConstVectorView& x,
               const SensorObsel& v,
               const SensorKeyType& key) {
  switch (key) {
    using enum SensorKeyType;
    case f:   set_frq(v, sensor, x); break;
    case za:  set_zag(v, sensor, x); break;
    case aa:  set_aag(v, sensor, x); break;
    case alt: set_alt(v, sensor, x); break;
    case lat: set_lat(v, sensor, x); break;
    case lon: set_lon(v, sensor, x); break;
  }
}
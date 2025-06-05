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
bool SparseStokvec::operator==(const SparseStokvec& other) const {
  return (irow == other.irow) and (icol == other.icol);
}

bool SparseStokvec::operator!=(const SparseStokvec& other) const {
  return not(*this == other);
}

bool SparseStokvec::operator<(const SparseStokvec& other) const {
  if (irow < other.irow) return true;
  if (irow > other.irow) return false;
  return icol < other.icol;
}

bool SparseStokvec::operator<=(const SparseStokvec& other) const {
  return (*this < other) or (*this == other);
}

bool SparseStokvec::operator>(const SparseStokvec& other) const {
  if (irow > other.irow) return true;
  if (irow < other.irow) return false;
  return icol > other.icol;
}

bool SparseStokvec::operator>=(const SparseStokvec& other) const {
  return (*this > other) or (*this == other);
}

bool SparseStokvecMatrix::operator==(const SparseStokvecMatrix& other) const {
  return rows == other.rows and cols == other.cols and
         sparse_data == other.sparse_data;
}

bool SparseStokvecMatrix::operator!=(const SparseStokvecMatrix& other) const {
  return not(*this == other);
}

SparseStokvecMatrix::SparseStokvecMatrix(const StokvecMatrix& m)
    : rows(m.nrows()), cols(m.ncols()) {
  //! NOTE: Retains sorting
  for (Size i = 0; i < rows; ++i) {
    for (Size j = 0; j < cols; ++j) {
      if (m[i, j].is_zero()) continue;
      sparse_data.push_back({i, j, m[i, j]});
    }
  }
  assert(stdr::is_sorted(sparse_data));
}

SparseStokvecMatrix& SparseStokvecMatrix::operator=(const StokvecMatrix& m) {
  return *this = SparseStokvecMatrix(m);
}

Index SparseStokvecMatrix::nrows() const { return rows; }

Index SparseStokvecMatrix::ncols() const { return cols; }

Size SparseStokvecMatrix::size() const { return sparse_data.size(); }

bool SparseStokvecMatrix::empty() const { return sparse_data.empty(); }

void SparseStokvecMatrix::resize(Size rows, Size cols, Size reserve) {
  this->rows = rows;
  this->cols = cols;
  sparse_data.clear();
  sparse_data.reserve(reserve);
}

std::array<Index, 2> SparseStokvecMatrix::shape() const {
  return {static_cast<Index>(rows), static_cast<Index>(cols)};
}

Stokvec& SparseStokvecMatrix::operator[](Size i, Size j) {
  const SparseStokvec newdata{.irow = i, .icol = j};
  auto it = stdr::lower_bound(sparse_data, newdata);

  if (it != sparse_data.end() and newdata == *it) {
    return it->data;
  }

  auto newit = sparse_data.insert(it, newdata);

  //! NOTE: Must retain sorting?
  assert(stdr::is_sorted(sparse_data));

  return newit->data;
}

Stokvec SparseStokvecMatrix::operator[](Size i, Size j) const {
  const SparseStokvec newdata{.irow = i, .icol = j};
  auto it = stdr::lower_bound(sparse_data, newdata);

  if (it != sparse_data.end() and newdata == *it) {
    return it->data;
  }

  return {};
}

std::vector<SparseStokvec>::iterator SparseStokvecMatrix::begin() {
  return sparse_data.begin();
}

std::vector<SparseStokvec>::iterator SparseStokvecMatrix::end() {
  return sparse_data.end();
}

std::vector<SparseStokvec>::const_iterator SparseStokvecMatrix::begin() const {
  return sparse_data.begin();
}

std::vector<SparseStokvec>::const_iterator SparseStokvecMatrix::end() const {
  return sparse_data.end();
}

SparseStokvecMatrix::operator StokvecMatrix() const {
  StokvecMatrix m(rows, cols, Stokvec{0.0, 0.0, 0.0, 0.0});
  for (const auto& elem : sparse_data) {
    m[elem.irow, elem.icol] = elem.data;
  }
  return m;
}

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

Obsel::Obsel(std::shared_ptr<const AscendingGrid> fs,
             std::shared_ptr<const PosLosVector> pl,
             SparseStokvecMatrix ws)
    : f{std::move(fs)}, poslos{std::move(pl)}, w{std::move(ws)} {
  check();
}

Obsel::Obsel(const AscendingGrid& fs,
             const PosLosVector& pl,
             SparseStokvecMatrix ws)
    : Obsel::Obsel(std::make_shared<const AscendingGrid>(fs),
                   std::make_shared<const PosLosVector>(pl),
                   std::move(ws)) {}

bool Obsel::same_freqs(const Obsel& other) const { return f == other.f; }

bool Obsel::same_freqs(
    const std::shared_ptr<const AscendingGrid>& other) const {
  return f == other;
}

bool Obsel::same_poslos(const Obsel& other) const {
  return poslos == other.poslos;
}

bool Obsel::same_poslos(
    const std::shared_ptr<const PosLosVector>& other) const {
  return poslos == other;
}

void Obsel::set_f_grid_ptr(std::shared_ptr<const AscendingGrid> n) {
  ARTS_USER_ERROR_IF(not n, "Must exist");
  ARTS_USER_ERROR_IF(n->size() != f->size(), "Mismatching size");
  f = std::move(n);
}

void Obsel::set_poslos_grid_ptr(std::shared_ptr<const PosLosVector> n) {
  ARTS_USER_ERROR_IF(not n, "Must exist");
  ARTS_USER_ERROR_IF(n->size() != poslos->size(), "Mismatching size");
  poslos = std::move(n);
}

void Obsel::set_weight_matrix(SparseStokvecMatrix n) {
  ARTS_USER_ERROR_IF(n.shape() != w.shape(), "Mismatching shape");
  w = std::move(n);
}

void Obsel::normalize(Stokvec pol) {
  const Stokvec x = std::transform_reduce(
      w.begin(), w.end(), Stokvec{}, std::plus<>(), [](const SparseStokvec& e) {
        return e.data;
      });

  if (x.I() != 0.0) pol.I() = std::abs(pol.I() / x.I());

  if (x.is_polarized()) {
    const Numeric hyp = std::hypot(x.Q(), x.U(), x.V());
    pol.Q()           = std::abs(pol.Q() / hyp);
    pol.U()           = std::abs(pol.U() / hyp);
    pol.V()           = std::abs(pol.V() / hyp);
  }

  for (auto& e : w) {
    e.data.I() *= pol.I();
    e.data.Q() *= pol.Q();
    e.data.U() *= pol.U();
    e.data.V() *= pol.V();
  }
}

Numeric Obsel::sumup(const StokvecVectorView& i, Index ip) const {
  ARTS_ASSERT(i.size() == f->size(), "Bad size");
  ARTS_ASSERT(ip < static_cast<Index>(poslos->size()) and ip >= 0, "Bad index");

  // w is a sparse sorted matrix of shape poslos->size() x f->size()
  const SparseStokvec v0{.irow = static_cast<Size>(ip),
                         .icol = static_cast<Size>(0)};
  const SparseStokvec vn{.irow = static_cast<Size>(ip),
                         .icol = std::numeric_limits<Size>::max()};
  std::span<const SparseStokvec> span{stdr::lower_bound(w, v0),
                                      stdr::upper_bound(w, vn)};

  Numeric sum = 0.0;
  for (const auto& ws : span) {
    ARTS_ASSERT(ws.icol < i.size(), "Bad index in sparse matrix");
    sum += dot(i[ws.icol], ws.data);
  }
  return sum;
}

void Obsel::sumup(VectorView out, const StokvecMatrixView& j, Index ip) const {
  // j is a matrix of shape JACS x f->size()
  // w is a sparse sorted matrix of shape poslos->size() x f->size()
  const SparseStokvec v0{.irow = static_cast<Size>(ip),
                         .icol = static_cast<Size>(0)};
  const SparseStokvec vn{.irow = static_cast<Size>(ip),
                         .icol = std::numeric_limits<Size>::max()};
  std::span<const SparseStokvec> span{stdr::lower_bound(w, v0),
                                      stdr::upper_bound(w, vn)};

  for (Index ij = 0; ij < j.nrows(); ij++) {
    auto jac = j[ij];
    for (auto& ws : span) {
      out[ij] += dot(jac[ws.icol], ws.data);
    }
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
  const auto&& first     = poslos->begin();
  const auto&& last      = poslos->end();
  const auto same_poslos = [&](const auto& p) {
    return &pos == &p.pos and &los == &p.los;
  };

  const auto&& it = std::find_if(first, last, same_poslos);

  return (it == last) ? dont_have : std::distance(first, it);
}

Index Obsel::find(const AscendingGrid& frequency_grid) const {
  return (f.get() != &frequency_grid) ? dont_have : 0;
}
}  // namespace sensor

SensorSimulations collect_simulations(
    const std::span<const SensorObsel>& obsels) {
  SensorSimulations out;

  for (const auto& obsel : obsels) {
    out[obsel.f_grid_ptr()].insert(obsel.poslos_grid_ptr());
  }

  return out;
}

void make_exhaustive(std::span<SensorObsel> obsels) {
  const SensorSimulations simuls = collect_simulations(obsels);

  // Early return on trivial cases
  if (simuls.size() == 0) return;
  if (simuls.size() == 1 and simuls.begin()->second.size() <= 1) return;

  std::vector<Numeric> all_freq;  // Collect all frequencies in the simulations
  std::vector<SensorPosLos> all_geom;  // Collect all poslos in the simulations

  for (auto& [freqs, geoms] : simuls) {
    all_freq.insert(all_freq.end(), freqs->begin(), freqs->end());

    for (auto& geom : geoms) {
      for (auto& poslos : *geom) {
        if (stdr::find(all_geom, poslos) == all_geom.end()) {
          all_geom.push_back(poslos);
        }
      }
    }
  }

  stdr::sort(all_freq);
  auto [first, last] = stdr::unique(all_freq);
  all_freq.erase(first, last);

  auto fptr            = std::make_shared<const AscendingGrid>(all_freq);
  auto poslos_grid_ptr = std::make_shared<const SensorPosLosVector>(all_geom);

  for (auto& elem : obsels) {
    sensor::SparseStokvecMatrix weights(all_geom.size(), all_freq.size());

    for (auto& w : elem.weight_matrix()) {
      const Size ig  = w.irow;  // Index in poslos_grid
      const Size iv  = w.icol;  // Index in f_grid
      const Size ign = stdr::distance(
          all_geom.begin(), stdr::find(all_geom, elem.poslos_grid()[ig]));
      const Size ivn = stdr::distance(
          all_freq.begin(), stdr::lower_bound(all_freq, elem.f_grid()[iv]));
      weights[ign, ivn] = w.data;
    }

    elem = SensorObsel(fptr, poslos_grid_ptr, std::move(weights));
  }
}

void make_exclusive(std::span<SensorObsel> obsels) {
  Vector fs{};
  SensorPosLosVector gs{};
  std::vector<Size> freq_indices{};
  std::vector<Size> poslos_indices{};

  for (auto& elem : obsels) {
    freq_indices.resize(0);
    poslos_indices.resize(0);

    for (auto& w : elem.weight_matrix()) {
      poslos_indices.push_back(w.irow);
      freq_indices.push_back(w.icol);
    }

    stdr::sort(poslos_indices);
    stdr::sort(freq_indices);

    {
      auto [first, last] = stdr::unique(poslos_indices);
      poslos_indices.erase(first, last);
    }
    {
      auto [first, last] = stdr::unique(freq_indices);
      freq_indices.erase(first, last);
    }

    const auto nposlos = static_cast<Index>(poslos_indices.size());
    gs.resize(nposlos);
    auto& old_gs = elem.poslos_grid();
    for (Index i = 0; i < nposlos; i++) {
      gs[i] = old_gs[poslos_indices[i]];
    }

    const auto nfreqs = static_cast<Index>(freq_indices.size());
    fs.resize(nfreqs);
    auto& old_fs = elem.f_grid();
    for (Index j = 0; j < nfreqs; j++) {
      fs[j] = fs[freq_indices[j]];
    }

    sensor::SparseStokvecMatrix ws(nposlos, nfreqs);
    for (const auto& w : elem.weight_matrix()) {
      const Index iposlos =
          stdr::distance(stdr::begin(gs), stdr::find(gs, old_gs[w.irow]));
      const Index ifreq  = stdr::distance(stdr::begin(fs),
                                         stdr::lower_bound(fs, old_fs[w.icol]));
      ws[iposlos, ifreq] = w.data;
    }

    elem = SensorObsel(AscendingGrid{fs}, gs, std::move(ws));
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
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
Obsel::Obsel()                            = default;
Obsel::Obsel(const Obsel&)                = default;
Obsel::Obsel(Obsel&&) noexcept            = default;
Obsel& Obsel::operator=(const Obsel&)     = default;
Obsel& Obsel::operator=(Obsel&&) noexcept = default;

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

const std::vector<SparseStokvec>& SparseStokvecMatrix::vector() const {
  return sparse_data;
}

void SparseStokvecMatrix::resize(Size nrows, Size ncols, Size reserve) {
  rows = nrows;
  cols = ncols;
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
  assert(f);

  assert(poslos);

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
  assert(i.size() == f->size());
  assert(ip < static_cast<Index>(poslos->size()) and ip >= 0);

  // w is a sparse sorted matrix of shape poslos->size() x f->size()
  const SparseStokvec v0{.irow = static_cast<Size>(ip),
                         .icol = static_cast<Size>(0)};
  const SparseStokvec vn{.irow = static_cast<Size>(ip),
                         .icol = std::numeric_limits<Size>::max()};
  std::span<const SparseStokvec> span{stdr::lower_bound(w, v0),
                                      stdr::upper_bound(w, vn)};

  Numeric sum = 0.0;
  for (const auto& ws : span) {
    assert(ws.icol < i.size());
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

void xml_io_stream<SensorPosLos>::write(std::ostream& os,
                                        const SensorPosLos& x,
                                        bofstream* pbofs,
                                        std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs) {
    put(&x, pbofs);
  } else {
    xml_write_to_stream(os, x.pos, pbofs, "POS");
    xml_write_to_stream(os, x.los, pbofs, "LOS");
  }

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SensorPosLos>::read(std::istream& is,
                                       SensorPosLos& x,
                                       bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&x, pbifs);
  } else {
    xml_read_from_stream(is, x.pos);
    xml_read_from_stream(is, x.los);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SensorPosLos>::put(const SensorPosLos* const x,
                                      bofstream* pbofs,
                                      Size n) {
  assert(pbofs);
  xml_io_stream<Numeric>::put(
      reinterpret_cast<const Numeric*>(x), pbofs, 5 * n);
}

void xml_io_stream<SensorPosLos>::get(SensorPosLos* x,
                                      bifstream* pbifs,
                                      Size n) {
  assert(pbifs);
  xml_io_stream<Numeric>::get(reinterpret_cast<Numeric*>(x), pbifs, 5 * n);
}

void xml_io_stream<SensorKey>::write(std::ostream& os,
                                     const SensorKey& x,
                                     bofstream* pbofs,
                                     std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.sensor_elem, pbofs);
  xml_write_to_stream(os, x.measurement_elem, pbofs);
  xml_write_to_stream(os, x.model, pbofs);
  xml_write_to_stream(os, x.polyorder, pbofs);
  xml_write_to_stream(os, x.original_grid, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SensorKey>::read(std::istream& is,
                                    SensorKey& x,
                                    bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.sensor_elem, pbifs);
  xml_read_from_stream(is, x.measurement_elem, pbifs);
  xml_read_from_stream(is, x.model, pbifs);
  xml_read_from_stream(is, x.polyorder, pbifs);
  xml_read_from_stream(is, x.original_grid, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<sensor::SparseStokvec>::put(
    const sensor::SparseStokvec* const x, bofstream* pbofs, Size n) {
  pbofs->putRaw(reinterpret_cast<const char*>(x),
                n * 2 * sizeof(Size) + sizeof(Stokvec));
}

void xml_io_stream<sensor::SparseStokvec>::write(std::ostream& os,
                                                 const sensor::SparseStokvec& x,
                                                 bofstream* pbofs,
                                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  if (pbofs) {
    put(&x, pbofs);
  } else {
    xml_write_to_stream(os, x.irow);
    xml_write_to_stream(os, x.icol);
    xml_write_to_stream(os, x.data);
  }

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<sensor::SparseStokvec>::get(sensor::SparseStokvec* x,
                                               bifstream* pbifs,
                                               Size n) {
  pbifs->getRaw(reinterpret_cast<char*>(x),
                n * 2 * sizeof(Size) + sizeof(Stokvec));
}

void xml_io_stream<sensor::SparseStokvec>::read(std::istream& is,
                                                sensor::SparseStokvec& x,
                                                bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get(&x, pbifs);
  } else {
    xml_read_from_stream(is, x.irow);
    xml_read_from_stream(is, x.icol);
    xml_read_from_stream(is, x.data);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<sensor::SparseStokvecMatrix>::write(
    std::ostream& os,
    const sensor::SparseStokvecMatrix& x,
    bofstream* pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.nrows(), pbofs);
  xml_write_to_stream(os, x.ncols(), pbofs);
  xml_write_to_stream(os, x.vector(), pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<sensor::SparseStokvecMatrix>::read(
    std::istream& is, sensor::SparseStokvecMatrix& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  Size nr, nc;
  std::vector<sensor::SparseStokvec> vector;

  xml_read_from_stream(is, nr, pbifs);
  xml_read_from_stream(is, nc, pbifs);
  xml_read_from_stream(is, vector, pbifs);

  x.resize(nr, nc, vector.size());

  stdr::move(stdv::as_rvalue(vector), x.begin());

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SensorObsel>::write(std::ostream& os,
                                       const SensorObsel& g,
                                       bofstream* pbofs,
                                       std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, g.f_grid(), pbofs, "f_grid");
  xml_write_to_stream(os, g.poslos_grid(), pbofs, "poslos");
  xml_write_to_stream(os, g.weight_matrix(), pbofs, "weight_matrix");

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<SensorObsel>::read(std::istream& is,
                                      SensorObsel& x,
                                      bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  AscendingGrid f_grid;
  SensorPosLosVector poslos_grid;
  sensor::SparseStokvecMatrix weight_matrix;

  xml_read_from_stream(is, f_grid, pbifs);
  xml_read_from_stream(is, poslos_grid, pbifs);
  xml_read_from_stream(is, weight_matrix, pbifs);

  x = SensorObsel{f_grid, poslos_grid, weight_matrix};

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<ArrayOfSensorObsel>::write(std::ostream& os,
                                              const ArrayOfSensorObsel& g,
                                              bofstream* pbofs,
                                              std::string_view name) {
  const auto sen = collect_simulations(g);

  std::vector<std::shared_ptr<const SensorPosLosVector>> plos;
  plos.reserve(sen.size());
  for (const auto& i : sen | stdv::values | stdv::join) {
    if (not stdr::contains(plos, i)) {
      plos.push_back(i);
    }
  }

  std::vector<std::shared_ptr<const AscendingGrid>> freqs;
  freqs.reserve(sen.size());
  for (const auto& i : sen | stdv::keys) freqs.push_back(i);

  std::println(os,
               R"(<{0} name="{1}" nelem="{2}" nfreq="{3}" nposlos="{4}">)",
               type_name,
               name,
               static_cast<Index>(g.size()),
               static_cast<Index>(freqs.size()),
               static_cast<Index>(plos.size()));

  if (not sen.empty()) {
    for (auto& f : freqs) xml_write_to_stream(os, *f, pbofs, "f_grid");
    for (auto& p : plos) xml_write_to_stream(os, *p, pbofs, "poslos");

    for (auto& elem : g) {
      const Index ifreq =
          std::distance(freqs.begin(), stdr::find(freqs, elem.f_grid_ptr()));
      const Index iplos =
          std::distance(plos.begin(), stdr::find(plos, elem.poslos_grid_ptr()));

      xml_write_to_stream(os, ifreq, pbofs, "f_grid index");
      xml_write_to_stream(os, iplos, pbofs, "poslos_grid index");
      xml_write_to_stream(os, elem.weight_matrix(), pbofs, "weight_matrix");
    }
  }

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<ArrayOfSensorObsel>::read(std::istream& is,
                                             ArrayOfSensorObsel& g,
                                             bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  g.resize(0);

  Index nelem, nfreq, nposlos;

  tag.get_attribute_value("nelem", nelem);
  tag.get_attribute_value("nfreq", nfreq);
  tag.get_attribute_value("nposlos", nposlos);

  std::vector<std::shared_ptr<const AscendingGrid>> freqs;
  freqs.reserve(nfreq);
  for (Index i = 0; i < nfreq; i++) {
    AscendingGrid f;
    xml_read_from_stream(is, f, pbifs);
    freqs.push_back(std::make_shared<AscendingGrid>(std::move(f)));
  }

  std::vector<std::shared_ptr<const SensorPosLosVector>> plos;
  plos.reserve(nposlos);
  for (Index i = 0; i < nposlos; i++) {
    SensorPosLosVector p;
    xml_read_from_stream(is, p, pbifs);
    plos.push_back(std::make_shared<SensorPosLosVector>(std::move(p)));
  }

  g.reserve(nelem);
  for (Index i = 0; i < nelem; i++) {
    Index ifreq;
    Index iplos;
    sensor::SparseStokvecMatrix weight_matrix;
    try {
      xml_read_from_stream(is, ifreq, pbifs);
      xml_read_from_stream(is, iplos, pbifs);
      xml_read_from_stream(is, weight_matrix, pbifs);

      ARTS_USER_ERROR_IF(ifreq >= nfreq, "Frequency index out of range")
      ARTS_USER_ERROR_IF(iplos >= nposlos, "Position index out of range")

      g.emplace_back(freqs.at(ifreq), plos.at(iplos), weight_matrix);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          std::format("Error reading ArrayOfSensorObsel element {}/{}:\n{}",
                      i,
                      nelem,
                      e.what()));
    }
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

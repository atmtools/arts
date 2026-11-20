#include "obsel.h"

#include <debug.h>

#include <algorithm>
#include <array>
#include <format>
#include <limits>
#include <numeric>
#include <ranges>
#include <utility>

namespace sensor {
bool SparseStokvec::operator==(const SparseStokvec& other) const {
  return (irow == other.irow) and (icol == other.icol);
}

bool SparseStokvec::operator!=(const SparseStokvec& other) const { return not(*this == other); }

bool SparseStokvec::operator<(const SparseStokvec& other) const {
  if (irow < other.irow) return true;
  if (irow > other.irow) return false;
  return icol < other.icol;
}

bool SparseStokvec::operator<=(const SparseStokvec& other) const { return (*this < other) or (*this == other); }

bool SparseStokvec::operator>(const SparseStokvec& other) const {
  if (irow > other.irow) return true;
  if (irow < other.irow) return false;
  return icol > other.icol;
}

bool SparseStokvec::operator>=(const SparseStokvec& other) const { return (*this > other) or (*this == other); }

bool SparseStokvecMatrix::operator==(const SparseStokvecMatrix& other) const {
  return rows == other.rows and cols == other.cols and sparse_data == other.sparse_data;
}

bool SparseStokvecMatrix::operator!=(const SparseStokvecMatrix& other) const { return not(*this == other); }

SparseStokvecMatrix::SparseStokvecMatrix(const StokvecMatrix& m) : rows(m.nrows()), cols(m.ncols()) {
  //! NOTE: Retains sorting
  for (Size i = 0; i < rows; ++i) {
    for (Size j = 0; j < cols; ++j) {
      if (m[i, j].is_zero()) continue;
      sparse_data.push_back({i, j, m[i, j]});
    }
  }
  assert(stdr::is_sorted(sparse_data));
}

SparseStokvecMatrix& SparseStokvecMatrix::operator=(const StokvecMatrix& m) { return *this = SparseStokvecMatrix(m); }

Index SparseStokvecMatrix::nrows() const { return rows; }

Index SparseStokvecMatrix::ncols() const { return cols; }

Size SparseStokvecMatrix::size() const { return sparse_data.size(); }

bool SparseStokvecMatrix::empty() const { return sparse_data.empty(); }

const std::vector<SparseStokvec>& SparseStokvecMatrix::vector() const { return sparse_data; }

void SparseStokvecMatrix::resize(Size nrows, Size ncols, Size size) {
  rows = nrows;
  cols = ncols;
  sparse_data.clear();
  sparse_data.resize(size);
}

std::array<Index, 2> SparseStokvecMatrix::shape() const { return {static_cast<Index>(rows), static_cast<Index>(cols)}; }

Stokvec& SparseStokvecMatrix::operator[](Size i, Size j) {
  const SparseStokvec newdata{.irow = i, .icol = j};
  auto                it = stdr::lower_bound(sparse_data, newdata);

  if (it != sparse_data.end() and newdata == *it) { return it->data; }

  auto newit = sparse_data.insert(it, newdata);

  //! NOTE: Must retain sorting?
  assert(stdr::is_sorted(sparse_data));

  return newit->data;
}

Stokvec SparseStokvecMatrix::operator[](Size i, Size j) const {
  const SparseStokvec newdata{.irow = i, .icol = j};
  auto                it = stdr::lower_bound(sparse_data, newdata);

  if (it != sparse_data.end() and newdata == *it) { return it->data; }

  return {};
}

std::vector<SparseStokvec>::iterator SparseStokvecMatrix::begin() { return sparse_data.begin(); }

std::vector<SparseStokvec>::iterator SparseStokvecMatrix::end() { return sparse_data.end(); }

std::vector<SparseStokvec>::const_iterator SparseStokvecMatrix::begin() const { return sparse_data.begin(); }

std::vector<SparseStokvec>::const_iterator SparseStokvecMatrix::end() const { return sparse_data.end(); }

SparseStokvecMatrix::operator StokvecMatrix() const {
  StokvecMatrix m(rows, cols, Stokvec{0.0, 0.0, 0.0, 0.0});
  for (const auto& elem : sparse_data) { m[elem.irow, elem.icol] = elem.data; }
  return m;
}

void normalize(SparseStokvecMatrix& w, Stokvec pol) {
  const Stokvec x = std::transform_reduce(
      w.begin(), w.end(), Stokvec{}, std::plus<>(), [](const SparseStokvec& e) { return e.data; });

  if (pol.is_polarized() and not x.is_polarized())
    arts::user_error(
        R"(The weights receive no polarization, so they cannot be normalized to the polarized {:B,}.

Normalization scales the weights and does not change their polarization.
Set the polarization of a channel when building the sensor instead.
)",
        pol);

  const Numeric polarized_x   = std::hypot(x.Q(), x.U(), x.V());
  const Numeric polarized_pol = std::hypot(pol.Q(), pol.U(), pol.V());

  Numeric scale = 0.0;
  if (x.I() != 0.0 and pol.I() != 0.0) {
    scale = std::abs(pol.I() / x.I());
  } else if (polarized_x != 0.0 and polarized_pol != 0.0) {
    scale = polarized_pol / polarized_x;
  } else {
    arts::user_error("Cannot normalize the summed weights {:B,} to {:B,}", x, pol);
  }

  for (auto& e : w) e.data *= scale;
}

namespace {
//! The weights of geometry ip, which are contiguous as the sparse matrix is sorted by row and then column
std::span<const SparseStokvec> geometry_weights(const SparseStokvecMatrix& w, Size ip) {
  assert(ip < static_cast<Size>(w.nrows()));

  const SparseStokvec v0{.irow = ip, .icol = static_cast<Size>(0)};
  const SparseStokvec vn{.irow = ip, .icol = std::numeric_limits<Size>::max()};
  return {stdr::lower_bound(w, v0), stdr::upper_bound(w, vn)};
}
}  // namespace

Numeric sumup(const SparseStokvecMatrix& w, const Stokvec& i, Size ip) {
  Numeric sum = 0.0;
  for (const auto& ws : geometry_weights(w, ip)) sum += dot(i, ws.data);
  return sum;
}

Numeric sumup(const SparseStokvecMatrix& w, const StokvecConstVectorView& i, Size ip) {
  assert(static_cast<Index>(i.size()) == w.ncols());

  return sumup(geometry_weights(w, ip), i);
}

void sumup(VectorView out, const SparseStokvecMatrix& w, const StokvecConstMatrixView& j, Size ip) {
  sumup(out, geometry_weights(w, ip), j);
}

Numeric sumup(std::span<const SparseStokvec> w, const StokvecConstVectorView& i) {
  Numeric sum = 0.0;
  for (const auto& ws : w) {
    assert(ws.icol < i.size());
    sum += dot(i[ws.icol], ws.data);
  }
  return sum;
}

void sumup(VectorView out, std::span<const SparseStokvec> w, const StokvecConstMatrixView& j) {
  // j is a matrix of shape JACS x frequencies
  for (Index ij = 0; ij < j.nrows(); ij++) {
    auto jac = j[ij];
    for (auto& ws : w) { out[ij] += dot(jac[ws.icol], ws.data); }
  }
}

SparseStokvecMatrix embed(
    const SparseStokvecMatrix& w, Size nrows, Size ncols, Size row_offset, std::span<const Size> col_map) {
  if (static_cast<Size>(w.nrows()) + row_offset > nrows or static_cast<Size>(w.ncols()) != col_map.size())
    arts::user_error("Cannot place weights of shape {:B,} at row {} of a {} by {} matrix with {} column indices",
                     w.shape(),
                     row_offset,
                     nrows,
                     ncols,
                     col_map.size());

  SparseStokvecMatrix out;
  out.resize(nrows, ncols, w.size());
  stdr::transform(w, out.begin(), [row_offset, col_map](SparseStokvec v) {
    v.irow += row_offset;
    v.icol  = col_map[v.icol];
    return v;
  });
  return out;
}

SparseStokvecMatrix embed(
    SparseStokvecMatrix&& w, Size nrows, Size ncols, Size row_offset, std::span<const Size> col_map) {
  if (static_cast<Size>(w.nrows()) + row_offset > nrows or static_cast<Size>(w.ncols()) != col_map.size())
    arts::user_error("Cannot place weights of shape {:B,} at row {} of a {} by {} matrix with {} column indices",
                     w.shape(),
                     row_offset,
                     nrows,
                     ncols,
                     col_map.size());

  SparseStokvecMatrix out(nrows, ncols);
  out.sparse_data = std::move(w.sparse_data);
  for (auto& v : out.sparse_data) {
    v.irow += row_offset;
    v.icol  = col_map[v.icol];
  }
  w = SparseStokvecMatrix{};
  return out;
}

}  // namespace sensor

void xml_io_stream<SensorPosLos>::write(std::ostream&       os,
                                        const SensorPosLos& x,
                                        bofstream*          pbofs,
                                        std::string_view    name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  if (pbofs) {
    put({&x, 1}, pbofs);
  } else {
    std::println(os, "{:IO}", x);
  }

  tag.write_to_end_stream(os);
}

void xml_io_stream<SensorPosLos>::read(std::istream& is, SensorPosLos& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get({&x, 1}, pbifs);
  } else {
    parse({&x, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<SensorPosLos>::put(std::span<const SensorPosLos> x, bofstream* pbofs) {
  xml_io_stream<Numeric>::put(std::span{reinterpret_cast<const Numeric*>(x.data()), x.size() * 5}, pbofs);
}

void xml_io_stream<SensorPosLos>::get(std::span<SensorPosLos> x, bifstream* pbifs) {
  xml_io_stream<Numeric>::get(std::span{reinterpret_cast<Numeric*>(x.data()), x.size() * 5}, pbifs);
}

void xml_io_stream<SensorPosLos>::parse(std::span<SensorPosLos> x, std::istream& is) {
  xml_io_stream<Numeric>::parse(std::span{reinterpret_cast<Numeric*>(x.data()), x.size() * 5}, is);
}

void xml_io_stream<sensor::SparseStokvec>::put(std::span<const sensor::SparseStokvec> x, bofstream* pbofs) {
  for (auto& elem : x) {
    xml_io_stream<Size>::put({&elem.irow, 1}, pbofs);
    xml_io_stream<Size>::put({&elem.icol, 1}, pbofs);
    xml_io_stream<Numeric>::put({elem.data.data_handle(), 4}, pbofs);
  }
}

void xml_io_stream<sensor::SparseStokvec>::write(std::ostream&                os,
                                                 const sensor::SparseStokvec& x,
                                                 bofstream*                   pbofs,
                                                 std::string_view             name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  if (pbofs) {
    put({&x, 1}, pbofs);
  } else {
    std::println(os, "{:IO}", x);
  }

  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::SparseStokvec>::get(std::span<sensor::SparseStokvec> x, bifstream* pbifs) {
  for (auto& elem : x) {
    xml_io_stream<Size>::get({&elem.irow, 1}, pbifs);
    xml_io_stream<Size>::get({&elem.icol, 1}, pbifs);
    xml_io_stream<Numeric>::get({elem.data.data_handle(), 4}, pbifs);
  }
}

void xml_io_stream<sensor::SparseStokvec>::parse(std::span<sensor::SparseStokvec> x, std::istream& is) {
  for (auto& v : x) { is >> v.irow >> v.icol >> v.data; }
}

void xml_io_stream<sensor::SparseStokvec>::read(std::istream& is, sensor::SparseStokvec& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  if (pbifs) {
    get({&x, 1}, pbifs);
  } else {
    parse({&x, 1}, is);
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<sensor::SparseStokvecMatrix>::write(std::ostream&                      os,
                                                       const sensor::SparseStokvecMatrix& x,
                                                       bofstream*                         pbofs,
                                                       std::string_view                   name) {
  XMLTag tag(type_name, "name", name);
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.nrows(), pbofs);
  xml_write_to_stream(os, x.ncols(), pbofs);
  xml_write_to_stream(os, x.vector(), pbofs);

  tag.write_to_end_stream(os);
}

void xml_io_stream<sensor::SparseStokvecMatrix>::read(std::istream&                is,
                                                      sensor::SparseStokvecMatrix& x,
                                                      bifstream*                   pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  Index                              nr, nc;
  std::vector<sensor::SparseStokvec> vector;

  xml_read_from_stream(is, nr, pbifs);
  xml_read_from_stream(is, nc, pbifs);
  xml_read_from_stream(is, vector, pbifs);

  x.resize(nr, nc, vector.size());

  stdr::move(stdv::as_rvalue(vector), x.begin());

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

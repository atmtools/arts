#pragma once

#include <matpack.h>

#include <variant>

namespace Atm::interp {
using namespace lagrange_interp;

using altlag1 = lag_t<1>;
using altlag0 = lag_t<0>;

using latlag1 = lag_t<1>;
using latlag0 = lag_t<0>;

using lonlag1 = lag_t<1, loncross>;
using lonlag0 = lag_t<0, loncross>;

using altlags = std::variant<altlag0, altlag1>;
using latlags = std::variant<latlag0, latlag1>;
using lonlags = std::variant<lonlag0, lonlag1>;

altlags altlag(const AscendingGrid& xs, Numeric x);
latlags latlag(const LatGrid& xs, Numeric x);
lonlags lonlag(const LonGrid& xs, Numeric x);

Tensor3 interpweights(const altlags&, const latlags&, const lonlags&);
Matrix interpweights(const latlags&, const lonlags&);

template <lagrange_interp::value_transformer transform>
Numeric get(const GeodeticField3& f,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [&data = f.data](auto&& al, auto&& la, auto&& lo) {
        return lagrange_interp::interp<transform>(data, al, la, lo);
      },
      altlag(f.grid<0>(), alt),
      latlag(f.grid<1>(), lat),
      lonlag(f.grid<2>(), lon));
}

template <lagrange_interp::value_transformer transform>
Numeric get(const GeodeticField3& f,
            const Tensor3& iw,
            const altlags& alt,
            const latlags& lat,
            const lonlags& lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [&data = f.data, &iw](auto&& al, auto&& la, auto&& lo) {
        return lagrange_interp::interp<transform>(data, iw, al, la, lo);
      },
      alt,
      lat,
      lon);
}

template <lagrange_interp::value_transformer transform>
Numeric get(const GeodeticField3& f,
            Index ialt,
            const Matrix& iw,
            const latlags& lat,
            const lonlags& lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [data = f.data[ialt], &iw](auto&& la, auto&& lo) {
        return lagrange_interp::interp<transform>(data, iw, la, lo);
      },
      lat,
      lon);
}

// Numeric interpolation
template <lagrange_interp::value_transformer transform>
Numeric get(const Numeric num, const Numeric, const Numeric, const Numeric) {
  return num;
}

// Functional interpolation
template <lagrange_interp::value_transformer transform>
Numeric get(const std::function<Numeric(Numeric, Numeric, Numeric)>& fd,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  return fd(alt, lat, lon);
}

std::vector<std::pair<Index, Numeric>> flat_weight(const GeodeticField3& gf3,
                                                   const Numeric alt,
                                                   const Numeric lat,
                                                   const Numeric lon);
}  // namespace Atm::interp

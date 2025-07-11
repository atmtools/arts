#include "atm_field_interp.h"

namespace Atm::interp {
altlags altlag(const AscendingGrid& xs, Numeric x) {
  return xs.size() == 1 ? altlags{altlag0(0, x, xs)}
                        : altlags{altlag1(0, x, xs)};
}

latlags latlag(const AscendingGrid& xs, Numeric x) {
  return xs.size() == 1 ? latlags{latlag0(0, x, xs)}
                        : latlags{latlag1(0, x, xs)};
}

lonlags lonlag(const AscendingGrid& xs, Numeric x) {
  return xs.size() == 1 ? lonlags{lonlag0(0, x, xs)}
                        : lonlags{lonlag1(0, x, xs)};
}

Tensor3 interpweights(const altlags& a, const latlags& b, const lonlags& c) {
  return std::visit(
      [](auto& alt, auto& lat, auto& lon) -> Tensor3 {
        auto d = my_interp::interpweights(alt, lat, lon);
        Tensor3 out(alt.size(), lat.size(), lon.size());
        for (Size i0 = 0; i0 < alt.size(); i0++)
          for (Size i1 = 0; i1 < lat.size(); i1++)
            for (Size i2 = 0; i2 < lon.size(); i2++)
              out[i0, i1, i2] = d[i0, i1, i2];
        return out;
      },
      a,
      b,
      c);
}

Matrix interpweights(const latlags& b, const lonlags& c) {
  return std::visit(
      [](auto& lat, auto& lon) -> Matrix {
        auto d = my_interp::interpweights(lat, lon);
        Matrix out(lat.size(), lon.size());
        for (Size i1 = 0; i1 < lat.size(); i1++)
          for (Size i2 = 0; i2 < lon.size(); i2++) out[i1, i2] = d[i1, i2];
        return out;
      },
      b,
      c);
}

Numeric get(const SortedGriddedField3& f,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [&data = f.data](auto&& al, auto&& la, auto&& lo) {
        return my_interp::interp(data, al, la, lo);
      },
      altlag(f.grid<0>(), alt),
      latlag(f.grid<1>(), lat),
      lonlag(f.grid<2>(), lon));
}

Numeric get(const SortedGriddedField3& f,
            const Tensor3& iw,
            const altlags& alt,
            const latlags& lat,
            const lonlags& lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [&data = f.data, &iw](auto&& al, auto&& la, auto&& lo) {
        return my_interp::interp(data, iw, al, la, lo);
      },
      alt,
      lat,
      lon);
}

Numeric get(const SortedGriddedField3& f,
            Index ialt,
            const Matrix& iw,
            const latlags& lat,
            const lonlags& lon) {
  if (not f.ok()) throw std::runtime_error("bad field");
  return std::visit(
      [data = f.data[ialt], &iw](auto&& la, auto&& lo) {
        return my_interp::interp(data, iw, la, lo);
      },
      lat,
      lon);
}


Numeric get(const Numeric num,
                      const Numeric,
                      const Numeric,
                      const Numeric) {
  return num;
}

Numeric get(const FunctionalData &fd,
            const Numeric alt,
            const Numeric lat,
            const Numeric lon) {
  return fd(alt, lat, lon);
}
}  // namespace Atm::interp

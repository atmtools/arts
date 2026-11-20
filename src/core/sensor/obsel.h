#pragma once

#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

#include <span>

#include "matpack_mdspan_helpers_grid_t.h"

namespace sensor {
struct PosLos {
  Vector3 pos;
  Vector2 los;

  constexpr bool operator==(const PosLos& other) const = default;
  constexpr bool operator!=(const PosLos& other) const = default;

  friend std::ostream& operator<<(std::ostream& os, const PosLos& poslos);

  template <typename Self> [[nodiscard]] constexpr decltype(auto) alt(this Self&& self) {
    return std::forward<Self>(self).pos[0];
  }

  template <typename Self> [[nodiscard]] constexpr decltype(auto) lat(this Self&& self) {
    return std::forward<Self>(self).pos[1];
  }

  template <typename Self> [[nodiscard]] constexpr decltype(auto) lon(this Self&& self) {
    return std::forward<Self>(self).pos[2];
  }

  template <typename Self> [[nodiscard]] constexpr decltype(auto) zen(this Self&& self) {
    return std::forward<Self>(self).los[0];
  }

  template <typename Self> [[nodiscard]] constexpr decltype(auto) azi(this Self&& self) {
    return std::forward<Self>(self).los[1];
  }
};

using PosLosVector = matpack::data_t<PosLos, 1>;

struct SparseStokvec {
  Size    irow;
  Size    icol;
  Stokvec data{0, 0, 0, 0};

  bool operator==(const SparseStokvec& other) const;
  bool operator!=(const SparseStokvec& other) const;
  bool operator<(const SparseStokvec& other) const;
  bool operator<=(const SparseStokvec& other) const;
  bool operator>(const SparseStokvec& other) const;
  bool operator>=(const SparseStokvec& other) const;
};

class SparseStokvecMatrix {
  Size rows;
  Size cols;

  std::vector<SparseStokvec> sparse_data{};

 public:
  SparseStokvecMatrix(Size r = 0, Size c = 0) : rows(r), cols(c) {}
  SparseStokvecMatrix(const SparseStokvecMatrix&)            = default;
  SparseStokvecMatrix(SparseStokvecMatrix&&)                 = default;
  SparseStokvecMatrix& operator=(const SparseStokvecMatrix&) = default;
  SparseStokvecMatrix& operator=(SparseStokvecMatrix&&)      = default;

  SparseStokvecMatrix(const StokvecMatrix& m);
  SparseStokvecMatrix& operator=(const StokvecMatrix& m);

  [[nodiscard]] bool operator==(const SparseStokvecMatrix& other) const;
  [[nodiscard]] bool operator!=(const SparseStokvecMatrix& other) const;

  [[nodiscard]] Index                nrows() const;
  [[nodiscard]] Index                ncols() const;
  [[nodiscard]] Size                 size() const;
  [[nodiscard]] bool                 empty() const;
  [[nodiscard]] std::array<Index, 2> shape() const;
  const std::vector<SparseStokvec>&  vector() const;

  void resize(Size rows, Size cols, Size size);

  Stokvec& operator[](Size i, Size j);
  Stokvec  operator[](Size i, Size j) const;

  [[nodiscard]] std::vector<SparseStokvec>::iterator       begin();
  [[nodiscard]] std::vector<SparseStokvec>::iterator       end();
  [[nodiscard]] std::vector<SparseStokvec>::const_iterator begin() const;
  [[nodiscard]] std::vector<SparseStokvec>::const_iterator end() const;

  explicit operator StokvecMatrix() const;

  friend SparseStokvecMatrix embed(
      SparseStokvecMatrix&& w, Size nrows, Size ncols, Size row_offset, std::span<const Size> col_map);
};

/** The weights are scaled by one common factor.
 *
 * After calling this method, the I of the summed weights is the I of pol.
 * When either has no I, the magnitudes of their summed Q, U, and V agree
 * instead.  The polarization of the weights is not changed, so it is an
 * error to ask for a polarized pol when the weights have no polarization.
 *
 * @param w The weights to scale
 * @param pol The summed weights to scale to. The default is [1, 0, 0, 0].
 */
void normalize(SparseStokvecMatrix& w, Stokvec pol = {1., 0., 0., 0.});

//! The weights of geometry ip dotted with i at every frequency, summed
[[nodiscard]] Numeric sumup(const SparseStokvecMatrix& w, const Stokvec& i, Size ip);

//! The weights of geometry ip dotted with the spectral radiance i, summed
[[nodiscard]] Numeric sumup(const SparseStokvecMatrix& w, const StokvecConstVectorView& i, Size ip);

//! Adds the weights of geometry ip dotted with each row of the spectral radiance Jacobian j to out
void sumup(VectorView out, const SparseStokvecMatrix& w, const StokvecConstMatrixView& j, Size ip);

//! The weights w, all of one geometry, dotted with the spectral radiance i, summed
[[nodiscard]] Numeric sumup(std::span<const SparseStokvec> w, const StokvecConstVectorView& i);

//! Adds the weights w, all of one geometry, dotted with each row of the spectral radiance Jacobian j to out
void sumup(VectorView out, std::span<const SparseStokvec> w, const StokvecConstMatrixView& j);

/** The weights placed in a larger matrix of nrows by ncols
 *
 *  Row i moves to row i + row_offset, and column j to column col_map[j], which
 *  must increase with j so the weights stay sorted.
 */
[[nodiscard]] SparseStokvecMatrix embed(
    const SparseStokvecMatrix& w, Size nrows, Size ncols, Size row_offset, std::span<const Size> col_map);

//! As embed above, but reusing the memory of w
[[nodiscard]] SparseStokvecMatrix embed(
    SparseStokvecMatrix&& w, Size nrows, Size ncols, Size row_offset, std::span<const Size> col_map);
}  // namespace sensor

using SensorPosLos       = sensor::PosLos;
using SensorPosLosVector = sensor::PosLosVector;

template <> struct std::formatter<SensorPosLos> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SensorPosLos& v, FmtContext& ctx) const {
    return tags.format(ctx, v.pos, tags.sep(), v.los);
  }
};

template <> struct std::formatter<sensor::SparseStokvec> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const sensor::SparseStokvec& w, FmtContext& ctx) const {
    const auto sep = tags.sep();

    tags.add_if_bracket(ctx, "["sv);
    tags.format(ctx, w.irow, sep, w.icol, sep, w.data);
    tags.add_if_bracket(ctx, "]"sv);

    return ctx.out();
  }
};

template <> struct std::formatter<sensor::SparseStokvecMatrix> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext> FmtContext::iterator format(const sensor::SparseStokvecMatrix& v, FmtContext& ctx) const {
    const auto sep = tags.sep();

    tags.add_if_bracket(ctx, "["sv);
    for (auto& w : v) tags.format(ctx, w, sep);
    tags.add_if_bracket(ctx, "]"sv);
    return ctx.out();
  }
};

template <> struct xml_io_stream<SensorPosLos> {
  static constexpr std::string_view type_name = "SensorPosLos"sv;

  static void write(std::ostream& os, const SensorPosLos& x, bofstream* pbofs = nullptr, std::string_view name = ""sv);

  static void read(std::istream& is, SensorPosLos& x, bifstream* pbifs = nullptr);

  static void put(std::span<const SensorPosLos> x, bofstream*);
  static void get(std::span<SensorPosLos> x, bifstream*);
  static void parse(std::span<SensorPosLos> x, std::istream&);
};

template <> struct xml_io_stream<sensor::SparseStokvec> {
  static constexpr std::string_view type_name = "SparseStokvec"sv;

  static void write(std::ostream&                os,
                    const sensor::SparseStokvec& x,
                    bofstream*                   pbofs = nullptr,
                    std::string_view             name  = ""sv);

  static void read(std::istream& is, sensor::SparseStokvec& x, bifstream* pbifs = nullptr);

  static void put(std::span<const sensor::SparseStokvec> x, bofstream*);
  static void get(std::span<sensor::SparseStokvec> x, bifstream*);
  static void parse(std::span<sensor::SparseStokvec> x, std::istream&);
};

template <> struct xml_io_stream<sensor::SparseStokvecMatrix> {
  static constexpr std::string_view type_name = "SparseStokvecMatrix"sv;

  static void write(std::ostream&                      os,
                    const sensor::SparseStokvecMatrix& x,
                    bofstream*                         pbofs = nullptr,
                    std::string_view                   name  = ""sv);

  static void read(std::istream& is, sensor::SparseStokvecMatrix& x, bifstream* pbifs = nullptr);
};

#pragma once

#include <enumsSensorJacobianModelType.h>
#include <enumsSensorKeyType.h>
#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

#include <boost/container_hash/hash.hpp>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace sensor {
struct PosLos {
  Vector3 pos;
  Vector2 los;

  constexpr bool operator==(const PosLos& other) const = default;
  constexpr bool operator!=(const PosLos& other) const = default;

  friend std::ostream& operator<<(std::ostream& os, const PosLos& poslos);

  [[nodiscard]] constexpr Numeric alt() const { return pos[0]; }
  [[nodiscard]] constexpr Numeric lat() const { return pos[1]; }
  [[nodiscard]] constexpr Numeric lon() const { return pos[2]; }
  [[nodiscard]] constexpr Numeric za() const { return los[0]; }
  [[nodiscard]] constexpr Numeric aa() const { return los[1]; }
};

using PosLosVector = matpack::data_t<PosLos, 1>;

struct SparseStokvec {
  Size irow;
  Size icol;
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

  [[nodiscard]] Index nrows() const;
  [[nodiscard]] Index ncols() const;
  [[nodiscard]] Size size() const;
  [[nodiscard]] bool empty() const;
  [[nodiscard]] std::array<Index, 2> shape() const;
  const std::vector<SparseStokvec>& vector() const;

  void resize(Size rows, Size cols, Size size);

  Stokvec& operator[](Size i, Size j);
  Stokvec operator[](Size i, Size j) const;

  [[nodiscard]] std::vector<SparseStokvec>::iterator begin();
  [[nodiscard]] std::vector<SparseStokvec>::iterator end();
  [[nodiscard]] std::vector<SparseStokvec>::const_iterator begin() const;
  [[nodiscard]] std::vector<SparseStokvec>::const_iterator end() const;

  explicit operator StokvecMatrix() const;
};

class Obsel {
  //! Frequency grid, must be ascending
  std::shared_ptr<const AscendingGrid> f{
      std::shared_ptr<const AscendingGrid>(new AscendingGrid{})};

  //! Position and line of sight grid of the sensor
  std::shared_ptr<const PosLosVector> poslos{
      std::shared_ptr<const PosLosVector>(new PosLosVector{})};

  //! FIXME: This should be made a variant of sparse/non-sparse!  Do this if/when we see an actual slowdown cf ARTS2.
  // (The type should be "std::variant<StokvecMatrix, std::array<Sparse, 4>>", where the "4" is for the 4 Stokes components.)
  // A matrix size of poslos_grid.size() x f_grid.size() with the polarized weight of the sensor
  SparseStokvecMatrix w{};

 public:
  Obsel();
  Obsel(const Obsel&);
  Obsel(Obsel&&) noexcept;
  Obsel& operator=(const Obsel&);
  Obsel& operator=(Obsel&&) noexcept;

  Obsel(std::shared_ptr<const AscendingGrid> fs,
        std::shared_ptr<const PosLosVector> pl,
        SparseStokvecMatrix ws);
  Obsel(const AscendingGrid& fs,
        const PosLosVector& pl,
        SparseStokvecMatrix ws);

  void check() const;

  [[nodiscard]] bool same_freqs(const Obsel& other) const;

  [[nodiscard]] bool same_freqs(
      const std::shared_ptr<const AscendingGrid>& other) const;

  [[nodiscard]] bool same_poslos(const Obsel& other) const;

  [[nodiscard]] bool same_poslos(
      const std::shared_ptr<const PosLosVector>& other) const;

  [[nodiscard]] const auto& f_grid_ptr() const { return f; }
  [[nodiscard]] const auto& poslos_grid_ptr() const { return poslos; }

  [[nodiscard]] const AscendingGrid& f_grid() const { return *f; }
  [[nodiscard]] const PosLosVector& poslos_grid() const { return *poslos; }
  [[nodiscard]] const SparseStokvecMatrix& weight_matrix() const { return w; }

  void set_f_grid_ptr(std::shared_ptr<const AscendingGrid> n);

  void set_poslos_grid_ptr(std::shared_ptr<const PosLosVector> n);

  void set_weight_matrix(SparseStokvecMatrix n);

  //! Constant indicating that the frequency or poslos is not found in the grid
  constexpr static Index dont_have = -1;

  friend std::ostream& operator<<(std::ostream& os, const Obsel& obsel);

  /** The weights are renormalized to the new value.
   * 
   * After calling this method, the sum of all weights dot-producted with the
   * polarization will be equal to the new value.
   *
   * @param pol The polarization that is sampled for the sum. The default is [1, 0, 0, 0].
   */
  void normalize(Stokvec pol = {1., 0., 0., 0.});

  [[nodiscard]] Numeric sumup(const StokvecVectorView& i, Index ip) const;
  void sumup(VectorView out, const StokvecMatrixView& j, Index ip) const;

  [[nodiscard]] Size flat_size(const SensorKeyType& key) const;
  void flat(VectorView x, const SensorKeyType& key) const;
  [[nodiscard]] Vector flat(const SensorKeyType& key) const;

  [[nodiscard]] Index find(const Vector3&, const Vector2&) const;
  [[nodiscard]] Index find(const AscendingGrid& frequency_grid) const;
};

std::ostream& operator<<(std::ostream& os, const Array<Obsel>& obsel);
}  // namespace sensor

struct SensorKey {
  SensorKeyType type{};

  Index sensor_elem{};

  Index measurement_elem{};

  bool operator==(const SensorKey& other) const;
};

template <>
struct std::hash<SensorKey> {
  std::size_t operator()(const SensorKey& g) const {
    std::size_t seed = 0;
    boost::hash_combine(seed, std::hash<SensorKeyType>{}(g.type));
    boost::hash_combine(seed, std::hash<Index>{}(g.sensor_elem));
    boost::hash_combine(seed, std::hash<Index>{}(g.measurement_elem));
    return seed;
  }
};

template <>
struct std::formatter<SensorKey> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SensorKey& v, FmtContext& ctx) const {
    return tags.format(
        ctx, v.type, tags.sep(), v.sensor_elem, tags.sep(), v.measurement_elem);
  }
};

using SensorPosLos       = sensor::PosLos;
using SensorPosLosVector = sensor::PosLosVector;
using SensorObsel        = sensor::Obsel;
using ArrayOfSensorObsel = Array<SensorObsel>;
using SensorSimulations  = std::unordered_map<
     std::shared_ptr<const AscendingGrid>,
     std::unordered_set<std::shared_ptr<const SensorPosLosVector>>>;

void unflatten(ArrayOfSensorObsel& sensor,
               const ConstVectorView& x,
               const SensorObsel& v,
               const SensorKeyType& key);

/** Make the sensor obsels exhaustive.
 * 
 * Exhaustive means that all poslos and frequency points are present
 * in every obsel.  This generally makes sensors that have a small
 * continous grid of frequencies and mostly overlapping observation
 * geometries.
 * 
 * @param obsels An existing list of observation elements to make exhaustive
 */
void make_exhaustive(std::span<SensorObsel> obsels);

/** Make the sensor obsels exclusive.
 * 
 * Exclusive means that all observation elements are completely 
 * independent of each other.  This can improve simulation speed
 * when there are sparse number of grids or when the elements have
 * very different observation geometries.
 *  
 * @param obsels An existing list of observation elements to make exclusive
 */
void make_exclusive(std::span<SensorObsel> obsels);

SensorSimulations collect_simulations(
    const std::span<const SensorObsel>& obsels);

template <>
struct std::formatter<SensorPosLos> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SensorPosLos& v, FmtContext& ctx) const {
    return tags.format(ctx, v.pos, tags.sep(), v.los);
  }
};

template <>
struct std::formatter<sensor::SparseStokvec> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::SparseStokvec& w,
                              FmtContext& ctx) const {
    const auto sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    tags.format(ctx, w.irow, sep, w.icol, sep, w.data);
    tags.add_if_bracket(ctx, ']');

    return ctx.out();
  }
};

template <>
struct std::formatter<sensor::SparseStokvecMatrix> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    std::format_parse_context::iterator v = parse_format_tags(tags, ctx);
    tags.newline                          = not tags.newline;
    return v;
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::SparseStokvecMatrix& v,
                              FmtContext& ctx) const {
    const auto sep = tags.sep();

    tags.add_if_bracket(ctx, '[');
    for (auto& w : v) tags.format(ctx, w, sep);
    tags.add_if_bracket(ctx, ']');
    return ctx.out();
  }
};

template <>
struct std::formatter<SensorObsel> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const SensorObsel& v, FmtContext& ctx) const {
    return tags.format(ctx,
                       "Obsel:\n  frequency grid:           "sv,
                       v.f_grid(),
                       "\n  pos-los grid:                   "sv,
                       v.poslos_grid(),
                       "\n  weights:                        "sv,
                       v.weight_matrix());
  }
};

template <>
struct xml_io_stream<SensorPosLos> {
  static constexpr std::string_view type_name = "SensorPosLos"sv;

  static void write(std::ostream& os,
                    const SensorPosLos& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SensorPosLos& x,
                   bifstream* pbifs = nullptr);

  static void put(std::span<const SensorPosLos> x, bofstream*);
  static void get(std::span<SensorPosLos> x, bifstream*);
  static void parse(std::span<SensorPosLos> x, std::istream&);
};

template <>
struct xml_io_stream_name<SensorKey> {
  static constexpr std::string_view name = "SensorKey"sv;
};

template <>
struct xml_io_stream_aggregate<SensorKey> {
  static constexpr bool value = true;
};

template <>
struct xml_io_stream<sensor::SparseStokvec> {
  static constexpr std::string_view type_name = "SparseStokvec"sv;

  static void write(std::ostream& os,
                    const sensor::SparseStokvec& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   sensor::SparseStokvec& x,
                   bifstream* pbifs = nullptr);

  static void put(std::span<const sensor::SparseStokvec> x, bofstream*);
  static void get(std::span<sensor::SparseStokvec> x, bifstream*);
  static void parse(std::span<sensor::SparseStokvec> x, std::istream&);
};

template <>
struct xml_io_stream<sensor::SparseStokvecMatrix> {
  static constexpr std::string_view type_name = "SparseStokvecMatrix"sv;

  static void write(std::ostream& os,
                    const sensor::SparseStokvecMatrix& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   sensor::SparseStokvecMatrix& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<SensorObsel> {
  static constexpr std::string_view type_name = "SensorObsel"sv;

  static void write(std::ostream& os,
                    const SensorObsel& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   SensorObsel& x,
                   bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<ArrayOfSensorObsel> {
  static constexpr std::string_view type_name = "ArrayOfSensorObsel"sv;

  static void write(std::ostream& os,
                    const ArrayOfSensorObsel& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   ArrayOfSensorObsel& x,
                   bifstream* pbifs = nullptr);
};

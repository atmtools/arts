#pragma once

#include <enumsSensorJacobianModelType.h>
#include <enumsSensorKeyType.h>
#include <matpack.h>
#include <rtepack.h>

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
  
  void resize(Size rows, Size cols, Size reserve);

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
      std::make_shared<const AscendingGrid>()};

  //! Position and line of sight grid of the sensor
  std::shared_ptr<const PosLosVector> poslos{
      std::make_shared<const PosLosVector>()};

  //! FIXME: This should be made a variant of sparse/non-sparse!  Do this if/when we see an actual slowdown cf ARTS2.
  // (The type should be "std::variant<StokvecMatrix, std::array<Sparse, 4>>", where the "4" is for the 4 Stokes components.)
  // A matrix size of poslos_grid.size() x f_grid.size() with the polarized weight of the sensor
  SparseStokvecMatrix w{};

 public:
  Obsel() = default;

  Obsel(std::shared_ptr<const AscendingGrid> fs,
        std::shared_ptr<const PosLosVector> pl,
        SparseStokvecMatrix ws);
  Obsel(const AscendingGrid& fs,
        const PosLosVector& pl,
        SparseStokvecMatrix ws);

  Obsel(const Obsel&)            = default;
  Obsel(Obsel&&)                 = default;
  Obsel& operator=(const Obsel&) = default;
  Obsel& operator=(Obsel&&)      = default;

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

  SensorJacobianModelType model{};

  Index polyorder{-1};

  Vector original_grid{};

  bool operator==(const SensorKey& other) const;
};

template <>
struct std::hash<SensorKey> {
  std::size_t operator()(const SensorKey& g) const {
    return std::hash<SensorKeyType>{}(g.type) ^
           (std::hash<SensorJacobianModelType>{}(g.model)
            << (8 * sizeof(SensorKeyType))) ^
           (std::hash<Index>{}(g.measurement_elem)
            << (8 * (sizeof(SensorKeyType) + sizeof(SensorJacobianModelType))));
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
    return tags.format(ctx,
                       v.type,
                       tags.sep(),
                       v.sensor_elem,
                       tags.sep(),
                       v.measurement_elem,
                       tags.sep(),
                       v.model);
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
struct std::formatter<sensor::SparseStokvecMatrix> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(
      std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext>
  FmtContext::iterator format(const sensor::SparseStokvecMatrix& v,
                              FmtContext& ctx) const {
    for (auto& w : v) {
      tags.format(ctx, w.irow, ' ', w.icol, ' ', w.data);
    }
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

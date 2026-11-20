#pragma once

#include <enumsSensorKeyType.h>
#include <matpack.h>
#include <rtepack.h>
#include <xml.h>

#include <boost/container_hash/hash.hpp>
#include <span>
#include <utility>
#include <vector>

#include "obsel.h"
#include "sensor_meta_info.h"

namespace sensor {
/** A sensor: observation elements on one frequency grid and one list of geometries
 *
 *  The radiance is simulated at every geometry over the whole frequency grid.
 *  Each element is one value of the measurement: its Stokes weights of each
 *  geometry (rows) and frequency (columns), which describe the channel response
 *  and the antenna response together.  The meta describes the elements as a
 *  structured gridded field, one value per element in their order, and its
 *  data name is the name of the sensor.
 */
struct Sensor {
  //! The frequencies at which the radiance is simulated
  AscendingGrid freq_grid{};

  //! The positions and lines of sight at which the radiance is simulated
  PosLosVector geometry{};

  //! The observation elements, each a geometry times frequency sparse matrix of Stokvec weights
  std::vector<SparseStokvecMatrix> elements{};

  //! The structure of the elements, named as the sensor
  MetaInfo meta{};

  //! The name of the sensor, the data name of its meta: a label for messages, not necessarily unique
  [[nodiscard]] const String& name() const;
  [[nodiscard]] String&       name();

  //! Throws if the weights or the meta do not match the grids and the number of elements
  void check() const;

  //! Normalizes each element, see sensor::normalize
  void normalize(const Stokvec& pol = {1., 0., 0., 0.});

  /** One value per element: its weights times the spectral radiance, summed
   *
   *  @param spectral_rad The spectral radiance of each geometry (rows) at each frequency (columns)
   */
  [[nodiscard]] Vector extract(const StokvecMatrix& spectral_rad) const;

  [[nodiscard]] Size   flat_size(const SensorKeyType& key) const;
  void                 flat(VectorView x, const SensorKeyType& key) const;
  [[nodiscard]] Vector flat(const SensorKeyType& key) const;

  //! Sets the frequency grid or a coordinate of the geometries from x, the inverse of flat
  void unflatten(const ConstVectorView& x, const SensorKeyType& key);
};

//! The weights of one element at one geometry
struct GeometryWeights {
  //! The index of the element in the sensor
  Size element;

  //! The weights of the element at the geometry
  std::span<const SparseStokvec> weights;
};

/** The elements of a sensor that have weights at each geometry
 *
 *  An element keeps its weights sorted by geometry and then by frequency, so
 *  its weights at one geometry are contiguous.  The index lists these by
 *  geometry instead: for each geometry, the elements that have weights there,
 *  in increasing order, each with its weights there.  Only the pairs of element
 *  and geometry that have weights are listed.
 *
 *  The weights are views of those of the sensor, so the index is only valid
 *  while the sensor and its elements are neither changed nor destroyed.
 */
class GeometryIndex {
  //! The entries of geometry ip are entries[starts[ip]] to entries[starts[ip + 1]]
  std::vector<Size> starts{};

  std::vector<GeometryWeights> entries{};

  //! The number of geometries at which each element has weights
  std::vector<Size> counts{};

 public:
  //! Throws if a weight is outside the geometries and frequencies of the sensor
  explicit GeometryIndex(const Sensor& sensor);

  //! The elements with weights at geometry ip, in increasing order
  [[nodiscard]] std::span<const GeometryWeights> operator[](Size ip) const;

  //! The number of geometries at which element ie has weights
  [[nodiscard]] Size geometry_count(Size ie) const;
};
}  // namespace sensor

struct SensorKey {
  SensorKeyType type{};

  //! The index of the sensor in the measurement sensor
  Index sensor_elem{};

  bool operator==(const SensorKey& other) const;
};

template <> struct std::hash<SensorKey> {
  static std::size_t operator()(const SensorKey& g) {
    std::size_t seed = 0;
    boost::hash_combine(seed, std::hash<SensorKeyType>{}(g.type));
    boost::hash_combine(seed, std::hash<Index>{}(g.sensor_elem));
    return seed;
  }
};

template <> struct std::formatter<SensorKey> {
  format_tags tags{};

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SensorKey& v, FmtContext& ctx) const {
    return tags.format(ctx, v.type, tags.sep(), v.sensor_elem);
  }
};

using Sensor        = sensor::Sensor;
using ArrayOfSensor = Array<Sensor>;

//! The number of measurement elements of all the sensors
[[nodiscard]] Size measurement_size(std::span<const Sensor> sensors);

//! The index of the first measurement element of sensor isensor
[[nodiscard]] Size measurement_start(std::span<const Sensor> sensors, Size isensor);

//! The sensor of measurement element iy, and the index of the element in it
[[nodiscard]] std::pair<Size, Size> measurement_element(std::span<const Sensor> sensors, Size iy);

/** What simulating the sensors is expected to cost
 *
 *  A table with, per sensor, its elements, the shape of its meta, its
 *  geometries (one agenda call each), its frequencies, the spectral points
 *  those calls simulate, and its non-zero weights; then the totals, and the
 *  geometries that more than one sensor simulates.
 */
[[nodiscard]] String cost_report(std::span<const Sensor> sensors);

template <> struct format_tag_aggregate<sensor::Sensor> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_name<SensorKey> {
  static constexpr std::string_view name = "SensorKey"sv;
};

template <> struct xml_io_stream_aggregate<SensorKey> {
  static constexpr bool value = true;
};

template <> struct xml_io_stream_name<sensor::Sensor> {
  static constexpr std::string_view name = "Sensor"sv;
};

template <> struct xml_io_stream_aggregate<sensor::Sensor> {
  static constexpr bool value = true;
};

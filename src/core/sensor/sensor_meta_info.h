#pragma once

#include <configtypes.h>
#include <format_tags.h>
#include <matpack.h>
#include <xml.h>

#include <variant>

namespace sensor {

/// A camera image gridded field:
///   grid<0> = AscendingGrid  (row angular offsets from boresight, degrees)
///   grid<1> = AscendingGrid  (column angular offsets from boresight, degrees)
///   grid<2> = Vector         (color / frequency values)
using CameraGriddedField = matpack::gridded_data_t<Numeric, AscendingGrid, AscendingGrid, Vector>;

/// Metadata describing the mapping from a contiguous block of
/// measurement_vec / measurement_sensor to a structured gridded field.
///
/// Each variant alternative corresponds to a different sensor type.
/// The gridded field's data_name serves as a human-readable label (e.g.,
/// "camera", "ffts").  Grids are set at sensor construction time; data is
/// zeroed and filled later by measurement_sensor_metaFromMeasurementVec.
///
/// The start index into measurement_vec is not stored — it is derived from
/// the cumulative count() of preceding elements in an ArrayOfSensorMetaInfo.
struct MetaInfo {
  using variant_t = std::variant<SortedGriddedField1, CameraGriddedField>;

  variant_t data{SortedGriddedField1{}};

  /// Total number of scalar elements (product of all grid sizes).
  [[nodiscard]] Size count() const;

  /// Validate that data dimensions match grid sizes.
  void check() const;
};

}  // namespace sensor

using SensorMetaInfo        = sensor::MetaInfo;
using ArrayOfSensorMetaInfo = Array<SensorMetaInfo>;

//
// Formatter
//

template <> struct std::formatter<SensorMetaInfo> {
  format_tags tags;

  [[nodiscard]] constexpr auto& inner_fmt() { return *this; }
  [[nodiscard]] constexpr auto& inner_fmt() const { return *this; }

  constexpr std::format_parse_context::iterator parse(std::format_parse_context& ctx) {
    return parse_format_tags(tags, ctx);
  }

  template <class FmtContext> FmtContext::iterator format(const SensorMetaInfo& x, FmtContext& ctx) const {
    return std::visit([this, &ctx](const auto& gf) { return tags.format(ctx, gf); }, x.data);
  }
};

//
// XML I/O
//

template <> struct xml_io_stream<SensorMetaInfo> {
  static constexpr std::string_view type_name = "SensorMetaInfo"sv;

  static void write(std::ostream&         os,
                    const SensorMetaInfo& x,
                    bofstream*            pbofs = nullptr,
                    std::string_view      name  = ""sv);

  static void read(std::istream& is, SensorMetaInfo& x, bifstream* pbifs = nullptr);
};

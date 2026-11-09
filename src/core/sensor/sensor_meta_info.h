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

/// Metadata describing the elements of a sensor as a structured gridded field.
///
/// The row-major data of the gridded field has one value per element of the
/// sensor, in the order of its elements.  Each variant alternative corresponds
/// to a different sensor type.  The gridded field's data_name serves as a
/// human-readable label (e.g., "camera", "ffts").  Grids are set at sensor
/// construction time; data is zeroed.
using MetaInfo = std::variant<SortedGriddedField1, SortedGriddedField2, CameraGriddedField, SortedGriddedField4>;

/// Total number of scalar elements (product of all grid sizes).
[[nodiscard]] Size count(const MetaInfo& meta);
}  // namespace sensor

using SensorMetaInfo = sensor::MetaInfo;

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

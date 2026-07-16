#include "sensor_meta_info.h"

#include <xml.h>

#include <stdexcept>

namespace sensor {

Size MetaInfo::count() const {
  return std::visit(
      [](const auto& gf) -> Size {
        Size n = 1;
        for (auto s : gf.shape()) n *= s;
        return n;
      },
      data);
}

void MetaInfo::check() const {
  std::visit(
      [](const auto& gf) {
        ARTS_USER_ERROR_IF(not gf.ok(),
                           "SensorMetaInfo gridded field shape mismatch: data shape "
                           "does not match grid sizes");
      },
      data);
}

}  // namespace sensor

//
// XML I/O — custom dispatch since all gridded_data_t share the same
// xml_io_stream type_name ("GriddedField"), which prevents using the
// generic variant XML machinery.
//

// String tags for discriminating the variant alternative in XML.
static constexpr std::string_view tag_sorted_gf1 = "SortedGriddedField1";
static constexpr std::string_view tag_camera     = "CameraGriddedField";

void xml_io_stream<SensorMetaInfo>::write(std::ostream&         os,
                                          const SensorMetaInfo& x,
                                          bofstream*            pbofs,
                                          std::string_view      name) {
  std::string_view type_tag;
  if (std::holds_alternative<SortedGriddedField1>(x.data))
    type_tag = tag_sorted_gf1;
  else if (std::holds_alternative<sensor::CameraGriddedField>(x.data))
    type_tag = tag_camera;

  XMLTag tag{type_name, "name", name, "type", type_tag};
  tag.write_to_stream(os);

  std::visit([&](const auto& gf) { xml_write_to_stream(os, gf, pbofs); }, x.data);

  tag.write_to_end_stream(os);
}

void xml_io_stream<SensorMetaInfo>::read(std::istream& is, SensorMetaInfo& x, bifstream* pbifs) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String type_tag;
  tag.get_attribute_value("type", type_tag);

  if (type_tag == tag_sorted_gf1) {
    SortedGriddedField1 gf;
    xml_read_from_stream(is, gf, pbifs);
    x.data = std::move(gf);
  } else if (type_tag == tag_camera) {
    sensor::CameraGriddedField gf;
    xml_read_from_stream(is, gf, pbifs);
    x.data = std::move(gf);
  } else {
    throw std::runtime_error(std::format("Unknown SensorMetaInfo type: \"{}\"", type_tag));
  }

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("Error reading SensorMetaInfo:\n{}", e.what()));
}

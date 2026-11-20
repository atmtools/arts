#include "sensor_meta_info.h"

#include <xml.h>

#include <array>
#include <stdexcept>

namespace sensor {
Size count(const MetaInfo& meta) {
  return std::visit(
      [](const auto& gf) -> Size {
        Size n = 1;
        for (auto s : gf.shape()) n *= s;
        return n;
      },
      meta);
}
}  // namespace sensor

//
// XML I/O — custom dispatch since all gridded_data_t share the same
// xml_io_stream type_name ("GriddedField"), which prevents using the
// generic variant XML machinery.
//

namespace {
// String tags for discriminating the variant alternative in XML, in variant order.
constexpr std::array<std::string_view, std::variant_size_v<SensorMetaInfo>> type_tags{
    "SortedGriddedField1", "SortedGriddedField2", "CameraGriddedField", "SortedGriddedField4"};

template <Size I = 0>
void read_alternative(std::string_view type_tag, std::istream& is, SensorMetaInfo& x, bifstream* pbifs) {
  if constexpr (I == std::variant_size_v<SensorMetaInfo>) {
    throw std::runtime_error(std::format("Unknown SensorMetaInfo type: \"{}\"", type_tag));
  } else if (type_tag == type_tags[I]) {
    std::variant_alternative_t<I, SensorMetaInfo> gf;
    xml_read_from_stream(is, gf, pbifs);
    x = std::move(gf);
  } else {
    read_alternative<I + 1>(type_tag, is, x, pbifs);
  }
}
}  // namespace

void xml_io_stream<SensorMetaInfo>::write(std::ostream&         os,
                                          const SensorMetaInfo& x,
                                          bofstream*            pbofs,
                                          std::string_view      name) {
  XMLTag tag{type_name, "name", name, "type", type_tags[x.index()]};
  tag.write_to_stream(os);

  std::visit([&](const auto& gf) { xml_write_to_stream(os, gf, pbofs); }, x);

  tag.write_to_end_stream(os);
}

void xml_io_stream<SensorMetaInfo>::read(std::istream& is, SensorMetaInfo& x, bifstream* pbifs) try {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  String type_tag;
  tag.get_attribute_value("type", type_tag);

  read_alternative(type_tag, is, x, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
} catch (const std::exception& e) {
  throw std::runtime_error(std::format("Error reading SensorMetaInfo:\n{}", e.what()));
}

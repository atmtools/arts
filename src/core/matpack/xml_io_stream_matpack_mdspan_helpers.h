#pragma once

#include <xml.h>

#include "matpack_mdspan_helpers.h"
#include "xml_io_stream_matpack_mdspan.h"

template <class Compare>
struct xml_io_stream<matpack::grid_t<Compare>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<matpack::data_t<Numeric, 1>>;

  static void write(std::ostream& os,
                    const matpack::grid_t<Compare>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    xml_io_stream<Vector>::write(os, x.vec(), pbofs, name);
  }

  static void read(std::istream& is,
                   matpack::grid_t<Compare>& x,
                   bifstream* pbifs = nullptr) {
    Vector y;
    xml_io_stream<Vector>::read(is, y, pbifs);
    x = y;
  }
};

template <typename T, typename... Grids>
struct xml_io_stream_name<matpack::gridded_data_t<T, Grids...>> {
  static constexpr std::string_view name = "GriddedField"sv;
};

template <arts_xml_ioable T, arts_xml_ioable... Grids>
struct xml_io_stream_gridded_field {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<matpack::gridded_data_t<T, Grids...>>;
  constexpr static Size codeversion = 2;

  static void write(std::ostream& os,
                    const matpack::gridded_data_t<T, Grids...>& gf,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    constexpr Size M = sizeof...(Grids);
    XMLTag tag{type_name, "name", name, "N", M, "version", codeversion};
    tag.write_to_stream(os);

    xml_write_to_stream(os, gf.data_name, pbofs, "Name of Data");
    xml_write_to_stream(os, gf.grid_names, pbofs, "Gridnames");
    xml_write_to_stream(os, gf.grids, pbofs, "Grids");
    xml_write_to_stream(os, gf.data, pbofs, "Data");

    tag.write_to_end_stream(os);
  }

  static void parse(std::istream& is,
                    matpack::gridded_data_t<T, Grids...>& gf,
                    bifstream* pbifs,
                    XMLTag& tag) {
    constexpr Size M = sizeof...(Grids);
    tag.check_name(type_name);
    tag.check_attribute("N", M);

    xml_read_from_stream(is, gf.data_name, pbifs);
    xml_read_from_stream(is, gf.grid_names, pbifs);
    xml_read_from_stream(is, gf.grids, pbifs);
    xml_read_from_stream(is, gf.data, pbifs);

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }

  static void read(std::istream& is,
                   matpack::gridded_data_t<T, Grids...>& gf,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);

    Index version = 0;
    if (tag.has_attribute("version"))
      tag.get_attribute_value("version", version);

    if (version <= 0) throw std::runtime_error("Old version, please update");

    parse(is, gf, pbifs, tag);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error reading {}:\n{}", type_name.data(), e.what()));
  }
};

template <arts_xml_ioable T, arts_xml_ioable... Grids>
struct xml_io_stream<matpack::gridded_data_t<T, Grids...>> {
  using GF = xml_io_stream_gridded_field<T, Grids...>;

  static constexpr std::string_view type_name = GF::type_name;

  static void write(std::ostream& os,
                    const matpack::gridded_data_t<T, Grids...>& gf,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    GF::write(os, gf, pbofs, name);
  }

  static void read(std::istream& is,
                   matpack::gridded_data_t<T, Grids...>& gf,
                   bifstream* pbifs = nullptr) {
    GF::read(is, gf, pbifs);
  }
};

#define GF_IO(GF)               \
  template <>                   \
  void xml_io_stream<GF>::read( \
      std::istream& is_xml, GF& gf, bifstream* pbifs);

GF_IO(GriddedField1)
GF_IO(GriddedField2)
GF_IO(GriddedField3)
GF_IO(GriddedField4)
GF_IO(GriddedField5)
GF_IO(GriddedField6)
GF_IO(ComplexGriddedField2)
GF_IO(SortedGriddedField1)
GF_IO(SortedGriddedField2)
GF_IO(SortedGriddedField3)
GF_IO(SortedGriddedField4)
GF_IO(SortedGriddedField5)
GF_IO(SortedGriddedField6)
GF_IO(CartesianSubsurfaceGriddedField3)
#undef GF_IO

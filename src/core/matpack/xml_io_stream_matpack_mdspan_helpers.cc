#include "xml_io_stream_matpack_mdspan_helpers.h"

namespace {
template <Size N, typename T, typename... Grids, Size M = sizeof...(Grids)> void xml_read_from_stream_recursive_old(
    std::istream& is_xml, matpack::gridded_data_t<T, Grids...>& gfield, bifstream* pbifs, XMLTag& tag) {
  if constexpr (M > N) {
    xml_read_from_stream(is_xml, gfield.template grid<N>(), pbifs);
    xml_read_from_stream_recursive_old<N + 1>(is_xml, gfield, pbifs, tag);
  }
}

template <typename T, typename... Grids> void xml_read_from_stream_gf_old(std::istream&                         is_xml,
                                                                          matpack::gridded_data_t<T, Grids...>& gfield,
                                                                          bifstream*                            pbifs,
                                                                          Index version) {
  using GF = matpack::gridded_data_t<T, Grids...>;

  XMLTag tag;

  if (version == 0) {
    const auto reader = [&](auto& grid, String& name) {
      using U = std::decay_t<decltype(grid)>;

      tag.get_attribute_value("name", name);

      if constexpr (std::same_as<U, Vector>) {
        if (tag.name != "Vector") arts::user_error("Must be Vector, is {}", tag.name);
        old_xml_io_parse(is_xml, grid, pbifs, tag);
        tag.read_from_stream(is_xml);
        tag.check_name("/Vector"sv);
      } else if constexpr (std::same_as<U, ArrayOfString>) {
        if (tag.name != "Array") arts::user_error("Must be Array, is {}", tag.name);
        String s;
        tag.get_attribute_value("type", s);
        if (s != "String") arts::user_error("Must be Array<String>, is Array<{}>", s);
        xml_parse_from_stream(is_xml, grid, pbifs, tag);
        tag.read_from_stream(is_xml);
        tag.check_name("/Array"sv);
      } else {
        arts::user_error("Unknown grid type: {}", tag.name);
      }
    };

    if constexpr (constexpr Size N = 0; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 1; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 2; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 3; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 4; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 5; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 6; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 7; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 8; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    if constexpr (constexpr Size N = 9; GF::dim > N) {
      tag.read_from_stream(is_xml);
      reader(gfield.template grid<N>(), gfield.template gridname<N>());
    }

    xml_io_stream<matpack::data_t<T, sizeof...(Grids)>>::read(is_xml, gfield.data, pbifs);
  } else if (version == 1) {
    ArrayOfString gridnames;
    xml_read_from_stream(is_xml, gridnames, pbifs);
    if (gridnames.size() != gfield.grid_names.size()) arts::user_error("Bad number of grid names:\n{}", gridnames);
    stdr::move(gridnames, gfield.grid_names.begin());

    xml_read_from_stream_recursive_old<0>(is_xml, gfield, pbifs, tag);

    xml_read_from_stream(is_xml, gfield.data, pbifs);
  } else {
    std::unreachable();
  }
}

template <typename T, typename... Grids> void xml_read_from_stream_tmpl(std::istream&                         is_xml,
                                                                        matpack::gridded_data_t<T, Grids...>& gfield,
                                                                        bifstream*                            pbifs,
                                                                        const char* old_name) {
  XMLTag tag;
  tag.read_from_stream(is_xml);

  Index version = 0;
  if (tag.has_attribute("version")) tag.get_attribute_value("version", version);

  if (version <= 1) {
    tag.check_name(std::string_view{old_name});

    tag.get_attribute_value("name", gfield.data_name);

    xml_read_from_stream_gf_old(is_xml, gfield, pbifs, version);

    tag.read_from_stream(is_xml);
    tag.check_end_name(std::string_view{old_name});
  } else {
    xml_io_stream_gridded_field<T, Grids...>::parse(is_xml, gfield, pbifs, tag);
  }

  if (not gfield.ok()) arts::user_error("Bad gridded field:\n{}", gfield);
}
}  // namespace

#define GF_IO(GF)                                                                                    \
  template <> void xml_io_stream<GF>::read(std::istream& is_xml, GF& gfield, bifstream* pbifs) try { \
    xml_read_from_stream_tmpl(is_xml, gfield, pbifs, #GF);                                           \
  }                                                                                                  \
  ARTS_METHOD_ERROR_CATCH

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

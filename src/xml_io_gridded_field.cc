#include <workspace.h>

#include <utility>

#include "xml_io.h"
#include "xml_io_arts_types.h"
#include "xml_io_base.h"
#include "xml_io_general_types.h"

inline constexpr Index CurrentVersion = 2;

template <typename T, typename... Grids>
void xml_write_to_stream_tmpl(
    std::ostream& os_xml,
    const matpack::gridded_data_t<T, Grids...>& gfield,
    bofstream* pbofs,
    const String& name) {
  ArtsXMLTag open_tag;
  open_tag.set_name("GriddedField");
  open_tag.add_attribute("N", gfield.dim);
  if (not name.empty()) open_tag.add_attribute("Name", name);
  open_tag.add_attribute("version", CurrentVersion);
  open_tag.write_to_stream(os_xml);
  std::println(os_xml);

  xml_write_to_stream(os_xml, gfield.data_name, pbofs, "Name of Data");

  for (Size i = 0; i < gfield.dim; i++) {
    xml_write_to_stream(os_xml, gfield.grid_names[i], pbofs, "Gridname");
  }

  std::apply(
      [&os_xml, pbofs](auto&... vals) {
        (xml_write_to_stream(os_xml, vals, pbofs, "Grid"), ...);
      },
      gfield.grids);

  xml_write_to_stream(os_xml, gfield.data, pbofs, "Data");

  ArtsXMLTag close_tag;
  close_tag.set_name("/GriddedField");
  close_tag.write_to_stream(os_xml);
  std::println(os_xml);
}

template <Size N, typename T, typename... Grids, Size M = sizeof...(Grids)>
void xml_read_from_stream_recursive_old(
    std::istream& is_xml,
    matpack::gridded_data_t<T, Grids...>& gfield,
    bifstream* pbifs,
    XMLTag& tag) {
  if constexpr (M > N) {
    xml_read_from_stream(is_xml, gfield.template grid<N>(), pbifs);
    xml_read_from_stream_recursive_old<N + 1>(is_xml, gfield, pbifs, tag);
  }
}

template <typename T, typename... Grids>
void xml_read_from_stream_gf_old(std::istream& is_xml,
                                 matpack::gridded_data_t<T, Grids...>& gfield,
                                 bifstream* pbifs,
                                 Index version) {
  using GF = matpack::gridded_data_t<T, Grids...>;

  XMLTag tag;

  if (version == 0) {
    const auto reader = [&](auto& grid, String& name) {
      using U = std::decay_t<decltype(grid)>;

      tag.get_attribute_value("name", name);

      if constexpr (std::same_as<U, Vector>) {
        ARTS_USER_ERROR_IF(tag.get_name() != "Vector",
                           "Must be Vector, is {}",
                           tag.get_name());
        xml_parse_from_stream(is_xml, grid, pbifs, tag);
        tag.read_from_stream(is_xml);
        tag.check_name("/Vector");
      } else if constexpr (std::same_as<U, ArrayOfString>) {
        ARTS_USER_ERROR_IF(
            tag.get_name() != "Array", "Must be Array, is {}", tag.get_name());
        String s;
        tag.get_attribute_value("type", s);
        ARTS_USER_ERROR_IF(
            s != "String", "Must be Array<String>, is Array<{}>", s);
        xml_parse_from_stream(is_xml, grid, pbifs, tag);
        tag.read_from_stream(is_xml);
        tag.check_name("/Array");
      } else {
        ARTS_USER_ERROR("Unknown grid type: {}", tag.get_name());
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

    xml_read_from_stream(is_xml, gfield.data, pbifs);
  } else if (version == 1) {
    ArrayOfString gridnames;
    xml_read_from_stream(is_xml, gridnames, pbifs);
    ARTS_USER_ERROR_IF(gridnames.size() != gfield.grid_names.size(),
                       "Bad number of grid names:\n{}",
                       gridnames);
    std::ranges::move(gridnames, gfield.grid_names.begin());

    xml_read_from_stream_recursive_old<0>(is_xml, gfield, pbifs, tag);

    xml_read_from_stream(is_xml, gfield.data, pbifs);
  } else {
    std::unreachable();
  }
}

template <typename T, typename... Grids>
void xml_read_from_stream_tmpl(std::istream& is_xml,
                               matpack::gridded_data_t<T, Grids...>& gfield,
                               bifstream* pbifs,
                               const char* old_name) {
  ArtsXMLTag tag;
  Index N{};
  tag.read_from_stream(is_xml);

  Index version = 0;
  if (tag.has_attribute("version")) tag.get_attribute_value("version", version);

  if (version <= 1) {
    tag.check_name(old_name);

    tag.get_attribute_value("name", gfield.data_name);

    xml_read_from_stream_gf_old(is_xml, gfield, pbifs, version);

    tag.read_from_stream(is_xml);
    tag.check_name(std::string{"/"} + old_name);
  } else {
    tag.check_name("GriddedField");
    tag.get_attribute_value("N", N);
    const Size MyN = static_cast<Size>(N);

    ARTS_USER_ERROR_IF(MyN != gfield.dim,
                       "Bad dimension N := {}, expected N := {}",
                       MyN,
                       gfield.dim);

    xml_read_from_stream(is_xml, gfield.data_name, pbifs);

    for (Size i = 0; i < gfield.dim; i++) {
      xml_read_from_stream(is_xml, gfield.grid_names[i], pbifs);
    }

    std::apply(
        [&is_xml, pbifs](auto&... vals) {
          (xml_read_from_stream(is_xml, vals, pbifs), ...);
        },
        gfield.grids);

    xml_read_from_stream(is_xml, gfield.data, pbifs);

    tag.read_from_stream(is_xml);
    tag.check_name("/GriddedField");
  }

  ARTS_USER_ERROR_IF(not gfield.ok(), "Bad gridded field:\n{}", gfield);
}

#define GF_IO(GF)                                               \
  void xml_read_from_stream(                                    \
      std::istream& is_xml, GF& gfield, bifstream* pbifs) try { \
    xml_read_from_stream_tmpl(is_xml, gfield, pbifs, #GF);      \
  }                                                             \
  ARTS_METHOD_ERROR_CATCH                                       \
                                                                \
  void xml_write_to_stream(std::ostream& os_xml,                \
                           const GF& gfield,                    \
                           bofstream* pbofs,                    \
                           const String& name) try {            \
    xml_write_to_stream_tmpl(os_xml, gfield, pbofs, name);      \
  }                                                             \
  ARTS_METHOD_ERROR_CATCH

GF_IO(GriddedField1)
GF_IO(GriddedField2)
GF_IO(GriddedField3)
GF_IO(GriddedField4)
GF_IO(GriddedField5)
GF_IO(GriddedField6)
GF_IO(NamedGriddedField2)
GF_IO(NamedGriddedField3)
GF_IO(GriddedField1Named)
GF_IO(ComplexGriddedField2)
GF_IO(SortedGriddedField1)
GF_IO(SortedGriddedField2)
GF_IO(SortedGriddedField3)
GF_IO(SortedGriddedField4)
GF_IO(SortedGriddedField5)
GF_IO(SortedGriddedField6)
GF_IO(StokvecSortedGriddedField1)
GF_IO(StokvecSortedGriddedField2)
GF_IO(StokvecSortedGriddedField3)
GF_IO(StokvecSortedGriddedField4)
GF_IO(StokvecSortedGriddedField5)
GF_IO(StokvecSortedGriddedField6)
GF_IO(CartesianSubsurfaceGriddedField3)

#undef GF_IO

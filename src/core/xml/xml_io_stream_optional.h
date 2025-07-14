#pragma once

#include <optional>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename T>
struct xml_io_stream_name<std::optional<T>> {
  static constexpr std::string_view name = "Optional"sv;
};

template <arts_xml_ioable T>
  requires(std::is_default_constructible_v<T>)
struct xml_io_stream<std::optional<T>> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<std::optional<T>>;

  static void write(std::ostream& os,
                    const std::optional<T>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    std::println(os,
                 R"(<{0} name="{1}" type="{2}" null="{3}">)",
                 type_name,
                 name,
                 xml_io_stream<T>::type_name,
                 Index{bool{x}});

    if (x) xml_io_stream<T>::write(os, *x, pbofs, name);

    std::println(os, R"(</{0}>)", type_name);
  }

  static void read(std::istream& is,
                   std::optional<T>& x,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("type", xml_io_stream<T>::type_name);

    Index null;
    tag.get_attribute_value("null", null);

    if (null) {
      x = std::make_optional<T>();
      xml_io_stream<T>::read(is, *x, pbifs);
    } else {
      x = std::nullopt;
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::format("Error reading {}<{}>:\n{}",
                                         type_name,
                                         xml_io_stream<T>::type_name,
                                         e.what()));
  }
};

#pragma once

#include <enums.h>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename T>
concept enum_option = requires(T a) {
  to<T>(std::string{});
  toString(a);
  enumdocs<T>::name;
};

template <enum_option T>
struct xml_io_stream_name<T> {
  static constexpr std::string_view name = enumdocs<T>::name;
};

template <enum_option T>
struct xml_io_stream<T> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<T>;

  static void write(std::ostream& os,
                    const T& x,
                    bofstream*       = nullptr,
                    std::string_view = ""sv) {
    std::println(os, R"(<{0} value="{1}"> </{0}>)", type_name, toString(x));
  }

  static void read(std::istream& is, T& x, bifstream* = nullptr) {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    std::string value;
    tag.get_attribute_value("value", value);
    x = to<T>(value);

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
};

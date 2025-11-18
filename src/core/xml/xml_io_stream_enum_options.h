#pragma once

#include <enums.h>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename T>
concept xml_enum_option = requires(T a) {
  to<T>(std::string{});
  toString(a);
  enumdocs<T>::name;
};

template <xml_enum_option T>
struct xml_io_stream_name<T> {
  static constexpr std::string_view name = enumdocs<T>::name;
};

template <xml_enum_option T>
struct xml_io_stream<T> {
  static constexpr std::string_view type_name = xml_io_stream_name_v<T>;

  static void write(std::ostream& os,
                    const T& x,
                    bofstream*            = nullptr,
                    std::string_view name = ""sv) {
    XMLTag tag{type_name, "value", toString(x), "name", name};
    tag.write_to_stream(os);
    tag.write_to_end_stream(os);
  }

  static void read(std::istream& is, T& x, bifstream* = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);

    std::string value;
    tag.get_attribute_value("value", value);
    try {
      x = to<T>(value);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          std::format("Failed to convert '{}' to enum {}:\n{}",
                      value,
                      type_name,
                      e.what()));
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error reading enum {}:\n{}", type_name, e.what()));
  }
};

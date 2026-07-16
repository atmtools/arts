#pragma once

#include <concepts>
#include <string_view>

#include "xml_io_base.h"
#include "xml_io_stream.h"

using namespace std::literals;

template <xml_io_writable T, std::derived_from<T> U> requires(not std::same_as<T, U>) struct xml_io_stream_inherit {
  static constexpr std::string_view type_name = xml_io_stream_name_v<T>;

  struct _no_name_ {};
  static_assert(type_name != xml_io_stream_name_v<_no_name_>,
                "xml_io_stream_inherit requires that the base type has a name");

  static void write(std::ostream& os, const U& n, bofstream* pbofs = nullptr, std::string_view name = ""sv) {
    if constexpr (type_name != xml_io_stream_name_v<T>) {
      XMLTag tag(type_name, "name", name);
      tag.write_to_stream(os);
      xml_io_stream<T>::write(os, reinterpret_cast<const T&>(n), pbofs);
      tag.write_to_end_stream(os);
    } else {
      xml_io_stream<T>::write(os, reinterpret_cast<const T&>(n), pbofs, name);
    }
  }

  static void read(std::istream& is, U& n, bifstream* pbifs = nullptr) {
    if constexpr (type_name != xml_io_stream_name_v<T>) {
      XMLTag tag{};
      tag.read_from_stream(is);
      tag.check_name(type_name);
    }

    xml_io_stream<T>::read(is, reinterpret_cast<T&>(n), pbifs);

    if constexpr (type_name != xml_io_stream_name_v<T>) {
      XMLTag tag{};
      tag.read_from_stream(is);
      tag.check_end_name(type_name);
    }
  }

  static void get(std::span<U> b, bifstream* pbifs) requires(xml_io_binary<T>) {
    xml_io_stream<T>::get(std::span{reinterpret_cast<T*>(b.data()), b.size()}, pbifs);
  }

  static void put(std::span<const U> b, bofstream* pbofs) requires(xml_io_binary<T>) {
    xml_io_stream<T>::put(std::span{reinterpret_cast<const T*>(b.data()), b.size()}, pbofs);
  }

  static void parse(std::span<U> b, std::istream& is) requires(xml_io_parseable<T>) {
    xml_io_stream<T>::parse(std::span{reinterpret_cast<T*>(b.data()), b.size()}, is);
  }
};

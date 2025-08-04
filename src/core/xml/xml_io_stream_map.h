#pragma once

#include <unordered_map>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename Key, typename T>
struct xml_io_stream_name<std::unordered_map<Key, T>> {
  static constexpr std::string_view name = "Map"sv;
};

template <arts_xml_ioable Key, arts_xml_ioable T>
  requires(std::is_default_constructible_v<Key> and
           std::is_default_constructible_v<T>)
struct xml_io_stream<std::unordered_map<Key, T>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::unordered_map<Key, T>>;

  static void write(std::ostream& os,
                    const std::unordered_map<Key, T>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) {
    XMLTag tag{type_name,
               "name",
               name,
               "type",
               xml_io_stream<T>::type_name,
               "key",
               xml_io_stream<Key>::type_name,
               "nelem",
               n.size()};
    tag.write_to_stream(os);

    for (const auto& [key, elem] : n) {
      xml_io_stream<Key>::write(os, key, pbofs);
      xml_io_stream<T>::write(os, elem, pbofs);
    }

    tag.write_to_end_stream(os);
  }

  static void read(std::istream& is,
                   std::unordered_map<Key, T>& n,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("key", xml_io_stream<Key>::type_name);
    tag.check_attribute("type", xml_io_stream<T>::type_name);

    Size nelem = 0;
    tag.get_attribute_value("nelem", nelem);
    n.clear();
    n.reserve(nelem);

    for (Size i = 0; i < nelem; ++i) {
      Key key{};
      xml_io_stream<Key>::read(is, key, pbifs);
      xml_io_stream<T>::read(is, n[key], pbifs);
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::format("Error reading {}<{}, {}>:\n{}",
                                         type_name,
                                         xml_io_stream<Key>::type_name,
                                         xml_io_stream<T>::type_name,
                                         e.what()));
  }
};

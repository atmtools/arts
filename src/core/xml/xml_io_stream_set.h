#pragma once

#include <set>
#include <string_view>

#include "xml_io_base.h"
#include "xml_io_stream.h"

using namespace std::literals;

template <typename T> struct xml_io_stream_name<std::set<T>> {
  static constexpr std::string_view name = "Set"sv;
};

template <arts_xml_ioable T> requires(std::is_default_constructible_v<T>) struct xml_io_stream<std::set<T>> {
  constexpr static std::string_view type_name = xml_io_stream_name_v<std::set<T>>;

  static void write(std::ostream& os, const std::set<T>& n, bofstream* pbofs = nullptr, std::string_view name = ""sv) {
    XMLTag tag{type_name, "name", name, "type", xml_io_stream<T>::type_name, "nelem", n.size()};
    tag.write_to_stream(os);

    for (const auto& elem : n) { xml_io_stream<T>::write(os, elem, pbofs); }

    tag.write_to_end_stream(os);
  }

  static void extend(std::istream& is, std::set<T>& n, bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("type", xml_io_stream<T>::type_name);

    Size nelem = 0;
    tag.get_attribute_value("nelem", nelem);

    for (Size i = 0; i < nelem; ++i) {
      T elem{};
      xml_io_stream<T>::read(is, elem, pbifs);
      n.insert(std::move(elem));
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Cannot extend {}<{}>:\n{}", type_name, xml_io_stream<T>::type_name, e.what()));
  }

  static void read(std::istream& is, std::set<T>& n, bifstream* pbifs = nullptr) try {
    n.clear();
    extend(is, n, pbifs);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::format("Cannot read {}<{}>:\n{}", type_name, xml_io_stream<T>::type_name, e.what()));
  }
};

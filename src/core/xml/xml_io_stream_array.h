#pragma once

#include <array.h>
#include <double_imanip.h>

#include "xml_io_base.h"
#include "xml_io_stream.h"
#include "xml_io_stream_core.h"

template <typename T>
struct xml_io_stream_name<Array<T>> {
  static constexpr std::string_view name = "Array"sv;
};

template <>
struct xml_io_stream_name<ArrayOfString> {
  static constexpr std::string_view name = "ArrayOfString"sv;
};

template <arts_xml_ioable T>
  requires(std::is_default_constructible_v<T>)
struct xml_io_stream<Array<T>> {
  constexpr static std::string_view type_name = xml_io_stream_name_v<Array<T>>;

  using inner = xml_io_stream<T>;

  static void write(std::ostream& os,
                    const Array<T>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(os,
                 R"(<{0} type="{1}" name="{2}" nelem="{3}">)",
                 type_name,
                 inner::type_name,
                 name,
                 n.size());

    if (pbofs) {
      if constexpr (xml_io_binary<T>)
        inner::put(std::span{n.data(), n.size()}, pbofs);
      else
        for (auto& v : n) inner::write(os, v, pbofs);
    } else {
      if constexpr (xml_io_parseable<T>)
        std::println(os, "{:IO}", n);
      else
        for (auto& v : n) inner::write(os, v, pbofs);
    }

    std::println(os, R"(</{0}>)", type_name);
  }
  ARTS_METHOD_ERROR_CATCH

  static void read(std::istream& is,
                   Array<T>& n,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("type", inner::type_name);

    Size nelem = 0;
    tag.get_attribute_value("nelem", nelem);
    n.resize(nelem);

    if (pbifs) {
      if constexpr (xml_io_binary<T>)
        inner::get(std::span{n.data(), n.size()}, pbifs);
      else
        for (auto& v : n) inner::read(is, v, pbifs);
    } else {
      if constexpr (xml_io_parseable<T>)
        inner::parse(std::span{n.data(), n.size()}, is);
      else
        for (auto& v : n) inner::read(is, v, pbifs);
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
  ARTS_METHOD_ERROR_CATCH
};

template <typename T, Size N>
struct xml_io_stream_name<std::array<T, N>> {
  static constexpr std::string_view name = xml_io_stream_name_v<Array<T>>;
};

template <arts_xml_ioable T, Size N>
struct xml_io_stream<std::array<T, N>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<std::array<T, N>>;

  using inner = xml_io_stream<T>;

  static void write(std::ostream& os,
                    const std::array<T, N>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(os,
                 R"(<{0} type="{1}" name="{2}" nelem="{3}">)",
                 type_name,
                 inner::type_name,
                 name,
                 N);

    if (pbofs) {
      if constexpr (xml_io_binary<T>)
        put(std::span{&n, 1}, pbofs);
      else
        for (auto& v : n) inner::write(os, v, pbofs);
    } else {
      if constexpr (xml_io_parseable<T>)
        std::println(os, "{:IO}", n);
      else
        for (auto& v : n) inner::write(os, v, pbofs);
    }

    std::println(os, R"(</{0}>)", type_name);
  }
  ARTS_METHOD_ERROR_CATCH

  static void read(std::istream& is,
                   std::array<T, N>& n,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("type", inner::type_name);

    Size nelem = 0;
    tag.get_attribute_value("nelem", nelem);
    if (N != nelem) throw std::runtime_error("Size-mismatch for constant type");

    if (pbifs) {
      if constexpr (xml_io_binary<T>)
        get(std::span{&n, 1}, pbifs);
      else
        for (auto& v : n) inner::read(is, v, pbifs);
    } else {
      if constexpr (xml_io_parseable<T>)
        parse(std::span{&n, 1}, is);
      else
        for (auto& v : n) inner::read(is, v, pbifs);
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
  ARTS_METHOD_ERROR_CATCH

  static void put(std::span<const std::array<T, N>> v, bifstream* pbifs)
    requires(xml_io_binary<T>)
  {
    inner::put(std::span(reinterpret_cast<const T*>(v.data()), N), pbifs);
  }

  static void get(std::span<std::array<T, N>> v, bofstream* pbofs)
    requires(xml_io_binary<T>)
  {
    inner::get(std::span(reinterpret_cast<T*>(v.data()), N), pbofs);
  }

  static void parse(std::span<std::array<T, N>> v, std::istream& is)
    requires(xml_io_parseable<T>)
  {
    inner::parse(std::span(reinterpret_cast<T*>(v.data()), v.size()), is);
  }
};

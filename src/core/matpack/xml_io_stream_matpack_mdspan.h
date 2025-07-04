#pragma once

#include <xml.h>

#include "matpack_mdspan.h"

void old_xml_io_read(XMLTag&, std::istream&, ComplexVector&, bifstream*);

void old_xml_io_parse(std::istream&, Vector&, bifstream*, XMLTag&);
void old_xml_io_read(XMLTag&, std::istream&, Vector&, bifstream*);
void old_xml_io_read(XMLTag&, std::istream&, Matrix&, bifstream*);
void old_xml_io_read(XMLTag&, std::istream&, Tensor3&, bifstream*);
void old_xml_io_read(XMLTag&, std::istream&, Tensor4&, bifstream*);
void old_xml_io_read(XMLTag&, std::istream&, Tensor5&, bifstream*);
void old_xml_io_read(XMLTag&, std::istream&, Tensor6&, bifstream*);
void old_xml_io_read(XMLTag&, std::istream&, Tensor7&, bifstream*);

template <typename T, Size N>
struct xml_io_stream_name<matpack::data_t<T, N>> {
  static constexpr std::string_view name = "Matpack"sv;
};

template <>
struct xml_io_stream_name<std::shared_ptr<Matrix>> {
  static constexpr std::string_view name = "SharedMatrix"sv;
};

template <arts_xml_ioable T, Size N>
struct xml_io_stream<matpack::data_t<T, N>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<matpack::data_t<T, N>>;
  static constexpr Size codeversion = 2;

  using inner = xml_io_stream<T>;

  static void write(std::ostream& os,
                    const matpack::data_t<T, N>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(
        os,
        R"(<{0} rank="{1}" type="{2}" shape="{3}" name="{4}" version="{5}">)",
        type_name,
        N,
        inner::type_name,
        x.shape(),
        name,
        codeversion);

    if (pbofs) {
      if constexpr (xml_io_binary<T>)
        inner::put(std::span{x.data_handle(), x.size()}, pbofs);
      else
        for (auto& v : elemwise_range(x)) inner::write(os, v, pbofs);
    } else {
      if constexpr (xml_io_parseable<T>)
        std::println(os, "{:IO}", x);
      else
        for (auto& v : elemwise_range(x)) inner::write(os, v, pbofs);
    }

    std::println(os, R"(</{0}>)", type_name);
  }
  ARTS_METHOD_ERROR_CATCH

  static void read(std::istream& is,
                   matpack::data_t<T, N>& x,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);

    Size fileversion = 1;  // Orig gets version number 1
    if (tag.has_attribute("version"))
      tag.get_attribute_value("version", fileversion);

    if (fileversion > codeversion or fileversion == 0) {
      throw std::runtime_error("Unknown version");
    } else if (fileversion == codeversion) {
      tag.check_name(type_name);

      String str_shape{};
      tag.get_attribute_value("shape", str_shape);

      tag.check_attribute("type", inner::type_name);
      tag.check_attribute("rank", N);

      std::istringstream iss(str_shape);
      std::array<Index, N> shape{};
      for (Size i = 0; i < N; ++i) iss >> shape[i];

      x.resize(shape);

      if (pbifs) {
        if constexpr (xml_io_binary<T>)
          inner::get(std::span{x.data_handle(), x.size()}, pbifs);
        else
          for (auto& v : elemwise_range(x)) inner::read(is, v, pbifs);
      } else {
        if constexpr (xml_io_parseable<T>)
          inner::parse(std::span{x.data_handle(), x.size()}, is);
        else
          for (auto& v : elemwise_range(x)) inner::read(is, v, pbifs);
      }

      tag.read_from_stream(is);
      tag.check_end_name(type_name);
    } else if (fileversion == 1) {
      if constexpr (requires { old_xml_io_read(tag, is, x, pbifs); }) {
        old_xml_io_read(tag, is, x, pbifs);
      } else {
        throw std::runtime_error("Unsupported old file format.");
      }
    }
  }
  ARTS_METHOD_ERROR_CATCH
};

template <typename T, Size... N>
struct xml_io_stream_name<matpack::cdata_t<T, N...>> {
  static constexpr std::string_view name =
      xml_io_stream_name_v<matpack::data_t<T, sizeof...(N)>>;
};

template <arts_xml_ioable T, Size... N>
struct xml_io_stream<matpack::cdata_t<T, N...>> {
  constexpr static std::string_view type_name =
      xml_io_stream_name_v<matpack::cdata_t<T, N...>>;

  using inner = xml_io_stream<T>;

  static void parse(std::span<matpack::cdata_t<T, N...>> v, std::istream& is)
    requires(xml_io_parseable<T>)
  {
    inner::parse(
        std::span{reinterpret_cast<T*>(v.data()), v.size() * (N * ...)}, is);
  }

  static void get(std::span<matpack::cdata_t<T, N...>> v, bifstream* pbifs)
    requires(xml_io_binary<T>)
  {
    inner::get(std::span{reinterpret_cast<T*>(v.data()), v.size() * (N * ...)},
               pbifs);
  }

  static void put(std::span<const matpack::cdata_t<T, N...>> v,
                  bofstream* pbofs)
    requires(xml_io_binary<T>)
  {
    inner::put(
        std::span{reinterpret_cast<const T*>(v.data()), v.size() * (N * ...)},
        pbofs);
  }

  static void write(std::ostream& os,
                    const matpack::cdata_t<T, N...>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(os,
                 R"(<{0} rank="{1}" type="{2}" shape="{3}" name="{4}">)",
                 type_name,
                 sizeof...(N),
                 inner::type_name,
                 x.shape(),
                 name);

    if (pbofs) {
      if constexpr (xml_io_binary<T>)
        put(std::span{&x, 1}, pbofs);
      else
        for (auto& v : x.data) inner::write(os, v, pbofs);
    } else {
      if constexpr (xml_io_parseable<T>)
        std::println(os, R"({:IO})", x);
      else
        for (auto& v : x.data) inner::write(os, v, pbofs);
    }

    std::println(os, R"(</{0}>)", type_name);
  }
  ARTS_METHOD_ERROR_CATCH

  static void read(std::istream& is,
                   matpack::cdata_t<T, N...>& x,
                   bifstream* pbifs = nullptr) try {
    XMLTag tag;
    tag.read_from_stream(is);
    tag.check_name(type_name);
    tag.check_attribute("type", inner::type_name);
    tag.check_attribute("rank", sizeof...(N));
    tag.check_attribute("shape", std::format("{}", x.shape()));

    if (pbifs) {
      if constexpr (xml_io_binary<T>)
        get(std::span{&x, 1}, pbifs);
      else
        for (auto& elem : x.data) inner::read(is, elem, pbifs);
    } else {
      if constexpr (xml_io_parseable<T>)
        parse(std::span{&x, 1}, is);
      else
        for (auto& elem : x.data) inner::read(is, elem, pbifs);
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
  ARTS_METHOD_ERROR_CATCH
};

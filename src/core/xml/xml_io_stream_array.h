#pragma once

#include <array.h>
#include <double_imanip.h>
#include <xml_io_base.h>
#include <xml_io_stream.h>
#include <xml_io_stream_core.h>

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

  static void write(std::ostream& os,
                    const Array<T>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(os,
                 R"(<{0} type="{1}" name="{2}" nelem="{3}">)",
                 type_name,
                 xml_io_stream<T>::type_name,
                 name,
                 n.size());

    if constexpr (xml_io_binary<T>) {
      if (pbofs) xml_io_stream<T>::put(n.data(), pbofs, n.size());
    } else {
      pbofs = nullptr;
    }

    if (not pbofs) {
      if constexpr (bitshift_readable<T> or xml_coretype<T>) {
        std::println(os, "{:IO}", n);
      } else {
        for (auto& v : n) xml_io_stream<T>::write(os, v);
      }
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

    Size nelem = 0;
    tag.get_attribute_value("nelem", nelem);
    n.resize(nelem);

    if constexpr (xml_io_binary<T>) {
      if (pbifs) xml_io_stream<T>::get(n.data(), pbifs, n.size());
    } else {
      pbifs = nullptr;
    }

    if (not pbifs) {
      if constexpr (bitshift_readable<T> or xml_coretype<T>) {
        for (auto& v : n) is >> v;
      } else {
        for (auto& v : n) xml_io_stream<T>::read(is, v);
      }
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

  static void write(std::ostream& os,
                    const std::array<T, N>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(os,
                 R"(<{0} type="{1}" name="{2}" nelem="{3}">)",
                 type_name,
                 xml_io_stream<T>::type_name,
                 name,
                 N);

    if constexpr (xml_io_binary<T>) {
      if (pbofs) xml_io_stream<T>::put(n.data(), pbofs, n.size());
    } else {
      pbofs = nullptr;
    }

    if (not pbofs) {
      if constexpr (bitshift_readable<T> or xml_coretype<T>) {
        std::println(os, "{:IO}", n);
      } else {
        for (auto& v : n) xml_io_stream<T>::write(os, v);
      }
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

    Size nelem = 0;
    tag.get_attribute_value("nelem", nelem);
    if (N != nelem) throw std::runtime_error("Size-mismatch for constant type");

    if constexpr (xml_io_binary<T>) {
      if (pbifs) xml_io_stream<T>::get(n.data(), pbifs, n.size());
    } else {
      pbifs = nullptr;
    }

    if (not pbifs) {
      if constexpr (bitshift_readable<T> or xml_coretype<T>) {
        for (auto& v : n) is >> v;
      } else {
        for (auto& v : n) xml_io_stream<T>::read(is, v);
      }
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
  ARTS_METHOD_ERROR_CATCH
};

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

namespace {
template <arts_xml_ioable T, Size N>
consteval std::string_view matpack_names() {
  if constexpr (std::same_as<matpack::data_t<T, N>, Vector>) {
    return "Vector"sv;
  } else if constexpr (std::same_as<matpack::data_t<T, N>, Matrix>) {
    return "Matrix"sv;
  } else if constexpr (std::same_as<matpack::data_t<T, N>, Tensor3>) {
    return "Tensor3"sv;
  } else if constexpr (std::same_as<matpack::data_t<T, N>, Tensor4>) {
    return "Tensor4"sv;
  } else if constexpr (std::same_as<matpack::data_t<T, N>, Tensor5>) {
    return "Tensor5"sv;
  } else if constexpr (std::same_as<matpack::data_t<T, N>, Tensor6>) {
    return "Tensor6"sv;
  } else if constexpr (std::same_as<matpack::data_t<T, N>, Tensor7>) {
    return "Tensor7"sv;
  } else {
    return "Matpack"sv;
  }
}
}  // namespace

template <arts_xml_ioable T, Size N>
struct xml_io_stream<matpack::data_t<T, N>> {
  constexpr static std::string_view type_name = matpack_names<T, N>();
  static constexpr Size codeversion           = 2;

  static void write(std::ostream& os,
                    const matpack::data_t<T, N>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(
        os,
        R"(<{0} rank="{1}" type="{2}" shape="{3}" name="{4}" version="{5}">)",
        type_name,
        N,
        xml_io_stream<T>::type_name,
        x.shape(),
        name,
        codeversion);

    if constexpr (xml_io_binary<T>) {
      if (pbofs) xml_io_stream<T>::put(x.data_handle(), pbofs, x.size());
    } else {
      pbofs = nullptr;
    }

    if (not pbofs) {
      if constexpr (bitshift_readable<T> or xml_coretype<T>) {
        std::println(os, "{:IO}", x);
      } else {
        for (auto& v : elemwise_range(x)) xml_io_stream<T>::write(os, v);
      }
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
      throw std::runtime_error(std::format(
          "Got version {}, compiled version is {}, minimum version is 1",
          fileversion,
          codeversion));
    } else if (fileversion == codeversion) {
      tag.check_name(type_name);

      String str_shape{};
      tag.get_attribute_value("shape", str_shape);

      String type;
      tag.get_attribute_value("type", type);
      if (type != xml_io_stream<T>::type_name)
        throw std::runtime_error("Type mismatch");

      Size rank = 0;
      tag.get_attribute_value("rank", rank);
      if (rank != N) throw std::runtime_error("Rank mismatch");

      std::istringstream iss(str_shape);
      std::array<Index, N> shape{};
      for (Size i = 0; i < N; ++i) iss >> shape[i];

      x.resize(shape);

      if constexpr (xml_io_binary<T>) {
        if (pbifs) xml_io_stream<T>::get(x.data_handle(), pbifs, x.size());
      } else {
        pbifs = nullptr;
      }

      if (not pbifs) {
        if constexpr (bitshift_readable<T> or xml_coretype<T>) {
          if constexpr (std::same_as<T, Numeric>) {
            for (auto& v : elemwise_range(x)) is >> double_imanip() >> v;
          } else if constexpr (std::same_as<T, Complex>) {
            for (auto& v : elemwise_range(x))
              is >> double_imanip() >> real_val(v) >> imag_val(v);
          } else {
            for (auto& v : elemwise_range(x)) is >> v;
          }
        } else {
          for (auto& v : elemwise_range(x)) xml_io_stream<T>::read(is, v);
        }
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

namespace {
// Make IO-compatible with data_t for many types
template <arts_xml_ioable T, Size... N>
consteval std::string_view matpack_cnames() {
  if constexpr (std::same_as<T, Numeric> and sizeof...(N) == 1) {
    return "Vector"sv;
  } else if constexpr (std::same_as<T, Numeric> and sizeof...(N) == 2) {
    return "Matrix"sv;
  } else if constexpr (std::same_as<T, Complex> and sizeof...(N) == 1) {
    return "ComplexVector"sv;
  } else if constexpr (std::same_as<T, Complex> and sizeof...(N) == 2) {
    return "ComplexMatrix"sv;
  } else {
    return "Matpack"sv;
  }
}
}  // namespace

template <arts_xml_ioable T, Size... N>
struct xml_io_stream<matpack::cdata_t<T, N...>> {
  constexpr static std::string_view type_name = matpack_cnames<T, N...>();

  static void get(matpack::cdata_t<T, N...>* v, bifstream* pbifs, Size n = 1)
    requires(xml_io_binary<T>)
  {
    assert(pbifs);
    xml_io_stream<T>::get(v->data.data(), pbifs, n * (N * ...));
  }

  static void put(const matpack::cdata_t<T, N...>* const v,
                  bofstream* pbofs,
                  Size n = 1)
    requires(xml_io_binary<T>)
  {
    assert(pbofs);
    xml_io_stream<T>::put(v->data.data(), pbofs, n * (N * ...));
  }

  static void write(std::ostream& os,
                    const matpack::cdata_t<T, N...>& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv) try {
    std::println(os,
                 R"(<{0} rank="{1}" type="{2}" shape="{3}" name="{4}">)",
                 type_name,
                 sizeof...(N),
                 xml_io_stream<T>::type_name,
                 x.shape(),
                 name);

    if constexpr (xml_io_binary<T>) {
      if (pbofs) put(&x, pbofs);
    } else {
      pbofs = nullptr;
    }

    if (not pbofs) {
      if constexpr (bitshift_readable<T> or xml_coretype<T>) {
        std::println(os, R"({:IO})", x);
      } else {
        for (auto& v : x.data) xml_io_stream<T>::write(os, v);
      }
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

    String str{};
    tag.get_attribute_value("shape", str);

    std::istringstream iss{str};
    std::array<Index, sizeof...(N)> shape{};
    for (Size i = 0; i < sizeof...(N); ++i) iss >> shape[i];
    if (shape != x.shape()) throw std::runtime_error("Mismatch size");

    tag.get_attribute_value("type", str);
    if (xml_io_stream<T>::type_name != str)
      throw std::runtime_error("Mismatch type");

    if constexpr (xml_io_binary<T>) {
      if (pbifs) get(&x, pbifs);
    } else {
      pbifs = nullptr;
    }

    if (not pbifs) {
      if constexpr (std::same_as<T, Numeric>) {
        for (auto& elem : elemwise_range(x)) is >> double_imanip() >> elem;
      } else if constexpr (std::same_as<T, Complex>) {
        for (auto& e : elemwise_range(x))
          is >> double_imanip() >> real_val(e) >> imag_val(e);
      } else {
        if constexpr (bitshift_readable<T> or xml_coretype<T>) {
          for (auto& elem : elemwise_range(x)) is >> elem;
        } else {
          for (auto& elem : elemwise_range(x)) xml_io_stream<T>::read(is, elem);
        }
      }
    }

    tag.read_from_stream(is);
    tag.check_end_name(type_name);
  }
  ARTS_METHOD_ERROR_CATCH
};

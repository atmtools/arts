#pragma once

#include <configtypes.h>
#include <mystring.h>
#include <supergeneric.h>

#include <complex>

#include "xml_io_stream.h"

template <>
struct xml_io_stream<Numeric> {
  constexpr static std::string_view type_name = "Numeric"sv;

  static void write(std::ostream& os,
                    const Numeric& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream& is, Numeric& n, bifstream* pbifs = nullptr);
  static void get(Numeric*, bifstream*, Size = 1);
  static void put(const Numeric* const, bofstream*, Size = 1);
};

template <>
struct xml_io_stream<Index> {
  constexpr static std::string_view type_name = "Index"sv;

  static void write(std::ostream& os,
                    const Index& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream& is, Index& n, bifstream* pbifs = nullptr);
  static void get(Index*, bifstream*, Size = 1);
  static void put(const Index* const, bofstream*, Size = 1);
};

template <>
struct xml_io_stream<Size> {
  constexpr static std::string_view type_name = "Size"sv;

  static void write(std::ostream& os,
                    const Size& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream& is, Size& n, bifstream* pbifs = nullptr);
  static void get(Size*, bifstream*, Size = 1);
  static void put(const Size* const, bofstream*, Size = 1);
};

template <>
struct xml_io_stream<std::complex<Numeric>> {
  constexpr static std::string_view type_name = "Complex"sv;

  static void write(std::ostream& os,
                    const std::complex<Numeric>& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is,
                   std::complex<Numeric>& n,
                   bifstream* pbifs = nullptr);
  static void get(std::complex<Numeric>*, bifstream*, Size = 1);
  static void put(const std::complex<Numeric>* const, bofstream*, Size = 1);
};

template <>
struct xml_io_stream<String> {
  constexpr static std::string_view type_name = "String"sv;

  static void write(std::ostream& os,
                    const String& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, String& n, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Any> {
  static constexpr std::string_view type_name = "Any"sv;

  static void write(std::ostream& os,
                    const Any& x,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream& is, Any& x, bifstream* pbifs = nullptr);
};

// No string or Any because they are not putable/getable
template <typename T>
concept xml_coretype =
    std::same_as<T, Numeric> or std::same_as<T, std::complex<Numeric>> or
    std::same_as<T, Index> or std::same_as<T, Size>;

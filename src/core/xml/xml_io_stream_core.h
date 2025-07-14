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
  static void get(std::span<Numeric>, bifstream*);
  static void put(std::span<const Numeric>, bofstream*);
  static void parse(std::span<Numeric>, std::istream&);
};

template <>
struct xml_io_stream<bool> {
  constexpr static std::string_view type_name = "Bool"sv;

  static void write(std::ostream& os,
                    const bool& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream& is, bool& n, bifstream* pbifs = nullptr);
};

template <>
struct xml_io_stream<Index> {
  constexpr static std::string_view type_name = "Index"sv;

  static void write(std::ostream& os,
                    const Index& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream& is, Index& n, bifstream* pbifs = nullptr);
  static void get(std::span<Index>, bifstream*);
  static void put(std::span<const Index>, bofstream*);
  static void parse(std::span<Index>, std::istream&);
};

template <>
struct xml_io_stream<Size> {
  constexpr static std::string_view type_name = "Size"sv;

  static void write(std::ostream& os,
                    const Size& n,
                    bofstream* pbofs      = nullptr,
                    std::string_view name = ""sv);
  static void read(std::istream& is, Size& n, bifstream* pbifs = nullptr);
  static void get(std::span<Size>, bifstream*);
  static void put(std::span<const Size>, bofstream*);
  static void parse(std::span<Size>, std::istream&);
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
  static void get(std::span<std::complex<Numeric>>, bifstream*);
  static void put(std::span<const std::complex<Numeric>>, bofstream*);
  static void parse(std::span<std::complex<Numeric>>, std::istream&);
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

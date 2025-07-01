#pragma once

#include <functional>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename... Ts>
struct xml_io_stream_name<std::function<Ts...>> {
  static constexpr std::string_view name = "Function"sv;
};

template <typename... Ts>
struct xml_io_stream<std::function<Ts...>> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<std::function<Ts...>>;

  static void write(std::ostream&,
                    const std::function<Ts...>&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv) {
    throw std::runtime_error("No XML IO for pure functional types");
  }

  static void read(std::istream&, std::function<Ts...>&, bifstream* = nullptr) {
    throw std::runtime_error("No XML IO for pure functional types");
  }
};

#pragma once

#include <functional>

#include "xml_io_base.h"
#include "xml_io_stream.h"

template <typename R, typename... Ts>
struct xml_io_stream_name<std::function<R(Ts...)>> {
  static constexpr std::string_view name = "Function"sv;
};

template <typename R, typename... Ts>
struct xml_io_stream<std::function<R(Ts...)>> {
  static constexpr std::string_view type_name =
      xml_io_stream_name_v<std::function<R(Ts...)>>;

  static void write(std::ostream&,
                    const std::function<R(Ts...)>&,
                    bofstream*       = nullptr,
                    std::string_view = ""sv) {
    throw std::runtime_error("No XML IO for pure functional types");
  }

  static void read(std::istream&,
                   std::function<R(Ts...)>&,
                   bifstream* = nullptr) try {
    throw std::runtime_error("No XML IO for pure functional types");
  } catch (const std::exception& e) {
    throw std::runtime_error(
        std::format("Error reading {}:\n{}", type_name, e.what()));
  }
};

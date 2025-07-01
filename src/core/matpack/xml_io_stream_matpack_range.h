#pragma once

#include <xml.h>

#include "matpack_mdspan_common.h"

template <>
struct xml_io_stream<Range> {
  static constexpr std::string_view type_name = "Range"sv;

  static void write(std::ostream &os,
                    const Range &x,
                    bofstream *pbofs      = nullptr,
                    std::string_view name = ""sv);

  static void read(std::istream &is, Range &x, bifstream *pbifs = nullptr);
};

#pragma once

#include <matpack.h>

struct Gravity {
  Numeric GM;
  Numeric a;
  Numeric e;

  Numeric operator()(Numeric h, Numeric lat, Numeric lon) const;

  static Gravity Mercury();
  static Gravity Venus();
  static Gravity Earth();
  static Gravity Moon();
  static Gravity Mars();
  static Gravity Jupiter();
  static Gravity Saturn();
};

template <>
struct xml_io_stream_name<Gravity> {
  static constexpr std::string_view name = "Gravity";
};

template <>
struct xml_io_stream_aggregate<Gravity> {
  static constexpr bool value = true;
};

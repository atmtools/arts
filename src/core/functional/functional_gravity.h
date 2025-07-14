#pragma once

#include <matpack.h>

struct EllipsoidGravity {
  Numeric GM;
  Numeric a;
  Numeric e;

  Numeric operator()(Numeric h, Numeric lat, Numeric lon) const;

  static EllipsoidGravity Mercury();
  static EllipsoidGravity Venus();
  static EllipsoidGravity Earth();
  static EllipsoidGravity Moon();
  static EllipsoidGravity Mars();
  static EllipsoidGravity Jupiter();
  static EllipsoidGravity Saturn();
};

template <>
struct xml_io_stream_name<EllipsoidGravity> {
  static constexpr std::string_view name = "EllipsoidGravity";
};

template <>
struct xml_io_stream_aggregate<EllipsoidGravity> {
  static constexpr bool value = true;
};

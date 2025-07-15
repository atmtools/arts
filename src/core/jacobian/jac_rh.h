#pragma once

#include <atm.h>
#include <matpack.h>

#include "jacobian_names.h"

struct rhfwd {
  NumericUnaryOperator psat;
  AtmKeyVal self_key;
  bool self_fix;

  Vector operator()(ConstVectorView x, const AtmField& y) const;
};

struct rhinv {
  NumericUnaryOperator psat;
  AtmKeyVal self_key;
  bool self_fix;

  Vector operator()(ConstVectorView x, const AtmField& y) const;

  Matrix operator()(ConstMatrixView dy,
                    ConstVectorView x,
                    const AtmField& y) const;
};

void make_rhfit(Jacobian::AtmTarget&,
                const AtmField&,
                const NumericUnaryOperator&,
                const bool);

template <>
struct xml_io_stream_name<rhfwd> {
  constexpr static std::string_view name = "rhfwd";
};

template <>
struct xml_io_stream_name<rhinv> {
  constexpr static std::string_view name = "rhinv";
};

template <>
struct xml_io_stream_aggregate<rhfwd> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_aggregate<rhinv> {
  constexpr static bool value = true;
};

#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurf.h>
#include <surf.h>

#include "jacobian_names.h"

struct logfwd {
  Size N;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const Jacobian::transformable_field auto&) const { return (*this)(x); }
};

struct loginv {
  Size N;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const Jacobian::transformable_field auto&) const { return (*this)(x); }

  Matrix operator()(ConstMatrixView dy, ConstVectorView x) const;
  Matrix operator()(ConstMatrixView dy, ConstVectorView x, const Jacobian::transformable_field auto&) const {
    return (*this)(dy, x);
  }
};

//! Makes the target the logarithm of its field
template <typename Target> void make_logfit(Target& x, const Jacobian::transformable_field auto& field) {
  const Size N       = field[x.type].flat_view().size();
  x.inverse_state    = loginv{.N = N};
  x.inverse_jacobian = loginv{.N = N};
  x.transform_state  = logfwd{.N = N};
}

template <> struct xml_io_stream_name<loginv> {
  constexpr static std::string_view name = "loginv";
};

template <> struct xml_io_stream_name<logfwd> {
  constexpr static std::string_view name = "logfwd";
};

template <> struct xml_io_stream_aggregate<loginv> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_aggregate<logfwd> {
  constexpr static bool value = true;
};

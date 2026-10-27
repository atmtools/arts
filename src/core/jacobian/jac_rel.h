#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurf.h>
#include <surf.h>

#include <memory>

#include "jacobian_names.h"

struct relfwd {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const Jacobian::transformable_field auto&) const { return (*this)(x); }
};

struct relinv {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const Jacobian::transformable_field auto&) const { return (*this)(x); }

  Matrix operator()(ConstMatrixView dy) const;
  Matrix operator()(ConstMatrixView dy, ConstVectorView, const Jacobian::transformable_field auto&) const {
    return (*this)(dy);
  }
};

//! Makes the target relative to the current values of its field
template <typename Target> void make_relfit(Target& x, const Jacobian::transformable_field auto& field) {
  const auto sorig   = std::make_shared<const Vector>(field[x.type].flat_view());
  x.inverse_state    = relinv{.sorig = sorig};
  x.inverse_jacobian = relinv{.sorig = sorig};
  x.transform_state  = relfwd{.sorig = sorig};
}

template <> struct xml_io_stream_name<relfwd> {
  constexpr static std::string_view name = "relfwd";
};

template <> struct xml_io_stream_name<relinv> {
  constexpr static std::string_view name = "relinv";
};

template <> struct xml_io_stream_aggregate<relfwd> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_aggregate<relinv> {
  constexpr static bool value = true;
};

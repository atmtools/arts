#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurf.h>
#include <surf.h>

#include <memory>

#include "jacobian_names.h"

struct logrelfwd {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const Jacobian::transformable_field auto&) const { return (*this)(x); }
};

struct logrelinv {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const Jacobian::transformable_field auto&) const { return (*this)(x); }

  Matrix operator()(ConstMatrixView dy, ConstVectorView x) const;
  Matrix operator()(ConstMatrixView dy, ConstVectorView x, const Jacobian::transformable_field auto&) const {
    return (*this)(dy, x);
  }
};

//! Makes the target the logarithm of its value relative to the current values of its field
template <typename Target> void make_logrelfit(Target& x, const Jacobian::transformable_field auto& field) {
  const auto sorig   = std::make_shared<const Vector>(field[x.type].flat_view());
  x.inverse_state    = logrelinv{.sorig = sorig};
  x.inverse_jacobian = logrelinv{.sorig = sorig};
  x.transform_state  = logrelfwd{.sorig = sorig};
}

template <> struct xml_io_stream_name<logrelfwd> {
  constexpr static std::string_view name = "logrelfwd";
};

template <> struct xml_io_stream_name<logrelinv> {
  constexpr static std::string_view name = "logrelinv";
};

template <> struct xml_io_stream_aggregate<logrelfwd> {
  constexpr static bool value = true;
};

template <> struct xml_io_stream_aggregate<logrelinv> {
  constexpr static bool value = true;
};

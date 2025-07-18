#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurface.h>
#include <surf.h>

#include <memory>

#include "jacobian_names.h"

struct logrelfwd {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const AtmField& y) const;
  Vector operator()(ConstVectorView x, const SurfaceField& y) const;
  Vector operator()(ConstVectorView x, const SubsurfaceField& y) const;
};

struct logrelinv {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const AtmField& y) const;
  Vector operator()(ConstVectorView x, const SurfaceField& y) const;
  Vector operator()(ConstVectorView x, const SubsurfaceField& y) const;

  Matrix operator()(ConstMatrixView dy, ConstVectorView x) const;
  Matrix operator()(ConstMatrixView dy,
                    ConstVectorView x,
                    const AtmField& y) const;
  Matrix operator()(ConstMatrixView dy,
                    ConstVectorView x,
                    const SurfaceField& y) const;
  Matrix operator()(ConstMatrixView dy,
                    ConstVectorView x,
                    const SubsurfaceField& y) const;
};

void make_logrelfit(Jacobian::AtmTarget&, const AtmField&);
void make_logrelfit(Jacobian::SurfaceTarget&, const SurfaceField&);
void make_logrelfit(Jacobian::SubsurfaceTarget&, const SubsurfaceField&);

template <>
struct xml_io_stream_name<logrelfwd> {
  constexpr static std::string_view name = "logrelfwd";
};

template <>
struct xml_io_stream_name<logrelinv> {
  constexpr static std::string_view name = "logrelinv";
};

template <>
struct xml_io_stream_aggregate<logrelfwd> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_aggregate<logrelinv> {
  constexpr static bool value = true;
};

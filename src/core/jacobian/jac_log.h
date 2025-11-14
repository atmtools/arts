#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurf.h>
#include <surf.h>

#include "jacobian_names.h"

struct logfwd {
  Size N;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const AtmField& y) const;
  Vector operator()(ConstVectorView x, const SurfaceField& y) const;
  Vector operator()(ConstVectorView x, const SubsurfaceField& y) const;
};

struct loginv {
  Size N;

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

void make_logfit(Jacobian::AtmTarget&, const AtmField&);
void make_logfit(Jacobian::SurfaceTarget&, const SurfaceField&);
void make_logfit(Jacobian::SubsurfaceTarget&, const SubsurfaceField&);

template <>
struct xml_io_stream_name<loginv> {
  constexpr static std::string_view name = "loginv";
};

template <>
struct xml_io_stream_name<logfwd> {
  constexpr static std::string_view name = "logfwd";
};

template <>
struct xml_io_stream_aggregate<loginv> {
  constexpr static bool value = true;
};

template <>
struct xml_io_stream_aggregate<logfwd> {
  constexpr static bool value = true;
};

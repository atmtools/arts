#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurface.h>
#include <surf.h>

#include <memory>

#include "jacobian.h"

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

#pragma once

#include <atm.h>
#include <matpack.h>
#include <subsurface.h>
#include <surf.h>

#include <memory>

#include "jacobian.h"

struct relfwd {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const AtmField& y) const;
  Vector operator()(ConstVectorView x, const SurfaceField& y) const;
  Vector operator()(ConstVectorView x, const SubsurfaceField& y) const;
};

struct relinv {
  std::shared_ptr<const Vector> sorig;

  Vector operator()(ConstVectorView x) const;
  Vector operator()(ConstVectorView x, const AtmField& y) const;
  Vector operator()(ConstVectorView x, const SurfaceField& y) const;
  Vector operator()(ConstVectorView x, const SubsurfaceField& y) const;

  Matrix operator()(ConstMatrixView dy) const;
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

void make_relfit(Jacobian::AtmTarget&, const AtmField&);
void make_relfit(Jacobian::SurfaceTarget&, const SurfaceField&);
void make_relfit(Jacobian::SubsurfaceTarget&, const SubsurfaceField&);

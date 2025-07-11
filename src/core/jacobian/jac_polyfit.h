#pragma once

#include <matpack.h>
#include <obsel.h>

#include <memory>

#include "jacobian.h"

void polyfit(VectorView param,
             const ConstVectorView x,
             const ConstVectorView y);

struct polyfit_t {
  std::shared_ptr<const Vector> st;
  Size polyorder;

  Vector operator()(ConstVectorView y) const;
  Vector operator()(ConstVectorView y, ConstVectorView) const;
};

struct polyfit_sensor_offset_t {
  polyfit_t fit;
  SensorKey key;

  Vector operator()(ConstVectorView y, const ArrayOfSensorObsel&) const;
};

struct polyinv_t {
  std::shared_ptr<const Vector> st;
  Size polyorder;

  Vector operator()(ConstVectorView x, ConstVectorView y) const;

  Matrix operator()(ConstMatrixView,
                    ConstVectorView x,
                    ConstVectorView y) const;
  Matrix operator()() const;
};

struct polyinv_sensor_offset_t {
  polyinv_t inv;
  SensorKey key;

  Vector operator()(ConstVectorView x, const ArrayOfSensorObsel& y) const;
};

void make_polyfit(Jacobian::ErrorTarget&, const Size, const Vector&);
void make_polyoffset(Jacobian::SensorTarget&,
                     const Size,
                     const ArrayOfSensorObsel&);

#include "jac_polyfit.h"

#include <minimize.h>

#include <utility>

#include "jacobian.h"

void polyfit(VectorView param,
             const ConstVectorView x,
             const ConstVectorView y) {
  ARTS_USER_ERROR_IF(param.size() < 1, "Must have atleast one fit-parameter.")

  using namespace Minimize;
  auto fit = polyfit(x, y, param.size() - 1);

  ARTS_USER_ERROR_IF(not fit,
                     R"(Could not fit polynomial:

  x:     {:B,},
  y:     {:B,},
  order: {}
)",
                     x,
                     y,
                     param.size() - 1)
  ARTS_USER_ERROR_IF(static_cast<Size>(fit->size()) != param.size(),
                     "Bad size. fit.size(): {}, param.size(): {}",
                     fit->size(),
                     param.size())

  param = *fit;
}

Vector polyfit_t::operator()(ConstVectorView y) const {
  ARTS_USER_ERROR_IF(not st, "No t-vector provided for polyinv.")
  const Vector& t = *st;

  ARTS_USER_ERROR_IF(y.size() != t.size(),
                     "Mismatched y and t sizes. y.size() vs t.size(): {} vs {}",
                     y.size(),
                     t.size())

  Vector p(polyorder + 1);

  polyfit(p, t, y);

  return p;
}

Vector polyfit_t::operator()(ConstVectorView y, ConstVectorView) const {
  return this->operator()(y);
}

Vector polyinv_t::operator()(ConstVectorView x, ConstVectorView y) const {
  ARTS_USER_ERROR_IF(not st, "No t-vector provided for polyfit.")
  const Vector& t = *st;

  ARTS_USER_ERROR_IF(y.size() != t.size(),
                     "Mismatched y and t sizes. y.size() vs t.size(): {} vs {}",
                     y.size(),
                     t.size())

  ARTS_USER_ERROR_IF(
      (polyorder + 1) != x.size(),
      "Mismatched polyorder and x sizes. polyorder + 1 vs x.size(): {} vs {}",
      polyorder + 1,
      x.size())

  Vector out{y};

  for (Size j = 0; j < t.size(); j++) {
    const Numeric tj = t[j];
    Numeric xn       = 1.0;
    for (Size i = 0; i <= polyorder; i++) {
      out[j] += x[i] * xn;
      xn     *= tj;
    }
  }

  return out;
}

Matrix polyinv_t::operator()() const {
  ARTS_USER_ERROR_IF(not st, "No t-vector provided for polyinv.")
  const Vector& t = *st;

  Matrix out(t.size(), polyorder + 1);

  for (Size j = 0; j < t.size(); j++) {
    const Numeric tj = t[j];
    Numeric xn       = 1.0;
    for (Size i = 0; i <= polyorder; i++) {
      out[j, i]  = xn;
      xn        *= tj;
    }
  }

  return out;
}

Matrix polyinv_t::operator()(ConstMatrixView,
                             ConstVectorView,
                             ConstVectorView) const {
  return this->operator()();
}

void make_polyfit(Jacobian::ErrorTarget& x,
                  const Size polyorder,
                  const Vector& t) {
  auto st = std::make_shared<const Vector>(t);
  const polyfit_t pfit{.st = st, .polyorder = polyorder};
  const polyinv_t pinv{.st = st, .polyorder = polyorder};
  x.x_size           = polyorder + 1;
  x.inverse_state    = pinv;
  x.inverse_jacobian = pinv;
  x.transform_state  = pfit;
}

namespace {
Vector rem_frq(const SensorObsel& v, const ConstVectorView x) {
  ARTS_USER_ERROR_IF(x.size() != v.f_grid().size(),
                     "Bad size. x.size(): {}, f_grid().size(): {}",
                     x.size(),
                     v.f_grid().size())

  const auto& fs = v.f_grid();

  Vector out(fs);
  out -= x;
  return out;
}

template <bool pos, Index k>
Vector rem_poslos(const SensorObsel& v, const ConstVectorView x) {
  ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                     "Bad size. x.size(): {}, poslos_grid().size(): {}",
                     x.size(),
                     v.poslos_grid().size())

  const SensorPosLosVector& xsv = v.poslos_grid();

  Vector out(xsv.size());
  std::transform(xsv.begin(),
                 xsv.end(),
                 x.begin(),
                 out.begin(),
                 [](auto poslos, Numeric val) {
                   if constexpr (pos) {
                     return poslos.pos[k] - val;
                   } else {
                     return poslos.los[k] - val;
                   }
                 });
  return out;
}

Vector rem_alt(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<true, 0>(v, x);
}

Vector rem_lat(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<true, 1>(v, x);
}

Vector rem_lon(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<true, 2>(v, x);
}

Vector rem_zen(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<false, 0>(v, x);
}

Vector rem_azi(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<false, 1>(v, x);
}
}  // namespace

Vector polyfit_sensor_offset_t::operator()(ConstVectorView,
                                           const ArrayOfSensorObsel& x) const {
  ARTS_USER_ERROR_IF(not fit.st, "No t-vector provided for polyinv.")
  const Vector& o      = *fit.st;
  const SensorObsel& v = x.at(key.measurement_elem);

  switch (key.type) {
    using enum SensorKeyType;
    case freq: return fit(rem_frq(v, o)); break;
    case zen:  return fit(rem_zen(v, o)); break;
    case azi:  return fit(rem_azi(v, o)); break;
    case alt:  return fit(rem_alt(v, o)); break;
    case lat:  return fit(rem_lat(v, o)); break;
    case lon:  return fit(rem_lon(v, o)); break;
  }

  std::unreachable();
}

Vector polyinv_sensor_offset_t::operator()(ConstVectorView x,
                                           const ArrayOfSensorObsel&) const {
  ARTS_USER_ERROR_IF(not inv.st, "No t-vector provided for polyinv.")
  const Vector& o = *inv.st;

  return inv(x, o);
}

void make_polyoffset(Jacobian::SensorTarget& x,
                     const Size polyorder,
                     const ArrayOfSensorObsel& sensor) {
  auto st = std::make_shared<const Vector>(
      sensor[x.type.measurement_elem].flat(x.type.type));
  const polyfit_sensor_offset_t pfit{
      .fit = polyfit_t{.st = st, .polyorder = polyorder}, .key = x.type};
  const polyinv_sensor_offset_t pinv{
      .inv = polyinv_t{.st = st, .polyorder = polyorder}, .key = x.type};
  x.x_size        = polyorder + 1;
  x.inverse_state = pinv;
  // x.inverse_jacobian = pinv;  Not neccessary, actively bad
  x.transform_state = pfit;
}

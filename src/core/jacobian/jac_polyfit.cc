#include "jac_polyfit.h"

#include <minimize.h>

#include <utility>

#include "jacobian.h"

void polyfit(VectorView param, const ConstVectorView x, const ConstVectorView y) {
  if (param.size() < 1) arts::user_error("Must have atleast one fit-parameter.");

  using namespace Minimize;
  auto fit = polyfit(x, y, param.size() - 1);

  if (not fit)
    arts::user_error(R"(Could not fit polynomial:

  x:     {:B,},
  y:     {:B,},
  order: {}
)",
                     x,
                     y,
                     param.size() - 1);
  if (static_cast<Size>(fit->size()) != param.size())
    arts::user_error("Bad size. fit.size(): {}, param.size(): {}", fit->size(), param.size());

  param = *fit;
}

Vector polyfit_t::operator()(ConstVectorView y) const {
  if (not st) arts::user_error("No t-vector provided for polyfit.");
  const Vector& t = *st;

  if (y.size() != t.size())
    arts::user_error("Mismatched y and t sizes. y.size() vs t.size(): {} vs {}", y.size(), t.size());

  Vector p(polyorder + 1);

  polyfit(p, t, y);

  return p;
}

Vector polyfit_t::operator()(ConstVectorView y, ConstVectorView) const { return this->operator()(y); }

Vector polyinv_t::operator()(ConstVectorView x, ConstVectorView y) const {
  if (not st) arts::user_error("No t-vector provided for polyinv.");
  const Vector& t = *st;

  if (y.size() != t.size())
    arts::user_error("Mismatched y and t sizes. y.size() vs t.size(): {} vs {}", y.size(), t.size());

  if ((polyorder + 1) != x.size())
    arts::user_error("Mismatched polyorder and x sizes. polyorder + 1 vs x.size(): {} vs {}", polyorder + 1, x.size());

  Vector out{y};

  for (Size j = 0; j < t.size(); j++) {
    const Numeric tj = t[j];
    Numeric       xn = 1.0;
    for (Size i = 0; i <= polyorder; i++) {
      out[j] += x[i] * xn;
      xn     *= tj;
    }
  }

  return out;
}

Matrix polyinv_t::operator()() const {
  if (not st) arts::user_error("No t-vector provided for polyinv.");
  const Vector& t = *st;

  Matrix out(t.size(), polyorder + 1);

  for (Size j = 0; j < t.size(); j++) {
    const Numeric tj = t[j];
    Numeric       xn = 1.0;
    for (Size i = 0; i <= polyorder; i++) {
      out[j, i]  = xn;
      xn        *= tj;
    }
  }

  return out;
}

Matrix polyinv_t::operator()(ConstMatrixView, ConstVectorView, ConstVectorView) const { return this->operator()(); }

void make_polyfit(Jacobian::ErrorTarget& x, const Size polyorder, const Vector& t) {
  auto            st = std::make_shared<const Vector>(t);
  const polyfit_t pfit{.st = st, .polyorder = polyorder};
  const polyinv_t pinv{.st = st, .polyorder = polyorder};
  x.x_size           = polyorder + 1;
  x.inverse_state    = pinv;
  x.inverse_jacobian = pinv;
  x.transform_state  = pfit;
}

namespace {
Vector rem_frq(const SensorObsel& v, const ConstVectorView x) {
  if (x.size() != v.f_grid().size())
    arts::user_error("Bad size. x.size(): {}, f_grid().size(): {}", x.size(), v.f_grid().size());

  const auto& fs = v.f_grid();

  Vector out(fs);
  out -= x;
  return out;
}

template <bool pos, Index k> Vector rem_poslos(const SensorObsel& v, const ConstVectorView x) {
  if (x.size() != v.poslos_grid().size())
    arts::user_error("Bad size. x.size(): {}, poslos_grid().size(): {}", x.size(), v.poslos_grid().size());

  const SensorPosLosVector& xsv = v.poslos_grid();

  Vector out(xsv.size());
  std::transform(xsv.begin(), xsv.end(), x.begin(), out.begin(), [](auto poslos, Numeric val) {
    if constexpr (pos) {
      return poslos.pos[k] - val;
    } else {
      return poslos.los[k] - val;
    }
  });
  return out;
}

Vector rem_alt(const SensorObsel& v, const ConstVectorView x) { return rem_poslos<true, 0>(v, x); }

Vector rem_lat(const SensorObsel& v, const ConstVectorView x) { return rem_poslos<true, 1>(v, x); }

Vector rem_lon(const SensorObsel& v, const ConstVectorView x) { return rem_poslos<true, 2>(v, x); }

Vector rem_zen(const SensorObsel& v, const ConstVectorView x) { return rem_poslos<false, 0>(v, x); }

Vector rem_azi(const SensorObsel& v, const ConstVectorView x) { return rem_poslos<false, 1>(v, x); }

//! The sensor grid mapped to [-1, 1], so the polynomial offset is well-conditioned and all its coefficients are shifts
std::shared_ptr<const Vector> normalized_grid(const std::shared_ptr<const Vector>& grid) {
  ARTS_USER_ERROR_IF(not grid or grid->size() == 0, "No sensor grid provided for the polynomial offset.")

  const auto [lo, hi] = stdr::minmax(*grid);
  const Numeric mid   = 0.5 * (hi + lo);
  const Numeric half  = hi > lo ? 0.5 * (hi - lo) : 1.0;

  auto t = std::make_shared<Vector>(*grid);
  for (auto& v : *t) v = (v - mid) / half;
  return t;
}
}  // namespace

polyfit_t polyfit_sensor_offset_t::normalized() const {
  return {.st = normalized_grid(fit.st), .polyorder = fit.polyorder};
}

polyinv_t polyinv_sensor_offset_t::normalized() const {
  return {.st = normalized_grid(inv.st), .polyorder = inv.polyorder};
}

Vector polyfit_sensor_offset_t::operator()(ConstVectorView, const ArrayOfSensorObsel& x) const {
  if (not fit.st) arts::user_error("No t-vector provided for polyfit.");
  const Vector&      o = *fit.st;
  const SensorObsel& v = x.at(key.measurement_elem);
  const polyfit_t    f = normalized();

  switch (key.type) {
    using enum SensorKeyType;
    case freq: return f(rem_frq(v, o)); break;
    case zen:  return f(rem_zen(v, o)); break;
    case azi:  return f(rem_azi(v, o)); break;
    case alt:  return f(rem_alt(v, o)); break;
    case lat:  return f(rem_lat(v, o)); break;
    case lon:  return f(rem_lon(v, o)); break;
  }

  std::unreachable();
}

Vector polyinv_sensor_offset_t::operator()(ConstVectorView x, const ArrayOfSensorObsel&) const {
  if (not inv.st) arts::user_error("No t-vector provided for polyinv.");
  const Vector& o = *inv.st;

  return normalized()(x, o);
}

void make_polyoffset(Jacobian::SensorTarget& x, const Size polyorder, const ArrayOfSensorObsel& sensor) {
  auto st = std::make_shared<const Vector>(sensor[x.type.measurement_elem].flat(x.type.type));

  ARTS_USER_ERROR_IF(polyorder >= st->size(),
                     "A polynomial offset of order {} needs more than {} sensor grid points",
                     polyorder,
                     st->size())

  const polyfit_sensor_offset_t pfit{.fit = polyfit_t{.st = st, .polyorder = polyorder}, .key = x.type};
  const polyinv_sensor_offset_t pinv{.inv = polyinv_t{.st = st, .polyorder = polyorder}, .key = x.type};
  x.x_size        = polyorder + 1;
  x.inverse_state = pinv;
  // No inverse_jacobian: spectral_rad_jacAddSensorJacobianPerturbations maps directly to the coefficients
  x.transform_state = pfit;
}

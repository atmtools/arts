#include "jacobian.h"

#include <compare.h>
#include <minimize.h>
#include <obsel.h>

#include <algorithm>
#include <utility>

namespace Jacobian {
void default_atm_x_set(VectorView x,
                       const AtmField& atm,
                       const AtmKeyVal& key) {
  ARTS_USER_ERROR_IF(
      not atm.contains(key), "Atmosphere does not contain key value {}", key)

  auto xn = atm[key].flat_view();

  ARTS_USER_ERROR_IF(
      atm[key].flat_view().size() not_eq x.size(),
      "Problem with sizes.  \n"
      "Did you change your atmosphere since you set the jacobian targets?\n"
      "Did you forget to finalize the JacobianTargets?")

  x = xn;
}

void default_x_atm_set(AtmField& atm,
                       const AtmKeyVal& key,
                       const ConstVectorView x) {
  ARTS_USER_ERROR_IF(
      not atm.contains(key), "Atmosphere does not contain key value {}", key)

  auto xn = atm[key].flat_view();

  ARTS_USER_ERROR_IF(
      atm[key].flat_view().size() not_eq x.size(),
      "Problem with sizes.  \n"
      "Did you change your atmosphere since you set the jacobian targets?\n"
      "Did you forget to finalize the JacobianTargets?")

  xn = x;
}

void default_surf_x_set(VectorView x,
                        const SurfaceField& surf,
                        const SurfaceKeyVal& key) {
  ARTS_USER_ERROR_IF(
      not surf.contains(key), "Surface does not contain key value {}", key)

  auto xn = surf[key].flat_view();

  ARTS_USER_ERROR_IF(
      xn.size() not_eq x.size(),
      "Problem with sizes.\n"
      "Did you change your surface since you set the jacobian targets?\n"
      "Did you forget to finalize the JacobianTargets?")

  x = xn;
}

void default_x_surf_set(SurfaceField& surf,
                        const SurfaceKeyVal& key,
                        const ConstVectorView x) {
  ARTS_USER_ERROR_IF(
      not surf.contains(key), "Surface does not contain key value {}", key)

  auto xn = surf[key].flat_view();

  ARTS_USER_ERROR_IF(
      xn.size() not_eq x.size(),
      "Problem with sizes.\n"
      "Did you change your surface since you set the jacobian targets?\n"
      "Did you forget to finalize the JacobianTargets?")

  xn = x;
}

void default_line_x_set(VectorView x,
                        const AbsorptionBands& bands,
                        const LblLineKey& key) {
  x = key.get_value(bands);
}

void default_x_line_set(AbsorptionBands& bands,
                        const LblLineKey& key,
                        const ConstVectorView x) {
  VectorView{key.get_value(bands)} = x;
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

Vector rem_zag(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<false, 0>(v, x);
}

Vector rem_aag(const SensorObsel& v, const ConstVectorView x) {
  return rem_poslos<false, 1>(v, x);
}
}  // namespace

void polyfit(VectorView param,
             const ConstVectorView x,
             const ConstVectorView y) {
  ARTS_USER_ERROR_IF(param.size() < 1, "Must have atleast one fit-parameter.")

  using namespace Minimize;
  auto [success, fit] = curve_fit<Polynom>(x, y, param.size() - 1);

  ARTS_USER_ERROR_IF(not success,
                     R"(Could not fit polynomial:
  x:   {:B,},
  y:   {:B,}
  fit: {:B,}
)",
                     x,
                     y,
                     Vector{fit})
  ARTS_USER_ERROR_IF(static_cast<Size>(fit.size()) != param.size(),
                     "Bad size. fit.size(): {}, param.size(): {}",
                     fit.size(),
                     param.size())

  param = fit;
}

void default_sensor_x_set(VectorView x,
                          const ArrayOfSensorObsel& sensor,
                          const SensorKey& key) {
  const SensorObsel& v = sensor.at(key.measurement_elem);
  const Vector& o      = key.original_grid;

  using enum SensorKeyType;
  switch (key.model) {
    case SensorJacobianModelType::None: v.flat(x, key.type); break;
    case SensorJacobianModelType::PolynomialOffset:
      switch (key.type) {
        case f:   polyfit(x, o, rem_frq(v, o)); break;
        case za:  polyfit(x, o, rem_zag(v, o)); break;
        case aa:  polyfit(x, o, rem_aag(v, o)); break;
        case alt: polyfit(x, o, rem_alt(v, o)); break;
        case lat: polyfit(x, o, rem_lat(v, o)); break;
        case lon: polyfit(x, o, rem_lon(v, o)); break;
      }
      break;
  }
}

namespace {
// Returns p + x[0] + x[1]*p + x[2]*p^2 + ...
Vector polynomial_offset_evaluate(const ConstVectorView x, const Vector& p) {
  ARTS_USER_ERROR_IF(x.empty(), "Must have some polynomial coefficients.")

  Vector out(p);

  for (Size i = 0; i < p.size(); ++i) {
    out[i]    += x[0];
    Numeric v  = 1.0;
    for (Size j = 1; j < x.size(); ++j) {
      out[i] += x[j] * (v *= p[i]);
    }
  }

  return out;
}
}  // namespace

void default_x_sensor_set(ArrayOfSensorObsel& sensor,
                          const SensorKey& key,
                          const ConstVectorView x) {
  auto& v = sensor.at(key.measurement_elem);

  using enum SensorKeyType;
  switch (key.model) {
    case SensorJacobianModelType::None:
      unflatten(sensor, x, v, key.type);
      break;
    case SensorJacobianModelType::PolynomialOffset: {
      auto& o        = key.original_grid;
      const Vector r = polynomial_offset_evaluate(x, o);
      unflatten(sensor, r, v, key.type);
    } break;
  }
}

void AtmTarget::update(AtmField& atm, const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(atm, type, x[Range(x_start, x_size)]);
}

void AtmTarget::update(Vector& x, const AtmField& atm) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x[Range(x_start, x_size)], atm, type);
}

bool AtmTarget::is_wind() const {
  return type == AtmKey::wind_u or type == AtmKey::wind_v or
         type == AtmKey::wind_w;
}

void SensorTarget::update(ArrayOfSensorObsel& sens, const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(sens, type, x[Range(x_start, x_size)]);
}

void SensorTarget::update(Vector& x, const ArrayOfSensorObsel& sens) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x[Range(x_start, x_size)], sens, type);
}

void SurfaceTarget::update(SurfaceField& surf, const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(surf, type, x[Range(x_start, x_size)]);
}

void SurfaceTarget::update(Vector& x, const SurfaceField& surf) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x[Range(x_start, x_size)], surf, type);
}

void LineTarget::update(AbsorptionBands& absorption_bands,
                        const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(absorption_bands, type, x[Range(x_start, x_size)]);
}

void LineTarget::update(Vector& x,
                        const AbsorptionBands& absorption_bands) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x[Range(x_start, x_size)], absorption_bands, type);
}

void ErrorTarget::update_y(Vector& y, Matrix& dy, const Vector& x) const {
  const Index szx = x.size();
  ARTS_USER_ERROR_IF(szx < static_cast<Index>(x_start + x_size),
                     "Got too small x.")

  const Index szy = y.size();
  ARTS_USER_ERROR_IF(szy < static_cast<Index>(type.y_start + type.y_size),
                     "Got too small y.")

  ARTS_USER_ERROR_IF((dy.shape() != std::array{szy, szx}),
                     "Mismatched dy shape. {:B,} vs [{}, {}]",
                     dy.shape(),
                     szy,
                     szx)

  set_y(y[Range(type.y_start, type.y_size)],
        dy[Range(type.y_start, type.y_size), Range(x_start, x_size)],
        x[Range(x_start, x_size)]);
}

void ErrorTarget::update_x(Vector& x, const Vector& y, const Vector& y0) const {
  ARTS_USER_ERROR_IF(y.size() != y0.size(), "Got different size y and y0.")

  const auto szx = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(szx < (x_start + x_size), "Got too small x.")

  const auto szy = static_cast<Size>(y.size());
  ARTS_USER_ERROR_IF(szy < (type.y_start + type.y_size), "Got too small y.")

  Vector e(y[Range(type.y_start, type.y_size)]);
  e -= y0[Range(type.y_start, type.y_size)];

  set_x(x[Range(x_start, x_size)], e);
}

const std::vector<AtmTarget>& Targets::atm() const {
  return target<AtmTarget>();
}

const std::vector<SurfaceTarget>& Targets::surf() const {
  return target<SurfaceTarget>();
}

const std::vector<LineTarget>& Targets::line() const {
  return target<LineTarget>();
}

const std::vector<SensorTarget>& Targets::sensor() const {
  return target<SensorTarget>();
}

const std::vector<ErrorTarget>& Targets::error() const {
  return target<ErrorTarget>();
}

std::vector<AtmTarget>& Targets::atm() { return target<AtmTarget>(); }

std::vector<SurfaceTarget>& Targets::surf() { return target<SurfaceTarget>(); }

std::vector<LineTarget>& Targets::line() { return target<LineTarget>(); }

std::vector<SensorTarget>& Targets::sensor() { return target<SensorTarget>(); }

std::vector<ErrorTarget>& Targets::error() { return target<ErrorTarget>(); }

void Targets::finalize(const AtmField& atmospheric_field,
                       const SurfaceField& surface_field,
                       const AbsorptionBands&,
                       const ArrayOfSensorObsel& measurement_sensor) {
  const Size natm    = atm().size();
  const Size nsurf   = surf().size();
  const Size nline   = line().size();
  const Size nsensor = sensor().size();

  Size last_size = 0;

  for (Size i = 0; i < natm; i++) {
    AtmTarget& t = atm()[i];
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(
            atm() | std::views::drop(i + 1), Cmp::eq(t.type), &AtmTarget::type),
        "Multiple targets of the same type: {}",
        t.type)
    t.x_start  = last_size;
    t.x_size   = atmospheric_field[t.type].flat_view().size();
    last_size += t.x_size;
  }

  for (Size i = 0; i < nsurf; i++) {
    SurfaceTarget& t = surf()[i];
    ARTS_USER_ERROR_IF(std::ranges::any_of(surf() | std::views::drop(i + 1),
                                           Cmp::eq(t.type),
                                           &SurfaceTarget::type),
                       "Multiple targets of the same type: {}",
                       t.type)
    t.x_start  = last_size;
    t.x_size   = surface_field[t.type].flat_view().size();
    last_size += t.x_size;
  }

  for (Size i = 0; i < nline; i++) {
    LineTarget& t = line()[i];
    ARTS_USER_ERROR_IF(std::ranges::any_of(line() | std::views::drop(i + 1),
                                           Cmp::eq(t.type),
                                           &LineTarget::type),
                       "Multiple targets of the same type: {}",
                       t.type)
    t.x_start  = last_size;
    t.x_size   = 1;
    last_size += t.x_size;
  }

  for (auto& elem : sensor()) {
    ARTS_USER_ERROR_IF(static_cast<Size>(elem.type.measurement_elem) >=
                           measurement_sensor.size(),
                       "Bad sensor elements {}, out-of-bounds",
                       elem);
  }

  for (Size i = 0; i < nsensor; i++) {
    SensorTarget& t = sensor()[i];
    ARTS_USER_ERROR_IF(
        std::ranges::any_of(
            sensor() | std::views::drop(i + 1),
            [&](const SensorKey& key) {
              if (t.type.type != key.type) return false;
              if (t.type.model != key.model) return false;

              const Index elem1 = t.type.measurement_elem;
              const Index elem2 = key.measurement_elem;

              switch (key.type) {
                using enum SensorKeyType;
                case f:
                  return measurement_sensor[elem1].f_grid_ptr() ==
                         measurement_sensor[elem2].f_grid_ptr();
                case za:
                case aa:
                case alt:
                case lat:
                case lon:
                  return measurement_sensor[elem1].poslos_grid_ptr() ==
                         measurement_sensor[elem2].poslos_grid_ptr();
              }

              std::unreachable();
            },
            &SensorTarget::type),
        "Multiple targets of the same type for target: {} "
        "- note that different sensor element targets may not share the same type, model, and relevant grids",
        t.type)

    t.x_start = last_size;
    switch (t.type.model) {
      case SensorJacobianModelType::None:
        t.x_size = measurement_sensor.at(t.type.measurement_elem)
                       .flat_size(t.type.type);
        break;
      case SensorJacobianModelType::PolynomialOffset:
        ARTS_USER_ERROR_IF(t.type.polyorder < 0,
                           "Must have a polynomial order.")
        t.x_size = t.type.polyorder + 1;
        break;
    }
    last_size += t.x_size;
  }

  for (auto& t : error()) {
    t.x_start = last_size;
    // t.x_size already known
    last_size += t.x_size;

    ARTS_USER_ERROR_IF(
        t.type.y_start + t.type.y_size > measurement_sensor.size(),
        R"(An error target is out of bounds of the y-vector.

The error target is defined to operator on the measurement_vector in the range {},
but the measurement_vector will only have {} element by the measurement_sensor.
)",
        t.type,
        measurement_sensor.size())
  }

  finalized = true;

  throwing_check(last_size);
}

AtmTarget& Targets::emplace_back(AtmKeyVal&& t, Numeric d) {
  return atm().emplace_back(std::move(t), d, target_count());
}

SurfaceTarget& Targets::emplace_back(SurfaceKeyVal&& t, Numeric d) {
  return surf().emplace_back(std::move(t), d, target_count());
}

LineTarget& Targets::emplace_back(LblLineKey&& t, Numeric d) {
  return line().emplace_back(std::move(t), d, target_count());
}

SensorTarget& Targets::emplace_back(SensorKey&& t, Numeric d) {
  return sensor().emplace_back(std::move(t), d, target_count());
}

AtmTarget& Targets::emplace_back(const AtmKeyVal& t, Numeric d) {
  return atm().emplace_back(t, d, target_count());
}

SurfaceTarget& Targets::emplace_back(const SurfaceKeyVal& t, Numeric d) {
  return surf().emplace_back(t, d, target_count());
}

LineTarget& Targets::emplace_back(const LblLineKey& t, Numeric d) {
  return line().emplace_back(t, d, target_count());
}

SensorTarget& Targets::emplace_back(const SensorKey& t, Numeric d) {
  return sensor().emplace_back(t, d, target_count());
}

ErrorTarget& Targets::emplace_back(const ErrorKey& target,
                                   Size x_size,
                                   ErrorTargetSetState&& set_y,
                                   ErrorTargetSetModel&& set_x) {
  return error().emplace_back(ErrorTarget{.type       = target,
                                          .target_pos = target_count(),
                                          .x_size     = x_size,
                                          .set_y      = std::move(set_y),
                                          .set_x      = std::move(set_x)});
}
}  // namespace Jacobian

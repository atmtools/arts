#include "jacobian.h"

#include <algorithm>

#include "compare.h"
#include "minimize.h"
#include "obsel.h"

namespace Jacobian {
void default_atm_x_set(ExhaustiveVectorView x,
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
                       const ExhaustiveConstVectorView x) {
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

void default_surf_x_set(ExhaustiveVectorView x,
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
                        const ExhaustiveConstVectorView x) {
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

void default_line_x_set(ExhaustiveVectorView x,
                        const AbsorptionBands& bands,
                        const LblLineKey& key) {
  x = key.get_value(bands);
}

void default_x_line_set(AbsorptionBands& bands,
                        const LblLineKey& key,
                        const ExhaustiveConstVectorView x) {
  ExhaustiveVectorView{key.get_value(bands)} = x;
}

Vector rem_frq(const SensorObsel& v, const ExhaustiveConstVectorView x) {
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
Vector rem_poslos(const SensorObsel& v, const ExhaustiveConstVectorView x) {
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

Vector rem_alt(const SensorObsel& v, const ExhaustiveConstVectorView x) {
  return rem_poslos<true, 0>(v, x);
}

Vector rem_lat(const SensorObsel& v, const ExhaustiveConstVectorView x) {
  return rem_poslos<true, 1>(v, x);
}

Vector rem_lon(const SensorObsel& v, const ExhaustiveConstVectorView x) {
  return rem_poslos<true, 2>(v, x);
}

Vector rem_zag(const SensorObsel& v, const ExhaustiveConstVectorView x) {
  return rem_poslos<false, 0>(v, x);
}

Vector rem_aag(const SensorObsel& v, const ExhaustiveConstVectorView x) {
  return rem_poslos<false, 1>(v, x);
}

void polyfit(ExhaustiveVectorView param, const Vector& x, const Vector& y) {
  ARTS_USER_ERROR_IF(param.size() < 1, "Must have atleast one fit-parameter.")

  using namespace Minimize;
  auto [success, fit] = curve_fit<Polynom>(x, y, param.size() - 1);

  ARTS_USER_ERROR_IF(not success, "Could not fit polynomial.")
  ARTS_USER_ERROR_IF(fit.size() != param.size(),
                     "Bad size. fit.size(): {}, param.size(): {}",
                     fit.size(),
                     param.size())

  param = fit;
}

void default_sensor_x_set(ExhaustiveVectorView x,
                          const ArrayOfSensorObsel& sensor,
                          const SensorKey& key) {
  const SensorObsel& v = sensor.at(key.elem);
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

// Returns p + x[0] + x[1]*p + x[2]*p^2 + ...
Vector polynomial_offset_evaluate(const ExhaustiveConstVectorView x,
                                  const Vector& p) {
  ARTS_USER_ERROR_IF(x.empty(), "Must have some polynomial coefficients.")

  Vector out(p);

  for (Index i = 0; i < p.size(); ++i) {
    out[i]    += x[0];
    Numeric v  = 1.0;
    for (Index j = 1; j < x.size(); ++j) {
      out[i] += x[j] * (v *= p[i]);
    }
  }

  return out;
}

void default_x_sensor_set(ArrayOfSensorObsel& sensor,
                          const SensorKey& key,
                          const ExhaustiveConstVectorView x) {
  auto& v = sensor.at(key.elem);

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
  set_model(atm, type, x.slice(x_start, x_size));
}

void AtmTarget::update(Vector& x, const AtmField& atm) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x.slice(x_start, x_size), atm, type);
}

bool AtmTarget::is_wind() const {
  return type == AtmKey::wind_u or type == AtmKey::wind_v or
         type == AtmKey::wind_w;
}

void SurfaceTarget::update(SurfaceField& surf, const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(surf, type, x.slice(x_start, x_size));
}

void SurfaceTarget::update(Vector& x, const SurfaceField& surf) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x.slice(x_start, x_size), surf, type);
}

void LineTarget::update(AbsorptionBands& absorption_bands,
                        const Vector& x) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_model(absorption_bands, type, x.slice(x_start, x_size));
}

void LineTarget::update(Vector& x,
                        const AbsorptionBands& absorption_bands) const {
  const auto sz = static_cast<Size>(x.size());
  ARTS_USER_ERROR_IF(sz < (x_start + x_size), "Got too small vector.")
  set_state(x.slice(x_start, x_size), absorption_bands, type);
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

std::vector<AtmTarget>& Targets::atm() { return target<AtmTarget>(); }

std::vector<SurfaceTarget>& Targets::surf() { return target<SurfaceTarget>(); }

std::vector<LineTarget>& Targets::line() { return target<LineTarget>(); }

std::vector<SensorTarget>& Targets::sensor() { return target<SensorTarget>(); }

void Targets::finalize(const AtmField& atmospheric_field,
                       const SurfaceField& surface_field,
                       const AbsorptionBands&,
                       const ArrayOfSensorObsel& measurement_sensor) {
  zero_out_x();

  const Size natm    = atm().size();
  const Size nsurf   = surf().size();
  const Size nline   = line().size();
  const Size nsensor = sensor().size();
  ARTS_ASSERT(target_count() == (natm + nsurf + nline + nsensor));

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

  for (Size i = 0; i < nsensor; i++) {
    SensorTarget& t = sensor()[i];
    ARTS_USER_ERROR_IF(std::ranges::any_of(sensor() | std::views::drop(i + 1),
                                           Cmp::eq(t.type),
                                           &SensorTarget::type),
                       "Multiple targets of the same type: {}",
                       t.type)
    t.x_start = last_size;
    switch (t.type.model) {
      case SensorJacobianModelType::None:
        t.x_size = measurement_sensor.at(t.type.elem).flat_size(t.type.type);
        break;
      case SensorJacobianModelType::PolynomialOffset:
        ARTS_USER_ERROR_IF(t.type.polyorder < 0,
                           "Must have a polynomial order.")
        t.x_size = t.type.polyorder + 1;
        break;
    }
    last_size += t.x_size;
  }

  finalized = true;

  throwing_check(last_size);
}
}  // namespace Jacobian

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

  using enum SensorKeyType;
  switch (key.model) {
    case SensorJacobianModelType::None: {
      switch (key.type) {
        case Frequency:
          ARTS_USER_ERROR_IF(x.size() != v.f_grid().size(),
                             "Bad size. x.size(): {}, f_grid().size(): {}",
                             x.size(),
                             v.f_grid().size())
          x = v.f_grid();
          break;
        case PointingZenith:
          ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                             "Bad size. x.size(): {}, poslos_grid().size(): {}",
                             x.size(),
                             v.poslos_grid().size())
          std::transform(v.poslos_grid().begin(),
                         v.poslos_grid().end(),
                         x.begin(),
                         [](auto& poslos) { return poslos.los[0]; });
          break;
        case PointingAzimuth:
          ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                             "Bad size. x.size(): {}, poslos_grid().size(): {}",
                             x.size(),
                             v.poslos_grid().size())
          std::transform(v.poslos_grid().begin(),
                         v.poslos_grid().end(),
                         x.begin(),
                         [](auto& poslos) { return poslos.los[1]; });
          break;
        case PointingAltitude:
          ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                             "Bad size. x.size(): {}, poslos_grid().size(): {}",
                             x.size(),
                             v.poslos_grid().size())
          std::transform(v.poslos_grid().begin(),
                         v.poslos_grid().end(),
                         x.begin(),
                         [](auto& poslos) { return poslos.pos[0]; });
          break;
        case PointingLatitude:
          ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                             "Bad size. x.size(): {}, poslos_grid().size(): {}",
                             x.size(),
                             v.poslos_grid().size())
          std::transform(v.poslos_grid().begin(),
                         v.poslos_grid().end(),
                         x.begin(),
                         [](auto& poslos) { return poslos.pos[1]; });
          break;
        case PointingLongitude:
          ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                             "Bad size. x.size(): {}, poslos_grid().size(): {}",
                             x.size(),
                             v.poslos_grid().size())
          std::transform(v.poslos_grid().begin(),
                         v.poslos_grid().end(),
                         x.begin(),
                         [](auto& poslos) { return poslos.pos[2]; });
          break;
        case Weights: {
          auto& w = v.weight_matrix();
          ARTS_USER_ERROR_IF(w.size() * 4 != x.size(),
                             "Bad size. x.size(): {} for w.size()*4: {}",
                             x.size(),
                             w.size() * 4)
          auto X = x.reshape_as(w.nrows(), w.ncols(), 4);
          for (Index i = 0; i < X.npages(); ++i) {
            for (Index j = 0; j < X.nrows(); ++j) {
              X(i, j, 0) = w(i, j).I();
              X(i, j, 1) = w(i, j).Q();
              X(i, j, 2) = w(i, j).U();
              X(i, j, 3) = w(i, j).V();
            }
          }
        } break;
      }
    } break;
    case SensorJacobianModelType::PolynomialOffset: {
      const Vector& o = key.original_grid;
      switch (key.type) {
        case Frequency:         polyfit(x, o, rem_frq(v, o)); break;
        case PointingZenith:    polyfit(x, o, rem_zag(v, o)); break;
        case PointingAzimuth:   polyfit(x, o, rem_aag(v, o)); break;
        case PointingAltitude:  polyfit(x, o, rem_alt(v, o)); break;
        case PointingLatitude:  polyfit(x, o, rem_lat(v, o)); break;
        case PointingLongitude: polyfit(x, o, rem_lon(v, o)); break;
        case Weights:           ARTS_USER_ERROR("Not available for polynomial model.");
      }
    } break;
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

void set_frq(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ExhaustiveConstVectorView x) {
  ARTS_USER_ERROR_IF(x.size() != v.f_grid().size(),
                     "Bad size. x.size(): {}, f_grid().size(): {}",
                     x.size(),
                     v.f_grid().size())

  const auto xs = std::make_shared<const AscendingGrid>(
      x.begin(), x.end(), [](auto& x) { return x; });

  // Must copy, as we may change the shared_ptr later
  const auto fs = v.f_grid_ptr();

  for (auto& elem : sensor) {
    if (elem.f_grid_ptr() == fs) {
      elem.set_f_grid_ptr(xs);  // may change here
    }
  }
}

template <bool pos, Index k>
void set_poslos(const SensorObsel& v,
                ArrayOfSensorObsel& sensor,
                const ExhaustiveConstVectorView x) {
  ARTS_USER_ERROR_IF(x.size() != v.poslos_grid().size(),
                     "Bad size. x.size(): {}, poslos_grid().size(): {}",
                     x.size(),
                     v.poslos_grid().size())

  SensorPosLosVector xsv = v.poslos_grid();

  std::transform(xsv.begin(),
                 xsv.end(),
                 x.begin(),
                 xsv.begin(),
                 [](auto poslos, Numeric val) {
                   if constexpr (pos) {
                     poslos.pos[k] = val;
                   } else {
                     poslos.los[k] = val;
                   }
                   return poslos;
                 });

  const auto xs = std::make_shared<const SensorPosLosVector>(std::move(xsv));

  // Must copy, as we may change the shared_ptr later
  const auto ps = v.poslos_grid_ptr();

  for (auto& elem : sensor) {
    if (elem.poslos_grid_ptr() == ps) {
      elem.set_poslos_grid_ptr(ps);
    }
  }
}

void set_alt(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ExhaustiveConstVectorView x) {
  set_poslos<true, 0>(v, sensor, x);
}

void set_lat(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ExhaustiveConstVectorView x) {
  set_poslos<true, 1>(v, sensor, x);
}

void set_lon(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ExhaustiveConstVectorView x) {
  set_poslos<true, 2>(v, sensor, x);
}

void set_zag(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ExhaustiveConstVectorView x) {
  set_poslos<false, 0>(v, sensor, x);
}

void set_aag(const SensorObsel& v,
             ArrayOfSensorObsel& sensor,
             const ExhaustiveConstVectorView x) {
  set_poslos<false, 1>(v, sensor, x);
}

void default_x_sensor_set(ArrayOfSensorObsel& sensor,
                          const SensorKey& key,
                          const ExhaustiveConstVectorView x) {
  auto& v = sensor.at(key.elem);

  using enum SensorKeyType;
  switch (key.model) {
    case SensorJacobianModelType::None: {
      switch (key.type) {
        case Frequency:         set_frq(v, sensor, x); break;
        case PointingZenith:    set_zag(v, sensor, x); break;
        case PointingAzimuth:   set_aag(v, sensor, x); break;
        case PointingAltitude:  set_alt(v, sensor, x); break;
        case PointingLatitude:  set_lat(v, sensor, x); break;
        case PointingLongitude: set_lon(v, sensor, x); break;
        case Weights:           {
          auto [r, c] = v.weight_matrix().shape();

          ARTS_USER_ERROR_IF(r * c * 4 != x.size(),
                             "Bad size. x.size(): {} for w.size()*4: {}",
                             x.size(),
                             r * c * 4)

          auto X = x.reshape_as(r, c, 4);
          StokvecMatrix ws{r, c};
          for (Index i = 0; i < r; ++i) {
            for (Index j = 0; j < c; ++j) {
              ws(i, j).I() = X(i, j, 0);
              ws(i, j).Q() = X(i, j, 1);
              ws(i, j).U() = X(i, j, 2);
              ws(i, j).V() = X(i, j, 3);
            }
          }

          v.set_weight_matrix(std::move(ws));
        } break;
      }
    } break;
    case SensorJacobianModelType::PolynomialOffset: {
      auto& o  = key.original_grid;
      Vector r = polynomial_offset_evaluate(x, o);

      switch (key.type) {
        case Frequency:         set_frq(v, sensor, r); break;
        case PointingZenith:    set_zag(v, sensor, r); break;
        case PointingAzimuth:   set_aag(v, sensor, r); break;
        case PointingAltitude:  set_alt(v, sensor, r); break;
        case PointingLatitude:  set_lat(v, sensor, r); break;
        case PointingLongitude: set_lon(v, sensor, r); break;
        case Weights:           ARTS_USER_ERROR("Not available for polynomial model.");
      }
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

std::vector<AtmTarget>& Targets::atm() { return target<AtmTarget>(); }

std::vector<SurfaceTarget>& Targets::surf() { return target<SurfaceTarget>(); }

std::vector<LineTarget>& Targets::line() { return target<LineTarget>(); }

void Targets::finalize(const AtmField& atmospheric_field,
                       const SurfaceField& surface_field,
                       const AbsorptionBands&,
                       const ArrayOfSensorObsel& measurement_sensor) {
  zero_out_x();

  const Size natm  = atm().size();
  const Size nsurf = surf().size();
  const Size nline = line().size();
  ARTS_ASSERT(target_count() == (natm + nsurf + nline));

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

  finalized = true;

  throwing_check(last_size);
}
}  // namespace Jacobian

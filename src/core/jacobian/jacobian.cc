#include "jacobian.h"

#include <compare.h>
#include <minimize.h>
#include <obsel.h>

#include <algorithm>
#include <utility>

namespace Jacobian {
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

////////////////////////////////////////////////////////////////////////
/// Templates for doing the common work of updating fields, model state vectors, and Jacobians
////////////////////////////////////////////////////////////////////////

namespace {
template <typename Func, typename Key, typename Field>
void update_x(VectorView x_state,
              const ConstVectorView x_field,
              Func&& transform_state,
              const Key& type,
              const Field& field) {
  if (transform_state) {
    const Vector xn_transformed = transform_state(Vector{x_field}, field);

    ARTS_USER_ERROR_IF(xn_transformed.size() not_eq x_state.size(),
                       R"(Size mismatch in transformation for target {}.

Cannot set the model state vector from the transformation of the field.

The size of the target field is       : {}.
The size of the transformed target is : {}.
The expected size of the target is    : {}.
)",
                       type,
                       x_field.size(),
                       xn_transformed.size(),
                       x_state.size())

    x_state = xn_transformed;
  } else {
    ARTS_USER_ERROR_IF(x_field.size() not_eq x_state.size(),
                       R"(Size mismatch in Jacobian target for {}.

Cannot set the model state vector from the field.

The size of the target field is    : {}.
The expected size of the target is : {}.
)",
                       type,
                       x_field.size(),
                       x_state.size())

    x_state = x_field;
  }
}
template <typename Func, typename Key, typename Field>
void update_dy(StridedMatrixView dy,
               Func&& transform_jacobian,
               const Key& type,
               const Field& field) {
  const Matrix dy_transformed = transform_jacobian(Matrix{dy}, field);

  ARTS_USER_ERROR_IF(dy_transformed.shape() != dy.shape(),
                     R"(Size mismatch in Jacobian transformation for target {}.

Cannot transform the Jacobian matrix.

The original shape is                 : {:B,}.
The size of the transformed target is : {:B,}.
)",
                     type,
                     dy.shape(),
                     dy_transformed.shape())

  dy = dy_transformed;
}

template <typename Func, typename Key, typename Field>
void update_field(VectorView x_field,
                  const ConstVectorView x_state,
                  Func&& inverse_state,
                  const Key& type,
                  const Field& field) {
  if (inverse_state) {
    const Vector x_inverse = inverse_state(Vector{x_state}, field);

    ARTS_USER_ERROR_IF(x_inverse.size() not_eq x_field.size(),
                       R"(Size mismatch in inverse transformation for target {}.

Cannot set the field from the inverse transformation of the model state vector.

The original size of the target is            : {}.
The size of the inverse transformed target is : {}.
The expected size of the target is            : {}.
)",
                       type,
                       x_field.size(),
                       x_inverse.size(),
                       x_state.size())

    x_field = x_inverse;
  } else {
    ARTS_USER_ERROR_IF(x_field.size() not_eq x_state.size(),
                       R"(Size mismatch for target {}.

Cannot set the field from the model state vector.

The size of the target is          : {}.
The expected size of the target is : {}.
)",
                       type,
                       x_field.size(),
                       x_state.size())

    x_field = x_state;
  }
}
}  // namespace

////////////////////////////////////////////////////////////////////////
/// Update the fields
////////////////////////////////////////////////////////////////////////

void ErrorTarget::update_model(VectorView meas, const ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  update_field(meas, x[Range(x_start, x_size)], inverse_state, type, meas);
}

void AtmTarget::update_model(AtmField& atm, const ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ARTS_USER_ERROR_IF(
      not atm.contains(type), "Atmosphere does not contain key value {}", type)

  update_field(atm[type].flat_view(),
               x[Range(x_start, x_size)],
               inverse_state,
               type,
               atm);
}

void SurfaceTarget::update_model(SurfaceField& surf,
                                 const ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ARTS_USER_ERROR_IF(
      not surf.contains(type), "Surface does not contain key value {}", type)

  update_field(surf[type].flat_view(),
               x[Range(x_start, x_size)],
               inverse_state,
               type,
               surf);
}

void SubsurfaceTarget::update_model(SubsurfaceField& subsurf,
                                    const ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ARTS_USER_ERROR_IF(not subsurf.contains(type),
                     "Subsurface does not contain key value {}",
                     type)

  update_field(subsurf[type].flat_view(),
               x[Range(x_start, x_size)],
               inverse_state,
               type,
               subsurf);
}

void LineTarget::update_model(AbsorptionBands& bands,
                              const ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  update_field(VectorView{type.get_value(bands)},
               x[Range(x_start, x_size)],
               inverse_state,
               type,
               bands);
}

void SensorTarget::update_model(ArrayOfSensorObsel& sens,
                                const ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  const ConstVectorView x_state = x[Range(x_start, x_size)];

  const SensorObsel& v = sens.at(type.measurement_elem);
  const Size N         = v.flat_size(type.type);

  if (inverse_state) {
    const Vector x_inverse = inverse_state(Vector{x_state}, sens);

    ARTS_USER_ERROR_IF(x_inverse.size() not_eq N,
                       R"(Size mismatch in inverse transformation for target {}.

Cannot set the field from the inverse transformation of the model state vector.

The original size of the target is            : {}.
The size of the inverse transformed target is : {}.
The expected size of the target is            : {}.
)",
                       type,
                       N,
                       x_inverse.size(),
                       x_state.size())
    unflatten(sens, x_inverse, v, type.type);
  } else {
    ARTS_USER_ERROR_IF(N not_eq x_state.size(),
                       R"(Size mismatch for target {}.

Cannot set the field from the model state vector.

The size of the target is          : {}.
The expected size of the target is : {}.
)",
                       type,
                       N,
                       x_state.size())

    unflatten(sens, x_state, v, type.type);
  }
}

////////////////////////////////////////////////////////////////////////
/// Update the model state vector
////////////////////////////////////////////////////////////////////////

void ErrorTarget::update_state(VectorView x, const ConstVectorView meas) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  update_x(x[Range(x_start, x_size)], meas, transform_state, type, meas);
}

void AtmTarget::update_state(VectorView x, const AtmField& atm) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ARTS_USER_ERROR_IF(
      not atm.contains(type), "Atmosphere does not contain key value {}", type)

  update_x(x[Range(x_start, x_size)],
           atm[type].flat_view(),
           transform_state,
           type,
           atm);
}

void SurfaceTarget::update_state(VectorView x, const SurfaceField& surf) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ARTS_USER_ERROR_IF(
      not surf.contains(type), "Surface does not contain key value {}", type)

  update_x(x[Range(x_start, x_size)],
           surf[type].flat_view(),
           transform_state,
           type,
           surf);
}

void SubsurfaceTarget::update_state(VectorView x,
                                    const SubsurfaceField& subsurf) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ARTS_USER_ERROR_IF(not subsurf.contains(type),
                     "Subsurface does not contain key value {}",
                     type)

  update_x(x[Range(x_start, x_size)],
           subsurf[type].flat_view(),
           transform_state,
           type,
           subsurf);
}

void LineTarget::update_state(VectorView x,
                              const AbsorptionBands& bands) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  update_x(x[Range(x_start, x_size)],
           ConstVectorView{type.get_value(bands)},
           transform_state,
           type,
           bands);
}

void SensorTarget::update_state(VectorView x,
                                const ArrayOfSensorObsel& sens) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  const SensorObsel& v = sens.at(type.measurement_elem);
  const Size N         = v.flat_size(type.type);

  const auto call = [&](VectorView x_field) {
    v.flat(x_field, type.type);

    update_x(x[Range(x_start, x_size)], x_field, transform_state, type, sens);
  };

  if (N == x_size) {
    call(x[Range(x_start, x_size)]);
  } else {
    Vector x_field(N);
    call(x_field);
  }
}

////////////////////////////////////////////////////////////////////////
/// Update the Jacobian
////////////////////////////////////////////////////////////////////////

void ErrorTarget::update_jac(MatrixView dy, const ConstVectorView meas) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (transform_jacobian) {
    update_dy(
        dy[joker, Range(x_start, x_size)], transform_jacobian, type, meas);
  }
}

void AtmTarget::update_jac(MatrixView dy, const AtmField& atm) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (transform_jacobian) {
    update_dy(dy[joker, Range(x_start, x_size)], transform_jacobian, type, atm);
  }
}

void SurfaceTarget::update_jac(MatrixView dy, const SurfaceField& surf) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (transform_jacobian) {
    update_dy(
        dy[joker, Range(x_start, x_size)], transform_jacobian, type, surf);
  }
}

void SubsurfaceTarget::update_jac(MatrixView dy,
                                  const SubsurfaceField& subsurf) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (transform_jacobian) {
    update_dy(
        dy[joker, Range(x_start, x_size)], transform_jacobian, type, subsurf);
  }
}

void LineTarget::update_jac(MatrixView dy, const AbsorptionBands& bands) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (transform_jacobian) {
    update_dy(
        dy[joker, Range(x_start, x_size)], transform_jacobian, type, bands);
  }
}

void SensorTarget::update_jac(MatrixView dy,
                              const ArrayOfSensorObsel& sens) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (transform_jacobian) {
    update_dy(
        dy[joker, Range(x_start, x_size)], transform_jacobian, type, sens);
  }
}

////////////////////////////////////////////////////////////////////////
/// ...
////////////////////////////////////////////////////////////////////////

bool AtmTarget::is_wind() const {
  return type == AtmKey::wind_u or type == AtmKey::wind_v or
         type == AtmKey::wind_w;
}

const std::vector<AtmTarget>& Targets::atm() const {
  return target<AtmTarget>();
}

const std::vector<SurfaceTarget>& Targets::surf() const {
  return target<SurfaceTarget>();
}

const std::vector<SubsurfaceTarget>& Targets::subsurf() const {
  return target<SubsurfaceTarget>();
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

std::vector<SubsurfaceTarget>& Targets::subsurf() {
  return target<SubsurfaceTarget>();
}

std::vector<LineTarget>& Targets::line() { return target<LineTarget>(); }

std::vector<SensorTarget>& Targets::sensor() { return target<SensorTarget>(); }

std::vector<ErrorTarget>& Targets::error() { return target<ErrorTarget>(); }

void Targets::finalize(const AtmField& atmospheric_field,
                       const SurfaceField& surface_field,
                       const SubsurfaceField& subsurface_field,
                       const AbsorptionBands&,
                       const ArrayOfSensorObsel& measurement_sensor) {
  const Size natm     = atm().size();
  const Size nsurf    = surf().size();
  const Size nsubsurf = subsurf().size();
  const Size nline    = line().size();
  const Size nsensor  = sensor().size();

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

  for (Size i = 0; i < nsubsurf; i++) {
    SubsurfaceTarget& t = subsurf()[i];
    ARTS_USER_ERROR_IF(std::ranges::any_of(subsurf() | std::views::drop(i + 1),
                                           Cmp::eq(t.type),
                                           &SubsurfaceTarget::type),
                       "Multiple targets of the same type: {}",
                       t.type)
    t.x_start  = last_size;
    t.x_size   = subsurface_field[t.type].flat_view().size();
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

    t.x_start  = last_size;
    // t.x_size already known
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

SubsurfaceTarget& Targets::emplace_back(SubsurfaceKeyVal&& t, Numeric d) {
  return subsurf().emplace_back(std::move(t), d, target_count());
}

LineTarget& Targets::emplace_back(LblLineKey&& t, Numeric d) {
  return line().emplace_back(std::move(t), d, target_count());
}

SensorTarget& Targets::emplace_back(SensorKey&& t, Numeric d) {
  return sensor().emplace_back(std::move(t), d, target_count());
}

ErrorTarget& Targets::emplace_back(ErrorKey&& t, Numeric d) {
  return error().emplace_back(std::move(t), d, target_count());
}

AtmTarget& Targets::emplace_back(const AtmKeyVal& t, Numeric d) {
  return atm().emplace_back(t, d, target_count());
}

SurfaceTarget& Targets::emplace_back(const SurfaceKeyVal& t, Numeric d) {
  return surf().emplace_back(t, d, target_count());
}

SubsurfaceTarget& Targets::emplace_back(const SubsurfaceKeyVal& t, Numeric d) {
  return subsurf().emplace_back(t, d, target_count());
}

LineTarget& Targets::emplace_back(const LblLineKey& t, Numeric d) {
  return line().emplace_back(t, d, target_count());
}

SensorTarget& Targets::emplace_back(const SensorKey& t, Numeric d) {
  return sensor().emplace_back(t, d, target_count());
}

ErrorTarget& Targets::emplace_back(const ErrorKey& t, Numeric d) {
  return error().emplace_back(t, d, target_count());
}
}  // namespace Jacobian

void xml_io_stream<JacobianTargetType>::write(std::ostream& os,
                                              const JacobianTargetType& x,
                                              bofstream* pbofs,
                                              std::string_view name) {
  xml_write_to_stream(os, x.target, pbofs, name);
}

void xml_io_stream<JacobianTargetType>::read(std::istream& is,
                                             JacobianTargetType& x,
                                             bifstream* pbifs) {
  xml_read_from_stream(is, x.target, pbifs);
}

void xml_io_stream<Jacobian::AtmTarget>::write(std::ostream& os,
                                               const Jacobian::AtmTarget& x,
                                               bofstream* pbofs,
                                               std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.d, pbofs);
  xml_write_to_stream(os, x.target_pos, pbofs);
  xml_write_to_stream(os, x.x_start, pbofs);
  xml_write_to_stream(os, x.x_size, pbofs);
  xml_write_to_stream(os, x.set_state, pbofs);
  xml_write_to_stream(os, x.set_model, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Jacobian::AtmTarget>::read(std::istream& is,
                                              Jacobian::AtmTarget& x,
                                              bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.d, pbifs);
  xml_read_from_stream(is, x.target_pos, pbifs);
  xml_read_from_stream(is, x.x_start, pbifs);
  xml_read_from_stream(is, x.x_size, pbifs);
  xml_read_from_stream(is, x.set_state, pbifs);
  xml_read_from_stream(is, x.set_model, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Jacobian::SurfaceTarget>::write(
    std::ostream& os,
    const Jacobian::SurfaceTarget& x,
    bofstream* pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.d, pbofs);
  xml_write_to_stream(os, x.target_pos, pbofs);
  xml_write_to_stream(os, x.x_start, pbofs);
  xml_write_to_stream(os, x.x_size, pbofs);
  xml_write_to_stream(os, x.set_state, pbofs);
  xml_write_to_stream(os, x.set_model, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Jacobian::SurfaceTarget>::read(std::istream& is,
                                                  Jacobian::SurfaceTarget& x,
                                                  bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.d, pbifs);
  xml_read_from_stream(is, x.target_pos, pbifs);
  xml_read_from_stream(is, x.x_start, pbifs);
  xml_read_from_stream(is, x.x_size, pbifs);
  xml_read_from_stream(is, x.set_state, pbifs);
  xml_read_from_stream(is, x.set_model, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Jacobian::SubsurfaceTarget>::write(
    std::ostream& os,
    const Jacobian::SubsurfaceTarget& x,
    bofstream* pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.d, pbofs);
  xml_write_to_stream(os, x.target_pos, pbofs);
  xml_write_to_stream(os, x.x_start, pbofs);
  xml_write_to_stream(os, x.x_size, pbofs);
  xml_write_to_stream(os, x.set_state, pbofs);
  xml_write_to_stream(os, x.set_model, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Jacobian::SubsurfaceTarget>::read(
    std::istream& is, Jacobian::SubsurfaceTarget& x, bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.d, pbifs);
  xml_read_from_stream(is, x.target_pos, pbifs);
  xml_read_from_stream(is, x.x_start, pbifs);
  xml_read_from_stream(is, x.x_size, pbifs);
  xml_read_from_stream(is, x.set_state, pbifs);
  xml_read_from_stream(is, x.set_model, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Jacobian::LineTarget>::write(std::ostream& os,
                                                const Jacobian::LineTarget& x,
                                                bofstream* pbofs,
                                                std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.d, pbofs);
  xml_write_to_stream(os, x.target_pos, pbofs);
  xml_write_to_stream(os, x.x_start, pbofs);
  xml_write_to_stream(os, x.x_size, pbofs);
  xml_write_to_stream(os, x.set_state, pbofs);
  xml_write_to_stream(os, x.set_model, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Jacobian::LineTarget>::read(std::istream& is,
                                               Jacobian::LineTarget& x,
                                               bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.d, pbifs);
  xml_read_from_stream(is, x.target_pos, pbifs);
  xml_read_from_stream(is, x.x_start, pbifs);
  xml_read_from_stream(is, x.x_size, pbifs);
  xml_read_from_stream(is, x.set_state, pbifs);
  xml_read_from_stream(is, x.set_model, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Jacobian::SensorTarget>::write(
    std::ostream& os,
    const Jacobian::SensorTarget& x,
    bofstream* pbofs,
    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.d, pbofs);
  xml_write_to_stream(os, x.target_pos, pbofs);
  xml_write_to_stream(os, x.x_start, pbofs);
  xml_write_to_stream(os, x.x_size, pbofs);
  xml_write_to_stream(os, x.set_state, pbofs);
  xml_write_to_stream(os, x.set_model, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Jacobian::SensorTarget>::read(std::istream& is,
                                                 Jacobian::SensorTarget& x,
                                                 bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.d, pbifs);
  xml_read_from_stream(is, x.target_pos, pbifs);
  xml_read_from_stream(is, x.x_start, pbifs);
  xml_read_from_stream(is, x.x_size, pbifs);
  xml_read_from_stream(is, x.set_state, pbifs);
  xml_read_from_stream(is, x.set_model, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<ErrorKey>::write(std::ostream& os,
                                    const ErrorKey& x,
                                    bofstream* pbofs,
                                    std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.y_start, pbofs);
  xml_write_to_stream(os, x.y_size, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<ErrorKey>::read(std::istream& is,
                                   ErrorKey& x,
                                   bifstream* pbifs) {
  XMLTag tag;

  xml_read_from_stream(is, x.y_start, pbifs);
  xml_read_from_stream(is, x.y_size, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<Jacobian::ErrorTarget>::write(std::ostream& os,
                                                 const Jacobian::ErrorTarget& x,
                                                 bofstream* pbofs,
                                                 std::string_view name) {
  std::println(os, R"(<{0} name="{1}">)", type_name, name);

  xml_write_to_stream(os, x.type, pbofs);
  xml_write_to_stream(os, x.target_pos, pbofs);
  xml_write_to_stream(os, x.x_start, pbofs);
  xml_write_to_stream(os, x.x_size, pbofs);
  xml_write_to_stream(os, x.set_y, pbofs);
  xml_write_to_stream(os, x.set_x, pbofs);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<Jacobian::ErrorTarget>::read(std::istream& is,
                                                Jacobian::ErrorTarget& x,
                                                bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  xml_read_from_stream(is, x.type, pbifs);
  xml_read_from_stream(is, x.target_pos, pbifs);
  xml_read_from_stream(is, x.x_start, pbifs);
  xml_read_from_stream(is, x.x_size, pbifs);
  xml_read_from_stream(is, x.set_y, pbifs);
  xml_read_from_stream(is, x.set_x, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}

void xml_io_stream<JacobianTargets>::write(std::ostream& os,
                                           const JacobianTargets& x,
                                           bofstream* pbofs,
                                           std::string_view name) {
  std::println(os,
               R"(<{0} name="{1}" final="{2}">)",
               type_name,
               name,
               Index{x.finalized});

  xml_write_to_stream(os, x.targets, pbofs, name);

  std::println(os, R"(</{0}>)", type_name);
}

void xml_io_stream<JacobianTargets>::read(std::istream& is,
                                          JacobianTargets& x,
                                          bifstream* pbifs) {
  XMLTag tag;
  tag.read_from_stream(is);
  tag.check_name(type_name);

  Index fin{};
  tag.get_attribute_value("final", fin);

  x.finalized = fin != 0;

  xml_read_from_stream(is, x.targets, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
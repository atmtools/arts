#include "jacobian.h"

#include <compare.h>
#include <obsel.h>

#include <algorithm>
#include <utility>

namespace Jacobian {
void Targets::clear() {
  atm.clear();
  surf.clear();
  subsurf.clear();
  line.clear();
  sensor.clear();
  error.clear();
}

Size Targets::x_size() const {
  if (target_count() != 0 and not finalized)
    throw std::runtime_error("Not finalized.");

  const auto sz = [](const auto& x) -> Size {
    return x.overlap ? 0 : x.x_size;
  };

  return std::transform_reduce(
             atm.begin(), atm.end(), Size{0}, std::plus<>{}, sz) +
         std::transform_reduce(
             surf.begin(), surf.end(), Size{0}, std::plus<>{}, sz) +
         std::transform_reduce(
             subsurf.begin(), subsurf.end(), Size{0}, std::plus<>{}, sz) +
         std::transform_reduce(
             line.begin(), line.end(), Size{0}, std::plus<>{}, sz) +
         std::transform_reduce(
             sensor.begin(), sensor.end(), Size{0}, std::plus<>{}, sz) +
         std::transform_reduce(
             error.begin(), error.end(), Size{0}, std::plus<>{}, sz);
}

Size Targets::target_count() const {
  return atm.size() + surf.size() + subsurf.size() + line.size() +
         sensor.size() + error.size();
}

bool Targets::empty() const { return target_count() == 0; }

void Targets::throwing_check(Size xsize) const {
  const auto t_size = target_count();

  if (xsize != x_size())
    throw std::runtime_error(
        "The size of the x-vector does not match the size of the targets.");

  for (auto& a : atm) {
    if ((a.x_start + a.x_size) > xsize)
      throw std::runtime_error("x-vector out-of-bounds");
    if (t_size <= a.target_pos)
      throw std::runtime_error("target-vector out-of-bounds.");
  }

  for (auto& a : surf) {
    if ((a.x_start + a.x_size) > xsize)
      throw std::runtime_error("x-vector out-of-bounds");
    if (t_size <= a.target_pos)
      throw std::runtime_error("target-vector out-of-bounds.");
  }

  for (auto& a : subsurf) {
    if ((a.x_start + a.x_size) > xsize)
      throw std::runtime_error("x-vector out-of-bounds");
    if (t_size <= a.target_pos)
      throw std::runtime_error("target-vector out-of-bounds.");
  }

  for (auto& a : line) {
    if ((a.x_start + a.x_size) > xsize)
      throw std::runtime_error("x-vector out-of-bounds");
    if (t_size <= a.target_pos)
      throw std::runtime_error("target-vector out-of-bounds.");
  }

  for (auto& a : sensor) {
    if ((a.x_start + a.x_size) > xsize)
      throw std::runtime_error("x-vector out-of-bounds");
    if (t_size <= a.target_pos)
      throw std::runtime_error("target-vector out-of-bounds.");
  }

  for (auto& a : error) {
    if ((a.x_start + a.x_size) > xsize)
      throw std::runtime_error("x-vector out-of-bounds");
    if (t_size <= a.target_pos)
      throw std::runtime_error("target-vector out-of-bounds.");
  }
}

std::vector<AtmTarget>::const_iterator Targets::find(const AtmKeyVal& t) const {
  return stdr::find_if(atm,
                       [&t](const auto& target) { return target.type == t; });
}

std::vector<SurfaceTarget>::const_iterator Targets::find(
    const SurfaceKeyVal& t) const {
  return stdr::find_if(surf,
                       [&t](const auto& target) { return target.type == t; });
}

std::vector<SubsurfaceTarget>::const_iterator Targets::find(
    const SubsurfaceKeyVal& t) const {
  return stdr::find_if(subsurf,
                       [&t](const auto& target) { return target.type == t; });
}

std::vector<LineTarget>::const_iterator Targets::find(
    const LblLineKey& t) const {
  return stdr::find_if(line,
                       [&t](const auto& target) { return target.type == t; });
}

std::vector<SensorTarget>::const_iterator Targets::find(
    const SensorKey& t) const {
  return stdr::find_if(sensor,
                       [&t](const auto& target) { return target.type == t; });
}

std::vector<ErrorTarget>::const_iterator Targets::find(
    const ErrorKey& t) const {
  return stdr::find_if(error,
                       [&t](const auto& target) { return target.type == t; });
}

Index Targets::target_position(const AtmKeyVal& t) const {
  const auto it = find(t);
  return it != atm.end() ? it->target_pos : -1;
}

Index Targets::target_position(const SurfaceKeyVal& t) const {
  const auto it = find(t);
  return it != surf.end() ? it->target_pos : -1;
}

Index Targets::target_position(const SubsurfaceKeyVal& t) const {
  const auto it = find(t);
  return it != subsurf.end() ? it->target_pos : -1;
}

Index Targets::target_position(const LblLineKey& t) const {
  const auto it = find(t);
  return it != line.end() ? it->target_pos : -1;
}

Index Targets::target_position(const SensorKey& t) const {
  const auto it = find(t);
  return it != sensor.end() ? it->target_pos : -1;
}

Index Targets::target_position(const ErrorKey& t) const {
  const auto it = find(t);
  return it != error.end() ? it->target_pos : -1;
}

bool TargetType::operator==(const TargetType&) const = default;

std::string TargetType::type() const {
  return apply([](auto&) { return "AtmKeyVal"s; },
               [](auto&) { return "SurfaceKeyVal"s; },
               [](auto&) { return "SubsurfaceKeyVal"s; },
               [](auto&) { return "LblLineKey"s; },
               [](auto&) { return "SensorKey"s; },
               [](auto&) { return "ErrorKey"s; });
}

////////////////////////////////////////////////////////////////////////
/// Templates for doing the common work of updating fields, model state vectors, and Jacobians
////////////////////////////////////////////////////////////////////////

namespace {
template <typename Func, typename Key, typename Field>
void update_x(VectorView x,
              ConstVectorView y,
              Func&& transform_state,
              const Key& type,
              const Field& field) {
  if (transform_state) {
    const Vector x_from_y = transform_state(y, field);

    ARTS_USER_ERROR_IF(x_from_y.size() not_eq x.size(),
                       R"(Size mismatch in transformation for target {}.

Cannot set the model state vector from the transformation of the field.

The size of the target field is       : {}.
The size of the transformed target is : {}.
The expected size of the target is    : {}.
)",
                       type,
                       y.size(),
                       x_from_y.size(),
                       x.size())

    x = x_from_y;
  } else {
    ARTS_USER_ERROR_IF(y.size() not_eq x.size(),
                       R"(Size mismatch in Jacobian target for {}.

Cannot set the model state vector from the field.

The size of the target field is    : {}.
The expected size of the target is : {}.
)",
                       type,
                       y.size(),
                       x.size())

    x = y;
  }
}

template <typename Func, typename Key, typename Field>
void update_dy(StridedMatrixView dy,
               ConstVectorView x,
               Func&& inverse_jacobian,
               const Key& type,
               const Field& y) {
  const Matrix dy_transformed = inverse_jacobian(Matrix{dy}, x, y);

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
void update_field(VectorView y,
                  ConstVectorView x,
                  Func&& inverse_state,
                  const Key& type,
                  const Field& field) {
  if (inverse_state) {
    const Vector y_of_x = inverse_state(x, field);

    ARTS_USER_ERROR_IF(y_of_x.size() != y.size(),
                       R"(Size mismatch in inverse transformation for target {}.

Cannot set the field from the inverse transformation of the model state vector.

The original size of the target is            : {}.
The size of the inverse transformed target is : {}.
The expected size of the target is            : {}.
)",
                       type,
                       y.size(),
                       y_of_x.size(),
                       x.size())

    y = y_of_x;
  } else {
    ARTS_USER_ERROR_IF(y.size() not_eq x.size(),
                       R"(Size mismatch for target {}.

Cannot set the field from the model state vector.

The size of the target is          : {}.
The expected size of the target is : {}.
)",
                       type,
                       y.size(),
                       x.size())

    y = x;
  }
}
}  // namespace

////////////////////////////////////////////////////////////////////////
/// Update the fields
////////////////////////////////////////////////////////////////////////

void ErrorTarget::update_model(VectorView y, ConstVectorView x) const {
  ARTS_USER_ERROR_IF(
      x.size() < (x_start + x_size) or y.size() < (type.y_start + type.y_size),
      R"(Measurement Jacobian too small.

Expected minimum shape: [{}, {}]
Got shape:              [{}, {}]
)",
      type.y_start + type.y_size,
      x_start + x_size,
      y.size(),
      x.size())

  update_field(y[Range(type.y_start, type.y_size)],
               x[Range(x_start, x_size)],
               inverse_state,
               type,
               y[Range(type.y_start, type.y_size)]);
}

void AtmTarget::update_model(AtmField& atm, ConstVectorView x) const {
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

void SurfaceTarget::update_model(SurfaceField& surf, ConstVectorView x) const {
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
                                    ConstVectorView x) const {
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

void LineTarget::update_model(AbsorptionBands& bands, ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  update_field(VectorView{type.get_value(bands)},
               x[Range(x_start, x_size)],
               inverse_state,
               type,
               bands);
}

void SensorTarget::update_model(ArrayOfSensorObsel& sens,
                                ConstVectorView x) const {
  ARTS_USER_ERROR_IF(x.size() < (x_start + x_size),
                     "Model state vector is too small.")

  ConstVectorView x_state = x[Range(x_start, x_size)];

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

void ErrorTarget::update_state(VectorView x, ConstVectorView y) const {
  ARTS_USER_ERROR_IF(
      x.size() < (x_start + x_size) or y.size() < (type.y_start + type.y_size),
      R"(Measurement Jacobian too small.

Expected minimum shape: [{}, {}]
Got shape:              [{}, {}]
)",
      type.y_start + type.y_size,
      x_start + x_size,
      y.size(),
      x.size())

  update_x(x[Range(x_start, x_size)],
           y[Range(type.y_start, type.y_size)],
           transform_state,
           type,
           y[Range(type.y_start, type.y_size)]);
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

void ErrorTarget::update_jac(MatrixView dy,
                             ConstVectorView x,
                             ConstVectorView y) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size) or
          static_cast<Size>(dy.nrows()) < (type.y_start + type.y_size),
      R"(Measurement Jacobian too small.

Expected minimum shape: [{}, {}]
Got shape:              {:B,}
)",
      type.y_start + type.y_size,
      x_start + x_size,
      dy.shape())

  if (inverse_jacobian) {
    update_dy(dy[Range(type.y_start, type.y_size), Range(x_start, x_size)],
              x[Range(x_start, x_size)],
              inverse_jacobian,
              type,
              y[Range(type.y_start, type.y_size)]);
  }
}

void AtmTarget::update_jac(MatrixView dy,
                           ConstVectorView x,
                           const AtmField& atm) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (inverse_jacobian) {
    update_dy(dy[joker, Range(x_start, x_size)],
              x[Range(x_start, x_size)],
              inverse_jacobian,
              type,
              atm);
  }
}

void SurfaceTarget::update_jac(MatrixView dy,
                               ConstVectorView x,
                               const SurfaceField& surf) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (inverse_jacobian) {
    update_dy(dy[joker, Range(x_start, x_size)],
              x[Range(x_start, x_size)],
              inverse_jacobian,
              type,
              surf);
  }
}

void SubsurfaceTarget::update_jac(MatrixView dy,
                                  ConstVectorView x,
                                  const SubsurfaceField& subsurf) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (inverse_jacobian) {
    update_dy(dy[joker, Range(x_start, x_size)],
              x[Range(x_start, x_size)],
              inverse_jacobian,
              type,
              subsurf);
  }
}

void LineTarget::update_jac(MatrixView dy,
                            ConstVectorView x,
                            const AbsorptionBands& bands) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (inverse_jacobian) {
    update_dy(dy[joker, Range(x_start, x_size)],
              x[Range(x_start, x_size)],
              inverse_jacobian,
              type,
              bands);
  }
}

void SensorTarget::update_jac(MatrixView dy,
                              ConstVectorView x,
                              const ArrayOfSensorObsel& sens) const {
  ARTS_USER_ERROR_IF(
      static_cast<Size>(dy.ncols()) < (x_start + x_size),
      "Model state vector dimension of the Jacobian is too small.")

  if (inverse_jacobian) {
    update_dy(dy[joker, Range(x_start, x_size)],
              x[Range(x_start, x_size)],
              inverse_jacobian,
              type,
              sens);
  }
}

////////////////////////////////////////////////////////////////////////
/// ...
////////////////////////////////////////////////////////////////////////

bool is_wind(const AtmTarget& t) {
  return t.type == AtmKey::wind_u or t.type == AtmKey::wind_v or
         t.type == AtmKey::wind_w;
}

bool is_mag(const AtmTarget& t) {
  return t.type == AtmKey::mag_u or t.type == AtmKey::mag_v or
         t.type == AtmKey::mag_w;
}

void Targets::finalize(const AtmField& atm_field,
                       const SurfaceField& surf_field,
                       const SubsurfaceField& subsurf_field,
                       const AbsorptionBands&,
                       const ArrayOfSensorObsel& measurement_sensor) {
  const Size natm     = atm.size();
  const Size nsurf    = surf.size();
  const Size nsubsurf = subsurf.size();
  const Size nline    = line.size();
  const Size nsensor  = sensor.size();
  const Size nerror   = error.size();

  Size last_size = 0;

  for (Size i = 0; i < natm; i++) {
    AtmTarget& t = atm[i];
    ARTS_USER_ERROR_IF(
        stdr::any_of(
            atm | stdv::drop(i + 1), Cmp::eq(t.type), &AtmTarget::type),
        "Multiple targets of the same type: {}",
        t.type)

    if (t.overlap) {
      const auto f = stdr::find(atm, t.overlap_key, &AtmTarget::type);
      ARTS_USER_ERROR_IF(
          static_cast<Index>(i) <= stdr::distance(atm.begin(), f),
          "Overlap target {} not found prior in targets.",
          t.overlap_key)
      t.x_start = f->x_start;
      t.x_size  = f->x_size;
    } else {
      t.x_start  = last_size;
      t.x_size   = atm_field[t.type].flat_view().size();
      last_size += t.x_size;
    }
  }

  for (Size i = 0; i < nsurf; i++) {
    SurfaceTarget& t = surf[i];
    ARTS_USER_ERROR_IF(
        stdr::any_of(
            surf | stdv::drop(i + 1), Cmp::eq(t.type), &SurfaceTarget::type),
        "Multiple targets of the same type: {}",
        t.type)

    if (t.overlap) {
      const auto f = stdr::find(surf, t.overlap_key, &SurfaceTarget::type);
      ARTS_USER_ERROR_IF(
          static_cast<Index>(i) <= stdr::distance(surf.begin(), f),
          "Overlap target {} not found prior in targets.",
          t.overlap_key)
      t.x_start = f->x_start;
      t.x_size  = f->x_size;
    } else {
      t.x_start  = last_size;
      t.x_size   = surf_field[t.type].flat_view().size();
      last_size += t.x_size;
    }
  }

  for (Size i = 0; i < nsubsurf; i++) {
    SubsurfaceTarget& t = subsurf[i];
    ARTS_USER_ERROR_IF(stdr::any_of(subsurf | stdv::drop(i + 1),
                                    Cmp::eq(t.type),
                                    &SubsurfaceTarget::type),
                       "Multiple targets of the same type: {}",
                       t.type)

    if (t.overlap) {
      const auto f =
          stdr::find(subsurf, t.overlap_key, &SubsurfaceTarget::type);
      ARTS_USER_ERROR_IF(
          static_cast<Index>(i) <= stdr::distance(subsurf.begin(), f),
          "Overlap target {} not found prior in targets.",
          t.overlap_key)
      t.x_start = f->x_start;
      t.x_size  = f->x_size;
    } else {
      t.x_start  = last_size;
      t.x_size   = subsurf_field[t.type].flat_view().size();
      last_size += t.x_size;
    }
  }

  for (Size i = 0; i < nline; i++) {
    LineTarget& t = line[i];
    ARTS_USER_ERROR_IF(
        stdr::any_of(
            line | stdv::drop(i + 1), Cmp::eq(t.type), &LineTarget::type),
        "Multiple targets of the same type: {}",
        t.type)

    if (t.overlap) {
      const auto f = stdr::find(line, t.overlap_key, &LineTarget::type);
      ARTS_USER_ERROR_IF(
          static_cast<Index>(i) <= stdr::distance(line.begin(), f),
          "Overlap target {} not found prior in targets.",
          t.overlap_key)
      t.x_start = f->x_start;
      t.x_size  = f->x_size;
    } else {
      t.x_start  = last_size;
      t.x_size   = 1;
      last_size += t.x_size;
    }
  }

  for (auto& elem : sensor) {
    ARTS_USER_ERROR_IF(static_cast<Size>(elem.type.measurement_elem) >=
                           measurement_sensor.size(),
                       "Bad sensor elements {}, out-of-bounds",
                       elem);
  }

  for (Size i = 0; i < nsensor; i++) {
    SensorTarget& t = sensor[i];
    ARTS_USER_ERROR_IF(
        stdr::any_of(
            sensor | stdv::drop(i + 1),
            [&](const SensorKey& key) {
              if (t.type.type != key.type) return false;

              const Index elem1 = t.type.measurement_elem;
              const Index elem2 = key.measurement_elem;

              switch (key.type) {
                using enum SensorKeyType;
                case freq:
                  return measurement_sensor[elem1].f_grid_ptr() ==
                         measurement_sensor[elem2].f_grid_ptr();
                case zen:
                case azi:
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

    if (t.overlap) {
      const auto f = stdr::find(sensor, t.overlap_key, &SensorTarget::type);
      ARTS_USER_ERROR_IF(
          static_cast<Index>(i) <= stdr::distance(sensor.begin(), f),
          "Overlap target {} not found prior in targets.",
          t.overlap_key)
      t.x_start = f->x_start;
      t.x_size  = f->x_size;
    } else {
      t.x_start = last_size;
      // t.x_size already known
      last_size += t.x_size;
    }
  }

  for (Size i = 0; i < nerror; i++) {
    ErrorTarget& t = error[i];

    if (t.overlap) {
      const auto f = stdr::find(error, t.overlap_key, &ErrorTarget::type);
      ARTS_USER_ERROR_IF(
          static_cast<Index>(i) <= stdr::distance(error.begin(), f),
          "Overlap target {} not found prior in targets.",
          t.overlap_key)
      t.x_start = f->x_start;
      t.x_size  = f->x_size;
    } else {
      t.x_start = last_size;
      // t.x_size already known
      last_size += t.x_size;
    }

    ARTS_USER_ERROR_IF(
        t.type.y_start + t.type.y_size > measurement_sensor.size(),
        R"(An error target is out of bounds of the y-vector.

The error target is defined to operator on the measurement_vec in the range {},
but the measurement_vec will only have {} element by the measurement_sensor.
)",
        t.type,
        measurement_sensor.size())
  }

  finalized = true;

  throwing_check(last_size);
}

AtmTarget& Targets::emplace_back(AtmKeyVal&& t, Numeric d) {
  const auto n = target_count();
  return atm.emplace_back(std::move(t), d, n);
}

SurfaceTarget& Targets::emplace_back(SurfaceKeyVal&& t, Numeric d) {
  const auto n = target_count();
  return surf.emplace_back(std::move(t), d, n);
}

SubsurfaceTarget& Targets::emplace_back(SubsurfaceKeyVal&& t, Numeric d) {
  const auto n = target_count();
  return subsurf.emplace_back(std::move(t), d, n);
}

LineTarget& Targets::emplace_back(LblLineKey&& t, Numeric d) {
  const auto n = target_count();
  return line.emplace_back(std::move(t), d, n);
}

SensorTarget& Targets::emplace_back(SensorKey&& t, Numeric d) {
  const auto n = target_count();
  return sensor.emplace_back(std::move(t), d, n);
}

ErrorTarget& Targets::emplace_back(ErrorKey&& t, Numeric d) {
  const auto n = target_count();
  return error.emplace_back(std::move(t), d, n);
}

AtmTarget& Targets::emplace_back(const AtmKeyVal& t, Numeric d) {
  return atm.emplace_back(t, d, target_count());
}

SurfaceTarget& Targets::emplace_back(const SurfaceKeyVal& t, Numeric d) {
  return surf.emplace_back(t, d, target_count());
}

SubsurfaceTarget& Targets::emplace_back(const SubsurfaceKeyVal& t, Numeric d) {
  return subsurf.emplace_back(t, d, target_count());
}

LineTarget& Targets::emplace_back(const LblLineKey& t, Numeric d) {
  return line.emplace_back(t, d, target_count());
}

SensorTarget& Targets::emplace_back(const SensorKey& t, Numeric d) {
  return sensor.emplace_back(t, d, target_count());
}

ErrorTarget& Targets::emplace_back(const ErrorKey& t, Numeric d) {
  return error.emplace_back(t, d, target_count());
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

void xml_io_stream<JacobianTargets>::write(std::ostream& os,
                                           const JacobianTargets& x,
                                           bofstream* pbofs,
                                           std::string_view name) {
  XMLTag tag(type_name, "name", name, "final", Index{x.finalized});
  tag.write_to_stream(os);

  xml_write_to_stream(os, x.atm, pbofs, "atm");
  xml_write_to_stream(os, x.surf, pbofs, "surf");
  xml_write_to_stream(os, x.subsurf, pbofs, "subsurf");
  xml_write_to_stream(os, x.line, pbofs, "line");
  xml_write_to_stream(os, x.sensor, pbofs, "sensor");
  xml_write_to_stream(os, x.error, pbofs, "error");

  tag.write_to_end_stream(os);
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

  xml_read_from_stream(is, x.atm, pbifs);
  xml_read_from_stream(is, x.surf, pbifs);
  xml_read_from_stream(is, x.subsurf, pbifs);
  xml_read_from_stream(is, x.line, pbifs);
  xml_read_from_stream(is, x.sensor, pbifs);
  xml_read_from_stream(is, x.error, pbifs);

  tag.read_from_stream(is);
  tag.check_end_name(type_name);
}
#pragma once

#include <atm.h>
#include <lbl.h>
#include <sensor.h>
#include <subsurf.h>
#include <surf.h>
#include <xml_io_stream_functional.h>

#include <concepts>
#include <functional>

struct ErrorKey;

//! The transformations of the targets, see jac_rel.h, jac_log.h, jac_logrel.h, jac_rh.h, and jac_polyfit.h
struct relfwd;
struct relinv;
struct logfwd;
struct loginv;
struct logrelfwd;
struct logrelinv;
struct rhfwd;
struct rhinv;
struct polyfit_t;
struct polyinv_t;
struct polyfit_sensor_offset_t;
struct polyinv_sensor_offset_t;

namespace Jacobian {
struct AtmTarget;
struct SurfaceTarget;
struct SubsurfaceTarget;
struct LineTarget;
struct SensorTarget;
struct ErrorTarget;

//! The fields whose targets can be relative or logarithmic
template <typename T>
concept transformable_field =
    std::same_as<T, AtmField> or std::same_as<T, SurfaceField> or std::same_as<T, SubsurfaceField>;

using cv          = ConstVectorView;
using cm          = ConstMatrixView;
using atm_vec     = std::function<Vector(cv, const AtmField&)>;
using atm_mat     = std::function<Matrix(cm, cv, const AtmField&)>;
using surf_vec    = std::function<Vector(cv, const SurfaceField&)>;
using surf_mat    = std::function<Matrix(cm, cv, const SurfaceField&)>;
using subsurf_vec = std::function<Vector(cv, const SubsurfaceField&)>;
using subsurf_mat = std::function<Matrix(cm, cv, const SubsurfaceField&)>;
using line_vec    = std::function<Vector(cv, const AbsorptionBands&)>;
using line_mat    = std::function<Matrix(cm, cv, const AbsorptionBands&)>;
using sensor_vec  = std::function<Vector(cv, const ArrayOfSensor&)>;
using sensor_mat  = std::function<Matrix(cm, cv, const ArrayOfSensor&)>;
using error_vec   = std::function<Vector(cv, cv)>;
using error_mat   = std::function<Matrix(cm, cv, cv)>;
}  // namespace Jacobian

/* The transformations that the functions of the targets can hold, by XML
 *
 * The XML IO of every function type is declared here, next to the type, and
 * defined in jacobian_names.cc, the only place where the transformation structs
 * must be complete.  Each function type lists its structs, and those that hold
 * no transformation yet list none, so the place to add one is clear.
 */

template <> struct xml_io_stream_functional<Jacobian::atm_vec> {
  using func_t    = Vector (*)(ConstVectorView, const AtmField&);
  using structs_t = std::variant<relinv, relfwd, loginv, logfwd, logrelfwd, logrelinv, rhfwd, rhinv>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::atm_mat> {
  using func_t    = Matrix (*)(ConstMatrixView, ConstVectorView, const AtmField&);
  using structs_t = std::variant<relinv, loginv, logrelinv, rhinv>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::surf_vec> {
  using func_t    = Vector (*)(ConstVectorView, const SurfaceField&);
  using structs_t = std::variant<relinv, relfwd, loginv, logfwd, logrelfwd, logrelinv>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::surf_mat> {
  using func_t    = Matrix (*)(ConstMatrixView, ConstVectorView, const SurfaceField&);
  using structs_t = std::variant<relinv, loginv, logrelinv>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::subsurf_vec> {
  using func_t    = Vector (*)(ConstVectorView, const SubsurfaceField&);
  using structs_t = std::variant<relinv, relfwd, loginv, logfwd, logrelfwd, logrelinv>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::subsurf_mat> {
  using func_t    = Matrix (*)(ConstMatrixView, ConstVectorView, const SubsurfaceField&);
  using structs_t = std::variant<relinv, loginv, logrelinv>;
  static constexpr std::array<func_t*, 0> funcs{};
};

//! Line targets hold no transformation yet, so their functions are always empty
template <> struct xml_io_stream_functional<Jacobian::line_vec> {
  using func_t    = Vector (*)(ConstVectorView, const AbsorptionBands&);
  using structs_t = void;
  static constexpr std::array<func_t*, 0> funcs{};
};

//! Line targets hold no transformation yet, so their functions are always empty
template <> struct xml_io_stream_functional<Jacobian::line_mat> {
  using func_t    = Matrix (*)(ConstMatrixView, ConstVectorView, const AbsorptionBands&);
  using structs_t = void;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::sensor_vec> {
  using func_t    = Vector (*)(ConstVectorView, const ArrayOfSensor&);
  using structs_t = std::variant<polyfit_sensor_offset_t, polyinv_sensor_offset_t>;
  static constexpr std::array<func_t*, 0> funcs{};
};

//! The sensor polynomial offset perturbs its coefficients directly, so no inverse Jacobian is held yet
template <> struct xml_io_stream_functional<Jacobian::sensor_mat> {
  using func_t    = Matrix (*)(ConstMatrixView, ConstVectorView, const ArrayOfSensor&);
  using structs_t = void;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::error_vec> {
  using func_t    = Vector (*)(ConstVectorView, ConstVectorView);
  using structs_t = std::variant<polyfit_t, polyinv_t>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream_functional<Jacobian::error_mat> {
  using func_t    = Matrix (*)(ConstMatrixView, ConstVectorView, ConstVectorView);
  using structs_t = std::variant<polyinv_t>;
  static constexpr std::array<func_t*, 0> funcs{};
};

template <> struct xml_io_stream<Jacobian::atm_vec> : xml_io_stream_function_declaration<Jacobian::atm_vec> {};
template <> struct xml_io_stream<Jacobian::atm_mat> : xml_io_stream_function_declaration<Jacobian::atm_mat> {};
template <> struct xml_io_stream<Jacobian::surf_vec> : xml_io_stream_function_declaration<Jacobian::surf_vec> {};
template <> struct xml_io_stream<Jacobian::surf_mat> : xml_io_stream_function_declaration<Jacobian::surf_mat> {};
template <> struct xml_io_stream<Jacobian::subsurf_vec> : xml_io_stream_function_declaration<Jacobian::subsurf_vec> {};
template <> struct xml_io_stream<Jacobian::subsurf_mat> : xml_io_stream_function_declaration<Jacobian::subsurf_mat> {};
template <> struct xml_io_stream<Jacobian::line_vec> : xml_io_stream_function_declaration<Jacobian::line_vec> {};
template <> struct xml_io_stream<Jacobian::line_mat> : xml_io_stream_function_declaration<Jacobian::line_mat> {};
template <> struct xml_io_stream<Jacobian::sensor_vec> : xml_io_stream_function_declaration<Jacobian::sensor_vec> {};
template <> struct xml_io_stream<Jacobian::sensor_mat> : xml_io_stream_function_declaration<Jacobian::sensor_mat> {};
template <> struct xml_io_stream<Jacobian::error_vec> : xml_io_stream_function_declaration<Jacobian::error_vec> {};
template <> struct xml_io_stream<Jacobian::error_mat> : xml_io_stream_function_declaration<Jacobian::error_mat> {};

extern template struct xml_io_stream_function_declaration<Jacobian::atm_vec>;
extern template struct xml_io_stream_function_declaration<Jacobian::atm_mat>;
extern template struct xml_io_stream_function_declaration<Jacobian::surf_vec>;
extern template struct xml_io_stream_function_declaration<Jacobian::surf_mat>;
extern template struct xml_io_stream_function_declaration<Jacobian::subsurf_vec>;
extern template struct xml_io_stream_function_declaration<Jacobian::subsurf_mat>;
extern template struct xml_io_stream_function_declaration<Jacobian::line_vec>;
extern template struct xml_io_stream_function_declaration<Jacobian::line_mat>;
extern template struct xml_io_stream_function_declaration<Jacobian::sensor_vec>;
extern template struct xml_io_stream_function_declaration<Jacobian::sensor_mat>;
extern template struct xml_io_stream_function_declaration<Jacobian::error_vec>;
extern template struct xml_io_stream_function_declaration<Jacobian::error_mat>;

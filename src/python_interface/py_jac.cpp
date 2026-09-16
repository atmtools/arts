#include <python_interface.h>

#include "atm.h"
#include "hpy_arts.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "nanobind/stl/bind_vector.h"
#include "nanobind/stl/function.h"
#include "nanobind/stl/variant.h"
#include "surf.h"

namespace Python {
void py_jac(py::module_& m) try {
  py::class_<SensorKey>(m, "SensorKey")
      .def(py::init<>(), "Create a sensor Jacobian key.")
      .def_rw("type", &SensorKey::type, "Sensor quantity.\n\n.. :class:`~pyarts3.arts.SensorKeyType`")
      .def_rw("sensor_elem", &SensorKey::sensor_elem, "Sensor element index.\n\n.. :class:`int`")
      .def_rw("measurement_elem", &SensorKey::measurement_elem, "Measurement element index.\n\n.. :class:`int`")
      .doc() = "Identifies a sensor quantity and its sensor and measurement elements.";

  py::class_<ErrorKey> errkey(m, "ErrorKey");
  generic_interface(errkey);
  errkey.doc() = "Error key";
  errkey.def_rw(
      "y_start", &ErrorKey::y_start, "Start index of target in measurement vector\n\n.. :class:`~pyarts3.arts.Index`");
  errkey.def_rw("y_size", &ErrorKey::y_size, "Size of target in measurement vector\n\n.. :class:`~pyarts3.arts.Index`");

  py::class_<Jacobian::AtmTarget> atm(m, "JacobianAtmTarget");
  generic_interface(atm);
  atm.doc() = "Atmospheric target";
  atm.def_ro(
      "type",
      &Jacobian::AtmTarget::type,
      "Type of target\n\n.. :class:`~pyarts3.arts.AtmKey | ~pyarts3.arts.SpeciesEnum | ~pyarts3.arts.SpeciesIsotope | ~pyarts3.arts.QuantumLevelIdentifier | ~pyarts3.arts.ScatteringSpeciesProperty`");
  atm.def_ro("d", &Jacobian::AtmTarget::d, "Perturbation magnitude\n\n.. :class:`~pyarts3.arts.Numeric`");
  atm.def_ro("target_pos", &Jacobian::AtmTarget::target_pos, "Target position\n\n.. :class:`~pyarts3.arts.Index`");
  atm.def_ro("x_start",
             &Jacobian::AtmTarget::x_start,
             "Start index of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");
  atm.def_ro(
      "x_size", &Jacobian::AtmTarget::x_size, "Size of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");

  py::class_<Jacobian::SurfaceTarget> surf(m, "JacobianSurfaceTarget");
  generic_interface(surf);
  surf.doc() = "Surface target";
  surf.def_ro("type",
              &Jacobian::SurfaceTarget::type,
              "Type of target\n\n.. :class:`~pyarts3.arts.SurfaceKey | ~pyarts3.arts.SurfacePropertyTag`");
  surf.def_ro("d", &Jacobian::SurfaceTarget::d, "Perturbation magnitude\n\n.. :class:`~pyarts3.arts.Numeric`");
  surf.def_ro("target_pos", &Jacobian::SurfaceTarget::target_pos, "Target position\n\n.. :class:`~pyarts3.arts.Index`");
  surf.def_ro("x_start",
              &Jacobian::SurfaceTarget::x_start,
              "Start index of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");
  surf.def_ro(
      "x_size", &Jacobian::SurfaceTarget::x_size, "Size of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");

  py::class_<Jacobian::SubsurfaceTarget> subsurf(m, "JacobianSubsurfaceTarget");
  generic_interface(subsurf);
  subsurf.doc() = "Subsurface target";
  subsurf.def_ro("type",
                 &Jacobian::SubsurfaceTarget::type,
                 "Type of target\n\n.. :class:`~pyarts3.arts.SubsurfaceKey | ~pyarts3.arts.SubsurfacePropertyTag`");
  subsurf.def_ro("d", &Jacobian::SubsurfaceTarget::d, "Perturbation magnitude\n\n.. :class:`~pyarts3.arts.Numeric`");
  subsurf.def_ro(
      "target_pos", &Jacobian::SubsurfaceTarget::target_pos, "Target position\n\n.. :class:`~pyarts3.arts.Index`");
  subsurf.def_ro("x_start",
                 &Jacobian::SubsurfaceTarget::x_start,
                 "Start index of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");
  subsurf.def_ro("x_size",
                 &Jacobian::SubsurfaceTarget::x_size,
                 "Size of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");

  py::class_<Jacobian::LineTarget> line(m, "JacobianLineTarget");
  generic_interface(line);
  line.doc() = "Line target";
  line.def_ro("type", &Jacobian::LineTarget::type, "Type of target\n\n.. :class:`~pyarts3.arts.lbl.line_key`");
  line.def_ro("d", &Jacobian::LineTarget::d, "Perturbation magnitude\n\n.. :class:`~pyarts3.arts.Numeric`");
  line.def_ro("target_pos", &Jacobian::LineTarget::target_pos, "Target position\n\n.. :class:`~pyarts3.arts.Index`");
  line.def_ro("x_start",
              &Jacobian::LineTarget::x_start,
              "Start index of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");
  line.def_ro(
      "x_size", &Jacobian::LineTarget::x_size, "Size of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");

  py::class_<Jacobian::SensorTarget> sensor(m, "JacobianSensorTarget");
  generic_interface(sensor);
  sensor.doc() = "Sensor target";
  sensor.def_ro("type", &Jacobian::SensorTarget::type, "Type of target\n\n.. :class:`~pyarts3.arts.SensorKey`");
  sensor.def_ro("d", &Jacobian::SensorTarget::d, "Perturbation magnitude\n\n.. :class:`~pyarts3.arts.Numeric`");
  sensor.def_ro(
      "target_pos", &Jacobian::SensorTarget::target_pos, "Target position\n\n.. :class:`~pyarts3.arts.Index`");
  sensor.def_ro("x_start",
                &Jacobian::SensorTarget::x_start,
                "Start index of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");
  sensor.def_ro(
      "x_size", &Jacobian::SensorTarget::x_size, "Size of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");

  py::class_<Jacobian::ErrorTarget> error(m, "JacobianErrorTarget");
  generic_interface(error);
  error.doc() = "Error target";
  error.def_ro("type", &Jacobian::ErrorTarget::type, "Type of target\n\n.. :class:`~pyarts3.arts.ErrorKey`");
  error.def_ro("target_pos", &Jacobian::ErrorTarget::target_pos, "Target position\n\n.. :class:`~pyarts3.arts.Index`");
  error.def_ro("x_start",
               &Jacobian::ErrorTarget::x_start,
               "Start index of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");
  error.def_ro(
      "x_size", &Jacobian::ErrorTarget::x_size, "Size of target in state vector\n\n.. :class:`~pyarts3.arts.Index`");

  const auto python_callable = [](py::object fn) {
    if (not PyCallable_Check(fn.ptr())) throw py::type_error("Jacobian transformation must be callable or None");
    return std::shared_ptr<py::object>(new py::object(std::move(fn)), [](py::object* ptr) {
      if (py::is_alive()) {
        py::gil_scoped_acquire guard;
        delete ptr;
      } else {
        ptr->release();
        delete ptr;
      }
    });
  };

  const auto bind_transformations = [python_callable]<typename Target, typename Context>(py::class_<Target>& target,
                                                                                         std::type_identity<Context>) {
    target.def_prop_rw(
        "transform_state",
        [](const Target& self) -> py::object {
          if (not self.transform_state) return py::none();
          if constexpr (std::same_as<Context, ConstVectorView>) {
            return py::cpp_function(
                [fn = self.transform_state](const Vector& state, const Vector& context) { return fn(state, context); });
          } else {
            return py::cpp_function([fn = self.transform_state](const Vector& state, const Context& context) {
              return fn(state, context);
            });
          }
        },
        [python_callable](Target& self, py::object fn) {
          if (fn.is_none()) {
            self.transform_state = {};
            return;
          }
          auto callable        = python_callable(std::move(fn));
          self.transform_state = [callable = std::move(callable)](ConstVectorView state,
                                                                  const Context&  context) -> Vector {
            py::gil_scoped_acquire guard;
            py::object             pycontext;
            if constexpr (std::same_as<Context, ConstVectorView>) {
              pycontext = py::cast(Vector{context});
            } else {
              pycontext = py::cast(&context, py::rv_policy::reference);
            }
            py::object result = (*callable)(Vector{state}, std::move(pycontext));
            return py::cast<Vector>(py::type<Vector>()(result));
          };
        },
        R"doc(Transform native target values into the model state.

The callable signature is ``(native_state, context) -> model_state``.  The
context is the complete field or data object that owns the target.

.. :class:`~collections.abc.Callable`)doc");
    target.def_prop_rw(
        "inverse_state",
        [](const Target& self) -> py::object {
          if (not self.inverse_state) return py::none();
          if constexpr (std::same_as<Context, ConstVectorView>) {
            return py::cpp_function(
                [fn = self.inverse_state](const Vector& state, const Vector& context) { return fn(state, context); });
          } else {
            return py::cpp_function(
                [fn = self.inverse_state](const Vector& state, const Context& context) { return fn(state, context); });
          }
        },
        [python_callable](Target& self, py::object fn) {
          if (fn.is_none()) {
            self.inverse_state = {};
            return;
          }
          auto callable      = python_callable(std::move(fn));
          self.inverse_state = [callable = std::move(callable)](ConstVectorView state,
                                                                const Context&  context) -> Vector {
            py::gil_scoped_acquire guard;
            py::object             pycontext;
            if constexpr (std::same_as<Context, ConstVectorView>) {
              pycontext = py::cast(Vector{context});
            } else {
              pycontext = py::cast(&context, py::rv_policy::reference);
            }
            py::object result = (*callable)(Vector{state}, std::move(pycontext));
            return py::cast<Vector>(py::type<Vector>()(result));
          };
        },
        R"doc(Transform model-state values back into native target values.

The callable signature is ``(model_state, context) -> native_state``.  The
context is the complete field or data object that owns the target.

.. :class:`~collections.abc.Callable`)doc");
    target.def_prop_rw(
        "inverse_jacobian",
        [](const Target& self) -> py::object {
          if (not self.inverse_jacobian) return py::none();
          if constexpr (std::same_as<Context, ConstVectorView>) {
            return py::cpp_function(
                [fn = self.inverse_jacobian](const Matrix& jacobian, const Vector& state, const Vector& context) {
                  return fn(jacobian, state, context);
                });
          } else {
            return py::cpp_function(
                [fn = self.inverse_jacobian](const Matrix& jacobian, const Vector& state, const Context& context) {
                  return fn(jacobian, state, context);
                });
          }
        },
        [python_callable](Target& self, py::object fn) {
          if (fn.is_none()) {
            self.inverse_jacobian = {};
            return;
          }
          auto callable         = python_callable(std::move(fn));
          self.inverse_jacobian = [callable = std::move(callable)](ConstMatrixView jacobian,
                                                                   ConstVectorView state,
                                                                   const Context&  context) -> Matrix {
            py::gil_scoped_acquire guard;
            py::object             pycontext;
            if constexpr (std::same_as<Context, ConstVectorView>) {
              pycontext = py::cast(Vector{context});
            } else {
              pycontext = py::cast(&context, py::rv_policy::reference);
            }
            py::object result = (*callable)(Matrix{jacobian}, Vector{state}, std::move(pycontext));
            return py::cast<Matrix>(py::type<Matrix>()(result));
          };
        },
        R"doc(Transform a native-unit Jacobian into model-state units.

The callable signature is ``(jacobian, model_state, context) -> jacobian`` and
must apply the derivative of ``inverse_state`` with respect to the model
state.

.. :class:`~collections.abc.Callable`)doc");
  };

  bind_transformations(atm, std::type_identity<AtmField>{});
  bind_transformations(surf, std::type_identity<SurfaceField>{});
  bind_transformations(subsurf, std::type_identity<SubsurfaceField>{});
  bind_transformations(line, std::type_identity<AbsorptionBands>{});
  bind_transformations(sensor, std::type_identity<ArrayOfSensorObsel>{});
  bind_transformations(error, std::type_identity<ConstVectorView>{});

  auto atm_targets =
      py::bind_vector<std::vector<Jacobian::AtmTarget>, py::rv_policy::reference_internal>(m, "ArrayOfAtmTargets");
  atm_targets.doc() = "List of atmospheric targets";
  auto surf_targets = py::bind_vector<std::vector<Jacobian::SurfaceTarget>, py::rv_policy::reference_internal>(
      m, "ArrayOfSurfaceTarget");
  surf_targets.doc()   = "List of surface targets";
  auto subsurf_targets = py::bind_vector<std::vector<Jacobian::SubsurfaceTarget>, py::rv_policy::reference_internal>(
      m, "ArrayOfSubsurfaceTarget");
  subsurf_targets.doc() = "List of subsurface targets";
  auto line_targets =
      py::bind_vector<std::vector<Jacobian::LineTarget>, py::rv_policy::reference_internal>(m, "ArrayOfLineTarget");
  line_targets.doc() = "List of line targets";
  auto sensor_targets =
      py::bind_vector<std::vector<Jacobian::SensorTarget>, py::rv_policy::reference_internal>(m, "ArrayOfSensorTarget");
  sensor_targets.doc() = "List of sensor targets";
  auto error_targets =
      py::bind_vector<std::vector<Jacobian::ErrorTarget>, py::rv_policy::reference_internal>(m, "ArrayOfErrorTarget");
  error_targets.doc() = "List of error targets";

  generic_interface(atm_targets);
  generic_interface(surf_targets);
  generic_interface(subsurf_targets);
  generic_interface(line_targets);
  generic_interface(sensor_targets);
  generic_interface(error_targets);

  py::class_<JacobianTargets> jacs(m, "JacobianTargets");
  generic_interface(jacs);
  jacs.def_rw(
      "atm", &JacobianTargets::atm, "List of atmospheric targets\n\n.. :class:`~pyarts3.arts.ArrayOfAtmTargets`");
  jacs.def_rw(
      "surf", &JacobianTargets::surf, "List of surface targets\n\n.. :class:`~pyarts3.arts.ArrayOfSurfaceTarget`");
  jacs.def_rw("subsurf",
              &JacobianTargets::subsurf,
              "List of subsurface targets\n\n.. :class:`~pyarts3.arts.ArrayOfSubsurfaceTarget`");
  jacs.def_rw("line", &JacobianTargets::line, "List of line targets\n\n.. :class:`~pyarts3.arts.ArrayOfLineTarget`");
  jacs.def_rw(
      "sensor", &JacobianTargets::sensor, "List of sensor targets\n\n.. :class:`~pyarts3.arts.ArrayOfSensorTarget`");
  jacs.def_rw(
      "error", &JacobianTargets::error, "List of error targets\n\n.. :class:`~pyarts3.arts.ArrayOfErrorTarget`");
  jacs.def("x_size", &JacobianTargets::x_size, "The size of the model state vector.");
  jacs.def("target_count", &JacobianTargets::target_count, "The number of targets added to the Jacobian.");
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize jac\n{}", e.what()));
}
}  // namespace Python

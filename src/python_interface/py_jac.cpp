#include <python_interface.h>

#include "atm.h"
#include "hpy_arts.h"
#include "jacobian.h"
#include "lbl_data.h"
#include "nanobind/stl/bind_vector.h"
#include "nanobind/stl/function.h"
#include "nanobind/stl/variant.h"
#include "surf.h"

NB_MAKE_OPAQUE(std::vector<Jacobian::AtmTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::SurfaceTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::LineTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::SensorTarget>)
NB_MAKE_OPAQUE(std::vector<Jacobian::ErrorTarget>)

namespace Python {
void py_jac(py::module_& m) try {
  py::class_<ErrorKey> errkey(m, "ErrorKey");
  errkey.doc() = "Error key";
  errkey.def_rw(
      "y_start",
      &ErrorKey::y_start,
      "Start index of target in measurement vector\n\n.. :class:`Index`");
  errkey.def_rw("y_size",
                &ErrorKey::y_size,
                "Size of target in measurement vector\n\n.. :class:`Index`");

  py::class_<Jacobian::AtmTarget> atm(m, "JacobianAtmTarget");
  atm.doc() = "Atmospheric target";
  atm.def_ro("type",
             &Jacobian::AtmTarget::type,
             "Type of target\n\n.. :class:`AtmKeyVal`");
  atm.def_ro("d",
             &Jacobian::AtmTarget::d,
             "Perturbation magnitude\n\n.. :class:`Numeric`");
  atm.def_ro("target_pos",
             &Jacobian::AtmTarget::target_pos,
             "Target position\n\n.. :class:`Index`");
  atm.def_ro("x_start",
             &Jacobian::AtmTarget::x_start,
             "Start index of target in state vector\n\n.. :class:`Index`");
  atm.def_ro("x_size",
             &Jacobian::AtmTarget::x_size,
             "Size of target in state vector\n\n.. :class:`Index`");
  str_interface(atm);

  py::class_<Jacobian::SurfaceTarget> surf(m, "JacobianSurfaceTarget");
  surf.doc() = "Surface target";
  surf.def_ro("type",
              &Jacobian::SurfaceTarget::type,
              "Type of target\n\n.. :class:`SurfaceKeyVal`");
  surf.def_ro("d",
              &Jacobian::SurfaceTarget::d,
              "Perturbation magnitude\n\n.. :class:`Numeric`");
  surf.def_ro("target_pos",
              &Jacobian::SurfaceTarget::target_pos,
              "Target position\n\n.. :class:`Index`");
  surf.def_ro("x_start",
              &Jacobian::SurfaceTarget::x_start,
              "Start index of target in state vector\n\n.. :class:`Index`");
  surf.def_ro("x_size",
              &Jacobian::SurfaceTarget::x_size,
              "Size of target in state vector\n\n.. :class:`Index`");
  str_interface(surf);

  py::class_<Jacobian::LineTarget> line(m, "JacobianLineTarget");
  line.doc() = "Line target";
  line.def_ro("type",
              &Jacobian::LineTarget::type,
              "Type of target\n\n.. :class:`LineKeyVal`");
  line.def_ro("d",
              &Jacobian::LineTarget::d,
              "Perturbation magnitude\n\n.. :class:`Numeric`");
  line.def_ro("target_pos",
              &Jacobian::LineTarget::target_pos,
              "Target position\n\n.. :class:`Index`");
  line.def_ro("x_start",
              &Jacobian::LineTarget::x_start,
              "Start index of target in state vector\n\n.. :class:`Index`");
  line.def_ro("x_size",
              &Jacobian::LineTarget::x_size,
              "Size of target in state vector\n\n.. :class:`Index`");
  str_interface(line);

  py::class_<Jacobian::SensorTarget> sensor(m, "JacobianSensorTarget");
  sensor.doc() = "Sensor target";
  sensor.def_ro("type",
                &Jacobian::SensorTarget::type,
                "Type of target\n\n.. :class:`SensorKeyVal`");
  sensor.def_ro("d",
                &Jacobian::SensorTarget::d,
                "Perturbation magnitude\n\n.. :class:`Numeric`");
  sensor.def_ro("target_pos",
                &Jacobian::SensorTarget::target_pos,
                "Target position\n\n.. :class:`Index`");
  sensor.def_ro("x_start",
                &Jacobian::SensorTarget::x_start,
                "Start index of target in state vector\n\n.. :class:`Index`");
  sensor.def_ro("x_size",
                &Jacobian::SensorTarget::x_size,
                "Size of target in state vector\n\n.. :class:`Index`");
  str_interface(sensor);

  py::class_<Jacobian::ErrorTarget> error(m, "JacobianErrorTarget");
  error.doc() = "Error target";
  error.def_ro("type",
               &Jacobian::ErrorTarget::type,
               "Type of target\n\n.. :class:`ErrorKeyVal`");
  error.def_ro("target_pos",
               &Jacobian::ErrorTarget::target_pos,
               "Target position\n\n.. :class:`Index`");
  error.def_ro("x_start",
               &Jacobian::ErrorTarget::x_start,
               "Start index of target in state vector\n\n.. :class:`Index`");
  error.def_ro("x_size",
               &Jacobian::ErrorTarget::x_size,
               "Size of target in state vector\n\n.. :class:`Index`");
  str_interface(error);

  py::bind_vector<std::vector<Jacobian::AtmTarget>>(m, "ArrayOfAtmTargets")
      .doc() = "List of atmospheric targets";
  py::bind_vector<std::vector<Jacobian::SurfaceTarget>>(m,
                                                        "ArrayOfSurfaceTarget")
      .doc() = "List of surface targets";
  py::bind_vector<std::vector<Jacobian::LineTarget>>(m, "ArrayOfLineTarget")
      .doc() = "List of line targets";
  py::bind_vector<std::vector<Jacobian::SensorTarget>>(m, "ArrayOfSensorTarget")
      .doc() = "List of sensor targets";
  py::bind_vector<std::vector<Jacobian::ErrorTarget>>(m, "ArrayOfErrorTarget")
      .doc() = "List of error targets";

  py::class_<JacobianTargets> jacs(m, "JacobianTargets");
  generic_interface(jacs);
  jacs.def_rw("atm",
              &JacobianTargets::atm,
              "List of atmospheric targets\n\n.. :class:`ArrayOfAtmTargets`");
  jacs.def_rw("surf",
              &JacobianTargets::surf,
              "List of surface targets\n\n.. :class:`ArrayOfSurfaceTarget`");
  jacs.def_rw(
      "subsurf",
      &JacobianTargets::subsurf,
      "List of subsurface targets\n\n.. :class:`ArrayOfSubsurfaceTarget`");
  jacs.def_rw("line",
              &JacobianTargets::line,
              "List of line targets\n\n.. :class:`ArrayOfLineTarget`");
  jacs.def_rw("sensor",
              &JacobianTargets::sensor,
              "List of sensor targets\n\n.. :class:`ArrayOfSensorTarget`");
  jacs.def_rw("error",
              &JacobianTargets::error,
              "List of error targets\n\n.. :class:`ArrayOfErrorTarget`");
  jacs.def("x_size",
           &JacobianTargets::x_size,
           "The size of the model state vector.");
  jacs.def("target_count",
           &JacobianTargets::target_count,
           "The number of targets added to the Jacobian.");
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize jac\n{}", e.what()));
}
}  // namespace Python
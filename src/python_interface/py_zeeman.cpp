#include <lbl.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_zeeman(py::module_& m) try {
  py::class_<lbl::zeeman::magnetic_angles>(m, "MagneticAngles")
      .def(py::init<Vector3, Vector2>())
      .PythonInterfaceBasicRepresentation(lbl::zeeman::magnetic_angles)
      .def_ro("u", &lbl::zeeman::magnetic_angles::u)
      .def_ro("v", &lbl::zeeman::magnetic_angles::v)
      .def_ro("w", &lbl::zeeman::magnetic_angles::w)
      .def_ro("sa", &lbl::zeeman::magnetic_angles::sa)
      .def_ro("ca", &lbl::zeeman::magnetic_angles::ca)
      .def_ro("sz", &lbl::zeeman::magnetic_angles::sz)
      .def_ro("cz", &lbl::zeeman::magnetic_angles::cz)
      .def_ro("H", &lbl::zeeman::magnetic_angles::H)
      .def_ro("uct", &lbl::zeeman::magnetic_angles::uct)
      .def_ro("duct", &lbl::zeeman::magnetic_angles::duct)
      .def_prop_ro("theta", &lbl::zeeman::magnetic_angles::theta)
      .def_prop_ro("eta", &lbl::zeeman::magnetic_angles::eta)
      .def_prop_ro("dtheta_du", &lbl::zeeman::magnetic_angles::dtheta_du)
      .def_prop_ro("dtheta_dv", &lbl::zeeman::magnetic_angles::dtheta_dv)
      .def_prop_ro("dtheta_dw", &lbl::zeeman::magnetic_angles::dtheta_dw)
      .def_prop_ro("deta_du", &lbl::zeeman::magnetic_angles::deta_du)
      .def_prop_ro("deta_dv", &lbl::zeeman::magnetic_angles::deta_dv)
      .def_prop_ro("deta_dw", &lbl::zeeman::magnetic_angles::deta_dw)
      .doc() = "Magnetic angles for the Zeeman effect, note that the numbers are in radians.";
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize IGRF\n", e.what()));
}
}  // namespace Python

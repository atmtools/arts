#include <lbl.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_zeeman(py::module_& m) try {
  artsclass<lbl::zeeman::magnetic_angles>(m, "MagneticAngles")
      .def(py::init<Vector3, Vector2>())
      .PythonInterfaceBasicRepresentation(lbl::zeeman::magnetic_angles)
      .def_readonly("u", &lbl::zeeman::magnetic_angles::u)
      .def_readonly("v", &lbl::zeeman::magnetic_angles::v)
      .def_readonly("w", &lbl::zeeman::magnetic_angles::w)
      .def_readonly("sa", &lbl::zeeman::magnetic_angles::sa)
      .def_readonly("ca", &lbl::zeeman::magnetic_angles::ca)
      .def_readonly("sz", &lbl::zeeman::magnetic_angles::sz)
      .def_readonly("cz", &lbl::zeeman::magnetic_angles::cz)
      .def_readonly("H", &lbl::zeeman::magnetic_angles::H)
      .def_readonly("uct", &lbl::zeeman::magnetic_angles::uct)
      .def_readonly("duct", &lbl::zeeman::magnetic_angles::duct)
      .def_property_readonly("theta", &lbl::zeeman::magnetic_angles::theta)
      .def_property_readonly("eta", &lbl::zeeman::magnetic_angles::eta)
      .def_property_readonly("dtheta_du", &lbl::zeeman::magnetic_angles::dtheta_du)
      .def_property_readonly("dtheta_dv", &lbl::zeeman::magnetic_angles::dtheta_dv)
      .def_property_readonly("dtheta_dw", &lbl::zeeman::magnetic_angles::dtheta_dw)
      .def_property_readonly("deta_du", &lbl::zeeman::magnetic_angles::deta_du)
      .def_property_readonly("deta_dv", &lbl::zeeman::magnetic_angles::deta_dv)
      .def_property_readonly("deta_dw", &lbl::zeeman::magnetic_angles::deta_dw)
      .doc() = "Magnetic angles for the Zeeman effect, note that the numbers are in radians.";
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize IGRF\n", e.what()));
}
}  // namespace Python

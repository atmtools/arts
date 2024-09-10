#include <lbl.h>

#include "py_macros.h"
#include "python_interface.h"

namespace Python {
void py_zeeman(py::module_& m) try {
  py::class_<lbl::zeeman::magnetic_angles>(m, "MagneticAngles")
      .def(py::init<Vector3, Vector2>())
      //.PythonInterfaceBasicRepresentation(lbl::zeeman::magnetic_angles)
      .def_ro("u", &lbl::zeeman::magnetic_angles::u, "Magnetic u")
      .def_ro("v", &lbl::zeeman::magnetic_angles::v, "Magnetic v")
      .def_ro("w", &lbl::zeeman::magnetic_angles::w, "Magnetic w")
      .def_ro("sa", &lbl::zeeman::magnetic_angles::sa, "Magnetic cos(asimuth)")
      .def_ro("ca", &lbl::zeeman::magnetic_angles::ca, "Magnetic sin(asimuth)")
      .def_ro("sz", &lbl::zeeman::magnetic_angles::sz, "Magnetic sin(zenith)")
      .def_ro("cz", &lbl::zeeman::magnetic_angles::cz, "Magnetic cos(zenith)")
      .def_ro("H",
              &lbl::zeeman::magnetic_angles::H,
              "Absolute magnetic field strength")
      .def_ro("uct", &lbl::zeeman::magnetic_angles::uct, "Angle constant")
      .def_ro("duct",
              &lbl::zeeman::magnetic_angles::duct,
              "Derivative of angle constant")
      .def_prop_ro(
          "theta", &lbl::zeeman::magnetic_angles::theta, "Magnetic angle theta")
      .def_prop_ro(
          "eta", &lbl::zeeman::magnetic_angles::eta, "Magnetic angle eta")
      .def_prop_ro("dtheta_du",
                   &lbl::zeeman::magnetic_angles::dtheta_du,
                   "Derivative of theta with respect to u")
      .def_prop_ro("dtheta_dv",
                   &lbl::zeeman::magnetic_angles::dtheta_dv,
                   "Derivative of theta with respect to v")
      .def_prop_ro("dtheta_dw",
                   &lbl::zeeman::magnetic_angles::dtheta_dw,
                   "Derivative of theta with respect to w")
      .def_prop_ro("deta_du",
                   &lbl::zeeman::magnetic_angles::deta_du,
                   "Derivative of eta with respect to u")
      .def_prop_ro("deta_dv",
                   &lbl::zeeman::magnetic_angles::deta_dv,
                   "Derivative of eta with respect to v")
      .def_prop_ro("deta_dw",
                   &lbl::zeeman::magnetic_angles::deta_dw,
                   "Derivative of eta with respect to w")
      .doc() =
      "Magnetic angles for the Zeeman effect, note that the numbers are in radians.";
} catch (std::exception& e) {
  throw std::runtime_error(
      var_string("DEV ERROR:\nCannot initialize IGRF\n", e.what()));
}
}  // namespace Python

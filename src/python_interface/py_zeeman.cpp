#include <arts_conversions.h>
#include <debug.h>
#include <lbl.h>
#include <nanobind/stl/array.h>
#include <nanobind/stl/string.h>

#include <format>

#include "hpy_arts.h"
#include "lbl_zeeman.h"
#include "matpack_mdspan_cdata_t.h"
#include "python_interface.h"

namespace Python {

void py_zeeman(py::module_& m) try {
  auto zee  = m.def_submodule("zeeman");
  zee.doc() = "Zeeman effect calculations";

  py::class_<lbl::zeeman::magnetic_angles> zzma(zee, "MagneticAngles");
  generic_interface(zzma);
  zzma.def(py::init<Vector3, Vector2>())
      .def_ro("u",
              &lbl::zeeman::magnetic_angles::u,
              "Magnetic u\n\n.. :class:`float`")
      .def_ro("v",
              &lbl::zeeman::magnetic_angles::v,
              "Magnetic v\n\n.. :class:`float`")
      .def_ro("w",
              &lbl::zeeman::magnetic_angles::w,
              "Magnetic w\n\n.. :class:`float`")
      .def_ro("sa",
              &lbl::zeeman::magnetic_angles::sa,
              "Magnetic cos(asimuth)\n\n.. :class:`float`")
      .def_ro("ca",
              &lbl::zeeman::magnetic_angles::ca,
              "Magnetic sin(asimuth)\n\n.. :class:`float`")
      .def_ro("sz",
              &lbl::zeeman::magnetic_angles::sz,
              "Magnetic sin(zenith)\n\n.. :class:`float`")
      .def_ro("cz",
              &lbl::zeeman::magnetic_angles::cz,
              "Magnetic cos(zenith)\n\n.. :class:`float`")
      .def_ro("H",
              &lbl::zeeman::magnetic_angles::H,
              "Absolute magnetic field strength\n\n.. :class:`float`")
      .def_ro("uct",
              &lbl::zeeman::magnetic_angles::uct,
              "Angle constant\n\n.. :class:`float`")
      .def_ro("duct",
              &lbl::zeeman::magnetic_angles::duct,
              "Derivative of angle constant\n\n.. :class:`float`")
      .def_prop_ro("theta",
                   &lbl::zeeman::magnetic_angles::theta,
                   "Magnetic angle theta\n\n.. :class:`float`")
      .def_prop_ro("eta",
                   &lbl::zeeman::magnetic_angles::eta,
                   "Magnetic angle eta\n\n.. :class:`float`")
      .def_prop_ro("dtheta_du",
                   &lbl::zeeman::magnetic_angles::dtheta_du,
                   "Derivative of theta with respect to u\n\n.. :class:`float`")
      .def_prop_ro("dtheta_dv",
                   &lbl::zeeman::magnetic_angles::dtheta_dv,
                   "Derivative of theta with respect to v\n\n.. :class:`float`")
      .def_prop_ro("dtheta_dw",
                   &lbl::zeeman::magnetic_angles::dtheta_dw,
                   "Derivative of theta with respect to w\n\n.. :class:`float`")
      .def_prop_ro("deta_du",
                   &lbl::zeeman::magnetic_angles::deta_du,
                   "Derivative of eta with respect to u\n\n.. :class:`float`")
      .def_prop_ro("deta_dv",
                   &lbl::zeeman::magnetic_angles::deta_dv,
                   "Derivative of eta with respect to v\n\n.. :class:`float`")
      .def_prop_ro("deta_dw",
                   &lbl::zeeman::magnetic_angles::deta_dw,
                   "Derivative of eta with respect to w\n\n.. :class:`float`")
      .def("__repr__",
           [](const lbl::zeeman::magnetic_angles& x) {
             return std::format(
                 "MagneticAngles(B = ({}, {}, {}), LOS = ({}, {})):\n"
                 "  theta: {} deg\n"
                 "  eta:   {} deg",
                 x.u,
                 x.v,
                 x.w,
                 Conversion::rad2deg(std::atan2(x.sz, x.cz)),
                 Conversion::rad2deg(std::atan2(x.sa, x.ca)),
                 Conversion::rad2deg(x.theta()),
                 Conversion::rad2deg(x.eta()));
           })
      .doc() =
      "Magnetic angles for the Zeeman effect, note that the numbers are in radians.";

  zzma.def_prop_ro(
      "k",
      [](const lbl::zeeman::magnetic_angles& x) -> Vector3 {
        return {x.sz * x.sa, x.sz * x.ca, x.cz};
      },
      "Propagation vector\n\n.. :class:`~pyarts3.arts.Vector3`");

  zzma.def_prop_ro(
      "e1",
      [](const lbl::zeeman::magnetic_angles& x) -> Vector3 {
        return {-x.cz * x.sa, -x.cz * x.ca, x.sz};
      },
      "Propagation vector upwards projection\n\n.. :class:`~pyarts3.arts.Vector3`");

  zzma.def_prop_ro(
      "e2",
      [](const lbl::zeeman::magnetic_angles& x) -> Vector3 {
        return {x.ca, -x.sa, 0};
      },
      "Propagation vector rightward projection\n\n.. :class:`~pyarts3.arts.Vector3`");

  zzma.def_prop_ro(
      "B_projected",
      [](const lbl::zeeman::magnetic_angles& x) -> Vector3 {
        Vector3 B{x.u, x.v, x.w};
        Vector3 k{x.sz * x.sa, x.sz * x.ca, x.cz};

        return B - proj(B, k);
      },
      "Projection of the magnetic field onto the plane perpendicular to the line of sight.\n\n.. :class:`~pyarts3.arts.Vector3`");

  zee.def("norm_view",
          &lbl::zeeman::norm_view,
          "Normalized view matrix for the Zeeman effect",
          py::arg("polarization"),
          py::arg("magnetic"),
          py::arg("line_of_sight"));

  zee.def(
      "dnorm_view_du",
      &lbl::zeeman::dnorm_view_du,
      "Derivative of the normalized view matrix for the Zeeman effect with respect to u",
      py::arg("polarization"),
      py::arg("magnetic"),
      py::arg("line_of_sight"));

  zee.def(
      "dnorm_view_dv",
      &lbl::zeeman::dnorm_view_dv,
      "Derivative of the normalized view matrix for the Zeeman effect with respect to v",
      py::arg("polarization"),
      py::arg("magnetic"),
      py::arg("line_of_sight"));

  zee.def(
      "dnorm_view_dw",
      &lbl::zeeman::dnorm_view_dw,
      "Derivative of the normalized view matrix for the Zeeman effect with respect to w",
      py::arg("polarization"),
      py::arg("magnetic"),
      py::arg("line_of_sight"));
} catch (std::exception& e) {
  throw std::runtime_error(
      std::format("DEV ERROR:\nCannot initialize Zeeman\n{}", e.what()));
}
}  // namespace Python

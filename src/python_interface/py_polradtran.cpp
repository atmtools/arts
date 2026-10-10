#include <polradtran.h>
#include <python_interface.h>

#include "hpy_arts.h"

namespace Python {
void py_polradtran(py::module_& m) try {
  auto pr = m.def_submodule("polradtran", R"(What Evans' RT3 and RT4 (polradtran) solvers share.

The quadrature of their streams and the grounds both have, used by
:mod:`~pyarts3.arts.rt3` and :mod:`~pyarts3.arts.rt4`.  See :doc:`dev.rt3`
and :doc:`dev.rt4` for their conventions.
)");

  py::enum_<polradtran::quadrature_type>(
      pr, "QuadratureType", "The quadrature rules of RT3's and RT4's streams on one hemisphere")
      .value("gauss",
             polradtran::quadrature_type::gauss,
             "'G' (RT3's 'E' with extra angles): positive half of a 2*nmu-point Gauss-Legendre rule on [-1, 1]")
      .value("double_gauss", polradtran::quadrature_type::double_gauss, "'D': nmu-point Gauss-Legendre rule on [0, 1]")
      .value("lobatto",
             polradtran::quadrature_type::lobatto,
             "'L': positive half of a 2*nmu-point Lobatto rule on [-1, 1]; includes mu = 1");

  py::class_<polradtran::quadrature>(pr, "Quadrature")
      .def_ro("mu", &polradtran::quadrature::mu, "Ascending nodes in (0, 1]\n\n.. :class:`~pyarts3.arts.Vector`")
      .def_ro("weights",
              &polradtran::quadrature::weights,
              "Weights for the integral over mu in [0, 1], summing to 1 (no 2 pi azimuth factor)\n\n.. "
              ":class:`~pyarts3.arts.Vector`")
      .doc() = "The streams of one hemisphere";

  pr.def("get_quadrature",
         &polradtran::get_quadrature,
         "nmu"_a,
         "type"_a,
         "The streams of RT3 and RT4, nmu >= 1 nodes per hemisphere: the positive half of ARTS's 2 nmu-point rule, "
         "as the solvers use them.",
         py::call_guard<py::gil_scoped_release>());

  py::class_<polradtran::lambertian_surface>(pr, "LambertianSurface")
      .def(
          "__init__",
          [](polradtran::lambertian_surface* s, Numeric albedo) {
            new (s) polradtran::lambertian_surface{.albedo = albedo};
          },
          "albedo"_a = 0.0)
      .def_rw("albedo", &polradtran::lambertian_surface::albedo, "Albedo A\n\n.. :class:`float`")
      .doc() =
      "'L': reflection 2 A mu_j w_j into every stream in the m = 0 mode, I to I only; emission (1 - A) B in I "
      "(in RT3 only with thermal set); RT3 also reflects the direct beam, A F / pi.  Energy is conserved on the "
      "streams only for double_gauss quadrature.";

  py::class_<polradtran::fresnel_surface>(pr, "FresnelSurface")
      .def(
          "__init__",
          [](polradtran::fresnel_surface* s, Complex n) { new (s) polradtran::fresnel_surface{.refractive_index = n}; },
          "refractive_index"_a)
      .def_rw("refractive_index",
              &polradtran::fresnel_surface::refractive_index,
              "Complex refractive index of the surface (medium above has index 1)\n\n.. :class:`complex`")
      .doc() =
      "'F': specular Fresnel reflection with :math:`R1 = (|r_v|^2 + |r_h|^2) / 2, R2 = (|r_v|^2 - |r_h|^2) / 2, "
      "R3 = Re(r_v r_h*), R4 = Im(r_v r_h*)` (RT4 has the [I, Q] block); emission "
      ":math:`[(1 - R1) B, -R2 B, 0, 0]`, also when RT3's thermal is false.  RT3 does not allow it with a direct "
      "beam.";
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize polradtran\n{}", e.what()));
}
}  // namespace Python

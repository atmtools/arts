#include <python_interface.h>

#include "core/functional/functional_gravity.h"
#include "core/util/planet_data.h"
#include "hpy_arts.h"
#include "hpy_numpy.h"

namespace Python {
#define ADD_PLANET(PLANET)                                                         \
  auto PLANET              = planets.def_submodule(#PLANET, #PLANET " constants"); \
  PLANET.attr("ellipsoid") = Vector2{Body::PLANET::a, Body::PLANET::b};            \
  PLANET.attr("GM")        = Body::PLANET::GM;                                     \
  PLANET.attr("day")       = Body::PLANET::day;                                    \
  PLANET.def("gravity_operator",                                                   \
             &EllipsoidGravity::PLANET,                                            \
             "Get " #PLANET "'s gravity object.\n\n.. :class:`~pyarts3.arts.EllipsoidGravity`");

void py_planets(py::module_& m) try {
  py::class_<EllipsoidGravity> eg(m, "EllipsoidGravity");
  eg.doc() = "Ellipsoid gravity operator - allows calculation of gravity on ellipsoidal bodies";
  generic_interface(eg);
  eg.def(py::init<>())
      .def_rw("GM",
              &EllipsoidGravity::GM,
              "Standard gravitational parameter in m^3/s^2.\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("a", &EllipsoidGravity::a, "Semi-major axis in m.\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def_rw("b", &EllipsoidGravity::b, "Semi-minor axis in m.\n\n.. :class:`~pyarts3.arts.Numeric`")
      .def(
          "__call__",
          [](const EllipsoidGravity& self, py::object h, py::object lat, py::object lon) {
            return vectorize([&self](Numeric h, Numeric lat, Numeric lon) { return self(h, lat, lon); }, h, lat, lon);
          },
          "h"_a,
          "lat"_a,
          "lon"_a,
          "Calculate gravity at a given position.")
      .def_static("Mercury", &EllipsoidGravity::Mercury, "Mercury ellipsoid gravity operator")
      .def_static("Venus", &EllipsoidGravity::Venus, "Venus ellipsoid gravity operator")
      .def_static("Earth", &EllipsoidGravity::Earth, "Earth ellipsoid gravity operator")
      .def_static("Moon", &EllipsoidGravity::Moon, "Moon ellipsoid gravity operator")
      .def_static("Mars", &EllipsoidGravity::Mars, "Mars ellipsoid gravity operator")
      .def_static("Jupiter", &EllipsoidGravity::Jupiter, "Jupiter ellipsoid gravity operator")
      .def_static("Saturn", &EllipsoidGravity::Saturn, "Saturn ellipsoid gravity operator");
  auto planets = m.def_submodule("planets", "Planetary constants");

  ADD_PLANET(Earth);
  ADD_PLANET(Jupiter);
  ADD_PLANET(Mars);
  ADD_PLANET(Moon);
  ADD_PLANET(Mercury);
  ADD_PLANET(Venus);
  ADD_PLANET(Saturn);
} catch (std::exception& e) {
  throw std::runtime_error(std::format("DEV ERROR:\nCannot initialize planets\n{}", e.what()));
}

#undef ADD_PLANET
}  // namespace Python